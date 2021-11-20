#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>
#include "mmio.h"
#include <pthread.h>

void coo2csc(
    uint32_t       * const row,       /*!< CSC row start indices */
    uint32_t       * const col,       /*!< CSC column indices */
    uint32_t const * const row_coo,   /*!< COO row indices */
    uint32_t const * const col_coo,   /*!< COO column indices */
    uint32_t const         nnz,       /*!< Number of nonzero elements */
    uint32_t const         n,         /*!< Number of rows/columns */
    uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
    );

//Function to count common elements in 2 arrays
int count_Common(uint32_t arr1[], int len1, uint32_t arr2[], int len2,uint32_t *c3, int nodnum,int col);

void *vertexWiseTriangleCounts(void *arg);

typedef struct {    /* Used as argument to thread_start() */
    pthread_t thread_id;        /* ID returned by pthread_create() */
    int       thread_num;       /* Application-defined thread # */
    int loop_start;		 /* Because we want to brake the 1 loop */
    int loop_finish;		 /* in as many as threads exist */
    uint32_t *rowI;		 // rowI and rowJ are used to temporary save the row indexes
    uint32_t *rowJ;
}thread_info;

//Global variables in order all of the threads to have access
int activeThreads = 0, maxthreads;
int arrSize;
uint32_t * vector;
uint32_t *c3;
uint32_t *csc_row;
uint32_t *csc_col;
int M;

pthread_attr_t attr;
pthread_mutex_t mux;

int main(int argc, char *argv[]){

    int ret_code;
    MM_typecode matcode;
    FILE *f,*f2;
    int M, N, nnz;
    int i;
    double *val;
    uint32_t num=0;
    int s;
    void *res;

    uint32_t isOneBased = 0;
    uint32_t *I,*J;

    int workers[6] = {1,2,5,8,10,20};
    //maxthreads = 1;
    //maxthreads = 2;
    //maxthreads = 5maxthreads;
    //maxehreads = 8;
    //maxthreads = 10;
    //maxthreads = 20;

    /* variables to hold execution time */
    struct timeval startwtime, endwtime;
    double par_time;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL){
            printf("Could not process Matrix Market banner.\n");
            exit(1);
        }
    }

    if ((f2 = fopen("pthreads_runtimes.csv","a")) == NULL){
        printf("Could not process data writing file.\n");
        exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nnz)) !=0)
        exit(1);

    printf("nnz = %d\n",nnz);
    printf("M = %d\n",M);
    printf("N = %d\n",N);


    /* reseve memory for COO format matrices */

    I = (uint32_t *) malloc(nnz * sizeof(int));
    J = (uint32_t *) malloc(nnz * sizeof(int));
    val = (double *) malloc(nnz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    /* Replace missing val column with 1s and change the fscanf to match patter matrices*/

    if (!mm_is_pattern(matcode))
    {
        for (i=0; i<nnz; i++)
        {
            fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
            I[i]--;  /* adjust from 1-based to 0-based */
            J[i]--;
        }
    }
    else
    {
        for (i=0; i<nnz; i++)
        {
            fscanf(f, "%d %d\n", &I[i], &J[i]);
            val[i]=1;
            I[i]--;  /* adjust from 1-based to 0-based */
            J[i]--;
        }
    }

    /* Allocate space for csc format */
    csc_row = (uint32_t *)malloc(nnz   * sizeof(uint32_t));
    csc_col = (uint32_t *)malloc((M+1) * sizeof(uint32_t));
    vector = (uint32_t *)malloc(M   * sizeof(uint32_t));
    for (uint32_t l = 0; l < M; l++) vector[l] = 0;
    c3 = (uint32_t *)malloc(M   * sizeof(uint32_t));
    for (uint32_t l = 0; l < M; l++) c3[l] = 0;

    /* initialize arguments to pass to pthread_create*/
    pthread_t t;
    pthread_attr_init( &attr );
    pthread_mutex_init (&mux, NULL);
    thread_info *tinfo;
    maxthreads = workers[0];
    tinfo = calloc(maxthreads, sizeof(*tinfo));

    coo2csc(csc_row, csc_col, I, J, nnz, M, isOneBased);
    //mm_write_banner(stdout, matcode);                                           
    //mm_write_mtx_crd_size(stdout, M, N, nnz);                                    

    for(int w=0; w<6; w++){
        maxthreads = workers[w];
        tinfo = realloc(tinfo, maxthreads*sizeof(*tinfo));
        //Calculate the nodes each thread will check for triangles
        int col_per_thread = M/maxthreads;
        //printf("col_per_thread = %d\n",col_per_thread);

        gettimeofday (&startwtime, NULL);
        for (int tnum = 0; tnum < maxthreads; tnum++) {
            //init thread info
            tinfo[tnum].thread_num = tnum + 1;
            tinfo[tnum].loop_start = tnum*col_per_thread;
            tinfo[tnum].loop_finish = (tnum+1)*col_per_thread;
            if(tnum == (maxthreads - 1))
                tinfo[tnum].loop_finish = M;

            /* The pthread_create() call stores the thread ID into
               corresponding element of tinfo[] */
            tinfo[tnum].rowI = (uint32_t *)malloc(M   * sizeof(uint32_t));
            tinfo[tnum].rowJ = (uint32_t *)malloc(M  * sizeof(uint32_t));

            // Create each thread
            s = pthread_create(&tinfo[tnum].thread_id, &attr,vertexWiseTriangleCounts, &tinfo[tnum]);
            if (s != 0){
                printf("error in pthread_create\n");
                exit(1);
            }
        }
        // Join all threads
        for (int tnum = 0; tnum < maxthreads; tnum++) {
            s = pthread_join(tinfo[tnum].thread_id, &res);
            if (s != 0)
                printf("Didn't Joined with thread \n");

            //printf("Joined with thread %d; returned value was %s\n",tinfo[tnum].thread_num, (char *) res);
        
        }
        gettimeofday (&endwtime, NULL);

        for(int col = 0; col < M; col++)
            num += vector[col];

        for(int col = 0; col < M; col++){
            //printf("%d\n",c3[col]);
            c3[col] += vector[col];
        }


        par_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);

        /* print execution time */
        printf("Counting triangles wall clock time: %f sec\n", par_time);
        printf("Workers : %d\n", workers[w]);

        printf("triang_number = %d\n",num);
        fprintf(f2,"%s,%s,%d,%d,%f,%d\n",argv[0],argv[1],M,maxthreads,par_time,num);
        num = 0;
        for (uint32_t l = 0; l < M; l++) vector[l] = 0;
        for (uint32_t l = 0; l < M; l++) c3[l] = 0;

    }
    free(res);      /* Free memory allocated by thread */
    free(tinfo);
    pthread_attr_destroy(&attr);

    /* cleanup variables */
    if (f != stdin) fclose(f);
    if (f != stdin) fclose(f2);
    free(I);
    free(J);
    free(val);
    free( csc_row );
    free( csc_col );

}


void *vertexWiseTriangleCounts(void *arg){

    thread_info *tinfo = arg;
    int col_start, col_start2;
    int col_end, col_end2;
    int temp;
    int col2;
    uint32_t num = 0;
    uint32_t triang_num = 0;


    // for every column-node
    for(int col=tinfo->loop_start; col < tinfo->loop_finish; col++){
        triang_num = 0;
        col_start = csc_col[col];
        col_end   = csc_col[col+1];

        temp = col_start;
        
        // save in rowI the row indexes which is the nodes that column-node has an edge with
        for(int l = 0; l < (col_end-col_start); l++){
            tinfo->rowI[l] =  csc_row[temp];
            temp++;
        }

        //Go through only nodes that have an edge with column-node
        /* Go and check how many common nodes column-node has with column2-node(rowI[col2])
           which is the number of triangles that both nodes participate */
        for(col2 = 0; col2 < (col_end-col_start); col2++){
            col_start2 = csc_col[tinfo->rowI[col2]];
            col_end2   = csc_col[tinfo->rowI[col2]+1];
            temp = col_start2;

            // save in rowJ the row indexes which is the nodes that column2-node has an edge with
            for(int l = 0; l < (col_end2-col_start2); l++){
                tinfo->rowJ[l] =  csc_row[temp];

                temp++;
            }
            /* cout_Common counts the number of common nodes between column-node and column2-node
               which is actually the number of triangles this 2 nodes participate
               !!!!!!IMPORTANT : we use c3 to count for each node diffently in how many triangles it participates */
            //num = number of common nodes between nod-col and node-col2
            //which is actually the number of triangles
            //num = number of common nodes between nod-col and node-col2
            num = count_Common(tinfo->rowI,(col_end-col_start),tinfo->rowJ,(col_end2-col_start2),c3,tinfo->rowI[col2], col);
            // sum up all the triangles for column-node
            triang_num += num;
        }
        // Only 1 thread can access at a time because vector[] is global
        pthread_mutex_lock (&mux);
        vector[col] = triang_num;
        pthread_mutex_unlock (&mux);
    }
    free(tinfo->rowJ);
    free(tinfo->rowI);

}

int count_Common(uint32_t arr1[], int len1, uint32_t arr2[], int len2,uint32_t *c3, int nodnum,int col) {

    int i=0,j=0,count=0;

    while(len1 > i && len2 > j){

        if (arr1[i] < arr2[j]) {
            i++;

        }else if(arr2[j] < arr1[i]){
            j++;

        } else {
            count++;
            // Only 1 thread can access at a time because vector[] is global
            pthread_mutex_lock (&mux);
            // Add 1 if the nodes exists in more than 1 triangle
            c3[arr1[i]]++;
            c3[nodnum]++;
            c3[col]++;
            pthread_mutex_unlock (&mux);

            i++;
            j++;
        }
    }
    return count;

}

void coo2csc(
        uint32_t       * const row,       /*!< CSC row start indices */
        uint32_t       * const col,       /*!< CSC column indices */
        uint32_t const * const row_coo,   /*!< COO row indices */
        uint32_t const * const col_coo,   /*!< COO column indices */
        uint32_t const         nnz,       /*!< Number of nonzero elements */
        uint32_t const         n,         /*!< Number of rows/columns */
        uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
        ) {

    // ----- cannot assume that input is already 0!
    for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


    // ----- find the correct column sizes
    for (uint32_t l = 0; l < nnz; l++)
        col[col_coo[l] - isOneBased]++;


    // ----- cumulative sum
    for (uint32_t i = 0, cumsum = 0; i < n; i++) {
        uint32_t temp = col[i];
        col[i] = cumsum;
        cumsum += temp;
    }
    col[n] = nnz;
    // ----- copy the row indices to the correct place
    for (uint32_t l = 0; l < nnz; l++) {
        uint32_t col_l;
        col_l = col_coo[l] - isOneBased;

        uint32_t dst = col[col_l];
        row[dst] = row_coo[l] - isOneBased;

        col[col_l]++;
    }
    // ----- revert the column pointers
    for (uint32_t i = 0, last = 0; i < n; i++) {
        uint32_t temp = col[i];
        col[i] = last;
        last = temp;
    }

}
