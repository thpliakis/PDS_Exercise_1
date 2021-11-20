#using Pkg

#Pkg.add("PlotlyJs")
#Pkg.add("CSV")

using PlotlyJS, CSV

s = CSV.File("Sequential-Triang-Count/seq_runtimes.csv")
p = CSV.File("Pthreads-Triang-Count/pthreads_runtimes.csv")
c = CSV.File("Cilk-Triang-Count/cilk_runtimes.csv")
o = CSV.File("OpenMP-Triang-Count/openMP_runtimes.csv")


graph = s.graph[1];
plot([
    scatter(x=s.threads_num[s.graph .== graph], y=s.times[s.graph.== graph],mode="markers+lines", name="sequential"),
    scatter(x=p.threads_num[p.graph .== graph], y=p.times[p.graph.== graph],mode="markers+lines", name="pthreads"),
    scatter(x=c.threads_num[c.graph .== graph], y=c.times[c.graph.== graph],mode="markers+lines", name="cilk"),
    scatter(x=o.threads_num[o.graph .== graph], y=o.times[o.graph.== graph],mode="markers+lines", name="openMP")
],  Layout(title=graph*" #Triangle"*string(unique(s.triang_num[s.graph .== graph]))*" #Nodes"*string(unique(s.node_num[s.graph .== graph])),yaxis_title="Times", xaxis_title="Number of Threads")
)
#=
graph = s.graph[2];
plot([
    scatter(x=s.threads_num[s.graph .== graph], y=s.times[s.graph.== graph],mode="markers+lines", name="sequential"),
    scatter(x=p.threads_num[p.graph .== graph], y=p.times[p.graph.== graph],mode="markers+lines", name="pthreads"),
    scatter(x=c.threads_num[c.graph .== graph], y=c.times[c.graph.== graph],mode="markers+lines", name="cilk"),
    scatter(x=o.threads_num[o.graph .== graph], y=o.times[o.graph.== graph],mode="markers+lines", name="openMP")
],Layout(title=graph*" #Triangle"*string(unique(s.triang_num[s.graph .== graph]))*" #Nodes"*string(unique(s.node_num[s.graph .== graph])),yaxis_title="Times", xaxis_title="Number of Threads")
)


graph = s.graph[3];
plot([
    scatter(x=s.threads_num[s.graph .== graph], y=s.times[s.graph.== graph],mode="markers+lines", name="sequential"),
    scatter(x=p.threads_num[p.graph .== graph], y=p.times[p.graph.== graph],mode="markers+lines", name="pthreads"),
    scatter(x=c.threads_num[c.graph .== graph], y=c.times[c.graph.== graph],mode="markers+lines", name="cilk"),
    scatter(x=o.threads_num[o.graph .== graph], y=o.times[o.graph.== graph],mode="markers+lines", name="openMP")
],Layout(title=graph*" #Triangle"*string(unique(s.triang_num[s.graph .== graph]))*" #Nodes"*string(unique(s.node_num[s.graph .== graph])),yaxis_title="Times", xaxis_title="Number of Threads")
)

graph = s.graph[4];
plot([
    scatter(x=s.threads_num[s.graph .== graph], y=s.times[s.graph.== graph],mode="markers+lines", name="sequential"),
    scatter(x=p.threads_num[p.graph .== graph], y=p.times[p.graph.== graph],mode="markers+lines", name="pthreads"),
    scatter(x=c.threads_num[c.graph .== graph], y=c.times[c.graph.== graph],mode="markers+lines", name="cilk"),
    scatter(x=o.threads_num[o.graph .== graph], y=o.times[o.graph.== graph],mode="markers+lines", name="openMP")
],Layout(title=graph*" #Triangle"*string(unique(s.triang_num[s.graph .== graph]))*" #Nodes"*string(unique(s.node_num[s.graph .== graph])),yaxis_title="Times", xaxis_title="Number of Threads")
)

graph = s.graph[5];
plot([
    scatter(x=s.threads_num[s.graph .== graph], y=s.times[s.graph.== graph],mode="markers+lines", name="sequential"),
    scatter(x=p.threads_num[p.graph .== graph], y=p.times[p.graph.== graph],mode="markers+lines", name="pthreads"),
    scatter(x=c.threads_num[c.graph .== graph], y=c.times[c.graph.== graph],mode="markers+lines", name="cilk"),
    scatter(x=o.threads_num[o.graph .== graph], y=o.times[o.graph.== graph],mode="markers+lines", name="openMP")
],Layout(title=graph*" #Triangle"*string(unique(s.triang_num[s.graph .== graph]))*" #Nodes"*string(unique(s.node_num[s.graph .== graph])),yaxis_title="Times", xaxis_title="Number of Threads")
)
=#
