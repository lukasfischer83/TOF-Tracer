
push!(LOAD_PATH, pwd())
import PyPlot
include("plotRunsCloud11.jl")

logscale = false

plotCompounds = [
"C10H16"
#"C10H16O2"
#"C10H16O3"
#"C8H12O6"
#"C10H14O5"
#"C10H16O5"
"C19H28O5"
"C19H30O8"
"C19H28O9"
#"C18H28O10"
]


(t1,h1) = readdlm("/data/CLOUD11/running1/sortedByPrecursor/1802-pure-apinene/results/ptr3traces.csv", '\t', header=true)
time1 = Dates.unix2datetime.(t1[:,1])


plotindices1 = findin(h1,plotCompounds)


#lbl1 = reshape(h1[plotindices1],(1,length(h1[plotindices1])))
lbl1 = h1[plotindices1]

PyPlot.figure()
if logscale
    PyPlot.semilogy(time1, t1[:,plotindices1])
else
    PyPlot.plot(time1, t1[:,plotindices1])
end
PyPlot.legend(lbl1)
plotStages(time1[1], time1[end])
