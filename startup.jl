if (nprocs() < 2)
    println("Adding process for file precaching!")
    addprocs(1)
end
@everywhere scriptspath =  pwd()

include("$(scriptspath)/manualMassLibrary.jl")
include("$(scriptspath)/combinedMassScaleAndExtractSumSpec.jl")
include("$(scriptspath)/peakShape.jl")
include("$(scriptspath)/deconvolutionMatrix.jl")
