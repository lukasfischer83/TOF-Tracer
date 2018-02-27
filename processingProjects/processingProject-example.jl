@everywhere scriptspath =  "/home/lukas/ownCloud/documents/UNI/Auswertung/JULIA/"
@everywhere include("$(scriptspath)/manualMassLibrary.jl")
@everywhere include("$(scriptspath)/combinedMassScaleAndExtractSumSpec.jl")
@everywhere include("$(scriptspath)/peakShape.jl")
@everywhere include("$(scriptspath)/deconvolutionMatrix.jl")


fp = abspath("/data/CLOUDX/running1/cleaned/") # All files in this path will be processed
filefilterRegexp = r"\.h5$"
rf = "/data/CLOUDX/running1/cleaned/2015-10-01-06h43m54s.h5"  # The mass scale from this file defines the mass scale of all
masslist = loadMasslist("/home/lukas/ownCloud/documents/UNI/Auswertung/CLOUDX/peaklists/apinene-low-nox-benjamin-reviewed.csv")
cr = [59 391]

#masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions = createMassList(C=0:20, O=0:20, N=0:1, allowRadicals=false) #
(masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
s = (masslistMasses.>0) & ( masslistMasses.<600)
masslistMasses = masslistMasses[s]
masslistCompositions = masslistCompositions[s,:]

####################### END OF SETTINGS ###############################################################

####################### Processing sequence ###########################################################

correctMassScaleAndExtractSumSpec(fp, masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions, rf, cr, filefilterRegexp=filefilterRegexp, onlyUseAverages = false, plotControlMass = false, firstNFiles=0, lastNFiles = 0, binWidth=4)
#include("$(scriptspath)/extractInstrumentParameters.jl")
baselineAndPeakshape(fp);
deconvolute(fp, binWidth=4);
