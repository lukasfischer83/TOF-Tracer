include("$(pwd())/startup.jl")

fp = "./ExampleFiles/TOFDATA/" # All files in this path will be processed
filefilterRegexp = r"\.h5$"
rf = "./ExampleFiles/TOFDATA/2016-10-02-19h15m05.h5"  # The mass scale from this file defines the mass scale of all
masslist = MasslistFunctions.loadMasslist("./ExampleFiles/MASSLISTS/exampleMassList.csv")
cr = [59 391]

# alternatively: use an auto generated masslist
# masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions = createMassList(C=0:20, O=0:20, N=0:1, allowRadicals=false) #
(masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
s = (masslistMasses.>0) & ( masslistMasses.<600)
masslistMasses = masslistMasses[s]
masslistCompositions = masslistCompositions[s,:]

####################### END OF SETTINGS ###############################################################

####################### Processing sequence ###########################################################

correctMassScaleAndExtractSumSpec(
    fp,
    masslistMasses,
    masslistElements,
    masslistElementsMasses,
    masslistCompositions,
    rf,
    cr,
    filefilterRegexp=filefilterRegexp,
    onlyUseAverages = false,
    plotControlMass = true,
    firstNFiles=0,
    lastNFiles = 0,
    )

baselineAndPeakshape(
    fp,
    peakshapeRegions=8,
    peakshapeQuantileValue = 0.1,
    peakfindingNoiseThresholdValue = 25
    )

mtrx = deconvolute(
    fp
    )
