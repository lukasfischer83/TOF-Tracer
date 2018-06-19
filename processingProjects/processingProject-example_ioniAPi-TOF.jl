cd("/home/markus/github/TOF-Tracer")
include("$(pwd())/startup.jl")

fp = "./ExampleFiles/ioniAPiTOFDATA/" # All files in this path will be processed
filefilterRegexp = r"\.h5$"
rf = "./ExampleFiles/ioniAPiTOFDATA/Data_2017_10_27-13_50_59.h5"  # The mass scale from this file defines the mass scale of all
masslist = MasslistFunctions.loadMasslist("./ExampleFiles/MASSLISTS/exampleMassList.csv")
cr = [37 55 73]

# alternatively: use an auto generated masslist
# masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions = createMassList(C=0:20, O=0:20, N=0:1, allowRadicals=false) #
(masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
s = (masslistMasses.>17) & ( masslistMasses.<1200)
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
    onlyUseAverages = true,
    plotControlMass = true,
    firstNFiles=0,
    lastNFiles = 0,
    testRangeStart = 52.0, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
    testRangeEnd = 56.5,
    )

baselineAndPeakshape(                   # recommended values for:
    fp,                                 # PTR3 # ioniAPi-TOF
    peakshapeRegions=6,                 # 8    # <10
    peakshapeQuantileValue = 0.3,       # 0.1  # 0.1, in case of low counting statistics use a higher value, but higher value increases the peak width
    peakfindingNoiseThresholdValue = 4, # 25   # [1<x<10], in APi-TOF mode (without CI-source) lower counting statistics
    peakWindowWidth = 80                # 200  # <90, ioniAPi-TOF has more than twice the mass range of the PTR3, e. g. m/z 15 to 2000 Th
    )

mtrx = deconvolute(
    fp,
    calcTransposed = true
    )
