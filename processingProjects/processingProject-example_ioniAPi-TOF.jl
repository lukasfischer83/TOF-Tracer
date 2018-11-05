#cd("home/user/path/TOF-Tracer")
include("$(pwd())/startupAPi.jl")

fp = "./ExampleFiles/ioniAPiTOFDATA/20171027/" # All files in this path will be processed

filefilterRegexp = r"\.h5$"

rf = "./ExampleFiles/ioniAPiTOFDATA/20171027/Data_2017_10_27-13_50_59.h5"  # The mass scale from this file defines the mass scale of all

masslist = MasslistFunctions.loadMasslist("./ExampleFiles/MASSLISTS/masslist_20171129.csv")#masslist_20170907.csv")
cr = [37 55 73] # compounds for mass axis calibration

# alternatively: use an auto generated masslist
# masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions = createMassList(C=0:20, O=0:20, N=0:1, allowRadicals=false) #
(masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
s = (masslistMasses.>26) & ( masslistMasses.<2000)
masslistMasses = masslistMasses[s]
masslistCompositions = masslistCompositions[s,:]

####################### INFO ##########################################################################
###       if the lower mass range of IONICON ioniTOF is set to values <25, the parameter:           ###
###           safetyMarginMasses = 10  (can be found in MultipeakFunctions)                         ###
###       should be set to smaller values, e.g. <4 in case of an "BoundsError".                     ###

####################### END OF SETTINGS ###############################################################

####################### Processing sequence ###########################################################

correctMassScaleAndExtractSumSpecAPi(
    fp,
    masslistMasses,
    masslistElements,
    masslistElementsMasses,
    masslistCompositions,
    rf,
    cr,
    filefilterRegexp=filefilterRegexp,
    recalibInterval = 600,
    onlyUseAverages = true,
    plotControlMass = true,
    firstNFiles=0,
    lastNFiles = 0,
    testRangeStart = 95.0, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
    testRangeEnd = 98.0,
    binWidth = 10
    )

baselineAndPeakshape(                   # recommended values for:
    fp,                                 # PTR3 # ioniAPi-TOF
    peakshapeRegions=4,                 # 8    # <10
    peakshapeQuantileValue = 0.1,       # 0.1  # 0.1, in case of low counting statistics use a higher value, but higher value increases the peak width
    peakfindingNoiseThresholdValue = 25,# 25   # [1<x<25], in APi-TOF mode (without CI-source) lower counting statistics, # @CLOUD: use 11 or more for BEAM experiments
    peakWindowWidth = 80,               # 200  # <90, ioniAPi-TOF has more than twice the mass range of the PTR3, e. g. m/z 15 to 2000 Th
    )

mtrx = deconvolute(
    fp,
    calcTransposed = false
    )
