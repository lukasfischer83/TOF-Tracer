@everywhere scriptspath =  "/media/wiebke/Samsung_T3/Tof-Tracer-18-03/"
@everywhere include("$(scriptspath)manualMassLibrary.jl")
@everywhere include("$(scriptspath)/CloudIncludes/CLOUD12Functions.jl")

import ResultFileFunctions
import CalibrationFunctions
import ExportFunctions
import MasslistFunctions
import InterpolationFunctions
import PyPlot

timedelay = Dates.Hour(0) # CLOUD12
#timedelay = Dates.Hour(1) # CLOUDX

measurementResultFile = "/media/wiebke/Samsung_T3/Tof-Tracer-18-03/ExampleFiles/resultFiles/_result.hdf5"
calbrationResultFile = "/media/wiebke/Samsung_T3/Tof-Tracer-18-03/ExampleFiles/resultFiles/_resultCalibs.hdf5"
calibrationMasses = [MasslistFunctions.massFromComposition(C=6, H=12, O=1) APINENE[1]]
calibrationMassesConcentrationsPPT = [1000 1000]
calibrationMassIdForUnknowns = 1 # With the calibration of this mass, concentrations of all masses not in the calgas are estimated.

includeLossCorr = true 	# then check additionally Fan100times and Fan12times in "correction" part

########################### Select changing Masses #############################
BGs=[DateTime(2017,10,14,13,40) DateTime(2017,10,14,14,40)] 

measResultsForMassSelection = ResultFileFunctions.loadResults(measurementResultFile, startTime=DateTime(2017,10,14,12,00), endTime=DateTime(2017,10,19,18,40), useAveragesOnly=false)
measResultsForMassSelection.Times = measResultsForMassSelection.Times .- timedelay
bgSelection = (measResultsForMassSelection.Times .> BGs[1,1]) .& (measResultsForMassSelection.Times .< BGs[1,2])
signalSelection = (measResultsForMassSelection.Times .> DateTime(2017,10,14,18,30)) .& (measResultsForMassSelection.Times .< DateTime(2017,10,14,20,00))
sortedIndices, means, stderror = ResultFileFunctions.findChangingMasses(
    measResultsForMassSelection.MasslistMasses,
    measResultsForMassSelection.MasslistCompositions,
    measResultsForMassSelection.Traces,
    measResultsForMassSelection.Times,
    bgSelection,
    signalSelection,
    sigmaThreshold = 5,
    noNitrogen = true,
    onlySaneMasses = true,
    filterCrosstalkMasses=true
    )
selectedMasses = measResultsForMassSelection.MasslistMasses[sortedIndices]
s=  (selectedMasses .> 17) .& (selectedMasses .< 450)
selectedMasses = selectedMasses[s]

################################################################################



########################## Load Selected Traces and calibrate ##################
measResults = ResultFileFunctions.loadResults(measurementResultFile, massesToLoad=selectedMasses, useAveragesOnly=false)
measResults.Times = measResults.Times .- timedelay

########### BG Correction ##############
BG = CalibrationFunctions.generateBgTraces(measResults.Times, InterpolationFunctions.smooth(measResults.Traces, 0.01), slices=1, quant=0.00)
measResults.Traces = measResults.Traces - BG


############### Calibration ############
calibResults =  ResultFileFunctions.loadResults(calbrationResultFile, massesToLoad=calibrationMasses, useAveragesOnly=true)
calibResults.Times = calibResults.Times .- timedelay
if length(calibResults.MasslistMasses) != length(calibrationMasses)
    println("Could not find all given calib masses in calib result. Found: $(calibResults.MasslistMasses)")
end
calibCountsTraces = Array{typeof(measResults.Traces[1,1])}(size(measResults.Traces,1), length(calibResults.MasslistMasses))
for i=1:length(calibResults.MasslistMasses)
    calibCountsTraces[:,i] = CalibrationFunctions.generateCalibFactorTrace(measResults.Times, calibResults.Times, calibResults.Traces[:,i], cloud12transitions())
end
meanCalibFactors = mean(calibCountsTraces,1)
PyPlot.figure()
PyPlot.plot(measResults.Times, calibCountsTraces)
PyPlot.plot(calibResults.Times, calibResults.Traces, "o")
legStrings = []
for i = 1:size(calibCountsTraces,2)
  push!(legStrings,"m/z $(calibrationMasses[i]) / $(calibResults.MasslistMasses[i])")
end
PyPlot.legend(legStrings)

for i=1:length(measResults.MasslistMasses)
    found = false
    for j=1:length(calibrationMasses)
        if measResults.MasslistMasses[i] == calibrationMasses[j]
            println("Found calib for m/z $(measResults.MasslistMasses[i]), calibrating directly.")
            found = true
            measResults.Traces[:,i] = calibrationMassesConcentrationsPPT[j] .* measResults.Traces[:,i] ./ calibCountsTraces[:,j]
        end
    end
    if !found
        mf = sqrt(calibrationMasses[calibrationMassIdForUnknowns]/measResults.MasslistMasses[i]) .* calibrationMassesConcentrationsPPT[calibrationMassIdForUnknowns] ./ meanCalibFactors[calibrationMassIdForUnknowns]
        println("Estimated Sensitivity for m/z $(measResults.MasslistMasses[i]), mean: $(1/mf)cps/ppt")
        measResults.Traces[:,i] = sqrt(calibrationMasses[calibrationMassIdForUnknowns]/measResults.MasslistMasses[i]) .* calibrationMassesConcentrationsPPT[calibrationMassIdForUnknowns] .* measResults.Traces[:,i] ./ calibCountsTraces[:,calibrationMassIdForUnknowns]
    end
end

########################## Export it ###########################################
ExportFunctions.exportTracesCSV("$(dirname(measurementResultFile))/", measResults.MasslistElements, measResults.MasslistCompositions, Dates.datetime2unix.(measResults.Times), measResults.Traces; average=120)

########################## find inlet loss correction factors and export corrected traces ###########################################
Fan100times = (measResults.Times .> DateTime(2017,10,19,13,50)) & (measResults.Times .< DateTime(2017,10,19,14,35))
Fan12times = (measResults.Times .> DateTime(2017,10,19,14,50)) & (measResults.Times .< DateTime(2017,10,19,16,40))

if includeLossCorr == true 
    InletCorrection = 4/0.66    #scale CorrFactor linearly with LossFactor from 1 to 5 for substances with "O"<=4, substances with "O">4: corrected by factor 5
    Fan100Conc = mean(measResults.Traces[Fan100times,:],1)[1,:] 
    Fan100stdErr = std(measResults.Traces[Fan100times,:],1)[1,:] ./ sqrt(size(measResults.Traces[Fan100times,:],2))
    Fan12Conc =  mean(measResults.Traces[Fan12times,:],1)[1,:] 
    Fan12stdErr = std(measResults.Traces[Fan12times,:],1)[1,:] ./ sqrt(size(measResults.Traces[Fan12times,:],2))
    LossFactor = 1 - Fan100Conc./Fan12Conc 
    LossFactorErr = abs.(LossFactor) .* sqrt((Fan100stdErr./Fan100Conc).*(Fan100stdErr./Fan100Conc) + (Fan12stdErr./Fan12Conc).*(Fan12stdErr./Fan12Conc))
    CorrFactor = ones(LossFactor)
    CorrFactorErr = zeros(LossFactor)
    CorrectionNotes = similar(CorrFactor, String)
    for i=1:length(LossFactor)
        if measResults.MasslistCompositions[6,i] <= 4
	    if LossFactor[i] >= 3*LossFactorErr[i] 
	        CorrFactor[i] = 1+InletCorrection*LossFactor[i]
	        CorrFactorErr[i] = InletCorrection*LossFactorErr[i]
	        CorrectionNotes[i] = "corr. lin. to lossfactor"
	    elseif LossFactor[i] < 3*LossFactorErr[i] && measResults.MasslistCompositions[6,i] <= 4
	        if LossFactor[i] < 0
		    CorrectionNotes[i] = "not corr, increased during fanstage"
		    CorrFactorErr[i] = LossFactorErr[i]
	        else
	       	    CorrFactorErr[i] = 5
	    	    CorrectionNotes[i] = "not corr., lower limit conc. Use large uncertainty of factor 5"
	        end
	        CorrFactor[i] = 1
	    end
        elseif  measResults.MasslistCompositions[6,i] > 4
	    CorrFactor[i] = 5
	    CorrFactorErr[i] = LossFactorErr[i]*(5/LossFactor[i])
    	    CorrectionNotes[i] = "max. corr."
        end
    end
    # export additional file with corrections
    ExportFunctions.exportTracesCSVLossCorr("$(dirname(measurementResultFile))/", measResults.MasslistElements, measResults.MasslistCompositions, Dates.datetime2unix.(measResults.Times), measResults.Traces, round.(LossFactor,5), round.(LossFactorErr,5), round.(CorrFactor,5), round.(CorrFactorErr,5), CorrectionNotes; average=120)
end

########################## Plot Some Masses ####################################
PyPlot.figure()
PyPlot.semilogy(InterpolationFunctions.averageSamples(measResults.Times,20), InterpolationFunctions.averageSamples(measResults.Traces[:,1:70],20), "-")
legStrings = []
for i = 1:length(measResults.MasslistMasses)
  push!(legStrings,"m/z $(round(measResults.MasslistMasses[i],3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResults.MasslistCompositions[:,i]))")
end
PyPlot.legend(legStrings)
