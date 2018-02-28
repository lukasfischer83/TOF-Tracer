push!(LOAD_PATH, pwd())
using HDF5
using ExportFunctions
import  MasslistFunctions
using ResultFileFunctions

using PyPlot

#file = "/data/CLOUD11/running1/sortedByPrecursor/1805-1806-bcary/results/_result.hdf5"
#file = "/data/CLOUD11/running2/comparisonEESI/results/_result.hdf5"
#file = "/media/lukas/Auswertung2/CLOUD12/H3O-hard-2/results/_result.hdf5"
file = "/media/lukas/Auswertung2/CLOUD12/H3O-soft-2/results/_result.hdf5"
#file = "/media/lukas/Auswertung2/CLOUD12/H3O-hard-2//results/_result-LIM.hdf5"
#file = "/media/lukas/Auswertung2/CLOUD12/H3O-hard-2/results/_result-NewAlgorithm.hdf5"
#file = "/media/lukas/Auswertung2/CLOUD12/H3O-soft-2/run1940/results/_result.hdf5"
#file = "/media/lukas/Auswertung2/CLOUD12/H3O-hard-2/run1940/results/_result.hdf5"
#file = "/data/CLOUD11/running1/Martin-AP-IP-Pinonaldehyd/results/_result.hdf5"
#filepath = abspath("/media/lukas/Auswertung1/LTOF-CIMS")
#file = "/data/HYYDE2016Data/running1/_result.h5"
plotHighTimeRes = true
noNitrogen = true
sortSignal = "mean"
#sortSignal = "supersaturation"
sortResultsByMass = false
signifficanceInSigma = 5

if (plotHighTimeRes == false)
masses = HDF5.h5read(file, "MassList")
traces = HDF5.h5read(file, "CorrAvgStickCps")
timesUnix = HDF5.h5read(file, "AvgStickCpsTimes")
else
  masses = HDF5.h5read(file, "MassList")
  traces = HDF5.h5read(file, "CorrStickCps")
  timesUnix = HDF5.h5read(file, "Times")
end
compositions = HDF5.h5read(file, "ElementalCompositions")
elementNames = HDF5.h5read(file, "ElementNames")

times = Array(DateTime, length(timesUnix))

for i=1:length(times)
  if (isnan(timesUnix[i]))
    if i>2
      timesUnix[i] = timesUnix[i-1] + timesUnix[i-1] - timesUnix[i-2]
    else
      timesUnix[i] = 0
    end
  end
  #times[i] = Dates.unix2datetime(timesUnix[i]) - Dates.Hour(1) #CLOUD11
  times[i] = Dates.unix2datetime(timesUnix[i]) - Dates.Hour(0) #CLUOD12
end

#Martin CLOUD11 AP IP OH masses
#selBGTimes = (times.>DateTime(2016,10,06,13,00)) & (times.<DateTime(2016,10,06,14,00))
#selTimes = (times.>DateTime(2016,10,03,21,00)) & (times.<DateTime(2016,10,03,23,00))
#selFanTimes = (times.>DateTime(2016,10,04,00,15)) & (times.<DateTime(2016,10,04,00,40))
#plotTitle = "Alpha Pinene run1802"

#selTimes = (times.>DateTime(2016,10,24,01)) & (times.<DateTime(2016,10,25,05))
#Isoprene 17xx
#selBGTimes = (times.>DateTime(2015,11,12,1)) & (times.<DateTime(2015,11,12,4))
#selTimes = (times.>DateTime(2015,11,12,16)) & (times.<DateTime(2015,11,12,19))

#TMB 1780
#selBGTimes = (times.>DateTime(2015,12,10,4)) & (times.<DateTime(2015,12,10,6))
#selTimes = (times.>DateTime(2015,12,10,18)) & (times.<DateTime(2015,12,10,20))

#BCY
#selBGTimes = (times.>DateTime(2016,10,07,20,00)) & (times.<DateTime(2016,10,07,22,00))
#selTimes = (times.>DateTime(2016,10,8,1,0)) & (times.<DateTime(2016,10,8,2,0))

#NAPHTHA
#selBGTimes = (times.>DateTime(2016,11,06,03,30)) & (times.<DateTime(2016,11,06,05,00))
#selTimes = (times.>DateTime(2016,11,06,08,00)) & (times.<DateTime(2016,11,06,10,00))

#Superbanana
#selBGTimes = (times.>DateTime(2015,09,30,12,30)) & (times.<DateTime(2015,09,30,13,40))
#selTimes = (times.>DateTime(2015,09,30,17,30)) & (times.<DateTime(2015,09,30,19,30))

#CLOUD12 AP 1933, BG vs Plateau
#selBGTimes = (times.>DateTime(2017,10,14,12,35)) & (times.<DateTime(2017,10,14,15,00))
#selTimes = (times.>DateTime(2017,10,15,07,00)) & (times.<DateTime(2017,10,15,08,00))
#plotTitle = "Oxidation Products Alpha Pinene + O3"
#plotTitle = "Alpha Pinene run1933 E/N=hard"

#CLOUD12 AP 1933, Fan vs Plateau
#selBGTimes = (times.>DateTime(2017,10,15,10,00)) & (times.<DateTime(2017,10,15,11,00))
#selTimes = (times.>DateTime(2017,10,15,07,00)) & (times.<DateTime(2017,10,15,08,00))

#CLOUD12 LiM 1935, pre O3 vs post O3
#selBGTimes = (times.>DateTime(2017,10,16,23,30)) & (times.<DateTime(2017,10,17,1,00))
#selTimes = (times.>DateTime(2017,10,17,04,00)) & (times.<DateTime(2017,10,17,04,45))
#plotTitle = "Oxidation Products Limonene + O3"

#CLOUD12 Mix AP LIM Soft 1938
#selBGTimes = (times.>DateTime(2017,10,17,23,30)) & (times.<DateTime(2017,10,18,0,40))
#selTimes = (times.>DateTime(2017,10,18,06,00)) & (times.<DateTime(2017,10,18,07,30))
#plotTitle = "Oxidation Products Mix + O3"

#CLOUD12 AP Soft 1939
#selBGTimes = (times.>DateTime(2017,10,18,19,30)) & (times.<DateTime(2017,10,18,20,30))
#selTimes = (times.>DateTime(2017,10,19,05,00)) & (times.<DateTime(2017,10,19,06,30))
#selFanTimes = (times.>DateTime(2017,10,19,04,00)) & (times.<DateTime(2017,10,19,04,20))
#plotTitle = "Alpha Pinene run1939 E/N=soft"

#CLOUD12 AP Soft 1937
#selBGTimes = (times.>DateTime(2017,10,17,12,30)) & (times.<DateTime(2017,10,17,13,30))
#selTimes = (times.>DateTime(2017,10,17,16,30)) & (times.<DateTime(2017,10,17,17,30))
#selFanTimes = (times.>DateTime(2017,10,17,17,50)) & (times.<DateTime(2017,10,17,18,20))
#plotTitle = "Alpha Pinene run1937 UVH on E/N=soft"

#CLOUD12 AP Soft 1937
#selBGTimes = (times.>DateTime(2017,10,17,11,30)) & (times.<DateTime(2017,10,17,13,30))
#selTimes = (times.>DateTime(2017,10,17,16,30)) & (times.<DateTime(2017,10,17,17,30))
#selFanTimes = (times.>DateTime(2017,10,19,17,50)) & (times.<DateTime(2017,10,19,18,20))
#plotTitle = "Alpha Pinene run1937 UVH off E/N=soft"
#a = 1.3

#CLOUD12 AP IP Hard 1940
#selBGTimes = (times.>DateTime(2017,10,19,20,00)) & (times.<DateTime(2017,10,19,21,00))
#selTimes = (times.>DateTime(2017,10,20,15,30)) & (times.<DateTime(2017,10,20,16,20))
#selFanTimes = (times.>DateTime(2017,10,20,17,00)) & (times.<DateTime(2017,10,20,17,45))
#plotTitle = "Alpha Pinene run1940 UVH on E/N=hard"


#CLOUD12 AP IP Soft 1940
selBGTimes = (times.>DateTime(2017,10,20,22,00)) & (times.<DateTime(2017,10,20,22,30))
selTimes = (times.>DateTime(2017,10,20,19,30)) & (times.<DateTime(2017,10,20,20,30))
selFanTimes = (times.>DateTime(2017,10,20,21,00)) & (times.<DateTime(2017,10,20,21,45))
plotTitle = "Alpha Pinene run1940 UVH on E/N=soft"

sortedIndices, allMeans, allStderrors = ResultFileFunctions.findChangingMasses(masses, compositions, traces, times, selBGTimes, selTimes, signifficanceInSigma, sorting="mean", noNitrogen = true)
println("selMasses: $(size(sortedIndices))\nmeans: $(size(allMeans))\nstderror: $(size(allStderrors))")

sortedIndices = sortedIndices[(masses[sortedIndices].>65) & (masses[sortedIndices].<450)]

for i = 1:length(sortedIndices)
  println("$(masses[sortedIndices[i]])\t$(sumFormulaStringFromCompositionArray(compositions[:,sortedIndices[i]]))\t$(allMeans[sortedIndices[i]])")
end

#exportTracesCSV("/home/lukas/ownCloud/documents/UNI/Auswertung/CLOUD12/WallLosses", elementNames, compositions20sorted, timesUnix, tracesSelected20sorted; average=100)
#plot(compositions20sorted[6,:], meansBG20sorted./meansSignal20sorted, "o")
#xlabel("Number of Oxygens")
#ylabel("Loss Ratio High Fan / Low Fan")


using PlotFunctions
OtoC = compositions[6,sortedIndices]./compositions[1,sortedIndices] # O/C
PlotFunctions.massDefectPlot(masses[sortedIndices], compositions[:,sortedIndices], allMeans[sortedIndices], OtoC, plotTitle, "O/C", minConc = 0.01, dotSize=20)


allMeansFan = mean(traces[selFanTimes,:],1)[1,:] - mean(traces[selBGTimes,:],1)[1,:]

volatility = 1- allMeansFan[sortedIndices]./allMeans[sortedIndices] # volatility
volatility[volatility.<0] = 0
volatility[volatility.>1] = 1
#volatility = 1-volatility
#saturation = saturationFromComposition(compositions[:,sortedIndices])
PlotFunctions.massDefectPlot(masses[sortedIndices], compositions[:,sortedIndices], allMeans[sortedIndices], volatility, plotTitle, "Concentration Drop Rate during high Fan period", minConc = 0.01, dotSize=20, sumformulas = true)

#=
figure()
scatter(volatility,OtoC, masses[sortedIndices].*masses[sortedIndices]/maximum( masses[sortedIndices]))
xlabel("Volatility")
ylabel("O/C")


figure()
scatter(compositions[6,sortedIndices], volatility, 20*log(allMeans[sortedIndices]/0.2),compositions[1,sortedIndices])
cb=colorbar()
cb["ax"]["set_ylabel"]("Carbon")
ylabel("Wall Loss Ratio")
xlabel("MyParameter")
ylim(0.01,0.99)
=#
#exportTracesCSV("/home/lukas/ownCloud/documents/UNI/Auswertung/CLOUD11/Martin_run1802/", elementNames, compositions[:,sortedIndices], timesUnix, traces[:,sortedIndices]; average=10)
