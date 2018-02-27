__precompile__()
module ResultFileFunctions
using HDF5
using MasslistFunctions

export MeasurementResult, joinResultsTime, joinResultsMasses, getTraces, getTimetraces, getNbrTraces, getTraceSamples, getNbrTraceSamples, findChangingMasses, saturationFromComposition, getIndicesInTimeframe

type MeasurementResult
    Times
    MasslistMasses
    MasslistElements
    MasslistElementsMasses
    MasslistCompositions
    Traces
end

function joinResultsTime(firstResult::MeasurementResult, secondResult::MeasurementResult)
    if (length(firstResult.Times) > 0) & (length(secondResult.Times) > 0)
        if (firstResult.MasslistMasses == secondResult.MasslistMasses)
            firstResult.Times = vcat(firstResult.Times, secondResult.Times)
            firstResult.MasslistCompositions = hcat(firstResult.MasslistCompositions, secondResult.MasslistCompositions)
            firstResult.Traces = hcat(firstResult.Traces, secondResult.Traces)
        else
        println("Masses did not match, could not merge results!")
        end
    end
return firstResult
end

function joinResultsMasses(firstResult::MeasurementResult, secondResult::MeasurementResult)
    if (length(firstResult.MasslistMasses) > 0) & (length(secondResult.MasslistMasses) > 0)
        if (firstResult.Times == secondResult.Times)
            firstResult.MasslistMasses = vcat(firstResult.MasslistMasses, secondResult.MasslistMasses)
            firstResult.MasslistCompositions = hcat(firstResult.MasslistCompositions, secondResult.MasslistCompositions)
            firstResult.Traces = hcat(firstResult.Traces, secondResult.Traces)
        else
        println("Times did not match, could not merge results!")
        end
    end
return firstResult
end

function loadResults(filename; useAveragesOnly = false, raw = false, startTime::DateTime=DateTime(0), endTime::DateTime=DateTime(3000), massesToLoad=Array{Float64}(0), massMatchTolerance = 0.00001)
    println("Loading Times:")
    tic()
  if useAveragesOnly
    timesUnix = h5read(filename, "AvgStickCpsTimes")
  else
      try
          timesUnix = h5read(filename, "Times")
      catch
          println("Could not load high res time, maybe this is an average-only result file?")
          return []
      end
  end
  toc()
println("Looking for errors in time axis")
tic()
  ############## Check for Errors in time axis ###################
#  if ! issorted(timesUnix)
#    println("Time in result file is not in accending order, selection might not be complete!")
#  end
  if (isnan(timesUnix[1])) | (timesUnix[1] < Dates.datetime2unix(DateTime(1975,1,1))) | (timesUnix[1] > Dates.datetime2unix(DateTime(2030,1,1)))
      timesUnix[1] = timesUnix[1] - (timesUnix[3]-timesUnix[2])
  end
  if (isnan(timesUnix[length(timesUnix)]))| (timesUnix[length(timesUnix)] < Dates.datetime2unix(DateTime(1975,1,1))) | (timesUnix[length(timesUnix)] > Dates.datetime2unix(DateTime(2030,1,1)))
      timesUnix[length(timesUnix)] = timesUnix[length(timesUnix)] + (timesUnix[length(timesUnix)-1]-timesUnix[length(timesUnix)-2])
  end

  for i=2:length(timesUnix)-1
      if (isnan(timesUnix[i])) | (timesUnix[i] < Dates.datetime2unix(DateTime(1975,1,1))) | (timesUnix[i] > Dates.datetime2unix(DateTime(2030,1,1)))
          timesUnix[i] = 0.5*(timesUnix[i-1]+timesUnix[i+1])
      end
  end
  toc()
  println("Selecting Times and Masses:")
  tic()

  masslistCompositions = h5read(filename, "ElementalCompositions")
  masslistMasses = h5read(filename, "MassList")
  masslistElements = h5read(filename, "ElementNames")
  masslistElementsMasses = h5read(filename, "ElementMasses")

  if length(massesToLoad) > 0
      selectionMassesIndices = Array{Int}(0)
      for i=1:length(masslistMasses)
        for j=1:length(massesToLoad)
          if (isapprox(masslistMasses[i], massesToLoad[j], atol=massMatchTolerance))
            push!(selectionMassesIndices, i)
          end
        end
      end
  else
      selectionMassesIndices=1:length(masslistMasses)
  end

  selectionTimeIndexStart = searchsortedfirst(timesUnix, Dates.datetime2unix(startTime))
  println("Start index: $selectionTimeIndexStart")
  selectionTimeIndexEnd = searchsortedlast(timesUnix, Dates.datetime2unix(endTime))
  println("End index: $selectionTimeIndexEnd")

  if (selectionTimeIndexStart > 0) & (selectionTimeIndexEnd > selectionTimeIndexStart)
    selectedTimesUnix = timesUnix[selectionTimeIndexStart:selectionTimeIndexEnd]
  else
    selectedTimesUnix = timesUnix
    selectionTimeIndexStart = 1
    selectionTimeIndexEnd = length(timesUnix)
  end
  toc()
  ############# Load Traces here ##################

  selectedTimes = Dates.unix2datetime.(selectedTimesUnix)

  selectedMasslistMasses = masslistMasses[selectionMassesIndices]
  selectedMassesCompositions = masslistCompositions[:,selectionMassesIndices]

  if length(selectionMassesIndices) == 0
    return [], selectedMasslistMasses, masslistElements, masslistElementsMasses, Array{Float64}(0,0), Array{Float64}(0,0)
  end

  println("Loading Traces:")
  tic()
  traces = getTraces(filename, timeIndexStart=selectionTimeIndexStart, timeIndexEnd=selectionTimeIndexEnd, massIndices=selectionMassesIndices, raw=raw, useAveragesOnly=useAveragesOnly)
  toc()
  println("Loaded $(size(traces)) traces")
  return MeasurementResult(selectedTimes, selectedMasslistMasses, masslistElements, masslistElementsMasses, selectedMassesCompositions, traces)
end


function getIndicesInTimeframe(filename, startTime::DateTime, endTime::DateTime)
  times = h5read(filename, "Times")
  return (1:length(times))[(times.>startTime) & (times.<endTime)]
end

function getTraces(filename; timeIndexStart=1, timeIndexEnd=0, massIndices=nothing, raw=false,  useAveragesOnly = false)
  fh = h5open(filename,"r")
  if useAveragesOnly
    if raw
      ds = fh["AvgStickCps"]
    else
      ds = fh["CorrAvgStickCps"]
    end
  else
    if raw
      ds = fh["StickCps"]
    else
      ds = fh["CorrStickCps"]
    end
  end
  result = 0

  if timeIndexEnd == 0
    timeIndexEnd = size(ds)[1]
    println("No end time given, using all $(size(ds)[1]) elements")
  end

  if massIndices == nothing
    massIndices = ( 1:size(ds)[2])
    println("No mass selection given, selecting all $(size(ds)[2]) masses.")
  end

  if typeof(massIndices) == BitArray{1}
    massIndices = ( 1:length(massIndices) )[massIndices]
  end

  if ((timeIndexStart>0) & (timeIndexEnd<=size(ds)[1])) & (minimum(massIndices)>0) & (maximum(massIndices)<=size(ds)[2])
    result = similar(ds[1,1],(timeIndexEnd-timeIndexStart+1),length(massIndices))
    println("Creating result Matrix with size $(size(result))")
    if length(massIndices) < 100
        for j=1:length(massIndices)
            result[:,j] = ds[timeIndexStart:timeIndexEnd,massIndices[j]]
            #println("Loaded $((timeIndexEnd-timeIndexStart)*j) samples")
        end
    else
        for j=timeIndexStart:timeIndexEnd
            tmp = ds[j,:]
            result[j-timeIndexStart+1,:] = tmp[massIndices]
            #println("Loaded $(length(massIndices)*j) samples")
        end
    end
  end

  close(fh)
  return result
end


function getTimetraces(filename, indices; raw=false)
  fh = h5open(filename,"r")
  if raw
    ds = fh["StickCps"]
  else
    ds = fh["CorrStickCps"]
  end

  result = 0

  if typeof(indices) == BitArray{1}
    indices = ( 1:length(indices) )[indices]
  end

  if (typeof(indices) == Array{Int,1}) | (typeof(indices) == Array{Integer})
    if (minimum(indices)>0) & (maximum(indices)<size(ds)[2])
      result = similar(ds[1,1],size(ds)[1],length(indices))
      for i=1:length(indices)
        result[:,i] = ds[:,i]
      end
    end
  end

  if typeof(indices) == UnitRange
    if (minimum(indices)>0) & (maximum(indices)<size(ds)[2])
      result = ds[:,indices]
    end
  end

  close(fh)
  return result
end

function getTraceSamples(filename, indices; raw=false)
  fh = h5open(filename,"r")
  if raw
    ds = fh["StickCps"]
  else
    ds = fh["CorrStickCps"]
  end
  result = 0

  if typeof(indices) == BitArray{1}
    indices = ( 1:length(indices) )[indices]
  end

  if (typeof(indices) == Array{Int,1}) | (typeof(indices) == Array{Integer})
    if (minimum(indices)>0) & (maximum(indices)<=size(ds)[1])
      result = similar(ds[1,1],length(indices),size(ds)[2])
      for i=1:length(indices)
        result[i,:] = ds[i,:]
      end
    end
  end

  if typeof(indices) == UnitRange{Int}
    if (minimum(indices)>0) & (maximum(indices)<=size(ds)[1])
      result = ds[indices,:]
    end
  end

  close(fh)
  if result == 0
    println("Failed getting Traces Subset!")
  end
  return result
end

function getNbrTraces(filename)
  fh = h5open(filename,"r")
  ds = fh["StickCps"]
  nbr = size(ds)[2]
  close(fh)
  return nbr
end

function getNbrTraceSamples(filename)
  fh = h5open(filename,"r")
  ds = fh["StickCps"]
  nbr = size(ds)[1]
  close(fh)
  return nbr
end

function findChangingMasses(masses, compositions, traces, times, bgTimesSelection, signalTimesSelection; minOxygen = 0, sigmaThreshold=3, sorting = "mean", noNitrogen = true, onlySaneMasses = true, filterCrosstalkMasses=true)
  bgTraces = traces[bgTimesSelection,:]
  sigTraces = traces[signalTimesSelection,:]
  meansBG = mean(bgTraces,1)[1,:]
  meansBG[meansBG.<=0] = 1e-99
  means = mean(sigTraces,1)[1,:] - meansBG
  println("Means received: $(size(means))")
  if size(bgTraces,1) < 3
    println("Not enough BG points for standard deviation! plotHighTimeRes = true ?")
  end

  stderror = std(bgTraces,1)[1,:] ./ sqrt(size(bgTraces,2))

  selMasses =  (means .> sigmaThreshold*stderror)
  if filterCrosstalkMasses
      # Filter out crosstalk masses
      s = filterMassListByContribution2(masses, means, 5000, 0.05)
      selMasses = selMasses & s
  end
  println("\nRemoving $(length(masses[!s])) masses")

  if noNitrogen == true
    selMasses = selMasses & (compositions[5,:] .== 0)
  end

  selMasses = selMasses & (compositions[6,:] .>= minOxygen)
  #selMasses = selMasses| ((masses.>137.1) & (masses.<137.15))

  ## Filter only sane masses
  if onlySaneMasses
    selMasses = selMasses & (compositions[3,:] .> 1*compositions[1,:]) # H:C > 1.3
    selMasses = selMasses & (compositions[3,:] .< 2.2*compositions[1,:]) # H:C < 2.2
    selMasses = selMasses & (compositions[6,:] .< 1.5*compositions[1,:]) # O:C < 1.5
  end
  println("Selected $(sum(ones(length(masses))[selMasses])) masses!")

  selIndices = linearindices(masses)[selMasses]
  if sorting == "mean"
      println("Sorting masses by mean value!")
      sorter = sortperm(means[selMasses], rev=true)
  elseif sorting == "mass"
      println("Sorting masses!")
    sorter = sortperm(masses[selMasses], rev=false)
  else
      println("Not sorting.")
    sorter = linearindices(masses)
  end

  return selIndices[sorter], means[selIndices[sorter]], stderror[selIndices[sorter]]
end


function saturationFromComposition(compositions)
  #return exp(10)*50000*exp.(-compositions[6,:]*5).*exp.(-compositions[1,:])
  n_C = float(compositions[1,:])
  n_O = float(compositions[6,:])
  b_add = n_C
  b_add[n_C .<= 15] = 0.92  # monomers
  b_add[n_C .> 15] = 1.15 # dimers

  return 10.^((25-n_C)*0.475-n_O.*(2.3-b_add)-2.*(n_C.*n_O./(n_C+n_O))*(-0.3))
end
end
