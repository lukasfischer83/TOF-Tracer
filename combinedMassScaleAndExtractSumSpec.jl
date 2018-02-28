push!(LOAD_PATH, pwd())

import HDF5
import PyPlot
import TOFFunctions
import InterpolationFunctions

function correctMassScaleAndExtractSumSpec(
  filepath,
  masslistMasses,
  masslistElements,
  masslistElementsMasses,
  masslistCompositions,
  referenceFile,
  calibRegions; # Regions which are pattern matched for mass shifts
  filefilterRegexp=r"\.h5$",
  outputfilename="results/_result.hdf5",
  firstNFiles=0, # only analyze the first N files, set to 0 for all files
  lastNFiles=0, # only analyze the last N files, set to 0 for all files
  debuglevel=3,
  searchWidth = 0.7, # Width of Mass Scale Search Regions in AMU
  dynamicMassScaleCorrection = true,
  recalibInterval = 61,
  peakshapeRegions = 10,
  createTotalAvg = true, # Needed for manualPeakFitter, not needed for timetraces
  onlyUseAverages = false, # Fast mode using only total file averages instead of individual spectra in file
  storeSubSpectra = true,
  plotControlMass = false, # plot the control mass peak for each file after correction
  testRangeStart = 137.0, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
  testRangeEnd = 137.5,
  massBorderCalculation = 2, # How to calculate borders? 0 = Cernter -0.1 to Center + 0.4,  1 = based on resolution, 2 = constant bin width
  binWidth = 6,
  resolution = 7500)



  tic()
  if (isfile(joinpath(filepath, outputfilename)))
    mv(joinpath(filepath, outputfilename), joinpath(filepath, "$outputfilename.bak"), remove_destination=true)
    if debuglevel > 0   println("Found existing result file, moving to '$outputfilename.bak'") end
  end

  #files = filter(r"2016-07-04.*\.h5", readdir(filepath)) #one day only
  files = filter(filefilterRegexp, readdir(filepath))




  nFiles = size(files,1)

  if firstNFiles != 0
    if length(files) > firstNFiles
      files = files[1:firstNFiles]
    end
  end

  if lastNFiles != 0
    if length(files) > lastNFiles
      files = files[length(files)-lastNFiles:end]
    end
  end

  nFiles = size(files,1)

  validFiles = TOFFunctions.validateHDF5Files(filepath, files)

  files = validFiles
  nFiles = size(files,1)


  if (referenceFile == "")
    referenceFile = joinpath(filepath, files[1])
    if (debuglevel > 1)
      println("Reference File Autodetected: $referenceFile")
    end
  end


  referenceSpectrum = TOFFunctions.getAvgSpectrumFromFile(referenceFile)
  referenceMassScaleMode, referenceMassScaleParameters = TOFFunctions.getMassCalibParametersFromFile(referenceFile)
  TOFFunctions.setMassScaleReferenceSpectrum(referenceSpectrum, calibRegions, searchWidth, referenceMassScaleMode, referenceMassScaleParameters, plotControlMass=plotControlMass, testRangeStart=testRangeStart, testRangeEnd=testRangeEnd)
  referenceMassAxis = []
  referenceMassAxis = HDF5.h5read(referenceFile, "FullSpectra/MassAxis")

  totalAvgSpectrum = SharedArray{Float64}(length(referenceSpectrum))
  totalAvgSubSpectrum = SharedArray{Float64}(length(referenceSpectrum))
  totalMinSpectrum = SharedArray{Float64}(length(referenceSpectrum))
  totalMaxSpectrum = SharedArray{Float64}(length(referenceSpectrum))

  @sync @parallel for i = 1 : length(referenceSpectrum)
    totalAvgSpectrum[i] = 0
    totalAvgSubSpectrum[i] = 0
    totalMinSpectrum[i] = 1e99
    totalMaxSpectrum[i] = 0
  end



  nMasses=length(masslistMasses)
  println("Calculating Stick CPS for $nMasses masses.")
  timebinshifts = zeros(Float64,(length(calibRegions),nFiles))
  intensities = zeros(Float64,(length(calibRegions),nFiles))
  monitorTimetrace = zeros(nFiles)

  time = SharedArray{DateTime}(nFiles)
  stickcps=SharedArray{Float64}(nMasses,nFiles)
  # Calculate Integration Borders
  if (massBorderCalculation == 0)
  mlow = masslistMasses - 0.1
  mhigh = masslistMasses + 0.4
  elseif  (massBorderCalculation == 1)
  mlow = masslistMasses .* resolution./(resolution+0.5)
  mhigh = masslistMasses .* resolution./(resolution-0.5)
  end


  if (plotControlMass == true)
    figure()
  end

  if (createTotalAvg == true)
    interpolatedSpectrum = SharedArray{Float64}(length(referenceSpectrum))
  end

  ############## Check output path and remove existing files #####################
  if (!isdir("$filepath/results"))
    mkdir("$filepath/results")
  end
  outfilepath = joinpath(filepath, outputfilename)
  if isfile(outfilepath)
    rm(outfilepath)
  end

  fid = HDF5.h5open(outfilepath, "w")

  ################################################################################

  ############## open outputfile and create extendable SumSpecs Array ############
  if (createTotalAvg == true)

    dsAvgSumWidth = length(referenceSpectrum)
    dspaceAvgSumSpecs = HDF5.dataspace((dsAvgSumWidth,1)::Dims, max_dims=(dsAvgSumWidth,typemax(Int64)))
    dtypeAvgSumSpecs = HDF5.datatype(Float32)
    dsetAvgSumSpecs = HDF5.d_create(fid, "SumSpecs", dtypeAvgSumSpecs, dspaceAvgSumSpecs, "chunk", (dsAvgSumWidth,1))
  end

  if (!onlyUseAverages)
    specCount = 0
    dsStickCpsWidth = nMasses
    dspaceStickCps = HDF5.dataspace((1,dsStickCpsWidth)::Dims, max_dims=(typemax(Int64),dsStickCpsWidth))
    dtypeStickCps = HDF5.datatype(Float32)
    dsetStickCps = HDF5.d_create(fid, "StickCps", dtypeStickCps, dspaceStickCps, "chunk", (1,dsStickCpsWidth))
    dsetStickCpsErr = HDF5.d_create(fid, "StickCpsErr", dtypeStickCps, dspaceStickCps, "chunk", (1,dsStickCpsWidth))

    dspaceTimes = HDF5.dataspace((1,)::Dims, max_dims=(typemax(Int64),))
    dtypeTimes = HDF5.datatype(Float64)
    dsetTimes = HDF5.d_create(fid, "Times", dtypeTimes, dspaceTimes, "chunk", (1,))

  end
  ################################################################################
  badFiles = Array{String}(0)

  for j=1:nFiles
    totalPath = joinpath(filepath, files[j])
    if debuglevel > 0   println("Processing File $j/$nFiles :  $totalPath") end

    fileIsBad = false;
    massAxis = []
    massAxis = HDF5.h5read(totalPath, "FullSpectra/MassAxis")
    time[j] = TOFFunctions.getTimeFromFile(totalPath)

    avgSpectrum = TOFFunctions.getAvgSpectrumFromFile(totalPath)
    newParams, success, tbs, ins = TOFFunctions.recalibrateMassScale(avgSpectrum, referenceSpectrum, calibRegions, searchWidth, referenceMassScaleMode, referenceMassScaleParameters)

    #if debuglevel > 0   println("Processing File $totalPath") end
    if (createTotalAvg == true && !fileIsBad)
      #println("153: rma: $(length(referenceMassAxis)), rmsm: $referenceMassScaleMode, p: $newParams")
      indexesExact = TOFFunctions.mass2timebin(referenceMassAxis, referenceMassScaleMode, newParams)
      interpolatedSpectrum = InterpolationFunctions.interpolate(indexesExact,avgSpectrum)
      totalAvgSpectrum += interpolatedSpectrum
      # calculate min and max spectra
      @sync @parallel for bin=1:length(interpolatedSpectrum)
        if interpolatedSpectrum[bin] > totalMaxSpectrum[bin]
          totalMaxSpectrum[bin] = interpolatedSpectrum[bin]
        end
        if interpolatedSpectrum[bin] < totalMinSpectrum[bin]
          totalMinSpectrum[bin] = interpolatedSpectrum[bin]
        end
      end
      # Write to hdf5, line by line, so there is no limit to number of files that can fit in RAM
      if (storeSubSpectra == true)
        currDims = size(dsetAvgSumSpecs)[2]
        dsetAvgSumSpecs[:,currDims] = interpolatedSpectrum
        if (j != nFiles)
          HDF5.set_dims!(dsetAvgSumSpecs, (dsAvgSumWidth,currDims+1)::Dims)
        end
      end
    end
    if debuglevel > 1   println() end
    for i=(1:nMasses)
      #if debuglevel > 0   println("Processing region mass($(mlow[i]):$(mhigh[i])) --> timebin($(round(TOFFunctions.mass2timebin(mlow[i],massCalibMode,newParams))):$(round(TOFFunctions.mass2timebin(mhigh[i],massCalibMode,newParams))))") end
      if (massBorderCalculation == 2)
        centerIndex = Int64(floor(TOFFunctions.mass2timebin(masslistMasses[i], referenceMassScaleMode, newParams)))
        if ((centerIndex-binWidth) > 0) && ((centerIndex+binWidth)<length(avgSpectrum))
          raw = sum(view(avgSpectrum, (centerIndex-binWidth : centerIndex + binWidth)))
        else
          raw = 0
        end
      else
        subIdxStartExact=TOFFunctions.mass2timebin(mlow[i],referenceMassScaleMode,newParams)
        subIdxEndExact = TOFFunctions.mass2timebin(mhigh[i],referenceMassScaleMode,newParams)
        print("summing $subIdxStartExact to $subIdxEndExact: ")
        raw = interpolatedSum(subIdxStartExact,subIdxEndExact,avgSpectrum)
        println(raw)
      end
      stickcps[i,j] = raw
    end
    if (!onlyUseAverages)
      spectrumMultFactor = TOFFunctions.getSpecMultiplicator(totalPath)
      subSpecStickCps = SharedArray{Float32}(nMasses)
      totalSubSpectra = TOFFunctions.getSubSpectraCount(totalPath)

      ################## prepare arrays for averaging #########################
      fileInternalLocalAvg = TOFFunctions.getSubSpectrumFromFile(totalPath,1)
      #println("referenceSpectrum is a $(summary(referenceSpectrum))")
      #println("SubSpec is a $(summary(fileInternalLocalAvg))")
      #println("spectrumMultFactor is a $(summary(spectrumMultFactor))")

      fill!(fileInternalLocalAvg,0)
      fileInternalLocalAvgCount = 0
      ################## Get First Set of Spectra for first mass scale calib
      if totalSubSpectra >= recalibInterval
        # Prepare first calib beforehead
        fileInternalLocalAvg = TOFFunctions.getSubSpectrumFromFile(totalPath,1)
        for avgIdx=2:minimum([recalibInterval totalSubSpectra])
          fileInternalLocalAvg += TOFFunctions.getSubSpectrumFromFile(totalPath,avgIdx)
        end
        newParams, success, tbs, ins = TOFFunctions.recalibrateMassScale(fileInternalLocalAvg*spectrumMultFactor, referenceSpectrum, calibRegions, searchWidth, referenceMassScaleMode, referenceMassScaleParameters)
        ######
        fill!(fileInternalLocalAvg,0)
        fileInternalLocalAvgCount = 0
      end
      #######################################################################

      ################### Loop over Subspectra Start ########################
      ################### Averaging and Mass Scale Calib ####################

      specCount += totalSubSpectra
      for subSpecIdx=1:totalSubSpectra
        subSpectrum = TOFFunctions.getSubSpectrumFromFile(totalPath,subSpecIdx)
        fileInternalLocalAvg += subSpectrum
        fileInternalLocalAvgCount += 1

        if (fileInternalLocalAvgCount > recalibInterval) || fileInternalLocalAvgCount == totalSubSpectra
          println("      Recalibrating at $subSpecIdx")
          newParams, success, tbs, ins = TOFFunctions.recalibrateMassScale(fileInternalLocalAvg, referenceSpectrum, calibRegions, searchWidth, referenceMassScaleMode, referenceMassScaleParameters)

          indexesExact = TOFFunctions.mass2timebin(referenceMassAxis, referenceMassScaleMode, newParams)
          interpolatedSpectrum = InterpolationFunctions.interpolate(indexesExact,fileInternalLocalAvg)

          totalAvgSubSpectrum += interpolatedSpectrum*spectrumMultFactor

          fill!(fileInternalLocalAvg,0)
          fileInternalLocalAvgCount = 0
        end
        ################## Peak Integration ################################
        @sync @parallel for i=(1:nMasses)
          if (mod(i,100) == 0)
            #println("Processing mass $(masslistMasses[i])")
          end
          if (massBorderCalculation == 2)
            centerIndex = Int64(floor(TOFFunctions.mass2timebin(masslistMasses[i], referenceMassScaleMode, newParams)))
            if ((centerIndex-binWidth) > 0) && ((centerIndex+binWidth)<length(avgSpectrum))
              raw = sum(view(subSpectrum, (centerIndex-binWidth : centerIndex + binWidth)))
            else
              raw = 0
            end
            subSpecStickCps[i]=raw
          else
            #if debuglevel > 0   println("Processing region mass($(mlow[i]):$(mhigh[i])) --> timebin($(round(TOFFunctions.mass2timebin(mlow[i],massCalibMode,newParams))):$(round(TOFFunctions.mass2timebin(mhigh[i],massCalibMode,newParams))))") end
            subIdxStartExact=TOFFunctions.mass2timebin(mlow[i],referenceMassScaleMode,newParams)
            subIdxStart::Int64 = ceil(subIdxStartExact)
            subIdxStartRoundError = subIdxStart - subIdxStartExact

            subIdxEndExact = TOFFunctions.mass2timebin(mhigh[i],referenceMassScaleMode,newParams)
            subIdxEnd::Int64 = floor(subIdxEndExact)
            subIdxEndRoundError = subIdxEndExact - subIdxEnd
            raw = sum(view(subSpectrum, subIdxStart:subIdxEnd)) + subSpectrum[subIdxStart-1]*subIdxStartRoundError + subSpectrum[subIdxEnd+1]*subIdxEndRoundError
            subSpecStickCps[i]=raw
          end
        end
        rawTime = TOFFunctions.getSubSpectrumTimeFromFile(totalPath,subSpecIdx)
        currDimsTime = size(dsetTimes)[1]
        absTimeOfCurrSample = Dates.datetime2unix(time[j]) + rawTime
        dsetTimes[currDimsTime] = absTimeOfCurrSample
        dsetStickCps[currDimsTime,:] = subSpecStickCps*spectrumMultFactor
        if (currDimsTime == 1)
          deltaT = 5 #rawTime # Not correct???
        else
          deltaT = absTimeOfCurrSample - dsetTimes[currDimsTime-1]
        end
        #println("DeltaT: $deltaT")
        dsetStickCpsErr[currDimsTime,:] = sqrt.(abs.(subSpecStickCps./deltaT))
        if !((j==nFiles) && (subSpecIdx== TOFFunctions.getSubSpectraCount(totalPath)))
          HDF5.set_dims!(dsetTimes, (currDimsTime+1,)::Dims)
          HDF5.set_dims!(dsetStickCps, (currDimsTime+1,dsStickCpsWidth)::Dims)
          HDF5.set_dims!(dsetStickCpsErr, (currDimsTime+1,dsStickCpsWidth)::Dims)
        end
      end
    end
    println("###################################################################\n\n")
  end
  close(fid)

  stickcpsFlat = zeros(nMasses,nFiles)
  stickcpsErrFlat = zeros(nMasses,nFiles)

  timeFlat = zeros(nFiles)
  timeFlat = Dates.datetime2unix.(time)
  # We need more than one timestamp to calculate a difference
  if length(timeFlat) > 1
      deltaT = median(diff(timeFlat))
  else
      deltaT = 1
  end

  totalAvgSpectrumFlat = zeros(length(totalAvgSpectrum))
  totalAvgSubSpectrumFlat = zeros(length(totalAvgSubSpectrum))


  for i=1:nMasses, j=1:nFiles
    stickcpsFlat[i,j] = stickcps[i,j]
    stickcpsErrFlat[i,j] = sqrt(abs(stickcps[i,j])/deltaT)
  end

  for i=1:length(totalAvgSpectrum)
    totalAvgSpectrumFlat[i] = totalAvgSpectrum[i]/nFiles
  end
  if !onlyUseAverages
    for i=1:length(totalAvgSpectrum)
      totalAvgSubSpectrumFlat[i] = totalAvgSubSpectrum[i]/specCount
    end
  end

  HDF5.h5write(outfilepath, "AvgStickCps", stickcpsFlat')
  HDF5.h5write(outfilepath, "AvgStickCpsErr", stickcpsErrFlat')
  HDF5.h5write(outfilepath, "MassList", masslistMasses)
  println("Wrote $(length(masslistMasses)) masses to file.")
  HDF5.h5write(outfilepath, "AvgStickCpsTimes", timeFlat)
  HDF5.h5write(outfilepath, "MassAxis", referenceMassAxis)
  if onlyUseAverages
    HDF5.h5write(outfilepath, "AvgSpectrum", totalAvgSpectrumFlat)
  else
    HDF5.h5write(outfilepath, "AvgSpectrum", totalAvgSubSpectrumFlat)
  end

  if (createTotalAvg == true)
    #HDF5.h5write(outfilepath, "SumSpecs", sdata(interpolatedSpectra))
    HDF5.h5write(outfilepath, "SumSpecMax", sdata(totalMaxSpectrum))
    HDF5.h5write(outfilepath, "SumSpecMin", sdata(totalMinSpectrum))
  end
  HDF5.h5write(outfilepath, "ElementNames", masslistElements)
  HDF5.h5write(outfilepath, "ElementMasses", masslistElementsMasses)
  compFlat = Array{Int64}(length(masslistCompositions[1]),length(masslistCompositions))
  [compFlat[:,i]=masslistCompositions[i] for i=1:length(masslistCompositions)]
  HDF5.h5write(outfilepath, "ElementalCompositions", compFlat)
  print("Processing took: ")
  toc()
end
