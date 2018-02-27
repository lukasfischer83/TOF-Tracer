__precompile__()
module TOFFunctions
using HDF5
using InterpolationFunctions
using PyPlot

SPEC_CACHE_SIZE_LIMIT = 1e9

#using PyCall
#@pyimport h5py

export mass2timebin, timebin2mass, getMassCalibParametersFromFile, getSubSpectraCount, getSubSpectrumFromFile, getSpecMultiplicator, getSubSpectrumTimeFromFile, getAvgSpectrumFromFile, getTimeFromFile, validateHDF5Files, setMassScaleReferenceSpectrum, recalibrateMassScale

debuglevel = 0

type h5cache
    filename
    content
end

timeCache = h5cache("",0)
spectraCache = h5cache("",0)

# HELPER @everywhere functionS #
function mass2timebin(mass::Number,mode,parameters)
  if mode == 0
    return parameters[1]*sqrt(mass) + parameters[2]
  end
  if mode == 1
    return parameters[1]/sqrt(mass) + parameters[2]
  end
  if mode == 2
    return parameters[1]*mass^parameters[3] + parameters[2]
  end
end
function mass2timebin(mass::AbstractArray,mode,parameters)
  ret = Array{Float64}(length(mass))
  for i = 1:length(mass)
    ret[i] = mass2timebin(mass[i],mode,parameters)
  end
  return ret
end

function timebin2mass(time::Number,mode,parameters)
  if mode == 0
    return ((time-parameters[2])/parameters[1])^2
  end
  if mode == 1
    return (parameters[1]/(time-parameters[2]))^2
  end
  if mode == 2
    return ((time-parameters[2])/parameters[1])^(1/parameters[3])
  end
end

function timebin2mass(time::AbstractArray,mode,parameters)
  ret = Array{Float64}(length(time))
  for i = 1:length(time)
    ret[i] = timebin2mass(time[i],mode,parameters)
  end
  return ret
end

function getMassCalibParametersFromFile(filename)
    attributesFullSpectra = h5readattr(filename, "/FullSpectra")

    mcm = attributesFullSpectra["MassCalibMode"]
    massCalibMode = mcm[1]
    massCalibParameters = []

    if (massCalibMode == 0)
      p1 = attributesFullSpectra["MassCalibration p1"]
      p2 = attributesFullSpectra["MassCalibration p2"]
      massCalibParameters = [p1[1] p2[1]]
    end

    if (massCalibMode == 2)
      p1 = attributesFullSpectra["MassCalibration p1"]
      p2 = attributesFullSpectra["MassCalibration p2"]
      p3 = attributesFullSpectra["MassCalibration p3"]
      massCalibParameters = [p1[1] p2[1] p3[1]]
    end

    return massCalibMode, massCalibParameters
end
function getAvgSpectrumFromFile(filename)
  #println("Getting avg spectrum from $filename")

  attributesRoot = h5readattr(filename, "/")
  attributesFullSpectra = h5readattr(filename, "/FullSpectra")

  H5NbrWrites = attributesRoot["NbrWrites"]
  H5NbrBufs = attributesRoot["NbrBufs"]
  H5NbrWaveForms = attributesRoot["NbrWaveforms"]
  H5TofPeriod = h5readattr(filename, "/TimingData")["TofPeriod"]
  H5NbrSegments = attributesRoot["NbrSegments"]
  H5NbrBlocks = attributesRoot["NbrBlocks"]

  H5SampleInterval = attributesFullSpectra["SampleInterval"] .* 1e9
  H5SingleIonSignal = attributesFullSpectra["Single Ion Signal"]

  H5inttime = H5NbrWaveForms.*H5TofPeriod.*H5NbrSegments.*H5NbrBlocks.*H5NbrBufs.*H5NbrWrites .* 1e-9

  avgSpectrum = h5read(filename, "FullSpectra/SumSpectrum").*H5SampleInterval./H5SingleIonSignal./H5inttime
  return avgSpectrum
end

function getSubSpectraCount(filename)
  #println("Getting sub spectrum count from $filename")

  fh = h5open(filename,"r")
  ds = fh["/FullSpectra/TofData"]
  l=size(ds)[3]
  m=size(ds)[4]
  close(fh)
  return l*m
end

function getSubSpectrumFromFile(filename, index)
  #allSpectra = h5read(totalPath, "/FullSpectra/TofData")
  #allSpectraTimes = h5read(totalPath, "/TimingData/BufTimes")
  #k=1:size(allSpectra)[2], m=1:size(allSpectra)[4], l=1:size(allSpectra)[3]
  #allSpectraTimes = h5read(filename, "/TimingData/BufTimes")
  if spectraCache.filename != filename
      fh = h5open(filename,"r")
      ds = fh["/FullSpectra/TofData"]
      if prod(size(ds)) < SPEC_CACHE_SIZE_LIMIT
          println("Caching spectra of $filename")
          spectraCache.content = h5read(filename, "/FullSpectra/TofData") #ds[:,:,:,:]
          spectraCache.filename = filename
          close(fh)
      else
          (l,m) = ind2sub((size(ds)[3],size(ds)[4]), index)
          subSpectrum = ds[:,1,l,m]
          close(fh)
          return subSpectrum
      end
  end
  (l,m) = ind2sub((size(spectraCache.content)[3],size(spectraCache.content)[4]), index)
  subSpectrum = spectraCache.content[:,1,l,m]
  return subSpectrum
end

function getSpecMultiplicator(filename)
  attributesRoot = h5readattr(filename, "/")
  attributesFullSpectra = h5readattr(filename, "/FullSpectra")

  H5NbrWrites = attributesRoot["NbrWrites"]
  H5NbrBufs = attributesRoot["NbrBufs"]
  H5NbrWaveForms = attributesRoot["NbrWaveforms"]
  H5TofPeriod = h5readattr(filename, "/TimingData")["TofPeriod"]
  H5NbrSegments = attributesRoot["NbrSegments"]
  H5NbrBlocks = attributesRoot["NbrBlocks"]

  H5SampleInterval = attributesFullSpectra["SampleInterval"] .* 1e9
  H5SingleIonSignal = attributesFullSpectra["Single Ion Signal"]

  H5inttime = H5NbrWaveForms.*H5TofPeriod.*H5NbrSegments.*H5NbrBlocks .* 1e-9

  return (H5SampleInterval./H5SingleIonSignal./H5inttime)[1]
end

function getSubSpectrumTimeFromFile(filename, index)
    if filename != timeCache.filename
        fh = h5open(filename,"r")
        allSpectraTimes = h5read(filename, "/TimingData/BufTimes")
        ds = fh["/TimingData/BufTimes"]
        timeCache.filename = filename
        timeCache.content = ds[:,:]
        close(fh)
    end
    (l,m) = ind2sub((size(timeCache.content)[1],size(timeCache.content)[2]), index)
    #println("Getting sub spectrum Time ($l,$m) from $filename")
    time = timeCache.content[l,m]
    return time
end


function getTimeFromFile(filename)

  # TODO: Replace with Aquisition log Start Time (Windows Timestamp)
  #=
  f = h5py.File(filename, "r")
  aqLog = f["AcquisitionLog"]["Log"]
  logEntryDescription = aqLog[1][3]
  time_windowsTimestamp = aqLog[1][1][1]
  println("Win Timestamp detected: $time_windowsTimestamp")
  t=DateTime()
  if logEntryDescription =="Acquisition started" # "Acquisition started"
    tUnix = time_windowsTimestamp/(10.0*1000.0*1000.0)-11644473600.0
    println("tUnix has size $(size(tUnix))")
    t = Dates.unix2datetime(tUnix)#[1]
  else
    println("Error: could not find start time in $filename/AcquisitionLog/Log\n...trying to use timestamp string from attributes")
    =#
    attributesRoot = h5readattr(filename, "/")
    time_s = attributesRoot["HDF5 File Creation Time"]
    acq_card = attributesRoot["DAQ Hardware"]
    println(acq_card)
    if (acq_card == "Cronologic HPTDC8-PCI")
      t = DateTime(time_s, "m/d/y H:M:S") # STOF
    else
      t = DateTime(time_s, "d/m/y H:M:S") # PTR3
    end

  #end
  #f[:close]()
  return t
end

function validateHDF5Files(filepath, files)
  if debuglevel > 0   println("$nFiles files found, checking if valid.") end
  validFiles = []
  badFiles = []
  nFiles = size(files,1)

  for j=1:nFiles
     try
        totalPath = joinpath(filepath, files[j])
        if debuglevel > 1 print("Checking $totalPath ... ") end
        if (!ishdf5(totalPath))
          if debuglevel > 0   println("Bad File: $totalPath") end
          push!(badFiles,files[j])
        else
          fh = h5open(totalPath,"r")
          ds = fh["TimingData/BufTimes"]
          if (length(ds) > 0)
            if (ds[end,end][end,end] > 1e-99) # Last timestamp seems to be very small on corrupted files
              if debuglevel > 1 println("OK") end
              push!(validFiles,files[j])
            end
          else
            if debuglevel > 0   println("Bad File: $totalPath") end
            push!(badFiles,files[j])
          end
          close(fh)
        end
    catch
        println("File seems corrupt: $(files[j])")
    end
  end
  if debuglevel > 0   println("$(length(files)-length(validFiles)) files removed.") end
  return validFiles
end



############# Mass recalibration stuff ##############################

m_regionsToMatch = []
m_regionMaxMatchCoeffs = []
m_calibRegions = []
m_searchWidth = 0
m_referenceMassScaleMode = 0
m_referenceMassScaleParameters = []
m_referenceSpectrum = []

#plot stuff
crIndStart = 0
crIndEnd = 0
crMassIndicesOriginal = []
crOriginalMasses = []
m_plotControlMass = false

function setMassScaleReferenceSpectrum(referenceSpectrum, calibRegions, searchWidth, referenceMassScaleMode, referenceMassScaleParameters; plotControlMass=false, testRangeStart=0, testRangeEnd=0)
  println("Setting mass scale reference spectrum.")
  global m_regionsToMatch = []
  global m_regionMaxMatchCoeffs = []
  global m_referenceSpectrum = referenceSpectrum
  global m_calibRegions = calibRegions
  println("calibRegions: $m_calibRegions")

  global m_searchWidth = searchWidth
  global m_referenceMassScaleMode = referenceMassScaleMode
  global m_referenceMassScaleParameters = referenceMassScaleParameters

  global m_plotControlMass = plotControlMass

  for region in m_calibRegions
    referenceIndexStart::Int64 = round(mass2timebin(region -  searchWidth, referenceMassScaleMode,referenceMassScaleParameters))
    referenceIndexEnd::Int64 = round(mass2timebin(region +  searchWidth, referenceMassScaleMode,referenceMassScaleParameters))
    regionToMatch = m_referenceSpectrum[referenceIndexStart:referenceIndexEnd]

    push!(m_regionsToMatch,  regionToMatch)

    # calc self correlation coeff max for later normalization of correlation quality
    correlation = (xcorr(regionToMatch,regionToMatch))
    maximumCorrelation = findmax(correlation)
    intensity = maximumCorrelation[1]
    push!(m_regionMaxMatchCoeffs,intensity)

  end
  if m_plotControlMass
    global crIndStart = round(mass2timebin(testRangeStart,referenceMassScaleMode,referenceMassScaleParameters))
    global crIndEnd = round(mass2timebin(testRangeEnd,referenceMassScaleMode,referenceMassScaleParameters))
    global crMassIndicesOriginal = collect(crIndStart:1:crIndEnd)
    global crOriginalMasses = zeros(size(crMassIndicesOriginal))
    for i = 1 : length(crOriginalMasses)
      crOriginalMasses[i] = timebin2mass(crMassIndicesOriginal[i],referenceMassScaleMode,referenceMassScaleParameters)
    end
  end
end

function recalibrateMassScale(spectrum, referenceSpectrum, calibRegions, searchWidth, massCalibMode, massCalibParameters)
#println("calibRegions: $m_calibRegions")
print(".")
success = true

A = SharedArray{Float64}(length(m_calibRegions),2)
A[:,2] = 1
B = SharedArray{Float64}(length(m_calibRegions),1)

timebinshifts = zeros(length(m_calibRegions))
intensities = zeros(length(m_calibRegions))

#@sync @parallel
for regionindex=1:length(m_calibRegions)
    region = m_calibRegions[regionindex]


    indexStart::Int64 = round(mass2timebin(region - m_searchWidth, m_referenceMassScaleMode,m_referenceMassScaleParameters))
    indexEnd::Int64 = round(mass2timebin(region + m_searchWidth, m_referenceMassScaleMode,m_referenceMassScaleParameters))
    regionToSearch = spectrum[indexStart:indexEnd]
    correlation = (xcorr(convert(Array{Float64,1}, regionToSearch),convert(Array{Float64,1}, m_regionsToMatch[regionindex])))
    maximumCorrelation = findmax(correlation)
    if debuglevel > 3   println("maximumCorrelation at $(maximumCorrelation)") end
    shift = 0;
    intensity = 0;
    if (maximumCorrelation[2]>1 && maximumCorrelation[2]<length(correlation))
      inta=correlation[maximumCorrelation[2]-1]
      intb=maximumCorrelation[1]
      intc=correlation[maximumCorrelation[2]+1]
      minAC = min(inta, intc)
      shiftDeltaInterpolated = (intc - inta)/(inta + intb + intc - 3*minAC) # Tricky: -0.5 if a=b, +0.5 if b=c, +0.0 if a=c, interpol in between
      #plot(correlation[maximumCorrelation[2]-10:maximumCorrelation[2]+10])
      if debuglevel > 3   println("Interpolated Shift: $shiftDeltaInterpolated") end
      shift = shiftDeltaInterpolated + maximumCorrelation[2] - (length(regionToSearch) + length(m_regionsToMatch[regionindex]))/2
      intensity = maximumCorrelation[1]/m_regionMaxMatchCoeffs[regionindex]
      timebinshifts[regionindex] = shift
      intensities[regionindex] = intensity
    else
      timebinshifts[regionindex] = 0
      intensities[regionindex] = 0

    end
    if debuglevel > 2   println("Mass $region found shifted by $shift timebins with correlation coeff $(intensity)") end
    if (intensity < 0.05)
      success = false
      println("Could not correctly match m$region")
    end
    A[regionindex,1] = sqrt(region)
    B[regionindex] = mass2timebin(region, m_referenceMassScaleMode,m_referenceMassScaleParameters)  + shift  + 0.5
  end

  if success
    if debuglevel > 3   println("A: $A") end
    if debuglevel > 3   println("B: $B") end
    newParams = \(A,B)
    if debuglevel > 3   println("New parameters: $(newParams)") end
  else
    newParams = m_referenceMassScaleParameters
    if debuglevel > 3   println("Using WRONG calib parameters.") end
  end

  if (m_plotControlMass == true)
    crNewInterpolatedValues = zeros(crOriginalMasses)
    indexesExact = mass2timebin(crOriginalMasses, m_referenceMassScaleMode, newParams)
    crNewInterpolatedValues = interpolate(indexesExact,spectrum)
    plot(crOriginalMasses, crNewInterpolatedValues/maximum(crNewInterpolatedValues),".-")
  end
return newParams, success, timebinshifts, intensities
end

end
