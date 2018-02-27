using HDF5
using PyPlot
using MultipeakFunctions
using MasslistFunctions
using ResultFileFunctions

#include("masslistFunctions.jl")

function deconvolute(
  filepath;
  outputfilename="results/_result.hdf5",
  binWidth = 6
  )
  #binWidth += 1 # Maybe crosstalk is calculated from "<" while summing is done from "<=" ?? Gives better crostalk removal.
  file = joinpath(filepath, outputfilename)

  peakShapesY = h5read(file, "MassDepPeakshape")
  peakShapesY = peakShapesY./sum(peakShapesY,1) # Normalize!

  peakShapesCenterMass = h5read(file, "MassDepPeakshapeCenterMasses")
  totalAvgSpectrum = h5read(file, "AvgSpectrum") - h5read(file,"AvgBaseline")
  massAxis = h5read(file, "MassAxis")
  masslistElements = h5read(file,"ElementNames")
  compositionsOrig = h5read(file, "ElementalCompositions")
  massesOrig = h5read(file, "MassList")
  if unique(massesOrig) != massesOrig
    println("Multiple entries of the same mass --> will produce singular matrix!!!")
  end

  figure()
  ax = subplot(111)

    selector = (massesOrig .> 0)
    masses = massesOrig[selector]
    compositions = compositionsOrig[:,selector]




    print("Populating matrix for inversion of linear system...")
    tic()
    mtrx = calculateCrossTalkMatrix(massAxis, binWidth, masses, masslistElements, compositions, peakShapesCenterMass, peakShapesY)
    stickRaw = sumBins(massAxis,binWidth,totalAvgSpectrum,masses)
    toc()
    println(" DONE")



    sparseMtrx = mtrx
    print("Inverting Matrix...")
    tic()
    deconvolutionMatrix = inv(sparseMtrx)
    toc()
    println(" DONE")

    print("Applying deconvolution kernel...")
    tic()
    counts = deconvolutionMatrix * stickRaw
    toc()
    println(" DONE")

  print("Reconstructing Spectrum for visual check...")
  tic()

  semilogy(massAxis,totalAvgSpectrum, "-o", label="Original", color="r")
  #reconstructedSpectrum = reconstructSpectrum(massAxis, masses[(masses.>158) & (masses.<162)], masslistElements, compositions[:,(masses.>158) & (masses.<162)], counts[(masses.>158) & (masses.<162)], peakShapesCenterMass, peakShapesY)
  reconstructedSpectrum = reconstructSpectrum(massAxis, masses, masslistElements, compositions, counts, peakShapesCenterMass, peakShapesY)

  plot(massAxis, reconstructedSpectrum, label="Fit", color="b")
  plot(massAxis, totalAvgSpectrum-reconstructedSpectrum, label="Residual", color="g")
  #=
  fittedPeaks = Array{Float64}(length(masses),2001)
  for i=1:length(masses)
    approxMassIndex = searchsortedfirst(massAxis,masses[i])
    fittedPeaks[i,:] = reconstructSpectrum(massAxis[approxMassIndex-300 : approxMassIndex+1700], masses[i], masslistElements, compositions[:,i], counts[i], peakShapesCenterMass, peakShapesY)
    plot(massAxis[approxMassIndex-300 : approxMassIndex+1700],fittedPeaks[i,:],"--", color="green")
  end
  =#

  legend()
  ax[:set_ylim]([minimum(totalAvgSpectrum),maximum(totalAvgSpectrum)])
  toc()
  println(" DONE")



  ############ WRITE OUTPUT TO FILE ##############################################
  haveStickCps = false
  # Correct timetraces

  fh = h5open(file,"r+")

  if exists(fh, "CorrStickCps")
  o_delete(fh,"CorrStickCps")
  end
  if exists(fh, "CorrStickCpsErrors")
  o_delete(fh,"CorrStickCpsErrors")
  end

  if exists(fh, "StickCps")
    haveStickCps = true

    # Create empty Dataspace
    nbrSpectra = getNbrTraceSamples(file)
    dset = d_create(fh, "CorrStickCps", datatype(Float32), dataspace(nbrSpectra, length(masses)), "chunk", (1,length(masses)))

    toProcessLow = 0
    toProcessHigh = 0
    while toProcessHigh < nbrSpectra
      toProcessLow = toProcessHigh + 1
      toProcessHigh = toProcessLow + 9999
      if toProcessHigh > nbrSpectra
        toProcessHigh = nbrSpectra
      end
      println("Correcting spectrum $toProcessLow to $toProcessHigh of $nbrSpectra")
      samplesSubRange = convert(SharedArray, getTraceSamples(file,toProcessLow:toProcessHigh, raw=true)[:,selector])

      for i=1:(toProcessHigh - toProcessLow + 1)
        #traces[i,:] = deconvolutionMatrix * traces[i,:]
        #tracesErrors[i,:] = abs(deconvolutionMatrix) * sqrt(abs(traces[i,:])/5)
        dset[toProcessLow - 1 + i,:] = deconvolutionMatrix *samplesSubRange[i,:]
      end
    end
    #h5write(file, "CorrStickCps", convert(Array,traces))
    #h5write(file, "CorrStickCpsErrors", tracesErrors)
  end
  close(fh)

  fh = h5open(file,"r+")
  if exists(fh, "AvgStickCps")

    traces = h5read(file, "AvgStickCps")[:,selector]
    tracesErrors = similar(traces)
    for i=1:size(traces,1)
      #traces[i,:] = deconvolutionMatrix \ traces[i,:]
      traces[i,:] = deconvolutionMatrix * traces[i,:]
      #tracesErrors[i,:] = abs(deconvolutionMatrix) * sqrt(abs(traces[i,:])/3600)
  end

  if exists(fh, "CorrAvgStickCps")
  o_delete(fh,"CorrAvgStickCps")
  end
  if exists(fh, "CorrAvgStickCpsErrors")
  o_delete(fh,"CorrAvgStickCpsErrors")
  end
  h5write(file, "CorrAvgStickCps", traces)
  h5write(file, "CorrAvgStickCpsErrors", tracesErrors)
  end
  close(fh)


  #if haveStickCps
  #  traces = h5read(file, "CorrStickCps")[:,selector]
  #else
  #  traces = h5read(file, "CorrAvgStickCps")[:,selector]
  #end
end
