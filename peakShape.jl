push!(LOAD_PATH, pwd())

using PyPlot
using HDF5
using InterpolationFunctions
using BaselineFunctions
using PeakshapeFunctions

function baselineAndPeakshape(
  filepath;
  outputfilename="results/_result.hdf5",
  peakshapeRegions=10,
  peakshapeQuantileValue = 0.05
  )

  file = joinpath(filepath, outputfilename)
  massAxis = h5read(file, "MassAxis")
  avgSpectrum = h5read(file, "AvgSpectrum")

  baselinePoints, baselineValues, baselineNoise = calculateBaseline(massAxis,avgSpectrum,baselinePointWidth = 0.8, threshold=0.2)
  baselineNoiseInterpolated = interpolate(massAxis, baselinePoints, baselineNoise)
  baselineInterpolated = interpolate(massAxis, baselinePoints, baselineValues)
  baselineCorrectedAvgSpec = avgSpectrum[:,1] - baselineInterpolated;


  peakIndices = findPeakIndices(massAxis, avgSpectrum, baselineInterpolated, baselineNoiseInterpolated)
  peakShapesCenterMass, peakShapesY = calculatePeakshapes(massAxis, baselineCorrectedAvgSpec, peakIndices, nbrMassRegions = peakshapeRegions, quantileValue = peakshapeQuantileValue)

  figure()
  semilogy(baselinePoints,baselineValues,".-")
  semilogy(baselinePoints,baselineValues + baselineNoise,".")
  #semilogy(peakMasses, peakValues, "x")
  semilogy(massAxis,avgSpectrum)


  ############ delete h5 data that will be overwritten ###########
  fh = h5open(file,"r+")
  if exists(fh, "AvgBaseline")
    o_delete(fh,"AvgBaseline")
  end
  if exists(fh, "MassDepPeakshape")
  o_delete(fh,"MassDepPeakshape")
  end
  if exists(fh, "MassDepPeakshapeCenterMasses")
  o_delete(fh,"MassDepPeakshapeCenterMasses")
  end
  if exists(fh, "BaseLines")
  o_delete(fh,"BaseLines")
  end
  close(fh)


  h5write(file, "AvgBaseline", baselineInterpolated)
  h5write(file, "MassDepPeakshape", peakShapesY)
  h5write(file, "MassDepPeakshapeCenterMasses", peakShapesCenterMass)

  fh = h5open(file,"r")
  ds = fh["SumSpecs"]
  spectra = size(ds)[2]
  close(fh)
  if (spectra > 0)
    SumSpecs = h5read(file, "SumSpecs")
    fill!(SumSpecs,0)
    h5write(file, "BaseLines", SumSpecs)

  else
    h5write(file, "SumSpecs", avgSpectrum)
    fill!(avgSpectrum,0)
    h5write(file, "BaseLines", avgSpectrum)
  end

  println("DONE")
end
