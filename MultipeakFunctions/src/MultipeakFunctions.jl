__precompile__()
module MultipeakFunctions
using PeakshapeFunctions
using MasslistFunctions
export calculateCrossTalkMatrix, reconstructSpectrum, sumBins

function calculateCrossTalkMatrix(massAxis, binWidth, masses, masslistElements, compositions, peakShapesCenterMass, peakShapesY)
  mtrx = Array(Float64,length(masses), length(masses))
  fill(mtrx,0)
  centerindices = Array(Float64,0)
  for i = 1:length(masses)
    peakShapeSize = size(peakShapesY,1)
    indexOffset::Int64 = round(peakShapeSize+1)/2
    centerIndex = searchsortedfirst(massAxis, masses[i])#-1
    push!(centerindices, centerIndex)
    nPointsBeforePeakMax = centerIndex - searchsortedfirst(massAxis,masses[i]-1)
    if sum(compositions[:,i]) > 0
      isotopeMasses, isotopeMasslistElements, isotopeCompositions, isotopeAbundances = isotopesFromCompositionArray(compositions[:,i])
    else
      println("Unidentified Peak at mass $(masses[i])")
      isotopeMasses = [masses[i]]
      isotopeMasslistElements = [masslistElements]
      isotopeCompositions = [compositions[:,i]]
      isotopeAbundances = [1]
    end

    nLocalPeakPatternPoints = searchsortedfirst(massAxis, masses[i]+3.5)+ - searchsortedfirst(massAxis, masses[i]-1)
    localPeakPattern = zeros(nLocalPeakPatternPoints)
    for isotope = 1:length(isotopeMasses)
      isotopeCenterIndex = searchsortedfirst(massAxis, isotopeMasses[isotope])-1
      nPointsBeforeLocalPeakMax = isotopeCenterIndex - searchsortedfirst(massAxis,masses[i]-1)
      localPeakshape = getLocalPeakshape(isotopeMasses[isotope], peakShapesCenterMass, peakShapesY)
      localPeakshape = localPeakshape * isotopeAbundances[isotope]
      localPeakPattern[nPointsBeforeLocalPeakMax-indexOffset+1:nPointsBeforeLocalPeakMax+length(localPeakshape)-indexOffset] += localPeakshape
    end
    for j = 1:length(masses)
        if ((masses[j] > masses[i]-1) && (masses[j] < masses[i]+3.5))
        dist = searchsortedfirst(massAxis,masses[j]) - centerIndex
        relativeStartIdx = (nPointsBeforePeakMax + dist - binWidth)
        relativeEndIdx = (nPointsBeforePeakMax + dist + binWidth)
        if (relativeStartIdx > 0) && (relativeEndIdx < length(localPeakPattern))
          overlap = sum(localPeakPattern[relativeStartIdx : relativeEndIdx])
          mtrx[j,i] = overlap
        else
            mtrx[j,i] = 0
        end
      else
          mtrx[j,i] = 0
    end
    end
  end
  return mtrx
end

function sumBins(massAxis, binWidth, spectrum, masses)
  stickRaw = Array(Float64,length(masses))
  fill(stickRaw,0)
  for i=1:length(masses)
    midx = searchsortedfirst(massAxis, masses[i])-1
    stickRaw[i] = sum(spectrum[midx-binWidth:midx+binWidth])
  end
  return stickRaw
end

function reconstructSpectrum(massAxis, masses, masslistElements, compositions, counts, peakShapesCenterMass, peakShapesY)

  reconstructedSpectrum = zeros(length(massAxis))
  fill(reconstructedSpectrum, 0)
  for i=1:length(masses)
    if sum(compositions[:,i]) > 0
      isotopeMasses, isotopeMasslistElements, isotopeCompositions, isotopeAbundances = isotopesFromCompositionArray(compositions[:,i])
    else
      println("Unidentified Peak at mass $(masses[i])")
      isotopeMasses = [masses[i]]
      isotopeMasslistElements = [masslistElements]
      isotopeCompositions = [compositions[:,i]]
      isotopeAbundances = [1]
    end
    for isotope = 1:length(isotopeMasses)
      centerIndex = searchsortedfirst(massAxis, isotopeMasses[isotope])
      lps = 0
      lps = getLocalPeakshape(isotopeMasses[isotope], peakShapesCenterMass, peakShapesY)
      reconstructedSpectrum[centerIndex-Int64((length(lps)-1)/2):centerIndex+Int64((length(lps)-1)/2)] += lps*counts[i]*isotopeAbundances[isotope]
    end
  end
  return reconstructedSpectrum
end


end
