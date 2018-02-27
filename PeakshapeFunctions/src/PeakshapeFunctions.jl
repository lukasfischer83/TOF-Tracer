__precompile__()
module PeakshapeFunctions
using InterpolationFunctions
using PyPlot

export findPeakIndices, calculatePeakshapes, getLocalPeakshape
"searches an average spectrum for peaks signifficantly higher than the baseline noise and returns their interpolated indices in the average spectrum"
  function findPeakIndices(massAxis, avgSpectrum, baseline, baselineNoise; noiseThreshold = 80, oddEven="both")
    totalMax = maximum(avgSpectrum)
    peakIndices = Array{Float64}(1)
    for j=searchsortedfirst(massAxis,10):length(massAxis)-1
      if ( avgSpectrum[j] > baseline[j] + baselineNoise[j] * noiseThreshold * 400/(massAxis[j]^1.3)) && (avgSpectrum[j] > avgSpectrum[j-1]) && (avgSpectrum[j] > avgSpectrum[j+1]) && (avgSpectrum[j] < totalMax*0.1)
          if (oddEven == "both") | ((oddEven == "even") & iseven(Int(round(massAxis[j],0)))) | ((oddEven == "odd") & !iseven(Int(round(massAxis[j],0))))
              ############ put interpolated exact peak position in peak list ###
              push!(peakIndices,interpolatedMax(j,avgSpectrum))
          end
      end
    end
    return peakIndices
  end

  function calculatePeakshapes(massAxis, baselineCorrectedAvgSpec, peakIndices; nbrMassRegions = 10, peakWindowWidth = 200, quantileValue = 0.05)
    peakMasses = interpolate(peakIndices, massAxis)
    peakValues = interpolate(peakIndices,baselineCorrectedAvgSpec)

    peakShapesY = Array{Float64}(2*peakWindowWidth+1, nbrMassRegions)
    peakShapesCenterMass = Array(Float64,nbrMassRegions)

    figure()

    for massRegion = 1:nbrMassRegions
      ax=subplot(ceil(nbrMassRegions/4),4,massRegion)

      #peakshapeRangeStart = (massRegion-1) *peakMasses[end] / nbrMassRegions
      #peakshapeRangeEnd = (massRegion) *peakMasses[end] / nbrMassRegions
      peakshapeRangeStart = (massRegion-1)^2 *peakMasses[end] / nbrMassRegions^2
      peakshapeRangeEnd = (massRegion)^2 *peakMasses[end] / nbrMassRegions^2

      peakShapesCenterMass[massRegion] = (peakshapeRangeEnd + peakshapeRangeStart)/2
      sel = (peakMasses.>peakshapeRangeStart) & (peakMasses.<peakshapeRangeEnd) & (peakIndices.> peakWindowWidth) & (peakIndices .< length(baselineCorrectedAvgSpec)-peakWindowWidth)

      peakIndicesInRegion = peakIndices[sel]
      peakMassesInRegion = peakMasses[sel]
      peakValuesInRegion = peakValues[sel]

      peakWood = SharedArray(Float64,2*peakWindowWidth+1,length(peakMassesInRegion))
      fill(peakWood, 0)
      for i=1:length(peakMassesInRegion)
          ############# interpolate exact peak position ###################
        lowIdx = Int64(floor(peakIndicesInRegion[i]))
        highIdx = lowIdx+1
        fract = peakIndicesInRegion[i]-lowIdx

        peakWood[:,i] = (
        (1-fract)*view(baselineCorrectedAvgSpec, lowIdx-peakWindowWidth : lowIdx+peakWindowWidth) +
        fract*view(baselineCorrectedAvgSpec, highIdx-peakWindowWidth : highIdx+peakWindowWidth))/peakValuesInRegion[i]
      end

      peakshapeY = Array(Float64,2*peakWindowWidth+1)
      for i=1:2*peakWindowWidth+1
        #peakshapeY[i] = minimum(peakWood[i,:])
        values = peakWood[i,:]
        #upperThreshold = quantile(values,quantileValue*peakshapeRangeEnd/100)
        ###### quantile based #############
        #upperThreshold = quantile(values,quantileValue*peakshapeRangeEnd/100)
        #peakshapeY[i] = mean(values[values.<upperThreshold])
        #lowest n perent based ############
        values = sort(values[values.>0])
        if length(values) > 0
          peakshapeY[i] = mean(values[1 : Int64(ceil(length(values)*quantileValue))])
        else
          peakshapeY[i] = 0
          if i == peakWindowWidth + 1
            peakshapeY[i] = 1
          end
        end
      end
      peakshapeY[peakshapeY.<0] = 0
      peakshapeY[peakWindowWidth+1] = 1
      peakShapesY[:,massRegion] = peakshapeY / sum(peakshapeY)
      if (size(peakWood,2) > 0)
          semilogy(peakWood)
      end
      semilogy(peakshapeY, "o-", label="$peakshapeRangeStart - $peakshapeRangeEnd")
      ax["set_title"]("$peakshapeRangeStart - $peakshapeRangeEnd")
    end
    return peakShapesCenterMass, peakShapesY
  end

  function getLocalPeakshape(mass, peakShapesCenterMass, peakShapesY)
    peakShapeCenterMassIndexHigh = searchsortedfirst(peakShapesCenterMass, mass)
    if (peakShapeCenterMassIndexHigh == 1)
      localPeakshape = peakShapesY[:,1]
    elseif (peakShapeCenterMassIndexHigh > length(peakShapesCenterMass))
      localPeakshape = peakShapesY[:,end]
    else
      peakShapeCenterMassIndexLow = peakShapeCenterMassIndexHigh - 1
      peakShapesCenterMassDistance = peakShapesCenterMass[peakShapeCenterMassIndexHigh] - peakShapesCenterMass[peakShapeCenterMassIndexLow]
      contribLow =  (peakShapesCenterMass[peakShapeCenterMassIndexHigh] - mass) / peakShapesCenterMassDistance
      contribHigh = 1 - contribLow
      localPeakshape = peakShapesY[:,peakShapeCenterMassIndexLow]*contribLow + peakShapesY[:,peakShapeCenterMassIndexHigh]*contribHigh
    end
    return localPeakshape
  end

end
