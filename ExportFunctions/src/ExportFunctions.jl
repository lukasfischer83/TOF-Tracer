__precompile__()
module ExportFunctions
using InterpolationFunctions
using MasslistFunctions

export exportTracesCSV

function exportTracesCSV(saveFolderPath, elementNames, compositions, times, traces; average=0)
  sumformulas = sumFormulaStringListFromCompositionArrayList(compositions)
  f = open("$saveFolderPath/ptr3compositions.txt", "w")
  writedlm(f, hcat(["Mass" "SumFormula"],reshape(elementNames,(1,length(elementNames)))))
  writedlm(f, hcat(massFromCompositionArrayList(compositions),sumformulas , compositions'))
  close(f)
  f = open("$saveFolderPath/ptr3traces.csv", "w")
  writedlm(f, hcat(["Time"], reshape(sumformulas,(1,length(sumformulas)))))
  if (average==0)
    writedlm(f, hcat(times ,traces))
  else
    writedlm(f, hcat(averageSamples(times,average) ,averageSamples(traces,average)))
  end
  close(f)
end

end
