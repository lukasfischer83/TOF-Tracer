__precompile__()
module ExportFunctions
using InterpolationFunctions
import  MasslistFunctions

export exportTracesCSV, exportTracesCSVLossCorr, toMatlabTime, fromMatlabTime

function exportTracesCSV(saveFolderPath, elementNames, compositions, times, traces; average=0)
  sumformulas = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(compositions)
  f = open("$saveFolderPath/ptr3compositions.txt", "w")
  writedlm(f, hcat(["Mass" "SumFormula"],reshape(elementNames,(1,length(elementNames)))))
  writedlm(f, hcat(MasslistFunctions.massFromCompositionArrayList(compositions),sumformulas , compositions'))
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

function exportTracesCSVLossCorr(saveFolderPath, elementNames, compositions, times, traces, lossfactor, lossfactorerr, corrfactor, corrfactorerr, corrnotes; average=0)
  sumformulas = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(compositions)
  f = open("$saveFolderPath/ptr3compositions.txt", "w")
  writedlm(f, hcat(["Mass"	"SumFormula"],reshape(elementNames,(1,length(elementNames))),["LossFactor"	"LossFactorError"	"CorrFactor"	"CorrFactorErr"	"CorrNotes"]))
  writedlm(f, hcat(MasslistFunctions.massFromCompositionArrayList(compositions),sumformulas , compositions', lossfactor, lossfactorerr, corrfactor, corrfactorerr, corrnotes))
  close(f)
  f = open("$saveFolderPath/ptr3tracesInletLossCorr.csv", "w")
  writedlm(f, hcat(["Time"], reshape(sumformulas,(1,length(sumformulas)))))
  if (average==0)
    writedlm(f, hcat(times , (corrfactor.*traces' )' ))
  else
    writedlm(f, hcat(averageSamples(times,average) ,(corrfactor.*(averageSamples(traces,average))' )' ))
  end
  close(f)
end


function toMatlabTime(t::DateTime)
    timespan = ((t+Dates.Day(1)) - DateTime(0,1,1,0,0,0))
    millis::Float64 = Dates.Millisecond(timespan)
    return millis/24/3600/1000
end

function fromMatlabTime(timestamp::Number)
    days=Int(floor(timestamp))
    millisecondsRemainder = Int(round((timestamp-days)*24*3600*1000))
    return DateTime(0,1,1,0,0,0)+Dates.Day(days-1)+Dates.Millisecond(millisecondsRemainder)+Dates.Year(1)
end

end
