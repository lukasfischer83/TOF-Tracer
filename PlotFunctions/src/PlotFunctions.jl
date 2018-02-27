__precompile__()
module PlotFunctions
using MasslistFunctions
using PyPlot

export massDefectPlot

function massDefectPlot(masses, compositions, concentrations, colors, plotTitle, colorCodeTitle; dotSize = 10, maxMass = 450, maxDefect = 0.25, minConc = 0.02, sumformulas = false)
  figure()

  h2o = createCompound(H=2,O=1, Hplus=0)
  o = createCompound(O=1, Hplus=0)

  for i=1:length(masses)
    m = masses[i]
    adduct = h2o
    if (inCompositions(compositions[:,i] + adduct[4], compositions))
      m1 = m + adduct[1]
      plot([m, m1], [m-round(m), m1 - round(m1)], color="lightblue", zorder=1)
    end
    adduct = o
    if (inCompositions(compositions[:,i] + adduct[4], compositions))
      m1 = m + adduct[1]
      plot([m, m1], [m-round(m), m1 - round(m1)], color="red", zorder=1)
    end
    if sumformulas text(masses[i],masses[i]-round(masses[i]),sumFormulaStringFromCompositionArray(compositions[:,i]), color="grey", clip_on=true, verticalalignment="center", size=10, zorder=100) end
  end

  scatter(masses, masses-round.(masses),dotSize*log.(concentrations./minConc), colors, zorder=10, linewidths=0.5)
  xlim(0,maxMass)
  ylim(0,maxDefect)
  cb=colorbar()
  cb["ax"]["set_ylabel"](colorCodeTitle)
  xlabel("Mass [amu]")
  ylabel("Kendrick Mass Defect")
  title(plotTitle)
  grid("on")
end

end
