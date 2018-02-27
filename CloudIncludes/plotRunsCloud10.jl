# Adds lines from the Runlog between startTime and endTime to the current plot
function plotStages(startTime, endTime)
  ax = PyPlot.gca()
  runlist_array = readdlm("/data/CLOUDX/runlist.txt",'\t', skipstart=1)
  runlist_times = [(Dates.unix2datetime(runlist_array[i,1]) - Dates.Hour(1)) for i = 1:size(runlist_array,1)]
  sel = (runlist_times .> startTime) & (runlist_times .< endTime)
  runlist_times = runlist_times[sel]
  runlist_array = runlist_array[sel,:]
  tform = PyPlot.matplotlib[:transforms][:blended_transform_factory](ax[:transData], ax[:transAxes])
  for i = 1:length(runlist_times)
    if (runlist_array[i,2] == "sor")
      PyPlot.axvline(x=runlist_times[i], linewidth=0.5, color="grey")
      PyPlot.text(runlist_times[i],0.95,"\n$(runlist_array[i,3]) $(runlist_array[i,5])",transform=tform, rotation=90, color="grey", clip_on=true)
    end
  end
end
