function cloud11humidity()

  humidityRaw = readdlm("/data/CLOUD11/external/Humidity-Temperature/humidity-temp-avg.txt", '\t', skipstart = 1)
  timeHumUnix = humidityRaw[:,1]/1000
  hum = humidityRaw[:,2]
  timeHum = Array{DateTime}(length(timeHumUnix))
  for i=1:length(timeHumUnix)
    timeHum[i] = Dates.unix2datetime(timeHumUnix[i])
  end
  return timeHum, hum
end

function cloud11daqapinene()
  daqbrokerapineneRaw = readdlm("/data/CLOUD11/external/apinene-daqbroker.txt", '\t', skipstart = 1)
  timeDaqbrokerapineneUnix = daqbrokerapineneRaw[:,1]/1000
  daqbrokerapinene = daqbrokerapineneRaw[:,2]
  timeDaqbrokerapinene = Array{DateTime}(length(timeDaqbrokerapineneUnix))
  for i=1:length(timeDaqbrokerapinene)
    timeDaqbrokerapinene[i] = Dates.unix2datetime(timeDaqbrokerapineneUnix[i])
  end
  return timeDaqbrokerapinene, daqbrokerapinene
end

function cloud11daqTMB()
  daqbrokerapineneRaw = readdlm("/data/CLOUD11/external/TMB-daqbroker.txt", '\t', skipstart = 1)
  timeDaqbrokerapineneUnix = daqbrokerapineneRaw[:,1]/1000
  daqbrokerapinene = daqbrokerapineneRaw[:,2]
  timeDaqbrokerapinene = Array{DateTime}(length(timeDaqbrokerapineneUnix))
  for i=1:length(timeDaqbrokerapinene)
    timeDaqbrokerapinene[i] = Dates.unix2datetime(timeDaqbrokerapineneUnix[i])
  end
  return timeDaqbrokerapinene, daqbrokerapinene
end

function cloud11daqNPT()
  daqbrokerapineneRaw = readdlm("/data/CLOUD11/external/NPT-daqbroker.txt", '\t', skipstart = 1)
  timeDaqbrokerapineneUnix = daqbrokerapineneRaw[:,1]/1000
  daqbrokerapinene = daqbrokerapineneRaw[:,2]
  timeDaqbrokerapinene = Array{DateTime}(length(timeDaqbrokerapineneUnix))
  for i=1:length(timeDaqbrokerapinene)
    timeDaqbrokerapinene[i] = Dates.unix2datetime(timeDaqbrokerapineneUnix[i])
  end
  return timeDaqbrokerapinene, daqbrokerapinene
end

function cloud11transitions()
    return [
    DateTime(2016,9,27,12,38),
    DateTime(2016,9,29,18,59),
    DateTime(2016,10,2,16,53),
    DateTime(2016,10,12,11,0),
    DateTime(2016,10,14,14,56),
    DateTime(2016,10,17,8,57),
    DateTime(2016,10,18,17,19),
    DateTime(2016,10,21,17,48),
    DateTime(2016,10,22,13,0),
    DateTime(2016,11,4,11,0),
    DateTime(2016,11,4,23,31),
    DateTime(2016,11,5,12,20),
    DateTime(2016,11,5,15,3),
    DateTime(2016,11,16,4,44),
    DateTime(2016,11,16,20,16),
    DateTime(2016,11,17,5,20),
    DateTime(2016,11,17,13,6),
    DateTime(2016,11,19,9,6),
    DateTime(2016,11,23,8,53),
    DateTime(2016,11,24,0,25),
    DateTime(2016,11,24,12,4),
    DateTime(2016,11,25,16,33),
    DateTime(2016,11,26,4,12),
    DateTime(2016,11,26,18,26),
    DateTime(2016,11,28,1,30),
    DateTime(2016,11,28,11,51),
    DateTime(2016,11,29,20,13),
    DateTime(2016,12,1,17,31)
    ]
end

function cloud11RHRampPeriods()
  return [
  DateTime(2016,10,13,22) DateTime(2016,10,14,22,22)
  DateTime(2016,10,18,14,41) DateTime(2016,10,19,4,44)
  DateTime(2016,10,21,15,30) DateTime(2016,10,22,5)
  DateTime(2016,11,04,10,30) DateTime(2016,11,04,15,30)
  DateTime(2016,11,04,22,30) DateTime(2016,11,05,03,00)
  DateTime(2016,11,05,06,00) DateTime(2016,11,05,08,30)
  DateTime(2016,11,05,11,30) DateTime(2016,11,05,15,00)
  DateTime(2016,11,15,21) DateTime(2016,11,16,21)
  DateTime(2016,11,17,01) DateTime(2016,11,17,07,50)
  DateTime(2016,11,17,11) DateTime(2016,11,17,14)
  DateTime(2016,11,18,7) DateTime(2016,11,19,13,40)
  DateTime(2016,11,23,7) DateTime(2016,11,23,22)
  DateTime(2016,11,24,11) DateTime(2016,11,24,16,20)
  ]
end

function shrinkSel!(sel)
  for i=length(sel):-1:2
    if (sel[i]==true) && (sel[i-1] == false)
        println("Setting $i to false")
        sel[i] = false
    end
  end
end
function growSel!(sel)
    for i=length(sel):-1:2
      if (sel[i]==true) && (sel[i-1] == false)
          println("Setting $i to false")
          sel[i] = true
      end
    end
end
