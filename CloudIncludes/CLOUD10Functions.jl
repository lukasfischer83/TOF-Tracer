using HDF5

function cloud10transitions()
    return [
    DateTime(2015,09,27,08,00),
    DateTime(2015,10,05,06,00),
    DateTime(2015,10,08,10,00),
    DateTime(2015,10,15,4,00)
    ]
end

function cloud10RHRampPeriods()
  return [
  DateTime(2016,10,13,22) DateTime(2016,10,14,22,22)
  DateTime(2016,10,18,14,41) DateTime(2016,10,19,4,44)
  DateTime(2016,10,21,15,30) DateTime(2016,10,22,5)
  ]
end

function readFixedColumnData(filename, indices; skipstart=0)
    nbrColumns = size(indices,1)
    results = Array{Float64,2}(0,nbrColumns)
    l=0
    open(filename) do io
        for i=0:skipstart
            line = readline(io)
            l+=1
        end
        while true
            line = readline(io)
            l+=1
            line == "" && break
            if length(line) >= indices[end][end]
                tmp = Array{Nullable{Float64},2}(1,nbrColumns)
                for i=1:nbrColumns
                    tmp[i] = tryparse(Float64, strip(line[indices[i][1]:indices[i][2]]))
                end
                try
                    results = vcat(results, get.(tmp))
                catch
                    println("Could not parse file $filename, line number $(l)!")
                end
            end
       end
    end
    results
end

function loadDEWPOINT(directoryPath; filefilterRegexp = r"\.data$")
      println("Loading DEWPOINT Data...")
      files = filter(filefilterRegexp, readdir(directoryPath))
      times = Array(Float64,0)
      DEW = Array(Float64,0)

      for i=1:length(files)
        println("Importing file $(files[i])")
        tmp = readFixedColumnData(directoryPath * "/"  * files[i], [(1,11),(12,19)], skipstart=1)
        append!(times,tmp[:,1])
        append!(DEW,tmp[:,2])
      end

      if !issorted(times)
        println("Times not in ascending order --> sorting data")
        sortedOrder = sortperm(times)
        times = times[sortedOrder]
        DEW=DEW[sortedOrder]
    end

    outpath=directoryPath * "/dew.hdf5"
    if (isfile(outpath))
      mv(outpath, directoryPath * "/dew-bak.hdf5", remove_destination=true)
    end
    HDF5.h5write(outpath, "times", times)
    HDF5.h5write(outpath, "dewpoints", DEW)

    times = Dates.unix2datetime.(times)
    return times, DEW
end

function loadDEWPOINT_cached(directoryPath)
    outpath=directoryPath * "/dew.hdf5"

    t = HDF5.h5read(outpath, "times")
    d = HDF5.h5read(outpath, "dewpoints")
    return Dates.unix2datetime.(t), d
end
