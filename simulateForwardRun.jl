using MAT
include("megamap.jl")


function convToJulia(arr)
    if length(arr) == 0
        return Vector{Float64}()
    elseif length(arr) == 1
        return [arr]
    else
        return Vector{Float64}(arr[:,1])
    end
end

predictedW = matread("fittedW_1dm_grid_noinput.mat");

vars = matread("place_field_data.mat")
placeCenters = convToJulia.(vars["PF_centers"])[:,1];


targetMap = TargetMap(placeCenters; placeWidth=50.0);

W = predictedW["fittedW"]
nc::Int64 = nCells(targetMap)

allAllF = Vector{Matrix{Float64}}()
for lap = 1:8
    println("Lap $(lap)")
    forwardMap = ForwardMap(targetMap.fPeak, targetMap.inhibThres, targetMap.wI, W, zeros(nc));
    allF = Matrix{Float32}(2000*length(50:10:3950), 800)#Vector{Matrix{Float64}}()
    for (i,x)=enumerate(50:10:3950)
        inp = input(x, targetMap);
        allF[i*2000 - 1999:i*2000, :] = simulate!(forwardMap, inp; timesteps=2000, noise_s=0.05);
    end
    push!(allAllF, allF)
end

meanF = mean(allAllF)

matwrite("noiseCorrelationSimulatedRun.mat", Dict("cor" => cor(allAllF[end] - meanF)))
