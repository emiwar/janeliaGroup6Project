import MAT
include("megamap.jl")


placeWidth = 50.0
learningRate = 1e-7


function convToJulia(arr)
    if length(arr) == 0
        return Vector{Float64}()
    elseif length(arr) == 1
        return [arr]
    else
        return Vector{Float64}(arr[:,1])
    end
end

vars = MAT.matread("place_field_data.mat")
placeCenters = convToJulia.(vars["PF_centers"])[:,1];

targetMap = TargetMap(placeCenters; placeWidth=placeWidth)
evalPoints = 50:10:3950
W = computeW(targetMap, evalPoints, -learningRate; maxIter=200000, initSigma=0.0)

MAT.matwrite("fittedW_1dm_grid_many_iters.mat", Dict("fittedW" => W, "evaluationPoints" => collect(evalPoints), "placeWidth" => placeWidth, "learningRate" => learningRate, "placeCenters" => placeCenters, "maxIter" => 200000))
