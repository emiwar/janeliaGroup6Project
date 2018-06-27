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

trialId = 20

vars = MAT.matread("20\ Poisson\ Trials\ with\ different\ n.mat")
placeCenters = convToJulia.(vars["poisson_data"][trialId])[:,1];

targetMap = TargetMap(placeCenters; placeWidth=placeWidth)
evalPoints = 50:10:3950
W = computeW(targetMap, evalPoints, -learningRate; maxIter=30000, initSigma=1e-3)

MAT.matwrite("poissonPlaceFields/trial20.mat", Dict("fittedW" => W, "evaluationPoints" => collect(evalPoints), "placeWidth" => placeWidth, "learningRate" => learningRate, "placeCenters" => placeCenters, "maxIter" => 30000))
