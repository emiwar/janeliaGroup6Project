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

predictedW = matread("./poissonPlaceFields/trial1.mat")
placeCenters = convToJulia.(predictedW["placeCenters"]);
targetMap = TargetMap(placeCenters; placeWidth=50.0);
W = predictedW["fittedW"]


targets = Vector{Vector{Float64}}()
fEnds = Vector{Vector{Float64}}()
fIntermediates = Vector{Vector{Float64}}()
W = predictedW["fittedW"]
for x=50:10:3950
    forwardMap = ForwardMap(targetMap.fPeak, targetMap.inhibThres, targetMap.wI, W, zeros(size(W, 1)));
    target = fTarget(x, targetMap);
    inp = input(x, targetMap);
    fIntermediate = simulate!(forwardMap, inp; timesteps=1500, noise_s=0.00)
    fEnd = simulate!(forwardMap, zeros(size(W, 1)); timesteps=1500, noise_s=0.00);
    push!(targets, target)
    push!(fEnds, fEnd[end,:])
    push!(fIntermediates, fIntermediate[end,:])
end

nPoints = length(targets)
intermediateCorrs = [cor(fIntermediates[i], targets[j]) for i=1:nPoints, j=1:nPoints];
finalCorrs = [cor(fEnds[i], targets[j]) for i=1:nPoints, j=1:nPoints];

predErr = 50.0*[indmax(finalCorrs[:,i]) - i for i=1:size(finalCorrs, 1)]

MAT.matwrite("./poissonPlaceFields/simTrial1.mat", Dict("W" => W, "placeCenters" => predictedW["placeCenters"],
                                    "placeWidth" => 50.0, "fEnds" => hcat(fEnds...), "fIntermediates" => hcat(fIntermediates...),
                                    "targets" => hcat(targets...), "xPos"=> collect(50:10:3950),
                                    "intermediateCorrs" => intermediateCorrs,
                                    "finalCorrs" => finalCorrs, "predErr" => predErr))
