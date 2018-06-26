using MAT
using ProgressMeter

function convToJulia(arr)
    if length(arr) == 0
        return Vector{Float64}()
    elseif length(arr) == 1
        return [arr]
    else
        return Vector{Float64}(arr[:,1])
    end
end

function loadMaps(weightMatrixFile::String, placeCentersFile::String="place_field_data.mat")
    predictedW = matread(weightMatrixFile);
    vars = matread(placeCentersFile)
    placeCenters = convToJulia.(vars["PF_centers"])[:,1]
    targetMap = TargetMap(placeCenters; placeWidth=50.0)
    W = predictedW["fittedW"]
    forwardMap = ForwardMap(targetMap.fPeak, targetMap.inhibThres, targetMap.wI, W, zeros(size(W, 1)))
    return targetMap, forwardMap
end