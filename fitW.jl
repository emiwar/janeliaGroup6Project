import MAT
include("megamap.jl")


placeWidth = 50.0
learningRate = 3e-7


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
W = computeW(targetMap, 50:500:3950, -learningRate; maxIter=20)

MAT.matwrite("fittedW_5m_grid.mat", Dict("fittedW" => W))
