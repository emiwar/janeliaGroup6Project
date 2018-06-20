mutable struct TargetMap
    placeWidth::Float64
    u0::Float64 #TODO: Find a better name...
    fPeak::Float64
    inhibThres::Float64
    IPeak::Float64
    wI::Float64
    placeCenters::Vector{Vector{Float64}}
end

TargetMap(nCells::Int64, placeWidth, u0, fPeak, inhibThres, IPeak, wI) = TargetMap(placeWidth, u0, fPeak, inhibThres, IPeak, wI, [Vector{Float64}() for i=1:nCells])


function TargetMap(nCells::Int64, allPlaceCenters::Vector{Float64}; placeWidth=5.0, u0=.2, fPeak=15.0, IPeak=0.3)
    targetMap = TargetMap(nCells, placeWidth, u0, fPeak, 0.0, IPeak, 0.0)
    for c in allPlaceCenters
        cell = rand(1:nCells)
        push!(targetMap.placeCenters[cell], c)
    end
    targetMap.inhibThres = mean(0.9*sum(fTarget(c, targetMap)) for c in allPlaceCenters)
    targetMap.wI = targetMap.u0 / (targetMap.inhibThres*(1/0.9 - 1))
    return targetMap
end

nCells(targetMap::TargetMap) = length(targetMap.placeCenters)

function uTune(d::Float64, targetMap::TargetMap)
    targetMap.fPeak * ((1+targetMap.u0) * exp(-d^2/(2*targetMap.placeWidth^2)) - targetMap.u0)
end

function fTarget(x::Number, targetMap::TargetMap)
    result = zeros(nCells(targetMap))
    for cell = 1:nCells(targetMap)
        for placeCenter in targetMap.placeCenters[cell]
            result[cell] += max(0, uTune(placeCenter - x, targetMap))
        end
    end
    return result
end

function input(x::Number, targetMap::TargetMap)
    result = zeros(nCells(targetMap))
    for cell = 1:nCells(targetMap)
        for placeCenter in targetMap.placeCenters[cell]
            result[cell] += exp(-(placeCenter - x)^2 / (2*targetMap.placeWidth^2))
        end
    end
    return result * targetMap.IPeak
end

function fInhibition(f::Vector, targetMap::TargetMap)
    return targetMap.wI * max(0, sum(f) - targetMap.inhibThres)
end

function fProjection(x::Number, fBar::Vector, W::Matrix, targetMap::TargetMap)
    #targetMap.fPeak * max.(0, W*fBar - fInhibition(fBar, targetMap) + input(x, targetMap))
    targetMap.fPeak * max.(0, W*fBar - targetMap.u0 + input(x, targetMap))
end


function computeW(targetMap::TargetMap, learningRegion, s::Number, tol::Number)
    N = nCells(targetMap)
    W = zeros(N,N)
    maxIter = 10000
    i = 0
    while i<maxIter
        i+=1
        deltaW = zeros(N, N)
        for x in learningRegion
            fBar = fTarget(x, targetMap)
            fProj = fProjection(x, fBar, W, targetMap)
            for j=1:N, k=1:N
                if j==k
                    continue
                end
                deltaW[j, k] += (fProj[j] - fBar[j])*fBar[k]
            end
        end
        W += s*deltaW
        println(sum(deltaW.^2))
        #break
        if sum(deltaW.^2)<tol
            break
        end
    end
    return W
end

