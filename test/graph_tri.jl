# idk the actual way to do this sorry
push!(LOAD_PATH, string(@__DIR__, "/../src/"))

using GeometryTypes
using FEA

include("plottess.jl")

# testing
divs = 100
degs = LinRange(0, 2*π, divs)
base1 = map(x -> Point2f0(1.0*cos(x), 1.0*sin(x)), degs)[1:divs-1]
base2 = map(x -> Point2f0(1.5*cos(x), 1.5*sin(x)), degs)[1:divs-1]
base3 = map(x -> Point2f0(1.25*cos(x), 1.25*sin(x)), degs)[1:divs-1]
base4 = map(x -> Point2f0(rand(Float32)*2-1, rand(Float32)*2-1), degs)[1:divs-1]
base = vcat(base1, base2, base3, base4)
tess = delaunay2D(base)
#plottess(tess)

idx = length(base)+1
push!(base, Point2f0(0.5, 0.5))
push!(base, Point2f0(-0.5, 0.5))
push!(base, Point2f0(-0.5, -0.5))
push!(base, Point2f0(0.5, -0.5))

pslg = FEA.PSLG([])
push!(pslg.segments, IndexedLineSegment(idx+0, idx+1))
push!(pslg.segments, IndexedLineSegment(idx+1, idx+2))
push!(pslg.segments, IndexedLineSegment(idx+2, idx+3))
push!(pslg.segments, IndexedLineSegment(idx+3, idx+0))

tess = conformingDelaunay2D(base, pslg)
plottess(tess)
