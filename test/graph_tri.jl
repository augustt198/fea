# idk the actual way to do this sorry
push!(LOAD_PATH, string(@__DIR__, "/../src/"))

using GeometryTypes
using FEA

using Makie

function plottess(tess::DelaunayTess2D{T}) where T <: Point2
    vertices = [0.0f0 0.0f0]
    connectivity = [1 1 1]
    for t in tess.faces
        if t.active
            idx = size(vertices)[1]
            vertices = vcat(vertices, [t.a[1] t.a[2]])
            vertices = vcat(vertices, [t.b[1] t.b[2]])
            vertices = vcat(vertices, [t.c[1] t.c[2]])
            connectivity = vcat(connectivity, [idx+1 idx+2 idx+3])
        end
    end

    scene = mesh(vertices, connectivity, color=:white, shading=false)
    wireframe!(scene[end][1], color=(:black, 0.6), linewidth=2)
end

# testing
divs = 25
degs = LinRange(0, 2*π, divs)
base1 = map(x -> Point2f0(cos(x), sin(x)), degs)[1:divs-1]
base2 = map(x -> Point2f0(2*cos(x), 2*sin(x)), degs)[1:divs-1]
base3 = map(x -> Point2f0(1.5*cos(x), 1.5*sin(x)), degs)[1:divs-1]
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
