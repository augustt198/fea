# idk the actual way to do this sorry
push!(LOAD_PATH, string(@__DIR__, "/../src/"))

using GeometryTypes
using FEA

using Gadfly

function plottess(tess::DelaunayTess2D{T}) where T <: Point2
    x0 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    y0 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    x1 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    y1 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    for i in 1:length(tess.faces)
        t = tess.faces[i]
        j = i - 1
        d = t.active
        # segment AB
        x0[j*3+1] = d ? t.a[1] : 0
        y0[j*3+1] = d ? t.a[2] : 0
        x1[j*3+1] = d ? t.b[1] : 0
        y1[j*3+1] = d ? t.b[2] : 0
        # segment BC
        x0[j*3+2] = d ? t.b[1] : 0
        y0[j*3+2] = d ? t.b[2] : 0
        x1[j*3+2] = d ? t.c[1] : 0
        y1[j*3+2] = d ? t.c[2] : 0
        # segment CA
        x0[j*3+3] = d ? t.c[1] : 0
        y0[j*3+3] = d ? t.c[2] : 0
        x1[j*3+3] = d ? t.a[1] : 0
        y1[j*3+3] = d ? t.a[2] : 0
    end
    #draw(SVG(24cm, 18cm), plot(x=x0, y=y0, xend=x1, yend=y1, Geom.segment))
    plot(x=x0, y=y0, xend=x1, yend=y1, Geom.segment, Coord.cartesian(fixed=true))
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
