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
        if !t.active
            continue
        end
        j = i - 1
        # segment AB
        x0[j*3+1] = t.a[1]
        y0[j*3+1] = t.a[2]
        x1[j*3+1] = t.b[1]
        y1[j*3+1] = t.b[2]
        # segment BC
        x0[j*3+2] = t.b[1]
        y0[j*3+2] = t.b[2]
        x1[j*3+2] = t.c[1]
        y1[j*3+2] = t.c[2]
        # segment CA
        x0[j*3+3] = t.c[1]
        y0[j*3+3] = t.c[2]
        x1[j*3+3] = t.a[1]
        y1[j*3+3] = t.a[2]
    end
    #draw(SVG(24cm, 18cm), plot(x=x0, y=y0, xend=x1, yend=y1, Geom.segment))
    plot(x=x0, y=y0, xend=x1, yend=y1, Geom.segment, Coord.cartesian(fixed=true))
end

# testing
#base = [Point2f0(1, 1.5), Point2f0(2, 2), Point2f0(0.5, 1.25)]
degs = LinRange(0, 2*π, 20)
base1 = map(x -> Point2f0(cos(x), sin(x)), degs)[1:15]
base2 = map(x -> Point2f0(2*cos(x), 2*sin(x)), degs)[1:15]
#base3 = map(x -> Point2f0(rand(Float32)*2-1, rand(Float32)*2-1), degs)[1:15]
#base3 = [Point2f0(-0.5, -0.5), Point2f0(0.5, -0.5), Point2f0(0.5, 0.5), Point2f0(-0.5, 0.5), Point2f0(0.1, -0.25)]
base = vcat(base1, base2)
#push!(base, Point2f0(0,0))
println(base)
tess = delaunay2D(base)
plottess(tess)

