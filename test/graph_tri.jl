# idk the actual way to do this sorry
push!(LOAD_PATH, string(@__DIR__, "/../src/"))

using GeometryTypes
using FEA

using Makie

function plottess(tess::DelaunayTess2D{T}) where T <: Point2
    vertices_x     = map(v -> v[1], tess.verts)
    vertices_y     = map(v -> v[2], tess.verts)
    vertices       = hcat(vertices_x, vertices_y)
    faces_active   = filter(t -> t.active, tess.faces)
    connectivity_a = map(t -> t.a, faces_active)
    connectivity_b = map(t -> t.b, faces_active)
    connectivity_c = map(t -> t.c, faces_active)
    connectivity   = hcat(connectivity_a, connectivity_b, connectivity_c)
    color          = sin.(collect(1:length(tess.verts)) / 50.0f0)

    # hide extreme points
    vertices[1:3, :]  .= vertices[4:4, :]

    scene = Scene()

    # enable to plot circumcircles
    if false
        circle_fn = (t) -> begin
            a = vertices[t[1], :]
            b = vertices[t[2], :]
            c = vertices[t[3], :]
            o, R = FEA.circumcenterwithradius(Point2f0(a), Point2f0(b), Point2f0(c))
            Circle(Point2f0(o), R)
        end
        circles = map(circle_fn, eachrow(connectivity))
        poly!(scene, circles, color=(:white, 0.0), strokecolor=:red, strokewidth=3)
    end

    scatter!(scene, vertices, color=:black, markersize=0.025)
    mesh!(scene, vertices, connectivity, color=color, shading=false, colormap=:plasma)
    cam = cam2d!(scene, panbutton=Mouse.left)
    wireframe!(scene[end][1], color=(:black, 0.6), linewidth=2)
end

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
