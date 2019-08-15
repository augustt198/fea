using GeometryTypes
using FEA

using Makie

function plottess(tess::DelaunayTess2D{T}, draw_circumcircles::Bool=false, ncolors=[]) where T <: Point2
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

    scene = Scene(scale_plot=false)

    # enable to plot circumcircles
    if draw_circumcircles
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
    mesh!(scene, vertices, connectivity, color=ncolors, shading=false, colormap=:plasma)
    cam = cam2d!(scene, panbutton=Mouse.left)
    wireframe!(scene[end][1], color=(:black, 0.6), linewidth=2)
    return scene
end
