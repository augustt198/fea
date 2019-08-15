# idk the actual way to do this sorry
push!(LOAD_PATH, string(@__DIR__, "/../src/"))

using GeometryTypes
using FEA

o = 10
tri = ACWTriangle(
    Point2f0(0.0+o, 1.0+o),
    Point2f0(0.0+o, 0.0+o),
    Point2f0(1.0+o, 0.0+o)
)

a = Point2f0(Float32[0.0, 0.0])
b = Point2f0(Float32[0.2, 0.0])
c = Point2f0(Float32[0.189163, 0.0649399])
#c = Point2f0(Float32[0.1, 0.2])
#f = (x) -> (barycentric(tri, Point2f0(x))[1] + barycentric(tri, Point2f0(x))[2] + barycentric(tri, Point2f0(x))[3])
f = (x) -> 1.0
curr_id = 1
load_integrand = (x) -> begin
    #FEA.barycentric(a, b, c, Point2f0(x))[3]
    detT = (b[2] - c[2])*(a[1] - c[1]) + (c[1] - b[1])*(a[2] - c[2])
    w1 = ( (b[2] - c[2])*(x[1] - c[1]) + (c[1] - b[1])*(x[2] - c[2]) ) / detT
    w2 = ( (c[2] - a[2])*(x[1] - c[1]) + (a[1] - c[1])*(x[2] - c[2]) ) / detT
    w3 = 1.0 - w1 - w2

    _w = FEA.barycentric2(a, b, c, x)[1]
    println(w1, " / ", _w)
    println(typeof(w1), " - ", typeof(_w))
    w1
end

println("Integrate:")
#println(FEA.integratetri(a, b, c, load_integrand))

_pa, _pb, _pc = Point2f0(a), Point2f0(b), Point2f0(c)
#load_contrib, err = FEA.integratetri(a, b, c, (x) -> begin
#    f(x) * FEA.barycentric(_pa, _pb, _pc, Point2f0(x))[1]
#end)

l, e = FEA.integratetri(a, b, c, load_integrand)

