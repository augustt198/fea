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

f = (x) -> x[1]

println("Integrate:")
println( integratetri(tri, f) )
