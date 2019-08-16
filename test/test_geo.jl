# idk the actual way to do this sorry
push!(LOAD_PATH, string(@__DIR__, "/../src/"))

using GeometryTypes
using FEA

pa = Point2f0(-0.1,0.1)
pb = Point2f0(0.5,0)
pc = Point2f0(0,1)
x = Point2f0(0, 0.25)
curr_id = 2

h = 0.00001
h_x, h_y = Float64[h ; 0], Float64[0 ; h]
w_0 = FEA.barycentric2(pa, pb, pc, x)[curr_id]
w_1 = FEA.barycentric2(pa, pb, pc, x+h_x)[curr_id]
w_2 = FEA.barycentric2(pa, pb, pc, x+h_y)[curr_id]

w_x = (w_1 - w_0) / h
w_y = (w_2 - w_0) / h

∇W = FEA.barycentric2_grad(pa, pb, pc, x)

println(">>> NUMERIC: ", (w_x, w_y))
println(">>> ANALYTICAL: ", ∇W[curr_id])
