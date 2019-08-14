using Cubature
using LinearAlgebra

function integratetri(tri::AbstractACWTriangle{T}, f::Function) where T <: Point2
    integratetri(tri.a, tri.b, tri.c, f)
end

# integrate f over tri (anticlockwise points)
function integratetri(a::T, b::T, c::T, f::Function) where T <: Point2
    u = b - a
    v = c - a
    # transform from i, j basis to u, v
    M = hcat(u, v)
    J = det(M)
    f′ = (x) -> begin
        # square to triangle transform
        x′ = [x[1] ; x[2] * (1.0 - x[1])]
        p = M*x′ + a
        return f(p) * J * (1.0 - x[1])
    end

    xmin = [ 0 ; 0 ]
    xmax = [ 1 ; 1 ]
    return hcubature(f′, xmin, xmax)
end
