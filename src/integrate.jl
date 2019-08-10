using Cubature
using LinearAlgebra

# integrate f over tri (anticlockwise points)
function integratetri(tri::AbstractACWTriangle{T}, f::Function) where T <: Point2
    u = tri.b - tri.a
    v = tri.c - tri.a
    # transform from i, j basis to u, v
    M = hcat(u, v)
    J = det(M)
    f′ = (x) -> begin
        # square to triangle transform
        x′ = [x[1] ; x[2] * (1.0 - x[1])]
        p = M*x′ + tri.a
        return f(p) * J * (1.0 - x[1])
    end

    xmin = [ 0 ; 0 ]
    xmax = [ 1 ; 1 ]
    return hcubature(f′, xmin, xmax)
end
