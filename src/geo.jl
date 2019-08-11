using GeometryTypes
using LinearAlgebra

const TRI_NEIGHBOR_A = 1
const TRI_NEIGHBOR_B = 2
const TRI_NEIGHBOR_C = 3

abstract type AbstractTriangle{T} end
# vertices stored clockwise
abstract type AbstractCWTriangle{T} <: AbstractTriangle{T} end
# vertices stored anticlockwise
abstract type AbstractACWTriangle{T} <: AbstractTriangle{T} end

struct ACWTriangle{T <: Point2} <: AbstractACWTriangle{T}
    a::T; b::T; c::T
end

# returns -1 if T not in tri
# returns 0 if T within tri
# returns TRI_NEIGHBOR_A/B/C if T is on edge A/B/C
function intriangle(tri::AbstractACWTriangle{T}, P::T) where T <: Point2
    a    = tri.c-tri.b; b = tri.a-tri.c; c = tri.b-tri.a
    ap   = P-tri.a; bp = P-tri.b; cp = P-tri.c
    a_bp = a[1]*bp[2] - a[2]*bp[1];
    c_ap = c[1]*ap[2] - c[2]*ap[1];
    b_cp = b[1]*cp[2] - b[2]*cp[1];

    t0        = zero(Float32)
    epsilon   = eps(Float32)
    a_bp_zero = isapprox(a_bp, t0, atol=epsilon)
    c_ap_zero = isapprox(c_ap, t0, atol=epsilon)
    b_cp_zero = isapprox(b_cp, t0, atol=epsilon)

    if (a_bp >= t0 || a_bp_zero) && (b_cp >= t0 || b_cp_zero) && (c_ap >= t0 || c_ap_zero)
        if a_bp_zero
            return TRI_NEIGHBOR_A
        elseif b_cp_zero
            return TRI_NEIGHBOR_B
        elseif c_ap_zero
            return TRI_NEIGHBOR_C
        else
            return 0
        end
    else
        return -1
    end
end

# assuming A, B, C are ccw
function incircumcircle(tri::AbstractACWTriangle{T}, P::T) where T<: Point2
    A = tri.a ; B = tri.b ; C = tri.c
    M11 = A[1] - P[1]
    M21 = B[1] - P[1]
    M31 = C[1] - P[1]
    M12 = A[2] - P[2]
    M22 = B[2] - P[2]
    M32 = C[2] - P[2]
    M13 = M11*M11 + M12*M12
    M23 = M21*M21 + M22*M22
    M33 = M31*M31 + M32*M32
    det = M11*(M22*M33 - M23*M32) - M12*(M21*M33 - M23*M31) + M13*(M21*M32 - M22*M31)
    return det > 0
end

# a, b, c, d are vertices (ACW order) of simple quadrilateral
# test method: http://mathworld.wolfram.com/ConvexPolygon.html
function isquadconvex(a::T, b::T, c::T, d::T) where T <: Point2
    vec_ab = b - a ; vec_ab_P = Point2(-vec_ab[2], vec_ab[1])
    vec_bc = c - b ; vec_bc_P = Point2(-vec_bc[2], vec_bc[1])
    vec_cd = d - c ; vec_cd_P = Point2(-vec_cd[2], vec_cd[1])
    vec_da = a - d ; vec_da_P = Point2(-vec_da[2], vec_da[1])

    sgn = sign(sum(vec_ab_P .*  vec_bc))

    if sign(sum(vec_bc_P .* vec_cd)) != sgn
        return false
    elseif sign(sum(vec_cd_P .* vec_da)) != sgn
        return false
    elseif sign(sum(vec_da_P .* vec_ab)) != sgn
        return false
    else
        return true
    end
end

# transform point P to barymetric coords by solving:
#   W_1*a_x + W_2*b_x + W_3*c_x = P_x
#   W_1*a_y + W_2*b_y + W_3*c_y = P_y
#   W_1     + W_2     + W_3     = 1
function barycentric(tri::AbstractTriangle{T}, P::T) where T <: Point2
    A = [
        tri.a[1] tri.b[1] tri.c[1]
        tri.a[2] tri.b[2] tri.c[2]
        1.0 1.0 1.0
    ]
    b = [P[1] ; P[2] ; 1.0]
    return A \ b
end

# solves
# a1_x + (a2_x - a1_x)*t = b1_x + (b2_x - b1_x)*s
# a1_y + (a2_y - a1_y)*t = b1_y + (b2_y - b1_y)*s
function lineintersection(a1::T, a2::T, b1::T, b2::T) where T <: Point2
    A = [
        (a2[1] - a1[1]) -(b2[1] - b1[1])
        (a2[2] - a1[2]) -(b2[2] - b1[2])
    ]
    if det(A) ≈ 0
        return Inf, Inf, T(0,0)
    end
    b = [
        b1[1] - a1[1]
        b1[2] - a1[2]
    ]
    t, s = A \ b
    # intersection point
    p = a1 + (a2 - a1)*t
    return t, s, p
end
