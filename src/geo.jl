using GeometryTypes
using LinearAlgebra

const TRI_NEIGHBOR_A = 1
const TRI_NEIGHBOR_B = 2
const TRI_NEIGHBOR_C = 3

const TRI_VERTEX_A = 1
const TRI_VERTEX_B = 2
const TRI_VERTEX_C = 3

abstract type AbstractTriangle{T} end
# vertices stored clockwise
abstract type AbstractCWTriangle{T} <: AbstractTriangle{T} end
# vertices stored anticlockwise
abstract type AbstractACWTriangle{T} <: AbstractTriangle{T} end

struct CWTriangle{T <: Point2} <: AbstractCWTriangle{T}
    a::T; b::T; c::T
end

struct ACWTriangle{T <: Point2} <: AbstractACWTriangle{T}
    a::T; b::T; c::T
end

struct IndexedLineSegment
    a::Int64
    b::Int64
    # which boundary this segment belongs to.
    boundary::Int64
end

IndexedLineSegment(a, b) = IndexedLineSegment(a, b, 1)

function intriangle(tri::AbstractACWTriangle{T}, P::T) where T <: Point2
    return intriangle(tri.a, tri.b, tri.c, P)
end

# returns -1 if T not in tri
# returns 0 if T within tri
# returns TRI_NEIGHBOR_A/B/C if T is on edge A/B/C
function intriangle(A::T, B::T, C::T, P::T) where T <: Point2
    a    = C-B; b = A-C; c = B-A
    ap   = P-A; bp = P-B; cp = P-C
    a_bp = a[1]*bp[2] - a[2]*bp[1];
    b_cp = b[1]*cp[2] - b[2]*cp[1];
    c_ap = c[1]*ap[2] - c[2]*ap[1];

    # TODO don't hardcode Float32
    t0        = zero(Float32)
    epsilon   = 1e-5 # eps(Float32)
    a_bp_zero = isapprox(a_bp, t0, atol=epsilon)
    b_cp_zero = isapprox(b_cp, t0, atol=epsilon)
    c_ap_zero = isapprox(c_ap, t0, atol=epsilon)

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
        if a_bp < 0
            return -TRI_NEIGHBOR_A
        elseif b_cp < 0
            return -TRI_NEIGHBOR_B
        else
            return -TRI_NEIGHBOR_C
        end
    end
end

function incircumcircle(tri::AbstractACWTriangle{T}, P::T) where T <: Point2
    return incircumcircle(tri.a, tri.b, tri.c, P)
end

# assuming A, B, C are ccw
function incircumcircle(A::T, B::T, C::T, P::T) where T <: Point2
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

    sgn = sign(vec_ab_P' *  vec_bc)

    if sign(vec_bc_P' * vec_cd) != sgn
        return false
    elseif sign(vec_cd_P' * vec_da) != sgn
        return false
    elseif sign(vec_da_P' * vec_ab) != sgn
        return false
    else
        return true
    end
end

function barycentric(tri::AbstractTriangle{T}, P::T) where T <: Point2
    return barycentric(tri.a, tri.b, tri.c, P)
end

# transform point P to barymetric coords by solving:
#   W_1*a_x + W_2*b_x + W_3*c_x = P_x
#   W_1*a_y + W_2*b_y + W_3*c_y = P_y
#   W_1     + W_2     + W_3     = 1
function barycentric(A::T, B::T, C::T, P::T) where T <: Point2
    M = [
        A[1] B[1] C[1]
        A[2] B[2] C[2]
        1.0 1.0 1.0
    ]
    b = [P[1] ; P[2] ; 1.0]
    return M \ b
end

function barycentric2(a, b, c, x)
    detT = (b[2] - c[2])*(a[1] - c[1]) + (c[1] - b[1])*(a[2] - c[2])
    w₁ = ( (b[2] - c[2])*(x[1] - c[1]) + (c[1] - b[1])*(x[2] - c[2]) ) / detT
    w₂ = ( (c[2] - a[2])*(x[1] - c[1]) + (a[1] - c[1])*(x[2] - c[2]) ) / detT
    w₃ = 1.0 - w₁ - w₂
    return w₁, w₂, w₃
end

# http://blog.marshalljiang.com/gradient-of-the-barycentric-coordinate-in-2d/
function barycentricgrad(a, b, c, x)
    det_T = Float64( (b[2] - c[2])*(a[1] - c[1]) + (c[1] - b[1])*(a[2] - c[2]) )
    bc = c - b
    ca = a - c
    ab = b - a
    L_bc = sqrt(bc' * bc) # side lengths
    L_ca = sqrt(ca' * ca)
    L_ab = sqrt(ab' * ab)
    bc_n = normalize([bc[2] ; -bc[1]]) # normals
    ca_n = normalize([ca[2] ; -ca[1]])
    ab_n = normalize([ab[2] ; -ab[1]])

    ∇w₁ = -L_bc/(det_T) * bc_n
    ∇w₂ = -L_ca/(det_T) * ca_n
    ∇w₃ = -L_ab/(det_T) * ab_n

    return ∇w₁, ∇w₂, ∇w₃
end

# solves
# a1_x + (a2_x - a1_x)*t = b1_x + (b2_x - b1_x)*s
# a1_y + (a2_y - a1_y)*t = b1_y + (b2_y - b1_y)*s
function lineintersection(a1::T, a2::T, b1::T, b2::T) where T <: Point2
    A = [
        (a2[1] - a1[1]) -(b2[1] - b1[1])
        (a2[2] - a1[2]) -(b2[2] - b1[2])
    ]
    ε = eps(eltype(T))
    if isapprox(det(A), 0, atol=ε)
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

# TODO change
function findregion(V::AbstractVector{T}, dR::AbstractVector{IndexedLineSegment},
    pt::T, off=3, nregions=0) where T <: Point2

    pt_proj = pt + T(1, 0)

    crossings = Vector{Int64}(undef, nregions)
    crossings .= 0
    for seg in dR
        # TODO make indexing more elegant
        pa, pb = V[seg.a+off], V[seg.b+off]
        l1, l2, _ = lineintersection(pa, pb, pt, pt_proj)

        # need to be careful about the projected line
        # crossing through a vertex in dR
        if 0 <= l1 < 1 && l2 >= 0
            if seg.boundary > length(crossings)
                s = length(crossings)
                resize!(crossings, seg.boundary)
                crossings[s+1:end] .= 0
            end
            crossings[seg.boundary] += 1
        end
    end

    for b in length(crossings):-1:1
        c = crossings[b]
        if c > 0 && isodd(c)
            return b
        end
    end
    return 0
end

# TODO change
# determines if pt is inside region with boundary dR by counting crossings
function regiontestinside(V::AbstractVector{T}, dR::AbstractVector{IndexedLineSegment},
    pt::T) where T <: Point2

    return isodd(regioncrossings(V, dR, pt))
end

#       + a
#      / \
#   d +   + e
#    /     \
# b +-------+ c
function circumcenterwithradius(a::T, b::T, c::T) where T <: Point2
    d = (a + b) / 2
    e = (a + c) / 2
    vec_ab = b - a
    vec_ac = c - a
    vec_d_normal = [-vec_ab[2] ; vec_ab[1]]
    vec_e_normal = [vec_ac[2] ; -vec_ac[1]]

    l1, l2, o = lineintersection(d, d + vec_d_normal, e, e + vec_e_normal)
    ao = o - a
    R = sqrt(ao' * ao)
    return o, R
end
