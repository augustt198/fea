using GeometryTypes

function intriangle(tri::DelaunayTriangle{T}, P::T) where T <: Point
    a = tri.c-tri.b; b = tri.a-tri.c; c = tri.b-tri.a
    ap = P-tri.a; bp = P-tri.b; cp = P-tri.c
    a_bp = a[1]*bp[2] - a[2]*bp[1];
    c_ap = c[1]*ap[2] - c[2]*ap[1];
    b_cp = b[1]*cp[2] - b[2]*cp[1];
    t0 = zero(eltype(T))
    return ((a_bp >= t0) && (b_cp >= t0) && (c_ap >= t0))
end

# assuming A, B, C are ccw
function incircumcircle(tri::DelaunayTriangle{T}, P::T) where T<: Point2
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