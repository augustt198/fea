using GeometryTypes
using Gadfly

const TRI_NEIGHBOR_A = 1
const TRI_NEIGHBOR_B = 2
const TRI_NEIGHBOR_C = 3

mutable struct DelaunayTriangle{T<:Point2} 
    # points in ccw order
    a::T; b::T; c::T
    # neighbor indices. 0 = no neighbor
    na::Int64; nb::Int64; nc::Int64
    active::Bool
end

mutable struct DelaunayTess2D{T<:Point2}
    faces::AbstractVector{DelaunayTriangle{T}}
    _last_tri_idx::Int
end

function delaunay2D(V::AbstractVector{T}) where T <: Point2
    minc, maxc = _coordinate_bounds(V)
    minc -= 1
    maxc += 1
    # base triangulation
    #     c   b
    # c + +---+
    #   |\ \  |
    #   | \ \ |
    #   |  \ \|
    #   +---+ +
    #   a   b a
    base_t1 = DelaunayTriangle(Point2(minc, minc), Point2(maxc, minc), Point2(minc, maxc), 0, 2, 0, true)
    base_t2 = DelaunayTriangle(Point2(maxc, minc), Point2(maxc, maxc), Point2(minc, maxc), 1, 0, 0, true)
    faces = [base_t1, base_t2]

    tess = DelaunayTess2D(faces, 1)

    for v_idx in 1:length(V)
        vert = V[v_idx]
        
        i = _find_tri_idx(tess, vert)
        t = tess.faces[i]
        #        a
        #       /|\
        #      / | \
        #     /  |  \
        #    / 1 x 2 \
        #   /  /   \  \
        #  / /   3   \ \
        # //           \\
        # +-------------+
        # b             c
        base_idx = length(tess.faces) + 1
        t₁ = DelaunayTriangle(t.a, t.b, vert, base_idx+2, base_idx+1, t.nc, true) # tri @ base_idx+0
        t₂ = DelaunayTriangle(t.c, t.a, vert, base_idx+0, base_idx+2, t.nb, true) # tri @ base_idx+1
        t₃ = DelaunayTriangle(t.b, t.c, vert, base_idx+1, base_idx+0, t.na, true) # tri @ base_idx+2
        push!(tess.faces, t₁)
        push!(tess.faces, t₂)
        push!(tess.faces, t₃)
        t.active = false

        _flip(tess)
        #break
    end
    return tess
end

# Lawson flip algorithm
function _flip(tess::DelaunayTess2D{T}) where T <: Point2
    # queue of edges. edge is represented by a pair a triangle and the neighbor that
    # forms the common edge
    queue::Vector{Tuple{DelaunayTriangle{T}, Int64}} = []
    # init queue
    for idx_offset in 0:2
        t = tess.faces[end - idx_offset]
        queue.push!((t, TRI_NEIGHBOR_C))
    end 

    while length(queue) > 0
        t, nbr = pop!(queue)
        if nbr == TRI_NEIGHBOR_A
            n = tess.faces[t.a]
            c1 = t.b ; c2 = t.c
        elseif nbr == TRI_NEIGHBOR_B
            n = tess.faces[t.b]
            c1 = t.a ; c2 = t.c
        else # nbr == TRI_NEIGHBOR_C
            n = tess.faces[t.c]
            c1 = t.a ; c2 = t.b
        end

        # find the vertex v in n not shared with t
        if t.a != c1 && t.a != c2
            v = t.a
        elseif t.b != c1 && t.b != c2
            v = t.b
        else
            v = t.c
        end

        if !isquadconvex(t.a, t.b, v, t.c)
            continue # nothing to do
        end

        _flip_tri(tess, t, n, v)
    end
end

function _flip_tri(tess::DelaunayTess2D{T}, t::DelaunayTriangle{T},
    n::DelaunayTriangle{T}, v::T) where T <: Point2
    # flip diagram:
    #        +(X)             +
    #       /|\              / \
    #      / | \     -->    /   \
    #     /  |  \          / ta  \
    #  u + t | n + v      +-------+
    #     \  |  /          \ tb  /
    #      \ | /            \   /
    #       \|/              \ /
    #        +(Y)             +

    X = n.a
    nX_a = n.na
    Y = n.b
    nY_a = n.nb
    if v == n.a
        X = n.b
        nX_a = n.nb
        Y = n.c
        nY_a = n.nc
    elseif v == n.b
        X = n.c
        nX_a = n.nc
        Y = n.a
        nY_a = n.na
    end

    u = t.c
    if t.a != X && t.a != Y
        u = t.a
        nX_b 
    elseif t.b != X && t.b != Y
        u = t.b
    end


    #ta = DelaunayTriangle()
end

function _find_tri_idx(tess::DelaunayTess2D{T}, pt::T) where T <: Point2
    for j in 1:length(tess.faces)
        tri = tess.faces[j]
        if !tri.active
            continue
        end
        if intriangle(tri, pt)
            return j
        end
    end
    
    return 0
end

function _coordinate_bounds(V::AbstractVector{T}) where T <: Point2
    mincoord = Inf; maxcoord = -Inf
    println("A")
    for v in V
        (v[1] > maxcoord) && (maxcoord = v[1])
        (v[1] < mincoord) && (mincoord = v[1])
        (v[2] > maxcoord) && (maxcoord = v[2])
        (v[2] < mincoord) && (mincoord = v[2])
    end
    println("C")
    return (mincoord, maxcoord)
end

function intriangle(tri::DelaunayTriangle{T}, P::T) where T <: Point2
    a = tri.c-tri.b; b = tri.a-tri.c; c = tri.b-tri.a
    ap = P-tri.a; bp = P-tri.b; cp = P-tri.c
    a_bp = a[1]*bp[2] - a[2]*bp[1];
    c_ap = c[1]*ap[2] - c[2]*ap[1];
    b_cp = b[1]*cp[2] - b[2]*cp[1];
    t0 = 0.0
    return ((a_bp >= t0) && (b_cp >= t0) && (c_ap >= t0))
end

# assuming A, B, C are ccw
function incircumcircle(tri::DelaunayTriangle{T}, P::T) where T <: Point2
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

# a, b, c, d are CCW vertices of simple quadrilateral
# test method: http://mathworld.wolfram.com/ConvexPolygon.html
function isquadconvex(a::T, b::T, c::T, d::T) where T<: Point2
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

function plottess(tess::DelaunayTess2D{T}) where T <: Point2
    x0 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    y0 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    x1 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    y1 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    println("\$", length(tess.faces))
    println("@", length(x0))
    for i in 1:length(tess.faces)
        t = tess.faces[i]
        if !t.active
            continue
        end
        j = i - 1
        # segment AB
        println(t.a[1])
        x0[j*3+1] = t.a[1]
        y0[j*3+1] = t.a[2]
        x1[j*3+1] = t.b[1]
        y1[j*3+1] = t.b[2]
        # segment BC
        x0[j*3+2] = t.b[1]
        y0[j*3+2] = t.b[2]
        x1[j*3+2] = t.c[1]
        y1[j*3+2] = t.c[2]
        # segment CA
        x0[j*3+3] = t.c[1]
        y0[j*3+3] = t.c[2]
        x1[j*3+3] = t.a[1]
        y1[j*3+3] = t.a[2]
    end
    println("XXX")
    #draw(SVG(24cm, 18cm), plot(x=x0, y=y0, xend=x1, yend=y1, Geom.segment))
    plot(x=x0, y=y0, xend=x1, yend=y1, Geom.segment)
end

# testing
base = [Point2f0(1, 1), Point2f0(2, 2), Point2f0(0.5, 1.25)]
println(base)
tess = delaunay2D(base)
plottess(tess)