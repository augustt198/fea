using GeometryTypes
using Gadfly

const TRI_NEIGHBOR_A = 1
const TRI_NEIGHBOR_B = 2
const TRI_NEIGHBOR_C = 3

# TODO: use hierarchical structure for faster triangle lookup
mutable struct DelaunayTriangle{T<:Point2} 
    # points in ccw order
    a::T; b::T; c::T
    # neighbor indices. 0 = no neighbor
    na::Int64; nb::Int64; nc::Int64
    active::Bool
    child1::Int64; child2::Int64; child3::Int64
end

DelaunayTriangle(a, b, c, na, nb, nc, active) = DelaunayTriangle(a, b, c, na, nb, nc, active, 0, 0, 0)
DelaunayTriangle(a, b, c, na, nb, nc) = DelaunayTriangle(a, b, c, na, nb, nc, false)
DelaunayTriangle(a, b, c) = DelaunayTriangle(a, b, c, 0, 0, 0)

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
    base_t1 = DelaunayTriangle(Point2(minc, minc), Point2(maxc, minc), Point2(minc, maxc), 2, 0, 0, true)
    base_t2 = DelaunayTriangle(Point2(maxc, minc), Point2(maxc, maxc), Point2(minc, maxc), 0, 1, 0, true)
    base_tri = DelaunayTriangle(Point2f0(0, 5), Point2f0(-5, -5), Point2f0(5, -5), 0, 0, 0, true)
    faces = [base_tri] # [base_t1, base_t2]

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
        
        (t.nc > 0) && _update_neighbor(tess.faces[t.nc], i, base_idx + 0)
        (t.nb > 0) && _update_neighbor(tess.faces[t.nb], i, base_idx + 1)
        (t.na > 0) && _update_neighbor(tess.faces[t.na], i, base_idx + 2)

        println("=== CHECK NBR PRE")
        _check_nbr(tess)

        _flip(tess)

        println("=== CHECK NBR POST")
        _check_nbr(tess)
    end
    return tess
end

# update neighbor indices for surrounding triangles
function _update_neighbor(t::DelaunayTriangle{T}, old_idx::Int, new_idx::Int) where T <: Point2
    if t.na == old_idx
        t.na = new_idx
    elseif t.nb == old_idx
        t.nb = new_idx
    else
        t.nc = new_idx
    end
end

function _shares_common_edge(t::DelaunayTriangle{T}, n::DelaunayTriangle{T}) where T <: Point2
    length(unique([t.a, t.b, t.c, n.a, n.b, n.c])) == 4
end

function _check_nbr_single(tess::DelaunayTess2D{T}, t::DelaunayTriangle{T}) where T<: Point2
    if t.na > 0
        if !_shares_common_edge(t, tess.faces[t.na])
            println("********************** BAD A")
            println("**********************-> T ", t)
            println("**********************-> N ", tess.faces[t.na])
        end
    end
    if t.nb > 0
        if !_shares_common_edge(t, tess.faces[t.nb])
            println("********************** BAD B")
            println("**********************-> T ", t)
            println("**********************-> N ", tess.faces[t.nb])
        end
    end
    if t.nc > 0
        if !_shares_common_edge(t, tess.faces[t.nc])
            println("********************** BAD C")
            println("**********************-> T ", t)
            println("**********************-> N ", tess.faces[t.nc])
        end
    end
end

function _check_nbr(tess::DelaunayTess2D{T}) where T <: Point2
    for t in tess.faces
        if !t.active
            continue
        end
        _check_nbr_single(tess, t)
    end
end

# Lawson flip algorithm
function _flip(tess::DelaunayTess2D{T}) where T <: Point2
    # queue of edges. edge is represented by a pair a triangle and the neighbor that
    # forms the common edge
    queue::Vector{Tuple{DelaunayTriangle{T}, Int64, Int64}} = []
    # init queue
    for offset in 0:2
        idx = length(tess.faces) - offset
        t = tess.faces[idx]
        push!(queue, (t, idx, TRI_NEIGHBOR_C))
    end 

    while length(queue) > 0
        t, it, nbr = pop!(queue)
        if nbr == TRI_NEIGHBOR_A
            if t.na < 1
                continue
            end
            n = tess.faces[t.na]
            in = t.na
            u = t.a ; c1 = t.b ; c2 = t.c
        elseif nbr == TRI_NEIGHBOR_B
            if t.nb < 1
                continue
            end
            n = tess.faces[t.nb]
            in = t.nb
            u = t.b ; c1 = t.a ; c2 = t.c
        else # nbr == TRI_NEIGHBOR_C
            if t.nc < 1
                continue
            end
            n = tess.faces[t.nc]
            in = t.nc
            u = t.c ; c1 = t.a ; c2 = t.b
        end

        # find the vertex v in n not shared with t
        if n.a != c1 && n.a != c2
            v = n.a
        elseif n.b != c1 && n.b != c2
            v = n.b
        else
            v = n.c
        end

        if !isquadconvex(u, c1, v, c2)
            continue # nothing to do
        end

        if incircumcircle(t, v)
            _flip_tri(tess, t, n, it, in, v)
        end

        break
    end
end

function _flip_tri(tess::DelaunayTess2D{T}, t::DelaunayTriangle{T},
    n::DelaunayTriangle{T}, it::Int64, in::Int64, v::T) where T <: Point2
    # flip diagram:
    #        +(X)               +
    #       /|\                / \
    #      / | \     --> nY_b /   \  nY_a
    #     /  |  \            / ta  \
    #  u + t | n + v        +-------+
    #     \  |  /            \ tb  /
    #      \ | /        nX_b  \   /  nX_a
    #       \|/                \ /
    #        +(Y)               +

    X = n.a
    nX_a = n.na
    Y = n.b
    nY_a = n.nb
    t_idx = n.nc
    if v == n.a
        X = n.b
        nX_a = n.nb
        Y = n.c
        nY_a = n.nc
        t_idx = n.na
    elseif v == n.b
        X = n.c
        nX_a = n.nc
        Y = n.a
        nY_a = n.na
        t_idx = n.nb
    end

    u = t.c
    nX_b = t.nb
    nY_b = t.na
    n_idx = t.nc
    if t.a != X && t.a != Y
        u = t.a
        nX_b = t.nc
        nY_b = t.nb
        n_idx = t.na
    elseif t.b != X && t.b != Y
        u = t.b
        nX_b = t.na
        nY_b = t.nc
        n_idx = t.nb
    end

    # t becomes ta, n becomes tb
    t.a = v ; t.b = X ; t.c = u
    t.na = nY_b ; t.nb = n_idx ; t.nc = nY_a

    n.a = u ; n.b = Y ; n.c = v
    n.na = nX_a ; n.nb = t_idx ; n.nc = nX_b

    (t.nc > 0) && _update_neighbor(tess.faces[t.nc], in, it)
    (n.nc > 0) && _update_neighbor(tess.faces[n.nc], it, in)
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

function plottess(tess::DelaunayTess2D{T}) where T <: Point2
    x0 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    y0 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    x1 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    y1 = Vector{eltype(T)}(undef, length(tess.faces)*3)
    for i in 1:length(tess.faces)
        t = tess.faces[i]
        if !t.active
            continue
        end
        j = i - 1
        # segment AB
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
    #draw(SVG(24cm, 18cm), plot(x=x0, y=y0, xend=x1, yend=y1, Geom.segment))
    plot(x=x0, y=y0, xend=x1, yend=y1, Geom.segment)
end

# testing
#base = [Point2f0(1, 1.5), Point2f0(2, 2), Point2f0(0.5, 1.25)]
degs = LinRange(0, 2*π, 10)
base = map(x -> Point2f0(cos(x), sin(x)), degs)[1:10]
push!(base, Point2f0(0,0))
println(base)
tess = delaunay2D(base)
plottess(tess)
