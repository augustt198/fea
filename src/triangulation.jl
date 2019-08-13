using GeometryTypes

const MAX_COORD = 10

# TODO: use hierarchical structure for faster triangle lookup
mutable struct DelaunayTriangle
    # point indices in ccw order
    a::Int64; b::Int64; c::Int64
    # neighbor indices. 0 = no neighbor
    na::Int64; nb::Int64; nc::Int64
    active::Bool
end

DelaunayTriangle(a, b, c, na, nb, nc) = DelaunayTriangle(a, b, c, na, nb, nc, false)
DelaunayTriangle(a, b, c) = DelaunayTriangle(a, b, c, 0, 0, 0)

struct TriEdge
    tri_idx::Int64
    # one of TRI_NEIGHBOR_A/B/C
    nbr::Int64
end

mutable struct DelaunayTess2D{T<:Point2}
    verts::AbstractVector{T}
    faces::AbstractVector{DelaunayTriangle}
end

struct PSLG
    segments::AbstractVector{IndexedLineSegment}
end

function delaunay2D(V::AbstractVector{T}) where T <: Point2
    maxc = _coordinate_bound(V)

    # base triangle large enough to enclose all points in V
    v₁ = T(0,        4*maxc)
    v₂ = T(-4*maxc, -4*maxc)
    v₃ = T(4*maxc,   -4*maxc)
    V = vcat([v₁, v₂, v₃], V)
    base_tri = DelaunayTriangle(1, 2, 3, 0, 0, 0, true)
    println(">>>>>>>>>> ", typeof(V))
    tess = DelaunayTess2D(V, [base_tri])

    for i in 4:length(tess.verts)
        _insert_point(tess, i)
    end

    _deactivate_extremal_triangles(tess)
    _check_nbr(tess)
 
    return tess
end

function _insert_point(tess::DelaunayTess2D{T}, vidx::Int64) where T <: Point2
    vert = tess.verts[vidx]
    i1, i2, ret1, ret2 = _find_tri_idx(tess, vert)

    if ret1 == 0 # interior
        i = i1
        t = tess.faces[i]
        # division diagram of t -> t1, t2, t3:
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
        t₁ = DelaunayTriangle(t.a, t.b, vidx, base_idx+2, base_idx+1, t.nc, true) # tri @ base_idx+0
        t₂ = DelaunayTriangle(t.c, t.a, vidx, base_idx+0, base_idx+2, t.nb, true) # tri @ base_idx+1
        t₃ = DelaunayTriangle(t.b, t.c, vidx, base_idx+1, base_idx+0, t.na, true) # tri @ base_idx+2
        push!(tess.faces, t₁)
        push!(tess.faces, t₂)
        push!(tess.faces, t₃)
        t.active = false
        
        (t.nc > 0) && _update_neighbor(tess.faces[t.nc], i, base_idx + 0)
        (t.nb > 0) && _update_neighbor(tess.faces[t.nb], i, base_idx + 1)
        (t.na > 0) && _update_neighbor(tess.faces[t.na], i, base_idx + 2)

        edges = [
            TriEdge(base_idx+0, TRI_NEIGHBOR_C),
            TriEdge(base_idx+1, TRI_NEIGHBOR_C),
            TriEdge(base_idx+2, TRI_NEIGHBOR_C)
        ]
        _flip(tess, edges)
        
        if !_check_nbr(tess)
            #break
        end 
    elseif ret1 > 0 # on edge
        #        + s
        #       /|\
        #      / | \
        #     /  |  \
        #    / 1 | 2 \
        #   +----+----+
        #  q    vert   r
        edges = Vector{TriEdge}()
        for (i, ret) in zip([i1, i2], [ret1, ret2])
            base_idx = length(tess.faces) + 1
            t = tess.faces[i]
            if ret == TRI_NEIGHBOR_A
                s, q, r = t.a, t.b, t.c
                ns, nq, nr = t.na, t.nb, t.nc
            elseif ret == TRI_NEIGHBOR_B
                s, q, r = t.b, t.c, t.a
                ns, nq, nr = t.nb, t.nc, t.na
            else # ret == TRI_NEIGHBOR_C
                s, q, r = t.c, t.a, t.b
                ns, nq, nr = t.nc, t.na, t.nb
            end
            t₁ = DelaunayTriangle(s, q, vidx, 0, base_idx+1, nr, true)
            t₂ = DelaunayTriangle(s, vidx, r, 0, nq, base_idx+0, true)
            (nr > 0) && _update_neighbor(tess.faces[nr], i, base_idx+0)
            (nq > 0) && _update_neighbor(tess.faces[nq], i, base_idx+1)
            push!(tess.faces, t₁)
            push!(tess.faces, t₂)
            push!(edges, TriEdge(base_idx+0, TRI_NEIGHBOR_C))
            push!(edges, TriEdge(base_idx+1, TRI_NEIGHBOR_B))
            t.active = false
        end
        # adjust neighboring triangles
        tess.faces[end - 3].na = length(tess.faces) - 0
        tess.faces[end - 2].na = length(tess.faces) - 1
        tess.faces[end - 1].na = length(tess.faces) - 2
        tess.faces[end - 0].na = length(tess.faces) - 3

        _flip(tess, edges)

        if !_check_nbr(tess)
            #break
        end
    end
end

# update neighbor indices for surrounding triangles
function _update_neighbor(t::DelaunayTriangle, old_idx::Int, new_idx::Int) where T <: Point2
    if t.na == old_idx
        t.na = new_idx
    elseif t.nb == old_idx
        t.nb = new_idx
    else
        t.nc = new_idx
    end
end

function _shares_common_edge(tess::DelaunayTess2D{T}, t::DelaunayTriangle,
    n::DelaunayTriangle) where T <: Point2

    p_ta, p_tb, p_tc = tess.verts[[t.a, t.b, t.c]]
    p_na, p_nb, p_nc = tess.verts[[n.a, n.b, n.c]]
    length(unique([p_ta, p_tb, p_tc, p_na, p_nb, p_nc])) == 4
end

function _check_nbr_single(tess::DelaunayTess2D{T}, t::DelaunayTriangle) where T<: Point2
    good = true
    if t.na > 0
        if !_shares_common_edge(tess, t, tess.faces[t.na])
            good = false
            println("********************** BAD A")
            println("**********************-> T ", t)
            println("**********************-> N ", tess.faces[t.na])
        end
    end
    if t.nb > 0
        if !_shares_common_edge(tess, t, tess.faces[t.nb])
            good = false
            println("********************** BAD B")
            println("**********************-> T ", t)
            println("**********************-> N ", tess.faces[t.nb])
        end
    end
    if t.nc > 0
        if !_shares_common_edge(tess, t, tess.faces[t.nc])
            good = false
            println("********************** BAD C")
            println("**********************-> T ", t)
            println("**********************-> N ", tess.faces[t.nc])
        end
    end
    return good
end

function _check_nbr(tess::DelaunayTess2D{T}) where T <: Point2
    good = true
    for t in tess.faces
        if !t.active
            continue
        end
        good = good && _check_nbr_single(tess, t)
    end
    return good
end

# Lawson flip algorithm
function _flip(tess::DelaunayTess2D{T}, queue_init::Vector{TriEdge}) where T <: Point2
    seen::Set{TriEdge} = Set()
    queue::Vector{TriEdge} = queue_init

    # init queue
    for offset in 0:2
        idx = length(tess.faces) - offset
        t = tess.faces[idx]
        push!(queue, TriEdge(idx, TRI_NEIGHBOR_C))
    end 

    while length(queue) > 0
        edge = pop!(queue)
        it = edge.tri_idx; nbr = edge.nbr
        t = tess.faces[it]
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

        if !isquadconvex(tess.verts[u], tess.verts[c1], tess.verts[v], tess.verts[c2])
            continue # nothing to do
        end

        pt_a, pt_b, pt_c, pt_v = tess.verts[[t.a, t.b, t.c, v]]
        if incircumcircle(pt_a, pt_b, pt_c, pt_v)
            edge1, edge2 = _flip_tri(tess, t, n, it, in, v)
            push!(queue, edge1)
            push!(queue, edge2)
        end
    end
end

function _flip_tri(tess::DelaunayTess2D{T}, t::DelaunayTriangle,
    n::DelaunayTriangle, it::Int64, in::Int64, v::Int64) where T <: Point2
    # flip diagram:
    #               ---->
    #        +(X)               +
    #       /|\                / \
    #      / | \         nY_b /   \  nY_a
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

    # new edges that are not incident to `u` need to be checked
    return (TriEdge(t_idx, TRI_NEIGHBOR_C), TriEdge(n_idx, TRI_NEIGHBOR_A))
end

function _deactivate_extremal_triangles(tess::DelaunayTess2D{T}) where T <: Point2
    for t in tess.faces
        # first three vertices are extremal pts
        extremal = (t.a == 1 || t.a == 2 || t.a == 3 ||
                    t.b == 1 || t.b == 2 || t.b == 3 ||
                    t.c == 1 || t.c == 2 || t.c == 3)
        t.active = t.active && !extremal
    end
end

function _find_tri_idx(tess::DelaunayTess2D{T}, pt::T) where T <: Point2
    t1 = 0
    ret1 = 0
    for (i, tri) in enumerate(tess.faces)
        if tri.active
            a, b, c = tess.verts[[tri.a, tri.b, tri.c]]
            ret = intriangle(a, b, c, pt)
            if ret == 0
                return (i, 0, ret, 0)
            elseif ret > 0
                if t1 == 0
                    t1 = i
                    ret1 = ret
                else
                    return (t1, i, ret1, ret)
                end
            end
        end
    end

    # means a triangle should be split but we only found one triangle
    if t1 != 0
        println("problem!")
        println("]]] ", pt)
        t = tess.faces[t1]
        println("]]] ", t)
        println("]]] NA -> ", tess.faces[t.na])
        println("]]] NB -> ", tess.faces[t.nb])
        println("]]] NC -> ", tess.faces[t.nc])
    end
    #@assert t1 == 0
    return (0, 0)
end

function _coordinate_bound(V::AbstractVector{T}) where T <: Point2
    maxcoord = -Inf
    for v in V
        x, y = abs(v[1]), abs(v[2])
        (x > maxcoord) && (maxcoord = x)
        (y > maxcoord) && (maxcoord = y)
    end
    return maxcoord
end

function conformingDelaunay2D(V::AbstractVector{T}, pslg::PSLG) where T <: Point2
    tess = delaunay2D(V)

    for seg in pslg.segments
        if !_edge_in_tess(V, tess, seg)
            println("Segment ", seg, " not in tess")
            _insert_segment(V, tess, seg)
        end
    end

    return tess
end

# FIX: use indexed vertices
function _edge_in_tess(V::AbstractVector{T}, tess::DelaunayTess2D{T},
    seg::IndexedLineSegment) where T <: Point2

    for t in tess.faces
        if !t.active
            continue
        end

        pa, pb = V[seg.a], V[seg.b]
        r1 = intriangle(t, pa)
        r2 = intriangle(t, pb)
        if r1 > 0 && r2 > 0
            if _edge_in_tri(t, pa, pb)
                return true
            end
        end
    end

    return false
end

# FIX: use indexed vertices
function _edge_in_tri(tri::DelaunayTriangle, p::T, q::T) where T <: Point2
    if tri.a == p
        return tri.b == q || tri.c == q
    elseif tri.b == p
        return tri.a == q || tri.c == q
    elseif tri.c == p
        return tri.a == q || tri.b == q
    end
    return false
end

# FIX: use indexed vertices
function _insert_segment(V::AbstractVector{T}, tess::DelaunayTess2D{T},
    seg::IndexedLineSegment) where T <: Point2

    pa, pb = V[seg.a], V[seg.b]
    pc = (pa + pb)/2

    i1, i2, ret1, ret2 = _find_tri_idx(tess, pc)
    @assert i1 != 0

    t = tess.faces[i1]
    if _edge_in_tri(t, pa, pb)
        println("Segment insertion ret!")
        return
    end

    l1_bc, l2_bc, pc_bc = lineintersection(t.b, t.c, pa, pb)
    l1_ac, l2_ac, pc_ac = lineintersection(t.a, t.c, pa, pb)
    l1_ab, l2_ab, pc_ab = lineintersection(t.a, t.b, pa, pb)

    if t.a == pa || t.a == pb
        @assert 0 <= l1_bc <= 1 && 0 <= l2_bc <= 1
        pc = pc_bc
    elseif t.b == pa || t.b == pb
        @assert 0 <= l1_ac <= 1 && 0 <= l2_ac <= 1
        pc = pc_ac
    elseif t.c == pa || t.c == pb
        @assert 0 <= l1_ab <= 1 && 0 <= l2_ab <= 1
        pc = pc_ab
    else
        if 0 <= l1_bc <= 1
            pc = pc_bc
        elseif 0 <= l1_ac <= 1
            pc = pc_ac
        elseif 0 <= l1_ab <= 1
            pc = pc_ab
        else
            @assert false
        end
    end

    push!(V, pc)
    pc_idx = length(V)
    _insert_point(tess, pc)

    seg1 = IndexedLineSegment(seg.a, pc_idx)
    seg2 = IndexedLineSegment(pc_idx, seg.b)

    _insert_segment(V, tess, seg1)
    _insert_segment(V, tess, seg2)
end
