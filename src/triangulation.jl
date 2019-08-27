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
    last_tri_idx::Int64
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
    tess = DelaunayTess2D(V, [base_tri], 1)

    for i in 4:length(tess.verts)
        _insert_point(tess, i)
    end

    #_deactivate_extremal_triangles(tess)
    #_check_nbr(tess)
 
    return tess
end

function _insert_point(tess::DelaunayTess2D{T}, vidx::Int64, flip=true) where T <: Point2
    vert = tess.verts[vidx]
    i1, i2, ret1, ret2 = _find_tri_idx(tess, vert) # _find_tri_idx_fast(tess, vert)

    # TODO FIX
    (ret1 > 0 && i2 < 1) && return

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
        flip && _flip(tess, edges)
        
        false && @assert(_check_nbr(tess))
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

        flip && _flip(tess, edges)

        false && @assert(_check_nbr(tess))
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
function _flip(tess::DelaunayTess2D{T}, queue::Vector{TriEdge}) where T <: Point2
    seen::Set{TriEdge} = Set()

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

        if isquadconvex(tess.verts[u], tess.verts[c1], tess.verts[v], tess.verts[c2])
            pt_a, pt_b, pt_c, pt_v = tess.verts[[t.a, t.b, t.c, v]]
            if incircumcircle(pt_a, pt_b, pt_c, pt_v)
                edge1, edge2 = _flip_tri(tess, t, n, it, in, v)
                push!(queue, edge1)
                push!(queue, edge2)
            end
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

    # TODO: consider handling case when we're supposed to find
    # two triangles but one was find (need to ignore boundary?)
    return (t1, 0, ret1, 0)
end

function _find_tri_idx_fast(tess::DelaunayTess2D{T}, pt::T) where T <: Point2
    idx = tess.last_tri_idx
    n = 0
    while true
        n += 1
        #println("> [$n]: $idx")
        #flush(stdout)
        #n > 100 && error("nop")
        t = tess.faces[idx]
        #println("thingy: $t / $(tess.last_tri_idx)")
        if n > 100
            #println("i slep $n")
            #sleep(10)
        end
        if !t.active
            idx = length(tess.faces)
            t = tess.faces[idx]
        end
        println("now [$idx] $t")
        
        a, b, c = tess.verts[[t.a, t.b, t.c]]
        ret = intriangle(a, b, c, pt)
        #println("ret $ret / $a $b $c")
        if ret == -TRI_NEIGHBOR_A
            idx = t.na
        elseif ret == -TRI_NEIGHBOR_B
            idx = t.nb
            #println(">> B: $(t.nb)")
        elseif ret == -TRI_NEIGHBOR_C
            idx = t.nc
        elseif ret == 0
            tess.last_tri_idx = idx
            return (idx, 0, ret, 0)
        else
            if ret == TRI_NEIGHBOR_A
                idx2 = t.na
            elseif ret == TRI_NEIGHBOR_B
                idx2 = t.nb
            else
                idx2 = t.nc
            end
            t2 = tess.faces[idx2]
            a2, b2, c2 = tess.verts[[t2.a, t2.b, t2.c]]
            ret2 = intriangle(a2, b2, c2, pt)

            tess.last_tri_idx = idx
            return (idx, idx2, ret, ret2)
        end
    end
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
        seg_offset_idx = IndexedLineSegment(seg.a+3, seg.b+3)

        if !_edge_in_tess(tess, seg_offset_idx)
            println("Segment ", seg, " not in tess")
            # to adjust for added triangle
            _insert_segment(tess, seg_offset_idx)
        end
    end

    return tess
end

function _edge_in_tess(tess::DelaunayTess2D{T}, seg::IndexedLineSegment) where T <: Point2
    for t in tess.faces
        if t.active && _edge_in_tri(t, seg)
            return true
        end
    end

    return false
end

function _edge_in_tri(tri::DelaunayTriangle, seg::IndexedLineSegment) where T <: Point2
    if tri.a == seg.a
        return tri.b == seg.b || tri.c == seg.b
    elseif tri.b == seg.a
        return tri.a == seg.b || tri.c == seg.b
    elseif tri.c == seg.a
        return tri.a == seg.b || tri.b == seg.b
    end
    return false
end

function _insert_segment(tess::DelaunayTess2D{T}, seg::IndexedLineSegment) where T <: Point2

    pa, pb = tess.verts[[seg.a, seg.b]]
    pc = (pa + pb)/2 # midpoint

    i1, i2, ret1, ret2 = _find_tri_idx(tess, pc)
    @assert i1 != 0

    t = tess.faces[i1]
    if _edge_in_tri(t, seg)
        return
    end

    p_ta, p_tb, p_tc = tess.verts[[t.a, t.b, t.c]]
    l1_bc, l2_bc, pc_bc = lineintersection(p_tb, p_tc, pa, pb)
    l1_ac, l2_ac, pc_ac = lineintersection(p_ta, p_tc, pa, pb)
    l1_ab, l2_ab, pc_ab = lineintersection(p_ta, p_tb, pa, pb)

    ϵ = 1e-6
    if t.a == seg.a || t.a == seg.b
        @assert 0-ϵ <= l1_bc <= 1+ϵ && 0-ϵ <= l2_bc <= 1+ϵ
        pc = pc_bc
    elseif t.b == seg.a || t.b == seg.b
        @assert 0-ϵ <= l1_ac <= 1+ϵ && 0-ϵ <= l2_ac <= 1+ϵ
        pc = pc_ac
    elseif t.c == seg.a || t.c == seg.b
        @assert 0-ϵ <= l1_ab <= 1+ϵ && 0-ϵ <= l2_ab <= 1+ϵ
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

    push!(tess.verts, pc)
    pc_idx = length(tess.verts)
    _insert_point(tess, pc_idx, false)

    seg1 = IndexedLineSegment(seg.a, pc_idx)
    seg2 = IndexedLineSegment(pc_idx, seg.b)

    _insert_segment(tess, seg1)
    _insert_segment(tess, seg2)
end

function deactivate_external(tess::DelaunayTess2D{T}, pslg::PSLG) where T <: Point2
    for t in tess.faces
        if t.active
            v = sum(tess.verts[[t.a, t.b, t.c]]) / 3 
            r = findregion(tess.verts, pslg.segments, v)
            if r == 0
                t.active = false
            end
        end
    end
end

function _longest_side(tess::DelaunayTess2D{T}, t::DelaunayTriangle) where T <: Point2
    pa, pb, pc = tess.verts[[t.a, t.b, t.c]]
    vec_a, vec_b, vec_c = pc - pb, pc - pa, pb - pa
    dists = sqrt(vec_a' * vec_a), sqrt(vec_b' * vec_b), sqrt(vec_c' * vec_c)
    idx = argmax(dists)
    return idx, dists[idx]
end

# Rivara refinement
function _refine_tri(tess::DelaunayTess2D{T}, t::DelaunayTriangle) where T <: Point2
    triangles = [t]
    side, maxdist = _longest_side(tess, t)
    while true
        if side == TRI_NEIGHBOR_A
            t.na < 1 && break
            t = tess.faces[t.na]
        elseif side == TRI_NEIGHBOR_B
            t.nb < 1 && break
            t = tess.faces[t.nb]
        else
            t.nc < 1 && break
            t = tess.faces[t.nc]
        end
        side, dist = _longest_side(tess, t)
        a,b,c = tess.verts[[t.a, t.b, t.c]]
        ax,ay = a
        bx,by = b
        cx,cy = c
        if dist <= maxdist
            break
        else
            maxdist = dist
        end
    end

    mids = (tess.verts[[t.b, t.c, t.a]] + tess.verts[[t.c, t.a, t.b]]) / 2
    midpoint = mids[side]

    push!(tess.verts, midpoint)
    midpoint_idx = length(tess.verts)

    pa, pb, pc = tess.verts[[t.a, t.b, t.c]]

    _insert_point(tess, midpoint_idx)
end

function minangle(tess::DelaunayTess2D{T}, t::DelaunayTriangle) where T <: Point2
    pa, pb, pc = tess.verts[[t.a, t.b, t.c]]
    ab, ac = pb - pa, pc - pa
    ba, bc = pa - pb, pc - pb
    
    c₁ = (ab' * ac) / (sqrt(ab' * ab) * sqrt(ac' * ac))
    c₁ = max(min(c₁, 1), -1)
    c₂ = (ba' * bc) / (sqrt(ba' * ba) * sqrt(bc' * bc))
    c₂ = max(min(c₂, 1), -1)
    α₁ = acos(c₁)
    α₂ = acos(c₂)
    α₃ = π - α₁ - α₂
    return minimum((α₁, α₂, α₃))
end

function refineMesh(tess::DelaunayTess2D{T}, ε) where T <: Point2
    triangle_mask = Vector{Bool}(undef, length(tess.faces))
    for (i, t) in enumerate(tess.faces)
        triangle_mask[i] = (t.active && minangle(tess, t) < ε)
    end
    while any(triangle_mask)
        c = count(identity, triangle_mask)
        println("refine pass: $c")
        for (i, b) in enumerate(triangle_mask)
            if b && tess.faces[i].active
                m = tess.faces[i]
                _refine_tri(tess, tess.faces[i])
            end
        end

        triangle_mask = Vector{Bool}(undef, length(tess.faces))
        for (i, t) in enumerate(tess.faces)
            triangle_mask[i] = (t.active && minangle(tess, t) < ε)
        end
        c2 = count(identity, triangle_mask)
        if c2 == c
            break
        else
            c = c2
        end
    end
end
