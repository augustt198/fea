using GeometryTypes
using SparseArrays
using LinearAlgebra

struct FEAMesh{T <: Point2} 
    tess::DelaunayTess2D{T}
    pslg::PSLG
    vertex_indexing::AbstractVector{Int64}
end

function createFEAMesh(V::AbstractVector{T}, pslg::PSLG) where T <: Point2
    tess = conformingDelaunay2D(V, pslg)

    boundary_vec = Vector{Bool}(undef, length(tess.verts))
    boundary_vec .= false
    for seg in pslg.segments
        if seg.boundary == 1
            boundary_vec[seg.a+3] = true
            boundary_vec[seg.b+3] = true
        end
    end

    vertex_indexing = Vector{Int64}(undef, length(boundary_vec))
    i = 0
    for (j, b) in enumerate(boundary_vec)
        vertex_indexing[j] = (j <= 3 || b) ? 0 : (i += 1)
    end

    FEAMesh(tess, pslg, vertex_indexing)
end

function _get_tri_vert_id(t::DelaunayTriangle, v::Int64)
    if t.a == v
        return TRI_VERTEX_A
    elseif t.b == v
        return TRI_VERTEX_B
    elseif t.c == v
        return TRI_VERTEX_C
    else
        return 0
    end
end

_next_vert_id_acw(i::Int64) = (i % 3) + 1

function _get_tri_vert_by_id(t::DelaunayTriangle, id::Int64)
    if id == TRI_VERTEX_A
        return t.a
    elseif id == TRI_VERTEX_B
        return t.b
    else
        return t.c
    end
end

function integrateVert(mesh::FEAMesh{T}, v::Int64, f::Function, A, F) where T <: Point2
    pv = mesh.tess.verts[v]
    start_tri_idx = 0
    for (i, t) in enumerate(mesh.tess.faces)
        if t.active
            id = _get_tri_vert_id(t, v)
            if id > 0
                start_tri_idx = i
                break
            end
        end
    end
    @assert start_tri_idx != 0

    load_val = 0.0
    self_val = 0.0
    curr_tri_idx = start_tri_idx
    mat_idx = mesh.vertex_indexing[v]

    while true
        t = mesh.tess.faces[curr_tri_idx]
        pa, pb, pc = mesh.tess.verts[[t.a, t.b, t.c]]
        curr_id = _get_tri_vert_id(t, v)

        load_integrand = (x) -> begin
            f(x) * barycentric2(pa, pb, pc, x)[curr_id]
        end
        load_contrib, err = integratetri(pa, pb, pc, load_integrand)
        load_val += load_contrib

        self_integrand = (x) -> begin
            ∇W = barycentricgrad(pa, pb, pc, x)
            return ∇W[curr_id]' * ∇W[curr_id]
        end
        self_contrib, err = integratetri(pa, pb, pc, self_integrand)
        self_val += self_contrib

        opp_id = curr_id
        for _iter in 1:2
            opp_id = _next_vert_id_acw(opp_id)
            stiffness_integrand = (x) -> begin
                ∇W = barycentricgrad(pa, pb, pc, x)
                return  ∇W[curr_id]' * ∇W[opp_id]
            end
            opp_vidx       = _get_tri_vert_by_id(t, opp_id)
            opp_mat_idx    = mesh.vertex_indexing[opp_vidx]
            if opp_mat_idx > 0
                stiffness_contrib, err = integratetri(pa, pb, pc, stiffness_integrand)
                A[mat_idx, opp_mat_idx] += stiffness_contrib
            end
        end

        next_tri_idx = 0
        if curr_id == TRI_VERTEX_A
            next_tri_idx = t.nb
        elseif curr_id == TRI_VERTEX_B
            next_tri_idx = t.nc
        elseif curr_id == TRI_VERTEX_C
            next_tri_idx = t.na
        end
        (next_tri_idx == start_tri_idx) && break
        curr_tri_idx = next_tri_idx
    end

    A[mat_idx, mat_idx] = self_val
    F[mat_idx] = load_val
end

function assemblePoisson(mesh::FEAMesh{T}, f::Function) where T <: Point2

    internal_verts   = filter(x -> x > 0, mesh.vertex_indexing)
    n_internal_verts = length(internal_verts)
    A    = spzeros(n_internal_verts, n_internal_verts)
    load = zeros(n_internal_verts)

    for (vert_i, mat_i) in enumerate(mesh.vertex_indexing)
        if mat_i > 0
            integrateVert(mesh, vert_i, f, A, load)
        end
    end

    println(">>>> IS POS DEF? ", isposdef(A))
    display(A)

    return A, load
end
