using GeometryTypes
using SparseArrays
using LinearAlgebra

# Stores a material function for each region, indexed
# by region ID.
struct MaterialFunctionsList
    funcs::AbstractVector{Function}
end

MaterialFunctionsList() = MaterialFunctionsList(Vector{Function}(undef,0))

NodeBoundaryType = Int64
const NODE_BOUNDARY_TYPE_NONE = 0
const NODE_BOUNDARY_TYPE_DIRICHLET = 1
const NODE_BOUNDARY_TYPE_NEUMANN = 2

struct NodeBoundaryCondition
    type::NodeBoundaryType
    value::Float64
end

# default node: not a boundary
NodeBoundaryCondition() = NodeBoundaryCondition(NODE_BOUNDARY_TYPE_NONE, 0.0)

# vertex_indexing: maps vertex index to load vector index
#                  (or zero) if not in load vector
# tri_region_id: maps triangle index to region ID
struct FEAMesh{T <: Point2} 
    tess::DelaunayTess2D{T}
    pslg::PSLG
    vertex_indexing::AbstractVector{Int64}
    tri_region_id::AbstractVector{Int64}
    materials::MaterialFunctionsList
    boundary_conditions::AbstractVector{NodeBoundaryCondition}
end

function createFEAMesh(V::AbstractVector{T}, pslg::PSLG,
    bdconds::AbstractVector{NodeBoundaryCondition},
    materials::MaterialFunctionsList) where T <: Point2

    tess = conformingDelaunay2D(V, pslg)

    for (i,b) in enumerate(bdconds)
        println(i, " -> ", b)
    end

    vertex_indexing = Vector{Int64}(undef, length(tess.verts))
    vertex_indexing[1:3] .= 0
    i = 0
    for j in 4:length(vertex_indexing)
        bd_idx = j - 3
        if bd_idx <= length(bdconds)
            bdc = bdconds[j - 3]
        else
            bdc = NodeBoundaryCondition()
        end
        if bdc.type == NODE_BOUNDARY_TYPE_DIRICHLET
            vertex_indexing[j] = 0
        else
            vertex_indexing[j] = (i += 1)
        end
    end

    tri_region_id = Vector{Int64}(undef, length(tess.faces))
    for (i, t) in enumerate(tess.faces)
        if t.active
            a, b, c = tess.verts[[t.a, t.b, t.c]]
            centroid = (a + b + c)/3
            R = findregion(tess.verts, pslg.segments, centroid)
            tri_region_id[i] = R
        else
            tri_region_id[i] = 0
        end
    end

    FEAMesh(tess, pslg, vertex_indexing, tri_region_id, materials, bdconds)
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


    first = true
    while true
        (!first && curr_tri_idx == start_tri_idx) && break

        t         = mesh.tess.faces[curr_tri_idx]
        region_id = mesh.tri_region_id[curr_tri_idx]
        curr_id   = _get_tri_vert_id(t, v)

        # next triangle in ACW cycle
        next_tri_idx = 0
        if curr_id == TRI_VERTEX_A
            next_tri_idx = t.nb
        elseif curr_id == TRI_VERTEX_B
            next_tri_idx = t.nc
        elseif curr_id == TRI_VERTEX_C
            next_tri_idx = t.na
        end
        curr_tri_idx = next_tri_idx

        if !t.active || region_id < 1
            println("skipping!")
            continue
        end

        pa, pb, pc = mesh.tess.verts[[t.a, t.b, t.c]]
         # material function
        k = mesh.materials.funcs[region_id]

        println(">>> ", t)

        load_integrand = (x) -> begin
            f(x) * barycentric2(pa, pb, pc, x)[curr_id]
        end
        load_contrib, err = integratetri(pa, pb, pc, load_integrand)
        load_val += load_contrib

        self_integrand = (x) -> begin
            ∇W = barycentricgrad(pa, pb, pc, x)
            return k(x) * ∇W[curr_id]' * ∇W[curr_id]
        end
        self_contrib, err = integratetri(pa, pb, pc, self_integrand)
        self_val += self_contrib

        opp_id = curr_id
        for _iter in 1:2
            opp_id      = _next_vert_id_acw(opp_id)
            opp_vidx    = _get_tri_vert_by_id(t, opp_id)
            opp_mat_idx = mesh.vertex_indexing[opp_vidx]

            stiffness_integrand = (x) -> begin
                ∇W = barycentricgrad(pa, pb, pc, x)
                return k(x) * ∇W[curr_id]' * ∇W[opp_id]
            end
            stiffness_contrib, err = integratetri(pa, pb, pc, stiffness_integrand)

            if opp_mat_idx > 0
                # neighboring node is free
                println("*** ", mat_idx, " --> ", opp_mat_idx, "   (", stiffness_contrib, ")")
                A[mat_idx, opp_mat_idx] += stiffness_contrib
            else
                println("*** DIRICHLET INT")
                opp_bdcond = mesh.boundary_conditions[opp_vidx - 3]
                @assert opp_bdcond.type == NODE_BOUNDARY_TYPE_DIRICHLET
                F[mat_idx] -= opp_bdcond.value * stiffness_contrib
            end
        end

        first = false
    end

    A[mat_idx, mat_idx] = self_val
    F[mat_idx]          += load_val

    if v - 3 < length(mesh.boundary_conditions)
        bdcond = mesh.boundary_conditions[v - 3]
        if bdcond.type == NODE_BOUNDARY_TYPE_NEUMANN
            println("Subtract neumann: ", bdcond.value)
            F[mat_idx] -= bdcond.value
        end
    end
end

function assemblePoisson(mesh::FEAMesh{T}, f::Function) where T <: Point2
    n_free_verts = count(x -> x > 0, mesh.vertex_indexing)
    A    = spzeros(n_free_verts, n_free_verts)
    load = zeros(n_free_verts)

    display(mesh.vertex_indexing)

    for (vert_i, mat_i) in enumerate(mesh.vertex_indexing)
        if mat_i > 0
            integrateVert(mesh, vert_i, f, A, load)
        end
    end

    println(">>>> IS POS DEF? ", isposdef(A))
    display(Matrix(A))
    println("LOAD: ")
    display(load)

    return A, load
end
