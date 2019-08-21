# idk the actual way to do this sorry
push!(LOAD_PATH, string(@__DIR__, "/../src/"))

using GeometryTypes
using FEA
import IterativeSolvers
import Random

include("plottess.jl")

function test()
    pts = Vector{Point2f0}(undef, 0)
    bdconds = Vector{FEA.NodeBoundaryCondition}(undef, 0)

    outer_tuples = [(1.0, 1.0), (1.0, 0.75), (1.0, 0.5), (1.0, 0.25), (1.0, 0.0), (0.5, 0.0),
        (0.0, 0.0), (0.0, -0.5), (0.0, -1.0), (-0.25,-1.0), (-0.5,-1.0), (-0.75,-1.0), (-1.0, -1.0), (-1.0, 0.0), (-1, 1.0), (0.0, 1.0)]
    for (x, y) in outer_tuples
        push!(pts, Point2f0(x, y))
        if x == 1.0
            type = FEA.NODE_BOUNDARY_TYPE_DIRICHLET
            val = 2.0
        elseif y == -1
            type = FEA.NODE_BOUNDARY_TYPE_DIRICHLET
            val = 1.0
        else
            type = FEA.NODE_BOUNDARY_TYPE_NEUMANN
            val = 0.0
        end
        bd = FEA.NodeBoundaryCondition(type, val)
        push!(bdconds, bd)
    end

    inner_tuples = [(0.75,0.75), (0.75,0.25),(-0.25,0.25),(-0.25,-0.75),(-0.75,-0.75),(-0.75,0.75)]
    for (x, y) in inner_tuples
        push!(pts, Point2f0(x, y))
        type = FEA.NODE_BOUNDARY_TYPE_NONE
        bd = FEA.NodeBoundaryCondition(type, 0.0)
        push!(bdconds, bd)
    end
    
    segments = Vector{FEA.IndexedLineSegment}(undef, 0)
    for i in 1:length(outer_tuples)-1
        j = (i % length(outer_tuples)) + 1
        seg = IndexedLineSegment(i, j, 1)
        push!(segments, seg)
    end

    for i in 1:length(inner_tuples)
        j = (i % length(inner_tuples)) + 1
        offset = length(outer_tuples)
        seg = IndexedLineSegment(i+offset, j+offset, 2)
        println("Create seg: ", pts[seg.a], " ", pts[seg.b])
        push!(segments, seg)
    end

    Random.seed!(4)
    for i in 1:100
        pt = Point2f0(rand(2)) * 2 - 1

        if 0.1 <= pt[1] <= 0.9 && 0.1 <= pt[2] <= 0.9
            push!(pts, pt)
        elseif -0.9 <= pt[1] <= -0.1 && 0.1 <= pt[2] <= 0.9
            push!(pts, pt)
        elseif -0.9 <= pt[1] <= -0.1 && -0.9 <= pt[2] <= -0.1
            push!(pts, pt)
        end
    end
    
    pslg = FEA.PSLG(segments)

    f = (x) -> 0.0
    k₁(x) = 1.0
    k₂(x) = 1.0
    materials = FEA.MaterialFunctionsList([k₁, k₂])
    mesh = FEA.createFEAMesh(pts, pslg, bdconds, materials)
    A, F = FEA.assemblePoisson(mesh, f)
    FEA.deactivate_external(mesh.tess, mesh.pslg)

    U = IterativeSolvers.cg(A, F)
    display(U)
    println("max min ", maximum(U), " ", minimum(U))

    node_colors = Vector{Float64}(undef, length(mesh.tess.verts))

    for (n_idx, m_idx) in enumerate(mesh.vertex_indexing)
        if m_idx > 0
            node_colors[n_idx] = U[m_idx]
        elseif n_idx > 3
            bd = bdconds[n_idx - 3]
            println("** BD: ", bd)
            node_colors[n_idx] = bd.value
        else
            node_colors[n_idx] = 0
        end
    end

    println(node_colors)

    p1 = plottess(mesh.tess, false, node_colors)
    scene = vbox(
        p1,
        colorlegend(
            p1[end-1],            # get Plot object from Scene
            camera = campixel!, # let vbox decide scene limits
            raw = true,          # no axes, other things as well
            width = (            # make the colorlegend longer so it looks nicer
                30,              # the width
                540              # the height
            )
        )
    )
end

test()
