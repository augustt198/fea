# idk the actual way to do this sorry
push!(LOAD_PATH, string(@__DIR__, "/../src/"))

using GeometryTypes
using FEA

import IterativeSolvers
import Random
using PyCall

include("plottess.jl")

function create_pslg(fname)
    ezdxf = pyimport("ezdxf")
    drawing = ezdxf.readfile(fname)

    region_start_pt  = nothing
    region_start_idx = nothing
    last_end  = nothing
    region_id = 0
    vertices  = Vector{Point2{Float64}}(undef, 0)
    segments  = Vector{FEA.IndexedLineSegment}(undef, 0)
    for e in drawing.entities
        if e.dxftype() == "LINE"
            if e.dxf.start != last_end
                region_id += 1
                push!(vertices, Point2{Float64}(e.dxf.start))
                region_start_pt = e.dxf.start
                region_start_idx = length(vertices)
            end
            
            start_idx = length(vertices)
            if e.dxf.end == region_start_pt
                end_idx = region_start_idx
            else
                push!(vertices, Point2{Float64}(e.dxf.end))
                end_idx = length(vertices)
            end
            seg = FEA.IndexedLineSegment(start_idx, end_idx, region_id)
            push!(segments, seg)

            last_end = e.dxf.end
        end
    end

    println("GOOD ?")
    println(length(unique(vertices)) == length(vertices))
    println(length(vertices))

    println("> REGION ID: ", region_id)
    println("> SEG LEN ", length(segments))

    fn = (x) -> print("($(x[1]), $(x[2])),")
    #map(fn, vertices)
    println()

    return vertices, segments
end

function test()
    vertices, segments = create_pslg("/Users/August/Documents/cad/6A01/rotor_hub.dxf")
    #vertices, segments = create_pslg("/Users/August/Documents/cad/side1_thing.dxf")
    #vertices, segments = create_pslg("/Users/August/Documents/cad/crane/dxf/support_1.dxf")

    k = (x) -> 1.0
    materials = FEA.MaterialFunctionsList([k, k])
    bdconds = Vector{FEA.NodeBoundaryCondition}(undef, 0)

    max_y = maximum(map((x) -> x[2], vertices))
    min_y = minimum(map((x) -> x[2], vertices))
    for v in vertices
        dist = sqrt(v' * v)
        if dist > 0.4
            th = atan(v[2], v[1])
            val = 4 * sin(5 * th)
            bd = FEA.NodeBoundaryCondition(FEA.NODE_BOUNDARY_TYPE_DIRICHLET, val)
        elseif abs(v[1]) > 0.1
            bd = FEA.NodeBoundaryCondition(FEA.NODE_BOUNDARY_TYPE_DIRICHLET, 0.0)
        elseif abs(v[2]) > 0.1
            bd = FEA.NodeBoundaryCondition(FEA.NODE_BOUNDARY_TYPE_DIRICHLET, 0.0)
        else
            bd = FEA.NodeBoundaryCondition(FEA.NODE_BOUNDARY_TYPE_NEUMANN, -1.0)
        end
        #if v[2] == max_y && v[1] < 5
        #    bd = FEA.NodeBoundaryCondition(FEA.NODE_BOUNDARY_TYPE_DIRICHLET, 1.0)
        #elseif v[2] == min_y && v[1] > 22.5
        #    bd = FEA.NodeBoundaryCondition(FEA.NODE_BOUNDARY_TYPE_DIRICHLET, -1.0)
        #else
        #    bd = FEA.NodeBoundaryCondition(FEA.NODE_BOUNDARY_TYPE_NEUMANN, 0.0)
        #end
        push!(bdconds, bd)
    end

    mapfn = (x) -> begin
        if x.boundary > 1
            return FEA.IndexedLineSegment(x.a, x.b, 1)
        else
            return x
        end
    end
    segments = map(mapfn, segments)

    pslg = FEA.PSLG(segments)

    mesh = FEA.createFEAMesh(vertices, pslg, bdconds, materials)
    α = 30.0 * π / 180.0
    FEA.refineMesh(mesh.tess, α)
    FEA.deactivate_external(mesh.tess, pslg)

    f = (x) -> 0.0
    #A, F = FEA.assemblePoisson(mesh, f)

    #=
    println()
    @time U = IterativeSolvers.cg(A, F)
    println()
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
    =#

    p1 = plottess(mesh.tess, false)
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
