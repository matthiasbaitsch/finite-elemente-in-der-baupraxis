import MMJMesh.Meshes: entities

function assembleKr(m::Mesh)
    N = nnodes(m)
    K = zeros(N, N)
    r = zeros(N)

    for d = 0:3
        for e ∈ entities(m, d)
            I = nodeindices(e)
            
            kef = e.data[:kefunc]
            if !isnothing(kef)
                K[I, I] += kef(e)
            end
            
            ref = e.data[:refunc]
            if !isnothing(ref)
                r[I] += ref(e)
            end
        end
    end

    return K, r
end

function applyDirichletBCs!(dofs, K, r)
    K[dofs, :] .= 0
    r[dofs] .= 0
    K[diagind(K)[dofs]] .= 1
end
