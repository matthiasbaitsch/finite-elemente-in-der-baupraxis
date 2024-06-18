using GmshJL
using LinearAlgebra

function assembleKr(s)
    N = size(s.nodes, 2)
    Ne = size(s.elements, 2)
    K = zeros(N, N)
    r = zeros(N)    
    for e ∈ 1:Ne
        I = s.elements[:, e]
        x = s.nodes[:, I]
        K[I, I] += s.keFunc(x)
        r[I] += s.reFunc(x)
    end    
    return K, r
end

function assembleKr(s::FEMeshGroups)
    N = size(s.nodes, 2)
    K = zeros(N, N)
    r = zeros(N)    
    for g ∈ values(m.groups)
        for e ∈ 1:g.Ne
            I = g.elements[:, e]
            x = s.nodes[:, I]
            if(haskey(g, :keFunc))
                K[I, I] += g.keFunc(x)
            end
            if(haskey(g, :reFunc))
                r[I] += g.reFunc(x)
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
