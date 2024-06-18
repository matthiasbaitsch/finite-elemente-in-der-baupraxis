import MMJMesh.Meshes: entities

dofs(idxs, nf) = collect(reshape([(i - 1) * nf + j for i = idxs, j = 1:nf]', :))


function assembleKr(m, n=1)
    N = n * nnodes(m)
    K = zeros(N, N)
    r = zeros(N)

    for d = 2:2
        for e ∈ entities(m, d)
            I = dofs(nodeindices(e), n)

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

function fixeddofs(dofs, fixed)
    cnt = 1
    nff = sum(fixed)
    nf = length(fixed)
    idxs = zeros(Int, nff * length(dofs))

    for i = dofs, j = 1:nf
        if fixed[j]
            idxs[cnt] = (i - 1) * nf + j
            cnt = cnt + 1
        end
    end

    return idxs
end

function applydirichletbcs!(dofs, K, r, fixed = [true])
    # pen = 1e8 * maximum(abs.(K))
    fdofs = fixeddofs(dofs, fixed)
    K[fdofs, :] .= 0
    r[fdofs] .= 0
    K[diagind(K)[fdofs]] .= 1 #.+= pen
    return nothing
end
