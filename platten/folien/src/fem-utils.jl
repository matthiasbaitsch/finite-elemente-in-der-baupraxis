function makenodalbasis(V, P, NP)
    nn = length(P)
    N = vcat([NP(p) for p in eachcol(V)]...)
    M = [n(p) for p in P, n in N]
    I = [i == j for i = 1:nn, j = 1:nn]
    return (M \ I) * P
end

function makevectorbasis(bs, nf)
    f1 = bs[1]
    n = length(bs)
    dt = domaintype(f1)
    ct = codomaintype(f1)
    d = domain(f1)
    null = MMJMesh.Mathematics.Zero{dt,ct,d}()
    return reshape([[j == k ? bs[i] : null for k = 1:nf] for j = 1:nf, i = 1:n], (:, 1))
end

function computeKe(ae, H)
    ne = length(H)
    Ke = zeros(Num, ne, ne)
    Threads.@threads for i = 1:ne
        for j = i:ne
            Ke[i, j] = ae(H[i], H[j])
        end
    end
    return Symmetric(Ke)
end

computeRe(be, H) = [be(h) for h in H]

function writecode(name, Ke, re, out::IO=stdout)
    ne = size(Ke, 1)

    println(out,
        "function $(name)Ke(h, E, ν)
    function keFunc(e)
        a, b = ab(e)
        Ke = zeros($ne, $ne)
        D = E*h^3 / (12*(1 - ν^2))")
    print(out, MatrixR("Ke", Ke))
    println(out,
        "        return D*Symmetric(Ke)
    end
    return keFunc
end")
    println(out)
    println(out,
        "function $(name)Re(q)
    function reFunc(e)
        a, b = ab(e)
        re = zeros($ne)"
    )
    for i = 1:length(re)
        println(out, "        re[$i] = $(simplifyx(re[i]))")
    end
    println(out,
        "        return q*re
    end
    return reFunc
end")
end
