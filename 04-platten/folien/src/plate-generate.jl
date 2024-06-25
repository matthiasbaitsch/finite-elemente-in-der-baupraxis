using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Plots
using MMJMesh.MMJBase
using MMJMesh.Utilities
using MMJMesh.Topologies
using MMJMesh.Geometries
using MMJMesh.Mathematics
using MMJMesh.Mathematics: domain, domaintype

using Symbolics
using IntervalSets
using StaticArrays
using LinearAlgebra
using OrderedCollections


include("matrixr.jl")


## Helper functions

# XXX
# struct DifferentialOperator
#     op
# end
# Base.:(*)(op::DifferentialOperator, f::AbstractMapping) = op.op(f)
# op(p) = DifferentialOperator(p)

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


## Shape functions

@variables a, b, E, ν, h;
points = [0 a a 0; 0 0 b b];

# Kirchhoff conforming
H4c = makenodalbasis(
    points,
    mmonomials(2, 3, QHat, type=BigInt),
    p -> [ValueAtLF(p), PDerivativeAtLF(p, [1, 0]), PDerivativeAtLF(p, [0, 1]), PDerivativeAtLF(p, [1, 1])]
)

# Kirchhoff non conforming
H4n = makenodalbasis(
    points,
    mmonomials(2, 3, QHat, (p1, p2) -> p1 + p2 <= 4 && p1 * p2 < 4, type=BigInt),
    p -> [ValueAtLF(p), PDerivativeAtLF(p, [1, 0]), PDerivativeAtLF(p, [0, 1])]
)

# Reissner Mindlin
L4 = makenodalbasis(
    points,
    mmonomials(2, 1, QHat),
    p -> [ValueAtLF(p)]
)
L43 = makevectorbasis(L4, 3)


## Reissner Mindlin plate (Notizen Vorlesung Kassel)

k = 5 // 6
c1 = (1 - ν) / 2
c2 = k * 6 * (1 - ν) / h^2

D = [
    0 ∂x 0
    0 0 ∂y
    0 ∂y ∂x
    ∂x 1 0
    ∂y 0 1
]

C = [
    1 ν 0 0 0
    ν 1 0 0 0
    0 0 c1 0 0
    0 0 0 c2 0
    0 0 0 0 c2
]

ae(w, δw) = simplifyx(integrate(((D * w) ⋅ (C * (D * δw))), 0 .. a, 0 .. b))
be(δw) = simplifyx(integrate(δw[1], 0 .. a, 0 .. b))

## Element matrix and vector
Kerm = computeKe(ae, L43)
rerm = computeRe(be, L43)
out = open("04-platten/folien/src/generated/plate-rrm-element.jl", "w")
# writecode("rrm", Symmetric(simplify.(expand.(h^2*8*9*a*b*Kerm))), rerm, out)
writecode("rrm", Kerm, rerm, out)
close(out)


## Kirchhoff plate

C = [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
B(w) = [∂xx(w), ∂yy(w), 2 * ∂xy(w)]
ae(w, δw) = simplifyx(integrate((B(w) ⋅ (C * B(δw))), 0 .. a, 0 .. b))
be(δw) = simplifyx(integrate(δw, 0 .. a, 0 .. b))

# Conforming
@time Kec = computeKe(ae, H4c)
@time rec = computeRe(be, H4c)
out = open("04-platten/folien/src/generated/plate-rkpc-element-conforming.jl", "w")
writecode("rkpc", Kec, rec, out)
close(out)

# Non conforming
@time Ken = computeKe(ae, H4n)
@time ren = computeRe(be, H4n)
out = open("04-platten/folien/src/generated/plate-rkpn-element-nonconforming.jl", "w")
writecode("rkpn", Ken, ren, out)
close(out)

