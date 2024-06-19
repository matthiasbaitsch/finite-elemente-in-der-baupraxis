using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Plots
using MMJMesh.Utilities
using MMJMesh.Geometries
using MMJMesh.Mathematics
using MMJMesh.Mathematics: domain, domaintype
using MMJMesh.Topologies
using MMJMesh.MMJBase

import MMJMesh.Mathematics.FixedPolynomials as FP
using LinearAlgebra
using IntervalSets

using Symbolics
using SymbolicUtils
using SymbolicUtils: Postwalk, Chain

ei(n, i) = [j == i ? 1 : 0 for j = 1:n]

##
@variables a, b, E, ν, h;

## Conforming shape functions 
V = [0 a a 0; 0 0 b b];
P = mmonomials(2, 3, QHat, type=BigInt)
N = vcat(
    [
        [
            ValueAtLF(p),
            PDerivativeAtLF(p, [1, 0]),
            PDerivativeAtLF(p, [0, 1]),
            PDerivativeAtLF(p, [1, 1])
        ]
        for p in eachcol(V)
    ]...
)
M = [n(p) for p in P, n in N]
invM = hcat([M \ ei(length(N), i) for i = 1:length(N)]...)
H4 = invM * P

## Non conforming shape functions
V = [0 a a 0; 0 0 b b];
P = mmonomials(2, 3, QHat, (p1, p2) -> p1 + p2 <= 4 && p1 * p2 < 4, type=BigInt)
N = vcat(
    [
        [
            ValueAtLF(p),
            PDerivativeAtLF(p, [1, 0]),
            PDerivativeAtLF(p, [0, 1])
        ]
        for p in eachcol(V)
    ]...
)
M = [n(p) for p in P, n in N]
invM = hcat([M \ ei(length(N), i) for i = 1:length(N)]...)
H4i = invM * P

##
D = E * h^3 / (12 * (1 - ν^2))
C = [1 ν 0; ν 1 0; 0 0 (1-ν)/2]

B(w) = [∂xx(w), ∂yy(w), 2 * ∂xy(w)]

ae(w, δw) = simplifyx(integrate((B(w) ⋅ (C * B(δw))), 0 .. a, 0 .. b))
be(δw) = simplifyx(integrate(δw, 0 .. a, 0 .. b))

function computeKe(H)
    ne = length(H)
    Ke = zeros(Num, ne, ne)
    Threads.@threads for i = 1:ne
        for j = i:ne
            Ke[i, j] = ae(H[i], H[j])
            if i != j
                Ke[j, i] = Ke[i, j]
            end
        end
    end
    return Ke
end

computeRe(H) = [be(h) for h in H]

## Incompatible case
@time Kei = computeKe(H4i)
@time rei = computeRe(H4i)


# Compare to reference solution (Vanam et al.)
@variables α, β
sol = 1 / (4 * a * b) * (4 * (α^2 + β^2) + 2 // 3 * (7 - 2 * ν))
substitute(sol, Dict(α => a / b, β => b / a)) |> simplify |> println
Kei[1, 1] |> println
println("Somethings wrong with my code or with Symbolics.jl")


## Compatible case
@time Ke = computeKe(H4)
@time re = computeRe(H4)


## Write code

function writecode(name, a, Ke, b, re, io::IO=stdout)

    ne = size(Ke, 1)

    println(io,
        "function $(name)Ke(h, E, ν)
    function keFunc(e)
        a, b = ab(e)
        Ke = zeros($ne, $ne)
        D = E*h^3 / (12*(1 - ν^2))"
    )

    for i = 1:size(Ke, 1), j = i:size(Ke, 2)
        if i == j
            println(io, "        Ke[$i, $j] = $(Ke[i, j])")
        else
            println(io, "        Ke[$i, $j] = Ke[$j, $i] = $(Ke[i, j])")
        end
    end

    println(io,
        "        return D*Ke
    end
    return keFunc
end
"
    )

    println(
        io,
        "
function $(name)Re(q)
    function reFunc(e)
        a, b = ab(e)
        re = zeros($ne)"
    )
    for i = 1:length(re)
        println(io, "        re[$i] = $(simplifyx(re[i]))")
    end
    println(
        io,
        "        return q / 144*re
    end
    return reFunc
end
"
    )
end

out = open("04-platten/folien/src/rkp-conforming.jl", "w")
writecode("rkp", 1, Ke, 1, 144 * re, out)
close(out)
