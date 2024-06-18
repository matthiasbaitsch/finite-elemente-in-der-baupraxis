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

# function isxnumeric(x)
#     for t = [Int32, Int64, BigInt, Float32, Float64, BigFloat, Rational]
#         if typeof(x) <: t
#             return true
#         end
#     end
#     return false
# end

# function toint(e)
#     @variables xone, xnull
#     isinteger(x) = isxnumeric(x) && round(x, digits=0) == x
#     e = simplify(e)
#     e = simplify(e, Postwalk(Chain([@rule ~x::isinteger => xone * (Int(~x) + xnull)])))
#     e = simplify(Symbolics.fixpoint_sub(e, Dict(xone => 1, xnull => 0)))
#     return e
# end

# function toint!(c::Vector{Num})
#     for i = eachindex(c)
#         c[i] = toint(c[i])
#     end
#     return c
# end

# function MMJMesh.Mathematics._simplify(exponents::Matrix, coefficients::Vector{Num})
#     return exponents, toint!(coefficients)
# end

ei(n, i) = [j == i ? 1 : 0 for j = 1:n]

function pmat(a)
    for i = 1:size(a, 1), j = 1:size(a, 2)
        println(i, ", ", j, ": ", a[i, j])
    end
end

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

B(w) = [-∂xx(w), -∂yy(w), 2 * ∂xy(w)]
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

function writecode(a, Ke, b, re)

    for i = 1:size(Ke, 1), j = i:size(Ke, 2)
        if i == j
            println("Ke[$i, $j] = $(Ke[i, j])")
        else
            println("Ke[$i, $j] = Ke[$j, $i] = $(Ke[i, j])")
        end
    end

    # for i = 1:length(re)
    #     println("re[$i] = $(re[i])")
    # end

end

