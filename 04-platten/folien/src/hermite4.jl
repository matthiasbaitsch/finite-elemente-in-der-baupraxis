using MMJMesh
using MMJMesh.Mathematics

V = [-1 1 1 -1; -1 -1 1 1]
P = mmonomials(2, 3, QHat)
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
H4 = round.(inv(M), digits=15) * P

# fplot3d(H4)

