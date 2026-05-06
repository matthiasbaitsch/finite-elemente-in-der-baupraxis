using LinearAlgebra

function ab(e)
    x = coordinates(e)
    u1 = x[:, 2] - x[:, 1]
    u2 = x[:, 3] - x[:, 4]
    v1 = x[:, 4] - x[:, 1]
    v2 = x[:, 3] - x[:, 2]
    a = norm(u1)
    b = norm(u2)

    @assert norm(u1 - u2) < 1e-10
    @assert norm(v1 - v2) < 1e-10
    @assert (u1 ⋅ v1) < 1e-10
    @assert a > 1e-10
    @assert b > 1e-10
    return a, b
end

include("generated/plate-rrm-element.jl")
include("generated/plate-rkpc-element-conforming.jl")
include("generated/plate-rkpn-element-nonconforming.jl")

