using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Plots
using MMJMesh.Utilities

using WGLMakie
using CairoMakie

using Revise
using VarStructs
using LinearAlgebra

include("rkp.jl")
include("fem.jl")
include("hermite4.jl")

set_theme!(theme_minimal())
update_theme!(colormap=Reverse(:blues))


# Make and solve plate model

params = @var Params()

params.q = 8
params.nu = 0
params.d = 0.2
params.E = 31000e6

function plate(lx, ly, nx, ny, p)
    m = makemeshonrectangle(lx, ly, nx, ny)
    m.data[:kefunc] = rkpKe(p.d, p.E, p.nu)
    m.data[:refunc] = rkpRe(p.q)
    K, r = assembleKr(m, 4)
    applydirichletbcs!(m.groups[:boundarynodes], K, r, [true, true, true, false])
    w = K \ r
    return m, w
end


## Analysis of one plate

m, w = plate(10, 6, 20, 12, params)

mplot(m, edgesvisible=true) |> mconf()
mplot(m, edgesvisible=true, w[1:4:end]) |> mconf()
mplot(m, edgesvisible=true, w[2:4:end]) |> mconf()
mplot(m, edgesvisible=true, w[3:4:end]) |> mconf()
mplot(m, edgesvisible=true, w[4:4:end]) |> mconf()


## Convergence study for quadratic plate with lx = ly = 10

l = 10
hh = [];
nn = [];
ww = [];
for n = 4:2:30
    print("$n ")
    mn, wn = plate(l, l, n, n, params)
    push!(hh, l / n)
    push!(nn, 4 * nnodes(mn))
    push!(ww, maximum(abs.(wn[1:4:end])))
end
println()
scatterlines(hh, ww)


## Comparison to Czerny

w_fe = ww[end]
w_czerny = params.q * l^4 / (params.E * params.d^3) * 0.0152
100 * abs(w_fe - w_czerny) / w_czerny


## Plots

function fsize(face)
    x = coordinates(face)
    p = x[:, 1]
    l1 = x[1, 2] - x[1, 1]
    l2 = x[2, 3] - x[2, 2]
    return p, l1, l2
end

function plotmeshsolution(m, mf, cr, s)
    fig = Figure()
    ax = Axis3(fig[1, 1], protrusions=0, aspect=:data)
    hidedecorations!(ax)
    for f = faces(m)
        fplot3d!(
            mf(f), gmap=_makegmap(f), mesh=0, npoints=50, colorrange=cr, zscale=s
        )
    end
    fig
end

function makefacefunction(w)
    function doit(face)
        idxs = dofs(nodeindices(face), 4)
        _, a, b = fsize(face)
        t = repeat([1, a / 2, b / 2, a * b / 2], 4)
        return sum(w[idxs] .* t .* H4) # TODO make dot work
    end
    return doit
end


##

include("mplot3d.jl")

WGLMakie.activate!()
m, w = plate(10, 5, 4, 2, params)

fig = Figure()
Axis3(fig[1, 1], aspect=:data)
mplot3d!(
    m, makefacefunction(w),
    zscale=2 / maximum(w),
    color=3,
    npoints=25
)
fig


##

mf = makefacefunction(w)
g = mf(face(m, 4))

fig = Figure()
Axis3(fig[1, 1])
fplot3d!(g, npoints=200)
fig

