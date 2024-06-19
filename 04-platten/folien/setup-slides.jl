using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Plots
using MMJMesh.Utilities

using Pkg
Pkg.add("Latexify")

using Revise
using Latexify
using CairoMakie
using GLMakie
using VarStructs
using LinearAlgebra

#WGLMakie.activate!()
GLMakie.activate!()

include("src/rkp.jl")
include("src/fem.jl")
include("src/hermite4.jl")
include("src/mplot3d.jl")

set_theme!(theme_minimal())
update_theme!(colormap=:blues)

function makefacefunction(w)
    function doit(face)
        idxs = dofs(nodeindices(face), 4)
        _, a, b = _fsize(face)
        t = repeat([1, a / 2, b / 2, a * b / 2], 4)
        return sum(w[idxs] .* t .* H4) # TODO make dot work
    end
    return doit
end

ei(n, i) = [j == i ? 1 : 0 for j = 1:n]

params = @var Params()
params.q = 5e3
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

function mkfig(
    ; a3d=true, w=600, h=400, title="",
    limits=(nothing, nothing, nothing)
)
    fig = Figure(size=(w, h))
    if a3d
        GLMakie.activate!()
        ax = Axis3(
            fig[1, 1],
            aspect=:data,
            title=title,
            viewmode=:stretch,
            perspectiveness=0.2,
            limits=limits,
            protrusions=0
        )
    else
        CairoMakie.activate!()
        ax = Axis(fig[1, 1], aspect=DataAspect(), title=title)
    end
    hidedecorations!(ax)
    hidespines!(ax)
    return fig
end

function plotw(
    m, ww;
    zs=2,
    edgesvisible=false,
    edgeslinewidth=2,
    mesh=3,
    a3d=true, w=600, h=400, title="",
    colorrange=MakieCore.automatic,
    colormap=MakieCore.theme(:colormap),
    limits=(nothing, nothing, nothing)
)
    fig = mkfig(a3d=a3d, w=w, h=h, title=title, limits=limits)
    mplot3d!(
        m, makefacefunction(ww),
        zscale=zs / maximum(ww),
        color=3,
        npoints=25,
        mesh=mesh,
        edgesvisible=edgesvisible,
        edgeslinewidth=edgeslinewidth,
        colorrange=colorrange,
        colormap=colormap
    )
    fig
end
