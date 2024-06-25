using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities

using MMJMesh.Plots
using MMJMesh.Plots: _fsize

using Revise
using Latexify
using GLMakie
using CairoMakie
using VarStructs
using SparseArrays
using LinearAlgebra

include("plate-elements.jl")
include("fem.jl")
include("hermite4.jl")
# include("mplot3d.jl")

set_theme!(theme_minimal())
update_theme!(faceplotmesh=5)
update_theme!(edgelinewidth=2.5)
update_theme!(colormap=:aquamarine)

GLMakie.activate!()

function makefacefunction(w)
    return face -> begin
        idxs = dofs(nodeindices(face), 4)
        _, a, b = _fsize(face)
        t = repeat([1, a / 2, b / 2, a * b / 4], 4)
        return sum(w[idxs] .* t .* H4) # TODO make dot work
    end
end

ei(n, i) = [j == i ? 1 : 0 for j = 1:n]

params = @var Params()
params.q = 5e3
params.nu = 0.0
params.d = 0.2
params.E = 31000e6

function plate(lx, ly, nx, ny, p, model=:kirchhoff_conforming)
    nf = 0
    bcs = []
    m = makemeshonrectangle(lx, ly, nx, ny)

    if model == :reissner_mindlin
        nf = 3
        bcs = [true, true, true]
        m.data[:kefunc] = rrmKe(p.d, p.E, p.nu)
        m.data[:refunc] = rrmRe(p.q)
    elseif model == :kirchhoff_conforming
        nf = 4
        bcs = [true, true, true, false]
        m.data[:kefunc] = rkpcKe(p.d, p.E, p.nu)
        m.data[:refunc] = rkpcRe(p.q)
    elseif model == :kirchhoff_nonconforming
        nf = 3
        bcs = [true, true, true]
        m.data[:kefunc] = rkpnKe(p.d, p.E, p.nu)
        m.data[:refunc] = rkpnRe(p.q)
    else
        error(p.model)
    end

    K, r = assembleKr(m, nf)
    applydirichletbcs!(m.groups[:boundarynodes], K, r, bcs)
    w = K \ r
    return m, w, nf
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
    a3d=true, w=600, h=400, title="",
    edgesvisible=false, nodesvisible=false, edgeslinewidth=2.5,
    mesh=5,
    colorrange=Makie.automatic,
    colormap=Makie.theme(:colormap),
    limits=(nothing, nothing, nothing)
)
    fig = mkfig(a3d=a3d, w=w, h=h, title=title, limits=limits)
    mplot!(
        m, makefacefunction(ww),
        faceplotzscale=zs / maximum(ww),
        faceplotmesh=mesh,
        edgesvisible=edgesvisible, edgeslinewidth=edgeslinewidth,
        nodesvisible=nodesvisible,
        color=3,
        colorrange=colorrange,
        colormap=colormap
    )
    fig
end
