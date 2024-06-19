using WGLMakie
using CairoMakie


using MMJMesh.Plots: sample2d, sample2dlines



function collectlines(cl)
    l1 = Float64[]
    l2 = Float64[]
    l3 = Float64[]
    for c in cl
        append!(l1, c[1])
        push!(l1, NaN)
        append!(l2, c[2])
        push!(l2, NaN)
        append!(l3, c[3])
        push!(l3, NaN)
    end
    return l1, l2, l3
end

function collectfaces(cf)
    xx = [[], [], []]
    tt = [[], [], []]
    for c in cf
        xf, tf = c
        pos = length(xx[1])
        for i = 1:3
            append!(xx[i], xf[i, :])
            append!(tt[i], pos .+ tf[:, i])
        end
    end
    x = Float64.(stack(xx))
    t = Int.(stack(tt))
    return x, t
end

function plotmeshsolution2(m, mf; npoints=20, zscale=1, mesh=4)

    cf = []
    cl1 = []
    cl2 = []

    for face = faces(m)
        f = mf(face)
        gmap = makegmap(face)
        push!(cf, sample2d(f, domain=QHat, npoints=2npoints, gmap=gmap, zscale=zscale))
        push!(cl1, sample2dlines(f, domain=QHat, npoints=npoints, mesh=0, gmap=gmap, zscale=zscale))
        push!(cl2, sample2dlines(f, domain=QHat, npoints=npoints, mesh=mesh, gmap=gmap, zscale=zscale))
    end

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data)
    hidedecorations!(ax)
    lines!(collectlines(cl1)..., color=:black, linewidth=2)
    lines!(collectlines(cl2)..., color=:black, linewidth=1)

    x, t = collectfaces(cf)
    mesh!(x, t, color=x[:, 3] / zscale)

    return fig
end
