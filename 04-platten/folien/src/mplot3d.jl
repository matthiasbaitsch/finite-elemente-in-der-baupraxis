using MakieCore

using MMJMesh.Plots: sample2d, sample2dlines

function _fsize(face)
    x = coordinates(face)
    p = x[:, 1]
    l1 = x[1, 2] - x[1, 1]
    l2 = x[2, 3] - x[2, 2]
    return p, l1, l2
end

function _makegmap(face) # TODO move to mmjmesh
    nn(x) = (1 + x) / 2
    p, a, b = _fsize(face)
    return x -> p + [nn(x[1]) * a, nn(x[2] * b)]
end


function _getcolor(x::Matrix, color, zscale)
    if typeof(color) == Int && 1 <= color <= 3
        return x[color, :] / zscale
    end
    return color
end

function _collectlines(cl)
    l1 = Float32[]
    l2 = Float32[]
    l3 = Float32[]
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

function _collectfaces(cf)
    xx = [Float32[], Float32[], Float32[]]
    tt = [Int[], Int[], Int[]]
    for c in cf
        xf, tf = c
        pos = length(xx[1])
        for i = 1:3
            append!(xx[i], xf[i, :])
            append!(tt[i], pos .+ tf[:, i])
        end
    end
    return stack(xx, dims=1), stack(tt)
end

function _sample(m, mf, npoints, zscale, mesh)
    cf = []
    cl1 = []
    cl2 = []
    for face = faces(m)
        f = mf(face)
        gmap = _makegmap(face)
        push!(cf, sample2d(f, domain=QHat, npoints=2 * npoints, gmap=gmap, zscale=zscale))
        push!(cl1, sample2dlines(f, domain=QHat, npoints=npoints, mesh=0, gmap=gmap, zscale=zscale))
        push!(cl2, sample2dlines(f, domain=QHat, npoints=npoints, mesh=mesh, gmap=gmap, zscale=zscale))
    end
    return cf, cl1, cl2
end

MakieCore.@recipe(MPlot3D, functions) do scene
    attr = Attributes(
        npoints=30,
        color=3,
        meshcolor=:black,
        mesh=5,
        edgesvisible=true,
        edgeslinewidth=2,
        colorrange=MakieCore.automatic,
        colormap=MakieCore.theme(scene, :colormap),
        zscale=1
    )
    MakieCore.generic_plot_attributes!(attr)
    MakieCore.colormap_attributes!(attr, MakieCore.theme(scene, :colormap))
    return attr
end

function MakieCore.plot!(plot::MPlot3D)
    attributes = plot.attributes
    npoints = attributes.npoints[]
    color = attributes.color[]
    meshcolor = attributes.meshcolor[]
    mesh = attributes.mesh[]
    edgesvisible = attributes.edgesvisible[]
    edgeslinewidth = attributes.edgeslinewidth[]
    colorrange = attributes.colorrange[]
    colormap = attributes.colormap[]
    zscale = attributes.zscale[]

    m = plot.args[1][]
    mf = plot.args[2][]
    cf, cl1, cl2 = _sample(m, mf, npoints, zscale, mesh)
    x, t = _collectfaces(cf)

    mesh!(plot, x, t, color=_getcolor(x, color, zscale), colormap=colormap, colorrange=colorrange)
    lines!(plot, _collectlines(cl2)..., color=meshcolor, linewidth=1.25)
    if edgesvisible
        lines!(plot, _collectlines(cl1)..., color=meshcolor, linewidth=edgeslinewidth)
    end

    return plot
end
