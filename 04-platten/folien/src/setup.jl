try
	using MMJMesh
catch e
	Pkg.add(url = "https://github.com/matthiasbaitsch/mmjmesh.git")
	using MMJMesh
end

using MMJMesh.Meshes
using MMJMesh.Utilities
using MMJMesh.Plots

using Revise
using Latexify
using GLMakie
using CairoMakie
using VarStructs
using SparseArrays
using LinearAlgebra
using FastGaussQuadrature

function _fsize(face)
	x = coordinates(face)
	p = x[:, 1]
	l1 = x[1, 2] - x[1, 1]
	l2 = x[2, 3] - x[2, 2]
	return p, l1, l2
end

function _makegmap(face) # TODO use geometry map for face
	nn(x) = (1 + x) / 2
	p, a, b = _fsize(face)
	return x -> p + [nn(x[1]) * a, nn(x[2]) * b]
end


include("fem-utils.jl")
include("plate-elements.jl")
include("fem.jl")
include("hermite4.jl")

set_theme!(theme_minimal())
update_theme!(faceplotmesh = 5)
update_theme!(edgelinewidth = 2.5)
update_theme!(colormap = :aquamarine)

params = @var Params()
params.lx = 8
params.ly = 8
params.q = 5e3
params.nu = 0.0
params.d = 0.2
params.E = 31000e6

## Gausspoint interpolation
xg, _ = gausslegendre(2)
VG = [xg[[1 2 2 1]]; xg[[1 1 2 2]]]
L4g = makenodalbasis(
	VG,
	mmonomials(2, 1, QHat),
	p -> [ValueAtLF(p)],
)

plot2d!() = CairoMakie.activate!()
plot3d!() = GLMakie.activate!()

function interpolateg(f, np = 2)
	if np == 2
		return sum(L4g .* f.([collect(p) for p ∈ eachcol(VG)]))
	else
		return MPolynomial([0; 0;;], [f(0, 0)], QHat)
	end
end

function valuerange(p)
	values = p.plots[1].kw[:color]
	min, max = extrema(values)
	return min, max
end

function makewe(wHat)
	return face -> begin
		idxs = dofs(nodeindices(face), 4)
		_, a, b = _fsize(face)
		t = repeat([1, a / 2, b / 2, a * b / 4], 4)
		return sum(wHat[idxs] .* t .* H4) # TODO make dot work
	end
end

ei(n, i) = [j == i ? 1 : 0 for j ∈ 1:n]

function nodalresult(m, result)
	sr = zeros(nnodes(m))
	VV = [[-1, -1], [1, -1], [1, 1], [-1, 1]]
	for f ∈ faces(m)
		sr[nodeindices(f)] .+= m.data[:post](f, result).(VV)
	end
	for (i, n) ∈ enumerate(nodes(m))
		sr[i] /= nfaces(n)
	end
	return sr
end

function postprocessor(params, wHat)
	return (face, name) -> begin

		# Element size and differential operators
		_, a, b = _fsize(face)
		∂X(we) = (2 / a) * ∂x(we)
		∂Y(we) = (2 / b) * ∂y(we)
		∂XX(we) = (2 / a)^2 * ∂xx(we)
		∂YY(we) = (2 / b)^2 * ∂yy(we)
		∂XY(we) = (2 / a) * (2 / b) * ∂xy(we)

		# Element displacement function
		idxs = dofs(nodeindices(face), 4)
		t = repeat([1, a / 2, b / 2, a * b / 4], 4)
		we = sum(wHat[idxs] .* t .* H4) # TODO make dot work

		# Quick return
		name == :w && return we

		# Derivatives
		wxx = ∂XX(we)
		wyy = ∂YY(we)
		wxy = ∂XY(we)
		Δw = wxx + wyy

		# Return
		name == :wx && return ∂X(we)
		name == :wy && return ∂Y(we)
		name == :wxx && return wxx
		name == :wyy && return wyy
		name == :wxy && return wxy
		name == :Δw && return Δw

		# Plate properties
		h = params.d
		E = params.E
		nu = params.nu
		D = E * h^3 / (12 * (1 - nu^2))

		# Section forces (Altenbach et al. p176)
		mx = -1e-3 * D * (wxx + nu * wyy)
		my = -1e-3 * D * (nu * wxx + wyy)
		mxy = -1e-3 * D * (1 - nu) * wxy
		qx = -1e-3 * D * ∂X(Δw)
		qy = -1e-3 * D * ∂Y(Δw)

		# From equilibrium
		qxe = ∂X(mx) + ∂Y(mxy)
		qye = ∂X(mxy) + ∂Y(my)

		# Return
		name == :mx && return mx
		name == :my && return my
		name == :mxy && return mxy
		name == :qx && return qx
		name == :qy && return qy
		name == :qxe && return qxe
		name == :qye && return qye
		name == :mxg && return interpolateg(mx, 2)
		name == :myg && return interpolateg(my, 2)
		name == :mxyg && return interpolateg(mxy, 2)
		name == :qxg && return interpolateg(qx, 2)
		name == :qyg && return interpolateg(qy, 2)
		name == :qxgg && return ∂X(interpolateg(mx, 2)) + ∂Y(interpolateg(mxy, 2))
		name == :qygg && return ∂X(interpolateg(mxy, 2)) + ∂Y(interpolateg(my, 2))

		# Unknown label
		error("Unkown function: ", name)
	end
end

function plate(p, nx, ny = nx, model = :kirchhoff_conforming)
	nf = 0
	bcs = []
	m = makemeshonrectangle(p.lx, p.ly, nx, ny)

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
	wHat = K \ r
	m.data[:post] = postprocessor(params, wHat)

	return m, wHat, nf
end

function mkfig(
	; a3d = true, w = 600, h = 400, title = "",
	limits = (nothing, nothing, nothing),
)
	fig = Figure(size = (w, h))
	if a3d
		plot3d!()
		ax = Axis3(
			fig[1, 1],
			aspect = :data,
			title = title,
			viewmode = :stretch,
			perspectiveness = 0.2,
			limits = limits,
			protrusions = 0,
		)
	else
		plot2d!()
		ax = Axis(fig[1, 1], aspect = DataAspect(), title = title)
	end
	hidedecorations!(ax)
	hidespines!(ax)
	return fig
end

function plotw(
	m, wHat;
	zs = 2,
	a3d = true, w = 600, h = 400, title = "",
	edgesvisible = false, nodesvisible = false, edgelinewidth = 2.5,
	mesh = 5,
	colorrange = Makie.automatic,
	colormap = Makie.theme(:colormap),
	limits = (nothing, nothing, nothing),
)
	fig = mkfig(a3d = a3d, w = w, h = h, title = title, limits = limits)
	mplot!(
		m, makewe(wHat),
		faceplotzscale = zs / maximum(wHat),
		faceplotmesh = mesh,
		edgesvisible = edgesvisible, edgelinewidth = edgelinewidth,
		nodesvisible = nodesvisible,
		color = 3,
		colorrange = colorrange,
		colormap = colormap,
	)
	fig
end

function maketitle(p, title)
	min, max = string.(round.(valuerange(p), digits = 2))
	return title * " | min: " * min * " | max: " * max
end

function plotr(m, result, title, cr, npoints = 15; nodal = false, a3d)
	if false
		fig = Figure(size = (1200, 800))
		if a3d
			ax = Axis3(fig[1, 1])
			zscale = 1
		else
			ax = Axis(fig[1, 1], aspect = DataAspect())
			zscale = 0
		end

		warpnodes = nothing

		if nodal
			r = nodalresult(m, result)
			if a3d
				warpnodes = r
			end
		else
			r = result
		end

		p = mplot!(
			m, r, edgesvisible = a3d,
			nodewarp = warpnodes,
			faceplotzscale = zscale, faceplotnpoints = npoints,
			colorrange = cr,
		)

		if !a3d
			mplot!(m, edgesvisible = false, facesvisible = false)
		end

		Colorbar(fig[1, 2], p)
		hidespines!(ax)
		hidedecorations!(ax)
		ax.title = maketitle(p, title)

		return fig
	else
		return L"Zeitweise außer Betrieb"
	end
end
