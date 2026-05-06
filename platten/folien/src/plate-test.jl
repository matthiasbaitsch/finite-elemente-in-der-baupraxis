using Test
using WGLMakie

include("setup.jl")

include("fem-utils.jl")

WGLMakie.activate!()


## Analysis of one plate
# params.nu = 0.3
m, wHat, nf = plate(params, 10, 10, :reissner_mindlin)
mplot(m, edgesvisible = true) |> mconf()
update_theme!(colormap = :acton)
mplot(m, edgesvisible = true, wHat[1:nf:end]) |> mconf()
update_theme!(colormap = :redblue)
mplot(m, edgesvisible = true, wHat[2:nf:end]) |> mconf()
mplot(m, edgesvisible = true, wHat[3:nf:end]) |> mconf()
#mplot(m, edgesvisible=true, w[4:nf:end]) |> mconf()

maximum(wHat[1:nf:end])

## Convergence study for quadratic plate - Czerny

# Compute
l = params.lx
params.h = 0.2
params.nu = 0

hh = []
nn = []
wwc = []
wwn = []
for n ∈ 4:2:30
	print("$n ")
	mn, wn, nf = plate(params, n, n, :kirchhoff_conforming)
	push!(hh, l / n)
	push!(nn, 4 * nnodes(mn))
	push!(wwc, maximum(abs.(wn[1:nf:end])))
	mn, wn, nf = plate(params, n, n, :kirchhoff_nonconforming)
	push!(wwn, maximum(abs.(wn[1:nf:end])))
end
println()

# Comparison to Czerny tables
w_fe = wwc[end]
w_czerny = params.q * l^4 / (params.E * params.d^3) * 0.0152
100 * abs(w_fe - w_czerny) / w_czerny

# Plot
f = Figure()
Axis(f[1, 1])
p3 = lines!(hh, 1000 * w_czerny * ones(length(hh)), color = :gray)
p1 = scatterlines!(hh, 1000 * wwc)
p2 = scatterlines!(hh, 1000 * wwn)
Legend(f[1, 2], [p1, p2, p3], ["Compatible", "Incompatible", "Czerny"])
f


## Convergence study for quadratic plate - comparison with Reissner Mindlin

# Compute
l = params.lx
params.d = 0.5
params.nu = 0.27

hh = []
nn = []
wwc = []
wwn = []
for n ∈ 4:2:30
	print("$n ")
	mn, wn, nf = plate(params, n, n, :kirchhoff_conforming)
	push!(hh, l / n)
	push!(nn, 4 * nnodes(mn))
	push!(wwc, maximum(abs.(wn[1:nf:end])))
	mn, wn, nf = plate(params, n, n, :kirchhoff_nonconforming)
	push!(wwn, maximum(abs.(wn[1:nf:end])))
end
println()

hhrm = []
nnrm = []
wwrm = []
for n ∈ 2 .^ collect((3:7))
	print("$n ")
	mn, wn, nf = plate(params, n, n, :reissner_mindlin)
	push!(hhrm, l / n)
	push!(nnrm, 4 * nnodes(mn))
	push!(wwrm, maximum(abs.(wn[1:nf:end])))
end
println()

# Plot
f = Figure()
Axis(f[1, 1])
p1 = scatterlines!(hh, 1000 * wwc)
p2 = scatterlines!(hh, 1000 * wwn)
p3 = scatterlines!(hhrm, 1000 * wwrm)
Legend(f[1, 2], [p1, p2, p3], ["Compatible", "Incompatible", "Reissner-Mindlin"])
f

## Plots

##

include("setup.jl")

using WGLMakie
WGLMakie.activate!()
set_theme!(theme_minimal())

update_theme!(
	colormap = :redblue,
	color = 3,
	faceplotzscale = 1,
	faceplotnpoints = 15,
	edgesvisible = true,
	featureedgelinewidth = 2.5,
)

params.E = 34000
params.nu = 0.2

m, wHat = plate(params, 20)
we = makewe(wHat)
post = m.data[:post]

println("Max w: ", 1000 * maximum(wHat[1:4:end]))

# Plot w
mplot(m, we)


##
# Test two versions of shear force
using Test
ff = face(m, 33)
qx = post(ff, :qx)
qxe = post(ff, :qxe)
@test qx.p.exponents == qxe.p.exponents
@test qx.p.coefficients ≈ qxe.p.coefficients


## Result to plot
result = :qxg

## Plot 3D
fig = Figure()
Axis3(fig[1, 1])
p = mplot!(m, result, faceplotzscale = 1)
Colorbar(fig[1, 2], p, label = string(result))
fig

## Plot 2D 
fig = Figure()
Axis(fig[1, 1], aspect = DataAspect())
(p = mplot!(m, result, faceplotzscale = 0, edgesvisible = false)) |> valuerange
mplot!(m, facesvisible = false, edgesvisible = false)
Colorbar(fig[1, 2], p, label = string(result))
fig

## Plot 2D Smooth result
fig = Figure()
Axis(fig[1, 1], aspect = DataAspect())
mplot!(m, nodalresult(m, result), edgesvisible = false) |> valuerange
Colorbar(fig[1, 2], p, label = string(result))
fig

## Plot one element
fig = Figure()
Axis3(fig[1, 1])
fplot3d!(qx, gmap = _makegmap(ff), npoints = 40)
fig

## Cerny
p = 1e-3 * params.q
lx = params.lx

mxerm = -p * lx^2 / 19.4
mxm = p * lx^2 / 56.8
qxerm = p * lx / 2.24

println("mxerm: $mxerm")
println("  mxm: $mxm")
println("qxerm: $qxerm")


