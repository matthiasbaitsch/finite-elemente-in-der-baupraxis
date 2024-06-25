include("setup.jl")

using WGLMakie
WGLMakie.activate!()


## Analysis of one plate
# params.nu = 0.3
m, w, nf = plate(10, 10, 50, 50, params, :reissner_mindlin)
mplot(m, edgesvisible=true) |> mconf()
update_theme!(colormap=:acton)
mplot(m, edgesvisible=true, w[1:nf:end]) |> mconf()
update_theme!(colormap=:redblue)
mplot(m, edgesvisible=true, w[2:nf:end]) |> mconf()
mplot(m, edgesvisible=true, w[3:nf:end]) |> mconf()

maximum(w[1:nf:end])

#mplot(m, edgesvisible=true, w[4:nf:end]) |> mconf()


## Convergence study for quadratic plate - Czerny

# Compute
l = 10
params.h = 0.2
params.nu = 0

hh = []
nn = []
wwc = []
wwn = []
for n = 4:2:30
    print("$n ")
    mn, wn, nf = plate(l, l, n, n, params, :kirchhoff_conforming)
    push!(hh, l / n)
    push!(nn, 4 * nnodes(mn))
    push!(wwc, maximum(abs.(wn[1:nf:end])))
    mn, wn, nf = plate(l, l, n, n, params, :kirchhoff_nonconforming)
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
p3 = lines!(hh, 1000 * w_czerny * ones(length(hh)), color=:gray)
p1 = scatterlines!(hh, 1000 * wwc)
p2 = scatterlines!(hh, 1000 * wwn)
Legend(f[1, 2], [p1, p2, p3], ["Compatible", "Incompatible", "Czerny"])
f


## Convergence study for quadratic plate - With Reissner Mindlin

# Compute
l = 10
params.d = 0.5
params.nu = 0.27

hh = []
nn = []
wwc = []
wwn = []
for n = 4:2:30
    print("$n ")
    mn, wn, nf = plate(l, l, n, n, params, :kirchhoff_conforming)
    push!(hh, l / n)
    push!(nn, 4 * nnodes(mn))
    push!(wwc, maximum(abs.(wn[1:nf:end])))
    mn, wn, nf = plate(l, l, n, n, params, :kirchhoff_nonconforming)
    push!(wwn, maximum(abs.(wn[1:nf:end])))
end
println()

hhrm = []
nnrm = []
wwrm = []
for n = 2 .^ collect((3:7))
    print("$n ")
    mn, wn, nf = plate(l, l, n, n, params, :reissner_mindlin)
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


# update_theme!(colormap=:acton)
#update_theme!(faceplotmesh=5)
#update_theme!(edgelinewidth=2.5)
# update_theme!(faceplotmeshcolor=:green)
# update_theme!(edgecolor=:black)


m, w = plate(10, 5, 12, 6, params)
fig = Figure()
Axis3(fig[1, 1], aspect=:data)
mplot!(
    m, makefacefunction(w),
    color=3,
    faceplotzscale=2 / maximum(w)
)
fig
