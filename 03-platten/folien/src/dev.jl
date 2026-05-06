include("setup.jl")

using WGLMakie
WGLMakie.activate!()
m, wHat = plate(params, 10)
we = makewe(wHat)

fig = Figure()
Axis3(fig[1, 1])
mplot!(m, we, faceplotzscale = 1000)
fig
