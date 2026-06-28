using CairoMakie

EA = 2400 * 10000
a = 1000
h0 = 50
l0 = hypot(a, h0)

function F(u)
    h = h0 - u
    l = hypot(a, h)
    alpha = atan(h, a)
    eps = (l - l0) / l0
    N = EA * eps
    return -2 * sin(alpha) * N
end

us = -20:1:120
Fs = F.(us)
lines(us, Fs)