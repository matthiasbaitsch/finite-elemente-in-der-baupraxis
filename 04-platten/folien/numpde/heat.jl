function heatKe(λ)
    function keFunc(e)
        x = coordinates(e)
        u = x[:, 2] - x[:, 1]
        v = x[:, 3] - x[:, 2]
        w = x[:, 1] - x[:, 3]
        A = 0.5 * (u[1] * v[2] - u[2] * v[1])
        B = 1 / (2 * A) * [
            -v[2] -w[2] -u[2]
             v[1]  w[1]  u[1]
        ]        
        return λ * A * B' * B
    end
    return keFunc
end

function heatRe(w)
    function reFunc(e)
        x = coordinates(e)
        u = x[:, 2] - x[:, 1]
        v = x[:, 3] - x[:, 2]
        A = 0.5 * (u[1] * v[2] - u[2] * v[1])
        return w * A / 3 * [1; 1; 1]
    end
    return reFunc
end

function robinRe(h, ts)
    function reFunc(e)
        l = length(e)
        return h * ts * l / 2 * [1; 1]
    end
    return reFunc
end

function robinKe(h)
    function keFunc(e)
        l = length(e)
        return h * l / 6 * [2 1; 1 2]
        return [2 1; 1 2]
    end
    return keFunc
end

