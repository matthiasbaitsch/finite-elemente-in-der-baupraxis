using Roots
using LinearAlgebra

# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

function arc_through_points(p_start, p_mid, p_end)
    # Center and radius from Gram matrix - with a little help from claude.ai
    D = hcat(p_mid - p_start, p_end - p_start)
    G = D' * D
    center = p_start + D * (G \ (diag(G) / 2))
    r = norm(p_start - center)

    # Basis vectors in plane of arc
    e1 = normalize(p_mid - p_start)
    e3 = normalize(cross(e1, p_end - p_start))
    e2 = cross(e3, e1)

    # Angles of start and endpoint
    theta_p(p) = atan(dot(p - center, e2), dot(p - center, e1))
    theta1 = theta_p(p_start)
    theta3 = theta_p(p_end)

    # Angle from parameter 0 <= t <= 1
    theta_t(t) = theta1 + t * (theta3 - theta1)

    return t -> center + r * (cos(theta_t(t)) * e1 + sin(theta_t(t)) * e2)
end

function plane_from_angle(alpha)
    p = [0, 0, 0]
    n = [-sin(deg2rad(alpha)), cos(deg2rad(alpha)), 0]
    return x -> dot(x - p, n)
end

# -------------------------------------------------------------------------------------------------
# Actual calculation
# -------------------------------------------------------------------------------------------------

curve = arc_through_points([108.5, -50, 0], [98, 0, 30], [108.5, 50, 0])

for alpha = -16:3.2:16
    plane = plane_from_angle(alpha)
    equation = t -> plane(curve(t))
    parameter = find_zero(equation, 0.5)
    print(round(100 * parameter, digits=2), "%, ")
end

