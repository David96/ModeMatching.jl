ε0 = 8.8541878128e-12
c = 299792458
μ0 = 1/(ε0*c^2)

struct RectangularWaveguide <: Waveguide
    a::AbstractFloat
    b::AbstractFloat
    f::AbstractFloat
    ω::AbstractFloat
    λ::AbstractFloat
    ε::AbstractFloat
    μ::AbstractFloat
    k::AbstractFloat
    x::AbstractFloat
    y::AbstractFloat
    z::AbstractFloat
    length::AbstractFloat

    function RectangularWaveguide(x, y, z, a, b, length, μ, ε, f)
        new(a, b, f, 2π * f, c / f, ε * ε0, μ * μ0, 2π * f * sqrt(μ * ε) / c, x, y, z, length)
    end
end

function mode_from_nr(_::RectangularWaveguide, nr::Integer,
        n_TE::Integer, n_TM::Integer, max_m::Integer)
    if nr > n_TE
        # For m == 0 || n == 0 => E_TM = 0
        # Therefore n >= 1 && 1 <= m <= max_m
        N = nr - n_TE - 1
        n = div(N, max_m) + 1
        m = (N) % max_m + 1
        TMMode(m, n)
    else
        n = div(nr, max_m + 1)
        m = nr % (max_m + 1)
        TEMode(m, n)
    end
end

k_c(g::RectangularWaveguide, mode::Mode) = sqrt((mode.m * π/g.a)^2 + (mode.n * π/g.b)^2)

"Returns (lower bound, higher bound, mask) of the intersection between two waveguides"
function intersect(g1::RectangularWaveguide, g2::RectangularWaveguide)
    lb = [maximum([g1.x, g2.x]), maximum([g1.y, g2.y])]
    hb = [minimum([g1.a+g1.x, g2.a+g2.x]), minimum([g1.b+g1.y, g2.b+g2.y])]
    (lb, hb, (_, _) -> 1)
end

function contains(g::RectangularWaveguide, x, y, z)
    return g.x <= x <= g.x + g.a && g.y <= y <= g.y + g.b && g.z <= z <= g.z + g.length
end

function integral_ExHy(a1, a2, b1, b2, x1, x2, y1, y2, m1, m2, n1, n2, ai, af, bi, bf)
    if a2 * m1 != a1 * m2 && b2 * n1 != b1 * n2
        ((a1*a2*b1*b2*(a2*m1*cos(((af - x2)*m2*π)/a2)*sin(((af - x1)*m1*π)/a1) - 
                      a2*m1*cos(((ai - x2)*m2*π)/a2)*sin(((ai - x1)*m1*π)/a1) + 
                      a1*m2*((-cos(((af - x1)*m1*π)/a1))*sin(((af - x2)*m2*π)/a2) + 
                             cos(((ai - x1)*m1*π)/a1)*sin(((ai - x2)*m2*π)/a2)))*
         (b1*n2*cos(((bf - y2)*n2*π)/b2)*sin(((bf - y1)*n1*π)/b1) - 
          b1*n2*cos(((bi - y2)*n2*π)/b2)*sin(((bi - y1)*n1*π)/b1) + 
          b2*n1*((-cos(((bf - y1)*n1*π)/b1))*sin(((bf - y2)*n2*π)/b2) + 
                 cos(((bi - y1)*n1*π)/b1)*sin(((bi - y2)*n2*π)/b2))))/
         ((a2*m1 - a1*m2)*(a2*m1 + a1*m2)*(b2^2*n1^2 - b1^2*n2^2)*π^2))
    elseif a2 * m1 == a1 * m2 && b2 * n1 != b1 * n2
        if m2 != 0
            (1/(4*m2*(b2*n1 - b1*n2)*(b2*n1 + b1*n2)*π^2))*(b1*b2*(2*(af - ai)*m2*π*cos((m2*π*(x1 - x2))/a2) + 
                a2*(-sin((m2*π*(-2*af + x1 + x2))/a2) + sin((m2*π*(-2*ai + x1 + x2))/a2)))*(b1*n2*cos((n2*π*(bf - y2))/b2)*sin((n1*π*(bf - y1))/b1) - 
                b1*n2*cos((n2*π*(bi - y2))/b2)*sin((n1*π*(bi - y1))/b1) + b2*n1*((-cos((n1*π*(bf - y1))/b1))*sin((n2*π*(bf - y2))/b2) + 
                cos((n1*π*(bi - y1))/b1)*sin((n2*π*(bi - y2))/b2))))
        elseif m1 != 0
            (a1*b1*b2*(sin(((af - x1)*m1*π)/a1) - sin(((ai - x1)*m1*π)/a1))*
               (b1*n2*cos(((bf - y2)*n2*π)/b2)*sin(((bf - y1)*n1*π)/b1) - 
                b1*n2*cos(((bi - y2)*n2*π)/b2)*sin(((bi - y1)*n1*π)/b1) + 
                b2*n1*((-cos(((bf - y1)*n1*π)/b1))*sin(((bf - y2)*n2*π)/b2) + 
                  cos(((bi - y1)*n1*π)/b1)*sin(((bi - y2)*n2*π)/b2))))/(m1*(b2^2*n1^2 - b1^2*n2^2)*π^2)
        else
            (1/((b2^2*n1^2 - b1^2*n2^2)*π))*((af - ai)*b1*b2*
               (b1*n2*cos(((bf - y2)*n2*π)/b2)*sin(((bf - y1)*n1*π)/b1) - 
                b1*n2*cos(((bi - y2)*n2*π)/b2)*sin(((bi - y1)*n1*π)/b1) + 
                b2*n1*((-cos(((bf - y1)*n1*π)/b1))*sin(((bf - y2)*n2*π)/b2) + 
                  cos(((bi - y1)*n1*π)/b1)*sin(((bi - y2)*n2*π)/b2))))
        end
    elseif a2 * m1 != a1 * m2 && b2 * n1 == b1 * n2
        if n2 != 0
            -((1/(4*(a2^2*m1^2 - a1^2*m2^2)*n2*π^2))*(a1*a2*((-a2)*m1*cos((m2*π*(af - x2))/a2)*sin((m1*π*(af - x1))/a1) + a2*m1*cos((m2*π*(ai - x2))/a2)*sin((m1*π*(ai - x1))/a1) + 
                a1*m2*(cos((m1*π*(af - x1))/a1)*sin((m2*π*(af - x2))/a2) - cos((m1*π*(ai - x1))/a1)*sin((m2*π*(ai - x2))/a2)))*
                (2*(bf - bi)*n2*π*cos((n2*π*(y1 - y2))/b2) + b2*(sin((n2*π*(-2*bf + y1 + y2))/b2) - sin((n2*π*(-2*bi + y1 + y2))/b2)))))
        else
            0
        end
    else
        if m2 != 0 && n2 != 0
            (1/16)*(2*(af - ai)*cos((m2*π*(x1 - x2))/a2) + (a2*(-sin((m2*π*(-2*af + x1 + x2))/a2) + sin((m2*π*(-2*ai + x1 + x2))/a2)))/(m2*π))*
                (2*(bf - bi)*cos((n2*π*(y1 - y2))/b2) + (b2*(sin((n2*π*(-2*bf + y1 + y2))/b2) - sin((n2*π*(-2*bi + y1 + y2))/b2)))/(n2*π))
        elseif m2 == 0 && m1 != 0 && n2 != 0
            (a1*(sin((m1*π*(af - x1))/a1) - sin((m1*π*(ai - x1))/a1))*(2*(bf - bi)*n2*π*cos((n2*π*(y1 - y2))/b2) + 
                b2*(sin((n2*π*(-2*bf + y1 + y2))/b2) - sin((n2*π*(-2*bi + y1 + y2))/b2))))/(4*m1*n2*π^2)
        elseif m2 == 0 && m1 == 0 && n2 != 0
            (1/4)*(af - ai)*(2*(bf - bi)*cos((n2*π*(y1 - y2))/b2) + (b2*(sin((n2*π*(-2*bf + y1 + y2))/b2) - sin((n2*π*(-2*bi + y1 + y2))/b2)))/(n2*π))
        elseif n2 == 0
            0
        else
            throw("Forgot something")
        end
    end
end

function Ey2Ex(g::RectangularWaveguide, mode::Mode, lb, hb)
    m = mode.n; a = g.b
    n = mode.m; b = g.a
    ai = lb[2]; af = hb[2]; bi = lb[1]; bf = hb[1]
    x0 = g.y; y0 = g.x
    (m, n, a, b, x0, y0, ai, af, bi, bf)
end

"Analytical modes - for rectangular waveguides, the spatial parts of TE and TM modes are identical"
function int_ExHy(g1::RectangularWaveguide, g2::RectangularWaveguide, lb, hb, _,
        z, mode1::Mode, mode2::Mode)
    a1 = g1.a; b1 = g1.b; a2 = g2.a; b2 = g2.b
    m1 = mode1.m; n1 = mode1.n; m2 = mode2.m; n2 = mode2.n
    integral_ExHy(a1, a2, b1, b2, g1.x, g2.x, g1.y, g2.y, m1, m2, n1, n2,
                  lb[1], hb[1], lb[2], hb[2])
end
function int_EyHx(g1::RectangularWaveguide, g2::RectangularWaveguide, lb, hb, _,
        z, mode1::Mode, mode2::Mode)
    # Switching the parameters around allows for a single analytical integral
    (m1, n1, a1, b1, x1, y1, ai1, af1, bi1, bf1) = Ey2Ex(g1, mode1, lb, hb)
    (m2, n2, a2, b2, x2, y2, ai2, af2, bi2, bf2) = Ey2Ex(g2, mode2, lb, hb)
    @assert ai1 == ai2 && af1 == af2 && bi1 == bi2 && bf1 == bf2
    integral_ExHy(a1, a2, b1, b2, x1, x2, y1, y2, m1, m2, n1, n2, ai1, af1, bi1, bf1)
end

j = -1im

#=
 Split modes into orthogonal (frequency independent), frequency dependent and z dependent
 part.
 This enables us to calculate a frequency independet scalar product as the frequency factor is
 independent of the coordinates.
=#

#
# TE Modes
#
"Spatial components of E field of the TE modes in a rectangular waveguide"
E_spatial(g::RectangularWaveguide, x, y, z, mode::TEMode) = [
    cos(mode.m * π * (x - g.x) / g.a) * sin(mode.n * π * (y - g.y) / g.b),
    sin(mode.m * π * (x - g.x) / g.a) * cos(mode.n * π * (y - g.y) / g.b),
    0
]
"Frequency dependent factor of E field of the TE modes in a rectangular waveguide"
E_freq(g::RectangularWaveguide, mode::TEMode) = [
    +j * g.μ * mode.n * π / (k_c(g, mode)^2 * g.b) * g.ω,
    -j * g.μ * mode.m * π / (k_c(g, mode)^2 * g.a) * g.ω,
    0
]
"Spatial components of H field of the TE modes in a rectangular waveguide"
H_spatial(g::RectangularWaveguide, x, y, z, mode::TEMode) = [
    sin(mode.m * π * (x - g.x) / g.a) * cos(mode.n * π * (y - g.y) / g.b),
    cos(mode.m * π * (x - g.x) / g.a) * sin(mode.n * π * (y - g.y) / g.b),
    cos(mode.m * π * (x - g.x) / g.a) * cos(mode.n * π * (y - g.y) / g.b)
]
"Frequency dependent factor of H field of the TE modes in a rectangular waveguide"
H_freq(g::RectangularWaveguide, mode::TEMode) = [
    j * mode.m * π / (k_c(g, mode)^2 * g.a) * β(g, mode),
    j * mode.n * π / (k_c(g, mode)^2 * g.b) * β(g, mode),
    1
]

#
# TM Modes
#
"Spatial components of E field of the TM modes in a rectangular waveguide"
E_spatial(g::RectangularWaveguide, x, y, z, mode::TMMode) = [
    cos(mode.m * π * (x - g.x) / g.a) * sin(mode.n * π * (y - g.y) / g.b),
    sin(mode.m * π * (x - g.x) / g.a) * cos(mode.n * π * (y - g.y) / g.b),
    sin(mode.m * π * (x - g.x) / g.a) * sin(mode.n * π * (y - g.y) / g.b)
]
"Frequency dependent factor of E field of the TM modes in a rectangular waveguide"
E_freq(g::RectangularWaveguide, mode::TMMode) = [
    -j * mode.m * π / (g.a * k_c(g, mode)^2) * β(g, mode),
    -j * mode.n * π / (g.b * k_c(g, mode)^2) * β(g, mode),
    1
]
"Spatial components of H field of the TM modes in a rectangular waveguide"
H_spatial(g::RectangularWaveguide, x, y, z, mode::TMMode) = [
    sin(mode.m * π * (x - g.x) / g.a) * cos(mode.n * π * (y - g.y) / g.b),
    cos(mode.m * π * (x - g.x) / g.a) * sin(mode.n * π * (y - g.y) / g.b),
    0
]
"Frequency dependent factor of H field of the TM modes in a rectangular waveguide"
H_freq(g::RectangularWaveguide, mode::TMMode) = [
    +j * g.ε * mode.n * π / (g.b * k_c(g, mode)^2) * g.ω,
    -j * g.ε * mode.m * π / (g.a * k_c(g, mode)^2) * g.ω,
    0
]
