using HCubature

"""
    Waveguide API:
    Every waveguide implementation is expected to provide:
    - `E_spatial(g::WaveguideImpl, x, y, z, mode::Mode)::Vector`
    - `E_freq(g::WaveguideImpl, mode::Mode)::Vector`
    - `H_spatial(g::WaveguideImpl, x, y, z, mode::Mode)::Vector`
    - `E_freq(g::WaveguideImpl, mode::Mode)::Vector`
    returning three dimensional vectors of the spatial and frequency dependent parts
    of the E and H fields of the Eigenmodes.
    - `k_c(g::WaveguideImpl, mode::Mode)::AbstractFloat` returning the cutoff wavevector
    - `intersect(g1::WaveguideImpl, g2::WaveguideImpl)::Tuple`
    returning (lowerbound, higherbound, mask) of the two dimensional intersection
    of the two waveguides.
    - `contains(g::WaveuideImpl, x, y, z)::Bool` returning true if the point lies within the
    waveguide, false otherwise.
    - `mode_from_nr(g::WaveguideImpl, nr::Integer, n_TE, n_TM, max_m)::Mode`
    has to provide a unique mapping of a nr to a waveguide mode with 1 <= nr <= n_TE + n_TM.

    - `int_ExHy(g1::WaveguideImpl, g2::WaveuideImpl, mode1::Mode, mode2::Mode)` and
    `int_EyHx(g1::WaveguideImpl, g2::WaveuideImpl, mode1::Mode, mode2::Mode)`
    can be used to provide analytical solutions to the overlap integrals.
    - `jacobi_det(g::WaveguideImpl, x, y, z)::Number` can be used in case non-cartesian
    coordinates are used.

    Each Waveguide struct is expected to provide:
    - `k::AbstractFloat`: The wavevector.
    - `A::AbstractFloat`: The cross section of the waveguide.
"""
abstract type Waveguide end
abstract type Mode end

struct TEMode <: Mode
    m::Integer
    n::Integer
end
struct TMMode <: Mode
    m::Integer
    n::Integer
end

"In case of non-cartesian coordinate system, supply a jacobi determinant"
jacobi_det(_::Waveguide, _, _, _) = 1

function integral(g1::Waveguide, g2::Waveguide, f)
    (lb, hb, mask) = intersect(g1, g2)
    println("WARNING: numerical integral")
    (I, _) = hcubature(r -> f(r[1], r[2]) * mask(r[1], r[2]),
                       lb, hb;
                       rtol=1e-4,
                       #atol=1e-12,
                       maxevals=Integer(1e6))
    I
end
function int_ExHy(g1::Waveguide, g2::Waveguide, z, mode1::Mode, mode2::Mode)
    E1(x, y) = E_spatial(g1, x, y, z, mode1)
    H2(x, y) = H_spatial(g2, x, y, z, mode2)
    integral(g1, g2, (x, y) -> E1(x, y)[1] * H2(x, y)[2] * jacobi_det(g1, x, y, z))
end
function int_EyHx(g1::Waveguide, g2::Waveguide, z, mode1::Mode, mode2::Mode)
    E1(x, y) = E_spatial(g1, x, y, z, mode1)
    H2(x, y) = H_spatial(g2, x, y, z, mode2)
    integral(g1, g2, (x, y) -> E1(x, y)[2] * H2(x, y)[1] * jacobi_det(g1, x, y, z))
end
function scalar(g1::Waveguide, g2::Waveguide, z, mode1::Mode, mode2::Mode; norm=true)
    # Remove if in doubt of orthonormal system!
    if norm && false
        if g1 == g2 && mode1 == mode2
            return Complex(1)
        elseif g1 == g2
            return Complex(0)
        end
    end
    I1 = int_ExHy(g1, g2, z, mode1, mode2)
    I2 = int_EyHx(g1, g2, z, mode1, mode2)

    E1_f = E_freq(g1, mode1)
    H2_f = H_freq(g2, mode2)
    ((norm ? C(g1, mode1) * C(g2, mode2) : 1) *
     (E1_f[1] * H2_f[2] * I1 - E1_f[2] * H2_f[1] * I2))
end
Cs = Dict()
lockCs = ReentrantLock()
function C(g::Waveguide, mode::Mode)
    hash = "$g;$mode"
    lock(lockCs)
    if !(hash in keys(Cs))
        Cm = scalar(g, g, 0, mode, mode; norm=false)
        Cs[hash] = 1 / sqrt(Complex(Cm))
    end
    unlock(lockCs)
    Cs[hash]
end

β(g::Waveguide, mode::Mode) = -sqrt(Complex(g.k^2 - k_c(g, mode)^2))
f_c(g::Waveguide, mode::Mode) = k_c(g, mode) / (2π * sqrt(g.μ * g.ε))

propagation(g::Waveguide, mode::Mode, z) = exp(-1im * β(g, mode) * z)

#=
 General mode functions combining normalisation, frequency dependence, orthogonal component and
 z dependence
=#
"""
    E(g::Waveguide, x, y, z, mode::Mode[; norm=true])

    Return the vector of the E field of a certain `mode` at position `x`, `y`, `z`.
    `norm` specifies whether the modes should be normalized in regards to `scalar` (default: true)
"""
function E(g::Waveguide, x, y, z, mode::Mode; norm=true)
    (norm ? C(g, mode) : 1) *
        E_freq(g, mode) .* E_spatial(g, x, y, z, mode) .*
        propagation(g, mode, z)
end

"""
    H(g::Waveguide, x, y, z, mode::Mode[; norm=true])

    Return the vector of the H field of a certain `mode` at position `x`, `y`, `z`.
    `norm` specifies whether the modes should be normalized in regards to `scalar` (default: true)
"""
function H(g::Waveguide, x, y, z, mode::Mode; norm=true)
    (norm ? C(g, mode) : 1) * H_freq(g, mode) .*
        H_spatial(g, x, y, z, mode) .* 
        propagation(g, mode, z)
end
