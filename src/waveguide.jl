using HCubature

"""
    Waveguide API:
    Every waveguide implementation is expected to provide:
    - `E_spatial(g::WaveguideImpl, x, y, z, mode::Mode)::SVector`
    - `E_freq(g::WaveguideImpl, mode::Mode)::SVector`
    - `H_spatial(g::WaveguideImpl, x, y, z, mode::Mode)::SVector`
    - `E_freq(g::WaveguideImpl, mode::Mode)::SVector`
    returning three dimensional vectors of the spatial and frequency dependent parts
    of the E and H fields of the Eigenmodes.
    - `k_c(g::WaveguideImpl, mode::Mode) <: AbstractFloat` returning the cutoff wavevector
    - `intersect(g1::WaveguideImpl, g2::WaveguideImpl)::Tuple`
    returning (lowerbound, higherbound, mask) of the two dimensional intersection
    of the two waveguides.
    - `contains(g::WaveuideImpl, x, y, z)::Bool` returning true if the point lies within the
    waveguide, false otherwise.
    - `mode_from_nr(g::WaveguideImpl, nr::T, n_TE, n_TM, max_m)::Mode where T<:Integer`
    has to provide a unique mapping of a nr to a waveguide mode with 1 <= nr <= n_TE + n_TM.

    - `int_ExHy(g1::WaveguideImpl, g2::WaveuideImpl, mode1::Mode, mode2::Mode)` and
    `int_EyHx(g1::WaveguideImpl, g2::WaveuideImpl, mode1::Mode, mode2::Mode)`
    can be used to provide analytical solutions to the overlap integrals.
    - `jacobi_det(g::WaveguideImpl, x, y, z)` can be used in case non-cartesian
    coordinates are used.

    Each Waveguide struct is expected to provide:
    - `k <: AbstractFloat`: The wavevector.
    - `A <: AbstractFloat`: The cross section of the waveguide.
"""
abstract type Waveguide end
abstract type Mode end

struct TEMode{T<:Integer} <: Mode
    m::T
    n::T
end
struct TMMode{T<:Integer} <: Mode
    m::T
    n::T
end
struct TEMMode{T<:Integer} <: Mode
    l::T
    p::T
end

@enum Direction fwd=1 bck=-1

"In case of non-cartesian coordinate system, supply a jacobi determinant"
jacobi_det(_::Waveguide, _, _, _) = 1

function integral(lb, hb, mask, f)
    #println("WARNING: numerical integral")
    (I, _) = hcubature(r -> f(r[1], r[2]) * mask(r[1], r[2]),
                       lb, hb;
                       rtol=1e-4,
                       #atol=1e-12,
                       maxevals=Integer(1e5))
    I
end
function int_ExHy(g1::Waveguide, g2::Waveguide, lb, hb, mask, z, mode1::Mode, mode2::Mode)
    E1(x, y) = E_spatial(g1, x, y, z, mode1, fwd)
    H2(x, y) = H_spatial(g2, x, y, z, mode2, fwd)
    integral(lb, hb, mask, (x, y) -> E1(x, y)[1] * H2(x, y)[2] * jacobi_det(g1, x, y, z))
end
function int_EyHx(g1::Waveguide, g2::Waveguide, lb, hb, mask, z, mode1::Mode, mode2::Mode)
    E1(x, y) = E_spatial(g1, x, y, z, mode1, fwd)
    H2(x, y) = H_spatial(g2, x, y, z, mode2, fwd)
    integral(lb, hb, mask, (x, y) -> E1(x, y)[2] * H2(x, y)[1] * jacobi_det(g1, x, y, z))
end
const integrals = Dict{UInt64, Tuple{ComplexF64, ComplexF64}}()
const Slock = ReentrantLock()
scalar(g1, g2, z, mode1, mode2; norm=true) = scalar(g1, g2, z, intersect(g1, g2)..., mode1, mode2; norm)
function scalar(g1::Waveguide, g2::Waveguide, z, lb, hb, mask, mode1::Mode, mode2::Mode; norm=true)
    # Remove if in doubt of orthonormal system!
    if norm# && false
        if g1 == g2 && mode1 == mode2
            return Complex(1)
        elseif g1 == g2
            return Complex(0)
        end
    end
    h = hash((integral_deps(g1, z), integral_deps(g2, z)), hash((mode1, mode2)))
    if !(h in keys(integrals))
        lock(Slock)
        if !(h in keys(integrals))
            I1 = int_ExHy(g1, g2, lb, hb, mask, z, mode1, mode2)
            I2 = int_EyHx(g1, g2, lb, hb, mask, z, mode1, mode2)
            integrals[h] = (I1, I2)
        end
        unlock(Slock)
    end
    I1, I2 = integrals[h]

    E1_f = E_freq(g1, mode1, fwd)
    H2_f = H_freq(g2, mode2, fwd)
    0.5 * (norm ? (C(g1, mode1) * C(g2, mode2)) : 1) *
            (E1_f[1] * H2_f[2] * I1 - E1_f[2] * H2_f[1] * I2)
end
const Cs = Dict()
const lockCs = ReentrantLock()
function C(g::Waveguide, mode::Mode)
    h = hash(g, hash(mode))
    lock(lockCs)
    if !(h in keys(Cs))
        Cm = scalar(g, g, 0, mode, mode; norm=false)
        Cs[h] = 1 / sqrt(Complex(Cm))
    end
    unlock(lockCs)
    Cs[h]
end

β(g::Waveguide, mode::Mode, dir::Direction) = Int(dir) * conj(sqrt(Complex(g.k^2 - k_c(g, mode)^2)))
f_c(g::Waveguide, mode::Mode) = k_c(g, mode) / (2π * sqrt(g.μ * g.ε))
is_propagating(g::Waveguide, mode::Mode) = isreal(β(g, mode, fwd))

function propagation(g::Waveguide, mode::Mode, dir::Direction, z)
    z_rel = dir == fwd ? g.z : (g.z + g.length)
    exp(-1im * β(g, mode, dir) * (z - z_rel))
end

#=
 General mode functions combining normalisation, frequency dependence, orthogonal component and
 z dependence
=#
"""
    E(g::Waveguide, x, y, z, mode::Mode[; norm=true])

    Return the vector of the E field of a certain `mode` at position `x`, `y`, `z`.
    `norm` specifies whether the modes should be normalized in regards to `scalar` (default: true)
"""
function E(g::Waveguide, x, y, z, mode::Mode, dir::Direction; norm=true)
    (norm ? C(g, mode) : 1) .*
        E_freq(g, mode, dir) .* E_spatial(g, x, y, z, mode, dir) .*
        propagation(g, mode, dir, z)
end

"""
    H(g::Waveguide, x, y, z, mode::Mode[; norm=true])

    Return the vector of the H field of a certain `mode` at position `x`, `y`, `z`.
    `norm` specifies whether the modes should be normalized in regards to `scalar` (default: true)
"""
function H(g::Waveguide, x, y, z, mode::Mode, dir::Direction; norm=true)
    (norm ? C(g, mode) : 1) .* H_freq(g, mode, dir) .*
        H_spatial(g, x, y, z, mode, dir) .* 
        propagation(g, mode, dir, z)
end
