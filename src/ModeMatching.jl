module ModeMatching

include("waveguide.jl")
include("waveguides/rectangular_waveguide.jl")
include("waveguides/cylindrical_waveguide.jl")
include("smatrix.jl")

export Waveguide, Mode, TEMode, TMMode, RectangularWaveguide, CylindricalWaveguide,
       scalar, E_freq, E_spatial, propagation, E, H, β, f_c, k_c, int_ExHy, int_EyHx,
       besselj_prime_zero, intersect, jacobi_det, C,
       WaveguideSetup, t_r_ab, power, mode_from_nr, calc_a_i, calc_b_i, Es_right, Es_left,
       Hs_right, Hs_left, get_region

end # module ModeMatching
