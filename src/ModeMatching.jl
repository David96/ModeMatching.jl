module ModeMatching

include("waveguide.jl")
include("waveguides/rectangular_waveguide.jl")
include("waveguides/cylindrical_waveguide.jl")

export Waveguide, Mode, TEMode, TMMode, RectangularWaveguide, CylindricalWaveguide,
       scalar, E_freq, E_spatial, propagation, E, H, β, f_c, k_c, int_ExHy, int_EyHx,
       besselj_prime_zero, intersect, jacobi_det

end # module ModeMatching
