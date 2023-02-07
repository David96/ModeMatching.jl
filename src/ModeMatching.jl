module ModeMatching

include("waveguide.jl")
include("waveguides/rectangular_waveguide.jl")

export Waveguide, Mode, TEMode, TMMode, RectangularWaveguide,
       scalar, E_freq, E_spatial, z_dep, E, H

end # module ModeMatching
