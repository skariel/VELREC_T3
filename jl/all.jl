using Logging
if !(isdefined(:ICORE) && ICORE==true)
    using Winston
    using NearestNeighbors
end

Logging.configure(level=INFO)

include("utils.jl")
include("realization.jl")
include("grid_grav.jl")
include("cic.jl")
include("2lpt.jl")
include("cosmology.jl")
include("kd_tools.jl")
include("correlations.jl")
include("dark_sky_halos_1600.jl")
include("simulation.jl")
include("optimization.jl")
