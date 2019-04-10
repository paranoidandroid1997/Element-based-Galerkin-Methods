#---------------------------------------------------------------------#
#This module stores everything needed to build different qudrature points
#in 1D.
#Written by F.X. Giraldo on 4/10/2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#---------------------------------------------------------------------#
module QuadraturePoints

#export chebysev_gauss, legendre_poly, legendre_gauss, lobatto_gauss

include("chebyshev_gauss.jl")
include("legendre_poly.jl")
include("legendre_gauss.jl")
include("lobatto_gauss.jl")

end
