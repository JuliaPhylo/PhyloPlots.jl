__precompile__()

module PhyloPlots

# standard libraries
using Markdown
using Printf: @printf, @sprintf

# dependencies that need explicit declaration
using DataFrames
using ColorTypes # used by Gadfly already. To resolve data type names (Colorant)
using Gadfly
import Gadfly.plot
using RCall
import RCall.sexp
using PhyloNetworks

export plot
export sexp

include("phylonetworksPlots.jl")
include("plotGadfly.jl")
include("plotRCall.jl")
include("substitutionmodels.jl")
include("rexport.jl")

end # of module
