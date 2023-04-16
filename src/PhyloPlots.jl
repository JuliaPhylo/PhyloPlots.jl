__precompile__()

module PhyloPlots

# standard libraries
using DataFrames: _broadcast_unalias_helper, LatexTableFormat
using Markdown
using Printf: @printf, @sprintf

# dependencies that need explicit declaration
using DataFrames
using RCall
import RCall.sexp
using PhyloNetworks

export plot
export sexp

include("phylonetworksPlots.jl")
include("plotRCall.jl")
include("substitutionmodels.jl")
include("rexport.jl")

end # of module
