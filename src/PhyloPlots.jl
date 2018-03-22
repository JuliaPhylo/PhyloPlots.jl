__precompile__()

module PhyloPlots

using DataFrames
using ColorTypes # used by Gadfly already. To resolve data type names (Colorant)
using Gadfly
import Gadfly.plot
using RCall
using PhyloNetworks

export
    # network comparative methods
    expectationsPlot,
    predintPlot, # prediction intervals at ancestral nodes
    plot

include("phylonetworksPlots.jl")
include("plotGadfly.jl")
include("plotRCall.jl")
include("substitutionmodels.jl")

end # of module
