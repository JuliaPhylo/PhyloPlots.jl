using RCall
using PhyloPlots
using Test

using DataFrames
using PhyloNetworks # to read topologies

@testset "PhyloPlots Tests" begin
  include("test_phylonetworkPlots.jl")
  include("test_plotRCall.jl")  # uses DataFrames to test annotations
  include("test_substitutionmodels.jl")
  include("test_rexport.jl")    # uses RCall, but NOT library(ape)
end
