using RCall # before using PhyloNetworks and PhyloPlots
using PhyloPlots
using Base.Test

using DataFrames
using PhyloNetworks # to read topologies
using Colors # is REQUIREd by Gadfly

@testset "PhyloPlots Tests" begin
  # each file should be its own testset, to run them all even if one has a failure
  include("test_phylonetworkPlots.jl")
  #include("test_plotGadfly.jl") # uses Colors
  include("test_plotRCall.jl")  # uses DataFrames to test annotations
  include("test_substitutionmodels.jl")
  include("test_rexport.jl")    # uses RCall, but NOT library(ape)
end
