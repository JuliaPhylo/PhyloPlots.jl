@testset "RCall-based plots" begin
  # testing for absence of errors, not for correctness

  # network rooted at a leaf: test for no error in warning message
  net = readnewick("(A:1,B:1);"); net.rooti = 2
  @test_logs (:warn, r"rootonedge!\(network_name, 2\)") plot(net)
  #net = readnewick("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
  net = readnewick("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");

  @test_logs plot(net);
  @test_logs plot(net, edgewidth=Dict(1=>4, 28=>28))
  @test_logs (:warn, "At least one non-missing edge length: plotting any missing length as 1.0") plot(net, useedgelength=true);
  @test_logs plot(net, showtiplabel=false);
  @test_logs plot(net, shownodenumber=true, shownodelabel=true);
  @test_logs plot(net, tipoffset=1, showgamma=true);
  @test_logs plot(net, showedgelength=true, showedgenumber=true);
  @test_logs plot(net, edgecolor="tomato4", minorhybridedgecolor="skyblue",
          majorhybridedgecolor="tan");
  dat = DataFrame(node=[-5,-10,-1],bs=["90","95","100"],edge=[11,22,26]);
  @test_logs (:warn, "Some node numbers in the nodelabel data frame are not found in the network:\n -1") plot(net, nodelabel=dat);
  @test_logs plot(net, edgelabel=dat[!,[:edge,:bs]]);

  @test_logs plot(net, style=:majortree, arrowlen=0.1);
  @test_logs (:warn, "Style bogus is unknown. Defaulted to :fulltree.") plot(net, style=:bogus);


  # coverage for the nomajorchild
  net2 = readnewick("((((B)#H1)#H2,((D,C,#H2)S1,(#H1:::.8,A)S2)S3)S4);")
  res = (@test_logs plot(net2, style=:majortree))
  @test keys(res) == (:xmin, :xmax, :ymin, :ymax, :node_x, :node_y,
    :node_y_lo, :node_y_hi, :edge_x_lo, :edge_x_hi, :edge_y_lo, :edge_y_hi,
    :arrow_x_lo, :arrow_x_hi, :arrow_y_lo, :arrow_y_hi,
    :node_data, :edge_data)
  @test res[:node_y_lo] == [3,3,2.9,1,2,1,4,3,1.5,2.5]
  @test res[:node_y_hi] == [3,3,2.9,1,2,2,4,4,3.5,2.9]
  @test res[:edge_x_lo] == [5,4,  1,  3,3,3,  2,  4,4,2,  1]
  @test res[:edge_x_hi] == [6,4,  4,  6,6,3,  3,  5,6,4,  2]
  @test res[:edge_y_lo] == [3,2.9,2.9,1,2,1.5,1.5,3,4,3.5,2.5]
  @test res[:edge_y_hi] == res[:edge_y_lo]
  @test res[:arrow_x_lo] == [4,3]
  @test res[:arrow_x_hi] == [5,4]
  @test res[:arrow_y_lo] == [2.9,1.5]
  @test res[:arrow_y_hi] == [3,  2.9]
  res = (@test_logs plot(net2, preorder=false))
  @test res[:node_y_lo] == [5,5,1,2,3,2,6,5,3,  1]
  @test res[:node_y_hi] == [5,5,1,2,3,4,6,6,5.5,4.25]
  @test res[:edge_x_lo] == [5,4,1,3,3,3,2,4,4,2,1]
  @test res[:edge_x_hi] == [6,5,4,6,6,4,3,5,6,4,2]
  @test res[:edge_y_lo] == [5,1,1,2,3,4,3,5,6,5.5,4.25]
  @test res[:edge_y_hi] == res[:edge_y_lo]
  @test res[:arrow_x_lo] == [5,4]
  @test res[:arrow_x_hi] == res[:arrow_x_lo] # bc :fulltree style
  @test res[:arrow_y_lo] == [1,4]
  @test res[:arrow_y_hi] == [5,1]
end
