@testset "RCall-based plots" begin

  @testset "basic case and warnings" begin
  # testing for absence of errors, not for correctness

  # network rooted at a leaf: test for no error in warning message
  net = readnewick("(A:1,B:1);"); net.rooti = 2
  @test_logs (:warn, r"rootonedge!\(network_name, 2\)") plot(net)
  #net = readnewick("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
  net = readnewick("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");

  @test_logs plot(net, style=:fulltree, curved=:none);
  @test_logs plot(net, style=:fulltree, curved=:both,
    edgewidth=Dict(1=>4, 9=>4, 2=>4, 28=>28), preorder=false)
  @test_logs (:warn, "At least one non-missing edge length: plotting any missing length as 1.0") plot(net,
    curved=:minor, preorder=false, useedgelength=true);
  @test_logs plot(net, style=:majortree, curved=:minor, showtiplabel=false, preorder=false);
  @test_logs plot(net, style=:majortree, curved=:both, shownodenumber=true, shownodelabel=true, preorder=false);
  @test_logs plot(net, style=:majortree, curved=:none, tipoffset=1,
    showgamma=true, minorlinetype=3, preorder=false); # 3=dotted
  @test_logs plot(net, style=:fulltree, curved=:none,
    showedgelength=true, showedgenumber=true, preorder=false);
  @test_logs plot(net, style=:fulltree, curved=:minor, preorder=false,
    edgecolor="tomato4", minorhybridedgecolor="skyblue", majorhybridedgecolor="tan");
  dat = DataFrame(node=[-5,-10,-1],bs=["90","95","100"],edge=[11,22,26]);
  @test_logs (:warn, "Some node numbers in the nodelabel data frame are not found in the network:\n -1") plot(net,
    nodelabel=dat, preorder=false);
  @test_logs plot(net, style=:majortree, curved=:both,
    edgelabel=dat[!,[:edge,:bs]], preorder=false);
  @test_logs plot(net, style=:majortree, curved=:none, arrowlen=0.1, preorder=false);
  @test_logs (:warn, "Style bogus is unknown. Defaulted to :majortree.") plot(net,
    style=:bogus, preorder=false);
  end # of basic case

end
