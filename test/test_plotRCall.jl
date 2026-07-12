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

  @testset "node with nomajorchild" begin
  # limited testing for correctness
  net2 = readnewick("((((B)#H1)#H2,((D,C,#H2)S1,(#H1:::.8,A)S2)S3)S4);")
  ecd = Dict(1 => SubString("grey50"), 2 => 2, 3 => "violet") # color 2 = red in R
  res = (@test_logs plot(net2, style=:majortree, curved=:none, edgecolor=ecd,
    defaultedgecolor="turquoise"))
  @test keys(res) == (:xmin, :xmax, :ymin, :ymax, :node_x, :node_y,
    :node_y_lo, :node_y_hi, :edge_x_lo, :edge_x_hi, :edge_y_lo, :edge_y_hi,
    :arrow_x_lo, :arrow_x_hi, :arrow_y_lo, :arrow_y_hi,
    :node_data, :edge_data)
  @test res[:node_y_lo] == [3,3,2.9,1,2,1,4,3,1.5,2.5]
  @test res[:node_y_hi] == [3,3,2.9,1,2,2,4,4,3.5,2.9]
  @test res[:edge_x_lo] == [4,3,  0,  2,2,2,  1,  3,3,1,  0]
  @test res[:edge_x_hi] == [5,3,  3,  5,5,2,  2,  4,5,3,  1]
  @test res[:edge_y_lo] == [3,2.9,2.9,1,2,1.5,1.5,3,4,3.5,2.5]
  @test res[:edge_y_hi] == res[:edge_y_lo]
  @test res[:arrow_x_lo] == [3,2]
  @test res[:arrow_x_hi] == [4,3]
  @test res[:arrow_y_lo] == [2.9,1.5]
  @test res[:arrow_y_hi] == [3,  2.9]
  @test res[:node_data].lab[[1,3]] == ["",""]
  @test res[:edge_data][[2,5],2:4] == DataFrame(
      gam=["0.2","1"], num=["2","5"], lab=["",""])
  res = (@test_logs plot(net2, style=:fulltree, curved=:none, preorder=false,
    edgelabel=DataFrame(num=[2,6], annotate=[85.0001,90])))
  @test res[:node_y_lo] == [5,5,1,2,3,2,6,5,3,  1]
  @test res[:node_y_hi] == [5,5,1,2,3,4,6,6,5.5,4]
  @test res[:edge_x_lo] == [4,3,0,2,2,2,1,3,3,1,0]
  @test res[:edge_x_hi] == [5,4,3,5,5,3,2,4,5,3,1]
  @test res[:edge_y_lo] == [5,1,1,2,3,4,3,5,6,5.5,4]
  @test res[:edge_y_hi] == res[:edge_y_lo]
  @test res[:arrow_x_lo] == [4,3]
  @test res[:arrow_x_hi] == res[:arrow_x_lo] # bc :fulltree style
  @test res[:arrow_y_lo] == [1,4]
  @test res[:arrow_y_hi] == [5,1]
  @test res[:node_data].lab[[1,3]] == ["",""]
  @test res[:edge_data][[2,5],2:4] == DataFrame(
      gam=["0.2","1"], num=["2","5"], lab=["85",""])
  end # of nomajorchild

  @testset "curved edges" begin
  # tree-child, level-2. not used for testing so far
  #net_lsa2 = readnewick("((((B)#H1:::0.6,C),((#H1:::0.4,D))#H2:::0.8),(#H2:::0.2,E));")

  # next net: overlapping segments when style=:majortree,
  #           and edge overlapping minor hybrid != major partner
  net = readnewick("((#H3:1.3,(A:3,((B:1.0)#H3:1,b1:2.0):0.5):1.5):0.3,D:3);")
  res = plot(net, curved=:both, style=:majortree,
    shownodenumber=true, showedgenumber=true)
  @test res[:node_data][!,4:6] == DataFrame(
    lea = Bool[1,1,0,1,0,0,0,1,0],
    x = [5,5,4,5,3,2,1,5,0.0],
    y = [1,2,2,3,2.5,2,2,4,2.5],
  )
  @test res[:edge_data][!,5:8] == DataFrame(
    hyb = Bool[1,0,0,1,0,0,0,0,0],
    min = Bool[1,0,0,0,0,0,0,0,0],
    x = [2.5,3.5,4.5,3.5,4,2.5,1.5,.5,2.5],
    y = [1.85,1,2,2.0428932188134525,3,2.5,2,2,4]
  )
  #
  # next net: not time-consistent, minor "corner" younger than youngest tip,
  # overlapping hybrid edges with style=:majortree (without curving edges)
  net = readnewick("((((b1:2,(B:1)#H3:1::0.3,b0:2):.5,#H3:1.3):.5,#H1:2):1,(((C:1)#H2:1::0.8)#H1:1::0.9,(((D:0.1)#H4:1::0.6,#H4:1.5):1,#H2:0.2):1):1);")
  #
  res = plot(net, curved=:both, style=:majortree, showedgelength=true)
  @test res[[:xmin,:xmax,:ymin,:ymax]] == (xmin=-.5, xmax=5.5, ymin=.5, ymax=5.5)
  @test res[[:node_y,:node_y_lo,:node_y_hi]] == (
    node_y    = [1,3,3,2,1.5,2,2,4,4,4,5,5,5,5,4.5,3],
    node_y_lo = [1,3,3,2,1,1.5,2,4,4,4,5,5,5,5,4.5,2],
    node_y_hi = [1,3,3,2,2,  2,2,4,4,4,5,5,5,5,5,4.5])
  @test res[[:edge_x_lo,:edge_x_hi, :edge_y_lo,:edge_y_hi]] == (
    edge_x_lo = [3,4,3,  3,2,  2,1,1,0,4,3,1,  4,3,3,2,2,1,0],
    edge_x_hi = [5,5,3,  5,3,  4,2,1,1,5,4,3,  5,4,3,3,2,2,1],
    edge_y_lo = [1,3,1.5,2,1.5,2,2,2,2,4,4,4.5,5,5,5,5,5,5,4.5],
    edge_y_hi = [1,3,1.5,2,1.5,3,2,2,2,4,4,4,  5,5,5,5,5,5,4.5]
  )
  @test res[[:arrow_x_lo, :arrow_x_hi, :arrow_y_lo, :arrow_y_hi]] == (
    arrow_x_lo = [3,  1,3,2],
    arrow_x_hi = [4,  3,4,4],
    arrow_y_lo = [1.5,2,5,5],
    arrow_y_hi = [3,  4,5,4]
  )
  @test res[:node_data] == DataFrame(
    name = ["b1","B","b0","C","D"],
    num = ["1","2","4","6","8"],
    lab = ["","","","",""],
    lea = [true,true,true,true,true],
    x = [5,5,5,5,5],
    y = [1,3,2,4,5],
  )
  @test res[:edge_data][!,[3,5,6,7]] == DataFrame(
    num = string.(1:19),
    hyb = Bool[0,0,1,0,0,1,0,1,0,0,1,1,0,1,1,0,1,0,0],
    min = Bool[0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0],
    x = [4,4.5,3.5,4,2.5,3,1.5,2,.5,4.5,3.5,2,4.5,3.5,3.5,2.5,3,1.5,.5],
  )
  @test res[:edge_data][!,:y] ≈
    [1,3,2.871320343559643,2,1.5,2.914213562373095,2,3.8284271247461903,2,4,4,4.042893218813453,5,5,4.85,5,4.085786437626905,5,4.5]
  #
  # same net, new options
  res = plot(net, curved=:minor, style=:majortree, showgamma=true,
    shownodenumber=true, preorder=false)
  @test res[:node_data] == DataFrame(
    name = ["b1","B","H3","b0","","","","C","H2","H1","D","H4","","","",""],
    num = string.([1,2,3,4,-5,-4,-3,6,7,5,8,9,-11,-10,-7,-2]),
    lab = repeat([""],16),
    lea = Bool[1,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0],
    x = [5,5,4,5,3,2,1,5,4,3,5,4,3,2,1,0],
    y = [1,3,3,2,1.5,2,2,4,4,4,5,5,5,5,4.5,3],
  )
  @test res[:edge_data][!,5:7] == DataFrame(
    hyb = Bool[0,0,1,0,0,1,0,1,0,0,1,1,0,1,1,0,1,0,0],
    min = Bool[0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0],
    x = [4,4.5,3.5,4,2.5,3,1.5,2,.5,4.5,3.5,2,4.5,3.5,3.5,2.5,3,1.5,.5],
  )
  @test res[:edge_data][!,:y] ≈
    [1,3,2.871320343559643,2,1.5,3,2,3.8284271247461903,2,4,4,4,5,5,4.85,5,4.085786437626905,5,4.5]
  #
  res = plot(net, curved=:both, style=:fulltree,
    shownodenumber=true, showedgenumber=true, preorder=false)
  @test res[:node_data][!,[:x,:y]] == DataFrame(
    x = [5,5,4,5,3,2,1,5,4,3,5,4,3,2,1,0],
    y = [1,4,4,3,2,2.5,3,6,6,6,7,7,7.5,8,7.5,5],
  )
  @test res[:edge_data][!,:x] ==
    [4,4.5,3.5,4,2.5,3,1.5,2,0.5,4.5,3.5,2,4.5,3.5,3.5,2.5,3,1.5,0.5]
  @test res[:edge_data][!,:y] ≈
    [1,4,2,3,2,3.871320343559643,2.5,5,3,6,6,6.128679656440358,7,7.042893218813453,8,7.5,9,8,7.5]
  #
  res = plot(net, curved=:minor, style=:fulltree, useedgelength=true,
    shownodenumber=true, showedgelength=true, preorder=false)
  @test res[:node_data][!,[:x,:y]] == DataFrame(
    x = [4,3.8,2.8,4,2,1.5,1,4,3,2,4.1,4,3,2,1,0],
    y = [1,4,4,3,2,2.5,3,6,6,6,7,7,7.5,8,7.5,5],
  )
  @test res[:edge_data][!,[:x,:y]] == DataFrame(
    x = [3,3.3,2.5,3,1.75,2.15,1.25,2,.5,3.5,2.5,1.5,4.05,3.5,3.75,2.5,2.1,1.5,.5],
    y = [1,4,2,3,2,4,2.5,5,3,6,6,6,7,7,8,7.5,9,8,7.5],
  )
  end
end
