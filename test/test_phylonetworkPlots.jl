@testset "correct coordinates" begin

@testset "basic, edgenode_coordinate" begin
  net = readnewick("(A:2.5,((B:1,#H1:0.5::0.1):1,(C:1,(D:0.5)#H1:0.5::0.9):1):0.5);")
  # in this order:
  # edge_xB, edge_xE, edge_yB, edge_yE,
  # node_x, node_y, node_yB, node_yE,
  # minoredge_xB, minoredge_xE, minoredge_yB, minoredge_yE,
  # xmin, xmax, ymin, ymax, rdisplacement
  @test PhyloPlots.edgenode_coordinates(net, true, :fulltree, false, false) == (
    [0.0, 1.5, 1.5, 0.5, 1.5, 2.0, 1.5, 0.5, 0.0],
    [2.5, 2.5, 2.0, 1.5, 2.5, 2.5, 2.0, 1.5, 0.5],
    [1.0, 2.0, 3.0, 2.5, 4.0, 5.0, 5.0, 4.5, 3.5],
    [1.0, 2.0, 3.0, 2.5, 4.0, 5.0, 5.0, 4.5, 3.5],
    [2.5, 2.5, 1.5, 2.5, 2.5, 2.0, 1.5, 0.5, 0.0],
    [1.0, 2.0, 2.5, 4.0, 5.0, 5.0, 4.5, 3.5, 3.0],
    [1.0, 2.0, 2.0, 4.0, 5.0, 5.0, 4.0, 2.5, 1.0],
    [1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 4.5, 3.5],
    [2.0], [2.0], [3.0], [5.0], 0.0, 2.5, 1.0, 5.0, 2.0)
  @test PhyloPlots.edgenode_coordinates(net, true, :majortree, false, false, false) == (
    [0.0, 1.5, 1.5, 0.5, 1.5, 2.0, 1.5, 0.5, 0.0],
    [2.5, 2.5, 1.5, 1.5, 2.5, 2.5, 2.0, 1.5, 0.5],
    [1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 3.0],
    [1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 3.0],
    [2.5, 2.5, 1.5, 2.5, 2.5, 2.0, 1.5, 0.5, 0.0],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 3.0, 2.5],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.5, 3.0],
    [1.5], [2.0], [2.0], [4.0], 0.0, 2.5, 1.0, 4, 2.0)
  @test PhyloPlots.edgenode_coordinates(net, false, :majortree, false, false, false) == (
    [0.0, 2.0, 2.0, 1.0, 2.0, 3.0, 2.0, 1.0, 0.0],
    [4.0, 4.0, 2.0, 2.0, 4.0, 4.0, 3.0, 2.0, 1.0],
    [1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 3.0],
    [1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 3.0],
    [4.0, 4.0, 2.0, 4.0, 4.0, 3.0, 2.0, 1.0, 0.0],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 3.0, 2.5],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.5, 3.0],
    [2.0], [3.0], [2.0], [4.0], 0.0, 4.0, 1.0, 4, 2.0)
  dat = DataFrame(node=[-5,-3,-4,5,100],bs=["90","95","99","mytip","bogus"],edge=[8,9,4,6,200]);
  @test_logs (:warn, "Some node numbers in the nodelabel data frame are not found in the network:\n 100") PhyloPlots.check_nodedataframe(net, dat);
  @test_logs (:warn, "nodelabel should have 2+ columns, the first one giving the node numbers (Integer)") PhyloPlots.check_nodedataframe(net, dat[!,2:3])
  dat = DataFrame(node=[-5,-3,missing,5],
                  bs=["90","95","99","mytip"],edge=[8,9,4,6]);
  @test PhyloPlots.check_nodedataframe(net, dat) == (true,
    DataFrame(node=[-5,-3,5],bs=["90","95","mytip"],edge=[8,9,6]))
  dat = DataFrame(node=[-5,-3,-4,5],bs=["90","95","99","mytips"]);
  @test PhyloPlots.prepare_nodedataframe(net,dat,true,true,true,collect(1.:9),collect(10.:18)) ==
    DataFrame(name=["A","B","","C","D","H1","","",""],
    num=["1","2","-4","4","5","3","-5","-3","-2"],
    lab=["","","99","","mytips","","90","95",""],
    lea=[true,true,false,true,true,false,false,false,false],
    x=collect(1.:9), y=collect(10.:18))
  dat = DataFrame(edge=[8,9,4,6,200],bs=["90","95","99","mytips","bogus"]);
  @test_logs (:warn, "Some edge numbers in the edgelabel data frame are not found in the network:\n 200") PhyloPlots.prepare_edgedataframe(
    net,dat,:fulltree,collect(1.:9),collect(10.:18),collect(19.:27),collect(28.:36),
    [6.5],[8.5],[24.5],[26.5], :minor, 0.3);
  dat = DataFrame(edge=[8,9,4,6],bs=[missing,"95","99","mytips"]);
  @test PhyloPlots.prepare_edgedataframe(net,dat,:majortree,collect(1.:9),collect(11.:19),
    collect(21.:29),collect(31.:39),[7.],[9.],[27.],[29.],:none,0.3) == (true, DataFrame(
    len=["2.5","1","0.5","1","1","0.5","0.5","1","0.5"],
    gam=["1","1","0.1","1","1","1","0.9","1","1"],
    num=["1","2","3","4","5","6","7","8","9"],
    lab=["","","","99","","mytips","","","95"], # if not missing: "90" second to last (row 8)
    hyb=[false,false,true,false,false,false,true,false,false],
    min=[false,false,true,false,false,false,false,false,false],
    x=collect(6.:14),y=collect(26.:34)))

  # example with level-2 network, non-tree child:
  # one hybrid node ends up as a leaf in the major tree.
  # no major child edge to follow to set coordinates
  net = readnewick("((((B)#H1:::0.2)#H2,((D,C,#H2:::0.8)S1,(#H1,A)S2)S3)S4);")
  @test_logs plot(net, shownodenumber=true, showgamma=true, style=:fulltree, curved=:none);
  @test PhyloPlots.edgenode_coordinates(net, false, :fulltree, false, false, false) == (
    [4.0, 3.0, 0.0, 2.0, 2.0, 2.0, 1.0, 3.0, 3.0, 1.0, 0.0],
    [5.0, 4.0, 3.0, 5.0, 5.0, 3.0, 2.0, 4.0, 5.0, 3.0, 1.0],
    [5.0, 4.0, 1.0, 2.0, 3.0, 4.0, 3.0, 5.0, 6.0, 5.5, 4.0],
    [5.0, 4.0, 1.0, 2.0, 3.0, 4.0, 3.0, 5.0, 6.0, 5.5, 4.0],
    [5.0, 4.0, 3.0, 5.0, 5.0, 2.0, 5.0, 3.0, 1.0, 0.0],
    [5.0, 5.0, 4.0, 2.0, 3.0, 3.0, 6.0, 5.5, 4.0, 3.5],
    [5.0, 5.0, 4.0, 2.0, 3.0, 2.0, 6.0, 5.0, 3.0, 1.0],
    [5.0, 5.0, 4.0, 2.0, 3.0, 4.0, 6.0, 6.0, 5.5, 4.0],
    [4.0, 3.0], [4.0, 3.0], [4.0, 1.0], [5.0, 4.0], 0.0, 5.0, 1.0, 6, 4.0)
  net = readnewick("((((B)#H1:::0.2)#H2,((D,C,#H2)S1,(#H1,A)S2)S3)S4);")
  @test_logs plot(net, shownodenumber=true, showgamma=true, style=:fulltree, curved=:none);
  @test PhyloPlots.edgenode_coordinates(net, false, :fulltree, false, false, false) == (
    [4.0, 3.0, 0.0, 2.0, 2.0, 2.0, 1.0, 3.0, 3.0, 1.0, 0.0],
    [5.0, 4.0, 3.0, 5.0, 5.0, 3.0, 2.0, 4.0, 5.0, 3.0, 1.0],
    [5.0, 1.0, 1.0, 2.0, 3.0, 4.0, 3.0, 5.0, 6.0, 5.5, 4.0],
    [5.0, 1.0, 1.0, 2.0, 3.0, 4.0, 3.0, 5.0, 6.0, 5.5, 4.0],
    [5.0, 4.0, 3.0, 5.0, 5.0, 2.0, 5.0, 3.0, 1.0, 0.0],
    [5.0, 5.0, 1.0, 2.0, 3.0, 3.0, 6.0, 5.5, 4.0, 3.5],
    [5.0, 5.0, 1.0, 2.0, 3.0, 2.0, 6.0, 5.0, 3.0, 1.0],
    [5.0, 5.0, 1.0, 2.0, 3.0, 4.0, 6.0, 6.0, 5.5, 4.0],
    [4.0, 3.0], [4.0, 3.0], [1.0, 4.0], [5.0, 1.0], 0.0, 5.0, 1.0, 6, 7.0)
end # basic

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
end # of curved edges

@testset "lsatree style" begin
  net = readnewick("((((a2,(a3)#H1),((#H1,(a5)#H2),#H2)),((a1)#H4)#H3),(((#H4,b1),#H3),(c2,c1)));")
  # to visualize the node indices:
  df = DataFrame(number = [n.number for n in net.node], index = [i for i in eachindex(net.node)])

  res = plot(net, showedgenumber=true, style=:lsatree) # optimizeRD=false
  @test res[:edge_data][[3,5,9,16,19],:y] ≈ # placement of (minor) hybrid edge annotation
    [3.871320343559643, 3.8284271247461903, 2.914213562373095, 7.742640687119286, 8.65685424949238]
  res = plot(net, nodelabel=df, style=:lsatree, curved=:minor);
  @test res[:node_y_lo][12] == 8 # vertical bar for hybrid node
  @test res[:node_y_hi][12] == 9

end # of lsatree style

@testset "reticulate displacement (rdisplacement)" begin
  # without optimization: rd=3 with majortree, rd=6 with fulltree, rd=18.5 with lsatree
  net = readnewick("((((a2,(a3)#H1),((#H1,(a5)#H2),#H2)),((a1)#H4)#H3),(((#H4,b1),#H3),(c2,c1)));")
  for (style, expected_rd) in ((:majortree, 3.0), (:fulltree, 6.0), (:lsatree, 18.5))
    rd_minor = PhyloPlots.edgenode_coordinates(net, false, style, false, false, true)[end] # curved=:minor / :none
    rd_both  = PhyloPlots.edgenode_coordinates(net, false, style, true,  false, true)[end] # curved=:both
    @test rd_minor ≈ expected_rd
    # majorcurved (curved=:both vs :minor) only bends how major edges are drawn;
    # it should not change the reticulate displacement cost itself
    @test rd_both == rd_minor
  end
end # of reticulate displacement

end
