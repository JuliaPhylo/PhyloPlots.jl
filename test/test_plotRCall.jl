@testset "RCall-based plots" begin

  @testset "basic case and warnings" begin
  # testing for absence of errors, not for correctness

  # network rooted at a leaf: test for no error in warning message
  net = readnewick("(A:1,B:1);"); net.rooti = 2
  @test_logs (:warn, r"rootonedge!\(network_name, 2\)") plot(net)
  #net = readnewick("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
  net = readnewick("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");

  @test_logs plot(net);
  @test_logs plot(net, curved=:both, edgewidth=Dict(1=>4, 9=>4, 2=>4, 28=>28))
  @test_logs (:warn, "At least one non-missing edge length: plotting any missing length as 1.0") plot(net,
    useedgelength=true, curved=:minor);
  @test_logs plot(net, showtiplabel=false, style=:majortree, curved=:minor);
  @test_logs plot(net, shownodenumber=true, shownodelabel=true, style=:majortree, curved=:both);
  @test_logs plot(net, tipoffset=1, showgamma=true, minorlinetype=3, style=:majortree); # 3=dotted
  @test_logs plot(net, showedgelength=true, showedgenumber=true);
  @test_logs plot(net, edgecolor="tomato4", minorhybridedgecolor="skyblue",
          majorhybridedgecolor="tan", curved=:minor);
  dat = DataFrame(node=[-5,-10,-1],bs=["90","95","100"],edge=[11,22,26]);
  @test_logs (:warn, "Some node numbers in the nodelabel data frame are not found in the network:\n -1") plot(net, nodelabel=dat);
  @test_logs plot(net, edgelabel=dat[!,[:edge,:bs]], style=:majortree, curved=:both);

  @test_logs plot(net, style=:majortree, arrowlen=0.1);
  @test_logs (:warn, "Style bogus is unknown. Defaulted to :fulltree.") plot(net, style=:bogus);
  end # of basic case

  @testset "node with nomajorchild" begin
  # limited testing for correctness
  net2 = readnewick("((((B)#H1)#H2,((D,C,#H2)S1,(#H1:::.8,A)S2)S3)S4);")
  ecd = Dict(1 => SubString("grey50"), 2 => 2, 3 => "violet") # color 2 = red in R
  res = (@test_logs plot(net2, style=:majortree, edgecolor=ecd,
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
  res = (@test_logs plot(net2, preorder=false,
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
  # next net: not time-consistent, minor "corner" younger than youngest tip,
  # overlapping hybrid edges with style=:majortree (without curving edges)
  net = readnewick("((((b1:2,(B:1)#H3:1::0.3,b0:2):.5,#H3:1.3):.5,#H1:2):1,(((C:1)#H2:1::0.8)#H1:1::0.9,(((D:0.1)#H4:1::0.6,#H4:1.5):1,#H2:0.2):1):1);")
  # next: bend forced bc overlapping segments when style=:majortree
  # edge overlapping minor hybrid != major partner
  net = readnewick("((#H3:1.3,(A:3,((B:1.0)#H3:1.0,b1:2.0):0.5):1.5),D:3);")

  end
end

@testset "curved hybrid bow direction" begin
  function h1_hybrid_controls(net, offset)
    c = PhyloPlots.edgenode_coordinates(net, false, true, true)
    edge_xB, edge_xE, _, edge_yE = c[1], c[2], c[3], c[4]
    node_y = c[6]
    hxB, hxE, hyB, hyE = c[9], c[10], c[11], c[12]
    h1_node = findfirst(n -> n.hybrid && n.name == "H1", net.node)
    i = findfirst(e -> e.hybrid && e.ismajor &&
                  PhyloNetworks.getchild(e) === net.node[h1_node], net.edge)
    pni = findfirst(x -> x === PhyloNetworks.getparent(net.edge[i]), net.node)
    major_seg = (edge_xB[i], node_y[pni], edge_xE[i], edge_yE[i])
    _, _, x2, y2 = major_seg
    minor_j = findfirst(j -> abs(hxE[j] - x2) < 1e-10 && abs(hyE[j] - y2) < 1e-10,
                        1:length(hxB))
    minor_seg = (hxB[minor_j], hyB[minor_j], hxE[minor_j], hyE[minor_j])
    major_bow = PhyloPlots._hybrid_bow_sign(major_seg; offset=offset, partner=minor_seg)
    minor_bow = PhyloPlots._hybrid_bow_sign(minor_seg; offset=offset, partner=major_seg)
    major_cx, major_cy, _ = PhyloPlots._quadbez_control(major_seg...;
        offset_override=offset, force_bow_sign=major_bow)
    minor_cx, minor_cy, _ = PhyloPlots._quadbez_control(minor_seg...;
        offset_override=offset, force_bow_sign=minor_bow)
    return (
        major=(cx=major_cx, cy=major_cy, mx=(major_seg[1] + major_seg[3]) / 2,
               y0=major_seg[2], y2=major_seg[4]),
        minor=(cx=minor_cx, cy=minor_cy, mx=(minor_seg[1] + minor_seg[3]) / 2,
               y0=minor_seg[2], y2=minor_seg[4]),
    )
  end

  net_op = readnewick("(A,((((B,(C)#H1:::0.7),(#H1:::0.3,D)))#H0,#H0),E);")
  net_lsa2 = readnewick("((((B)#H1:::0.6,C),((#H1:::0.4,D))#H2:::0.8),(#H2:::0.2,E));")
  offset_op = 0.3 * 5.0 / net_op.numtaxa
  offset_lsa2 = 0.3 * 4.0 / net_lsa2.numtaxa

  h1_op = h1_hybrid_controls(net_op, offset_op)
  h1_lsa2 = h1_hybrid_controls(net_lsa2, offset_lsa2)

  @test h1_op.major.cx < h1_op.major.mx
  @test h1_lsa2.major.cx < h1_lsa2.major.mx
  @test h1_op.major.cy > (h1_op.major.y0 + h1_op.major.y2) / 2
  @test h1_lsa2.major.cy < (h1_lsa2.major.y0 + h1_lsa2.major.y2) / 2

  @test_logs plot(net_op, curved=:both, style=:majortree)
  @test_logs plot(net_lsa2, curved=:both, style=:majortree)
  @test_logs plot(net_op, curved=:minor, style=:majortree)
end
