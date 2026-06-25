# try_curved.jl
# Run from the package root (always start a fresh Julia session after code changes):
#   julia --project=. try_curved.jl


using PhyloNetworks
using PhyloPlots
using RCall


# # ── Network 1: simple network with one hybrid node ────────────────────────
# # This is the standard test network — one minor hybrid edge
# net1 = readnewick("(((A:.2,(B:.1)#H1:.1::0.9):.1,(C:.11,#H1:.01::0.1):.19):.1,D:.4);")

# # ── Network 2: bigger network with two hybrid nodes ───────────────────────
# net2 = readnewick("((((B)#H1:::0.6,C),((#H1:::0.4,D))#H2:::0.8),(#H2:::0.2,E));")

# # ── Network 3: a plain tree (no hybrids) — curved should do nothing ───────
# net3 = readnewick("((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);")

# # ─────────────────────────────────────────────────────────────────────────
# # TEST GROUP 1 — compare straight vs curved on the same network
# # Open three windows side by side
# # ─────────────────────────────────────────────────────────────────────────

# R"dev.new()"
# plot(net1)
# R"title(main='net1   curved=:none   (straight, default)')"

# R"dev.new()"
# plot(net1, curved=:minor)
# R"title(main='net1   curved=:minor')"

# R"dev.new()"
# plot(net1, curved=:both)
# R"title(main='net1   curved=:both')"

# # ─────────────────────────────────────────────────────────────────────────
# # TEST GROUP 2 — try different bend values
# # See how pronounced the curve gets
# # ─────────────────────────────────────────────────────────────────────────

# R"dev.new()"
# plot(net1, curved=:minor)
# R"title(main='curved=:minor  (single bend value, passes through midpoint)')"

# # ─────────────────────────────────────────────────────────────────────────
# # TEST GROUP 3 — with actual edge lengths turned on
# # This tests that the Bezier still works when x coordinates are
# # proportional to real branch lengths instead of uniform spacing
# # ─────────────────────────────────────────────────────────────────────────

# R"dev.new()"
# plot(net1, curved=:none, useedgelength=true)
# R"title(main='useedgelength=true   curved=:none')"

# R"dev.new()"
# plot(net1, curved=:minor, useedgelength=true)
# R"title(main='useedgelength=true   curved=:minor')"

# R"dev.new()"
# plot(net1, curved=:both, useedgelength=true)
# R"title(main='useedgelength=true   curved=:both')"

# # ─────────────────────────────────────────────────────────────────────────
# # TEST GROUP 4 — different styles
# # ─────────────────────────────────────────────────────────────────────────

# R"dev.new()"
# plot(net1, curved=:minor, style=:fulltree)
# R"title(main='style=:fulltree   curved=:minor')"

# R"dev.new()"
# plot(net1, curved=:minor, style=:majortree)
# R"title(main='style=:majortree   curved=:minor')"

# # ─────────────────────────────────────────────────────────────────────────
# # TEST GROUP 5 — with colors and edge widths
# # Confirm curved still respects color and width settings
# # ─────────────────────────────────────────────────────────────────────────

# R"dev.new()"
# plot(net1, curved=:minor,
#      minorhybridedgecolor="red",
#      majorhybridedgecolor="blue")
# R"title(main='curved=:minor   custom hybrid colors')"

# R"dev.new()"
# plot(net1, curved=:both,
#      minorhybridedgecolor="red",
#      majorhybridedgecolor="blue",
#      edgewidth=3)
# R"title(main='curved=:both   thick edges')"

# # ─────────────────────────────────────────────────────────────────────────
# # TEST GROUP 6 — bigger network with two hybrids
# # ─────────────────────────────────────────────────────────────────────────

# R"dev.new()"
# plot(net2, curved=:none)
# R"title(main='net2 (2 hybrids)   curved=:none')"

# R"dev.new()"
# plot(net2, curved=:minor)
# R"title(main='net2 (2 hybrids)   curved=:minor')"

# R"dev.new()"
# plot(net2, curved=:both)
# R"title(main='net2 (2 hybrids)   curved=:both')"

# # ─────────────────────────────────────────────────────────────────────────
# # TEST GROUP 7 — plain tree, no hybrids
# # curved should make no visual difference here
# # ─────────────────────────────────────────────────────────────────────────

# R"dev.new()"
# plot(net3, curved=:none)
# R"title(main='plain tree   curved=:none')"

# R"dev.new()"
# plot(net3, curved=:minor)
# R"title(main='plain tree   curved=:minor  (should look identical)')"

# # ─────────────────────────────────────────────────────────────────────────
# # TEST GROUP 8 — test error handling
# # These should print error messages, NOT crash Julia
# # ─────────────────────────────────────────────────────────────────────────

# println("\n── error handling tests ──")

# try
#     plot(net1, curved=:diagonal)
#     println("WRONG: should have thrown an error")
# catch e
#     println("CORRECT: bad curved value caught → ", e.msg)
# end

# try
#     plot(net1, curved=:minor, bend=-0.1)
#     println("WRONG: should have thrown an error")
# catch e
#     println("CORRECT: negative bend caught → ", e.msg)
# end

# try
#     plot(net1, curved=:minor, bend=0.0)
#     println("WRONG: should have thrown an error")
# catch e
#     println("CORRECT: zero bend caught → ", e.msg)
# end

# # ─────────────────────────────────────────────────────────────────────────
# # TEST GROUP 9 — test _quadbez_control directly
# # Print the control poin


# net1 = readnewick("(A:3.3,((B:1.5,#H1:1.2):1.5,((C:1.8)#H1:1,D:1.1):.2):0.3);");
# net2 = readnewick("(A:3.3,((B:1.5,#H1:0.2):1.5,((C:1)#H1:1.8,D:1.1):.2):0.3);");

# for e in net1.edge 
#     e.length *= 1000
# end
# for e in net2.edge 
#     e.length /= 1000
# end

# R"layout"([1 2])
# plot(net1, useedgelength=true, style = :majortree, showedgelength=true, arrowlen=0.1,curved = :both);
# plot(net2, useedgelength=true, style = :majortree, showedgelength=true, arrowlen=0.1,curved = :both);




# ─────────────────────────────────────────────────────────────────────────
# Plain tree (no hybrids):
#   root -> (parent of A,B,C)  and  sibling (parent of D,E)
# curved=:both should match curved=:none (nothing to curve)
# ─────────────────────────────────────────────────────────────────────────

net_abc = readnewick("((A:1.0,B:1.0,C:1.0):1.0,(D:1.0,E:1.0):1.0);")

R"dev.new()"
plot(net_abc, curved=:none)
R"title(main='net_abc: plain tree, curved=:none')"

R"dev.new()"
plot(net_abc, curved=:both)
R"title(main='net_abc: plain tree, curved=:both')"

net_op = readnewick("(A,((((B,(C)#H1:::0.7),(#H1:::0.3,D)))#H0,#H0),E);")

R"dev.new()"
plot(net_op, curved=:none, style=:majortree)
R"title(main='net_op majortree curved=:none baseline')"

R"dev.new()"
plot(net_op, curved=:both, style=:majortree)
R"title(main='net_op majortree curved=:both: H0 must show 2 arcs')"

R"dev.new()"
plot(net_op, curved=:both, style=:fulltree)
R"title(main='net_op fulltree curved=:both: H0 must show 2 arcs')"

R"dev.new()"
plot(net_op, curved=:minor, style=:majortree)
R"title(main='net_op majortree curved=:minor')"

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 10 — minor hybrid edges bowing in opposite directions
#
# Three cases that together demonstrate the bow-direction variety:
#
#   A. Same parent (overlap) — forced opposite bows by the fix.
#      Both major and minor share identical start+end coordinates.
#      Before the fix both bowed the same way and overlapped.
#      After the fix major bows one direction, minor bows the other.
#
#   B. Different parents, minor comes from BELOW the hybrid (chord goes
#      up-right). Minor bows leftward relative to the chord.
#
#   C. Different parents, minor comes from ABOVE the hybrid (chord goes
#      down-right). Minor bows leftward relative to the chord — which
#      is visually the opposite direction from case B.
#
# In cases B and C the bow direction is determined purely by _bow_left_sign,
# which picks whichever bow_sign puts cx further to the left (smaller x).
# ─────────────────────────────────────────────────────────────────────────

println("\n── TEST GROUP 10: bow direction variety ──")

# A — same parent: major and minor share identical chord → forced apart
#     node N has two children: H1 (major, gamma=0.6) and H1 (minor, gamma=0.4)
#     N is therefore BOTH the major and minor parent of H1.
# net10a = readnewick("((A,((B)#H1:::0.6,#H1:::0.4)),C);")

# R"dev.new()"
# plot(net10a, curved=:none, style=:majortree)
# R"title(main='10A :none  same-parent baseline')"

# R"dev.new()"
# plot(net10a, curved=:both, style=:majortree)
# R"title(main='10A :both  same-parent → 2 arcs in OPPOSITE directions')"

# R"dev.new()"
# plot(net10a, curved=:both, style=:fulltree)
# R"title(main='10A :both fulltree  same-parent → OPPOSITE arcs')"

# # B — minor parent is BELOW the hybrid (minor diagonal goes up-right)
# #     H1 major parent: node containing A and B (top of tree)
# #     H1 minor parent: node containing C and D (bottom of tree)
# #     → minor chord has positive dy → bows upper-left
# net10b = readnewick("(((A,(B)#H1:::0.7),C),(#H1:::0.3,D));")

# R"dev.new()"
# plot(net10b, curved=:none, style=:majortree)
# R"title(main='10B :none  minor-below baseline')"

# R"dev.new()"
# plot(net10b, curved=:both, style=:majortree)
# R"title(main='10B :both  minor from BELOW: chord goes up, bows upper-left')"

# # C — minor parent is ABOVE the hybrid (minor diagonal goes down-right)
# #     H1 major parent: node in the lower subtree (B, C, D)
# #     H1 minor parent: node in the upper subtree (A side)
# #     → minor chord has negative dy → bows lower-left (opposite to B)
# net10c = readnewick("((A,(#H1:::0.3,B)),(((C,(D)#H1:::0.7),E),F));")

# R"dev.new()"
# plot(net10c, curved=:none, style=:majortree)
# R"title(main='10C :none  minor-above baseline')"

# R"dev.new()"
# plot(net10c, curved=:both, style=:majortree)
# R"title(main='10C :both  minor from ABOVE: chord goes down, bows lower-left')"

# Side-by-side: 10B vs 10C shows the minor arc flipping direction
# R"layout"([1 2])
# plot(net10b, curved=:both, style=:majortree)
# R"title(main='10B: minor bows upper-left')"
# plot(net10c, curved=:both, style=:majortree)
# R"title(main='10C: minor bows lower-left (opposite)')"
# R"layout"(1)

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 11 — arrowhead pointing UP vs DOWN
#
# The direction the minor hybrid arrowhead points is set purely by
# TOPOLOGY (not edge lengths). The minor diagonal goes from the minor
# parent's node_y to the hybrid child's node_y. If the minor parent is
# placed LOWER in the tree (lower y) than the hybrid child, the
# arrowhead points UP; if higher, it points DOWN.
#
# Rule of thumb for :majortree layout:
#   - Minor parent LOWER than hybrid → arrowhead UP ↑
#     Achieved when: the minor-ref "companion" leaf appears LEFT (earlier)
#     in the Newick than the hybrid's major-parent subtree. LIFO traversal
#     assigns lower y to leaves that appear earlier.
#   - Minor parent HIGHER than hybrid → arrowhead DOWN ↓
#     Standard case: minor parent's subtree is on the right (higher y).
# ─────────────────────────────────────────────────────────────────────────

println("\n── TEST GROUP 11: minor arrowhead UP vs DOWN ──")

# ── Upward arrowheads ─────────────────────────────────────────────────

# U1: simplest upward case — 3 taxa
#   minor parent (B-side, left) gets y=1 < H1's y=2 → arrowhead UP
# net_u1 = readnewick("((B,#H1:::0.3),((C)#H1:::0.7,D));")
# R"dev.new()"
# plot(net_u1, curved=:none, style=:majortree)
# R"title(main='U1 :none  minor arrowhead UP (y_from < y_to)')"
# R"dev.new()"
# plot(net_u1, curved=:both, style=:majortree)
# R"title(main='U1 :both  arc bows upper-left, arrowhead points UP')"

# # U2: 5-taxon upward — same pattern, more taxa
# net_u2 = readnewick("(A,((B,#H1:::0.3),((C)#H1:::0.7,(D,E))));")
# R"dev.new()"
# plot(net_u2, curved=:none, style=:majortree)
# R"title(main='U2 :none  5-taxon, arrowhead UP')"
# R"dev.new()"
# plot(net_u2, curved=:both, style=:majortree)
# R"title(main='U2 :both  5-taxon arrowhead UP')"

# # ── Downward arrowheads (reference) ──────────────────────────────────

# # D1: standard layout — minor parent higher → arrowhead DOWN
# net_d1 = readnewick("(((A,(B)#H1:::0.7),C),(#H1:::0.3,D));")
# R"dev.new()"
# plot(net_d1, curved=:none, style=:majortree)
# R"title(main='D1 :none  minor arrowhead DOWN (y_from > y_to)')"
# R"dev.new()"
# plot(net_d1, curved=:both, style=:majortree)
# R"title(main='D1 :both  arc bows lower-left, arrowhead points DOWN')"

# # ── Side-by-side contrast ────────────────────────────────────────────
# R"layout"([1 2])
# plot(net_u1, curved=:both, style=:majortree)
# R"title(main='UP: minor parent below hybrid')"
# plot(net_d1, curved=:both, style=:majortree)
# R"title(main='DOWN: minor parent above hybrid')"
# R"layout"(1)

# ── Two hybrids: one UP, one DOWN in the same network ────────────────
# H1 minor → UP ↑  (B-side minor parent is below H1)
# H2 minor → DOWN ↓ (H2's minor parent is above H2)
# net_mixed = readnewick("((A,((B,#H1:::0.3),((C)#H1:::0.7,D))),(((E,(F)#H2:::0.8),G),(#H2:::0.2,H)));")
# R"dev.new()"
# plot(net_mixed, curved=:none, style=:majortree)
# R"title(main='mixed :none  H1-minor UP, H2-minor DOWN')"
# R"dev.new()"
# plot(net_mixed, curved=:both, style=:majortree)
# R"title(main='mixed :both  H1 arc UP, H2 arc DOWN — both visible')"
# R"dev.new()"
# plot(net_mixed, curved=:minor, style=:majortree)
# R"title(main='mixed :minor  minor curves only — H1 UP, H2 DOWN')"



# # -------------------------------------------------------------------------
# # TEST CASE 1: overlapping major/minor edge at hybrid node #H3
# # -------------------------------------------------------------------------

# net_overlap1 = readnewick(
#     "((#H3:1.3,(A:3,((B:1.0)#H3:1.0,b1:2.0,b0:2.0):0.5):1.5),D:3);"
# )

# R"dev.new()"
# plot(net_overlap1, curved=:none, style=:majortree)
# R"title(main='net_overlap1 majortree curved=:none baseline')"

# R"dev.new()"
# plot(net_overlap1, curved=:both, style=:majortree)
# R"title(main='net_overlap1 majortree curved=:both: #H3 major/minor edges must show 2 arcs')"

# # -------------------------------------------------------------------------
# # TEST CASE 2: overlapping major/minor edges involving #H1, #H2, and #H4
# # -------------------------------------------------------------------------

# net_overlap2 = readnewick(
#     "((((b1:2,(B:1)#H3:1,b0:2):.5,#H3:1.3):.5,#H1:2):1,(((C:1)#H2:1)#H1:1,(((D:0.1)#H4:1,#H4:1.5):1,#H2:0.2):1):1);"
# )

# R"dev.new()"
# plot(net_overlap2, curved=:none, style=:majortree)
# R"title(main='net_overlap2 majortree curved=:none baseline')"

# R"dev.new()"
# plot(net_overlap2, curved=:both, style=:majortree)
# R"title(main='net_overlap2 majortree curved=:both: overlapping hybrid edges must show distinct arcs')"

# ============================================================
# Test 1
# Hybrids after straight children
# ============================================================

net_test1 = readnewick(
    "(((a:1,b:1,(h1:1)#H1:1,(h2:1)#H2:1):1,#H1:1):1,#H2:1);"
)

R"dev.new()"
plot(net_test1, curved=:none, style=:majortree)
R"title(main='test1 majortree curved=:none')"

R"dev.new()"
plot(net_test1, curved=:minor, style=:majortree)
R"title(main='test1 majortree curved=:minor')"

R"dev.new()"
plot(net_test1, curved=:both, style=:majortree)
R"title(main='test1 majortree curved=:both')"

R"dev.new()"
plot(net_test1, curved=:none, style=:fulltree)
R"title(main='test1 fulltree curved=:none')"

R"dev.new()"
plot(net_test1, curved=:minor, style=:fulltree)
R"title(main='test1 fulltree curved=:minor')"

R"dev.new()"
plot(net_test1, curved=:both, style=:fulltree)
R"title(main='test1 fulltree curved=:both')"


# ============================================================
# Test 2
# Hybrids interleaved with straight children
# ============================================================

net_test2 = readnewick(
    "(((a:1,(h1:1)#H1:1,b:1,(h2:1)#H2:1):1,#H1:1):1,#H2:1);"
)

R"dev.new()"
plot(net_test2, curved=:none, style=:majortree)
R"title(main='test2 majortree curved=:none')"

R"dev.new()"
plot(net_test2, curved=:minor, style=:majortree)
R"title(main='test2 majortree curved=:minor')"

R"dev.new()"
plot(net_test2, curved=:both, style=:majortree)
R"title(main='test2 majortree curved=:both')"

R"dev.new()"
plot(net_test2, curved=:none, style=:fulltree)
R"title(main='test2 fulltree curved=:none')"

R"dev.new()"
plot(net_test2, curved=:minor, style=:fulltree)
R"title(main='test2 fulltree curved=:minor')"

R"dev.new()"
plot(net_test2, curved=:both, style=:fulltree)
R"title(main='test2 fulltree curved=:both')"

# ============================================================
# Node-bar Ymin/Ymax tests
# Verify: for nodes with Bézier-major children, the bar spans
#   [min(node_y, min(normal_children_y)), max(node_y, max(normal_children_y))]
# so that ordinary tree children are never visually disconnected.
# ============================================================

# --- net_test2: interleaved hybrids and tree children ----------
# Parent at y=2.5 has children: a(y=1), H1-Bézier(y=2), b(y=3), H2-Bézier(y=4)
# Ymin = min(2.5, 1.0) = 1.0   Ymax = max(2.5, 3.0) = 3.0
# Expected bar: single segment from y=1 to y=3 (b must be connected)

R"dev.new()"
R"layout"([1 2])
plot(net_test2, curved=:none, style=:majortree)
R"title(main='test2 :none  full bar 1→4  (reference)')"
plot(net_test2, curved=:both, style=:majortree)
R"title(main='test2 :both  bar 1→3 (b connected, H1/H2 as arcs)')"
R"layout"(1)

# --- net_op: single Bézier child above parent ------------------
# Parent at y=2.5 has children: B-tree(y=2), H0-Bézier(y=3)
# Ymin = min(2.5, 2.0) = 2.0   Ymax = max(2.5, 2.0) = 2.5
# Expected bar: y=2 to y=2.5

net_op = readnewick("(A,((((B,(C)#H1:::0.7),(#H1:::0.3,D)))#H0,#H0),E);")
R"dev.new()"
R"layout"([1 2])
plot(net_op, curved=:none, style=:majortree)
R"title(main='net_op :none  reference')"
plot(net_op, curved=:both, style=:majortree)
R"title(main='net_op :both  bar 2→2.5  (H0 as arc)')"
R"layout"(1)

# --- net_lsa2: single Bézier child below parent ---------------
# Two affected nodes:
#   ni=4: parent y=1.5, normal child C(y=2), Bézier H1(y=1)  → bar 1.5→2
#   ni=8: parent y=2.25, normal child (y=1.5), Bézier H2(y=3) → bar 1.5→2.25

net_lsa2 = readnewick("((((B)#H1:::0.6,C),((#H1:::0.4,D))#H2:::0.8),(#H2:::0.2,E));")
R"dev.new()"
R"layout"([1 2])
plot(net_lsa2, curved=:none, style=:majortree)
R"title(main='net_lsa2 :none  reference')"
plot(net_lsa2, curved=:both, style=:majortree)
R"title(main='net_lsa2 :both  Ymin/Ymax bars')"
R"layout"(1)