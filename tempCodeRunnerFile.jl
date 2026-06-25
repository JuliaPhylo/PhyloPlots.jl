# try_curved.jl
# Run this from the package root with:w
#   julia --project=. try_curved.jl

using PhyloNetworks
using PhyloPlots
using RCall


# ── Network 1: simple network with one hybrid node ────────────────────────
# This is the standard test network — one minor hybrid edge
net1 = readnewick("(((A:.2,(B:.1)#H1:.1::0.9):.1,(C:.11,#H1:.01::0.1):.19):.1,D:.4);")

# ── Network 2: bigger network with two hybrid nodes ───────────────────────
net2 = readnewick("((((B)#H1:::0.6,C),((#H1:::0.4,D))#H2:::0.8),(#H2:::0.2,E));")

# ── Network 3: a plain tree (no hybrids) — curved should do nothing ───────
net3 = readnewick("((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);")

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 1 — compare straight vs curved on the same network
# Open three windows side by side
# ─────────────────────────────────────────────────────────────────────────

R"dev.new()"
plot(net1)
R"title(main='net1   curved=:none   (straight, default)')"

R"dev.new()"
plot(net1, curved=:minor)
R"title(main='net1   curved=:minor')"

R"dev.new()"
plot(net1, curved=:both)
R"title(main='net1   curved=:both')"

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 2 — try different bend values
# See how pronounced the curve gets
# ─────────────────────────────────────────────────────────────────────────

R"dev.new()"
plot(net1, curved=:minor, bend=0.1)
R"title(main='bend=0.1  (subtle)')"

R"dev.new()"
plot(net1, curved=:minor, bend=0.3)
R"title(main='bend=0.3  (default)')"

R"dev.new()"
plot(net1, curved=:minor, bend=0.6)
R"title(main='bend=0.6  (pronounced)')"

R"dev.new()"
plot(net1, curved=:minor, bend=1.0)
R"title(main='bend=1.0  (extreme)')"

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 3 — with actual edge lengths turned on
# This tests that the Bezier still works when x coordinates are
# proportional to real branch lengths instead of uniform spacing
# ─────────────────────────────────────────────────────────────────────────

R"dev.new()"
plot(net1, curved=:none, useedgelength=true)
R"title(main='useedgelength=true   curved=:none')"

R"dev.new()"
plot(net1, curved=:minor, useedgelength=true)
R"title(main='useedgelength=true   curved=:minor')"

R"dev.new()"
plot(net1, curved=:both, useedgelength=true)
R"title(main='useedgelength=true   curved=:both')"

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 4 — different styles
# ─────────────────────────────────────────────────────────────────────────

R"dev.new()"
plot(net1, curved=:minor, style=:fulltree)
R"title(main='style=:fulltree   curved=:minor')"

R"dev.new()"
plot(net1, curved=:minor, style=:majortree)
R"title(main='style=:majortree   curved=:minor')"

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 5 — with colors and edge widths
# Confirm curved still respects color and width settings
# ─────────────────────────────────────────────────────────────────────────

R"dev.new()"
plot(net1, curved=:minor,
     minorhybridedgecolor="red",
     majorhybridedgecolor="blue")
R"title(main='curved=:minor   custom hybrid colors')"

R"dev.new()"
plot(net1, curved=:both,
     minorhybridedgecolor="red",
     majorhybridedgecolor="blue",
     edgewidth=3)
R"title(main='curved=:both   thick edges')"

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 6 — bigger network with two hybrids
# ─────────────────────────────────────────────────────────────────────────

R"dev.new()"
plot(net2, curved=:none)
R"title(main='net2 (2 hybrids)   curved=:none')"

R"dev.new()"
plot(net2, curved=:minor)
R"title(main='net2 (2 hybrids)   curved=:minor')"

R"dev.new()"
plot(net2, curved=:both)
R"title(main='net2 (2 hybrids)   curved=:both')"

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 7 — plain tree, no hybrids
# curved should make no visual difference here
# ─────────────────────────────────────────────────────────────────────────

R"dev.new()"
plot(net3, curved=:none)
R"title(main='plain tree   curved=:none')"

R"dev.new()"
plot(net3, curved=:minor)
R"title(main='plain tree   curved=:minor  (should look identical)')"

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 8 — test error handling
# These should print error messages, NOT crash Julia
# ─────────────────────────────────────────────────────────────────────────

println("\n── error handling tests ──")

try
    plot(net1, curved=:diagonal)
    println("WRONG: should have thrown an error")
catch e
    println("CORRECT: bad curved value caught → ", e.msg)
end

try
    plot(net1, curved=:minor, bend=-0.1)
    println("WRONG: should have thrown an error")
catch e
    println("CORRECT: negative bend caught → ", e.msg)
end

try
    plot(net1, curved=:minor, bend=0.0)
    println("WRONG: should have thrown an error")
catch e
    println("CORRECT: zero bend caught → ", e.msg)
end

# ─────────────────────────────────────────────────────────────────────────
# TEST GROUP 9 — test _quadbez_control directly
# Print the control poin


net1 = readnewick("(A:3.3,((B:1.5,#H1:1.2):1.5,((C:1.8)#H1:1,D:1.1):.2):0.3);");
net2 = readnewick("(A:3.3,((B:1.5,#H1:0.2):1.5,((C:1)#H1:1.8,D:1.1):.2):0.3);");

for e in net1.edge 
    e.length *= 1000
end
for e in net2.edge 
    e.length /= 1000
end

R"layout"([1 2])
plot(net1, useedgelength=true, style = :majortree, showedgelength=true, arrowlen=0.1,curved = :both);
plot(net2, useedgelength=true, style = :majortree, showedgelength=true, arrowlen=0.1,curved = :both);

# more bending if less species less benidng if more species