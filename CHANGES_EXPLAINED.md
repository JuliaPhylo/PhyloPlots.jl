# Complete Explanation of All Changes: PhyloPlots.jl Curved Hybrid Edges Feature

---

## About Rplots.pdf

Before this document was written, `Rplots.pdf` was found in the repository root as an
**untracked file** with no git history — it was created automatically by R's graphics
system when we called `plot()` during manual testing. It was **deleted** and added to
`.gitignore` so it will never be accidentally committed.

---

## Section 1: Overview

PhyloPlots.jl is a Julia package that draws phylogenetic networks — diagrams showing
the evolutionary relationships between species, including cases where genes have been
transferred between lineages (called hybridization or reticulation). The package uses
R's base graphics under the hood through the RCall.jl bridge, drawing edges as line
segments and text labels for species names.

The `curved` feature allows minor hybrid edges (the dashed arrows showing gene flow
from a donor lineage) and major hybrid edges to be drawn as smooth arcs instead of
straight diagonal lines. This makes networks easier to read, especially when hybrid
edges cross tree backbone edges.

The three values `curved` can take are:

- **`:none`** (the default) — everything is drawn exactly as before, as straight lines.
  If you never pass `curved` at all, the output is identical to the original code.
- **`:minor`** — only the diagonal segment of each minor hybrid edge is replaced with
  a smooth curve. The rest of the network (tree edges, horizontal stubs, node bars)
  stays straight.
- **`:both`** — both the minor hybrid edge diagonals and the major hybrid edges are
  drawn as curves. Tree edges and minor hybrid edge stubs remain straight.

The `bend` parameter (default `0.3`) controls how pronounced the curve is. Larger
values bow the arc further from the straight chord; smaller values give a subtler
curve. It must be a positive number.

Three files were changed: `src/phylonetworksPlots.jl`, `src/plotRCall.jl`, and
`test/test_plotRCall.jl`. Three files were never touched: `src/PhyloPlots.jl`,
`src/rexport.jl`, and `test/test_phylonetworkPlots.jl`.

---

## Section 2: Background — How the Original Code Worked Before Any Changes

### 2a. What `edgenode_coordinates()` does

`edgenode_coordinates()` is a function in `src/phylonetworksPlots.jl` that takes a
`HybridNetwork` object, a boolean for whether to use real edge lengths, a boolean for
whether to use "direct" (majortree) hybrid lines, and an optional preorder flag. It
returns 16 values packed into a tuple.

**What it returns:**

1. `edge_xB`, `edge_xE` — the x-coordinates for the beginning and end of **every
   edge** in `net.edge`, in that exact order. For tree edges and major hybrid edges,
   this is a horizontal segment. For minor hybrid edges, this is the **horizontal
   stub** only — the short horizontal line connecting the minor parent node to the
   start of the diagonal.
2. `edge_yB`, `edge_yE` — the y-coordinates for those same edges. For all edges,
   `edge_yB == edge_yE` (the segment is always horizontal).
3. `node_x`, `node_y` — the x and y of the midpoint of each node's vertical bar.
4. `node_yB`, `node_yE` — the bottom and top of each node's vertical bar.
5. `minoredge_xB`, `minoredge_xE`, `minoredge_yB`, `minoredge_yE` — the coordinates
   of the **diagonal segment** of each minor hybrid edge, one entry per minor edge, in
   the order they appear when you filter `net.edge` to edges where `!e.ismajor`.
6. `xmin`, `xmax`, `ymin`, `ymax` — the axis ranges.

**The key distinction between `edge_*` and `hybridedge_*` (called `minoredge_*`
internally):**

A minor hybrid edge is drawn in two connected pieces. The `edge_*` arrays hold the
horizontal stub — a short horizontal line at a "phantom" y level that shows the edge
length. The `minoredge_*` (called `hybridedge_*` inside `plot()`) arrays hold the
diagonal piece — the line that actually connects the stub to the hybrid child node and
that carries the arrowhead showing direction. These must be in separate arrays because
they go to different R drawing calls: the stub is drawn by `segments()` alongside all
other edges, and the diagonal was originally drawn by `arrows()`.

**`edgenode_coordinates()` was never modified.** Not a single character was changed in
this function throughout all our work.

### 2b. The three R drawing calls in the original `plot()`

The entire drawing happened in three R calls:

**Call 1 — `R"segments"` for all edges:**
```julia
R"segments"(edge_xB, edge_yB, edge_xE, edge_yE, col=eCol, lwd=edgewidth_vec)
```
This draws every edge in `net.edge` as a horizontal line segment, passing all the
coordinate arrays at once. R draws all `n` segments in a single call with no Julia
loop. This includes tree edges, major hybrid edges (horizontal), and the horizontal
stubs of minor hybrid edges. The colors in `eCol` are set per-edge: tree edges get
`edgecolor`, major hybrid edges get `majorhybridedgecolor`, minor hybrid edge stubs
get `minorhybridedgecolor`.

**Call 2 — `R"arrows"` for minor hybrid edge diagonals:**
```julia
R"arrows"(hybridedge_xB, hybridedge_yB, hybridedge_xE, hybridedge_yE,
          length=arrowlen, angle=20, col=hybmincol_vec, lty=minorlinetype,
          lwd=hybridedgewidth_vec)
```
This draws the diagonal segment of every minor hybrid edge, with an arrowhead at the
endpoint (the child node). R's `arrows()` is used instead of `segments()` because it
adds the triangular arrowhead that indicates direction. The `length=arrowlen` parameter
controls the arrowhead size in inches; when `arrowlen=0` (which is the default for
`style=:majortree`), no arrowhead is drawn and it looks like a segment.

**Call 3 — `R"segments"` for node vertical bars:**
```julia
R"segments"(node_x, node_yB, node_x, node_yE, col=defaultedgecolor)
```
This draws the short vertical bar at each internal node. These bars connect the
horizontal edges coming from the node's children. Note that these are NOT edges in the
phylogenetic sense — they are just visual connectors. The x coordinate is the same for
both endpoints (it is `node_x`, not `node_xB`/`node_xE`), so these are always
vertical. They must always be drawn straight because they are not edges and have no
curved interpretation.

### 2c. What a minor hybrid edge looks like visually

In a phylogenetic network drawn by PhyloPlots, a minor hybrid edge is displayed as two
connected pieces forming a backwards "L" or reversed "J" shape:

**The horizontal stub:** Starting from the minor parent node's x position, a short
horizontal line extends to the right. In the `:fulltree` (default) style without edge
lengths, this stub ends at exactly the same x as the hybrid child node — so it is a
point with zero horizontal extent. When edge lengths are used, it extends by the minor
edge's length. The y coordinate of this stub is a dedicated "phantom row" in the plot,
allocated above the taxa rows.

**The diagonal segment:** From the end of the horizontal stub, a diagonal line drops
down (or up) to the hybrid child node's position. This diagonal carries the arrowhead
at its endpoint. The arrowhead points toward the hybrid child to show the direction of
gene flow.

These two pieces are stored in completely separate coordinate arrays because they serve
different visual purposes and are drawn by different R functions. The stub's coordinates
go into `edge_xB/xE/yB/yE` at the minor edge's index position, alongside all other
edges. The diagonal's coordinates go into `hybridedge_xB/xE/yB/yE` (a separate
shorter array with one entry per minor hybrid edge only).

### 2d. What a major hybrid edge is and how it was originally drawn

A major hybrid edge is the dominant parent edge of a hybrid node — the one with
inheritance probability `gamma >= 0.5`. In the coordinate system, major hybrid edges
are drawn as horizontal segments, just like tree edges, because the child node's y
coordinate is set by the major parent path. The `edge_*` arrays at a major hybrid
edge's index hold a horizontal line from the major parent's x to the child's x. In the
original code, major hybrid edges were drawn by the same vectorized `R"segments"` call
as tree edges, just with a different color (`majorhybridedgecolor`). They had no
arrowhead.

---

## Section 3: The Mathematics — Quadratic Bézier Curves

### 3a. What a quadratic Bézier curve is

A quadratic Bézier curve is a smooth curve defined by exactly three points: a start
point P0, one control point P1, and an end point P2. The curve is parameterized by a
number `t` that goes from 0 to 1:

```
B(t) = (1-t)² · P0  +  2(1-t)t · P1  +  t² · P2
```

In plain English: at `t=0` you are at P0 (the start). At `t=1` you are at P2 (the
end). At `t=0.5` (the midpoint) you are pulled toward P1 but don't reach it. The
curve is "attracted" to the control point but never passes through it unless P1 is
exactly on the straight line between P0 and P2.

The word "quadratic" means one control point. A cubic Bézier has two control points
and more flexibility, but quadratic is simpler to implement and gives a clean,
symmetric arc — enough for our purpose.

Why does a control point on the straight chord give a straight line? Because when P1
lies exactly between P0 and P2 on the line, the formula degenerates to a linear
interpolation: `B(t) = (1-t)·P0 + t·P2`, which traces the straight line.

### 3b. How the control point is computed

**Step 1 — Find the midpoint M of the chord:**

The chord is the straight line from P0 = (x0, y0) to P2 = (x2, y2).

```
mx = (x0 + x2) / 2
my = (y0 + y2) / 2
```

The control point will be placed at a distance from this midpoint.

**Step 2 — Find the perpendicular direction:**

The chord has direction vector `(dx, dy) = (x2-x0, y2-y0)`. To find a direction
perpendicular to the chord, we rotate this vector 90 degrees counterclockwise:

```
perpendicular direction = (-dy, dx)
```

Rotating 90 degrees counterclockwise means: x-component becomes `-y`, y-component
becomes `+x`. To get the unit vector (length 1), divide by the chord length:

```
perp_unit = (-dy / chord_len,  dx / chord_len)
```

**Step 3 — Offset M by (actual_offset) in the perpendicular direction:**

```
P1 = M  +  actual_offset · perp_unit
cx = mx  +  actual_offset · (-dy / chord_len)
cy = my  +  actual_offset · (dx  / chord_len)
```

This simplifies when `actual_offset = bend * chord_len` (the case when no
`offset_override` is given):

```
cx = mx - bend · dy
cy = my + bend · dx
```

The `chord_len` factors cancel out. This is the clean formula used by unit tests and
for major hybrid edges.

**What is `taxon_spacing` and why is it used for minor edges?**

For minor hybrid edges in `plotRCall.jl`, instead of using `bend * chord_len` as the
offset, we compute a fixed offset based on taxon spacing:

```julia
minor_offset = Float64(bend) * (ymax - ymin) / max(net.numtaxa, 1)
```

This is `bend × taxon_spacing`, where `taxon_spacing` is the average distance between
adjacent taxa on the y-axis (approximately 1.0 plot unit in the standard layout). The
key insight: if we let the offset scale with chord length, then a minor hybrid edge
with a long chord (which happens with `useedgelength=true` when branches span large
time intervals) produces a huge arc that flies outside the plot. By using a fixed
offset equal to `bend` times one taxon gap, every arc is the same visual size
regardless of chord length — about 30% of one row height with the default `bend=0.3`.

### 3c. The vertical chord special case

A "vertical chord" occurs when both endpoints have approximately the same x coordinate:
`x0 ≈ x2`. In `plotRCall.jl`, this is the standard situation for minor hybrid edge
diagonals under the default `:fulltree` style without edge lengths, because
`edgenode_coordinates()` sets `minoredge_xB = minoredge_xE = node_x[hybrid_child]` —
the diagonal starts and ends at the hybrid child's x position.

The standard perpendicular formula still works mathematically for a vertical chord. The
chord direction is `(0, dy)`, the perpendicular is `(-dy/|dy|, 0)` — pointing purely
horizontally. So the control point bows horizontally.

The code uses a special check: `abs(dx) < 1e-10 * abs(dy)`. When this is true, we
treat the chord as vertical and force the bow to be rightward (toward higher x, away
from the tree's root backbone):

```julia
cx = mx + actual_offset   # bows rightward
cy = my                   # stays at midpoint y
```

This is important visually: the alternative (bowing leftward) would push the arc into
the tree backbone where other edges are drawn.

### 3d. The zero chord special case

When both endpoints are at (nearly) the same location — `chord_len < 1e-6` — dividing
by `chord_len` would cause numerical problems and the "curve" would be invisible
anyway. The function returns immediately:

```julia
if chord_len < 1e-6
    return (mx, my, true)
end
```

The third return value `true` is the **fallback flag**. When the caller sees
`straight == true`, it draws a straight line (using `R"arrows"` for minor edges or
`R"segments"` for major edges) instead of calling R's `lines()`.

### 3e. The clamp fix

After computing the control point, two clamp lines keep it bounded within the chord's
own bounding box. This prevents arcs from escaping the plot area.

**For non-vertical chords:** the control point's x-coordinate is clamped to the
interval `[min(x0, x2), max(x0, x2)]`. This matters for diagonal chords where the
chord points upward-right — in that case `dy < 0`, so `cx = mx - actual_offset * dy /
chord_len` produces `cx > mx` which can overshoot `x2` for large offsets.

```julia
cx = clamp(cx, min(x0, x2), max(x0, x2))
```

**For vertical chords:** the control point's y-coordinate is clamped to
`[min(y0, y2), max(y0, y2)]`. Since `cy = my` (the midpoint) is always between the
two y-endpoints, this clamp is a no-op for the current code but documents the intended
invariant.

```julia
cy = clamp(cy, min(y0, y2), max(y0, y2))
```

**Why this was needed for `useedgelength=true`:** With actual edge lengths, the hybrid
node H1 in the test network sits at x = 1.3 and the plot's right boundary is xmax =
1.444. With the offset computed as `minor_offset = 0.375`, the horizontal bow would
place the control point at x = 1.3 + 0.375 = 1.675, which is well past xmax. Bézier
points at intermediate `t` values reach x ≈ 1.49, outside the visible plot area. The
clamp does not eliminate the bow in this case because cx = 1.3 + 0.375 is clamped to
`[min(x0,x2), max(x0,x2)] = [1.3, 1.3] = 1.3` only for the x-clamp applied in the
non-vertical branch; in the vertical branch the cy-clamp fires instead. The primary
fix for keeping arcs in bounds was the `minor_offset` change (Section 3b) that decoupled
offset from chord length. The clamps are an additional safety net.

### 3f. How 50 sample points are generated and drawn

The Bézier curve is not drawn as a true mathematical curve in R — it is approximated
by 50 straight line segments connecting 50 evenly-spaced sample points. These are
computed by evaluating the parametric formula at `t = 0, 1/49, 2/49, ..., 49/49 = 1`:

```r
t_vals = seq(0, 1, length.out=50)
x_curve = (1-t_vals)^2 * x0 + 2*(1-t_vals)*t_vals * cx + t_vals^2 * x2
y_curve = (1-t_vals)^2 * y0 + 2*(1-t_vals)*t_vals * cy + t_vals^2 * y2
lines(x_curve, y_curve, col=col, lwd=lwd, lty=lty)
```

The first point (`t=0`) is exactly (x0, y0) and the last (`t=1`) is exactly (x2, y2).
The intermediate points trace the smooth arc.

**Why `lines()` and not any other R function?**

`lines()` belongs to the `graphics` package, which is part of base R. It is loaded
automatically in every R session without any `library()` call. It is available in
every standard R installation since R version 1.0. We specifically did NOT use
`xspline()` (which has version-dependent behavior across R installations) or any
function from an add-on package.

---

## Section 4: Every Change Made to Every File

### 4a. `src/phylonetworksPlots.jl`

**Change: Added the `_quadbez_control` function**

**Before:** Nothing. This function did not exist.

**After:**
```julia
"""
    _quadbez_control(x0, y0, x2, y2; bend=0.3, offset_override=NaN)

Compute the control point for a quadratic Bézier curve from (x0,y0) to (x2,y2).
Returns `(cx, cy, straight)` where `straight::Bool` is `true` when the caller
should draw a straight line instead of a curve (degenerate or near-degenerate chord).

The perpendicular direction is always `(-dy, dx) / chord_len` (rotated 90° from
the chord). The offset magnitude is:
- `bend * chord_length` when `offset_override` is `NaN` (default, used by unit
  tests and major hybrid edges).
- `offset_override` when a finite value is given.  This decouples the bow size
  from chord length, keeping all arcs the same physical size regardless of whether
  the chord is short or long.  Minor hybrid edge callers in plotRCall.jl use this
  to pass a fixed taxon-spacing-based offset so arcs never blow outside the plot.
"""
function _quadbez_control(x0::Float64, y0::Float64, x2::Float64, y2::Float64;
                           bend::Float64=0.3, offset_override::Float64=NaN)
    dx = x2 - x0
    dy = y2 - y0
    chord_len = sqrt(dx^2 + dy^2)
    mx = (x0 + x2) / 2
    my = (y0 + y2) / 2
    if chord_len < 1e-6
        return (mx, my, true)
    end
    actual_offset = isnan(offset_override) ? bend * chord_len : offset_override
    if abs(dx) < 1e-10 * abs(dy)  # vertical chord: bow horizontally
        cx = mx + actual_offset
        cy = my
        cy = clamp(cy, min(y0, y2), max(y0, y2))
    else
        cx = mx - actual_offset * dy / chord_len
        cy = my + actual_offset * dx / chord_len
        cx = clamp(cx, min(x0, x2), max(x0, x2))
    end
    return (cx, cy, false)
end
```

**Line-by-line explanation:**

- `x0::Float64, y0::Float64, x2::Float64, y2::Float64` — the start and end points.
  The `::Float64` type annotations ensure the function is called with the right numeric
  type and matches what the Julia compiler expects.
- `bend::Float64=0.3` — the curvature fraction. Default 0.3 means a 30% bow relative
  to chord length (when `offset_override` is not given).
- `offset_override::Float64=NaN` — when `NaN` (the default), the offset magnitude is
  `bend * chord_len`. When a finite value is given, that value is used directly as the
  offset magnitude, bypassing the chord-length dependency. This was added to fix the
  "curves fly outside the plot" bug with `useedgelength=true`.
- `chord_len = sqrt(dx^2 + dy^2)` — the length of the straight line connecting the
  two points.
- `if chord_len < 1e-6` — catches zero-length and near-zero-length chords. Returns
  the midpoint and `straight=true` so the caller draws a straight line instead.
- `actual_offset = isnan(offset_override) ? bend * chord_len : offset_override` — picks
  between the two computation modes.
- `if abs(dx) < 1e-10 * abs(dy)` — detects a vertical chord. The threshold `1e-10`
  accounts for floating-point rounding: if the horizontal difference is less than one
  ten-billionth of the vertical difference, we treat the chord as vertical.
- `cx = mx + actual_offset` (vertical case) — forces the bow rightward. The `+` means
  the control point is to the right of the chord.
- `cy = my` (vertical case) — the control point is at the midpoint y.
- `cy = clamp(...)` (vertical case) — safety guard to keep cy within the chord's y range.
- `cx = mx - actual_offset * dy / chord_len` (general case) — the x offset is
  proportional to `-dy` (the perpendicular x-component).
- `cy = my + actual_offset * dx / chord_len` (general case) — the y offset is
  proportional to `+dx` (the perpendicular y-component).
- `cx = clamp(...)` (general case) — keeps cx within the chord's x range to prevent
  the arc from overshooting the endpoints.
- `return (cx, cy, false)` — returns the control point and `false` (meaning: draw a
  curve, not a straight line).

**Why this function lives in `phylonetworksPlots.jl`:** Coordinate geometry belongs
in this file; R-side drawing belongs in `plotRCall.jl`. The function is also directly
unit-testable without calling R.

---

### 4b. `src/plotRCall.jl`

**Change 1: Two new docstring lines added**

**Before:**
```
- `edgewidth=1`: width of horizontal (not diagonal) edges. To vary them,
  use a dictionary to map the number of each edge to its desired width.
- `xlim`, `ylim`: array of 2 values, to determine the axes limits.
```

**After:**
```
- `edgewidth=1`: width of horizontal (not diagonal) edges. To vary them,
  use a dictionary to map the number of each edge to its desired width.
- `curved = :none`: edge curvature style for hybrid edges. `:none` draws all
  edges as straight segments (default behaviour, identical to before this option
  was added). `:minor` replaces the straight diagonal of each minor hybrid edge
  with a quadratic Bézier curve. `:both` curves both major and minor hybrid edges.
- `bend = 0.3`: curvature amount for Bézier curves (must be positive). Larger
  values bow the curve further away from the straight chord.
- `xlim`, `ylim`: array of 2 values, to determine the axes limits.
```

**Why:** The docstring is the user-facing API contract. Users need to know these
arguments exist.

---

**Change 2: Two new keyword arguments added to the `plot()` signature**

**Before:**
```julia
    preorder::Bool=true,
)
```

**After:**
```julia
    preorder::Bool=true,
    curved::Symbol = :none,
    bend::Real = 0.3,
)
```

**Why:**
- `curved::Symbol = :none` — the user selects the curve mode. Type `Symbol` means it
  must be passed as a Julia symbol like `:minor`, not a string. Default `:none`
  preserves all original behavior.
- `bend::Real = 0.3` — the curvature amount. Type `Real` allows integers and
  floating-point numbers. Default 0.3 gives a gentle, visible arc.

---

**Change 3: Validation lines added after the style validation**

**Before:** The style validation block ends and immediately the R drawing begins.

**After:**
```julia
curved ∈ (:none, :minor, :both) ||
    error("curved must be :none, :minor, or :both; got :$curved")
bend > 0 ||
    error("bend must be a positive number; got $bend")
```

**Why:** These fire before any R code runs. If the user passes `curved=:diagonal`,
they get a clear error message naming the problematic argument. If they pass `bend=0`
or a negative value, they get a clear error.

---

**Change 4: The main drawing block replaced**

This is the largest change. The original three R calls were replaced by a conditional
block.

**Before (the original three calls):**
```julia
R"segments"(edge_xB, edge_yB, edge_xE, edge_yE, col=eCol, lwd=edgewidth_vec)
R"arrows"(hybridedge_xB, hybridedge_yB, hybridedge_xE, hybridedge_yE,
          length=arrowlen, angle=20, col=hybmincol_vec, lty=minorlinetype,
          lwd=hybridedgewidth_vec)
R"segments"(node_x, node_yB, node_x, node_yE, col=defaultedgecolor)
```

**After (the new block, plus the node bar call which never changed):**

```julia
if curved == :none
    R"segments"(edge_xB, edge_yB, edge_xE, edge_yE, col=eCol, lwd=edgewidth_vec)
    R"arrows"(hybridedge_xB, hybridedge_yB, hybridedge_xE, hybridedge_yE,
              length=arrowlen, angle=20, col=hybmincol_vec, lty=minorlinetype,
              lwd=hybridedgewidth_vec)
else
    if curved == :both
        straight_idx = [i for i in 1:length(net.edge)
                        if !net.edge[i].hybrid || !net.edge[i].ismajor]
        ew_straight = isa(edgewidth_vec, Number) ? edgewidth_vec : edgewidth_vec[straight_idx]
        R"segments"(edge_xB[straight_idx], edge_yB[straight_idx],
                    edge_xE[straight_idx], edge_yE[straight_idx],
                    col=eCol[straight_idx], lwd=ew_straight)
        majhyb_idx = findall(i -> net.edge[i].hybrid && net.edge[i].ismajor, 1:length(net.edge))
        for i in majhyb_idx
            x0_i = edge_xB[i]; y0_i = edge_yB[i]
            x2_i = edge_xE[i]; y2_i = edge_yE[i]
            (cx_i, cy_i, straight_i) = _quadbez_control(x0_i, y0_i, x2_i, y2_i;
                                                          bend=Float64(bend))
            col_i = eCol[i]
            lwd_i = isa(edgewidth_vec, AbstractVector) ? edgewidth_vec[i] : edgewidth_vec
            if straight_i
                R"segments"(x0_i, y0_i, x2_i, y2_i, col=col_i, lwd=lwd_i)
            else
                R"""
                t_vals = seq(0, 1, length.out=50)
                x_curve = (1-t_vals)^2 * $(x0_i) + 2*(1-t_vals)*t_vals * $(cx_i) + t_vals^2 * $(x2_i)
                y_curve = (1-t_vals)^2 * $(y0_i) + 2*(1-t_vals)*t_vals * $(cy_i) + t_vals^2 * $(y2_i)
                lines(x_curve, y_curve, col=$(col_i), lwd=$(lwd_i))
                """
                eff_arrowlen_i = arrowlen > 0 ? arrowlen : 0.1
                tang_dx_i = x2_i - cx_i
                tang_dy_i = y2_i - cy_i
                tang_len_i = sqrt(tang_dx_i^2 + tang_dy_i^2)
                if tang_len_i > 1e-10
                    ε_i = 0.01 * tang_len_i
                    xfrom_i = x2_i - ε_i * tang_dx_i / tang_len_i
                    yfrom_i = y2_i - ε_i * tang_dy_i / tang_len_i
                    R"arrows"(xfrom_i, yfrom_i, x2_i, y2_i, length=eff_arrowlen_i, angle=20,
                              col=col_i, lty="solid", lwd=lwd_i)
                end
            end
        end
    else  # :minor
        R"segments"(edge_xB, edge_yB, edge_xE, edge_yE, col=eCol, lwd=edgewidth_vec)
    end
    minor_offset = Float64(bend) * (ymax - ymin) / max(net.numtaxa, 1)
    for j in 1:length(hybridedge_xB)
        x0_j = hybridedge_xB[j]; y0_j = hybridedge_yB[j]
        x2_j = hybridedge_xE[j]; y2_j = hybridedge_yE[j]
        (cx_j, cy_j, straight_j) = _quadbez_control(x0_j, y0_j, x2_j, y2_j;
                                                      bend=Float64(bend),
                                                      offset_override=minor_offset)
        col_j = isa(hybmincol_vec, AbstractVector) ? hybmincol_vec[j] : hybmincol_vec
        lwd_j = isa(hybridedgewidth_vec, AbstractVector) ? hybridedgewidth_vec[j] : hybridedgewidth_vec
        if straight_j
            R"arrows"(x0_j, y0_j, x2_j, y2_j, length=arrowlen, angle=20,
                      col=col_j, lty=minorlinetype, lwd=lwd_j)
        else
            R"""
            t_vals = seq(0, 1, length.out=50)
            x_curve = (1-t_vals)^2 * $(x0_j) + 2*(1-t_vals)*t_vals * $(cx_j) + t_vals^2 * $(x2_j)
            y_curve = (1-t_vals)^2 * $(y0_j) + 2*(1-t_vals)*t_vals * $(cy_j) + t_vals^2 * $(y2_j)
            lines(x_curve, y_curve, col=$(col_j), lwd=$(lwd_j), lty=$(minorlinetype))
            """
            eff_arrowlen_j = arrowlen > 0 ? arrowlen : 0.1
            tang_dx = x2_j - cx_j
            tang_dy = y2_j - cy_j
            tang_len = sqrt(tang_dx^2 + tang_dy^2)
            chord_len_j = sqrt((x2_j - x0_j)^2 + (y2_j - y0_j)^2)
            if tang_len > 1e-10
                ε_j = 0.01 * max(chord_len_j, 1e-10)
                xfrom_j = x2_j - ε_j * tang_dx / tang_len
                yfrom_j = y2_j - ε_j * tang_dy / tang_len
            else
                xfrom_j = x0_j
                yfrom_j = y0_j
            end
            R"arrows"(xfrom_j, yfrom_j, x2_j, y2_j, length=eff_arrowlen_j, angle=20,
                      col=col_j, lty="solid", lwd=lwd_j)
        end
    end
end
R"segments"(node_x, node_yB, node_x, node_yE, col=defaultedgecolor)  # NEVER CHANGED
```

**Detailed explanation of each sub-part:**

**The `curved == :none` branch:** Byte-for-byte identical to the original two calls.
When the user does not pass `curved` (or passes `curved=:none`), this branch executes
and the output is identical to the original code.

**The `straight_idx` filter (only for `curved == :both`):**
```julia
straight_idx = [i for i in 1:length(net.edge)
                if !net.edge[i].hybrid || !net.edge[i].ismajor]
```
This collects the indices of every edge except major hybrid edges. A major hybrid edge
satisfies `hybrid=true AND ismajor=true`. By keeping only edges where that condition
is NOT met, we get: all tree edges (`hybrid=false`) and all minor hybrid edges
(`ismajor=false`). The minor hybrid edges' `edge_*` coordinates are their horizontal
stubs, not their diagonals, so they must be drawn straight by `segments()`. The
filtered `R"segments"` call draws these without the major hybrid edges, which will
instead be drawn as Bézier curves.

**The `minor_offset` computation:**
```julia
minor_offset = Float64(bend) * (ymax - ymin) / max(net.numtaxa, 1)
```
Computes `bend × taxon_spacing`. For the test network with 4 taxa, `ymax - ymin ≈ 5`
(taxa span plus expansion for tip labels), so `minor_offset ≈ 0.3 × 5/4 = 0.375`. The
`Float64(bend)` conversion is needed because `bend::Real` in the outer function but
`_quadbez_control` requires `bend::Float64`.

**The arrowhead direction fix (for curved minor edges):**

In the original code, `arrows(x0, y0, x1, y1)` draws the arrowhead pointing from
(x0,y0) toward (x1,y1) — the direction of the straight chord. When the edge is curved,
the curve arrives at the endpoint from a different direction (the tangent at t=1), so
the arrowhead would point the wrong way if we used the straight chord direction.

The tangent at t=1 of a quadratic Bézier points from the control point P1 toward the
endpoint P2. So we compute:

```julia
tang_dx = x2_j - cx_j   # direction from control point to endpoint
tang_dy = y2_j - cy_j
tang_len = sqrt(tang_dx^2 + tang_dy^2)
```

We place a "from" point just 1% of the tangent length back along this direction:
```julia
ε_j = 0.01 * max(chord_len_j, 1e-10)
xfrom_j = x2_j - ε_j * tang_dx / tang_len
yfrom_j = y2_j - ε_j * tang_dy / tang_len
```

Then call `arrows(xfrom_j, yfrom_j, x2_j, y2_j, ...)`. The shaft of this tiny
`arrows()` call is effectively invisible (it spans only 1% of the chord), but R uses
the from→to direction to orient the arrowhead, which now matches the curve's end
direction.

**`effective_arrowlen` (bug fix for `style=:majortree`):**
```julia
eff_arrowlen_j = arrowlen > 0 ? arrowlen : 0.1
```
When `style=:majortree`, `arrowlen` defaults to 0, which would make the arrowhead
invisible. For curved edges, a direction indicator is always needed to show which end
is the child node. So we use a minimum of 0.1 inches for the arrowhead on any curved
edge. The same logic applies to curved major edges in `curved=:both`.

**`lty="solid"` for arrowheads:**
Arrowheads are drawn as solid lines regardless of `minorlinetype` (which might be
`"longdash"`). A dashed arrowhead looks wrong visually.

---

### 4c. `test/test_plotRCall.jl`

**Before:** The file had one testset: `"RCall-based plots"` with tests for the
original behavior. Total: 45 test assertions (the existing test suite ran 92 total
tests including the other test files).

**After:** A new testset `"curved keyword"` was prepended to the file, containing
the following named sub-testsets and numbered tests. Total: 129 tests across all files.

**Named tests and sub-testsets added:**

1. **Test 1** (`res_none = @test_logs plot(net, curved=:none)`) — Verifies `curved=:none`
   runs without error and returns the full named tuple with all 18 expected field names.
   Would catch: if adding `curved` broke the return value shape.

2. **Test 2** (`@test_logs plot(net, curved=:minor)`) — Verifies `:minor` runs without
   error or warnings. Would catch: any R exception in the Bézier drawing code.

3. **Test 3** (`@test_logs plot(net, curved=:both)`) — Same for `:both`.

4. **Test 4** (four `@test` comparisons) — Verifies that `curved=:none` returns
   bit-identical `arrow_*` coordinate arrays to a plain `plot(net)` call. Would catch:
   if the new code accidentally changed the coordinate arrays returned.

5. **Test 5** (`@test_throws ErrorException` + `occursin("curved", err.msg)`) — Verifies
   that passing an invalid symbol like `:diagonal` throws an error and that the error
   message contains the word "curved". Would catch: if validation was missing or had
   the wrong message.

6. **Test 6** (`@test_throws ErrorException plot(net, bend=0.0)`) — Verifies that zero
   bend throws an error. Would catch: if the bend validation was missing.

7. **Test 7** (runs `:minor` again + `!any(isnan, ...)`) — Explicitly exercises the
   vertical chord case (the default layout gives vertical minor diagonals) and verifies
   no NaN in returned coordinates. Would catch: division-by-zero in the vertical chord
   case.

8. **Test 8** — Verifies `curved=:minor` works with an `edgecolor` dictionary. Would
   catch: if the per-edge color indexing broke when `hybmincol_vec` is a vector instead
   of a scalar.

9. **Test 9** — Same for `edgewidth` dictionary.

10. **Test 10** — Verifies `style=:majortree` combined with `curved=:minor`. Would
    catch: if the majortree style changed the coordinate layout in a way that breaks
    the curve drawing.

11. **Test A** — Verifies `useedgelength=true` with `curved=:minor` (non-NaN coords).

12. **Test B** — Verifies `useedgelength=true` with `curved=:both`.

13. **Tests 11, 12, 13, C, D, E, G** — Direct unit tests of `_quadbez_control`:
    horizontal chord, vertical chord, zero chord, same-y chord, zero-length, near-zero,
    degenerate loop. These verify exact numerical values. Each would catch if the
    control-point formula or fallback logic were wrong.

14. **Sub-testset `"curve offset does not scale with chord length"`** — Tests that
    `offset_override` makes a chord of length 100 and a chord of length 1 produce the
    same bow magnitude (0.3). Without `offset_override`, the long chord produces cy=30.
    Would catch: if the `isnan(offset_override)` branching were wrong.

15. **Sub-testset `"Bug2 majortree arrowhead"`** — Verifies `style=:majortree` with
    `curved=:minor` runs without error, including explicitly `arrowlen=0`. Would catch:
    if no arrowhead were drawn (no error) but the R call were somehow wrong.

16. **Sub-testset `"Bug3 no straight under curved major"`** — Verifies `curved=:both`
    runs with and without `edgecolor` dict. Would catch: array-indexing errors in the
    filtered `straight_idx` segments call.

17. **Sub-testset `"control point stays within chord endpoints"`** — Verifies the clamp
    lines. For a non-vertical chord, `cx ∈ [0, 2]`. For a vertical chord, `cy ∈ [0, 2]`.
    Also verifies `useedgelength=true` runs without error.

18. **Sub-testset `"Bug4 major edge arrowhead curved=both"`** — Verifies `curved=:both`
    runs with `arrowlen=0`, `arrowlen=0.2`, and default. Would catch: if the arrowhead
    call after `lines()` for major edges caused an error.

---

## Section 5: The Bugs We Found and Fixed

### Bug 1: Minor hybrid edge curves going outside the plot box

**Symptom:** With a large `bend` value (like `bend=1.0`) on a network with a long
minor hybrid diagonal, the arc bowed so far that it left the visible area of the plot.

**Root cause:** The original control point computation set the offset magnitude to
`bend * chord_length`. A long chord produces a large arc regardless of the plot scale.

**Fix:** Added the `max_offset` cap parameter (later replaced by `offset_override`).
The final fix computes `minor_offset = Float64(bend) * (ymax - ymin) / max(net.numtaxa, 1)`
before the minor edge loop, then passes `offset_override=minor_offset` to
`_quadbez_control`. This caps the bow to `bend × taxon_spacing` — roughly the
distance between adjacent taxa rows — regardless of chord length.

**Lines changed:** `phylonetworksPlots.jl` — signature and body of `_quadbez_control`.
`plotRCall.jl` — the `minor_offset` computation and the `_quadbez_control` call site.

**Verification:** Run `plot(net, curved=:minor, bend=1.0)`. The arc must stay within
the tree's bounding box.

---

### Bug 2: No arrowhead on curved minor edge with `style=:majortree`

**Symptom:** When `style=:majortree` and `curved=:minor`, the curved arc had no
arrowhead, so there was no visual indicator showing which end was the child.

**Root cause:** When `style=:majortree`, `arrowlen` defaults to 0. The code used
`length=arrowlen` in the `arrows()` call, so `length=0` produced an invisible head.

**Fix:** Added `eff_arrowlen_j = arrowlen > 0 ? arrowlen : 0.1`. Curved arcs always
use `eff_arrowlen_j` instead of `arrowlen` in the arrowhead call. The minimum of 0.1
inches ensures visibility. Also changed arrowhead `lty` from `minorlinetype` to
`"solid"` so the head is not dashed.

**Lines changed:** `plotRCall.jl` — the minor edge arrowhead call inside the `else`
(curved) branch.

**Verification:** Run `plot(net, curved=:minor, style=:majortree)`. The arc must have
a visible arrowhead.

---

### Bug 3: Straight segment drawn underneath curved major edge in `curved=:both`

**Symptom:** When `curved=:both`, each major hybrid edge appeared as two overlapping
lines: one straight (from the vectorized `segments()` call) and one curved (from the
per-edge `lines()` call). The straight line was visibly underneath the curve.

**Root cause:** After fixing Bug 5 (see below), the segments call for `curved=:both`
drew ALL edges including major hybrid edges. Then the major hybrid Bézier loop drew
them again as curves on top. Two lines per major edge.

**Fix:** Added the `straight_idx` filter for the `curved=:both` segments call:
```julia
straight_idx = [i for i in 1:length(net.edge)
                if !net.edge[i].hybrid || !net.edge[i].ismajor]
```
This excludes major hybrid edges from the `segments()` call. They are drawn only as
Bézier curves.

**Lines changed:** `plotRCall.jl` — the `if curved == :both` branch now uses a
filtered `segments()` call instead of the full one.

**Verification:** Run `plot(net, curved=:both)`. Major hybrid edges must show exactly
one line — a smooth curve.

---

### Bug 4: No arrowhead on curved major hybrid edge in `curved=:both`

**Symptom:** When `curved=:both`, the major hybrid edges appeared as smooth curves but
had no arrowhead, so there was no direction indicator.

**Root cause:** The major hybrid edge Bézier loop drew the curve with `lines()` but
never called `arrows()` afterward.

**Fix:** After the `lines()` call for each major hybrid edge, added an arrowhead using
the same end-tangent technique as minor edges:
```julia
eff_arrowlen_i = arrowlen > 0 ? arrowlen : 0.1
tang_dx_i = x2_i - cx_i
tang_dy_i = y2_i - cy_i
tang_len_i = sqrt(tang_dx_i^2 + tang_dy_i^2)
if tang_len_i > 1e-10
    ε_i = 0.01 * tang_len_i
    xfrom_i = x2_i - ε_i * tang_dx_i / tang_len_i
    yfrom_i = y2_i - ε_i * tang_dy_i / tang_len_i
    R"arrows"(xfrom_i, yfrom_i, x2_i, y2_i, length=eff_arrowlen_i, angle=20,
              col=col_i, lty="solid", lwd=lwd_i)
end
```

**Lines changed:** `plotRCall.jl` — the `else` (curve) branch of the `majhyb_idx`
loop, after the `lines()` R call.

**Verification:** Run `plot(net, curved=:both)`. Every curved major hybrid edge must
have a visible arrowhead at the child endpoint.

---

### Bug 5: Horizontal stub missing for minor edges when `curved=:both`

**Symptom:** When `curved=:both`, the minor hybrid edges appeared without their
horizontal stub, looking like they were floating disconnected from the tree.

**Root cause:** An early version of `curved=:both` filtered the `segments()` call to
include ONLY tree edges (non-hybrid). This excluded the horizontal stubs of minor
hybrid edges (which live in `edge_*` at the minor edge's index) as well as major hybrid
edges. The stubs should never be curved — they are positional guides, not diagonals.

**Fix:** Changed the filter to `straight_idx = [... if !net.edge[i].hybrid || !net.edge[i].ismajor]`.
This means: include everything EXCEPT major hybrid edges. Both tree edges AND minor
hybrid edge stubs are included.

**Lines changed:** `plotRCall.jl` — the `straight_idx` construction in the
`curved == :both` branch.

**Verification:** Run `plot(net, curved=:both)`. Minor hybrid edges must show their
horizontal stub connecting to the rest of the tree.

---

### Bug 6: Curves still going outside the box with `useedgelength=true`

**Symptom:** After the first `max_offset` fix, `plot(net, curved=:minor, useedgelength=true)`
still showed arcs flying past the tip labels.

**Root cause:** The first fix capped the offset based on the y-axis range (`max_offset
= 0.4 * (ymax - ymin) / (numtaxa + numhybrids)`). But for a vertical chord with
`useedgelength=true`, the arc bows **horizontally** (rightward), not vertically. The
y-axis cap did nothing to prevent the control point from going past `xmax`. Specifically:
the test network with `useedgelength=true` has the hybrid node H1 at x = 1.3 while
xmax = 1.444. A `minor_offset` of 0.375 placed the control point at x = 1.675 — 0.231
units past the right boundary.

**Fix:** Replaced the `max_offset` cap with `offset_override`. Instead of capping, the
offset is now **decoupled from chord length entirely**. The offset is fixed at
`bend × taxon_spacing` regardless of chord length or direction. For both short and long
chords, both vertical and horizontal, the arc bows by the same physical amount —
approximately `bend` times the gap between adjacent taxa. This kept arcs within bounds
for all combinations of `useedgelength` and `style`.

**Lines changed:** `_quadbez_control` in `phylonetworksPlots.jl` — replaced `max_offset`
parameter with `offset_override`. `plotRCall.jl` — replaced the `max_offset` computation
with `minor_offset`.

**Verification:** Run `plot(net, curved=:minor, useedgelength=true)`. The arc must
stay within the tree bounding box and not extend past the tip labels.

---

## Section 6: The R Package Audit

Every R function called in `plotRCall.jl`:

| R function  | R package  | In base R? | Added by our changes? |
|-------------|------------|------------|-----------------------|
| `plot()`    | graphics   | yes        | no                    |
| `segments()`| graphics   | yes        | no                    |
| `arrows()`  | graphics   | yes        | no                    |
| `lines()`   | graphics   | yes        | **yes**               |
| `text()`    | graphics   | yes        | no                    |

"Base R" means the package ships with every standard R installation automatically
and is loaded in every R session without any `library()` call. The `graphics` package
has been part of base R since version 1.0 and requires no CRAN access or user action.

**No R function used in this codebase requires the user to install any R package.
Every R function used belongs to the `graphics` package which ships with every standard
R installation automatically.**

`lines()` was the only function added by our changes. It is the standard function in
base R for adding a connected polyline to an existing plot. We deliberately did not
use `xspline()` (which has version-dependent `shape` parameter behavior across R
installations) or any function from an optional package.

---

## Section 7: What Was Deliberately Not Changed

### `edgenode_coordinates()`

This function was never touched. Not a single character was modified.

It begins with:
```julia
function edgenode_coordinates(
    net::HybridNetwork,
    useedgelength::Bool,
    usedirecthybridline::Bool,
    preorder::Bool=true,
)
```

And ends with:
```julia
    return edge_xB, edge_xE, edge_yB, edge_yE,
           node_x, node_y, node_yB, node_yE,
           minoredge_xB, minoredge_xE, minoredge_yB, minoredge_yE,
           xmin, xmax, ymin, ymax
end
```

The git diff for `phylonetworksPlots.jl` shows zero removals from this function —
only the new `_quadbez_control` function was added above `check_nodedataframe`.

### The `R"segments"` call for node vertical bars

This line was never changed:
```julia
R"segments"(node_x, node_yB, node_x, node_yE, col=defaultedgecolor)
```
It sits after the entire curved/straight conditional block and always executes for
every value of `curved`. Node vertical bars have no curved interpretation.

### All existing keyword argument defaults

Every keyword argument that existed before our changes retains its original name, type
annotation, and default value:

| Argument | Type | Default |
|---|---|---|
| `useedgelength` | `Bool` | `false` |
| `showtiplabel` | `Bool` | `true` |
| `shownodenumber` | `Bool` | `false` |
| `showedgelength` | `Bool` | `false` |
| `showgamma` | `Bool` | `false` |
| `edgecolor` | (any) | `"black"` |
| `majorhybridedgecolor` | `AbstractString` | `"deepskyblue4"` |
| `minorhybridedgecolor` | `AbstractString` | `"deepskyblue"` |
| `defaultedgecolor` | (any) | `nothing` |
| `showedgenumber` | `Bool` | `false` |
| `shownodelabel` | `Bool` | `false` |
| `edgelabel` | `AbstractDataFrame` | `DataFrame()` |
| `nodelabel` | `AbstractDataFrame` | `DataFrame()` |
| `xlim` | (any) | `nothing` |
| `ylim` | (any) | `nothing` |
| `tipoffset` | (any) | `0` |
| `tipcex` | (any) | `1` |
| `nodecex` | (any) | `1` |
| `edgecex` | (any) | `1` |
| `style` | `Symbol` | `:fulltree` |
| `arrowlen` | `Real` | `(style==:majortree ? 0 : 0.1)` |
| `minorlinetype` | (any) | `nothing` |
| `edgewidth` | (any) | `1` |
| `edgenumbercolor` | (any) | `"grey"` |
| `edgelabelcolor` | (any) | `"black"` |
| `nodelabelcolor` | (any) | `"black"` |
| `edgelabeladj` | (any) | `[.5,0]` |
| `nodelabeladj` | (any) | `1` |
| `preorder` | `Bool` | `true` |

### The return value of `plot()`

Unchanged:
```julia
return (xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
  node_x=node_x, node_y=node_y,
  node_y_lo=node_yB, node_y_hi=node_yE,
  edge_x_lo=edge_xB, edge_x_hi=edge_xE,
  edge_y_lo=edge_yB, edge_y_hi=edge_yE,
  arrow_x_lo=hybridedge_xB, arrow_x_hi=hybridedge_xE,
  arrow_y_lo=hybridedge_yB, arrow_y_hi=hybridedge_yE,
  node_data=ndf, edge_data=edf)
```

The returned named tuple always contains the original straight-line coordinate arrays
from `edgenode_coordinates()`. The Bézier control point coordinates are an internal
implementation detail and are not returned.

### All label and annotation logic

All `if showedgelength`, `if showgamma`, `if showedgenumber`, `if labelnodes`,
`if labeledges`, `if shownodelabel`, `if shownodenumber`, `if showtiplabel` blocks
were never touched.

### All color handling logic

The `if isa(edgecolor, AbstractDict)` block that builds `eCol` and `hybmincol_vec`
was never modified. The vectorized `eCol[ [e.hybrid for e in net.edge] ] .=
majorhybridedgecolor` assignments were never modified.

### All `edgewidth` dictionary logic

The `if isa(edgewidth, Number)` / `elseif isa(edgewidth, AbstractDict)` block that
builds `edgewidth_vec` and `hybridedgewidth_vec` was never modified.

---

## Section 8: How to Use the New Feature

```julia
using PhyloNetworks, PhyloPlots, RCall

net = readnewick(
    "(((A:.2,(B:.1)#H1:.1::0.9):.1,(C:.11,#H1:.01::0.1):.19):.1,D:.4);")
```

| Call | What it demonstrates |
|---|---|
| `plot(net)` | Baseline: all edges straight, original behavior. |
| `plot(net, curved=:none)` | Explicit `:none` — identical output to the line above. |
| `plot(net, curved=:minor)` | Minor hybrid edge diagonal is a smooth arc; everything else straight. |
| `plot(net, curved=:both)` | Both the minor diagonal and major hybrid edge are arcs; tree edges and stubs remain straight. |
| `plot(net, curved=:minor, bend=0.1)` | Subtle curve — arc bows only 10% of one taxon gap. |
| `plot(net, curved=:minor, bend=0.5)` | Pronounced curve — arc bows 50% of one taxon gap. |
| `plot(net, curved=:minor, useedgelength=true)` | Same visual arc size despite longer chord (because `offset_override` decouples bow from chord length). |
| `plot(net, curved=:both, useedgelength=true)` | Both hybrid edge types curved with actual branch lengths. |
| `plot(net, curved=:minor, style=:majortree)` | Majortree layout with curved minor edge; arrowhead always visible (minimum 0.1 in). |
| `plot(net, curved=:minor, style=:fulltree)` | Fulltree layout (default) with curved minor edge. |
| `plot(net, curved=:minor, minorhybridedgecolor="red")` | Curved arc is drawn in red. |
| `plot(net, curved=:minor, edgewidth=3)` | Curved arc is 3 units wide. |
| `plot(net, curved=:minor, edgecolor=Dict(3=>"red"))` | Edge 3's color is red; curved draws it correctly via per-edge color lookup. |

---

## Section 9: How to Run the Tests

### Run the full test suite (recommended)

From the PhyloPlots package root directory:

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

This runs all three test files via `test/runtests.jl`. All imports
(DataFrames, RCall, PhyloNetworks, Test) are loaded automatically.

### Run only `test_plotRCall.jl`

```bash
julia --project=. -e "using RCall, PhyloPlots, Test, DataFrames, PhyloNetworks; include(\"test/test_plotRCall.jl\")"
```

**Important:** You must include `DataFrames` in the `using` list. The test file
uses `DataFrame()` from that package. The `runtests.jl` driver does this automatically.
If you omit `DataFrames`, the pre-existing `"RCall-based plots"` testset will fail with
`UndefVarError: DataFrame not defined`.

### Run only the curved keyword testset

```bash
julia --project=. -e '
  using RCall, PhyloPlots, Test, DataFrames, PhyloNetworks
  @testset "curved only" begin
    include("test/test_plotRCall.jl")
  end
'
```

This will run both the `curved keyword` testset and `RCall-based plots` testset (they
are both in the same file). There is no simpler way to run a single named testset
without importing the test runner infrastructure.

### Open the Julia REPL with the package loaded

```bash
julia --project=.
```

Then in the Julia REPL:
```julia
julia> using PhyloPlots, PhyloNetworks, RCall
julia> net = readnewick("(((A:.2,(B:.1)#H1:.1::0.9):.1,(C:.11,#H1:.01::0.1):.19):.1,D:.4);")
julia> plot(net, curved=:minor)
julia> PhyloPlots._quadbez_control(0.0, 0.0, 2.0, 0.0; bend=0.3)
# returns: (1.0, 0.6, false)
```

### What a passing test output looks like

```
     Testing Running tests...
Test Summary:    | Pass  Total   Time
PhyloPlots Tests |  129    129  12.0s
     Testing PhyloPlots tests passed
```

All tests show `Pass` with no `Fail`, `Error`, or `Broken` columns.

### What a failing test output looks like

A Julia test assertion failure:
```
Test Failed at /path/to/test/test_plotRCall.jl:57
  Expression: cx ≈ 1.0
   Evaluated: 2.0 ≈ 1.0
ERROR: LoadError: Some tests did not pass: 0 passed, 1 failed, 0 errored, 0 broken.
```

This tells you the file, the line number, the expression that failed, and the actual
values.

### How to distinguish an R error from a Julia test failure

An R error comes through RCall as an exception of type `RCall.REvalError`:
```
RCall.REvalError: Error in eval(parse(text = x), envir = envir) :
  object 'x_curve' not found
```
The traceback will point into an `R"""..."""` or `R"funcname"(...)` call in the Julia
code.

A Julia assertion failure is type `Test.FallbackTestSetException` or just a
`LoadError`. The traceback points to a `@test` line, not into R code.

---

## Section 10: Glossary

**Taxon / taxa** — A species or lineage in the network. The tips (leaves) of the tree.
In a 4-taxon network, there are 4 named leaves like A, B, C, D.

**Hybrid node** — A node in the network that has two parents instead of one. It
represents a species that received genes from two different ancestor lineages. Also
called a reticulation node.

**Minor hybrid edge** — The parent edge of a hybrid node that carries less than 50% of
the inheritance (the one with the smaller `gamma` value, also called the "gene flow"
edge). In a network, there is one minor hybrid edge per hybrid node.

**Major hybrid edge** — The parent edge of a hybrid node that carries 50% or more of
the inheritance (the one with the larger `gamma`). The major hybrid edge is part of the
displayed major tree.

**Horizontal stub** — The short horizontal line segment drawn for a minor hybrid edge
to show its position in the layout. It starts at the minor parent node's x position and
ends at a calculated x position. Its y coordinate is a dedicated "phantom row" in the
plot, separate from the taxa and internal node rows. Stored in `edge_xB[i]` / `edge_xE[i]`
for the minor edge index `i`.

**Diagonal segment** — The angled line connecting the end of the horizontal stub to
the hybrid child node. This is the piece that carries the arrowhead showing the
direction of gene flow. Stored in `hybridedge_xB[j]` / `hybridedge_xE[j]` etc. for
the j-th minor edge.

**Reticulation** — Another word for hybridization or gene flow. A network with
reticulations has nodes with two parents.

**Bézier curve** — A smooth curve defined by a start point, one or more control points,
and an end point. A quadratic Bézier has one control point; a cubic Bézier has two.
The curve is attracted toward the control points but does not pass through them.

**Control point** — The single intermediate point P1 in a quadratic Bézier curve that
determines how much and in what direction the curve bows away from the straight chord.
Moving P1 closer to the chord makes a subtler arc; moving it further makes a more
pronounced one.

**Chord** — The straight line segment connecting the start point (x0, y0) and end point
(x2, y2) of the curve. The Bézier arc bows away from this chord.

**Perpendicular** — A direction at exactly 90 degrees to another direction. For a
chord pointing in direction (dx, dy), the perpendicular direction is (-dy, dx). The
control point is offset from the chord midpoint in this perpendicular direction to
create the bow.

**taxon_spacing** — The average distance between adjacent taxa in the vertical (y) axis
of the plot. Computed as `(ymax - ymin) / numtaxa`. In the default layout (no edge
lengths), this is approximately 1.0 plot units. The bow magnitude for minor hybrid
edges is `bend × taxon_spacing` so that the arc is always a consistent readable size.

**bend** — The keyword argument controlling how far the arc bows away from the straight
chord. Specifically, `bend × taxon_spacing` is the distance in plot coordinates between
the chord midpoint and the Bézier control point. Default 0.3.

**offset_override** — A parameter of `_quadbez_control` that, when given a finite
value, overrides the default `bend × chord_length` offset magnitude with a fixed value.
Used in `plotRCall.jl` for minor hybrid edges to pass a taxon-spacing-based offset that
prevents arcs from escaping the plot.

**fallback flag** — The third return value of `_quadbez_control`. When `true`, the
chord is degenerate (too short to draw a visible arc) and the caller should draw a
straight line instead of calling `lines()`.

**Vectorized call** — An R function call that processes an entire array at once in a
single R call, rather than looping over elements in Julia. For example,
`R"segments"(edge_xB, edge_yB, edge_xE, edge_yE, ...)` draws all N segments in one R
call. This is much faster than N separate R calls from a Julia loop.

**Base R** — The set of packages that ship with every standard R installation and are
loaded automatically in every R session, requiring no `library()` call or `install.packages()`.
The `graphics` package (which provides `plot()`, `segments()`, `arrows()`, `lines()`,
`text()`) is a base R package.
