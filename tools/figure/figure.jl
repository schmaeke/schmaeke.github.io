# Header figure of index.html, computed with Unfitted.jl.
#
# A plate with a five-lobed hole, Ω = {φ ≤ 0}, is immersed in a 6 × 6 Cartesian
# base mesh. One overlay, four times finer, is switched on only in a band that
# tracks ∂Ω. Cells outside Ω are dropped; every cut integration region carries
# the non-negative moment-fitted (NNMF) quadrature rule the package assembles
# with. The figure draws the live base cells, the live overlay cells, ∂Ω and
# every quadrature point of the regions that couple base and overlay (tensor
# Gauss where full, NNMF where cut) as squares whose area is proportional to
# their weight.
#
#     julia --project=tools/figure -e 'using Pkg; Pkg.instantiate()'
#     julia --project=tools/figure tools/figure/figure.jl
#
# The SVG replaces everything between the figure markers in index.html.

using Unfitted

# ── Geometry ──────────────────────────────────────────────────────────────────

const CENTER = (0.5, 0.5)
const R0, AMPLITUDE, LOBES, PHASE = 0.335, 0.13, 5, deg2rad(35)

radius_at(θ) = R0 * (1 + AMPLITUDE * cos(LOBES * (θ - PHASE)))

function φ(x)
    dx, dy = x[1] - CENTER[1], x[2] - CENTER[2]
    return radius_at(atan(dy, dx)) - sqrt(dx^2 + dy^2 + 1e-30)     # Ω lies outside the lobes
end

# ── Discretisation ────────────────────────────────────────────────────────────

const NBASE, FACTOR, ORDER = 6, 4, 1
const NFINE = NBASE * FACTOR
const H = 1 / NFINE

cellbox(i, j, h) = box(((i - 1) * h, (j - 1) * h), (i * h, j * h))

# Overlay cells that come within one fine cell of ∂Ω.
function near_boundary(b, _)
    lo, hi = b.lower, b.upper
    values = (φ(lo), φ((hi[1], lo[2])), φ((lo[1], hi[2])), φ(hi), φ((lo .+ hi) ./ 2))
    return minimum(values) < H && maximum(values) > -H
end

omega = box((0.0, 0.0), (1.0, 1.0))
domain = physical_domain(leaf(φ); subcell_length_scale=H / 16, max_depth=4)
V = overlay(space(omega; cells=NBASE, order=ORDER, physical=domain), omega; cells=NFINE,
            active=near_boundary)
plan = Unfitted.integration_plan(V)

# Quadrature points of the coupling regions: those covered by both levels.
points = Tuple{Float64,Float64,Float64}[]
for region in plan.regions
    any(p -> p.level == 2, region.parents) || continue
    jacobian = Unfitted._region_jacobian(region)
    for (η, w) in zip(region.quadrature.points, region.quadrature.weights)
        x = Unfitted.reference_to_physical(region.box, η)
        push!(points, (x[1], x[2], w * jacobian))
    end
end

live(b) = classify_cell(domain, b) !== :fictitious
base_live = [(i, j) for i in 1:NBASE, j in 1:NBASE if live(cellbox(i, j, 1 / NBASE))]
base_void = [(i, j) for i in 1:NBASE, j in 1:NBASE if !live(cellbox(i, j, 1 / NBASE))]
fine_live = [(i, j) for i in 1:NFINE, j in 1:NFINE
             if near_boundary(cellbox(i, j, H), nothing) && live(cellbox(i, j, H))]

# ── SVG ───────────────────────────────────────────────────────────────────────

const S = 320.0                                  # viewBox edge
const UNIT = S / NFINE                           # one fine cell in viewBox units

num(v) = (s = string(round(v; digits=1) + 0.0); endswith(s, ".0") ? s[1:end-2] : s)

# Unit edges on the fine lattice, as (orientation, line, position).
function cell_edges(cells, size)
    edges = Set{Tuple{Char,Int,Int}}()
    for (i, j) in cells, k in 0:size-1
        i0, j0 = (i - 1) * size, (j - 1) * size
        push!(edges, ('h', j0, i0 + k), ('h', j0 + size, i0 + k), ('v', i0, j0 + k),
              ('v', i0 + size, j0 + k))
    end
    return edges
end

# Merge collinear unit edges into maximal runs; optionally skip the frame.
function edge_path(edges; skip_frame=false)
    lines = Dict{Tuple{Char,Int},Vector{Int}}()
    for (o, a, b) in edges
        skip_frame && a in (0, NFINE) && continue
        push!(get!(lines, (o, a), Int[]), b)
    end
    io = IOBuffer()
    for ((o, a), bs) in sort!(collect(lines))
        sort!(bs)
        start = prev = bs[1]
        for b in [bs[2:end]; typemax(Int)]
            b == prev + 1 && (prev = b; continue)
            if o == 'h'
                print(io, "M", num(start * UNIT), " ", num(S - a * UNIT), "H", num((prev + 1) * UNIT))
            else
                print(io, "M", num(a * UNIT), " ", num(S - start * UNIT), "V", num(S - (prev + 1) * UNIT))
            end
            start = prev = b
        end
    end
    return String(take!(io))
end

base_edges = cell_edges(base_live, FACTOR)
void_edges = setdiff(cell_edges(base_void, FACTOR), base_edges)
fine_edges = setdiff(cell_edges(fine_live, 1), base_edges)

boundary = "M" * join((let r = radius_at(θ)
                           num(S * (CENTER[1] + r * cos(θ))) * " " * num(S * (1 - CENTER[2] - r * sin(θ)))
                       end for θ in range(0, 2π; length=361)[1:end-1]), "L") * "Z"

wmax = maximum(last, points)
squares = join(begin
                   s = max(0.8, 3.0 * sqrt(w / wmax))
                   "M$(num(S * x - s / 2)) $(num(S * (1 - y) - s / 2))h$(num(s))v$(num(s))h-$(num(s))z"
               end for (x, y, w) in points)

svg = """
        <svg viewBox="0 0 320 320" role="img" aria-labelledby="fig-title">
          <title id="fig-title">Computed with Unfitted.jl: a Cartesian base mesh cut by a five-lobed hole, one finer overlay active only along the boundary, and the quadrature points where the two couple: tensor Gauss in full cells, non-negative moment-fitted in cut cells.</title>
          <path class="fill" fill-rule="evenodd" d="M0 0H320V320H0Z$boundary"/>
          <path class="void" d="$(edge_path(void_edges; skip_frame=true))"/>
          <path class="fine" d="$(edge_path(fine_edges))"/>
          <path class="mesh" d="$(edge_path(base_edges))"/>
          <path class="boundary" d="$boundary"/>
          <path class="points" d="$squares"/>
          <rect class="frame" width="320" height="320"/>
        </svg>
"""

page = joinpath(@__DIR__, "..", "..", "index.html")
html = read(page, String)
head, rest = split(html, "<!-- figure:start -->\n"; limit=2)
_, tail = split(rest, "        <!-- figure:end -->"; limit=2)
write(page, head * "<!-- figure:start -->\n" * svg * "        <!-- figure:end -->" * tail)

println(length(base_live), " base cells, ", length(fine_live), " overlay cells, ",
        count(r -> any(p -> p.level == 2, r.parents), plan.regions), " coupling regions, ",
        length(points), " quadrature points")
println("max moment-fit residual ", plan.moment_fit_residual_max, ", svg ", sizeof(svg), " bytes")
