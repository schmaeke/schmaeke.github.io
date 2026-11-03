# Rayleigh–Bénard convection in a box heated from below and cooled from above: hot plumes rise, cold ones sink. Equal-
# order P1 elements with SUPG/PSPG stabilisation, backward Euler in time, and a mesh that is refined and coarsened as the
# plumes move. It makes the banner of schmaeke.dev.
#
#     julia --project=$HOME/Software/MicroFE.jl tools/plumes.jl
#
# Boussinesq equations in the box (0, 8) × (0, 2) with no-slip walls, a hot floor (T = 1), a cold lid (T = -1) and
# insulated side walls:
#     ∂u/∂t + (∇u) u + ∇p = ν Δu + T e₂,    div u = 0,    ∂T/∂t + u ⋅ ∇T = κ ΔT.
# With ν = κ = 4·10⁻³ the Rayleigh number is 10⁶. The fluid starts at rest at T = 0 with a little noise, which decides
# where the boundary layers at the floor and the lid break up into plumes. Output: banner.svg next to index.html, the
# temperature at t = 12 sampled per pixel into a WebP, with the mesh edges drawn over it. It is generated: change this
# script, not the SVG. Needs ImageMagick's `magick` for the WebP; takes about a minute on one thread.
using MicroFE, StaticArrays, LinearAlgebra, Base64, Random

out = joinpath(@__DIR__, "..", "banner.svg")

const W, H = 8.0, 2.0                                           # the box, 4:1 like the banner
const ν, κ = 4e-3, 4e-3                                         # viscosity, diffusivity
const Δt, steps = 0.04, 300                                     # the banner shows t = 12
const e₂ = SA[0.0, 1.0]

# velocity (no slip), pressure (fixed at the corner (0, 0): it is determined up to a constant only) and temperature
# (fixed on the floor "x2-" and the lid "x2+"); walls! sets their Dirichlet values
fields(mesh) = (Field(mesh, P1, 2; dirichlet = Boundary()), Field(mesh, P1; dirichlet = x -> x == SA[0.0, 0.0]),
                Field(mesh, P1; dirichlet = ["x2-", "x2+"]))
walls!(uh) = interpolate!(uh, x -> (SA[0.0, 0.0], 0.0, x[2] < H / 2 ? 1.0 : -1.0); dofs = fixeddofs(uh.fields))

# stabilisation parameter from the advecting velocity w, the cell size h and the viscosity or diffusivity d
τ(w, h, d) = 1 / sqrt((2 / Δt)^2 + (2norm(w) / h)^2 + (4d / h^2)^2)

# One backward Euler step, linearised about the velocity w of the previous step: a bilinear form in (u, p, T), tested
# with (v, q, s), and a linear form with the previous state.
#   Galerkin:   u ⋅ v / Δt + (∇u w) ⋅ v + ν ∇u : ∇v - p div v - q div u - T e₂ ⋅ v + T s / Δt + (w ⋅ ∇T) s + κ ∇T ⋅ ∇s
#   SUPG/PSPG:  the strong momentum residual (u / Δt + ∇u w + ∇p - T e₂) times τ (∇v w - ∇q), and the strong temperature
#               residual (T / Δt + w ⋅ ∇T) times τ w ⋅ ∇s; ν Δu and κ ΔT vanish for P1
function lhs((w, _, _), (u, p, T), (v, q, s), (; h))
    τm, τt = τ(w, h, ν), τ(w, h, κ)
    momentum = u / Δt + ∇(u) * w + ∇(p) - T * e₂
    temperature = T / Δt + w ⋅ ∇(T)
    u ⋅ v / Δt + (∇(u) * w) ⋅ v + ν * ∇(u) ⋅ ∇(v) - p * div(v) - q * div(u) - T * (e₂ ⋅ v) +
        τm * momentum ⋅ (∇(v) * w - ∇(q)) + T * s / Δt + (w ⋅ ∇(T)) * s + κ * ∇(T) ⋅ ∇(s) + τt * temperature * (w ⋅ ∇(s))
end
rhs((w, _, T), (v, q, s), (; h)) =
    (w / Δt) ⋅ (v + τ(w, h, ν) * (∇(v) * w - ∇(q))) + T / Δt * (s + τ(w, h, κ) * (w ⋅ ∇(s)))

# Adaptivity on the temperature: coarsen the cells across which it changes by less than 0.04, then refine those across
# which it changes by more than 0.16, down to a diameter of 0.046 (six bisections below the initial cells). CG P1 dofs
# are numbered by vertex, so uh[3].u[c] are the temperatures at the vertices of cell c.
variation(uh) = [maximum(uh[3].u[c]) - minimum(uh[3].u[c]) for c in uh.fields[1].mesh.cells]
function adapt(uh)
    uh = walls!(transfer(uh, last(coarsen(uh.fields[1].mesh, variation(uh) .< 0.04))))
    mesh = uh.fields[1].mesh
    diameter(c) = maximum(norm(mesh.x[c[i]] - mesh.x[c[j]]) for (i, j) in ((1, 2), (2, 3), (1, 3)))
    marked = findall(variation(uh) .> 0.16 .&& diameter.(mesh.cells) .> 0.046)
    isempty(marked) ? uh : walls!(transfer(uh, last(refine(mesh, marked))))
end

# initial state: at rest, T = 0 plus modes in x with random amplitudes (standard deviation 0.005) and phases, vanishing
# at the floor and the lid; the initial 32 × 8 mesh is adapted to it a few times
const noise = (rng = MersenneTwister(1); [(k, 0.005 * randn(rng), 2π * rand(rng)) for k in 4:40])
T₀(x) = sum(a * sin(k * π * x[1] / W + φ) for (k, a, φ) in noise) * sin(π * x[2] / H)
start!(uh) = walls!(interpolate!(uh, x -> (SA[0.0, 0.0], 0.0, T₀(x))))
uh = start!(FEFunction(fields(cubemesh(range(0, W, 33), range(0, H, 9)))))
for _ in 1:6
    global uh = start!(adapt(uh))
end

# time steps; the mesh is adapted after every second one
for step in 1:steps
    X = uh.fields
    A = assemble(lhs, uh, X, X; degree = 2)
    b = assemble(rhs, uh, X; degree = 2)
    global uh = solve!(walls!(FEFunction(X)), A, b)
    step % 2 == 0 && (global uh = adapt(uh))
end

# The picture, 1200 × 300 SVG units: the temperature at the pixel centres of a 2400 × 600 image (linear within each cell,
# as P1), coloured by a rainbow of Tango colours, linear in between: Sky Blue (from dark to light) for cold, Chameleon in
# the middle, Butter, Orange, Scarlet Red (from light to dark) and dark Plum for hot, so both ends darken. The stops crowd
# towards the middle, where most of the fluid is, so the plumes stand out. ImageMagick turns the image into a WebP, which
# is embedded in the SVG; the mesh edges, each once, are drawn over it as white lines.
const S, nx, ny = 150, 2400, 600                                # SVG units per length unit, pixels
const stops = [(-1.0, "#204a87"), (-0.45, "#3465a4"), (-0.18, "#729fcf"), (0.0, "#73d216"),
               (0.15, "#edd400"), (0.3, "#f57900"), (0.45, "#cc0000"), (0.65, "#a40000"), (1.0, "#5c3566")]
rgb(c) = SA[parse(Int, c[2:3]; base = 16), parse(Int, c[4:5]; base = 16), parse(Int, c[6:7]; base = 16)] / 255
function colour(t)
    k = clamp(something(findfirst(s -> s[1] ≥ t, stops), length(stops)), 2, length(stops))
    (a, ca), (b, cb) = stops[k-1], stops[k]
    s = clamp((t - a) / (b - a), 0, 1)
    round.(UInt8, 255 * ((1 - s) * rgb(ca) + s * rgb(cb)))
end

mesh = uh.fields[1].mesh
T = evaluate(uh[3], [SA[(i - 0.5) / nx * W, H - (j - 0.5) / ny * H] for i in 1:nx for j in 1:ny]; outside = 0.0)
ppm, webp = joinpath.(mktempdir(), ("banner.ppm", "banner.webp"))
open(ppm, "w") do io                                            # binary PPM: a header, then RGB row by row
    write(io, "P6\n$nx $ny\n255\n")
    foreach(t -> write(io, colour(t)), permutedims(reshape(T, ny, nx)))
end
run(`magick $ppm -quality 88 $webp`)

edges = unique(sort([minmax(c[i], c[j]) for c in mesh.cells for (i, j) in ((1, 2), (2, 3), (1, 3))]))
point(x) = string(round(S * x[1]; digits = 1), " ", round(S * (H - x[2]); digits = 1))
title = "Hot and cold plumes in Rayleigh–Bénard convection on an adaptively refined finite element mesh, " *
        "computed with MicroFE.jl"
image = "data:image/webp;base64," * base64encode(read(webp))
lines = join("M" * point(mesh.x[a]) * "L" * point(mesh.x[b]) for (a, b) in edges)
write(out, """
    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $(Int(S * W)) $(Int(S * H))">
    <title>$title</title>
    <image width="$(Int(S * W))" height="$(Int(S * H))" preserveAspectRatio="none" href="$image"/>
    <path fill="none" stroke="#fff" stroke-opacity="0.6" stroke-width="0.8" stroke-linecap="round" d="$lines"/>
    </svg>
    """)
println("wrote ", normpath(out), ": ", length(mesh.cells), " cells")
