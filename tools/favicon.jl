# A single convection plume in a square box heated from below and cooled from above: the hot fluid rises near the left
# wall and its head curls over to the right, like the hooks in the banner. The same model and picture as
# tools/plumes.jl, on a coarser mesh, for the favicon of schmaeke.dev.
#
#     julia --project=$HOME/Software/MicroFE.jl tools/favicon.jl
#
# Boussinesq equations in the unit square with no-slip walls, a hot floor (T = 1), a cold lid (T = -1) and insulated side
# walls (see tools/plumes.jl). With ν = κ = 3·10⁻³ the Rayleigh number is 2·10⁵. The initial temperature is 0 plus a
# small perturbation that is warmest at x = 0.3, so a single plume rises left of the middle; the wall on its left side
# makes it lean and curl to the right. Output: favicon.svg next to index.html, the temperature at t = 10 sampled into a
# 128 × 128 WebP, with the mesh edges drawn over it. It is generated: change this script, not the SVG. Needs
# ImageMagick's `magick`; takes about 20 s.
using MicroFE, StaticArrays, LinearAlgebra, Base64

out = joinpath(@__DIR__, "..", "favicon.svg")

const ν, κ = 3e-3, 3e-3                                         # viscosity, diffusivity
const Δt, steps = 0.02, 500                                     # the icon shows t = 10
const e₂ = SA[0.0, 1.0]

# velocity (no slip), pressure (fixed at the corner (0, 0)) and temperature (fixed on the floor and the lid)
fields(mesh) = (Field(mesh, P1, 2; dirichlet = Boundary()), Field(mesh, P1; dirichlet = x -> x == SA[0.0, 0.0]),
                Field(mesh, P1; dirichlet = ["x2-", "x2+"]))
walls!(uh) = interpolate!(uh, x -> (SA[0.0, 0.0], 0.0, x[2] < 1 / 2 ? 1.0 : -1.0); dofs = fixeddofs(uh.fields))

# one backward Euler step, linearised about the previous velocity w, with SUPG/PSPG stabilisation (see tools/plumes.jl)
τ(w, h, d) = 1 / sqrt((2 / Δt)^2 + (2norm(w) / h)^2 + (4d / h^2)^2)
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
# which it changes by more than 0.16, down to a diameter of 0.07 (two bisections below the initial cells): coarse
# enough for the triangles to show at icon size. uh[3].u[c] are the temperatures at the vertices of cell c.
variation(uh) = [maximum(uh[3].u[c]) - minimum(uh[3].u[c]) for c in uh.fields[1].mesh.cells]
function adapt(uh)
    uh = walls!(transfer(uh, last(coarsen(uh.fields[1].mesh, variation(uh) .< 0.04))))
    mesh = uh.fields[1].mesh
    diameter(c) = maximum(norm(mesh.x[c[i]] - mesh.x[c[j]]) for (i, j) in ((1, 2), (2, 3), (1, 3)))
    marked = findall(variation(uh) .> 0.16 .&& diameter.(mesh.cells) .> 0.07)
    isempty(marked) ? uh : walls!(transfer(uh, last(refine(mesh, marked))))
end

# initial state: at rest, T = 0 plus a perturbation of amplitude 0.02, a bump around x = 0.3 on a slightly cool
# background; the initial 8 × 8 mesh is adapted to it a few times
T₀(x) = 0.02 * (exp(-((x[1] - 0.3) / 0.12)^2) - 0.3) * sin(π * x[2])
start!(uh) = walls!(interpolate!(uh, x -> (SA[0.0, 0.0], 0.0, T₀(x))))
uh = start!(FEFunction(fields(cubemesh(range(0, 1, 9), range(0, 1, 9)))))
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

# The picture, 64 × 64 SVG units: the temperature at the pixel centres of a 128 × 128 image, in the colour map of
# tools/plumes.jl, as an embedded WebP; the mesh edges over it as white lines.
const S, n = 64, 128                                            # SVG units per length unit, pixels
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
T = evaluate(uh[3], [SA[(i - 0.5) / n, 1 - (j - 0.5) / n] for i in 1:n for j in 1:n]; outside = 0.0)
ppm, webp = joinpath.(mktempdir(), ("favicon.ppm", "favicon.webp"))
open(ppm, "w") do io                                            # binary PPM: a header, then RGB row by row
    write(io, "P6\n$n $n\n255\n")
    foreach(t -> write(io, colour(t)), permutedims(reshape(T, n, n)))
end
run(`magick $ppm -quality 90 $webp`)

edges = unique(sort([minmax(c[i], c[j]) for c in mesh.cells for (i, j) in ((1, 2), (2, 3), (1, 3))]))
point(x) = string(round(S * x[1]; digits = 2), " ", round(S * (1 - x[2]); digits = 2))
image = "data:image/webp;base64," * base64encode(read(webp))
lines = join("M" * point(mesh.x[a]) * "L" * point(mesh.x[b]) for (a, b) in edges)
write(out, """
    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $S $S">
    <image width="$S" height="$S" preserveAspectRatio="none" href="$image"/>
    <path fill="none" stroke="#fff" stroke-opacity="0.6" stroke-width="0.35" stroke-linecap="round" d="$lines"/>
    </svg>
    """)
println("wrote ", normpath(out), ": ", length(mesh.cells), " cells")
