using Documenter
using NMRHamiltonian

makedocs(
    sitename = "NMRHamiltonian",
    format = Documenter.HTML(),
    modules = [NMRHamiltonian]
)

makedocs(
    sitename="NMRHamiltonian.jl",
    modules=[NMRHamiltonian],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
    pages=[
        "Overview" => "index.md",
        "Public API" => "api.md",
        "Demo: code walk-through" => "demo_code.md",
        "Demo: return variables" => "return_vars.md",
    ],
    #strict=true,
)

deploydocs(
    repo = "github.com/AI4DBiological-Systems/NMRHamiltonian.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)