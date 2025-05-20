```@meta
CurrentModule = Photodynamics
```

# Photodynamics
[![GitHub](https://img.shields.io/badge/GitHub-grey?logo=github)](https://github.com/langfzac/photodynamics.jl)
[![license](https://img.shields.io/badge/license-MIT-green)](https://github.com/langfzac/Photodynamics.jl/blob/main/LICENSE)
[![arxiv](https://img.shields.io/badge/arXiv-2410.03874-orange?logo=arxiv&logoColor=%23B31B1B)](https://arxiv.org/abs/2410.03874)

A differentiable photodynamical model for multiplanet transit light curves, described in [Langford & Agol 2025](https://arxiv.org/abs/2410.03874). 

Photodynamical models are composed of a dynamical model and a transit (i.e. photometric) model. The currently implemented methods are:
- Dynamical Models:
    - [NbodyGradient.jl](https://github.com/ericagol/NbodyGradient.jl) + PK20: Differentiable 4th-order symplectic N-body integrator combined with an implementation of [Parviainen & Korth 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.3356P/abstract) (PK20). 
- Transit Models:
    - [Limbdark.jl](https://github.com/rodluger/limbdark.jl): Differentiable quadratic limbdarkened transits.

## Getting Started
First, you'll need to add the package. Using the Julia REPL, type `]` to enter the package manager. 

To add the latest tagged release, run:
```julia
pkg> add Photodynamics
```

If you'd like to use a developement version of the code, add `#branch-name`:
```julia
pkg> add Photodynamics#branch-name
```

See the [Tutorials](@ref) page for basic usage.
