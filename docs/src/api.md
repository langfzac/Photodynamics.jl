# API

## Functions
```@docs
compute_lightcurve!
```

## Types

### `TransitSeries`
```@docs
TransitSeries
```

### `Integrator`
```@docs
Integrator
```

### `Lightcurve`
```@docs
Lightcurve
```
#### Constructors
```@docs
Lightcurve(dt::T, duration::T, u_n::Vector{T}, k::Vector{T}, rstar::T) where T<:Real
```

### `dLightcurve`
```@docs
dLightcurve
```
#### Constructors
```@docs
dLightcurve(dt::T, duration::T, u_n::Vector{T}, k::Vector{T}, rstar::T) where T<:Real
```
