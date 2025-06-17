# Basic Usage
Here we walk through an example of generating a light curve for a 3-body system.

##### Units
A quick note on units.
- Distance: AU
- Time: Days
- Mass: Solar Masses
- Angles: Radians

### Transit Model

We start by defining the stellar/transit parameters:
```@example 1
using Photodynamics

rstar = 0.00465047 * 0.1192 # Radius of the star (AU)
k = [0.085, 0.083]  # Planet-star radius ratios
u_n = [0.28, 0.11]  # Limbdarkening coefficients
nothing # hide
```

Next, we set up the simulated light curve. We assume uniform cadence over length of the observations:
```@example 1
exp_time = 2 / 60 / 24 # 2 min cadence in days
obs_duration = 2.0 # 2 days total observations
lc = Lightcurve(exp_time, obs_duration, u_n, k, rstar)
```

### Dynamical Model
We re-export the `NbodyGradient.jl` API for initial conditions. See the docs for more details ([NbodyGradient.jl docs](https://ericagol.github.io/NbodyGradient.jl/dev/))

Start by defining the dynamical model initial conditions:
```@example 1
star = Elements(m = 1.0) # Mass of the star (solar masses)
planet_b = Elements(
    m = 3e-5,      # Planet mass (solar masses)
    P = 1.5,       # Period (days)
    t0 = 0.32,     # Initial time of transit (days from start of simulation)
    ecosω = 0.03,  # Eccentricity vector component (eccentricity * cos(argument of periastron))
    I = π/2)       # Inclination (radians; edge-on)
planet_c = Elements(
    m = 6e-5, 
    P = 24,
    t0 = 0.35,
    ecosω = 0.02,
    I = π/2)

t0 = 0.0 # Time of initial conditions 
N = 3 # Number of bodies (or a hierarchy matrix; see NbodyGradient.jl docs)
ic = ElementsIC(t0, N, star, planet_b, planet_c)
nothing # hide
```

Set up the `NbodyGradient.jl` integrator:
```@example 1
h = 0.05 # Time step [days]
intr = Integrator(h, obs_duration)
nothing # hide 
```

Now we integrate the N-body system, record the transit times, and compute the PK20 expansion points about each transit:
```@example 1
s = State(ic) # Get initial Cartesian coordinates from orbital elements
tt = TransitTiming(obs_duration, ic) # Pre-allocate transit timing arrays
ts = TransitSeries(obs_duration, ic) # Pre-allocate expansion point arrays
intr(s, ts, tt) # Dispatch on the desired outputs 
```

### Generating a light curve

Finally, we combine the photometric and dynamical model output and compute a light curve:
```@example 1
compute_lightcurve!(lc, ts)
nothing # hide
```

Plot a section of the light curve:
```@example 1
using Plots
mask = 0.3 .< lc.tobs .< 0.375
plot!(lc.tobs[mask], lc.flux[mask] .+ 1, legend=false, marker=true)
xlabel!("Time [Days]"); ylabel!("Relative Flux")
```

### Modifying light curve parameters
Since the transit and dynamical models are computed independently, we can modify the light curve parameters without needing to re-compute the N-body integration. 

As an example, let's say we want to better resolve the ingress/egress. We can simply decrease the exposure times and recompute the light curve:
```@example 1
exp_time = exp_time / 20 # reduce the exposure times
lc = Lightcurve(exp_time, obs_duration, u_n, k, rstar)
compute_lightcurve!(lc, ts)

mask = 0.3 .< lc.tobs .< 0.375
plot(lc.tobs[mask], lc.flux[mask] .+ 1, legend=false, marker=true)
xlabel!("Time [Days]"); ylabel!("Relative Flux")
```

We can also make changes to the dynamical parameters. Let's change the inclination of one of the planets:
```@example 1
planet_c = Elements(
    m = 6e-5, 
    P = 24,
    t0 = 0.35,
    ecosω = 0.02,
    I = π/2+0.002) # Vary the inclination
ic = ElementsIC(t0, N, star, planet_b, planet_c) # Re-initialize
s = State(ic) 
tt = TransitTiming(obs_duration, ic) 
ts = TransitSeries(obs_duration, ic) 
intr(s, ts, tt) 

# Plot original and new light curve
mask = 0.3 .< lc.tobs .< 0.375
plot(lc.tobs[mask], lc.flux[mask] .+ 1, lw=3, label="Original")
xlabel!("Time [Days]"); ylabel!("Relative Flux")

# Compute the new light curve and plot
compute_lightcurve!(lc, ts)
plot!(lc.tobs[mask], lc.flux[mask] .+ 1, lw=2, label="Modified")
```

### Multi-band light curves
We can model photometry from multiple instruments with different wavelength bands, which will have different limbdarkening coefficients. We use the same N-body model for both:
```@example 1
u_n_1 = copy(u_n)
u_n_2 = [0.1, 0.22]
lc_1 = Lightcurve(exp_time, obs_duration, u_n_1, k, rstar)
lc_2 = Lightcurve(exp_time, obs_duration, u_n_2, k, rstar)
compute_lightcurve!.([lc_1, lc_2], ts) # Vectorize over the light curve objects 

mask_1 = 0.3 .< lc_1.tobs .< 0.375
mask_2 = 0.3 .< lc_2.tobs .< 0.375
plot(lc_1.tobs[mask_1], lc_1.flux[mask_1] .+ 1, lw=3, label="Original")
plot!(lc_2.tobs[mask_2], lc_2.flux[mask_2] .+ 1, lw=2, label="Modified")
xlabel!("Time [Days]"); ylabel!("Relative Flux")
```

### Derivative Lightcurves
To generate derivatives of the flux with respect to the model parameters we first need to
get the derivatives for the dynamical model. We simply re-initialize the N-body simulation
and run the integration with the `grad=true` keyword.
```@example 1
s = State(ic)
tt = TransitTiming(obs_duration, ic)
ts = TransitSeries(obs_duration, ic)
intr(s, ts, tt; grad=true) 
``` 

Then, we get a differentiable light curve (`dLightcurve` tells `compute_lightcurve!` to compute derivatives):
```@example 1
dlc = dLightcurve(exp_time, obs_duration, u_n, k, rstar)
compute_lightcurve!(dlc, ts)
```

We now have the derivatives of the light curve with respect to the initial stellar parameters, masses, and initial *Cartesian* coordinates. 
If we want the derivatives with respect to the orbital elements, we run:
```@example 1
transform_to_elements!(s, dlc)
nothing # hide
```

Now we visualize the same transits, but plot the derivatives of the flux with respect to the inital mass of the planets. The `dLightcurve` type has a field `dfdelements` which is a `Matrix` of dimension number-of-exposure-times X orbital-element-index. The orbital elements are ordered (P, t0, ecosω, esinω, I, Ω, m). So, to get the mass of planet b, we want orbital element index `14` (for flexibility, the indices include the "orbital elements of the star". So, the planets start at index `8`). 
```@example 1
mask = 0.3 .< dlc.tobs .< 0.375
p1 = plot(dlc.tobs[mask], dlc.dfdelements[mask,14], legend=false, lw=3)
ylabel!("dF / dm_b")
p2 = plot(dlc.tobs[mask], dlc.dfdelements[mask,21], legend=false, lw=3, color=2)
xlabel!("Time [Days]"); ylabel!("dF / dm_c")
plot(p1,p2,layout=(2,1))
```

Similarly, we can access the other parameter derivative light curves. Below shows these fields and how the arrays are indexed:
```julia
dlc.dfdr  # Stellar radius [times]
dlc.dfdu  # Limbdark coefficients [times X u_n]
dlc.dfdk  # Radius ratios [times X k]
dlc.dfdq0 # Cartesian coordinates [times X cartesian coordinates]
```