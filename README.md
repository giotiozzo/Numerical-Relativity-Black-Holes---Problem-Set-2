# ADM 1+1D Gauge Wave Simulator

## Description

This code is a numerical simulation of Einstein's general relativity in 1+1 dimensions (one spatial dimension + time) using the ADM (Arnowitt-Deser-Misner) formulation. 

### What It Does

The simulator evolves a **gauge wave** through flat spacetime. Unlike gravitational waves in our 3D universe, this wave is purely a coordinate effect—it's essentially watching how different ways of slicing up spacetime create apparent "waves" that aren't physically real.

Think of it like this: if you're filming a still object but moving your camera up and down, the video looks like the object is moving. Similarly, this code shows how choosing different coordinate systems in spacetime creates apparent wave-like behavior, even though the underlying spacetime is completely flat (Minkowski).

### Key Features

- **ADM formulation**: Uses the standard 3+1 split of spacetime used in numerical relativity
- **Bona-Massó slicing**: Implements a stable gauge condition (coordinate choice) that prevents numerical instabilities
- **4th-order Runge-Kutta**: High-accuracy time integration method
- **Periodic boundaries**: Allows the wave to wrap around the domain
- **Pure gauge dynamics**: Demonstrates coordinate effects without physical gravitational radiation

### Physical Significance

In 1+1 dimensions, there are no physical gravitational waves—the Riemann curvature tensor vanishes identically. This makes it an ideal testbed for:

- Testing numerical relativity algorithms
- Understanding gauge freedom in general relativity
- Validating evolution schemes without the complexity of 3+1D
- Learning how coordinate choices affect spacetime evolution

### What the Code Evolves

The simulation tracks the time evolution of:

1. **Lapse function (α)**: Controls how fast proper time advances
2. **Spatial metric (g_xx)**: Describes spatial distances
3. **Extrinsic curvature (K)**: Measures how the spatial slice is embedded in spacetime

All three quantities show characteristic wave-like propagation patterns as the gauge wave moves through the domain at approximately unit speed.

### Initial Conditions

The code starts with a Gaussian-shaped "bump" in the initial time slice:
```
h(x) = 0.1 * exp(-x²/6)
```

This creates a localized disturbance that propagates to the right, wraps around due to periodic boundaries, and re-enters from the left.

### Technical Implementation

- **Language**: Python (typically using NumPy)
- **Spatial domain**: x ∈ [-10, 40] with 400 grid points
- **Time evolution**: t = 0 to t = 10
- **Numerical methods**: 
  - 2nd-order centered finite differences for spatial derivatives
  - 4th-order Runge-Kutta for time stepping
  - CFL condition: Δt = 0.05 × Δx for stability

### Educational Value

This code is particularly valuable for:
- Students learning numerical relativity
- Understanding the difference between gauge and physical degrees of freedom
- Seeing how the ADM formalism works in its simplest setting
- Testing new numerical schemes before applying them to full 3+1D problems

### Output

The code produces plots showing the evolution of α(x,t), g_xx(x,t), and K(x,t) at different time snapshots, clearly demonstrating the gauge wave propagation, dipole structure, and periodic wrapping behavior.

---

**In Summary**: This is a pedagogical numerical relativity code that simulates pure coordinate effects (gauge waves) in the simplest possible setting of general relativity, serving as both a learning tool and a numerical methods testbed.
