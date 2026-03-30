# Numerical Time Evolution of Quantum States in Superposed Gravitational Fields
### Project 1 — Branch-Resolved Quantum Dynamics
This repository provides a numerical framework for studying quantum dynamics in gravitational fields, with particular emphasis on branch-dependent gravitational potentials arising from superposed configurations.

The code implements a spectral split-operator solver for the time-dependent Schrödinger equation and provides diagnostics for verifying the numerical consistency of simulations involving multiple gravitational branches.


## Research Program

This project is part of a three-repository research program investigating quantum superpositions of gravitational fields:
- [Project 1: quantum-gravity-superposition-dynamics](https://github.com/kiossesv/quantum-gravity-superposition-dynamics)

Numerically stable solver for branch-resolved quantum dynamics in superposed gravitational fields.

- [Project 2: quantum-gravity-superposition-stability](https://github.com/kiossesv/quantum-gravity-superposition-stability) (this repository)

Statistical stability analysis of effective gravitational models.

- [Project 3: quantum-gravity-superposition-sciml](https://github.com/kiossesv/quantum-gravity-superposition-sciml) (planned)

Physics-informed machine learning for predicting stability regimes of effective gravity.

## Overview
The repository implements a solver for the time-dependent Schrödinger equation in one spatial dimension.
![Schrödinger equation](docs/schrodinger_equation.svg)

The solver is applied to systems evolving under gravitational potentials, allowing comparison between:
- branch-dependent gravitational fields
- effective symmetric gravitational descriptions

The purpose of the project is not to test a specific theory but to provide a **controlled computational environment** for exploring quantum evolution in such scenarios.

## Scientific Motivation
When a mass is placed in a quantum superposition of positions, the gravitational field associated with it may also be considered branch-dependent.

This raises a conceptual question:
`How should a quantum system evolve in the presence of such branch-dependent gravitational fields?`

This repository provides a numerical framework for studying this question through:
- branch-resolved Schrödinger evolution
- comparison with effective gravitational descriptions
- verification of numerical consistency and stability.

## Numerical Method
### Spectral Split-Operator Propagation
The solver uses a Fourier spectral discretization combined with a second-order Strang split-operator scheme.

Time evolution over a small time step $\Delta t$ is approximated as:
![Time evolution](docs/time_evolution.svg)

where
- $T$ is the kinetic operator
- $V$ is the gravitational potential.

This method provides
- unitary time evolution
- high spectral accuracy
- efficient FFT-based implementation.

## Physical Model
### Quantum System
- single non-relativistic particle
- one spatial dimension
- external gravitational potential

### Gravitational Potentials
Two scenarios are considered:

**Branch-Dependent Fields**

Each branch of a gravitational configuration produces its own potential, $V_i(x) = mg_i x$.
The quantum state evolves independently within each branch.

**Effective Symmetric Field**

A single effective gravitational field represents the averaged configuration.
This allows comparison between branch-resolved dynamics and effective descriptions.

## Numerical Verification
The repository includes diagnostics ensuring the numerical reliability of simulations.

### Norm Conservation
The split-operator scheme preserves the wavefunction norm to machine precision, confirming discrete unitarity.

### Ehrenfest Trajectory Test
For linear gravitational potentials, the expectation value of position follows the analytical uniformly accelerated trajectory.
Simulations reproduce this behaviour with high accuracy.

### Spectral Time-of-Validity Bound
Uniform gravitational acceleration produces a drift of the wavepacket in momentum space.

Because the solver uses a periodic Fourier grid, this introduces a maximum reliable simulation time before spectral aliasing occurs.

The repository includes numerical tests confirming this bound.

### Branch Independence
Multi-branch simulations verify that independent gravitational sectors evolve without spurious numerical coupling.

## Example Output
Typical simulations compute
- wavefunction evolution
- probability density dynamics
- expectation values
- comparison between gravitational branches.

These outputs allow direct inspection of branch-resolved quantum dynamics.


## Code Structure
The implementation is modular and separates numerical infrastructure from physical modeling.
Main components include
```
solvers.py              Split-operator propagation engine
potentials.py          Gravitational potential definitions
initial_states.py      Gaussian wavepacket initialization
observables.py       Expectation value calculations
```
This structure allows new potentials and models to be added without modifying the core solver.

## Reproducibility
All simulations are fully deterministic.

Each run is defined by explicit parameters:
- domain size $L$
- grid resolution $N$ 
- timestep $\Delta t$
- total simulation time $T$
- gravitational accelerations $g_i$
- initial wavepacket parameters

No stochastic elements are present, ensuring exact reproducibility for a given parameter set.

## Related Publication
The numerical framework implemented in this repository is described in the paper

**Branch-Resolved Spectral Dynamics for Quantum Gravitational Fields:
A Split-Operator Approach**

The paper presents
- the theoretical formulation
- numerical verification tests
- performance analysis
- applications to gravitational branch dynamics.

## License
MIT License.
