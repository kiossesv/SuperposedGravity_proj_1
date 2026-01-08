# SuperposedGravity_proj_1

# Numerical Time Evolution of Quantum States in Superposed Gravitational Fields


## Overview
This project implements a numerical solver for the time-dependent Schrödinger equation in one spatial dimension using a split-operator Fourier method.
The solver is applied to quantum systems evolving under gravitational potentials, with particular emphasis on branch-dependent gravitational fields versus an effective symmetric gravitational field.
The project is motivated by foundational questions in quantum mechanics and gravity, and focuses on clean numerical experimentation, stability, and physical interpretability.

## Scientific Motivation
In standard non-relativistic quantum mechanics, time evolution is defined with respect to a single global frame.
However, when gravitational fields are associated with superposed configurations, the consistency of quantum evolution becomes nontrivial.
This project explores, at a numerical level:
- how quantum states evolve under branch-dependent gravitational potentials
- how this evolution compares to an effective symmetric gravitational field
The goal is not to test a specific theory, but to provide a robust computational framework for comparing different gravitational scenarios.

## Physical Model
### Quantum System
- Single non-relativistic quantum particle
- One spatial dimension
- External gravitational potential
### Governing Equation
The time evolution is governed by the time-dependent Schrödinger equation:
$$ i \hbar \frac{\partial \psi(x,t)}{\partial t} = \left[ -\frac{\hbar^2}{2m} \frac{\partial^2}{\partial x^2} + V_g(x)\right] \psi(x,t) $$
where:
- $\hat{T} = -\frac{\hbar^2}{2m} \frac{d^2}{d x^2}$ is the kinetic operator
- $V_g(x)$ is the gravitational potential

### Gravitational Scenarios
#### Case A — Branch-dependent gravitational field
Two independent gravitational potentials are considered:
$$ V_g^{(1)}(x) = m g_1 x, \qquad  V_g^{(2)}(x) = m g_2 x $$
Each branch evolves under its own Hamiltonian.
#### Case B — Effective symmetric gravitational field
A single effective gravitational potential:
$$
V_g^{(eff)}(x) = m g_{eff} x,
$$
This provides a reference evolution with a common time parameter.

## Numerical Method


## Computational Implementation


## Results


## Outlook
