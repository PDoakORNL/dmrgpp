# cincuenta Design Document

## Overview

cincuenta is a C++ implementation of the Dynamical Mean-Field Theory (DMFT)
self-consistency algorithm for strongly correlated electron models.  It is
distributed as part of the DMRG++ repository and compiles alongside it, using
DMRG++ as its primary impurity solver.  A secondary, exact-diagonalization
solver based on LanczosPlusPlus is also provided for validation and small
system tests.

The name "cincuenta" is Spanish for fifty; cincuenta is developed around
version 50 (0.50) of the DMRG++ codebase.

UML artifacts for this design live alongside this file:

* `classes.mmd` — class diagram
* `sequence_dmft.mmd` — DMFT loop sequence diagram
* `component_arch.mmd` — component/dependency diagram

---

## Physical Background

### Dynamical Mean-Field Theory

DMFT maps a lattice model of strongly correlated electrons onto a single
impurity embedded in a self-consistent bath.  The key approximation is that the
self-energy Σ(iωₙ) is local (momentum-independent).  Given that assumption the
problem reduces to:

1. Choose a self-energy Σ(iωₙ).
2. Compute the local lattice Green's function G_latt from Σ.
3. Determine the bath that the impurity must sit in so that solving the
   impurity problem reproduces the same local Green's function.
4. Solve the impurity problem to get G_imp and the new self-energy.
5. Repeat until Σ converges.

The "bath" in step 3 is parameterised as a finite set of non-interacting
orbitals (bath sites) coupled to the interacting impurity site — the Anderson
impurity model.

### Anderson Impurity Model

cincuenta uses a star geometry: one correlated impurity site connected to
`nBath` non-interacting bath sites with no hopping between bath sites.  The
Hamiltonian is

```
H = U n↑ n↓  +  Σ_α [ V_α (c†_imp c_α + h.c.)  +  ε_α c†_α c_α ]
```

The hybridisation function

```
Δ(iωₙ) = Σ_α  V_α² / (iωₙ + μ − ε_α)
```

encodes the effect of the bath on the impurity.  The 2·nBath real numbers
`{V_α, ε_α}` are the "bath parameters" that must be fitted at each DMFT
iteration.  They are stored as a flat vector with hoppings first:

```
bathParams = [ V₀, V₁, …, V_{nBath−1},  ε₀, ε₁, …, ε_{nBath−1} ]
```

---

## Component Architecture

Four external components collaborate (see `component_arch.mmd`):

| Component | Role |
|-----------|------|
| **cincuenta** | DMFT driver: self-consistency loop, bath fitting, lattice G |
| **DMRG++** | Impurity solver (ground state + spectral functions via DMRG) |
| **LanczosPlusPlus** | Impurity solver (exact diagonalisation + continued fractions) |
| **PsimagLite** | Shared utilities: I/O, linear algebra, frequency grids, minimiser |

cincuenta owns the self-consistency loop and delegates the heavy numerics:
impurity solving to DMRG++ or LanczosPlusPlus, and bath fitting to
`PsimagLite::Minimizer`.

---

## The DMFT Self-Consistency Loop

`DmftSolver::selfConsistencyLoop()` is the top-level driver.  Each iteration
performs four steps in sequence.

### Step 1 — Update the lattice Green's function

`LatticeGf::update()` recomputes G_latt(iωₙ) from the current self-energy
Σ(iωₙ) for each Matsubara frequency.  Two modes are supported, selected by
the `LatticeGf=` input keyword:

**Momentum-space** (`momentum,1D,Nk`):
```
G_latt(iωₙ) = (1/Nk) Σ_k  1 / (iωₙ − εk + μ − Σ(iωₙ))
```
`Dispersion` provides the band energies εk.  Currently only the 1D
tight-binding dispersion εk = −2 cos(k) is implemented.

**Energy-space** (`energy,semicircular,W`):
```
G_latt(iωₙ) = ∫ dε  ρ(ε) / (iωₙ − ε + μ − Σ(iωₙ))
```
`DensityOfStates` provides ρ(ε); currently only the semicircular DOS is
implemented.  `PsimagLite::Integrator` (GSL adaptive quadrature) evaluates
the integral numerically.

After G_latt is known, the γ-function is formed:
```
Γ(iωₙ) = iωₙ − 1/G_latt(iωₙ) − Σ(iωₙ)
```
Γ is the target hybridisation function that the Anderson bath must reproduce.

### Step 2 — Fit the bath parameters

`Fit::fit()` minimises the least-squares objective
```
F(V, ε) = (1/N) Σₙ | Δ(V, ε, iωₙ) − Γ(iωₙ) |²
```
where Δ is the Anderson hybridisation function (see above) and the sum runs
over all Matsubara frequencies.

`FitFunction` implements both `operator()` (the objective) and `df` (its
gradient).  `df` uses analytic derivatives via `AndersonFunction::andersonPrime`:
```
∂F/∂V_α  =  2 V_α / (iωₙ + μ − ε_α)
∂F/∂ε_α  =  V_α² / (iωₙ + μ − ε_α)²
```

`PsimagLite::Minimizer` wraps GSL's conjugate-gradient or Nelder–Mead simplex
algorithm, selected by `FitMethod=` in the input.

**Particle-hole symmetry option**: when `FitOptions=ParticleHoleSymmetric` is
set, the bath energies are constrained to come in ±ε pairs (plus an optional
zero-energy site for odd `nBath`).  This halves the number of free energy
unknowns, making the fit better conditioned for half-filled, particle-hole
symmetric models.

### Step 3 — Solve the impurity problem

`ImpuritySolverBase::solve(bathParams, freqEnum, iter)` is the interface.  The
concrete solver is chosen at construction time from the `ImpuritySolver=` input
keyword and held by a raw `ImpuritySolverBase*` pointer in `DmftSolver`.

Both solvers build a `ModelParams` object from `bathParams`.  `ModelParams`
maps the flat bath-parameter vector onto the star-geometry Hamiltonian that
DMRG++ and LanczosPlusPlus expect: hoppings `V_α` connect the impurity
(centre) to bath sites, and on-site energies `ε_α` go on the bath sites.

#### ImpuritySolverDmrg

The DMRG solver proceeds in two stages:

1. **Ground state**: Constructs a DMRG++ input string embedding the current
   bath parameters and calls `Dmrg::DmrgRunner::doOneRun()`.

2. **Spectral function**: Calls `ManyOmegas` to loop over the requested
   Matsubara or real frequencies, running a correction-vector DMRG calculation
   at each.  `ProcOmegas` parses the resulting output files and assembles
   `gimp_`.

Because different frequencies are independent, the MPI parallelism in DMRG++
can distribute them across nodes.

#### ImpuritySolverExactDiag

The exact-diagonalisation solver is intended for small `nBath` (typically ≤ 5)
and serves primarily as a reference and validation tool.

1. **Hamiltonian build**: `LanczosPlusPlus::ModelSelector` constructs the full
   many-body Hamiltonian matrix in the occupation-number basis using the star
   geometry and `HubbardOneBand` model.

2. **Lanczos diagonalisation**: `LanczosPlusPlus::Engine` finds the ground
   state via the Lanczos algorithm.

3. **Green's function**: The resolvent
   G_imp(ω) = ⟨GS| c  (ω − H + E_GS)⁻¹ c† |GS⟩  + (ω → −ω) term
   is evaluated as a continued fraction.
   `PsimagLite::ContinuedFractionCollection` accumulates the Lanczos
   coefficients and evaluates the fraction at the requested frequencies.

### Step 4 — Update the self-energy

`DmftSolver::computeNewSelfEnergy()` uses the Dyson equation to extract the
new self-energy from the solved impurity Green's function:
```
Σ_new(iωₙ) = iωₙ + μ − Δ(iωₙ) − 1/G_imp(iωₙ)
```
The convergence error is the mean-squared change:
```
error = (1/N) Σₙ |Σ_new(iωₙ) − Σ_old(iωₙ)|²
```
If `error < DmftTolerance` the loop exits; otherwise Σ is overwritten with
Σ_new and the next iteration begins.

After the Matsubara-frequency loop converges (or exhausts `DmftNumberOfIterations`)
a final impurity solve is run at real frequencies to produce the spectral
function A(ω) = −Im[G_imp(ω + iη)] / π.

---

## Supporting Data Structures

### FunctionOfFrequency

A frequency-indexed vector of complex values, parameterised by `ficticiousBeta`
and `nMatsubaras`.  The underlying `PsimagLite::Matsubaras` object owns the
grid ωₙ = (2n + 1)π/β.  The container stores 2·nMatsubaras elements so that
both positive and negative Matsubara frequencies are available.

`sigma_`, `latticeG_`, and `gammaG_` in `DmftSolver` are all
`FunctionOfFrequency` objects sharing the same frequency grid.

### ModelParams

`ModelParams` is a thin adapter between the flat `bathParams` vector and the
geometry/Hamiltonian structures that the impurity solvers need.  It constructs
the star geometry (`PsimagLite::Star`) and sets up hopping and on-site
potential arrays in the format expected by DMRG++ input files.

### ParamsDmftSolver

Read-once configuration struct populated directly from the `PsimagLite::InputNg`
reader at startup.  All solver components receive a const reference to it.

---

## Input File Reference

| Keyword | Type | Description |
|---------|------|-------------|
| `FicticiousBeta` | real | Inverse fictitious temperature β; sets Matsubara grid spacing |
| `ChemicalPotential` | real | μ, chemical potential |
| `Matsubaras` | int | Number of positive Matsubara frequencies |
| `LatticeGf` | string | Lattice G mode, e.g. `momentum,1D,64` or `energy,semicircular,2.0` |
| `NumberOfBathPoints` | int | nBath, number of bath sites |
| `DmftNumberOfIterations` | int | Maximum self-consistency iterations |
| `DmftTolerance` | real | Convergence threshold for Σ error |
| `ImpuritySolver` | string | `dmrg` or `exactdiag` |
| `FitMethod` | string | `conjugate_gradient` or `simplex` |
| `FitOptions` | string | Optional: `ParticleHoleSymmetric` |
| `MinParamsMaxIter` | int | Maximum iterations for bath parameter minimiser |

---

## Output Files

After convergence `DmftSolver::print()` writes the following sections to
standard output (redirected to a file by convention):

| Section | Contents |
|---------|----------|
| `Sigma` | Self-energy Σ(iωₙ): frequency, Re, Im |
| `BathParams` | Optimised V_α and ε_α vectors |
| `LatticeG` | G_latt(iωₙ) |
| `Gamma` | Γ(iωₙ), the bath hybridisation target |
| `SiteExcludedG` | Local Green's function with impurity site removed |
| `AndersonFunction` | Δ(iωₙ) evaluated at converged bath parameters |

In addition, `logDebug()` writes `gimp_dmrg.txt` / `gimp_exactdiag.txt`
(Matsubara) and `gimp_dmrg_real.txt` (real-frequency) after every impurity
solve so intermediate results are available even if the run is interrupted.

The utility executable `freeGimp` reads a set of bath parameters from an input
file and outputs the corresponding free (U = 0) impurity Green's function.
This is useful for verifying bath fitting in isolation.
