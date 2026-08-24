# Project #15: Density-Functional Theory (DFT)

This project adds a guided introduction to **Kohn–Sham density-functional
theory**. It is intended to follow [the Hartree–Fock SCF project](../Project%2303):
the matrix machinery is deliberately familiar, while the exchange–correlation
energy and its numerical integration are new. The three examples below are
progressive. Complete them in order rather than attempting the molecular SCF
immediately.

> **Scope.** Here, DFT means *density-functional theory*, not the discrete
> Fourier transform. We use a closed-shell, spin-unpolarized density and the
> Dirac (Slater) local-density approximation (LDA) for exchange. Correlation,
> analytic gradients, unrestricted spin, and production-quality molecular
> quadratures are useful extensions, but are not required.

## Learning objectives

By the end of the project you should be able to:

1. explain the Hohenberg–Kohn and Kohn–Sham ideas without treating the
   Kohn–Sham orbitals as the exact many-electron wave function;
2. evaluate an LDA exchange energy and potential on a real-space grid;
3. transform a grid potential into an AO-basis matrix;
4. assemble and iterate the closed-shell Kohn–Sham equations; and
5. distinguish the total energy from the sum of Kohn–Sham orbital energies.

## Theory: from the density to working equations

The first Hohenberg–Kohn theorem says that, for a non-degenerate ground state,
the ground-state density determines the external potential up to a constant.
The second supplies a variational principle: the exact energy functional is
minimized by the exact ground-state density. These theorems establish that the
density is sufficient in principle; they do **not** provide the unknown
functional explicitly.

Kohn and Sham introduce a non-interacting reference system having the same
density as the interacting system. In atomic units, its energy is partitioned
as

\[
E[\rho]=T_s[\rho]+\int v_{\mathrm{ext}}(\mathbf r)\rho(\mathbf r)d\mathbf r
+J[\rho]+E_{\mathrm{xc}}[\rho]+E_{\mathrm{NN}},
\]

where \(T_s\) is the non-interacting kinetic energy, \(J\) is the classical
Coulomb energy, and everything not already represented exactly is assigned to
the exchange–correlation functional \(E_{\mathrm{xc}}\). Functional
differentiation gives the local potential

\[
v_{\mathrm{xc}}(\mathbf r)=\frac{\delta E_{\mathrm{xc}}}{\delta\rho(\mathbf r)}.
\]

Expanding the orbitals in an AO basis produces the generalized eigenproblem

\[
\mathbf F^{\mathrm{KS}}\mathbf C=\mathbf S\mathbf C\boldsymbol\epsilon,
\qquad
\mathbf F^{\mathrm{KS}}=\mathbf H+\mathbf J+\mathbf V^{\mathrm{xc}}.
\]

We use the **spin-summed** closed-shell density matrix
\(P_{\mu\nu}=2\sum_i^{\mathrm{occ}}C_{\mu i}C_{\nu i}\). With this convention,
\(\rho(\mathbf r)=\sum_{\mu\nu}P_{\mu\nu}\phi_\mu(\mathbf r)
\phi_\nu(\mathbf r)\), and
\(J_{\mu\nu}=\sum_{\lambda\sigma}P_{\lambda\sigma}
(\mu\nu|\lambda\sigma)\). Do not silently combine these formulas with the
half-density convention used by some Hartree–Fock texts.

For exchange-only LDA,

\[
E_x[\rho]=-C_x\int\rho(\mathbf r)^{4/3}d\mathbf r,
\quad C_x=\frac34\left(\frac3\pi\right)^{1/3},
\quad v_x(\mathbf r)=-\left(\frac3\pi\right)^{1/3}\rho(\mathbf r)^{1/3}.
\]

The factor \(4/3\) connecting the energy density and potential is an excellent
first debugging check: \(v_x\rho=(4/3)e_x\) at every point with nonzero
density.

## Example 1: evaluate LDA exchange on a tiny grid

### Why this example comes first

This example separates functional evaluation from SCF bookkeeping. A grid is
an integration rule, not merely a collection of plotting points:
\(\int f(\mathbf r)d\mathbf r\approx\sum_g w_gf(\mathbf r_g)\). The supplied
density values may be thought of as samples from any normalized model density.

### Input

Use these deliberately small values (atomic units):

| grid point \(g\) | \(w_g\) | \(\rho_g\) |
|---:|---:|---:|
| 0 | 0.50 | 0.80 |
| 1 | 1.00 | 0.20 |
| 2 | 0.75 | 0.05 |
| 3 | 1.25 | 0.00 |

### Instructions

1. Store weights and densities in `std::vector<double>` and verify that their
   lengths match.
2. For every point, calculate \(e_{x,g}=-C_x\rho_g^{4/3}\) and
   \(v_{x,g}=-(3/\pi)^{1/3}\rho_g^{1/3}\). Treat exactly zero density as zero;
   never evaluate a negative density to a fractional power.
3. Accumulate \(E_x=\sum_gw_ge_{x,g}\). Also accumulate the diagnostic
   \(I=\sum_gw_g\rho_gv_{x,g}\).
4. Print one row per point and at least ten digits for both totals.
5. Assert that \(I\approx(4/3)E_x\). A tolerance of \(10^{-12}\) is reasonable
   for this tiny data set.

Reference values are

```text
E_x                  = -0.370832523781
sum_g w_g rho_g v_xg = -0.494443365041
```

**Theory question:** Why is the last point harmless to the exchange energy,
even though low-density regions can still require many grid points for an
accurate molecular integral?

## Example 2: build an AO exchange–correlation matrix

### Theory

Because \(v_{\mathrm{xc}}\) is multiplicative in real space, its AO matrix is

\[
V^{\mathrm{xc}}_{\mu\nu}\approx\sum_gw_g
\phi_\mu(\mathbf r_g)v_{\mathrm{xc}}(\mathbf r_g)
\phi_\nu(\mathbf r_g).
\]

This equation connects a density functional evaluated point-by-point to the
matrix eigenproblem. The result must be symmetric for real basis functions.

### Input and instructions

Use two AO values at each point:

| \(g\) | \(w_g\) | \(\phi_0(\mathbf r_g)\) | \(\phi_1(\mathbf r_g)\) |
|---:|---:|---:|---:|
| 0 | 0.6 | 0.8 | 0.1 |
| 1 | 0.9 | 0.5 | 0.4 |
| 2 | 0.7 | 0.2 | 0.7 |

and the trial spin-summed density matrix

\[
\mathbf P=\begin{pmatrix}1.20&0.30\\0.30&0.80\end{pmatrix}.
\]

1. At each point form
   \(\rho_g=\sum_{\mu\nu}P_{\mu\nu}\phi_{\mu g}\phi_{\nu g}\).
   Sum over **all** matrix elements; if you exploit symmetry, remember the
   factor of two on off-diagonal terms.
2. Reject a density below \(-10^{-12}\) as an error. Clamp values in
   \([-10^{-12},0)\) to zero to tolerate roundoff only.
3. Evaluate \(e_{x,g}\) and \(v_{x,g}\), add \(w_ge_{x,g}\) to \(E_x\), and
   add \(w_g\phi_{\mu g}v_{x,g}\phi_{\nu g}\) to every matrix element.
4. Check \(\max_{\mu\nu}|V^x_{\mu\nu}-V^x_{\nu\mu}|<10^{-12}\).

Reference results are

```text
rho = [0.8240000000, 0.5480000000, 0.5240000000]
Ex  = -0.8588097459
Vx  = [[-0.5580559013, -0.2671683053],
       [-0.2671683053, -0.3938894955]]
```

**Theory question:** The trace of \(\mathbf P\) is not generally the electron
number in a non-orthogonal AO basis. Which contraction with the overlap matrix
should equal \(N\)?

## Example 3: a closed-shell molecular LDA SCF

This example turns the previous kernels into a working Kohn–Sham program. Use
the H\(_2\)O/STO-3G nuclear repulsion, overlap, kinetic, nuclear-attraction,
and electron-repulsion integral files from [Project #3](../Project%2303/input/h2o/STO-3G).
You must additionally supply or generate an atom-centered quadrature and AO
values on that grid; a production grid normally combines a radial rule,
angular points, and molecular partition weights.

### SCF algorithm

1. Read the integrals, form \(\mathbf H=\mathbf T+\mathbf V\), verify matrix
   symmetry, and build \(\mathbf X=\mathbf S^{-1/2}\).
2. Obtain an initial \(\mathbf P\) from the core Hamiltonian or a converged
   Project #3 density converted consistently to the spin-summed convention.
3. Evaluate \(\rho_g\), \(E_x\), and \(\mathbf V^x\) using Example 2.
4. Build \(\mathbf J\) from the two-electron integrals and form
   \(\mathbf F^{KS}=\mathbf H+\mathbf J+\mathbf V^x\). There is no exact
   Hartree–Fock exchange matrix in this exchange-only LDA example.
5. Transform \(\mathbf F'=\mathbf X^T\mathbf F^{KS}\mathbf X\), diagonalize
   it, back-transform the coefficients, occupy the lowest \(N/2\) orbitals,
   and form a new spin-summed density.
6. Mix the density, for example
   \(\mathbf P\leftarrow0.75\mathbf P_{new}+0.25\mathbf P_{old}\). DIIS from
   [Project #8](../Project%2308) is a stronger extension once simple damping
   works.
7. Evaluate

   \[
   E=\mathrm{Tr}[\mathbf P\mathbf H]+\tfrac12\mathrm{Tr}[\mathbf P\mathbf J]
   +E_x+E_{NN}.
   \]

   Do not use the Hartree–Fock energy formula: \(E_x\) is already an integrated
   functional and is not \(\tfrac12\mathrm{Tr}[\mathbf P\mathbf V^x]\).
8. Iterate until both \(|E_n-E_{n-1}|<10^{-10}\ E_h\) and the RMS density
   change is below \(10^{-8}\). Report the electron-count check
   \(\mathrm{Tr}[\mathbf P\mathbf S]\), grid-integrated electron count, energy
   components, and SCF iteration count.

### Pseudocode

```cpp
P = initial_density(H, S, nelectron);
for (int iter = 0; iter < max_iter; ++iter) {
    auto [rho, Ex, Vx, Ngrid] = lda_exchange(P, ao_values, weights);
    Matrix J = coulomb_matrix(P, eri);
    Matrix F = H + J + Vx;
    Matrix Pnew = occupied_density(diagonalize(X.transpose() * F * X), X);
    double E = trace(P * H) + 0.5 * trace(P * J) + Ex + Enuc;
    report(iter, E, rms(Pnew - P), trace(P * S), Ngrid);
    if (energy_and_density_converged()) break;
    P = 0.75 * Pnew + 0.25 * P;
}
```

The pseudocode emphasizes data flow, not a required API. Keep grid evaluation
in a separate function so it can be unit-tested using Examples 1 and 2.

### Interpretation questions

1. Why does refining the quadrature change the energy even when the AO basis
   and integral files are unchanged?
2. Why is the highest occupied Kohn–Sham eigenvalue special in exact DFT, while
   arbitrary orbital-energy differences are not generally exact excitation
   energies?
3. Which pieces of \(E[\rho]\) are approximated in this project, and which are
   evaluated exactly within the chosen AO basis and quadrature?

## Validation and debugging checklist

- [ ] `P`, `J`, `Vx`, and the Kohn–Sham matrix are symmetric.
- [ ] `trace(P * S)` agrees with the number of electrons.
- [ ] The grid integral \(\sum_gw_g\rho_g\) approaches the electron number as
      the grid is refined.
- [ ] All occupied orbitals satisfy \(\mathbf C_{occ}^T\mathbf S
      \mathbf C_{occ}=\mathbf I\).
- [ ] Example 1 passes before Example 2, and Example 2 passes before the SCF.
- [ ] Tightening the grid changes the converged energy less than the requested
      accuracy.
- [ ] Convergence requires both an energy criterion and a density criterion.

## Further examples and extensions

After the three required examples work, try one change at a time:

* **Grid convergence:** run coarse, medium, and fine quadratures and tabulate
  the electron-count error and total-energy change. This distinguishes SCF
  convergence from integration convergence.
* **Local correlation:** add an LDA correlation parameterization behind the
  same `evaluate(rho)` interface, returning both energy per volume and its
  functional derivative. Validate derivatives with finite differences.
* **Spin DFT:** maintain \(\rho_\alpha\) and \(\rho_\beta\) separately and test
  the closed-shell limit before attempting an open-shell atom.
* **A GGA:** evaluate density gradients and include the integration-by-parts
  terms needed for the potential. Merely replacing the LDA energy expression
  without changing the potential is not a self-consistent GGA implementation.
* **Hybrid DFT:** mix a specified fraction of exact exchange into the
  Kohn–Sham matrix and correct the energy expression consistently. Compare its
  additional cost with LDA.

## Common conceptual mistakes

* DFT is exact only with the exact functional; LDA is an approximation.
* Kohn–Sham orbitals reproduce the density of the auxiliary non-interacting
  system. They are not a Slater-determinant representation of the exact
  interacting wave function.
* A converged SCF on an inadequate grid is precisely converged to an inaccurate
  discretized problem.
* `Ex`, `trace(P * Vx)`, and exact Hartree–Fock exchange are three different
  quantities.
* Density-matrix conventions determine factors of two. Document the convention
  at every public function boundary rather than repairing factors empirically.

## Suggested reading

* P. Hohenberg and W. Kohn, *Physical Review* **136**, B864 (1964).
* W. Kohn and L. J. Sham, *Physical Review* **140**, A1133 (1965).
* R. G. Parr and W. Yang, *Density-Functional Theory of Atoms and Molecules*,
  Oxford University Press (1989).

## Supplied molecular test cases

The repository includes complete H\(_2\)O and CH\(_4\) test cases following the
same molecule/basis directory convention as Project #3:

| Molecule | Input | Reference output |
|---|---|---|
| H\(_2\)O/STO-3G | [`input/h2o/STO-3G`](./input/h2o/STO-3G) | [`output/h2o/STO-3G/output.txt`](./output/h2o/STO-3G/output.txt) |
| CH\(_4\)/STO-3G | [`input/ch4/STO-3G`](./input/ch4/STO-3G) | [`output/ch4/STO-3G/output.txt`](./output/ch4/STO-3G/output.txt) |

Each input directory contains `geom.dat`, `enuc.dat`, `s.dat`, `t.dat`,
`v.dat`, and `eri.dat` in exactly the formats described in Project #3. The
additional `grid.dat` records the Cartesian midpoint quadrature used to make
the reference output. Grid limits are half-open: generate points
`min + (i + 0.5) * spacing` while the result is less than `max`; the weight of
every point is `spacing^3`.

These intentionally simple Cartesian grids make independent implementations
reproducible, but they are **teaching grids**, not recommended molecular DFT
quadratures. In particular, tight core functions converge slowly on a uniform
Cartesian mesh. Reproduce the supplied files first, then replace this grid with
an atom-centered radial/angular quadrature and demonstrate convergence of both
the total energy and integrated electron count.

### Running the two examples

1. Copy or select one molecule's `STO-3G` directory; never mix matrices from
   different directories.
2. Determine `nao` from the largest one-electron-integral index and determine
   the electron count from `geom.dat`. Both examples are neutral, closed-shell
   ten-electron systems, so five spatial orbitals are occupied.
3. Read the six Project #3 files and `grid.dat`, construct AO values at each
   midpoint, and follow the SCF algorithm above using the stated spin-summed
   density convention.
4. Print an iteration table containing total energy, RMS density change, and
   grid-integrated electron count. Print final orbital energies and
   `Tr[P S]`. The supplied output uses these fields so results can be compared
   line-by-line.
5. First compare `Tr[P S]`, which should be ten to numerical precision. Then
   compare the grid electron count and energy. Differences in either usually
   indicate a basis normalization, AO ordering, grid-origin, or endpoint error.
6. After matching the reference, halve the spacing. Treat the result as a grid
   convergence study rather than expecting the coarse-grid energy to be a
   chemically accurate literature value.

CH\(_4\) also provides a useful symmetry check: its three occupied valence
orbital energies and three virtual valence orbital energies occur in degenerate
sets. A grid or AO ordering that breaks tetrahedral symmetry will spuriously
split these values.
