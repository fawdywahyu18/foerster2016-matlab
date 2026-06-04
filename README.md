# MATLAB Replication: Foerster et al. (2016)

Unofficial MATLAB translation of the replication package for:

> Foerster, A., Rubio-Ramirez, J.F., Waggoner, D.F., and Zha, T. (2016).
> **Perturbation Methods for Markov-Switching dynamic stochastic general equilibrium models.**
> *Quantitative Economics*, 7, 637–669.
> https://doi.org/10.3982/QE596


The original replication package by the authors is written in Mathematica. This repository provides a MATLAB translation that replicates Table 1 (Section 5) and Table 3 (Section 6) of the paper.

---

## Replication Status

**All partition perturbation results match the paper exactly.**

| Result | Method | Order | Computed | Paper | Status |
|--------|--------|-------|----------|-------|--------|
| Table 1 EE error | Partition | 1st | -Inf | -Inf | Pass |
| Table 1 EE error | Naive | 1st | -0.5564 | -0.5564 | Pass |
| Table 1 EE error | Naive | 2nd | -0.3171 | -1.3691 | Fail |
| Table 3 EE error | Partition | 1st | -3.01 | -3.01 | Pass |
| Table 3 EE error | Partition | 2nd | -3.59 | -3.59 | Pass |
| Table 3 EE error | Naive | 1st | -2.48 | -2.48 | Pass |
| Table 3 EE error | Naive | 2nd | -3.07 | -3.07 | Pass |

The only discrepancy is in the **naive perturbation 2nd order** result for the simple inflation model (Table 1). All first-order coefficients (G, H, Psi, Phi) for both partition and naive methods pass numerical verification. See `verify_results.m` for the full check.

---

## Repository Structure

```
foerster2016-matlab/
│
├── README.md
├── setup.m                       # Add all folders to MATLAB path
│
├── section5_simple_model/
│   └── main_simple_model.m       # Replicates Table 1
│
├── section6_rbc_model/
│   ├── main_rbc_model.m          # Replicates Table 3
│   ├── rbc_steady_state.m        # Steady state: c_ss, k_ss, z_ss
│   ├── rbc_eval_f.m              # 3 equilibrium equations
│   ├── rbc_jacobians.m           # Numerical Jacobians (central differences)
│   ├── solve_rbc_first_order.m   # First-order solver (slope, impact, chi)
│   ├── eval_rbc_F.m              # Evaluate F_i(z) via Gauss-Hermite quadrature
│   └── solve_rbc_second_order.m  # Second-order solver (numerical Hessian)
│
├── shared/
│   ├── ergodic_probs.m           # Ergodic probabilities from transition matrix
│   ├── check_mss.m               # Mean Square Stability criterion
│   └── compute_ee_rbc.m          # Euler equation errors via simulation
│
└── verify_results.m              # Automated check against paper values
```

---

## How to Run

First, run `setup.m` to add all folders to the MATLAB path:

```matlab
% Run once at the start of each MATLAB session
setup
```

Then run:

```matlab
% Replicate Table 1 (Section 5 - Simple Inflation Model)
main_simple_model

% Replicate Table 3 (Section 6 - RBC Model)
main_rbc_model

% Verify all results against the paper
verify_results
```

**Requirements:** MATLAB with Optimization Toolbox (`fsolve`). No other toolboxes needed.

---

## Model Overview (Section 6, RBC)

Three equilibrium equations: Euler equation, resource constraint, technology process.

| Variable | Description | Dimension |
|----------|-------------|-----------|
| y_t = c_t | Consumption (control) | ny = 1 |
| x_t = [k_t, z_t] | Capital, technology (states) | nx = 2 |
| eps_t | Technology shock | ne = 1 |
| s_t | Markov regime | ns = 2 |

**Parameters (Table 2):**

| Parameter | Value |
|-----------|-------|
| alpha | 0.33 |
| beta | 0.9976 |
| nu | -1 |
| delta | 0.025 |

**Regime parameters:**

| | Regime 1 | Regime 2 |
|-|----------|----------|
| mu | 0.0274 | -0.0337 |
| rho | 0.10 | 0.00 |
| sigma | 0.0072 | 0.0216 |

Transition matrix: P = [0.75, 0.25; 0.50, 0.50]

---

## Notes

- This is an **unofficial translation** and is not affiliated with the original authors.
- The original Mathematica replication package is available from the authors via the journal website.
- The partition perturbation method does **not** require manual log-linearization — equilibrium conditions are entered in nonlinear form and Jacobians are computed numerically.

---

## License

MIT License. If you use this code, please cite the original paper by Foerster et al. (2016).
