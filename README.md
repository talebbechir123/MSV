# Within-Host Evolutionary Dynamics of Plant Viruses — ODE Model in R

Mathematical model of competition and mutation between four variants of **Potato virus Y (PVY)** inside a single pepper plant. The variants carry one or two amino-acid substitutions in the VPg region, which affects pathogenicity. The model reproduces how selection and mutation shape variant frequencies over time, implementing a simplified version of the approach of Fabre *et al.*

Project for *Introduction à la Modélisation en Sciences du Vivant et en Santé* (IMSV), M2 MSCHP, Université Paris-Saclay. Pair work by **Ahmed Taleb Bechir** and **Abdelghani Ameziane**, supervised by Lulla Opatowski. Report: [`Projett_MSV.pdf`](Projett_MSV.pdf) · slides: [`MSV_slides.pdf`](MSV_slides.pdf).

## Model

For variants `V₁…V₄` (NN, DN, NH, DH), the population of each variant follows a Lotka-Volterra competition model with mutation:

```
dVᵢ/dt = rᵢ Vᵢ (1 − (Vᵢ + Σⱼ≠ᵢ bᵢⱼ Vⱼ) / K)  +  Σⱼ mᵢⱼ (Vⱼ − Vᵢ)
```

- `rᵢ` — intrinsic growth rate (fitness) of variant *i*
- `K` — within-host carrying capacity
- `bᵢⱼ = rⱼ / rᵢ` — competition coefficients
- `mᵢⱼ` — mutation rate: `μ` for variants one substitution apart, `μ²` for two

The system is integrated over 35 days from the inoculum frequencies (0.32 / 0.22 / 0.22 / 0.24) with `deSolve`, and the resulting variant frequencies are plotted over time.

## Run

```r
install.packages("deSolve")
source("projet.R")
```

Parameters (`r1…r4`, `K`, `mut`) and initial conditions are set at the top of the script section *Set up parameters*.

## Repository contents

| File | Description |
|---|---|
| `projet.R` | ODE model, numerical integration and frequency plot |
| `Projett_MSV.pdf` | Project report |
| `MSV_slides.pdf` | Presentation slides |
| `J1_…`, `J2_…`, `J3_…`, `J4_…` | Course lecture notes (epidemiology, ODE methods, within-host microbiology, cellular automata, parameter estimation) |
| `Simpson et al - 2010.pdf` | Reference article from the course material |

## Authors

**Ahmed Taleb Bechir** — [atbechir.com](https://atbechir.com) · [LinkedIn](https://www.linkedin.com/in/ahmed-taleb-bechir-97a7301a4)
**Abdelghani Ameziane**
