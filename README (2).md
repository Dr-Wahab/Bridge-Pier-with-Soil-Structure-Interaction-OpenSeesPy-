# Bridge Pier with Soil–Structure Interaction (OpenSeesPy)

A reproducible **OpenSeesPy** model of a realistically sized reinforced-concrete
bridge pier founded on a large-diameter drilled shaft, with soil–structure
interaction (SSI) represented by a **Beam-on-Nonlinear-Winkler-Foundation
(BNWF)** spring system. The script builds the model, runs gravity, modal,
pushover, and nonlinear time-history analyses, and produces publication-style
figures of the model, deformed shape, and results.

![Model schematic](fig1_model.png)

## Features

- Single-column RC pier continuous with a drilled-shaft foundation
- Distributed-plasticity `forceBeamColumn` elements with **fiber sections**
  (confined core + unconfined cover `Concrete02`, `Steel02` reinforcement)
- SSI via depth-dependent `PySimple1` (lateral), `TzSimple1` (skin friction),
  and `QzSimple1` (tip bearing) springs, calibrated with the API RP2A sand method
- Four analysis stages: **gravity → modal → pushover → nonlinear time-history**
- Self-contained synthetic ground motion (Kanai–Tajimi filtered noise with a
  Jennings envelope), so no external record is required
- Consistent, publication-ready plots generated automatically

## Model summary

| Component | Properties |
|---|---|
| RC column | Ø1.5 m, H = 8 m, C35 concrete, 28 Ø32 longitudinal bars |
| Drilled shaft | Ø1.8 m, L = 15 m embedded, 32 Ø36 longitudinal bars |
| Reinforcing steel | Grade 420 (fy = 420 MPa), `Steel02` |
| Soil profile | Medium-dense sand, φ = 35°, γ′ = 10 kN/m³, water table at surface |
| Superstructure mass | 600 Mg (≈ 5.9 MN) lumped at the pier top |
| Units | kN, m, s, Mg, kPa |

## Requirements

```bash
pip install openseespy numpy matplotlib
```

Tested with Python 3.10+.

## Usage

```bash
python bridge_pier_ssi.py
```

All model, material, and soil parameters are exposed as named variables at the
top of the script for direct modification. Running the script regenerates the
five figures listed below.

## Outputs

Running the script produces the following five figures.

**Model schematic** — labelled components and BNWF soil springs:

![Model schematic](fig1_model.png)

**Deformed shape** at the end of the pushover analysis:

![Deformed shape](fig2_deformed.png)

**Pushover capacity curve** — base shear vs. top displacement, with a drift-ratio axis:

![Pushover capacity curve](fig3_pushover.png)

**Time-history response** — input motion, pier-top displacement, and base shear:

![Time-history response](fig4_timehistory.png)

**Hysteretic response** — base shear vs. top displacement:

![Hysteretic response](fig5_hysteresis.png)

## Representative results

| Quantity | Value |
|---|---|
| Fundamental period, T₁ | ≈ 1.66 s |
| Second / third periods | ≈ 0.77 s / 0.08 s |
| Pushover peak base shear | ≈ 944 kN (≈ 3% column drift) |
| Time-history peak top drift (PGA = 0.40 g) | ≈ 173 mm |
| Time-history peak base shear | ≈ 1117 kN |

## Assumptions and limitations

- 2D plane model (ndm = 2, ndf = 3); out-of-plane behaviour is not represented.
- Single homogeneous sand layer; soil-spring backbone parameters are
  **illustrative API-based calibrations**, not site-specific design values.
- `UniformExcitation` imposes a rigid free-field motion at the spring supports,
  so **kinematic** interaction is neglected while **inertial** SSI is captured.
- The synthetic accelerogram is for demonstration; substitute a recorded or
  spectrum-matched motion for design studies.

This repository is intended for research, teaching, and demonstration. It is
**not** a substitute for code-compliant design and has not been independently
verified for any specific structure or site.

## License

Released under the MIT License. Attribution is appreciated.
