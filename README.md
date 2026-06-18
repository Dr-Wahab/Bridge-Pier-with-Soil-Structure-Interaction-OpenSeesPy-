# Bridge Pier with Soil-Structure Interaction (OpenSeesPy)

Authors: Dr. Wahab (and contributors)  
Repository: Bridge-Pier-with-Soil-Structure-Interaction-OpenSeesPy  
Date: 2026-06-18

## Abstract
This repository contains models and scripts to study the dynamic and static response of bridge piers interacting with surrounding soil using OpenSeesPy. The work focuses on soil-structure interaction (SSI) effects on the seismic performance of piers, including nonlinear material behavior, radiation damping approximations, and parametric studies of soil stiffness and pier geometry. All models are implemented with OpenSeesPy to promote reproducible research and easy extension.

Keywords: soil-structure interaction, bridge pier, OpenSeesPy, nonlinear analysis, seismic response, reproducibility

## Repository overview
This repository provides:
- OpenSeesPy scripts to build and analyze pier-soil systems.
- Example input files for parametric and time-history analyses.
- Post-processing utilities to extract displacements, accelerations, bending moments, and foundation forces.
- Documentation on model assumptions, validation checks, and expected outputs.

Intended audience: researchers and engineers performing computational SSI studies with OpenSeesPy and those wishing to reproduce or extend the provided simulations.

## Background and motivation
Soil-structure interaction influences the seismic demands on bridge piers, altering displacement and force distributions compared to fixed-base assumptions. Accurate SSI modelling is essential when foundation flexibility, embedment, or wave propagation effects are non-negligible. This repository demonstrates practical OpenSeesPy implementations of commonly used SSI modeling approaches, emphasizing reproducibility and documentation of model choices.

## Methodology and model description
- Structural model:
  - Bridge pier modeled as a 2D/3D frame element assembly (beam-column elements) with geometric nonlinearity where required.
  - Material models: use of OpenSees material modules (e.g., `Concrete01`, `Steel02`) — replace with the material models used in your scripts.
  - Boundary conditions: pier connected to foundation nodes representing pile/footing.

- Soil model and SSI:
  - Surrounding soil represented using:
    - Lumped springs and dashpots for radiation damping (Winkler/elastic springs + dashpots).
    - Discrete continuum models (3D/2D solid elements) where higher fidelity is necessary.
  - Interface behavior: contact or connection elements as needed (gap/friction models can be added).

- Loading and analysis:
  - Static gravity and lateral pushover analyses.
  - Linear and nonlinear time-history analyses using user-provided ground motions.
  - Parametric studies over soil stiffness, pier height, and foundation properties.

## Software and dependencies
Minimum requirements:
- Python 3.8+ (recommended 3.9+)
- OpenSeesPy (install via pip)
- NumPy, SciPy, pandas, matplotlib for post-processing
- (Optional) Jupyter Notebook for interactive exploration

Recommended quick install:
```bash
python -m venv venv
source venv/bin/activate       # On Windows use `venv\Scripts\activate`
pip install --upgrade pip
pip install openseespy numpy scipy pandas matplotlib
```

If you maintain a `requirements.txt`:
```text
openseespy
numpy
scipy
pandas
matplotlib
```

## Installation and running examples
1. Clone the repository:
```bash
git clone https://github.com/Dr-Wahab/Bridge-Pier-with-Soil-Structure-Interaction-OpenSeesPy-.git
cd Bridge-Pier-with-Soil-Structure-Interaction-OpenSeesPy-
```

2. Prepare virtual environment and install dependencies (see previous section).

3. Running an example static analysis:
```bash
python examples/run_static_pier.py
```
(Replace `examples/run_static_pier.py` with the actual script in the repository. The example script builds the model, applies gravity and lateral loads, and writes results to `output/`.)

4. Running a time-history analysis:
```bash
python examples/run_timehistory_pier.py --ground-motion data/ground_motion.txt --dt 0.01
```
(Replace script name and ground motion path accordingly. The script should accept input file and time-step options.)

5. Post-processing:
```bash
python postprocess/plot_results.py --results output/results.h5 --variable displacement
```
(Adjust arguments to match your output format; tools here illustrate a recommended workflow.)

## Reproducibility and recommended workflow
- For every simulation, include:
  - Random seeds (if stochastic processes used).
  - Exact OpenSeesPy version (report `openseespy.__version__`).
  - Input ground motion files and any processing scripts used to scale or filter the motions.
  - Parameter sets used in parametric studies (e.g., soil stiffness, damping ratios).
- Use a single script (or Jupyter notebook) to run the entire pipeline from model creation to figures for reproducible figures and tables.
- Store numerical results in human-readable formats (CSV, HDF5) and provide plotting scripts.

## Expected outputs and validation
- Typical outputs:
  - Time histories of pier top displacement and acceleration.
  - Base shear, bending moment diagrams, and foundation forces.
  - Mode shapes and eigenvalues (for modal analyses).
- Validation checks (suggested):
  - Compare fundamental period and mode shapes with analytical estimates (e.g., simple cantilever formulas).
  - For small-amplitude excitations, compare linear response with fixed-base linear elastic solutions.
  - Convergence checks: mesh sensitivity for continuum soil models; step-size sensitivity for time integration.

## File structure (example)
Update this section to reflect exact repository structure.
- examples/
  - run_static_pier.py
  - run_timehistory_pier.py
- models/
  - pier_model.py
  - soil_model.py
- postprocess/
  - plot_results.py
  - extract_signals.py
- data/
  - ground_motion.txt
- output/
  - results.h5
- README.md
- LICENSE

## How to extend or reuse the models
- To change geometry: edit `models/pier_model.py` geometry parameters.
- To change soil properties: edit `models/soil_model.py` or the input parameter file used by your script.
- To add new ground motions: place files into `data/` and reference them in the runner scripts.
- To integrate more advanced soil models (e.g., dynamic 3D continuum), see OpenSeesPy examples and libraries; ensure time-step and mesh resolution adhere to stability criteria.

## Citation and attribution
If you use this repository in academic work, cite it as:
- Dr. Wahab (2026). Bridge Pier with Soil-Structure Interaction (OpenSeesPy). GitHub repository. https://github.com/Dr-Wahab/Bridge-Pier-with-Soil-Structure-Interaction-OpenSeesPy-

Additionally cite OpenSeesPy:
- McKenna, F., Fenves, G. L., & Scott, M. H. (2000). OpenSees — Open System for Earthquake Engineering Simulation. [OpenSeesPy documentation / relevant citation].

Please add DOIs or formal paper references if available.

## Licensing
This repository is provided under the [LICENSE NAME] license. (Replace with your actual license, e.g., MIT, CC-BY, GPLv3.)

## Acknowledgements
Acknowledge funding, collaborators, or data sources here.

## Contact
For questions or contributions, open an issue or contact: Dr-Wahab (GitHub: @Dr-Wahab). Provide an email address if you wish to receive direct correspondence.

## Appendix: Common OpenSeesPy snippets
Model initialization and basic commands:
```python
import openseespy.opensees as ops

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Example: define nodes and fixities
ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 3.0)
ops.fix(1, 1, 1, 1)  # fully fixed base

# Example: define materials
ops.uniaxialMaterial('Concrete01', 1, -30.0, -0.002, -0.01, -0.002)

# Add elements, loads, analysis options, then run
```

Notes:
- Replace the snippets above with the exact sequences used in your scripts.
- Document non-default solver options, integration schemes, and convergence tolerances used in nonlinear analyses.
