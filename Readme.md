# TiO2 Sensor Tool Project

## Overview
This project aims to develop tools for the synthesis, characterization, and methane gas detection using pristine TiO2 nanoparticles at low temperatures.

## Plans

### Plan 1: Computational Simulation Tool (Extended)
**Objective:** Simulate TiO2 nanoparticle synthesis, structure, doping mechanisms, various formations, and methane adsorption/sensing to optimize for low-temperature methane detection.

**Execution Process:**
1. Research and select simulation software (e.g., LAMMPS for MD, Quantum ESPRESSO for DFT).
2. Set up computational environment (install software, configure HPC if needed).
3. Model pristine TiO2 crystal structure and nanoparticle formation.
4. Implement doping mechanisms (e.g., metal dopants like Co, Pd).
5. Explore different formations (sizes, shapes, phases: anatase/rutile).
6. Simulate methane adsorption and sensing response (binding energies, charge transfer, resistance modeling).
7. Compare with existing best performers (literature review and benchmarking).
8. Validate models against experimental data.
9. Create visualization and analysis scripts for optimization.

**Technologies:** Python, ASE, Matplotlib, NumPy, LAMMPS.

**Execution Status:** Extended with doping and formations.
- Pristine TiO2 nanoparticle (181 atoms) and Co-doped (2 Co atoms) generated.
- Nanoparticles of different sizes (22-606 atoms) created.
- Adsorption simulations for pristine and doped systems.
- Binding energies computed (~ -103 eV, similar for both with simple potential).
- Literature review: Best performers include Pd-doped TiO2 thin films (high sensitivity at 80-120°C), noble metal single atoms on TiO2 for low-temp activation.
- Next: Implement DFT/ReaxFF for accurate electronic properties, charge transfer modeling, resistance simulation.

### Plan 2: Data Analysis Tool
**Objective:** Process and analyze experimental characterization data (XRD, SEM, TEM) and sensor response data.

**Execution Process:**
1. Identify data formats and sources.
2. Develop data import and preprocessing modules.
3. Implement analysis algorithms (peak fitting, particle size distribution, response curve fitting).
4. Create visualization dashboards.
5. Add statistical analysis and error estimation.
6. Integrate with database for data storage.

**Technologies:** Python, Pandas, SciPy, Plotly, SQLite.

### Plan 3: Sensor Modeling Tool
**Objective:** Model the gas sensor behavior for methane detection at low temperatures.

**Execution Process:**
1. Research sensor physics and TiO2 gas sensing mechanisms.
2. Develop mathematical models for resistance change vs gas concentration.
3. Implement temperature-dependent models.
4. Create simulation interface for different nanoparticle properties.
5. Validate against experimental sensor data.
6. Generate prediction tools for sensor design.

**Technologies:** Python, NumPy, SciPy, Jupyter.

### Plan 4: Machine Learning Optimization Tool
**Objective:** Use ML to optimize synthesis parameters for best methane detection performance.

**Execution Process:**
1. Collect or generate training data (synthesis params vs sensor performance).
2. Preprocess and feature engineer data.
3. Train ML models (regression for performance prediction).
4. Implement optimization algorithms (Bayesian optimization).
5. Create user interface for parameter suggestion.
6. Validate predictions experimentally.

**Technologies:** Python, Scikit-learn, TensorFlow, Optuna.

### Plan 5: ZnO Nanoparticle Simulation Tool (New Extension)
**Objective:** Simulate ZnO nanoparticle synthesis, structure, doping mechanisms, various formations, and methane adsorption/sensing to optimize for low-temperature methane detection, comparing with TiO2.

**Execution Process:**
1. Research and select simulation software (e.g., LAMMPS for MD, Quantum ESPRESSO for DFT).
2. Set up computational environment (install software, configure HPC if needed).
3. Model pristine ZnO crystal structure (wurtzite phase) and nanoparticle formation.
4. Implement doping mechanisms (e.g., metal dopants like Ag, Cu).
5. Explore different formations (sizes, shapes, phases).
6. Simulate methane adsorption and sensing response (binding energies, charge transfer, resistance modeling).
7. Compare with TiO2 and existing best performers (literature review and benchmarking).
8. Validate models against experimental data.
9. Create visualization and analysis scripts for optimization.

**Technologies:** Python, ASE, Matplotlib, NumPy, LAMMPS.

**Execution Status:** Completed basic MD adsorption simulation; binding energy -1.83 eV.

### Comparison: ZnO vs TiO2 for Methane Sensing

| Aspect | TiO2 | ZnO |
|--------|------|-----|
| **Crystal Structure** | Anatase (tetragonal)/Rutile | Wurtzite (hexagonal) |
| **Bandgap** | ~3.2 eV | ~3.4 eV |
| **Semiconductor Type** | n-type | n-type |
| **Sensing Mechanism** | Oxygen adsorption/desorption; resistance change on gas exposure | Similar, but higher sensitivity due to polar surfaces; better at room temp |
| **Doping Examples** | Pd, Co (catalytic enhancement) | Ag, Cu (improves conductivity and affinity) |
| **Methane Adsorption (Estimated)** | Binding ~ -103 eV (pristine); improved with Pd | Simulated -1.83 eV (334-atom nanoparticle); potentially stronger with real potentials |
| **Advantages** | Stable, tunable phases | Lower operating temp, higher selectivity for reducing gases |
| **Disadvantages** | Needs higher temp for activation | Less stable in humid environments |
| **Simulation Status** | Completed basic MD for TiO2; Pd-doped promising | Planned; expect better low-temp performance |

ZnO may outperform TiO2 in low-temperature methane sensing due to its polar nature, but TiO2 is more phase-diverse for optimization.

## Next Steps
Select a plan to implement or combine elements from multiple plans. Start ZnO simulations by adapting TiO2 scripts.
