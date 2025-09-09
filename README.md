# LHON Penetrance and Prevalence Modeling

## Mathematical Modeling of Leber Hereditary Optic Neuropathy (LHON) Penetrance and Prevalence: A Comprehensive Analysis Integrating Population Genomics and Environmental Factors

**Author:** Dmitrii Smirnov  
**Affiliation:** Institute of Human Genetics, School of Medicine, Technical University of Munich, Munich, Germany  
**Email:** dmitrii.smirnov@tum.de

---

## Project Overview

This repository contains a comprehensive mathematical modeling framework for analyzing the penetrance and prevalence of Leber Hereditary Optic Neuropathy (LHON), a maternally inherited mitochondrial disease. The project addresses the long-standing discrepancy between traditional, disease-cohort-based penetrance estimates and emerging population-genomics-based estimates.

### Key Findings

- **Revised Penetrance Estimates:** Overall penetrance of LHON is approximately 2.44%, significantly lower than the 40-50% often cited in traditional literature
- **Carrier Prevalence:** Approximately 1 in 910 individuals carry a primary LHON mutation
- **Patient Prevalence:** Approximately 1 in 37,308 individuals are clinically affected
- **Risk Stratification:** Individual risk varies dramatically based on genetic and environmental factors

## Repository Structure

```
├── README.md                    # This file
├── paper/                       # Scientific papers and manuscripts
│   ├── lhon_scientific_paper_final.md    # Main scientific paper
│   └── lhon_scientific_paper.md          # Draft version
├── scripts/                     # Python analysis scripts
│   ├── lhon_mathematical_models.py       # Core mathematical models
│   ├── lhon_bayesian_network.py          # Bayesian hierarchical model
│   ├── lhon_sensitivity_analysis.py      # Sensitivity analysis
│   ├── lhon_real_prevalence_analysis.py  # Real prevalence calculations
│   ├── lhon_model_validation.py          # Model validation and calibration
│   ├── lhon_visualizations.py            # Main visualization script
│   ├── lhon_key_findings_visualizations.py    # Key findings plots
│   ├── lhon_real_prevalence_visualizations.py # Prevalence visualizations
│   └── lhon_sensitivity_visualizations.py     # Sensitivity analysis plots
├── data/                        # Generated data and results
│   ├── lhon_liability_model_results.csv       # Liability model outputs
│   ├── lhon_bayesian_model_results.csv        # Bayesian model results
│   ├── lhon_monte_carlo_results.csv           # Monte Carlo simulation data
│   ├── lhon_sensitivity_indices.csv           # Sensitivity analysis results
│   ├── real_prevalence_summary.csv            # Real prevalence estimates
│   └── [additional CSV files]                 # Supporting data files
├── figures/                     # Generated visualizations
│   ├── figure1_penetrance_comparison.png      # Main penetrance comparison
│   ├── figure2_population_prevalence.png      # Population prevalence analysis
│   ├── figure3_environmental_risk.png         # Environmental risk factors
│   ├── figure4_model_validation.png           # Model validation results
│   ├── figure5_bayesian_network.png           # Bayesian network structure
│   ├── key_findings_*.png                     # Key findings visualizations
│   ├── real_prevalence_*.png                  # Real prevalence plots
│   └── sensitivity_*.png                      # Sensitivity analysis plots
└── reports/                     # Analysis reports and summaries
    ├── lhon_sensitivity_analysis_report.md    # Sensitivity analysis report
    ├── lhon_real_prevalence_report.md         # Real prevalence report
    ├── lhon_key_findings_summary.md           # Key findings summary
    └── [additional reports]                   # Supporting reports
```

## Requirements

### Python Dependencies

```bash
pip install numpy pandas matplotlib seaborn scipy networkx
```

### System Requirements

- Python 3.7+
- Minimum 8GB RAM (for Monte Carlo simulations)
- ~2GB disk space for generated data and figures

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/mityasmirnov/LHON_analytics.git
cd LHON_analytics
```

### 2. Install Dependencies

```bash
pip install numpy pandas matplotlib seaborn scipy networkx
```

### 3. Run the Complete Analysis

Execute the scripts in the following order to reproduce the full analysis:

```bash
# 1. Core mathematical models
python scripts/lhon_mathematical_models.py

# 2. Bayesian hierarchical model
python scripts/lhon_bayesian_network.py

# 3. Model validation and calibration
python scripts/lhon_model_validation.py

# 4. Sensitivity analysis
python scripts/lhon_sensitivity_analysis.py

# 5. Real prevalence analysis
python scripts/lhon_real_prevalence_analysis.py

# 6. Generate all visualizations
python scripts/lhon_visualizations.py
python scripts/lhon_key_findings_visualizations.py
python scripts/lhon_real_prevalence_visualizations.py
python scripts/lhon_sensitivity_visualizations.py
```

## Detailed Analysis Workflow

### Step 1: Mathematical Models (`lhon_mathematical_models.py`)

This script implements three core modeling approaches:

1. **Liability-Threshold Model:** Models LHON as a dichotomous trait influenced by multiple continuous risk factors
2. **Bayesian Hierarchical Model:** Accounts for uncertainty in parameters using probabilistic distributions
3. **Monte Carlo Population Simulation:** Simulates 100,000-person populations to estimate epidemiological parameters

**Key Outputs:**
- `data/lhon_liability_model_results.csv`
- `data/lhon_monte_carlo_results.csv`
- `data/lhon_revised_penetrance_estimates.csv`

### Step 2: Bayesian Network Analysis (`lhon_bayesian_network.py`)

Implements a probabilistic graphical model representing conditional dependencies between LHON risk factors.

**Key Outputs:**
- `data/lhon_bayesian_model_results.csv`
- `data/lhon_bayesian_samples.csv`
- `data/lhon_bayesian_penetrance.csv`

### Step 3: Model Validation (`lhon_model_validation.py`)

Validates and calibrates models against empirical data from population-based studies.

**Key Outputs:**
- `data/lhon_model_validation_report.json`
- Calibrated model parameters

### Step 4: Sensitivity Analysis (`lhon_sensitivity_analysis.py`)

Identifies key parameters that most influence model outcomes through:
- One-at-a-time (OAT) analysis
- Monte Carlo sensitivity analysis
- Scenario analysis

**Key Outputs:**
- `data/lhon_sensitivity_indices.csv`
- `data/lhon_scenario_analysis.csv`
- `reports/lhon_sensitivity_analysis_report.md`

### Step 5: Real Prevalence Analysis (`lhon_real_prevalence_analysis.py`)

Calculates real-world prevalence estimates using calibrated models and gnomAD population data.

**Key Outputs:**
- `data/real_prevalence_summary.csv`
- `data/real_carrier_prevalence.csv`
- `data/real_patient_prevalence.csv`
- `reports/lhon_real_prevalence_report.md`

### Step 6: Visualization Generation

Multiple visualization scripts generate publication-quality figures:

- `lhon_visualizations.py`: Main figures for the scientific paper
- `lhon_key_findings_visualizations.py`: Summary of key findings
- `lhon_real_prevalence_visualizations.py`: Real prevalence analysis plots
- `lhon_sensitivity_visualizations.py`: Sensitivity analysis visualizations

## Key Data Sources

### Population Genomics Data
- **gnomAD v4.0:** Carrier frequencies for primary LHON mutations (m.11778G>A, m.14484T>C, m.3460G>A)
- **Population size:** >800,000 individuals

### Clinical Data
- **Kirkman et al. (2009):** Environmental risk factors (smoking, alcohol)
- **Watson et al. (2022):** Population-based penetrance estimates
- **Recent epidemiological studies:** Population prevalence data

### Genetic Modifiers
- **mtDNA haplogroups:** Particularly haplogroup J effects
- **Nuclear modifiers:** DNAJC30 and other nuclear genes

## Model Parameters

### Core Parameters
- **Liability threshold:** 5.0 (calibrated)
- **Male effect (OR):** 7.11
- **Heavy smoking (OR):** 3.16
- **Haplogroup J (OR):** 2.00

### Mutation-Specific Parameters
- **m.11778G>A:** Base liability = -2.5, Average penetrance = 4.56%
- **m.14484T>C:** Base liability = -4.0, Average penetrance = 0.67%
- **m.3460G>A:** Base liability = -1.5, Average penetrance = 17.15%

## Reproducing Specific Figures

### Figure 1: Penetrance Comparison
```bash
python scripts/lhon_visualizations.py
# Generates: figures/figure1_penetrance_comparison.png
```

### Figure 2: Population Prevalence
```bash
python scripts/lhon_visualizations.py
# Generates: figures/figure2_population_prevalence.png
```

### Key Findings Summary
```bash
python scripts/lhon_key_findings_visualizations.py
# Generates: figures/key_findings_*.png
```

### Real Prevalence Analysis
```bash
python scripts/lhon_real_prevalence_visualizations.py
# Generates: figures/real_prevalence_*.png
```

### Sensitivity Analysis
```bash
python scripts/lhon_sensitivity_visualizations.py
# Generates: figures/sensitivity_*.png
```

## Validation and Quality Control

### Model Validation Checks
1. **Carrier frequency validation:** Model predictions vs. gnomAD data
2. **Population prevalence validation:** Model outputs vs. epidemiological studies
3. **Penetrance validation:** Calibrated estimates vs. literature values

### Expected Results
- **Total carrier prevalence:** ~109.9 per 100,000 (matches gnomAD)
- **Total patient prevalence:** ~2.68 per 100,000 (within literature range: 0.79-3.23)
- **Overall penetrance:** ~2.44% (calibrated to population data)

## Troubleshooting

### Common Issues

1. **Memory errors during Monte Carlo simulation:**
   - Reduce population size in simulation parameters
   - Increase available RAM or use smaller batch sizes

2. **Missing dependencies:**
   ```bash
   pip install --upgrade numpy pandas matplotlib seaborn scipy networkx
   ```

3. **Figure generation errors:**
   - Ensure all data files are generated before running visualization scripts
   - Check that the `figures/` directory exists and is writable

4. **Inconsistent results:**
   - Monte Carlo simulations include random components
   - Set random seeds for reproducible results
   - Run multiple iterations and average results

### Performance Optimization

- **For faster execution:** Reduce Monte Carlo iterations (default: 1,000)
- **For higher precision:** Increase simulation population size (default: 100,000)
- **For memory efficiency:** Process data in smaller batches

## Citation

If you use this work in your research, please cite:

```
Smirnov, D. (2025). Insights into prevalence and penetrance bias estimations for 
Leber's hereditary optic neuropathy LHON. Institute of Human Genetics, 
School of Medicine, Technical University of Munich.
```

## License

This project is released under the MIT License. See LICENSE file for details.

## Contact

**Dmitrii Smirnov**  
Institute of Human Genetics  
School of Medicine  
Technical University of Munich  
Munich, Germany  
Email: dmitrii.smirnov@tum.de

## Acknowledgments

- gnomAD Consortium for population genomics data
- LHON research community for clinical and epidemiological data
- Technical University of Munich for computational resources

---

*Last updated: September 9, 2025*

