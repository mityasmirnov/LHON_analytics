#!/usr/bin/env python3
"""
LHON Visualization and Figure Generation
Creates publication-quality figures for the scientific paper
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle, FancyBboxPatch
import matplotlib.patches as mpatches
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set style for publication-quality figures
plt.style.use('default')
sns.set_palette("husl")

class LHONVisualizer:
    """
    Creates comprehensive visualizations for LHON modeling results
    """
    
    def __init__(self):
        """Initialize with data and styling parameters"""
        
        # Load data
        self.load_data()
        
        # Set figure parameters
        self.fig_params = {
            'figure.figsize': (12, 8),
            'font.size': 12,
            'axes.titlesize': 14,
            'axes.labelsize': 12,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 10,
            'figure.dpi': 300
        }
        
        plt.rcParams.update(self.fig_params)
        
        # Color scheme
        self.colors = {
            'traditional': '#d62728',  # Red
            'population': '#2ca02c',   # Green
            'modeled': '#1f77b4',      # Blue
            'gnomad': '#ff7f0e',       # Orange
            'environmental': '#9467bd', # Purple
            'male': '#1f77b4',         # Blue
            'female': '#e377c2'        # Pink
        }
    
    def load_data(self):
        """Load modeling results"""
        
        try:
            self.liability_results = pd.read_csv('/home/ubuntu/lhon_liability_model_results.csv')
            self.bayesian_penetrance = pd.read_csv('/home/ubuntu/lhon_bayesian_penetrance.csv')
            self.monte_carlo_results = pd.read_csv('/home/ubuntu/lhon_monte_carlo_results.csv')
            self.revised_estimates = pd.read_csv('/home/ubuntu/lhon_revised_penetrance_estimates.csv')
            
            print("Successfully loaded visualization data")
            
        except FileNotFoundError as e:
            print(f"Warning: Could not load some data files: {e}")
    
    def create_penetrance_comparison_figure(self):
        """Create figure comparing traditional vs population-based penetrance"""
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # Data for comparison
        mutations = ['11778G>A', '14484T>C', '3460G>A']
        traditional_penetrance = [43, 65, 78]  # Literature values
        gnomad_calculated = [2.86, 0.57, 10.56]  # From our analysis
        watson_overall = 1.11  # Population-based estimate
        
        # Panel A: Penetrance by mutation
        x = np.arange(len(mutations))
        width = 0.35
        
        bars1 = ax1.bar(x - width/2, traditional_penetrance, width, 
                       label='Traditional (Disease Cohorts)', color=self.colors['traditional'], alpha=0.8)
        bars2 = ax1.bar(x + width/2, gnomad_calculated, width, 
                       label='gnomAD-informed (Population)', color=self.colors['population'], alpha=0.8)
        
        ax1.axhline(y=watson_overall, color=self.colors['gnomad'], linestyle='--', 
                   label='Watson et al. (Overall)', linewidth=2)
        
        ax1.set_xlabel('mtDNA Mutation')
        ax1.set_ylabel('Penetrance (%)')
        ax1.set_title('A. Penetrance Estimates: Traditional vs Population-Based')
        ax1.set_xticks(x)
        ax1.set_xticklabels(mutations)
        ax1.legend()
        ax1.set_ylim(0, 85)
        
        # Add value labels on bars
        for bar in bars1:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{height:.0f}%', ha='center', va='bottom', fontweight='bold')
        
        for bar in bars2:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{height:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        # Panel B: Overestimation ratios
        overestimation_ratios = [15.0, 114.0, 7.4]  # From our analysis
        
        bars3 = ax2.bar(mutations, overestimation_ratios, color=self.colors['traditional'], alpha=0.7)
        ax2.axhline(y=1, color='black', linestyle='-', alpha=0.5, label='Perfect Agreement')
        
        ax2.set_xlabel('mtDNA Mutation')
        ax2.set_ylabel('Overestimation Ratio\n(Traditional / Population-based)')
        ax2.set_title('B. Traditional Penetrance Overestimation')
        ax2.set_yscale('log')
        ax2.legend()
        
        # Add value labels
        for i, (bar, ratio) in enumerate(zip(bars3, overestimation_ratios)):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{ratio:.0f}x', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/figure1_penetrance_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created Figure 1: Penetrance Comparison")
    
    def create_population_prevalence_figure(self):
        """Create figure showing population prevalence estimates"""
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # Panel A: Prevalence by study
        studies = ['Madrid\n2024', 'Watson\nAustralia', 'Finland\nHistorical', 'UK\nHistorical']
        prevalences = [0.79, 1.46, 2.0, 3.23]
        colors_studies = [self.colors['population']] * len(studies)
        
        bars = ax1.bar(studies, prevalences, color=colors_studies, alpha=0.8)
        ax1.axhline(y=np.mean(prevalences), color='red', linestyle='--', 
                   label=f'Average: {np.mean(prevalences):.2f}', linewidth=2)
        
        ax1.set_ylabel('Prevalence (per 100,000)')
        ax1.set_title('A. Population Prevalence by Study')
        ax1.legend()
        ax1.set_ylim(0, 3.5)
        
        # Add value labels
        for bar, prev in zip(bars, prevalences):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                    f'{prev:.2f}', ha='center', va='bottom', fontweight='bold')
        
        # Panel B: Carrier frequency vs prevalence
        mutations = ['11778G>A', '14484T>C', '3460G>A']
        carrier_freqs = [42.54, 65.58, 1.77]  # per 100,000 from gnomAD
        expected_cases_traditional = [18.3, 42.6, 1.4]  # If traditional penetrance
        observed_cases = [1.2, 0.4, 0.2]  # Estimated from prevalence studies
        
        x = np.arange(len(mutations))
        width = 0.25
        
        ax2.bar(x - width, carrier_freqs, width, label='Carriers (gnomAD)', 
               color=self.colors['gnomad'], alpha=0.8)
        ax2.bar(x, expected_cases_traditional, width, label='Expected Cases (Traditional)', 
               color=self.colors['traditional'], alpha=0.8)
        ax2.bar(x + width, observed_cases, width, label='Observed Cases (Population)', 
               color=self.colors['population'], alpha=0.8)
        
        ax2.set_xlabel('mtDNA Mutation')
        ax2.set_ylabel('Frequency (per 100,000)')
        ax2.set_title('B. Carrier Frequency vs Disease Prevalence')
        ax2.set_xticks(x)
        ax2.set_xticklabels(mutations)
        ax2.legend()
        ax2.set_yscale('log')
        ax2.set_ylim(0.1, 100)
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/figure2_population_prevalence.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created Figure 2: Population Prevalence")
    
    def create_environmental_risk_figure(self):
        """Create figure showing environmental risk factors"""
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # Panel A: Odds ratios for environmental factors
        factors = ['Male Sex', 'Heavy\nSmoking', 'Light\nSmoking', 'Heavy\nAlcohol', 'Light\nAlcohol']
        odds_ratios = [7.11, 3.16, 1.54, 3.27, 1.01]
        ci_lower = [4.58, 1.82, 0.99, 1.44, 0.45]
        ci_upper = [11.03, 5.50, 2.41, 7.46, 2.24]
        
        colors_factors = [self.colors['male'], self.colors['environmental'], 
                         self.colors['environmental'], self.colors['environmental'], 
                         self.colors['environmental']]
        
        bars = ax1.bar(factors, odds_ratios, color=colors_factors, alpha=0.8)
        
        # Add confidence intervals
        for i, (bar, lower, upper) in enumerate(zip(bars, ci_lower, ci_upper)):
            ax1.errorbar(bar.get_x() + bar.get_width()/2., bar.get_height(),
                        yerr=[[bar.get_height() - lower], [upper - bar.get_height()]],
                        fmt='none', color='black', capsize=5, capthick=2)
        
        ax1.axhline(y=1, color='black', linestyle='-', alpha=0.5, label='No Effect')
        ax1.set_ylabel('Odds Ratio (95% CI)')
        ax1.set_title('A. Environmental Risk Factors (Kirkman et al.)')
        ax1.set_yscale('log')
        ax1.set_ylim(0.3, 15)
        ax1.legend()
        
        # Add value labels
        for bar, or_val in zip(bars, odds_ratios):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{or_val:.2f}', ha='center', va='bottom', fontweight='bold')
        
        # Panel B: Penetrance by risk profile
        risk_profiles = ['Female\nNon-smoker', 'Male\nNon-smoker', 'Male\nLight smoker', 'Male\nHeavy smoker']
        penetrance_11778 = [9.7, 74.6, 85.2, 98.1]  # From liability model
        
        colors_risk = [self.colors['female'], self.colors['male'], 
                      self.colors['male'], self.colors['male']]
        
        bars2 = ax2.bar(risk_profiles, penetrance_11778, color=colors_risk, alpha=0.8)
        
        ax2.set_ylabel('Penetrance (%)')
        ax2.set_title('B. 11778G>A Penetrance by Risk Profile')
        ax2.set_ylim(0, 105)
        
        # Add value labels
        for bar, pen in zip(bars2, penetrance_11778):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{pen:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/figure3_environmental_risk.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created Figure 3: Environmental Risk Factors")
    
    def create_bayesian_network_diagram(self):
        """Create Bayesian network structure diagram using manus-render-diagram"""
        
        # Create PlantUML diagram
        plantuml_content = """
@startuml
!theme plain
skinparam backgroundColor white
skinparam defaultFontSize 12

package "Demographic Factors" {
  [Population] as Pop
  [Sex] as Sex
  [Age] as Age
  [Haplogroup] as Hap
}

package "Genetic Factors" {
  [mtDNA Mutation] as mtDNA
  [Nuclear Variants] as Nuclear
}

package "Environmental Factors" {
  [Smoking] as Smoke
  [Alcohol] as Alcohol
  [Environmental Stress] as EnvStress
}

package "Pathophysiology" {
  [Mitochondrial Function] as MitoFunc
  [Oxidative Stress] as OxStress
  [Liability] as Liability
}

package "Clinical Outcomes" {
  [LHON Phenotype] as Phenotype
  [Recovery] as Recovery
}

Pop --> Sex
Pop --> Age
Pop --> Hap
Pop --> mtDNA
Pop --> Nuclear
Pop --> EnvStress

Sex --> Smoke
Sex --> Alcohol

mtDNA --> MitoFunc
Hap --> MitoFunc
Nuclear --> MitoFunc

Smoke --> OxStress
Alcohol --> OxStress
EnvStress --> OxStress

MitoFunc --> Liability
OxStress --> Liability
Sex --> Liability
Age --> Liability

Liability --> Phenotype
Phenotype --> Recovery
Age --> Recovery
mtDNA --> Recovery

@enduml
"""
        
        # Save PlantUML file
        with open('/home/ubuntu/bayesian_network.puml', 'w') as f:
            f.write(plantuml_content)
        
        print("Created Bayesian Network diagram source")
    
    def create_model_validation_figure(self):
        """Create figure showing model validation results"""
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        # Panel A: Penetrance validation
        models = ['Traditional\n(Literature)', 'Bayesian\n(Uncalibrated)', 'Watson et al.\n(Population)', 'Theoretical\n(gnomAD + Watson)']
        penetrance_values = [43, 44, 1.11, 1.21]  # Overall penetrance estimates
        colors_models = [self.colors['traditional'], self.colors['modeled'], 
                        self.colors['population'], self.colors['gnomad']]
        
        bars1 = ax1.bar(models, penetrance_values, color=colors_models, alpha=0.8)
        ax1.set_ylabel('Overall Penetrance (%)')
        ax1.set_title('A. Model Penetrance Validation')
        ax1.set_yscale('log')
        ax1.set_ylim(0.5, 50)
        
        # Add value labels
        for bar, val in zip(bars1, penetrance_values):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        # Panel B: Prevalence validation
        prevalence_models = ['Observed\n(Literature)', 'Monte Carlo\n(Uncalibrated)', 'Theoretical\n(gnomAD + Watson)']
        prevalence_values = [1.87, 32.15, 1.21]
        colors_prev = [self.colors['population'], self.colors['modeled'], self.colors['gnomad']]
        
        bars2 = ax2.bar(prevalence_models, prevalence_values, color=colors_prev, alpha=0.8)
        ax2.set_ylabel('Population Prevalence (per 100,000)')
        ax2.set_title('B. Population Prevalence Validation')
        ax2.set_yscale('log')
        ax2.set_ylim(0.5, 50)
        
        # Add value labels
        for bar, val in zip(bars2, prevalence_values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.1f}', ha='center', va='bottom', fontweight='bold')
        
        # Panel C: Sex ratio validation
        sex_ratios = ['Literature\n(Kirkman et al.)', 'Bayesian Model\n(Uncalibrated)']
        ratio_values = [7.11, 1.6]
        colors_sex = [self.colors['traditional'], self.colors['modeled']]
        
        bars3 = ax3.bar(sex_ratios, ratio_values, color=colors_sex, alpha=0.8)
        ax3.set_ylabel('Male:Female Penetrance Ratio')
        ax3.set_title('C. Sex Ratio Validation')
        ax3.set_ylim(0, 8)
        
        # Add value labels
        for bar, val in zip(bars3, ratio_values):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{val:.1f}', ha='center', va='bottom', fontweight='bold')
        
        # Panel D: Overestimation by mutation
        mutations = ['11778G>A', '14484T>C', '3460G>A']
        overestimation = [15.0, 114.0, 7.4]
        
        bars4 = ax4.bar(mutations, overestimation, color=self.colors['traditional'], alpha=0.8)
        ax4.set_ylabel('Overestimation Ratio\n(Traditional / Population-based)')
        ax4.set_title('D. Traditional Penetrance Overestimation by Mutation')
        ax4.set_yscale('log')
        ax4.set_ylim(1, 200)
        
        # Add value labels
        for bar, val in zip(bars4, overestimation):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.0f}x', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/figure4_model_validation.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created Figure 4: Model Validation")
    
    def create_summary_table(self):
        """Create comprehensive summary table"""
        
        # Create summary data
        summary_data = {
            'Parameter': [
                'Carrier Frequency (per 100,000)',
                'Population Prevalence (per 100,000)', 
                'Overall Penetrance (%)',
                'Male Penetrance (%)',
                'Female Penetrance (%)',
                'Male:Female Ratio',
                'Heavy Smoker OR',
                'Male Sex OR',
                'Recovery Rate 11778G>A (%)',
                'Recovery Rate 14484T>C (%)',
                'Recovery Rate 3460G>A (%)'
            ],
            'Traditional Literature': [
                '22-74', '1.5-4.6', '12.5-50', '50-78', '10-32', '5-8', 
                'Not quantified', 'Not quantified', '4', '37-65', '20-78'
            ],
            'Population Studies': [
                '110 (gnomAD)', '0.79-3.23', '1.11 (0.5-3.4)', 'Not specified', 
                'Not specified', '7.11', '3.16', '7.11', '4', '37', '20'
            ],
            'Our Models': [
                '111.6', '1.21 (theoretical)', '1.11-44 (calibrated)', 
                '51.7 (uncalibrated)', '33.3 (uncalibrated)', '1.6 (uncalibrated)',
                '3.16', '7.11', '0', '56.2', '0'
            ],
            'Discrepancy Ratio': [
                '1.5-5x higher', '1-17x higher', '7-114x higher',
                'Variable', 'Variable', '1-4x higher',
                'Consistent', 'Consistent', 'Consistent', 'Consistent', 'Variable'
            ]
        }
        
        summary_df = pd.DataFrame(summary_data)
        
        # Save as CSV
        summary_df.to_csv('/home/ubuntu/table1_summary_parameters.csv', index=False)
        
        print("Created Table 1: Summary Parameters")
        
        return summary_df
    
    def create_all_figures(self):
        """Create all figures and tables"""
        
        print("LHON VISUALIZATION GENERATION")
        print("=" * 35)
        
        # Create figures
        self.create_penetrance_comparison_figure()
        self.create_population_prevalence_figure()
        self.create_environmental_risk_figure()
        self.create_bayesian_network_diagram()
        self.create_model_validation_figure()
        
        # Create tables
        summary_table = self.create_summary_table()
        
        print("\nAll visualizations created successfully!")
        print("\nGenerated files:")
        print("- figure1_penetrance_comparison.png")
        print("- figure2_population_prevalence.png") 
        print("- figure3_environmental_risk.png")
        print("- figure4_model_validation.png")
        print("- bayesian_network.puml")
        print("- table1_summary_parameters.csv")
        
        return {
            'figures_created': 4,
            'tables_created': 1,
            'summary_table': summary_table
        }

def main():
    """Generate all visualizations"""
    
    visualizer = LHONVisualizer()
    results = visualizer.create_all_figures()
    
    return results

if __name__ == "__main__":
    results = main()

