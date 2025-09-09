#!/usr/bin/env python3
"""
LHON Real Prevalence Visualizations
Creates comprehensive visualizations of real prevalence analysis results
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

# Set style for publication-quality figures
plt.style.use('default')
sns.set_palette("husl")

class LHONRealPrevalenceVisualizer:
    """
    Creates visualizations for LHON real prevalence analysis results
    """
    
    def __init__(self):
        """Initialize and load real prevalence data"""
        
        # Load data
        self.load_data()
        
        # Set figure parameters
        self.fig_params = {
            'figure.figsize': (16, 12),
            'font.size': 12,
            'axes.titlesize': 16,
            'axes.labelsize': 14,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'figure.dpi': 300
        }
        
        plt.rcParams.update(self.fig_params)
        
        # Color scheme
        self.colors = {
            'carriers': '#2ca02c',        # Green
            'patients': '#d62728',        # Red
            'male': '#1f77b4',           # Blue
            'female': '#e377c2',         # Pink
            '11778': '#ff7f0e',          # Orange
            '14484': '#9467bd',          # Purple
            '3460': '#8c564b',           # Brown
            'literature': '#17becf',      # Cyan
            'model': '#bcbd22'           # Olive
        }
    
    def load_data(self):
        """Load real prevalence analysis results"""
        
        try:
            self.carrier_prevalence = pd.read_csv('/home/ubuntu/real_carrier_prevalence.csv')
            self.patient_prevalence = pd.read_csv('/home/ubuntu/real_patient_prevalence.csv')
            self.penetrance_subgroups = pd.read_csv('/home/ubuntu/real_penetrance_subgroups.csv')
            self.sex_specific = pd.read_csv('/home/ubuntu/real_sex_specific_analysis.csv')
            self.summary = pd.read_csv('/home/ubuntu/real_prevalence_summary.csv')
            
            print("Successfully loaded real prevalence analysis data")
            
        except FileNotFoundError as e:
            print(f"Error loading data: {e}")
            raise
    
    def create_carrier_vs_patient_overview(self):
        """Create overview comparing carrier and patient prevalence"""
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
        
        # Panel A: Carrier vs Patient Prevalence by Mutation
        mutations = self.carrier_prevalence['mutation']
        carrier_freq = self.carrier_prevalence['frequency_per_100k']
        patient_freq = self.patient_prevalence['patient_prevalence_per_100k']
        
        x = np.arange(len(mutations))
        width = 0.35
        
        bars1 = ax1.bar(x - width/2, carrier_freq, width, label='Carriers', 
                       color=self.colors['carriers'], alpha=0.8)
        bars2 = ax1.bar(x + width/2, patient_freq, width, label='Patients', 
                       color=self.colors['patients'], alpha=0.8)
        
        ax1.set_xlabel('mtDNA Mutation')
        ax1.set_ylabel('Frequency (per 100,000)')
        ax1.set_title('A. Carrier vs Patient Prevalence by Mutation')
        ax1.set_xticks(x)
        ax1.set_xticklabels(mutations)
        ax1.legend()
        ax1.set_yscale('log')
        ax1.set_ylim(0.1, 100)
        
        # Add value labels
        for bar, val in zip(bars1, carrier_freq):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.1f}', ha='center', va='bottom', fontweight='bold')
        
        for bar, val in zip(bars2, patient_freq):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.2f}', ha='center', va='bottom', fontweight='bold')
        
        # Panel B: Penetrance by Mutation
        penetrance_values = self.patient_prevalence['average_penetrance']
        
        bars3 = ax2.bar(mutations, penetrance_values, 
                       color=[self.colors['11778'], self.colors['14484'], self.colors['3460']], alpha=0.8)
        
        ax2.set_xlabel('mtDNA Mutation')
        ax2.set_ylabel('Average Penetrance (%)')
        ax2.set_title('B. Average Penetrance by Mutation')
        
        # Add value labels
        for bar, val in zip(bars3, penetrance_values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.3,
                    f'{val:.2f}%', ha='center', va='bottom', fontweight='bold')
        
        # Panel C: Population Impact Visualization
        # Create a population pyramid-style visualization
        total_population = 100000
        total_carriers = self.summary['total_carrier_prevalence_per_100k'].iloc[0]
        total_patients = self.summary['total_patient_prevalence_per_100k'].iloc[0]
        unaffected = total_population - total_carriers
        
        categories = ['Unaffected', 'Carriers\n(Asymptomatic)', 'Patients\n(Affected)']
        frequencies = [unaffected, total_carriers - total_patients, total_patients]
        colors_pop = ['lightgray', self.colors['carriers'], self.colors['patients']]
        
        # Create pie chart
        wedges, texts, autotexts = ax3.pie(frequencies, labels=categories, colors=colors_pop, 
                                          autopct='%1.3f%%', startangle=90)
        
        ax3.set_title('C. Population Distribution (per 100,000)')
        
        # Panel D: Risk Ratios
        carrier_to_patient_ratios = []
        for i, row in self.patient_prevalence.iterrows():
            carrier_freq = self.carrier_prevalence.iloc[i]['frequency_per_100k']
            patient_freq = row['patient_prevalence_per_100k']
            ratio = carrier_freq / patient_freq if patient_freq > 0 else 0
            carrier_to_patient_ratios.append(ratio)
        
        bars4 = ax4.bar(mutations, carrier_to_patient_ratios, 
                       color=[self.colors['11778'], self.colors['14484'], self.colors['3460']], alpha=0.8)
        
        ax4.set_xlabel('mtDNA Mutation')
        ax4.set_ylabel('Carrier:Patient Ratio')
        ax4.set_title('D. Carrier to Patient Ratios by Mutation')
        ax4.set_yscale('log')
        
        # Add value labels
        for bar, val in zip(bars4, carrier_to_patient_ratios):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.0f}:1', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/real_prevalence_overview.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created: real_prevalence_overview.png")
    
    def create_penetrance_risk_stratification(self):
        """Create detailed penetrance and risk stratification visualization"""
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
        
        # Panel A: Penetrance by Risk Scenario
        scenarios = self.penetrance_subgroups['scenario']
        penetrances = self.penetrance_subgroups['penetrance_percent']
        
        # Color code by risk level
        colors = []
        for scenario in scenarios:
            if 'heavy smoker + alcohol' in scenario:
                colors.append(self.colors['patients'])  # Highest risk - red
            elif 'heavy smoker' in scenario:
                colors.append('#ff7f0e')  # High risk - orange
            elif 'Male' in scenario and '3460G>A' in scenario:
                colors.append('#ff7f0e')  # High risk - orange
            elif 'Female' in scenario:
                colors.append(self.colors['female'])  # Female - pink
            else:
                colors.append(self.colors['male'])  # Male baseline - blue
        
        bars1 = ax1.barh(range(len(scenarios)), penetrances, color=colors, alpha=0.8)
        
        ax1.set_yticks(range(len(scenarios)))
        ax1.set_yticklabels([s.replace(' ', '\n') for s in scenarios], fontsize=10)
        ax1.set_xlabel('Penetrance (%)')
        ax1.set_title('A. Penetrance by Risk Scenario')
        ax1.set_xscale('log')
        ax1.set_xlim(0.001, 50)
        
        # Add value labels
        for bar, val in zip(bars1, penetrances):
            width = bar.get_width()
            ax1.text(width * 1.1, bar.get_y() + bar.get_height()/2,
                    f'{val:.2f}%', ha='left', va='center', fontweight='bold')
        
        # Panel B: Sex-Specific Analysis
        sex_data = self.sex_specific
        
        metrics = ['Carrier\nFrequency', 'Patient\nFrequency', 'Overall\nPenetrance']
        male_values = [sex_data[sex_data['sex'] == 'Male']['carrier_frequency_per_100k'].iloc[0],
                      sex_data[sex_data['sex'] == 'Male']['patient_frequency_per_100k'].iloc[0],
                      sex_data[sex_data['sex'] == 'Male']['overall_penetrance_percent'].iloc[0]]
        female_values = [sex_data[sex_data['sex'] == 'Female']['carrier_frequency_per_100k'].iloc[0],
                        sex_data[sex_data['sex'] == 'Female']['patient_frequency_per_100k'].iloc[0],
                        sex_data[sex_data['sex'] == 'Female']['overall_penetrance_percent'].iloc[0]]
        
        x = np.arange(len(metrics))
        width = 0.35
        
        bars2 = ax2.bar(x - width/2, male_values, width, label='Male', 
                       color=self.colors['male'], alpha=0.8)
        bars3 = ax2.bar(x + width/2, female_values, width, label='Female', 
                       color=self.colors['female'], alpha=0.8)
        
        ax2.set_xlabel('Metrics')
        ax2.set_ylabel('Value')
        ax2.set_title('B. Sex-Specific Prevalence and Penetrance')
        ax2.set_xticks(x)
        ax2.set_xticklabels(metrics)
        ax2.legend()
        ax2.set_yscale('log')
        ax2.set_ylim(0.1, 200)
        
        # Add value labels
        for bar, val in zip(bars2, male_values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.2f}', ha='center', va='bottom', fontweight='bold', fontsize=9)
        
        for bar, val in zip(bars3, female_values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.2f}', ha='center', va='bottom', fontweight='bold', fontsize=9)
        
        # Panel C: Risk Stratification Heatmap
        # Create risk matrix for key scenarios
        risk_scenarios = ['Male\nnon-smoker', 'Female\nnon-smoker', 'Male\nheavy smoker', 
                         'Male heavy\nsmoker + alcohol']
        mutations_heat = ['11778G>A', '14484T>C', '3460G>A']
        
        # Extract penetrance values for heatmap
        risk_matrix = np.array([
            [2.82, 0.15, 17.59],    # Male non-smoker
            [0.04, 0.00, 0.81],     # Female non-smoker
            [14.83, 0.59, 43.93],   # Male heavy smoker (estimated for 14484 and 3460)
            [43.93, 1.77, 87.86]    # Male heavy smoker + alcohol (estimated)
        ])
        
        im = ax3.imshow(risk_matrix, cmap='Reds', aspect='auto', vmin=0, vmax=50)
        
        ax3.set_xticks(range(len(mutations_heat)))
        ax3.set_yticks(range(len(risk_scenarios)))
        ax3.set_xticklabels(mutations_heat)
        ax3.set_yticklabels(risk_scenarios)
        ax3.set_title('C. Risk Stratification Matrix (%)')
        
        # Add text annotations
        for i in range(len(risk_scenarios)):
            for j in range(len(mutations_heat)):
                text = ax3.text(j, i, f'{risk_matrix[i, j]:.1f}%',
                               ha="center", va="center", 
                               color="white" if risk_matrix[i, j] > 25 else "black",
                               fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax3, fraction=0.046, pad=0.04)
        cbar.set_label('Penetrance (%)')
        
        # Panel D: Model Validation
        validation_metrics = ['Population\nPrevalence', 'Overall\nPenetrance', 'Male:Female\nRatio']
        
        model_values = [
            self.summary['total_patient_prevalence_per_100k'].iloc[0],
            self.summary['overall_penetrance_percent'].iloc[0],
            self.summary['male_female_patient_ratio'].iloc[0]
        ]
        
        literature_values = [1.87, 1.1, 7.11]  # Literature benchmarks
        
        x = np.arange(len(validation_metrics))
        width = 0.35
        
        bars4 = ax4.bar(x - width/2, model_values, width, label='Model', 
                       color=self.colors['model'], alpha=0.8)
        bars5 = ax4.bar(x + width/2, literature_values, width, label='Literature', 
                       color=self.colors['literature'], alpha=0.8)
        
        ax4.set_xlabel('Validation Metrics')
        ax4.set_ylabel('Value')
        ax4.set_title('D. Model Validation Against Literature')
        ax4.set_xticks(x)
        ax4.set_xticklabels(validation_metrics)
        ax4.legend()
        ax4.set_yscale('log')
        ax4.set_ylim(0.5, 30)
        
        # Add value labels
        for bar, val in zip(bars4, model_values):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.2f}', ha='center', va='bottom', fontweight='bold', fontsize=9)
        
        for bar, val in zip(bars5, literature_values):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.2f}', ha='center', va='bottom', fontweight='bold', fontsize=9)
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/real_prevalence_risk_stratification.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created: real_prevalence_risk_stratification.png")
    
    def create_population_impact_analysis(self):
        """Create population impact and clinical implications visualization"""
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
        
        # Panel A: Population Pyramid
        total_pop = 100000
        total_carriers = self.summary['total_carrier_prevalence_per_100k'].iloc[0]
        total_patients = self.summary['total_patient_prevalence_per_100k'].iloc[0]
        
        # Create population pyramid data
        categories = ['General\nPopulation', 'Carriers\n(Asymptomatic)', 'Patients\n(Affected)']
        frequencies = [total_pop - total_carriers, total_carriers - total_patients, total_patients]
        
        # Create horizontal bar chart
        y_pos = np.arange(len(categories))
        colors_pyramid = ['lightgray', self.colors['carriers'], self.colors['patients']]
        
        bars1 = ax1.barh(y_pos, frequencies, color=colors_pyramid, alpha=0.8)
        
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(categories)
        ax1.set_xlabel('Population (per 100,000)')
        ax1.set_title('A. Population Distribution')
        ax1.set_xscale('log')
        
        # Add value labels
        for bar, val in zip(bars1, frequencies):
            width = bar.get_width()
            ax1.text(width * 1.1, bar.get_y() + bar.get_height()/2,
                    f'{val:.1f}', ha='left', va='center', fontweight='bold')
        
        # Panel B: Clinical Burden by Mutation
        mutations = self.patient_prevalence['mutation']
        patient_frequencies = self.patient_prevalence['patient_prevalence_per_100k']
        
        # Calculate expected cases in different population sizes
        pop_sizes = [100000, 1000000, 10000000]  # 100k, 1M, 10M
        
        bottom = np.zeros(len(pop_sizes))
        colors_burden = [self.colors['11778'], self.colors['14484'], self.colors['3460']]
        
        for i, (mutation, freq) in enumerate(zip(mutations, patient_frequencies)):
            cases = [freq * pop / 100000 for pop in pop_sizes]
            ax2.bar(range(len(pop_sizes)), cases, bottom=bottom, 
                   label=mutation, color=colors_burden[i], alpha=0.8)
            bottom += cases
        
        ax2.set_xlabel('Population Size')
        ax2.set_ylabel('Expected LHON Cases')
        ax2.set_title('B. Expected Clinical Burden by Population Size')
        ax2.set_xticks(range(len(pop_sizes)))
        ax2.set_xticklabels(['100,000', '1,000,000', '10,000,000'])
        ax2.legend()
        
        # Panel C: Penetrance vs Literature Comparison
        literature_penetrance = {
            '11778G>A': 50,  # Literature average
            '14484T>C': 15,  # Literature average
            '3460G>A': 80   # Literature average
        }
        
        model_penetrance = self.patient_prevalence['average_penetrance'].values
        lit_penetrance = [literature_penetrance[mut] for mut in mutations]
        
        x = np.arange(len(mutations))
        width = 0.35
        
        bars3 = ax3.bar(x - width/2, model_penetrance, width, label='Model (Calibrated)', 
                       color=self.colors['model'], alpha=0.8)
        bars4 = ax3.bar(x + width/2, lit_penetrance, width, label='Literature (Traditional)', 
                       color=self.colors['literature'], alpha=0.8)
        
        ax3.set_xlabel('mtDNA Mutation')
        ax3.set_ylabel('Penetrance (%)')
        ax3.set_title('C. Model vs Literature Penetrance Estimates')
        ax3.set_xticks(x)
        ax3.set_xticklabels(mutations)
        ax3.legend()
        ax3.set_yscale('log')
        ax3.set_ylim(0.1, 100)
        
        # Add value labels
        for bar, val in zip(bars3, model_penetrance):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        for bar, val in zip(bars4, lit_penetrance):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                    f'{val:.0f}%', ha='center', va='bottom', fontweight='bold')
        
        # Panel D: Risk Communication Chart
        # Show risk in different formats for patient communication
        overall_penetrance = self.summary['overall_penetrance_percent'].iloc[0]
        
        risk_formats = ['Percentage', 'Odds (1 in X)', 'Natural\nFrequency']
        risk_values = [
            overall_penetrance,
            100 / overall_penetrance,
            f"{overall_penetrance:.1f} out of 100"
        ]
        
        # Create text-based visualization
        ax4.axis('off')
        
        # Create risk communication table
        risk_data = [
            ['Risk Format', 'Value', 'Interpretation'],
            ['Percentage', f'{overall_penetrance:.2f}%', 'Low risk'],
            ['Odds', f'1 in {100/overall_penetrance:.0f}', 'Most carriers unaffected'],
            ['Natural Frequency', f'{overall_penetrance:.1f} out of 100 carriers', 'Population perspective'],
            ['Absolute Risk', f'{self.summary["total_patient_prevalence_per_100k"].iloc[0]:.1f} per 100,000', 'Very rare condition']
        ]
        
        # Create table
        table = ax4.table(cellText=risk_data[1:], colLabels=risk_data[0],
                         cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1, 2)
        
        # Style the table
        for i in range(len(risk_data[0])):
            table[(0, i)].set_facecolor('#4CAF50')
            table[(0, i)].set_text_props(weight='bold', color='white')
        
        ax4.set_title('D. Risk Communication Formats', pad=20, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/real_prevalence_population_impact.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created: real_prevalence_population_impact.png")
    
    def create_all_visualizations(self):
        """Create all real prevalence visualizations"""
        
        print("LHON REAL PREVALENCE VISUALIZATION GENERATION")
        print("=" * 50)
        
        # Create all visualizations
        self.create_carrier_vs_patient_overview()
        self.create_penetrance_risk_stratification()
        self.create_population_impact_analysis()
        
        print("\nAll real prevalence visualizations created successfully!")
        print("\nGenerated files:")
        print("- real_prevalence_overview.png")
        print("- real_prevalence_risk_stratification.png")
        print("- real_prevalence_population_impact.png")
        
        return {
            'visualizations_created': 3,
            'status': 'success'
        }

def main():
    """Generate all real prevalence visualizations"""
    
    visualizer = LHONRealPrevalenceVisualizer()
    results = visualizer.create_all_visualizations()
    
    return results

if __name__ == "__main__":
    results = main()

