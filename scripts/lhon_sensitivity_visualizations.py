#!/usr/bin/env python3
"""
LHON Sensitivity Analysis Visualizations
Creates comprehensive visualizations of sensitivity analysis results
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

class LHONSensitivityVisualizer:
    """
    Creates visualizations for LHON sensitivity analysis results
    """
    
    def __init__(self):
        """Initialize and load sensitivity analysis data"""
        
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
            'high_sensitivity': '#d62728',    # Red
            'medium_sensitivity': '#ff7f0e',  # Orange
            'low_sensitivity': '#2ca02c',     # Green
            'correlation_pos': '#1f77b4',     # Blue
            'correlation_neg': '#e377c2',     # Pink
            'scenario_base': '#9467bd',       # Purple
            'scenario_opt': '#2ca02c',        # Green
            'scenario_pess': '#d62728'        # Red
        }
    
    def load_data(self):
        """Load sensitivity analysis results"""
        
        try:
            self.sensitivity_indices = pd.read_csv('/home/ubuntu/lhon_sensitivity_indices.csv', index_col=0)
            self.monte_carlo_correlations = pd.read_csv('/home/ubuntu/lhon_monte_carlo_correlations.csv', index_col=0)
            self.scenario_analysis = pd.read_csv('/home/ubuntu/lhon_scenario_analysis.csv')
            
            print("Successfully loaded sensitivity analysis data")
            
        except FileNotFoundError as e:
            print(f"Error loading data: {e}")
            raise
    
    def create_sensitivity_indices_plot(self):
        """Create plot showing sensitivity indices for all parameters"""
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
        
        # Panel A: Sensitivity Index Ranking
        sensitivity_data = self.sensitivity_indices.sort_values('sensitivity_index', ascending=True)
        
        # Color code by sensitivity level
        colors = []
        for idx in sensitivity_data['sensitivity_index']:
            if idx > 0.5:
                colors.append(self.colors['high_sensitivity'])
            elif idx > 0.2:
                colors.append(self.colors['medium_sensitivity'])
            else:
                colors.append(self.colors['low_sensitivity'])
        
        bars1 = ax1.barh(range(len(sensitivity_data)), sensitivity_data['sensitivity_index'], 
                        color=colors, alpha=0.8)
        
        ax1.set_yticks(range(len(sensitivity_data)))
        ax1.set_yticklabels([self.format_parameter_name(name) for name in sensitivity_data.index])
        ax1.set_xlabel('Sensitivity Index')
        ax1.set_title('A. Parameter Sensitivity Ranking')
        ax1.grid(axis='x', alpha=0.3)
        
        # Add value labels
        for i, (bar, val) in enumerate(zip(bars1, sensitivity_data['sensitivity_index'])):
            width = bar.get_width()
            ax1.text(width + 0.01, bar.get_y() + bar.get_height()/2,
                    f'{val:.3f}', ha='left', va='center', fontweight='bold')
        
        # Panel B: Coefficient of Variation
        cv_data = sensitivity_data.sort_values('coefficient_of_variation', ascending=True)
        
        bars2 = ax2.barh(range(len(cv_data)), cv_data['coefficient_of_variation'], 
                        color=self.colors['medium_sensitivity'], alpha=0.8)
        
        ax2.set_yticks(range(len(cv_data)))
        ax2.set_yticklabels([self.format_parameter_name(name) for name in cv_data.index])
        ax2.set_xlabel('Coefficient of Variation')
        ax2.set_title('B. Parameter Variability Impact')
        ax2.grid(axis='x', alpha=0.3)
        
        # Panel C: Prevalence Range by Parameter
        prevalence_ranges = sensitivity_data['max_prevalence'] - sensitivity_data['min_prevalence']
        prevalence_ranges = prevalence_ranges.sort_values(ascending=True)
        
        bars3 = ax3.barh(range(len(prevalence_ranges)), prevalence_ranges, 
                        color=self.colors['correlation_pos'], alpha=0.8)
        
        ax3.set_yticks(range(len(prevalence_ranges)))
        ax3.set_yticklabels([self.format_parameter_name(name) for name in prevalence_ranges.index])
        ax3.set_xlabel('Prevalence Range (per 100,000)')
        ax3.set_title('C. Prevalence Range by Parameter')
        ax3.grid(axis='x', alpha=0.3)
        
        # Panel D: Sensitivity vs Baseline Impact
        ax4.scatter(sensitivity_data['sensitivity_index'], sensitivity_data['coefficient_of_variation'],
                   c=[colors[i] for i in range(len(sensitivity_data))], s=100, alpha=0.7)
        
        # Add parameter labels
        for i, (idx, row) in enumerate(sensitivity_data.iterrows()):
            ax4.annotate(self.format_parameter_name(idx), 
                        (row['sensitivity_index'], row['coefficient_of_variation']),
                        xytext=(5, 5), textcoords='offset points', fontsize=10)
        
        ax4.set_xlabel('Sensitivity Index')
        ax4.set_ylabel('Coefficient of Variation')
        ax4.set_title('D. Sensitivity vs Variability')
        ax4.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/sensitivity_indices_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created: sensitivity_indices_analysis.png")
    
    def create_monte_carlo_correlations_plot(self):
        """Create plot showing Monte Carlo correlation analysis"""
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
        
        # Panel A: Correlations with Prevalence
        prevalence_corrs = self.monte_carlo_correlations['prevalence'].sort_values(ascending=True)
        
        colors = ['red' if x < 0 else 'blue' for x in prevalence_corrs]
        
        bars1 = ax1.barh(range(len(prevalence_corrs)), prevalence_corrs, color=colors, alpha=0.7)
        
        ax1.set_yticks(range(len(prevalence_corrs)))
        ax1.set_yticklabels([self.format_parameter_name(name) for name in prevalence_corrs.index])
        ax1.set_xlabel('Correlation with Prevalence')
        ax1.set_title('A. Parameter Correlations with Population Prevalence')
        ax1.axvline(x=0, color='black', linestyle='-', alpha=0.5)
        ax1.grid(axis='x', alpha=0.3)
        
        # Add value labels
        for bar, val in zip(bars1, prevalence_corrs):
            width = bar.get_width()
            label_x = width + 0.02 if width >= 0 else width - 0.02
            ha = 'left' if width >= 0 else 'right'
            ax1.text(label_x, bar.get_y() + bar.get_height()/2,
                    f'{val:.3f}', ha=ha, va='center', fontweight='bold')
        
        # Panel B: Correlations with Male 11778G>A Penetrance
        male_11778_corrs = self.monte_carlo_correlations['penetrance_male_11778'].sort_values(ascending=True)
        
        colors = ['red' if x < 0 else 'green' for x in male_11778_corrs]
        
        bars2 = ax2.barh(range(len(male_11778_corrs)), male_11778_corrs, color=colors, alpha=0.7)
        
        ax2.set_yticks(range(len(male_11778_corrs)))
        ax2.set_yticklabels([self.format_parameter_name(name) for name in male_11778_corrs.index])
        ax2.set_xlabel('Correlation with Male 11778G>A Penetrance')
        ax2.set_title('B. Parameter Correlations with Male 11778G>A Penetrance')
        ax2.axvline(x=0, color='black', linestyle='-', alpha=0.5)
        ax2.grid(axis='x', alpha=0.3)
        
        # Panel C: Correlation Heatmap
        corr_matrix = self.monte_carlo_correlations.T
        
        # Create custom colormap
        im = ax3.imshow(corr_matrix, cmap='RdBu_r', aspect='auto', vmin=-1, vmax=1)
        
        # Set ticks and labels
        ax3.set_xticks(range(len(corr_matrix.columns)))
        ax3.set_yticks(range(len(corr_matrix.index)))
        ax3.set_xticklabels([self.format_parameter_name(name) for name in corr_matrix.columns], 
                           rotation=45, ha='right')
        ax3.set_yticklabels(['Prevalence', 'Male 11778G>A', 'Female 11778G>A', 'Male 14484T>C', 'Male 3460G>A'])
        
        # Add correlation values
        for i in range(len(corr_matrix.index)):
            for j in range(len(corr_matrix.columns)):
                text = ax3.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}',
                               ha="center", va="center", 
                               color="white" if abs(corr_matrix.iloc[i, j]) > 0.5 else "black",
                               fontsize=8, fontweight='bold')
        
        ax3.set_title('C. Correlation Matrix: Parameters vs Outcomes')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax3, fraction=0.046, pad=0.04)
        cbar.set_label('Correlation Coefficient')
        
        # Panel D: Absolute Correlation Comparison
        abs_corr_prevalence = np.abs(prevalence_corrs)
        abs_corr_male_11778 = np.abs(male_11778_corrs)
        
        x = np.arange(len(abs_corr_prevalence))
        width = 0.35
        
        bars3 = ax4.bar(x - width/2, abs_corr_prevalence, width, 
                       label='Prevalence', color=self.colors['correlation_pos'], alpha=0.8)
        bars4 = ax4.bar(x + width/2, abs_corr_male_11778, width, 
                       label='Male 11778G>A', color=self.colors['correlation_neg'], alpha=0.8)
        
        ax4.set_xlabel('Parameters')
        ax4.set_ylabel('Absolute Correlation')
        ax4.set_title('D. Absolute Correlation Comparison')
        ax4.set_xticks(x)
        ax4.set_xticklabels([self.format_parameter_name(name) for name in abs_corr_prevalence.index], 
                           rotation=45, ha='right')
        ax4.legend()
        ax4.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/monte_carlo_correlations.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created: monte_carlo_correlations.png")
    
    def create_scenario_analysis_plot(self):
        """Create plot showing scenario analysis results"""
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
        
        # Panel A: Prevalence by Scenario
        scenarios = self.scenario_analysis['scenario']
        prevalences = self.scenario_analysis['prevalence']
        
        scenario_colors = [self.colors['scenario_opt'], self.colors['scenario_base'], self.colors['scenario_pess']]
        
        bars1 = ax1.bar(scenarios, prevalences, color=scenario_colors, alpha=0.8)
        
        ax1.set_ylabel('Population Prevalence (per 100,000)')
        ax1.set_title('A. Population Prevalence by Scenario')
        ax1.set_ylim(0, max(prevalences) * 1.1)
        
        # Add value labels
        for bar, val in zip(bars1, prevalences):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{val:.1f}', ha='center', va='bottom', fontweight='bold')
        
        # Add reference line for literature range
        ax1.axhline(y=1.87, color='red', linestyle='--', alpha=0.7, 
                   label='Literature Average (1.87)')
        ax1.legend()
        
        # Panel B: Penetrance Comparison Across Scenarios
        penetrance_columns = ['male_11778_nonsmoker', 'female_11778_nonsmoker', 
                             'male_11778_heavy_smoker', 'male_14484_j_haplotype', 'male_3460_nonsmoker']
        
        x = np.arange(len(penetrance_columns))
        width = 0.25
        
        for i, scenario in enumerate(scenarios):
            scenario_data = self.scenario_analysis[self.scenario_analysis['scenario'] == scenario]
            penetrances = [scenario_data[col].iloc[0] for col in penetrance_columns]
            
            bars = ax2.bar(x + i*width, penetrances, width, 
                          label=scenario.title(), color=scenario_colors[i], alpha=0.8)
        
        ax2.set_xlabel('Risk Scenarios')
        ax2.set_ylabel('Penetrance (%)')
        ax2.set_title('B. Penetrance by Risk Scenario and Model Scenario')
        ax2.set_xticks(x + width)
        ax2.set_xticklabels(['Male 11778\nNon-smoker', 'Female 11778\nNon-smoker', 
                            'Male 11778\nHeavy smoker', 'Male 14484\nJ haplotype', 'Male 3460\nNon-smoker'],
                           rotation=45, ha='right')
        ax2.legend()
        ax2.set_yscale('log')
        ax2.set_ylim(0.1, 100)
        
        # Panel C: Scenario Range Analysis
        penetrance_ranges = {}
        for col in penetrance_columns:
            values = self.scenario_analysis[col]
            penetrance_ranges[col] = max(values) - min(values)
        
        range_names = list(penetrance_ranges.keys())
        range_values = list(penetrance_ranges.values())
        
        bars3 = ax3.bar(range_names, range_values, color=self.colors['medium_sensitivity'], alpha=0.8)
        
        ax3.set_ylabel('Penetrance Range (%)')
        ax3.set_title('C. Penetrance Uncertainty by Risk Scenario')
        ax3.set_xticklabels(['Male 11778\nNon-smoker', 'Female 11778\nNon-smoker', 
                            'Male 11778\nHeavy smoker', 'Male 14484\nJ haplotype', 'Male 3460\nNon-smoker'],
                           rotation=45, ha='right')
        
        # Add value labels
        for bar, val in zip(bars3, range_values):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        # Panel D: Relative Risk Ratios
        base_case = self.scenario_analysis[self.scenario_analysis['scenario'] == 'base_case']
        
        risk_ratios = {}
        for scenario in ['optimistic', 'pessimistic']:
            scenario_data = self.scenario_analysis[self.scenario_analysis['scenario'] == scenario]
            
            ratio = scenario_data['prevalence'].iloc[0] / base_case['prevalence'].iloc[0]
            risk_ratios[scenario] = ratio
        
        scenario_names = list(risk_ratios.keys())
        ratio_values = list(risk_ratios.values())
        
        colors_ratio = [self.colors['scenario_opt'], self.colors['scenario_pess']]
        
        bars4 = ax4.bar(scenario_names, ratio_values, color=colors_ratio, alpha=0.8)
        
        ax4.set_ylabel('Prevalence Ratio (vs Base Case)')
        ax4.set_title('D. Scenario Risk Ratios')
        ax4.axhline(y=1, color='black', linestyle='-', alpha=0.5, label='Base Case')
        ax4.legend()
        
        # Add value labels
        for bar, val in zip(bars4, ratio_values):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                    f'{val:.2f}x', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/scenario_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created: scenario_analysis.png")
    
    def create_parameter_impact_summary(self):
        """Create comprehensive parameter impact summary"""
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
        
        # Panel A: Top 5 Most Sensitive Parameters
        top_5_params = self.sensitivity_indices.nlargest(5, 'sensitivity_index')
        
        bars1 = ax1.bar(range(len(top_5_params)), top_5_params['sensitivity_index'], 
                       color=self.colors['high_sensitivity'], alpha=0.8)
        
        ax1.set_xticks(range(len(top_5_params)))
        ax1.set_xticklabels([self.format_parameter_name(name) for name in top_5_params.index], 
                           rotation=45, ha='right')
        ax1.set_ylabel('Sensitivity Index')
        ax1.set_title('A. Top 5 Most Sensitive Parameters')
        
        # Add value labels
        for bar, val in zip(bars1, top_5_params['sensitivity_index']):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                    f'{val:.3f}', ha='center', va='bottom', fontweight='bold')
        
        # Panel B: Parameter Categories
        param_categories = {
            'Genetic': ['base_liability_11778', 'base_liability_14484', 'base_liability_3460', 'haplogroup_j_14484_effect'],
            'Environmental': ['smoking_heavy_effect', 'smoking_light_effect', 'alcohol_heavy_effect'],
            'Demographic': ['male_effect'],
            'Model': ['liability_threshold', 'liability_sigma']
        }
        
        category_sensitivities = {}
        for category, params in param_categories.items():
            category_sens = []
            for param in params:
                if param in self.sensitivity_indices.index:
                    category_sens.append(self.sensitivity_indices.loc[param, 'sensitivity_index'])
            category_sensitivities[category] = np.mean(category_sens) if category_sens else 0
        
        categories = list(category_sensitivities.keys())
        sensitivities = list(category_sensitivities.values())
        
        bars2 = ax2.bar(categories, sensitivities, 
                       color=[self.colors['high_sensitivity'], self.colors['medium_sensitivity'], 
                             self.colors['low_sensitivity'], self.colors['correlation_pos']], alpha=0.8)
        
        ax2.set_ylabel('Average Sensitivity Index')
        ax2.set_title('B. Sensitivity by Parameter Category')
        
        # Add value labels
        for bar, val in zip(bars2, sensitivities):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                    f'{val:.3f}', ha='center', va='bottom', fontweight='bold')
        
        # Panel C: Correlation Strength Distribution
        all_correlations = []
        for outcome in self.monte_carlo_correlations.columns:
            all_correlations.extend(np.abs(self.monte_carlo_correlations[outcome].values))
        
        ax3.hist(all_correlations, bins=20, color=self.colors['correlation_pos'], alpha=0.7, edgecolor='black')
        ax3.axvline(np.mean(all_correlations), color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {np.mean(all_correlations):.3f}')
        ax3.set_xlabel('Absolute Correlation Coefficient')
        ax3.set_ylabel('Frequency')
        ax3.set_title('C. Distribution of Parameter-Outcome Correlations')
        ax3.legend()
        
        # Panel D: Summary Statistics Table
        ax4.axis('off')
        
        summary_stats = [
            ['Metric', 'Value'],
            ['Most Sensitive Parameter', self.format_parameter_name(self.sensitivity_indices.idxmax()['sensitivity_index'])],
            ['Highest Sensitivity Index', f"{self.sensitivity_indices['sensitivity_index'].max():.3f}"],
            ['Average Sensitivity Index', f"{self.sensitivity_indices['sensitivity_index'].mean():.3f}"],
            ['Strongest Correlation', f"{np.abs(self.monte_carlo_correlations.values).max():.3f}"],
            ['Prevalence Range (Optimistic)', f"{self.scenario_analysis[self.scenario_analysis['scenario']=='optimistic']['prevalence'].iloc[0]:.1f} per 100k"],
            ['Prevalence Range (Pessimistic)', f"{self.scenario_analysis[self.scenario_analysis['scenario']=='pessimistic']['prevalence'].iloc[0]:.1f} per 100k"],
            ['Uncertainty Factor', f"{self.scenario_analysis[self.scenario_analysis['scenario']=='pessimistic']['prevalence'].iloc[0] / self.scenario_analysis[self.scenario_analysis['scenario']=='optimistic']['prevalence'].iloc[0]:.1f}x"]
        ]
        
        # Create table
        table = ax4.table(cellText=summary_stats[1:], colLabels=summary_stats[0],
                         cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        table.scale(1, 2)
        
        # Style the table
        for i in range(len(summary_stats[0])):
            table[(0, i)].set_facecolor('#4CAF50')
            table[(0, i)].set_text_props(weight='bold', color='white')
        
        ax4.set_title('D. Sensitivity Analysis Summary Statistics', pad=20, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('/home/ubuntu/parameter_impact_summary.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created: parameter_impact_summary.png")
    
    def format_parameter_name(self, param_name):
        """Format parameter names for display"""
        
        name_mapping = {
            'male_effect': 'Male Sex',
            'smoking_heavy_effect': 'Heavy Smoking',
            'smoking_light_effect': 'Light Smoking',
            'alcohol_heavy_effect': 'Heavy Alcohol',
            'haplogroup_j_14484_effect': 'Haplogroup J',
            'liability_threshold': 'Liability Threshold',
            'liability_sigma': 'Liability Sigma',
            'base_liability_11778': '11778G>A Base',
            'base_liability_14484': '14484T>C Base',
            'base_liability_3460': '3460G>A Base'
        }
        
        return name_mapping.get(param_name, param_name)
    
    def create_all_visualizations(self):
        """Create all sensitivity analysis visualizations"""
        
        print("LHON SENSITIVITY ANALYSIS VISUALIZATION GENERATION")
        print("=" * 50)
        
        # Create all visualizations
        self.create_sensitivity_indices_plot()
        self.create_monte_carlo_correlations_plot()
        self.create_scenario_analysis_plot()
        self.create_parameter_impact_summary()
        
        print("\nAll sensitivity analysis visualizations created successfully!")
        print("\nGenerated files:")
        print("- sensitivity_indices_analysis.png")
        print("- monte_carlo_correlations.png")
        print("- scenario_analysis.png")
        print("- parameter_impact_summary.png")
        
        return {
            'visualizations_created': 4,
            'status': 'success'
        }

def main():
    """Generate all sensitivity analysis visualizations"""
    
    visualizer = LHONSensitivityVisualizer()
    results = visualizer.create_all_visualizations()
    
    return results

if __name__ == "__main__":
    results = main()

