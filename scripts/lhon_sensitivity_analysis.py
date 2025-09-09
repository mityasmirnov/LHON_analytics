#!/usr/bin/env python3
"""
LHON Sensitivity Analysis Framework
Conducts comprehensive sensitivity analysis of key genetic and environmental factors
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from itertools import product
import warnings
warnings.filterwarnings('ignore')

class LHONSensitivityAnalyzer:
    """
    Comprehensive sensitivity analysis for LHON models
    """
    
    def __init__(self):
        """Initialize sensitivity analysis framework"""
        
        # Base parameters from our models
        self.base_parameters = {
            # Liability threshold model parameters
            'liability_threshold': 2.0,
            'liability_sigma': 1.0,
            
            # Base liability scores by mutation
            'base_liability': {
                '11778G>A': 0.5,
                '14484T>C': -0.5,
                '3460G>A': 1.8
            },
            
            # Effect sizes (log odds ratios)
            'male_effect': np.log(7.11),
            'smoking_heavy_effect': np.log(3.16),
            'smoking_light_effect': np.log(1.54),
            'alcohol_heavy_effect': np.log(3.27),
            'haplogroup_j_14484_effect': np.log(2.0),  # Estimated
            'age_peak_effect': 0.5,
            
            # Population frequencies
            'carrier_frequencies': {
                '11778G>A': 42.54 / 100000,
                '14484T>C': 65.58 / 100000,
                '3460G>A': 1.77 / 100000
            },
            
            'male_proportion': 0.5,
            'smoking_heavy_proportion': 0.25,
            'smoking_light_proportion': 0.35,
            'alcohol_heavy_proportion': 0.2,
            'haplogroup_j_proportion': 0.08,
            'age_peak_proportion': 0.4
        }
        
        # Define parameter ranges for sensitivity analysis
        self.parameter_ranges = {
            'male_effect': (np.log(3), np.log(15)),  # OR 3-15
            'smoking_heavy_effect': (np.log(1.5), np.log(6)),  # OR 1.5-6
            'smoking_light_effect': (np.log(1.0), np.log(3)),  # OR 1-3
            'alcohol_heavy_effect': (np.log(1.5), np.log(6)),  # OR 1.5-6
            'haplogroup_j_14484_effect': (np.log(1.0), np.log(4)),  # OR 1-4
            'liability_threshold': (1.0, 4.0),
            'liability_sigma': (0.5, 2.0),
            'base_liability_11778': (-1.0, 2.0),
            'base_liability_14484': (-2.0, 1.0),
            'base_liability_3460': (0.5, 3.0)
        }
        
        # Results storage
        self.sensitivity_results = {}
        
    def calculate_penetrance(self, mutation, sex='Male', smoking='None', alcohol='None', 
                           haplogroup='Other', age='Peak', parameters=None):
        """Calculate penetrance for given scenario and parameters"""
        
        if parameters is None:
            parameters = self.base_parameters
        
        # Base liability for mutation
        if mutation == '11778G>A':
            liability = parameters.get('base_liability_11778', parameters['base_liability']['11778G>A'])
        elif mutation == '14484T>C':
            liability = parameters.get('base_liability_14484', parameters['base_liability']['14484T>C'])
        elif mutation == '3460G>A':
            liability = parameters.get('base_liability_3460', parameters['base_liability']['3460G>A'])
        else:
            liability = 0
        
        # Add effects
        if sex == 'Male':
            liability += parameters['male_effect']
        
        if smoking == 'Heavy':
            liability += parameters['smoking_heavy_effect']
        elif smoking == 'Light':
            liability += parameters['smoking_light_effect']
        
        if alcohol == 'Heavy':
            liability += parameters['alcohol_heavy_effect']
        
        if haplogroup == 'J' and mutation == '14484T>C':
            liability += parameters['haplogroup_j_14484_effect']
        
        if age == 'Peak':
            liability += parameters.get('age_peak_effect', 0)
        
        # Calculate penetrance
        threshold = parameters['liability_threshold']
        sigma = parameters['liability_sigma']
        
        penetrance = stats.norm.cdf((liability - threshold) / sigma)
        
        return max(0, min(1, penetrance))
    
    def calculate_population_prevalence(self, parameters=None):
        """Calculate overall population prevalence"""
        
        if parameters is None:
            parameters = self.base_parameters
        
        total_prevalence = 0
        
        # Iterate through all combinations
        mutations = ['11778G>A', '14484T>C', '3460G>A']
        sexes = ['Male', 'Female']
        smoking_levels = ['None', 'Light', 'Heavy']
        alcohol_levels = ['None', 'Heavy']
        haplogroups = ['Other', 'J']
        ages = ['Young', 'Peak', 'Middle']
        
        for mutation in mutations:
            carrier_freq = parameters['carrier_frequencies'][mutation]
            
            for sex, smoking, alcohol, haplogroup, age in product(sexes, smoking_levels, alcohol_levels, haplogroups, ages):
                
                # Calculate subgroup frequency
                subgroup_freq = carrier_freq
                subgroup_freq *= parameters['male_proportion'] if sex == 'Male' else (1 - parameters['male_proportion'])
                
                if smoking == 'Heavy':
                    subgroup_freq *= parameters['smoking_heavy_proportion']
                elif smoking == 'Light':
                    subgroup_freq *= parameters['smoking_light_proportion']
                else:
                    subgroup_freq *= (1 - parameters['smoking_heavy_proportion'] - parameters['smoking_light_proportion'])
                
                subgroup_freq *= parameters['alcohol_heavy_proportion'] if alcohol == 'Heavy' else (1 - parameters['alcohol_heavy_proportion'])
                subgroup_freq *= parameters['haplogroup_j_proportion'] if haplogroup == 'J' else (1 - parameters['haplogroup_j_proportion'])
                subgroup_freq *= parameters['age_peak_proportion'] if age == 'Peak' else ((1 - parameters['age_peak_proportion']) / 2)
                
                # Calculate penetrance for this subgroup
                penetrance = self.calculate_penetrance(mutation, sex, smoking, alcohol, haplogroup, age, parameters)
                
                # Add to total prevalence
                total_prevalence += subgroup_freq * penetrance
        
        return total_prevalence * 100000  # Convert to per 100,000
    
    def one_at_a_time_sensitivity(self, n_points=20):
        """Perform one-at-a-time sensitivity analysis"""
        
        print("Performing one-at-a-time sensitivity analysis...")
        
        results = {}
        
        for param_name, (min_val, max_val) in self.parameter_ranges.items():
            
            param_values = np.linspace(min_val, max_val, n_points)
            prevalences = []
            penetrances_male = []
            penetrances_female = []
            
            for param_val in param_values:
                # Create modified parameters
                modified_params = self.base_parameters.copy()
                modified_params[param_name] = param_val
                
                # Calculate outcomes
                prevalence = self.calculate_population_prevalence(modified_params)
                pen_male = self.calculate_penetrance('11778G>A', 'Male', parameters=modified_params)
                pen_female = self.calculate_penetrance('11778G>A', 'Female', parameters=modified_params)
                
                prevalences.append(prevalence)
                penetrances_male.append(pen_male * 100)
                penetrances_female.append(pen_female * 100)
            
            results[param_name] = {
                'values': param_values,
                'prevalences': prevalences,
                'penetrances_male': penetrances_male,
                'penetrances_female': penetrances_female
            }
            
            print(f"Completed sensitivity analysis for {param_name}")
        
        self.sensitivity_results['one_at_a_time'] = results
        return results
    
    def calculate_sensitivity_indices(self):
        """Calculate sensitivity indices for each parameter"""
        
        print("Calculating sensitivity indices...")
        
        indices = {}
        
        # Base case
        base_prevalence = self.calculate_population_prevalence()
        
        for param_name, param_data in self.sensitivity_results['one_at_a_time'].items():
            
            prevalences = np.array(param_data['prevalences'])
            
            # Calculate range-based sensitivity index
            prevalence_range = np.max(prevalences) - np.min(prevalences)
            sensitivity_index = prevalence_range / base_prevalence
            
            # Calculate normalized sensitivity (coefficient of variation)
            cv = np.std(prevalences) / np.mean(prevalences)
            
            indices[param_name] = {
                'sensitivity_index': sensitivity_index,
                'coefficient_of_variation': cv,
                'min_prevalence': np.min(prevalences),
                'max_prevalence': np.max(prevalences),
                'base_prevalence': base_prevalence
            }
        
        self.sensitivity_results['indices'] = indices
        return indices
    
    def monte_carlo_sensitivity(self, n_samples=1000):
        """Perform Monte Carlo sensitivity analysis"""
        
        print("Performing Monte Carlo sensitivity analysis...")
        
        # Generate random parameter samples
        parameter_samples = {}
        
        for param_name, (min_val, max_val) in self.parameter_ranges.items():
            parameter_samples[param_name] = np.random.uniform(min_val, max_val, n_samples)
        
        # Calculate outcomes for each sample
        prevalences = []
        penetrances_male_11778 = []
        penetrances_female_11778 = []
        penetrances_male_14484 = []
        penetrances_male_3460 = []
        
        for i in range(n_samples):
            
            # Create parameter set for this sample
            sample_params = self.base_parameters.copy()
            for param_name in self.parameter_ranges.keys():
                sample_params[param_name] = parameter_samples[param_name][i]
            
            # Calculate outcomes
            prevalence = self.calculate_population_prevalence(sample_params)
            pen_male_11778 = self.calculate_penetrance('11778G>A', 'Male', parameters=sample_params)
            pen_female_11778 = self.calculate_penetrance('11778G>A', 'Female', parameters=sample_params)
            pen_male_14484 = self.calculate_penetrance('14484T>C', 'Male', parameters=sample_params)
            pen_male_3460 = self.calculate_penetrance('3460G>A', 'Male', parameters=sample_params)
            
            prevalences.append(prevalence)
            penetrances_male_11778.append(pen_male_11778 * 100)
            penetrances_female_11778.append(pen_female_11778 * 100)
            penetrances_male_14484.append(pen_male_14484 * 100)
            penetrances_male_3460.append(pen_male_3460 * 100)
        
        # Calculate correlations
        correlations = {}
        outcomes = {
            'prevalence': prevalences,
            'penetrance_male_11778': penetrances_male_11778,
            'penetrance_female_11778': penetrances_female_11778,
            'penetrance_male_14484': penetrances_male_14484,
            'penetrance_male_3460': penetrances_male_3460
        }
        
        for param_name in self.parameter_ranges.keys():
            correlations[param_name] = {}
            for outcome_name, outcome_values in outcomes.items():
                corr = np.corrcoef(parameter_samples[param_name], outcome_values)[0, 1]
                correlations[param_name][outcome_name] = corr
        
        self.sensitivity_results['monte_carlo'] = {
            'parameter_samples': parameter_samples,
            'outcomes': outcomes,
            'correlations': correlations
        }
        
        return correlations
    
    def scenario_analysis(self):
        """Perform scenario analysis for extreme cases"""
        
        print("Performing scenario analysis...")
        
        scenarios = {
            'optimistic': {
                'male_effect': np.log(3),  # Lower male effect
                'smoking_heavy_effect': np.log(1.5),  # Lower smoking effect
                'liability_threshold': 3.0,  # Higher threshold
                'base_liability_11778': -0.5,  # Lower base liability
                'base_liability_14484': -1.5,
                'base_liability_3460': 0.5
            },
            'pessimistic': {
                'male_effect': np.log(15),  # Higher male effect
                'smoking_heavy_effect': np.log(6),  # Higher smoking effect
                'liability_threshold': 1.0,  # Lower threshold
                'base_liability_11778': 1.5,  # Higher base liability
                'base_liability_14484': 0.5,
                'base_liability_3460': 2.5
            },
            'base_case': self.base_parameters
        }
        
        scenario_results = {}
        
        for scenario_name, scenario_params in scenarios.items():
            
            if scenario_name == 'base_case':
                params = scenario_params
            else:
                params = self.base_parameters.copy()
                params.update(scenario_params)
            
            # Calculate key outcomes
            prevalence = self.calculate_population_prevalence(params)
            
            # Calculate penetrances for key scenarios
            penetrances = {
                'male_11778_nonsmoker': self.calculate_penetrance('11778G>A', 'Male', parameters=params),
                'female_11778_nonsmoker': self.calculate_penetrance('11778G>A', 'Female', parameters=params),
                'male_11778_heavy_smoker': self.calculate_penetrance('11778G>A', 'Male', 'Heavy', parameters=params),
                'male_14484_j_haplotype': self.calculate_penetrance('14484T>C', 'Male', haplogroup='J', parameters=params),
                'male_3460_nonsmoker': self.calculate_penetrance('3460G>A', 'Male', parameters=params)
            }
            
            scenario_results[scenario_name] = {
                'prevalence': prevalence,
                'penetrances': {k: v * 100 for k, v in penetrances.items()}
            }
        
        self.sensitivity_results['scenarios'] = scenario_results
        return scenario_results
    
    def run_complete_analysis(self):
        """Run complete sensitivity analysis"""
        
        print("LHON SENSITIVITY ANALYSIS")
        print("=" * 30)
        
        # Run all analyses
        oat_results = self.one_at_a_time_sensitivity()
        indices = self.calculate_sensitivity_indices()
        mc_correlations = self.monte_carlo_sensitivity()
        scenarios = self.scenario_analysis()
        
        # Save results
        self.save_results()
        
        # Print summary
        self.print_summary()
        
        return {
            'one_at_a_time': oat_results,
            'sensitivity_indices': indices,
            'monte_carlo_correlations': mc_correlations,
            'scenarios': scenarios
        }
    
    def save_results(self):
        """Save sensitivity analysis results"""
        
        # Save sensitivity indices
        indices_df = pd.DataFrame(self.sensitivity_results['indices']).T
        indices_df.to_csv('/home/ubuntu/lhon_sensitivity_indices.csv')
        
        # Save Monte Carlo correlations
        correlations_df = pd.DataFrame(self.sensitivity_results['monte_carlo']['correlations']).T
        correlations_df.to_csv('/home/ubuntu/lhon_monte_carlo_correlations.csv')
        
        # Save scenario results
        scenario_data = []
        for scenario_name, scenario_data_dict in self.sensitivity_results['scenarios'].items():
            row = {'scenario': scenario_name, 'prevalence': scenario_data_dict['prevalence']}
            row.update(scenario_data_dict['penetrances'])
            scenario_data.append(row)
        
        scenarios_df = pd.DataFrame(scenario_data)
        scenarios_df.to_csv('/home/ubuntu/lhon_scenario_analysis.csv', index=False)
        
        print("Saved sensitivity analysis results to CSV files")
    
    def print_summary(self):
        """Print summary of sensitivity analysis"""
        
        print("\nSENSITIVITY ANALYSIS SUMMARY")
        print("-" * 30)
        
        # Most sensitive parameters
        indices = self.sensitivity_results['indices']
        sorted_params = sorted(indices.items(), key=lambda x: x[1]['sensitivity_index'], reverse=True)
        
        print("Most sensitive parameters (by sensitivity index):")
        for i, (param_name, data) in enumerate(sorted_params[:5]):
            print(f"{i+1}. {param_name}: {data['sensitivity_index']:.3f}")
        
        # Scenario analysis summary
        scenarios = self.sensitivity_results['scenarios']
        print(f"\nScenario Analysis:")
        print(f"Base case prevalence: {scenarios['base_case']['prevalence']:.2f} per 100,000")
        print(f"Optimistic prevalence: {scenarios['optimistic']['prevalence']:.2f} per 100,000")
        print(f"Pessimistic prevalence: {scenarios['pessimistic']['prevalence']:.2f} per 100,000")
        
        # Monte Carlo correlations
        correlations = self.sensitivity_results['monte_carlo']['correlations']
        print(f"\nStrongest correlations with prevalence:")
        
        prevalence_corrs = [(param, data['prevalence']) for param, data in correlations.items()]
        prevalence_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
        
        for param, corr in prevalence_corrs[:5]:
            print(f"{param}: {corr:.3f}")

def main():
    """Run sensitivity analysis"""
    
    analyzer = LHONSensitivityAnalyzer()
    results = analyzer.run_complete_analysis()
    
    return results

if __name__ == "__main__":
    results = main()

