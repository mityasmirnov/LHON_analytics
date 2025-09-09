#!/usr/bin/env python3
"""
Comprehensive Mathematical Models for LHON Penetrance and Prevalence
Based on literature review and gnomAD population data
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

# Set random seed for reproducibility
np.random.seed(42)

class LHONPenetranceModels:
    """
    Comprehensive mathematical models for LHON penetrance and prevalence
    """
    
    def __init__(self):
        """Initialize with literature-derived parameters"""
        
        # gnomAD-derived carrier frequencies (per 100,000)
        self.carrier_frequencies = {
            '11778G>A': 42.54,  # 1 in 2,351
            '14484T>C': 65.58,  # 1 in 1,524  
            '3460G>A': 1.77     # 1 in 56,426
        }
        
        # Population prevalence (cases per 100,000)
        self.population_prevalence = {
            'Madrid_2024': 0.79,
            'UK_historical': 3.23,
            'Finland': 2.0,
            'Australia': 1.46
        }
        
        # Environmental risk factors (odds ratios from Kirkman et al.)
        self.environmental_ors = {
            'smoking_heavy': 3.16,
            'smoking_light': 1.54,
            'alcohol_heavy': 3.27,
            'alcohol_light': 1.01,
            'male_sex': 7.11,
            'heteroplasmy_protective': 0.37
        }
        
        # Haplogroup effects
        self.haplogroup_ors = {
            '11778G>A': {'J': 1.31, 'H': 0.79},
            '14484T>C': {'J': 27.0, 'non_J': 0.037},  # Derived from gnomAD analysis
            '3460G>A': {'K': 2.0}  # Literature estimate
        }
        
        # Age-related parameters
        self.age_params = {
            'peak_onset_age': 25,
            'age_std': 10,
            'late_onset_proportion': 0.10  # 10% after age 50
        }
        
    def liability_threshold_model(self, mutation, sex, haplogroup=None, 
                                environmental_factors=None, age=25):
        """
        Liability threshold model for LHON penetrance
        P(affected) = Φ((β₀ + β₁X₁ + ... + βₙXₙ - T)/σ)
        """
        
        # Base liability scores by mutation (from literature)
        base_liability = {
            '11778G>A': 0.5,   # 3.8% baseline penetrance
            '14484T>C': 0.2,   # 0.8% baseline penetrance
            '3460G>A': 1.8     # 14.1% baseline penetrance
        }
        
        liability = base_liability.get(mutation, 0.5)
        
        # Sex effect (males have higher liability)
        if sex == 'male':
            liability += np.log(self.environmental_ors['male_sex'])
        
        # Haplogroup effects
        if haplogroup and mutation in self.haplogroup_ors:
            if haplogroup in self.haplogroup_ors[mutation]:
                or_value = self.haplogroup_ors[mutation][haplogroup]
                liability += np.log(or_value)
        
        # Environmental factors
        if environmental_factors:
            for factor, present in environmental_factors.items():
                if present and factor in self.environmental_ors:
                    liability += np.log(self.environmental_ors[factor])
        
        # Age effect (younger onset has slightly higher liability)
        age_factor = np.exp(-(age - self.age_params['peak_onset_age'])**2 / 
                           (2 * self.age_params['age_std']**2))
        liability += 0.2 * age_factor
        
        # Convert to penetrance using cumulative normal distribution
        threshold = 2.0  # Liability threshold
        sigma = 1.0      # Residual variance
        
        penetrance = stats.norm.cdf((liability - threshold) / sigma)
        
        return penetrance, liability
    
    def bayesian_hierarchical_model(self, n_simulations=10000):
        """
        Bayesian hierarchical model with uncertainty quantification
        """
        
        results = []
        
        for _ in range(n_simulations):
            # Sample from prior distributions
            base_penetrance = {
                '11778G>A': np.random.beta(4, 96),      # ~4% mean
                '14484T>C': np.random.beta(1, 124),     # ~0.8% mean  
                '3460G>A': np.random.beta(14, 86)       # ~14% mean
            }
            
            # Sample environmental effects
            male_or = np.random.lognormal(np.log(7.11), 0.2)
            smoking_or = np.random.lognormal(np.log(3.16), 0.3)
            
            # Sample haplogroup effects
            hap_j_11778 = np.random.lognormal(np.log(1.31), 0.1)
            hap_j_14484 = np.random.lognormal(np.log(27.0), 0.5)
            
            # Calculate penetrance for different scenarios
            scenarios = [
                ('11778G>A', 'female', None, {}),
                ('11778G>A', 'male', None, {}),
                ('11778G>A', 'male', 'J', {'smoking_heavy': True}),
                ('14484T>C', 'male', 'J', {}),
                ('14484T>C', 'male', 'non_J', {}),
                ('3460G>A', 'male', None, {})
            ]
            
            sim_results = {}
            for mutation, sex, hap, env in scenarios:
                base = base_penetrance[mutation]
                
                # Apply modifiers
                if sex == 'male':
                    base *= male_or
                if env.get('smoking_heavy'):
                    base *= smoking_or
                if mutation == '11778G>A' and hap == 'J':
                    base *= hap_j_11778
                elif mutation == '14484T>C' and hap == 'J':
                    base *= hap_j_14484
                elif mutation == '14484T>C' and hap == 'non_J':
                    base *= 0.037  # Very low penetrance on non-J
                
                # Cap at 100%
                penetrance = min(base, 1.0)
                
                scenario_name = f"{mutation}_{sex}_{hap}_{list(env.keys())}"
                sim_results[scenario_name] = penetrance
            
            results.append(sim_results)
        
        return pd.DataFrame(results)
    
    def monte_carlo_population_model(self, population_size=100000, n_simulations=1000):
        """
        Monte Carlo simulation of LHON in a population
        """
        
        results = []
        
        for sim in range(n_simulations):
            # Generate population
            population = []
            
            for i in range(population_size):
                # Assign mutation (based on gnomAD frequencies)
                rand = np.random.random() * 100000
                
                if rand < self.carrier_frequencies['3460G>A']:
                    mutation = '3460G>A'
                elif rand < self.carrier_frequencies['3460G>A'] + self.carrier_frequencies['11778G>A']:
                    mutation = '11778G>A'
                elif rand < sum(self.carrier_frequencies.values()):
                    mutation = '14484T>C'
                else:
                    mutation = None  # No LHON mutation
                
                if mutation:
                    # Assign demographics
                    sex = 'male' if np.random.random() < 0.5 else 'female'
                    age = np.random.normal(25, 10)
                    age = max(15, min(80, age))  # Constrain age range
                    
                    # Assign haplogroup (simplified)
                    if mutation == '14484T>C':
                        haplogroup = 'J' if np.random.random() < 0.08 else 'non_J'  # 8% J in gnomAD
                    elif mutation == '11778G>A':
                        haplogroup = 'J' if np.random.random() < 0.15 else 'other'
                    else:
                        haplogroup = 'L2'  # 3460G>A only on L2 in gnomAD
                    
                    # Assign environmental factors
                    smoking = np.random.random() < 0.6  # 60% smoking rate
                    heavy_smoking = smoking and np.random.random() < 0.3
                    alcohol = np.random.random() < 0.9   # 90% drink alcohol
                    heavy_alcohol = alcohol and np.random.random() < 0.2
                    
                    env_factors = {
                        'smoking_heavy': heavy_smoking,
                        'smoking_light': smoking and not heavy_smoking,
                        'alcohol_heavy': heavy_alcohol,
                        'alcohol_light': alcohol and not heavy_alcohol
                    }
                    
                    # Calculate penetrance
                    penetrance, liability = self.liability_threshold_model(
                        mutation, sex, haplogroup, env_factors, age
                    )
                    
                    # Determine if affected
                    affected = np.random.random() < penetrance
                    
                    population.append({
                        'id': i,
                        'mutation': mutation,
                        'sex': sex,
                        'age': age,
                        'haplogroup': haplogroup,
                        'smoking_heavy': heavy_smoking,
                        'alcohol_heavy': heavy_alcohol,
                        'penetrance': penetrance,
                        'liability': liability,
                        'affected': affected
                    })
            
            # Calculate population statistics
            df = pd.DataFrame(population)
            
            if len(df) > 0:
                total_carriers = len(df)
                total_affected = df['affected'].sum()
                
                # By mutation
                mutation_stats = df.groupby('mutation').agg({
                    'affected': ['count', 'sum', 'mean']
                }).round(4)
                
                # By sex
                sex_stats = df.groupby('sex').agg({
                    'affected': ['count', 'sum', 'mean']
                }).round(4)
                
                # Overall penetrance
                overall_penetrance = total_affected / total_carriers if total_carriers > 0 else 0
                
                # Population prevalence (per 100,000)
                pop_prevalence = (total_affected / population_size) * 100000
                
                results.append({
                    'simulation': sim,
                    'total_carriers': total_carriers,
                    'total_affected': total_affected,
                    'overall_penetrance': overall_penetrance,
                    'population_prevalence': pop_prevalence,
                    'carrier_frequency': (total_carriers / population_size) * 100000
                })
        
        return pd.DataFrame(results), df if len(results) > 0 else None
    
    def calculate_revised_penetrance_estimates(self):
        """
        Calculate revised penetrance estimates based on gnomAD data
        """
        
        # Known population prevalence (average from studies)
        avg_prevalence = np.mean(list(self.population_prevalence.values()))  # per 100,000
        
        results = {}
        
        for mutation in self.carrier_frequencies:
            carrier_freq = self.carrier_frequencies[mutation]  # per 100,000
            
            # Assume mutation accounts for its proportion of cases
            mutation_proportions = {
                '11778G>A': 0.65,
                '14484T>C': 0.20,
                '3460G>A': 0.10
            }
            
            mutation_prevalence = avg_prevalence * mutation_proportions[mutation]
            
            # Calculate penetrance
            penetrance = (mutation_prevalence / carrier_freq) * 100 if carrier_freq > 0 else 0
            
            results[mutation] = {
                'carrier_frequency_per_100k': carrier_freq,
                'estimated_prevalence_per_100k': mutation_prevalence,
                'calculated_penetrance_percent': penetrance,
                'literature_penetrance_percent': {
                    '11778G>A': 43,
                    '14484T>C': 65,
                    '3460G>A': 78
                }[mutation],
                'penetrance_ratio_literature_vs_calculated': {
                    '11778G>A': 43,
                    '14484T>C': 65,
                    '3460G>A': 78
                }[mutation] / penetrance if penetrance > 0 else np.inf
            }
        
        return results

def main():
    """Run comprehensive LHON modeling analysis"""
    
    print("LHON Mathematical Modeling Analysis")
    print("=" * 50)
    
    # Initialize models
    models = LHONPenetranceModels()
    
    # 1. Liability Threshold Model Examples
    print("\n1. LIABILITY THRESHOLD MODEL")
    print("-" * 30)
    
    scenarios = [
        ('11778G>A', 'female', None, {}, 25),
        ('11778G>A', 'male', None, {}, 25),
        ('11778G>A', 'male', 'J', {'smoking_heavy': True}, 25),
        ('14484T>C', 'male', 'J', {}, 25),
        ('14484T>C', 'male', 'non_J', {}, 25),
        ('3460G>A', 'male', None, {}, 25),
        ('11778G>A', 'male', None, {'heteroplasmy_protective': True}, 25)
    ]
    
    liability_results = []
    
    for mutation, sex, hap, env, age in scenarios:
        penetrance, liability = models.liability_threshold_model(
            mutation, sex, hap, env, age
        )
        
        scenario_desc = f"{mutation} {sex}"
        if hap:
            scenario_desc += f" {hap}"
        if env:
            scenario_desc += f" {list(env.keys())}"
        
        liability_results.append({
            'Scenario': scenario_desc,
            'Liability': round(liability, 3),
            'Penetrance': f"{penetrance*100:.1f}%"
        })
        
        print(f"{scenario_desc:<40} Liability: {liability:.3f} Penetrance: {penetrance*100:.1f}%")
    
    # 2. Revised Penetrance Estimates
    print("\n2. REVISED PENETRANCE ESTIMATES (gnomAD-informed)")
    print("-" * 50)
    
    revised_estimates = models.calculate_revised_penetrance_estimates()
    
    for mutation, data in revised_estimates.items():
        print(f"\n{mutation}:")
        print(f"  Carrier frequency: 1 in {100000/data['carrier_frequency_per_100k']:.0f}")
        print(f"  Calculated penetrance: {data['calculated_penetrance_percent']:.2f}%")
        print(f"  Literature penetrance: {data['literature_penetrance_percent']}%")
        print(f"  Literature overestimate: {data['penetrance_ratio_literature_vs_calculated']:.1f}x")
    
    # 3. Bayesian Hierarchical Model
    print("\n3. BAYESIAN HIERARCHICAL MODEL (1000 simulations)")
    print("-" * 50)
    
    bayesian_results = models.bayesian_hierarchical_model(1000)
    
    # Summary statistics
    for col in bayesian_results.columns:
        mean_val = bayesian_results[col].mean()
        ci_low = bayesian_results[col].quantile(0.025)
        ci_high = bayesian_results[col].quantile(0.975)
        print(f"{col:<40} Mean: {mean_val*100:.1f}% (95% CI: {ci_low*100:.1f}-{ci_high*100:.1f}%)")
    
    # 4. Monte Carlo Population Model
    print("\n4. MONTE CARLO POPULATION MODEL")
    print("-" * 40)
    
    mc_results, sample_pop = models.monte_carlo_population_model(
        population_size=100000, n_simulations=100
    )
    
    if mc_results is not None and len(mc_results) > 0:
        print(f"Average carrier frequency: {mc_results['carrier_frequency'].mean():.1f} per 100,000")
        print(f"Average population prevalence: {mc_results['population_prevalence'].mean():.2f} per 100,000")
        print(f"Average overall penetrance: {mc_results['overall_penetrance'].mean()*100:.2f}%")
        
        # Save results
        liability_df = pd.DataFrame(liability_results)
        
        # Create summary tables
        print("\n5. SAVING RESULTS")
        print("-" * 20)
        
        # Save liability threshold results
        liability_df.to_csv('/home/ubuntu/lhon_liability_model_results.csv', index=False)
        print("Saved: lhon_liability_model_results.csv")
        
        # Save Bayesian results
        bayesian_results.to_csv('/home/ubuntu/lhon_bayesian_model_results.csv', index=False)
        print("Saved: lhon_bayesian_model_results.csv")
        
        # Save Monte Carlo results
        mc_results.to_csv('/home/ubuntu/lhon_monte_carlo_results.csv', index=False)
        print("Saved: lhon_monte_carlo_results.csv")
        
        # Save revised estimates
        revised_df = pd.DataFrame(revised_estimates).T
        revised_df.to_csv('/home/ubuntu/lhon_revised_penetrance_estimates.csv')
        print("Saved: lhon_revised_penetrance_estimates.csv")
        
        print("\nModeling analysis complete!")
        
        return {
            'liability_results': liability_df,
            'bayesian_results': bayesian_results,
            'monte_carlo_results': mc_results,
            'revised_estimates': revised_df
        }
    
    else:
        print("Error: Monte Carlo simulation failed")
        return None

if __name__ == "__main__":
    results = main()

