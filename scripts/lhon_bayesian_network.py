#!/usr/bin/env python3
"""
Bayesian Network and Advanced Monte Carlo Models for LHON
Implements probabilistic graphical models for penetrance and prevalence
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import networkx as nx
from itertools import product
import warnings
warnings.filterwarnings('ignore')

class LHONBayesianNetwork:
    """
    Bayesian Network for LHON penetrance modeling
    """
    
    def __init__(self):
        """Initialize Bayesian Network structure and parameters"""
        
        # Network structure (parent -> children relationships)
        self.network_structure = {
            'Population': [],  # Root node
            'mtDNA_Mutation': ['Population'],
            'Haplogroup': ['Population'],
            'Sex': ['Population'],
            'Age': ['Population'],
            'Smoking': ['Population', 'Sex'],
            'Alcohol': ['Population', 'Sex'],
            'Nuclear_Variants': ['Population'],
            'Environmental_Stress': ['Population'],
            'Mitochondrial_Function': ['mtDNA_Mutation', 'Haplogroup', 'Nuclear_Variants'],
            'Oxidative_Stress': ['Smoking', 'Alcohol', 'Environmental_Stress'],
            'Liability': ['Mitochondrial_Function', 'Oxidative_Stress', 'Sex', 'Age'],
            'LHON_Phenotype': ['Liability'],
            'Recovery': ['LHON_Phenotype', 'Age', 'mtDNA_Mutation']
        }
        
        # Prior probabilities and conditional probability tables
        self.initialize_parameters()
        
    def initialize_parameters(self):
        """Initialize all conditional probability tables"""
        
        # Population base rates
        self.priors = {
            'Population': {'General': 1.0},
            
            # mtDNA mutations (per 100,000 from gnomAD)
            'mtDNA_Mutation': {
                'None': 0.9988,
                '11778G>A': 0.0004254,
                '14484T>C': 0.0006558,
                '3460G>A': 0.0000177
            },
            
            # Haplogroups (European population)
            'Haplogroup': {
                'H': 0.45,
                'J': 0.08,
                'K': 0.06,
                'T': 0.09,
                'U': 0.18,
                'Other': 0.14
            },
            
            # Sex distribution
            'Sex': {
                'Male': 0.5,
                'Female': 0.5
            },
            
            # Age distribution (normal around 25)
            'Age': {
                'Young': 0.4,      # <20
                'Peak': 0.4,       # 20-35
                'Middle': 0.15,    # 35-50
                'Late': 0.05       # >50
            }
        }
        
        # Conditional probabilities
        self.conditional_probs = {
            
            # Smoking depends on sex
            'Smoking': {
                ('Male',): {'None': 0.35, 'Light': 0.35, 'Heavy': 0.30},
                ('Female',): {'None': 0.45, 'Light': 0.35, 'Heavy': 0.20}
            },
            
            # Alcohol depends on sex
            'Alcohol': {
                ('Male',): {'None': 0.05, 'Light': 0.70, 'Heavy': 0.25},
                ('Female',): {'None': 0.10, 'Light': 0.75, 'Heavy': 0.15}
            },
            
            # Nuclear variants (rare)
            'Nuclear_Variants': {
                (): {'None': 0.995, 'DNAJC30': 0.003, 'Other_CI': 0.002}
            },
            
            # Environmental stress
            'Environmental_Stress': {
                (): {'Low': 0.6, 'Moderate': 0.3, 'High': 0.1}
            },
            
            # Mitochondrial function depends on mutation, haplogroup, nuclear variants
            'Mitochondrial_Function': self._create_mito_function_cpt(),
            
            # Oxidative stress depends on smoking, alcohol, environmental stress
            'Oxidative_Stress': self._create_oxidative_stress_cpt(),
            
            # Liability threshold depends on multiple factors
            'Liability': self._create_liability_cpt(),
            
            # LHON phenotype depends on liability
            'LHON_Phenotype': {
                ('Very_Low',): {'Unaffected': 0.999, 'Affected': 0.001},
                ('Low',): {'Unaffected': 0.98, 'Affected': 0.02},
                ('Moderate',): {'Unaffected': 0.90, 'Affected': 0.10},
                ('High',): {'Unaffected': 0.50, 'Affected': 0.50},
                ('Very_High',): {'Unaffected': 0.05, 'Affected': 0.95}
            },
            
            # Recovery depends on phenotype, age, mutation
            'Recovery': self._create_recovery_cpt()
        }
    
    def _create_mito_function_cpt(self):
        """Create conditional probability table for mitochondrial function"""
        
        cpt = {}
        
        mutations = ['None', '11778G>A', '14484T>C', '3460G>A']
        haplogroups = ['H', 'J', 'K', 'T', 'U', 'Other']
        nuclear_vars = ['None', 'DNAJC30', 'Other_CI']
        
        for mut, hap, nuc in product(mutations, haplogroups, nuclear_vars):
            
            # Base function level
            if mut == 'None':
                if nuc == 'None':
                    probs = {'Normal': 0.95, 'Mild_Impair': 0.04, 'Severe_Impair': 0.01}
                else:  # Nuclear variants cause impairment
                    probs = {'Normal': 0.1, 'Mild_Impair': 0.4, 'Severe_Impair': 0.5}
            else:
                # mtDNA mutations cause impairment
                base_impair = {
                    '11778G>A': 0.7,
                    '14484T>C': 0.5,
                    '3460G>A': 0.9
                }[mut]
                
                # Haplogroup modifies impairment
                if mut == '14484T>C' and hap == 'J':
                    base_impair *= 1.5  # J haplogroup worsens 14484T>C
                elif mut == '11778G>A' and hap == 'J':
                    base_impair *= 1.1  # Mild worsening for 11778G>A
                elif hap in ['H', 'U']:
                    base_impair *= 0.9  # Slight protection
                
                # Nuclear variants worsen impairment
                if nuc != 'None':
                    base_impair *= 1.3
                
                # Convert to probabilities
                base_impair = min(base_impair, 0.95)
                probs = {
                    'Normal': 1 - base_impair,
                    'Mild_Impair': base_impair * 0.4,
                    'Severe_Impair': base_impair * 0.6
                }
            
            cpt[(mut, hap, nuc)] = probs
        
        return cpt
    
    def _create_oxidative_stress_cpt(self):
        """Create conditional probability table for oxidative stress"""
        
        cpt = {}
        
        smoking_levels = ['None', 'Light', 'Heavy']
        alcohol_levels = ['None', 'Light', 'Heavy']
        env_stress = ['Low', 'Moderate', 'High']
        
        for smoke, alc, env in product(smoking_levels, alcohol_levels, env_stress):
            
            # Base stress level
            stress_score = 0
            
            # Smoking contribution
            stress_score += {'None': 0, 'Light': 1, 'Heavy': 3}[smoke]
            
            # Alcohol contribution
            stress_score += {'None': 0, 'Light': 0.5, 'Heavy': 2}[alc]
            
            # Environmental contribution
            stress_score += {'Low': 0, 'Moderate': 1, 'High': 2}[env]
            
            # Convert to categorical
            if stress_score <= 1:
                probs = {'Low': 0.8, 'Moderate': 0.15, 'High': 0.05}
            elif stress_score <= 3:
                probs = {'Low': 0.3, 'Moderate': 0.5, 'High': 0.2}
            else:
                probs = {'Low': 0.1, 'Moderate': 0.3, 'High': 0.6}
            
            cpt[(smoke, alc, env)] = probs
        
        return cpt
    
    def _create_liability_cpt(self):
        """Create conditional probability table for liability"""
        
        cpt = {}
        
        mito_func = ['Normal', 'Mild_Impair', 'Severe_Impair']
        ox_stress = ['Low', 'Moderate', 'High']
        sex = ['Male', 'Female']
        age = ['Young', 'Peak', 'Middle', 'Late']
        
        for mf, os, s, a in product(mito_func, ox_stress, sex, age):
            
            # Base liability
            liability_score = 0
            
            # Mitochondrial function (strongest predictor)
            liability_score += {'Normal': 0, 'Mild_Impair': 2, 'Severe_Impair': 4}[mf]
            
            # Oxidative stress
            liability_score += {'Low': 0, 'Moderate': 1, 'High': 2}[os]
            
            # Sex (males higher risk)
            liability_score += {'Male': 1, 'Female': 0}[s]
            
            # Age (peak age has highest risk)
            liability_score += {'Young': 0.5, 'Peak': 1, 'Middle': 0.5, 'Late': 0.2}[a]
            
            # Convert to liability categories
            if liability_score <= 1:
                probs = {'Very_Low': 0.7, 'Low': 0.25, 'Moderate': 0.04, 'High': 0.01, 'Very_High': 0.0}
            elif liability_score <= 2:
                probs = {'Very_Low': 0.4, 'Low': 0.4, 'Moderate': 0.15, 'High': 0.04, 'Very_High': 0.01}
            elif liability_score <= 4:
                probs = {'Very_Low': 0.1, 'Low': 0.3, 'Moderate': 0.4, 'High': 0.15, 'Very_High': 0.05}
            elif liability_score <= 6:
                probs = {'Very_Low': 0.02, 'Low': 0.08, 'Moderate': 0.3, 'High': 0.4, 'Very_High': 0.2}
            else:
                probs = {'Very_Low': 0.0, 'Low': 0.02, 'Moderate': 0.08, 'High': 0.3, 'Very_High': 0.6}
            
            cpt[(mf, os, s, a)] = probs
        
        return cpt
    
    def _create_recovery_cpt(self):
        """Create conditional probability table for recovery"""
        
        cpt = {}
        
        phenotypes = ['Unaffected', 'Affected']
        ages = ['Young', 'Peak', 'Middle', 'Late']
        mutations = ['None', '11778G>A', '14484T>C', '3460G>A']
        
        for pheno, age, mut in product(phenotypes, ages, mutations):
            
            if pheno == 'Unaffected':
                probs = {'No_Recovery': 1.0, 'Partial': 0.0, 'Complete': 0.0}
            else:  # Affected
                # Base recovery rates by mutation
                recovery_rates = {
                    'None': 0.0,
                    '11778G>A': 0.04,    # 4% recovery
                    '14484T>C': 0.37,    # 37% recovery
                    '3460G>A': 0.20      # 20% recovery
                }
                
                base_recovery = recovery_rates.get(mut, 0.0)
                
                # Age effect (younger patients recover better)
                age_multiplier = {'Young': 2.0, 'Peak': 1.0, 'Middle': 0.5, 'Late': 0.2}[age]
                
                total_recovery = min(base_recovery * age_multiplier, 0.8)
                
                # Split between partial and complete recovery
                complete_recovery = total_recovery * 0.3
                partial_recovery = total_recovery * 0.7
                no_recovery = 1 - total_recovery
                
                probs = {
                    'No_Recovery': no_recovery,
                    'Partial': partial_recovery,
                    'Complete': complete_recovery
                }
            
            cpt[(pheno, age, mut)] = probs
        
        return cpt
    
    def sample_from_network(self, n_samples=10000):
        """Sample from the Bayesian network"""
        
        samples = []
        
        for _ in range(n_samples):
            sample = {}
            
            # Sample in topological order
            for node in ['Population', 'mtDNA_Mutation', 'Haplogroup', 'Sex', 'Age',
                        'Smoking', 'Alcohol', 'Nuclear_Variants', 'Environmental_Stress',
                        'Mitochondrial_Function', 'Oxidative_Stress', 'Liability',
                        'LHON_Phenotype', 'Recovery']:
                
                sample[node] = self._sample_node(node, sample)
            
            samples.append(sample)
        
        return pd.DataFrame(samples)
    
    def _sample_node(self, node, current_sample):
        """Sample a single node given its parents"""
        
        if node in self.priors:
            # Root node - sample from prior
            probs = self.priors[node]
        else:
            # Get parent values
            parents = self.network_structure[node]
            parent_values = tuple(current_sample[p] for p in parents)
            
            # Get conditional probabilities
            if node in self.conditional_probs:
                if parent_values in self.conditional_probs[node]:
                    probs = self.conditional_probs[node][parent_values]
                else:
                    # Default probabilities if combination not found
                    if node == 'Smoking':
                        probs = {'None': 0.4, 'Light': 0.35, 'Heavy': 0.25}
                    elif node == 'Alcohol':
                        probs = {'None': 0.075, 'Light': 0.725, 'Heavy': 0.2}
                    else:
                        # Uniform distribution as fallback
                        states = list(next(iter(self.conditional_probs[node].values())).keys())
                        prob_val = 1.0 / len(states)
                        probs = {state: prob_val for state in states}
            else:
                raise ValueError(f"No probabilities defined for node {node}")
        
        # Sample from categorical distribution
        states = list(probs.keys())
        probabilities = list(probs.values())
        
        # Normalize probabilities to ensure they sum to 1
        prob_sum = sum(probabilities)
        if prob_sum > 0:
            probabilities = [p / prob_sum for p in probabilities]
        else:
            probabilities = [1.0 / len(states)] * len(states)
        
        return np.random.choice(states, p=probabilities)
    
    def calculate_penetrance_by_subgroup(self, samples_df):
        """Calculate penetrance for different subgroups"""
        
        results = {}
        
        # Overall penetrance
        total_carriers = samples_df[samples_df['mtDNA_Mutation'] != 'None']
        if len(total_carriers) > 0:
            overall_penetrance = (total_carriers['LHON_Phenotype'] == 'Affected').mean()
            results['Overall'] = overall_penetrance
        
        # By mutation
        for mutation in ['11778G>A', '14484T>C', '3460G>A']:
            mut_carriers = samples_df[samples_df['mtDNA_Mutation'] == mutation]
            if len(mut_carriers) > 0:
                penetrance = (mut_carriers['LHON_Phenotype'] == 'Affected').mean()
                results[f'{mutation}'] = penetrance
        
        # By sex
        for sex in ['Male', 'Female']:
            sex_carriers = total_carriers[total_carriers['Sex'] == sex]
            if len(sex_carriers) > 0:
                penetrance = (sex_carriers['LHON_Phenotype'] == 'Affected').mean()
                results[f'{sex}'] = penetrance
        
        # By mutation and sex
        for mutation in ['11778G>A', '14484T>C', '3460G>A']:
            for sex in ['Male', 'Female']:
                subgroup = samples_df[
                    (samples_df['mtDNA_Mutation'] == mutation) & 
                    (samples_df['Sex'] == sex)
                ]
                if len(subgroup) > 0:
                    penetrance = (subgroup['LHON_Phenotype'] == 'Affected').mean()
                    results[f'{mutation}_{sex}'] = penetrance
        
        # By environmental factors
        heavy_smokers = total_carriers[total_carriers['Smoking'] == 'Heavy']
        if len(heavy_smokers) > 0:
            penetrance = (heavy_smokers['LHON_Phenotype'] == 'Affected').mean()
            results['Heavy_Smokers'] = penetrance
        
        # High-risk combination (male, heavy smoker, 11778G>A)
        high_risk = samples_df[
            (samples_df['mtDNA_Mutation'] == '11778G>A') & 
            (samples_df['Sex'] == 'Male') &
            (samples_df['Smoking'] == 'Heavy')
        ]
        if len(high_risk) > 0:
            penetrance = (high_risk['LHON_Phenotype'] == 'Affected').mean()
            results['High_Risk_11778_Male_Heavy_Smoker'] = penetrance
        
        return results
    
    def calculate_recovery_rates(self, samples_df):
        """Calculate recovery rates by subgroup"""
        
        affected = samples_df[samples_df['LHON_Phenotype'] == 'Affected']
        
        if len(affected) == 0:
            return {}
        
        results = {}
        
        # Overall recovery
        recovery_rate = (affected['Recovery'] != 'No_Recovery').mean()
        complete_recovery = (affected['Recovery'] == 'Complete').mean()
        results['Overall_Any_Recovery'] = recovery_rate
        results['Overall_Complete_Recovery'] = complete_recovery
        
        # By mutation
        for mutation in ['11778G>A', '14484T>C', '3460G>A']:
            mut_affected = affected[affected['mtDNA_Mutation'] == mutation]
            if len(mut_affected) > 0:
                recovery_rate = (mut_affected['Recovery'] != 'No_Recovery').mean()
                complete_recovery = (mut_affected['Recovery'] == 'Complete').mean()
                results[f'{mutation}_Any_Recovery'] = recovery_rate
                results[f'{mutation}_Complete_Recovery'] = complete_recovery
        
        # By age
        for age in ['Young', 'Peak', 'Middle', 'Late']:
            age_affected = affected[affected['Age'] == age]
            if len(age_affected) > 0:
                recovery_rate = (age_affected['Recovery'] != 'No_Recovery').mean()
                results[f'{age}_Recovery'] = recovery_rate
        
        return results

def run_bayesian_analysis():
    """Run comprehensive Bayesian network analysis"""
    
    print("LHON Bayesian Network Analysis")
    print("=" * 40)
    
    # Initialize network
    bn = LHONBayesianNetwork()
    
    # Generate samples
    print("Generating samples from Bayesian network...")
    samples = bn.sample_from_network(50000)
    
    # Calculate basic statistics
    print("\nBasic Population Statistics:")
    print("-" * 30)
    
    # Carrier frequencies
    mutation_counts = samples['mtDNA_Mutation'].value_counts()
    total_samples = len(samples)
    
    for mutation in ['11778G>A', '14484T>C', '3460G>A']:
        if mutation in mutation_counts:
            frequency = (mutation_counts[mutation] / total_samples) * 100000
            print(f"{mutation}: {frequency:.1f} per 100,000 (1 in {100000/frequency:.0f})")
    
    # Penetrance analysis
    print("\nPenetrance Analysis:")
    print("-" * 20)
    
    penetrance_results = bn.calculate_penetrance_by_subgroup(samples)
    
    for subgroup, penetrance in penetrance_results.items():
        print(f"{subgroup:<35} {penetrance*100:.1f}%")
    
    # Recovery analysis
    print("\nRecovery Analysis:")
    print("-" * 18)
    
    recovery_results = bn.calculate_recovery_rates(samples)
    
    for subgroup, rate in recovery_results.items():
        print(f"{subgroup:<35} {rate*100:.1f}%")
    
    # Population prevalence
    print("\nPopulation Prevalence:")
    print("-" * 22)
    
    affected_count = (samples['LHON_Phenotype'] == 'Affected').sum()
    prevalence_per_100k = (affected_count / total_samples) * 100000
    
    print(f"Modeled prevalence: {prevalence_per_100k:.2f} per 100,000")
    print(f"Literature range: 0.79-3.23 per 100,000")
    
    # Save results
    print("\nSaving Results:")
    print("-" * 15)
    
    # Save samples
    samples.to_csv('/home/ubuntu/lhon_bayesian_samples.csv', index=False)
    print("Saved: lhon_bayesian_samples.csv")
    
    # Save penetrance results
    penetrance_df = pd.DataFrame(list(penetrance_results.items()), 
                                columns=['Subgroup', 'Penetrance'])
    penetrance_df.to_csv('/home/ubuntu/lhon_bayesian_penetrance.csv', index=False)
    print("Saved: lhon_bayesian_penetrance.csv")
    
    # Save recovery results
    recovery_df = pd.DataFrame(list(recovery_results.items()), 
                              columns=['Subgroup', 'Recovery_Rate'])
    recovery_df.to_csv('/home/ubuntu/lhon_bayesian_recovery.csv', index=False)
    print("Saved: lhon_bayesian_recovery.csv")
    
    print("\nBayesian network analysis complete!")
    
    return {
        'samples': samples,
        'penetrance_results': penetrance_results,
        'recovery_results': recovery_results,
        'prevalence_per_100k': prevalence_per_100k
    }

if __name__ == "__main__":
    results = run_bayesian_analysis()

