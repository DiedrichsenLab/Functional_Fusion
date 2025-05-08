"""
This module generates synthetic data to simulate the full variance decomposition
for a given set of datasets and subjects. It creates a covariance matrix with
specific variance components and then decomposes the variance using the 
`reliability.decompose_variance_from_SS_scaled` function.
It also compares the estimated parameters with the true parameters using Mean Absolute Percentage Error (MAPE).

@author: Ali Shahbazi
"""

import numpy as np
import pandas as pd
import Functional_Fusion.reliability as rel


def generate_vectors(subjects_per_dataset, N_part, N_SS):
    """
    Generate dataset, subject, and partition vectors.
    """
    dataset_vec = np.empty(N_SS, dtype=object)
    sub_vec = np.empty(N_SS, dtype=object)
    part_vec = np.empty(N_SS, dtype=object)
    subjects = []

    current_idx = 0
    for d, n_d in enumerate(subjects_per_dataset):
        for s in range(n_d):
            sub_id = f'sub-{s + 1:02d}'
            for p in range(N_part):
                dataset_vec[current_idx] = f'dataset_{d+1}'
                sub_vec[current_idx] = sub_id
                part_vec[current_idx] = f'partition_{p+1}'
                if p == 0:
                    subjects.append((f'dataset_{d+1}', sub_id))
                current_idx += 1

    return dataset_vec, sub_vec, part_vec, subjects


def generate_true_parameters(N_datasets, N_subj):
    """
    Generate true parameters for the simulation.
    """
    MIN_SC = 0.01
    MAX_SC = 0.1
    MIN_V_D = 0.5
    MAX_V_D = 5.0
    MIN_V_S = 1.0
    MAX_V_S = 8.0
    MIN_V_M = 4.0
    MAX_V_M = 14.0

    # Generate true parameters
    true_sc = np.random.uniform(MIN_SC, MAX_SC, N_subj)
    true_V_u = 1.0
    true_V_d = np.random.uniform(MIN_V_D, MAX_V_D, N_datasets)
    true_V_s = np.random.uniform(MIN_V_S, MAX_V_S, N_datasets)
    true_V_m = np.random.uniform(MIN_V_M, MAX_V_M, N_subj)

    return true_sc, true_V_u, true_V_d, true_V_s, true_V_m


def compare_resutls(true_df, Q_df):
    params = ['sc', 'v_u', 'v_d', 'v_s', 'v_m']
    # Compare Q_df with true_df
    comparison_results = {}
    for column in params:
        true_values = true_df[column].values
        estimated_values = Q_df[column].values
        # Calculate Mean Absolute Percentage Error (MAPE)
        mape = np.mean(np.abs(estimated_values - true_values) / true_values) * 100
        comparison_results[column] = {'MAPE': mape}

    # Print comparison results
    print('----------------')
    print(f"Mean Absolute Percentage Error:")
    for col, metrics in comparison_results.items():
        print(f"{col}: {metrics['MAPE']:.2f}%")


def generate_synthetic_data(subjects_per_dataset, N_part, data_size=1000):
    """
    Generate synthetic data to make covariance matrix with specific variance components.
    Args:
        subjects_per_dataset (list): List of integers representing the number of subjects in each dataset.
        N_part (int): Number of partitions.
        generete (str): Type of data to generate ('data' or 'cov').
        data_size (int): Size of the data to generate.
    Returns:
        covariance_matrix (np.ndarray): Covariance matrix.
        dataset_vec (np.ndarray): Vector of dataset names.
        sub_vec (np.ndarray): Vector of subject IDs.
        part_vec (np.ndarray): Vector of partition IDs.
        true_Q_df (pd.DataFrame): DataFrame containing true parameters.
    """

    N_datasets = len(subjects_per_dataset)
    N_subj = sum(subjects_per_dataset)
    N_SS = N_subj * N_part
    
    # Fill the vectors with dataset, subject, and partition information
    dataset_vec, sub_vec, part_vec, subjects = generate_vectors(subjects_per_dataset, N_part, N_SS)
    
    unique_datasets = list(dict.fromkeys(dataset_vec))

    # Generate true parameters
    true_sc, true_V_u, true_V_d, true_V_s, true_V_m = generate_true_parameters(N_datasets, N_subj)
    
    # Generate data and covariance matrix
    W_u = np.random.normal(0, np.sqrt(true_V_u), data_size)
    W_d = np.random.normal(0, np.sqrt(true_V_d), (data_size, N_datasets))
    W_hat = np.zeros((data_size, N_SS))
    
    current_idx = 0
    for d, n_d in enumerate(subjects_per_dataset):
        for _ in range(n_d):
            # Generate W_s for each subject in the dataset
            W_s = np.random.normal(0, np.sqrt(true_V_s[d]), data_size)
            for p in range(N_part):
                # Generate W_m for each partition
                W_m = np.random.normal(0, np.sqrt(true_V_m[current_idx]), data_size)

                # Generate W_hat for each subject and partition
                W_hat[:, current_idx * N_part + p] = true_sc[current_idx] * (W_u + W_d[:,d] + W_s + W_m)
            current_idx += 1
    
    covariance_matrix = np.cov(W_hat, rowvar=False)
    
    # True Q DataFrame
    true_Q_df = pd.DataFrame({
        'train_dataset': [sid[0] for sid in subjects],
        'subj_id': [sid[1] for sid in subjects],
        'sc': true_sc,
        'v_u': true_V_u,
        'v_d': true_V_d[[unique_datasets.index(ds) for ds, _ in subjects]],
        'v_s': true_V_s[[unique_datasets.index(ds) for ds, _ in subjects]],
        'v_m': true_V_m
    })
    
    return covariance_matrix, dataset_vec, sub_vec, part_vec, true_Q_df


if __name__ == "__main__":

    # ----- Parameters -----
    subjects_per_dataset = [15, 6, 26, 12, 35]
    N_part = 2
    
    # ----- Generate synthetic data -----
    covariance_matrix, dataset_vec, sub_vec, part_vec, true_df = generate_synthetic_data(
        subjects_per_dataset, N_part, data_size=10000)
    
    # ----- Decompose variance -----
    Q_df = rel.decompose_variance_scaled_from_SS(covariance_matrix, dataset_vec, sub_vec, part_vec)

    # ----- Compare results -----
    compare_resutls(true_df, Q_df)
    
