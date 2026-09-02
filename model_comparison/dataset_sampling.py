import numpy as np
from scipy.stats import norm


def empirical_sample(
    true_distribution,
    variable_number,
    subject_number
):

    true_effects = true_distribution.rvs(variable_number)

    noise = norm.rvs(
        loc=0,
        scale=1/np.sqrt(subject_number),
        size=variable_number)
    dataset_effects = true_effects + noise

    mean_fit, std_fit = norm.fit(dataset_effects)
    
    # Return true effects and emperical distribution
    return norm(
        loc=mean_fit,
        scale=std_fit,
    )
