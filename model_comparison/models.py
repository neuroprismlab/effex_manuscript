from scipy.stats import norm
import numpy as np
from parameters import (
    MEAN,
    VARIANCE
)


# Normal to normal comparisson
def normal_model(mean, std_deviation):

    return norm(
        loc=mean,
        scale=std_deviation
    )


true_normal_dist = normal_model(MEAN, np.sqrt(VARIANCE))
