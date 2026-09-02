import numpy as np
from scipy.stats import entropy
import matplotlib.pyplot as plt
from parameters import(
    VARIABLE_NUMBER,
    N_MAJOR_REPS,
    DATASET_SIZE
)
from models import(
    true_normal_dist
)
from dataset_sampling import empirical_sample

# Change this to change model
true_dist = true_normal_dist

# Code 
mean_divergence = 0
for _ in N_MAJOR_REPS:
    
    fitted_dist = empirical_sample(
        true_dist,
        VARIABLE_NUMBER,
        DATASET_SIZE
    )

    # Shared grid — use ppf to get a sensible range covering both distributions' mass
    lo = min(true_dist.ppf(0.0001), fitted_dist.ppf(0.0001))
    hi = max(true_dist.ppf(0.9999), fitted_dist.ppf(0.9999))
    x = np.linspace(lo, hi, 10000)

    p = true_dist.pdf(x)
    q = fitted_dist.pdf(x)

    # Avoid division by zero
    eps = 1e-10
    kl_div = entropy(p + eps, q + eps)

    mean_divergence += kl_div

    last_fitted_dist = fitted_dist

mean_divergence = mean_divergence / N_MAJOR_REPS
print(f"Mean KL divergence over {N_MAJOR_REPS} reps: {mean_divergence:.6f}")

# --- Representative plot (uses the last rep's fitted distribution) ---
lo = min(true_dist.ppf(0.0001), last_fitted_dist.ppf(0.0001))
hi = max(true_dist.ppf(0.9999), last_fitted_dist.ppf(0.9999))
x = np.linspace(lo, hi, 1000)

plt.figure(figsize=(8, 5))
plt.plot(
    x,
    true_dist.pdf(x),
    label="True distribution",
    linewidth=2
)
plt.plot(
    x,
    last_fitted_dist.pdf(x),
    label="Fitted (empirical) distribution",
    linewidth=2,
    linestyle="--"
)
plt.fill_between(x, true_dist.pdf(x), alpha=0.15)
plt.fill_between(x, last_fitted_dist.pdf(x), alpha=0.15)
plt.title(
    "True vs. Fitted Distribution"
    f"(mean KL divergence = {mean_divergence:.4f})"
)
plt.xlabel("Effect size")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.show()



