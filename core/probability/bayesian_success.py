import numpy as np
from scipy.stats import beta

class BayesianSuccessModel:
    def __init__(self, prior_success=2, prior_failure=2):
        """
        prior_success / prior_failure:
        Low values = weak prior
        High values = strong belief before seeing data
        """
        self.prior_success = prior_success
        self.prior_failure = prior_failure

    def update(self, recovery_score, trials=20, noise=0.1):
        """
        Convert recovery score into probabilistic evidence.

        trials: virtual experiments
        noise: uncertainty in mapping score → success
        """
        score = np.clip(recovery_score + np.random.normal(0, noise), 0, 1)

        successes = int(round(score * trials))
        failures = trials - successes

        alpha = self.prior_success + successes
        beta_param = self.prior_failure + failures

        return alpha, beta_param

    def probability(self, alpha, beta_param, samples=10000):
        dist = beta(alpha, beta_param)
        mean = dist.mean()
        ci = dist.interval(0.95)

        return {
            "mean": float(mean),
            "ci_low": float(ci[0]),
            "ci_high": float(ci[1])
        }
