import numpy as np


def rmse(reference: np.ndarray, estimate: np.ndarray, mask: np.ndarray) -> float:
    diff = estimate[mask] - reference[mask]
    return float(np.sqrt(np.mean(diff**2)))


def relative_l2(reference: np.ndarray, estimate: np.ndarray, mask: np.ndarray) -> float:
    diff = estimate[mask] - reference[mask]
    denom = np.linalg.norm(reference[mask])
    if denom == 0.0:
        raise ValueError("reference norm is zero on the metric mask")
    return float(np.linalg.norm(diff) / denom)
