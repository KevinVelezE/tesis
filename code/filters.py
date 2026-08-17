import numpy as np


def ramp(sigma: np.ndarray) -> np.ndarray:
    return np.abs(sigma)


def triangular_symbol(sigma: np.ndarray, sigma_c: float) -> np.ndarray:
    sigma_c_grid = np.fft.fftshift(sigma)
    sigma_abs = np.abs(sigma_c_grid)
    tau_c = np.zeros_like(sigma_c_grid)
    mask1 = sigma_abs <= sigma_c
    tau_c[mask1] = sigma_abs[mask1]
    mask2 = (sigma_abs > sigma_c) & (sigma_abs <= 2.0 * sigma_c)
    tau_c[mask2] = 2.0 * sigma_c - sigma_abs[mask2]
    return np.fft.ifftshift(tau_c)


def decay_symbol(sigma: np.ndarray, sigma_c: float, alpha: float) -> np.ndarray:
    sigma_c_grid = np.fft.fftshift(sigma)
    sigma_abs = np.abs(sigma_c_grid)
    tau_c = np.zeros_like(sigma_c_grid)
    mask1 = sigma_abs <= sigma_c
    tau_c[mask1] = sigma_abs[mask1]
    mask2 = sigma_abs > sigma_c
    tau_c[mask2] = sigma_c * np.exp(-alpha * (sigma_abs[mask2] - sigma_c))
    return np.fft.ifftshift(tau_c)


def symbols(sigma: np.ndarray):
    sigma_abs = np.abs(sigma)
    sigma_max = float(np.max(sigma_abs))
    normalized = sigma_abs / sigma_max
    cutoff = 0.6

    compact = np.zeros_like(sigma)
    mask = normalized < cutoff
    compact[mask] = sigma_abs[mask] * (1.0 - normalized[mask] / cutoff)

    return {
        "ramp": (ramp(sigma), r"$|\sigma|$"),
        "A": (sigma_abs * (1.0 + np.sin(sigma_abs) / (1.0 + sigma_abs)), r"$|\sigma|(1+\sin|\sigma|/(1+|\sigma|))$"),
        "B": (np.sqrt(0.3**2 + sigma_abs**2), r"$\sqrt{0.3^2+|\sigma|^2}$"),
        "C": (compact, r"corte compacto, $\sigma_c=0.6\sigma_{\max}$"),
        "D": (triangular_symbol(sigma, sigma_c=20.0), r"triangular, $\sigma_c=20$"),
        "E": (decay_symbol(sigma, sigma_c=20.0, alpha=0.08), r"exponencial, $\sigma_c=20$, $\alpha=0.08$"),
        "F": (decay_symbol(sigma, sigma_c=20.0, alpha=0.50), r"exponencial, $\sigma_c=20$, $\alpha=0.50$"),
    }
