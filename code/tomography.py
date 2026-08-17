import numpy as np
from scipy.fft import fft, fftfreq, ifft
from scipy.ndimage import center_of_mass
from skimage.transform import radon


def make_phantom(n: int, radius: float):
    grid = np.linspace(-1.0, 1.0, n)
    xv, yv = np.meshgrid(grid, grid, indexing="xy")
    image = np.zeros((n, n), dtype=float)
    image[(xv**2 + yv**2) <= radius**2] = 0.6
    image[(xv**2 + yv**2) <= (0.75 * radius) ** 2] = 1.0
    image[(xv**2 + yv**2) <= (0.55 * radius) ** 2] = 0.2
    image[(xv**2 + yv**2) <= (0.35 * radius) ** 2] = 1.6
    image[(xv**2 + yv**2) <= (0.18 * radius) ** 2] = 2.5
    for cx, cy, rad, val in [
        (0.25 * radius, 0.15 * radius, 0.20 * radius, 2.0),
        (-0.30 * radius, -0.10 * radius, 0.06 * radius, 2.8),
        (0.05 * radius, -0.35 * radius, 0.10 * radius, 3.2),
        (-0.10 * radius, 0.35 * radius, 0.10 * radius, 3.6),
    ]:
        image[((xv - cx) ** 2 + (yv - cy) ** 2) <= rad**2] = val
    ellipse = ((xv + 0.20 * radius) / (0.22 * radius)) ** 2 + ((yv - 0.25 * radius) / (0.10 * radius)) ** 2 <= 1.0
    image[ellipse] = 5.0
    return xv, yv, image


def angle_grid(theta_step_degrees: int):
    return np.arange(0.0, 180.0, theta_step_degrees, dtype=float)


def sinogram(image: np.ndarray, theta_degrees: np.ndarray, delta_x: float):
    return radon(image, theta=theta_degrees, circle=True) * delta_x


def detector_grid(n_detector: int, delta_x: float):
    center = n_detector // 2
    return (np.arange(n_detector, dtype=float) - center) * delta_x


def frequency_grid(n_detector: int, delta_t: float):
    return fftfreq(n_detector, d=delta_t)


def centered_padding_slices(n_t: int, n_fft: int):
    if n_fft < n_t:
        raise ValueError("n_fft must be at least n_t")
    left = (n_fft - n_t) // 2
    right = left + n_t
    return left, right


def filter_sinogram(g: np.ndarray, tau: np.ndarray, n_fft: int):
    if tau.shape != (n_fft,):
        raise ValueError("tau must have one value per padded FFT sample")
    left, right = centered_padding_slices(g.shape[0], n_fft)
    padded = np.zeros((n_fft, g.shape[1]), dtype=float)
    padded[left:right, :] = g
    spectrum = fft(padded, axis=0)
    filtered_padded = ifft(tau[:, None] * spectrum, axis=0)
    return np.real(filtered_padded[left:right, :]), spectrum


def backproject(xv: np.ndarray, yv: np.ndarray, g: np.ndarray, t: np.ndarray, theta_degrees: np.ndarray):
    theta = np.deg2rad(theta_degrees)
    delta_theta = float(np.deg2rad(theta_degrees[1] - theta_degrees[0])) if theta_degrees.size > 1 else np.pi
    acc = np.zeros_like(xv, dtype=float)
    for j, angle in enumerate(theta):
        coords = xv * np.cos(angle) - yv * np.sin(angle)
        values = np.interp(coords.ravel(), t, g[:, j], left=0.0, right=0.0)
        acc += values.reshape(xv.shape)
    return acc * delta_theta


def reconstruct_fbp(xv: np.ndarray, yv: np.ndarray, g: np.ndarray, t: np.ndarray, theta_degrees: np.ndarray, tau: np.ndarray, n_fft: int):
    filtered, spectrum = filter_sinogram(g, tau, n_fft)
    return backproject(xv, yv, filtered, t, theta_degrees), filtered, spectrum


def orientation_check(xv: np.ndarray, yv: np.ndarray, image: np.ndarray, reconstruction: np.ndarray):
    ref = np.clip(image - np.percentile(image, 90), 0.0, None)
    est = np.clip(reconstruction - np.percentile(reconstruction, 90), 0.0, None)
    ref_cy, ref_cx = center_of_mass(ref)
    est_cy, est_cx = center_of_mass(est)
    if not np.all(np.isfinite([ref_cx, ref_cy, est_cx, est_cy])):
        raise ValueError("orientation check failed: non-finite center of mass")
    ref_x = float(np.interp(ref_cx, np.arange(xv.shape[1]), xv[0, :]))
    ref_y = float(np.interp(ref_cy, np.arange(yv.shape[0]), yv[:, 0]))
    est_x = float(np.interp(est_cx, np.arange(xv.shape[1]), xv[0, :]))
    est_y = float(np.interp(est_cy, np.arange(yv.shape[0]), yv[:, 0]))
    if np.sign(ref_x) != np.sign(est_x) or np.sign(ref_y) != np.sign(est_y):
        raise ValueError("orientation check failed: likely reflection or rotation")
    return {"phantom_center": [ref_x, ref_y], "ramp_reconstruction_center": [est_x, est_y]}
