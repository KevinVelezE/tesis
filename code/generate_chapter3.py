import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scipy
import skimage
from matplotlib import colors

from config import CFG
from filters import symbols
from metrics import relative_l2, rmse
from tomography import (
    angle_grid,
    backproject,
    detector_grid,
    filter_sinogram,
    frequency_grid,
    make_phantom,
    orientation_check,
    reconstruct_fbp,
    sinogram,
)


def save_figure(fig, path: Path):
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=CFG.dpi_png, bbox_inches="tight")
    plt.close(fig)


def finite_check(name: str, array: np.ndarray):
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains NaN or infinite values")


def check_symbols(symbol_map, sigma):
    zero_index = int(np.where(sigma == 0.0)[0][0])
    paired = np.array([np.where(np.isclose(sigma, -value))[0][0] for value in sigma if np.any(np.isclose(sigma, -value))])
    zero_values = {}
    for name, (tau, _) in symbol_map.items():
        if tau.shape != sigma.shape:
            raise ValueError(f"{name} has invalid shape")
        finite_check(name, tau)
        for i in paired:
            j = int(np.where(np.isclose(sigma, -sigma[i]))[0][0])
            if not np.isclose(tau[i], tau[j]):
                raise ValueError(f"{name} is not even on paired fftfreq samples")
        zero_values[name] = float(tau[zero_index])
    if zero_values["ramp"] != 0.0 or zero_values["B"] != 0.3:
        raise ValueError("unexpected zero-frequency symbol values")
    return zero_values


def check_column_application(g, tau, n_fft):
    filtered, _ = filter_sinogram(g, tau, n_fft)
    first_column, _ = filter_sinogram(g[:, :1], tau, n_fft)
    if not np.allclose(filtered[:, 0], first_column[:, 0]):
        raise ValueError("filter is not applied independently by detector columns")


def compute_run():
    theta = angle_grid(CFG.theta_step_degrees)
    if theta.size != CFG.m:
        raise ValueError("theta_step_degrees does not produce M angles")

    xv, yv, phantom = make_phantom(CFG.n, CFG.radius)
    delta_x = float(xv[0, 1] - xv[0, 0])
    g = sinogram(phantom, theta, delta_x)
    n_t = int(g.shape[0])
    n_fft = 1024
    t = detector_grid(g.shape[0], delta_x)
    delta_t = float(t[1] - t[0])
    sigma = frequency_grid(n_fft, delta_t)
    delta_sigma = float(abs(sigma[1] - sigma[0]))
    delta_theta = float(np.deg2rad(CFG.theta_step_degrees))
    symbol_map = symbols(sigma)
    zero_values = check_symbols(symbol_map, sigma)

    expected_n_t = CFG.n
    if phantom.shape != (CFG.n, CFG.n) or g.shape != (expected_n_t, CFG.m):
        raise ValueError("invalid phantom or sinogram dimensions")
    if n_t != 256 or n_fft != 1024:
        raise ValueError("unexpected detector or padded FFT length")
    finite_check("phantom", phantom)
    finite_check("sinogram", g)

    rng = np.random.default_rng(CFG.seed)
    sigma_noise = CFG.noise_level * float(np.max(np.abs(g)))
    noise = rng.normal(0.0, sigma_noise, size=g.shape)
    g_noisy = g + noise
    finite_check("noisy sinogram", g_noisy)
    check_column_application(g, symbol_map["ramp"][0], n_fft)

    bp = backproject(xv, yv, g, t, theta)
    recon_noiseless = {}
    recon_noisy = {}
    filtered_examples = {}
    for name, (tau, _) in symbol_map.items():
        if tau.shape[0] != n_fft:
            raise ValueError(f"M_tau for {name} has invalid padded FFT dimension")
        recon, filtered, _ = reconstruct_fbp(xv, yv, g, t, theta, tau, n_fft)
        recon_n, _, _ = reconstruct_fbp(xv, yv, g_noisy, t, theta, tau, n_fft)
        finite_check(f"reconstruction {name}", recon)
        finite_check(f"noisy reconstruction {name}", recon_n)
        recon_noiseless[name] = recon
        recon_noisy[name] = recon_n
        filtered_examples[name] = filtered

    orient = orientation_check(xv, yv, phantom, recon_noiseless["ramp"])
    mask = xv**2 + yv**2 <= CFG.radius**2
    rows = []
    for noise_case, reconstructions in [("noiseless", recon_noiseless), ("noisy", recon_noisy)]:
        for name, recon in reconstructions.items():
            rows.append(
                {
                    "case": noise_case,
                    "filter": name,
                    "rmse": rmse(phantom, recon, mask),
                    "relative_l2": relative_l2(phantom, recon, mask),
                    "mean_ratio_on_mask": float(np.mean(recon[mask]) / np.mean(phantom[mask])),
                }
            )

    return {
        "xv": xv,
        "yv": yv,
        "phantom": phantom,
        "theta": theta,
        "sinogram": g,
        "sinogram_noisy": g_noisy,
        "sigma_noise": sigma_noise,
        "t": t,
        "sigma": sigma,
        "n_t": n_t,
        "n_fft": n_fft,
        "delta_x": delta_x,
        "delta_t": delta_t,
        "delta_sigma": delta_sigma,
        "delta_theta": delta_theta,
        "symbols": symbol_map,
        "zero_values": zero_values,
        "bp": bp,
        "recon_noiseless": recon_noiseless,
        "recon_noisy": recon_noisy,
        "filtered_examples": filtered_examples,
        "orientation": orient,
        "mask": mask,
        "metrics": rows,
    }


def plot_image(ax, image, title, cmap, vmin=None, vmax=None):
    im = ax.imshow(image, cmap=cmap, origin="lower", vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    return im


def display_filter_name(name: str) -> str:
    return {
        "ramp": "rampa",
        "A": "A",
        "B": "B",
        "C": "C: corte compacto",
        "D": "D: triangular",
        "E": "E: exponencial",
        "F": "F: exponencial",
    }[name]


def display_symbol_label(name: str, label: str) -> str:
    return f"{'rampa' if name == 'ramp' else name}: {label}"


def display_case_name(name: str) -> str:
    return {"noiseless": "sin ruido", "noisy": "con ruido"}[name]


def create_figures(data):
    cmap = CFG.cmap
    phantom = data["phantom"]
    g = data["sinogram"]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    im0 = plot_image(axes[0], phantom, "Fantoma asimétrico", cmap)
    plt.colorbar(im0, ax=axes[0], fraction=0.046)
    im1 = axes[1].imshow(g, cmap=cmap, aspect="auto", origin="lower", extent=[0, 178, data["t"][0], data["t"][-1]])
    axes[1].set_title("Sinograma")
    axes[1].set_xlabel(r"$\theta$ (grados)")
    axes[1].set_ylabel(r"$t$")
    plt.colorbar(im1, ax=axes[1], fraction=0.046)
    save_figure(fig, CFG.figure_dir / "phantom_and_sinogram")

    ramp = data["recon_noiseless"]["ramp"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    bp_vmin, bp_vmax = np.percentile(data["bp"], [1, 99])
    fbp_vmin, fbp_vmax = np.percentile(ramp, [1, 99])
    im0 = plot_image(axes[0], data["bp"], "Retroproyección sin filtrar", cmap, bp_vmin, bp_vmax)
    im1 = plot_image(axes[1], ramp, "FBP con rampa", cmap, fbp_vmin, fbp_vmax)
    plt.colorbar(im0, ax=axes[0], fraction=0.046, label="escala BP")
    plt.colorbar(im1, ax=axes[1], fraction=0.046, label="escala FBP")
    save_figure(fig, CFG.figure_dir / "bp_vs_fbp")

    fig, ax = plt.subplots(figsize=(8, 5))
    sigma_c = np.fft.fftshift(data["sigma"])
    for name, (tau, label) in data["symbols"].items():
        ax.plot(sigma_c, np.fft.fftshift(tau), label=display_symbol_label(name, label), linewidth=1.7)
    ax.set_title("Rampa y símbolos exploratorios A-F")
    ax.set_xlabel(r"$\sigma$")
    ax.set_ylabel(r"$\tau(\sigma)$")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)
    save_figure(fig, CFG.figure_dir / "symbols_ramp_A_F")

    names = list(data["symbols"].keys())
    for key, output_name, title in [
        ("recon_noiseless", "reconstructions_noiseless", "Reconstrucciones sin ruido"),
        ("recon_noisy", "reconstructions_noisy", "Reconstrucciones con ruido"),
    ]:
        stack = np.stack([data[key][name] for name in names])
        vmin, vmax = np.percentile(stack, [1, 99])
        fig, axes = plt.subplots(2, 4, figsize=(13, 7))
        for ax, name in zip(axes.ravel(), names):
            im = plot_image(ax, data[key][name], display_filter_name(name), cmap, vmin, vmax)
        im = plot_image(axes.ravel()[-1], phantom, "Fantoma", cmap, vmin, vmax)
        fig.suptitle(title)
        plt.colorbar(im, ax=axes.ravel().tolist(), fraction=0.025)
        save_figure(fig, CFG.figure_dir / output_name)


def write_metrics(rows):
    csv_path = CFG.result_dir / "metrics.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["case", "filter", "rmse", "relative_l2", "mean_ratio_on_mask"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    tex_path = CFG.table_dir / "metrics_table.tex"
    with tex_path.open("w", encoding="utf-8") as f:
        f.write("\\begin{tabular}{llrrr}\n")
        f.write("\\toprule\n")
        f.write("Caso & Filtro & RMSE & Error relativo $L^2$ & Razón media \\\\\n")
        f.write("\\midrule\n")
        for row in rows:
            f.write(f"{display_case_name(row['case'])} & {display_filter_name(row['filter'])} & {row['rmse']:.6f} & {row['relative_l2']:.6f} & {row['mean_ratio_on_mask']:.6f} \\\\\n")
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")


def write_report(data, metric_hash):
    report = CFG.result_dir / "normalization_report.md"
    ramp_ratio = [row for row in data["metrics"] if row["case"] == "noiseless" and row["filter"] == "ramp"][0]["mean_ratio_on_mask"]
    with report.open("w", encoding="utf-8") as f:
        f.write("# Auditoría de normalización\n\n")
        f.write("Rampa implementada: `tau_rampa(sigma)=|sigma|`; no se usa `scale=0.5`.\n\n")
        f.write("El sinograma se calcula con `skimage.transform.radon(..., circle=True)`.\n\n")
        f.write(f"N_t = {data['n_t']}\n\n")
        f.write(f"N_fft = {data['n_fft']}\n\n")
        f.write(f"Delta x = {data['delta_x']:.12g}\n\n")
        f.write(f"Delta t = {data['delta_t']:.12g}\n\n")
        f.write(f"Delta sigma = {data['delta_sigma']:.12g}\n\n")
        f.write(f"Delta theta = {data['delta_theta']:.12g} rad ({CFG.theta_step_degrees} grados)\n\n")
        f.write("El sinograma de longitud `N_t` se rellena de forma centrada hasta `N_fft`, se transforma con `fft(..., axis=0)`, se multiplica como `tau[:, None] * G`, se retorna con `ifft(..., axis=0)` y se recorta al intervalo detector original.\n\n")
        f.write("La grilla de frecuencias es `scipy.fft.fftfreq(N_fft, d=Delta t)`. Los símbolos se multiplican en el orden nativo de `fftfreq`; `fftshift` se usa solo para graficar.\n\n")
        f.write("La retroproyección evalúa `t=x cos(theta)-y sin(theta)` por interpolación unidimensional en detector y aplica cuadratura angular con `Delta theta`. El signo menos alinea `y` cartesiano con las filas de imagen usadas por `skimage.transform.radon`; la prueba de orientación con fantoma asimétrico verifica esta elección.\n\n")
        f.write(f"Ruido: `sigma_noise = 0.05 * max(abs(g)) = {data['sigma_noise']:.12g}`; la misma realización ruidosa se usa para todos los filtros.\n\n")
        f.write(f"Valores en frecuencia cero: `{json.dumps(data['zero_values'], sort_keys=True)}`. En particular, `tau_B(0)=0.3`.\n\n")
        f.write(f"Centros de la prueba de orientación: `{json.dumps(data['orientation'])}`. Los signos coinciden, sin rotación/reflexión detectada.\n\n")
        f.write(f"No se aplica corrección empírica ni normalización por reconstrucción. La razón media de la rampa sin ruido en el disco de evaluación es {ramp_ratio:.6f}; se conserva solo como diagnóstico, no como factor de escala.\n\n")
        f.write(f"Hash reproducible de métricas redondeadas: `{metric_hash}`.\n")


def write_metadata(data, metric_hash):
    metadata = {
        "parameters": CFG.metadata(),
        "versions": {
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "scikit-image": skimage.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "sigma_noise": data["sigma_noise"],
        "n_t": data["n_t"],
        "n_fft": data["n_fft"],
        "circle": CFG.circle,
        "delta_t": data["delta_t"],
        "delta_sigma": data["delta_sigma"],
        "delta_theta": data["delta_theta"],
        "metric_hash": metric_hash,
        "orientation_check": data["orientation"],
        "zero_frequency_values": data["zero_values"],
    }
    with (CFG.result_dir / "run_metadata.json").open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, sort_keys=True)


def metric_hash(rows):
    rounded = []
    for row in rows:
        rounded.append({key: (round(value, 12) if isinstance(value, float) else value) for key, value in row.items()})
    payload = json.dumps(rounded, sort_keys=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def main():
    for directory in (CFG.figure_dir, CFG.table_dir, CFG.result_dir):
        directory.mkdir(parents=True, exist_ok=True)

    data = compute_run()
    repeated = compute_run()
    first_hash = metric_hash(data["metrics"])
    second_hash = metric_hash(repeated["metrics"])
    if first_hash != second_hash:
        raise ValueError("metrics are not reproducible with the same seed")

    create_figures(data)
    write_metrics(data["metrics"])
    write_report(data, first_hash)
    write_metadata(data, first_hash)
    print(json.dumps({"status": "ok", "metric_hash": first_hash, "metrics": data["metrics"]}, indent=2))


if __name__ == "__main__":
    main()
