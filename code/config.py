from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class Config:
    n: int = 256
    radius: float = 0.9
    theta_step_degrees: int = 2
    m: int = 90
    circle: bool = True
    seed: int = 42
    noise_level: float = 0.05
    cmap: str = "bone"
    dpi_png: int = 300

    root: Path = Path(__file__).resolve().parents[1]
    figure_dir: Path = root / "Figures" / "Chapter3"
    table_dir: Path = root / "Tables" / "Chapter3"
    result_dir: Path = root / "Results" / "Chapter3"

    def metadata(self):
        data = asdict(self)
        data["figure_dir"] = "Figures/Chapter3"
        data["table_dir"] = "Tables/Chapter3"
        data["result_dir"] = "Results/Chapter3"
        data.pop("root")
        return data


CFG = Config()
