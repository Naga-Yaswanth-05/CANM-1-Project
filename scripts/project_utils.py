"""Shared path and analytical-function helpers for the project scripts."""
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
ACTIVE_DATASET_FILE = OUTPUT_DIR / "active_dataset.txt"

DATASETS = {
    "data.csv": {
        "key": "sin",
        "label": "sin(5x)",
    },
    "group1.csv": {
        "key": "group1",
        "label": "1/(1+25x²)",
    },
    "tanh_data.csv": {
        "key": "tanh",
        "label": "tanh(10x)",
    },
}


def active_dataset_filename() -> str:
    """Return the dataset last used by the C++ pipeline."""
    if not ACTIVE_DATASET_FILE.is_file():
        raise FileNotFoundError(
            f"{ACTIVE_DATASET_FILE} is missing. Run the C++ main program first."
        )
    filename = ACTIVE_DATASET_FILE.read_text(encoding="utf-8").strip()
    if filename not in DATASETS:
        raise ValueError(
            f"Unknown active dataset {filename!r} in {ACTIVE_DATASET_FILE}. "
            f"Expected one of: {', '.join(DATASETS)}"
        )
    return filename


def dataset_key(filename: str) -> str:
    try:
        return DATASETS[filename]["key"]
    except KeyError as exc:
        raise ValueError(f"Unsupported dataset: {filename!r}") from exc


def dataset_label(filename: str) -> str:
    try:
        return DATASETS[filename]["label"]
    except KeyError as exc:
        raise ValueError(f"Unsupported dataset: {filename!r}") from exc


def exact_function(filename: str, x):
    """Evaluate the analytical function represented by a supported dataset."""
    x = np.asarray(x, dtype=float)
    if filename == "data.csv":
        return np.sin(5.0 * x)
    if filename == "group1.csv":
        return 1.0 / (1.0 + 25.0 * x**2)
    if filename == "tanh_data.csv":
        return np.tanh(10.0 * x)
    raise ValueError(f"Unsupported dataset: {filename!r}")


def exact_derivative(filename: str, x):
    """Evaluate the analytical derivative for a supported dataset."""
    x = np.asarray(x, dtype=float)
    if filename == "data.csv":
        return 5.0 * np.cos(5.0 * x)
    if filename == "group1.csv":
        return -50.0 * x / (1.0 + 25.0 * x**2) ** 2
    if filename == "tanh_data.csv":
        y = np.tanh(10.0 * x)
        return 10.0 * (1.0 - y**2)
    raise ValueError(f"Unsupported dataset: {filename!r}")
