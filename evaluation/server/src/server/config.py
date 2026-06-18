import os
from pathlib import Path

PARALLAX_URL = os.environ.get("PARALLAX_URL", "http://localhost:8080")
_default_results = Path(__file__).parent.parent.parent / "results"
RESULTS_DIR = Path(os.environ.get("RESULTS_DIR", _default_results))
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
