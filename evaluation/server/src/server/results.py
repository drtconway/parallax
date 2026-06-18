import uuid
from pathlib import Path
from .config import RESULTS_DIR


def new_result_dir() -> tuple[str, Path]:
    result_id = str(uuid.uuid4())
    path = RESULTS_DIR / result_id
    path.mkdir(parents=True)
    return result_id, path
