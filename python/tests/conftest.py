import sys
from pathlib import Path

# Make the package importable without installation.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SYNTHETIC = Path(__file__).resolve().parents[2] / "inst" / "extdata" / "synthetic"
