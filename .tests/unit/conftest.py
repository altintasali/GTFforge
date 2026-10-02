"""Make workflow/scripts importable as plain modules.

The scripts are standalone CLIs invoked by path from the rules; each inserts
its own directory into sys.path to import gtf_utils, so adding that directory
here is enough to import them in tests.
"""
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
