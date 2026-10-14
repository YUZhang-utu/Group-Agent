"""Entry point for the matched block-screening benchmark."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from aidd_agent.matched_block_benchmark import main


if __name__ == '__main__':
    main()
