from fetch_sequence import main as fetch
from smith_waterman import main as align
from run_mafft import main as mafft

if __name__ == "__main__":
    print("=== Step 1: Fetch Sequences ===")
    fetch()
    print("=== Step 2: Pairwise Alignment ===")
    align()
    print("=== Step 3: MAFFT ===")
    mafft()