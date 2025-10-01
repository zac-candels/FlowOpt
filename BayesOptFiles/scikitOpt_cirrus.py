"""
What happens in each iteration of optimization procedure:
 1. Write `input.txt` with current parameters.
 2. Run C++ simulation (e.g. `./run.exe`).
 3. Read in value output by the c++ function
 4. Parse and return objective.
"""
import subprocess
import numpy as np
import os
import re
import time

# --- Utility: run one simulation & analysis ---

def run_simulation(a1: float, b1: float, c1: float) -> float:
    """
    Run one simulation via SLURM, wait until it finishes, then run analysis.
    Returns: negative velocity (for minimization in gp_minimize).
    """
    full_path = "../ratchetGeom/examples/binary/superhydrophobic_wellbalanced"

    # --- Step 1: Write input.txt ---
    with open(os.path.join(full_path, 'input.txt'), 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        if line.strip().startswith("a1="):
            new_lines.append(f"a1={a1:.2f} \n")
        elif line.strip().startswith("b1="):
            new_lines.append(f"b1={b1:.2f} \n")
        elif line.strip().startswith("c1="):
            new_lines.append(f"c1={c1:.2f} \n")
        else:
            new_lines.append(line)

    with open(os.path.join(full_path, 'input.txt'), 'w') as f:
        f.writelines(new_lines)

    # --- Step 2: Submit job with sbatch ---
    submit = subprocess.run(
        ["sbatch", "submit.slurm"],
        cwd=full_path,
        capture_output=True,
        text=True
    )
    submit.check_returncode()

    # Example output: "Submitted batch job 12345"
    m = re.search(r"Submitted batch job (\d+)", submit.stdout)
    if not m:
        raise RuntimeError(f"Couldn't parse job ID from sbatch output: {submit.stdout}")
    job_id = m.group(1)
    print(f"[run_simulation] Submitted job {job_id} for a1={a1:.2f}, b1={b1:.2f}, c1={c1:.2f}")

    # --- Step 3: Wait until job completes ---
    while True:
        check = subprocess.run(
            ["squeue", "-j", job_id],
            capture_output=True,
            text=True
        )
        if job_id not in check.stdout:
            break  # job finished
        time.sleep(5)  # wait 5 seconds before checking again

    print(f"[run_simulation] Job {job_id} finished. Running analysis...")

    # --- Step 4: Run analysis script ---
    analysis = subprocess.run(
        ["python", "Analysis.py"],
        cwd=full_path,
        capture_output=True,
        text=True,
        timeout=300
    )
    analysis.check_returncode()

    lines = analysis.stdout.strip().splitlines()
    if not lines:
        raise RuntimeError("Analysis.py produced no output")

    last = lines[-1]
    m = re.search(r"spread_dist\s*=\s*:? *([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)", last)
    if not m:
        raise RuntimeError(f"Couldn't parse a number from Analysis.py output: {last!r}")

    val = float(m.group(1))
    print(f"[run_simulation] Parsed velocity = {val:.4f}")
    return -val  # minus sign because gp_minimize minimizes



from skopt import gp_minimize
from skopt.space import Real
from skopt.callbacks import VerboseCallback

def print_best_so_far(res):
    current_best_idx = int(np.argmin(res.func_vals))
    current_best_val = -res.func_vals[current_best_idx]
    current_best_params = res.x_iters[current_best_idx]
    a1, b1, c1 = current_best_params
    print(f"[Iter {len(res.func_vals)}] Best objective so far: {current_best_val:.4f} at (a1={a1:.2f}, b1={b1:.2f}, c1={c1:.2f})")


# Search space
search_space = [
    Real(3, 20, name='a1'),
    Real(3, 20, name='b1'),
    Real(3, 20, name='c1')
]

def objective_sk(params):
    a1, b1, c1 = params
    try:
        val = run_simulation(a1, b1, c1)
    except Exception as e:
        print(f"Error at {params}: {e}")
        return 1e6
    return val  



def run_skopt():
    result = gp_minimize(
        func=objective_sk,
        dimensions=search_space,
        acq_func='EI',
        n_initial_points=5,
        n_calls=100,
        random_state=42,
        callback=[print_best_so_far]  # <-- Use your callback here
    )
    a1, b1, c1 = result.x
    print("== skopt best ==")
    print(f"a1={a1:.2e}, b1={b1:.2e}, c1={c1:.2e}")
    print(f"Max value={-result.fun:.4f}")


run_skopt()
