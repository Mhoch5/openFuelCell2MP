import os
import subprocess
from of_modifier import update_all_openfuelcell_files

def prepare_openfoam_case(case_dir):
    """
    Checks if the mesh and parallel decomposition have already been created.
    If not, it executes the necessary Makefile commands.
    """
    # 1. Check for Meshing
    # The 'make mesh' command pipes preprocessing output to 'log.pre'.
    mesh_log = os.path.join(case_dir, "log.pre")
    if os.path.exists(mesh_log):
        print("[-] Mesh already exists (found log.pre). Skipping 'make mesh'.")
    else:
        print("[+] Mesh not found. Executing 'make mesh'...")
        subprocess.run(["make", "mesh"], cwd=case_dir, check=True)

    # 2. Check for Decomposition and Parallel Setup
    # The 'make parallel' command pipes output to 'log.parallel'[cite: 15].
    parallel_log = os.path.join(case_dir, "log.parallel")
    if os.path.exists(parallel_log):
        print("[-] Parallel decomposition already exists (found log.parallel). Skipping setup.")
    else:
        print("[+] Decomposition not found. Executing 'make decompose'...")
        # 'make decompose' generates the cellID for multiple regions[cite: 16].
        subprocess.run(["make", "decompose"], cwd=case_dir, check=True)
        
        print("[+] Executing 'make parallel'...")
        subprocess.run(["make", "parallel"], cwd=case_dir, check=True)


def run_openfuelcell(case_dir):
    """
    Executes the OpenFOAM solver in parallel and scans for masked crash codes.
    """
    print("[+] Starting openFuelCell solver in parallel...")
    
    result = subprocess.run(
        ["make", "run"], 
        cwd=case_dir,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True
    )
    
    # Combine standard output and error to scan for crash signatures
    full_output = result.stdout + result.stderr
    
    # Check for OpenFOAM fatal errors, MPI aborts, or segfaults
    crash_signatures = ["FOAM FATAL", "Process received signal", "FOAM exiting"]
    if any(sig in full_output for sig in crash_signatures) or result.returncode != 0:
        print("Error: The OpenFOAM solver crashed.")
        # Print the last 15 lines of the output to diagnose the crash directly in the loop
        print("\n".join(full_output.splitlines()[-15:]))
        return False
        
    print("[-] Solver completed successfully.")
    return True

def get_latest_ibar(file_path):
    latest_ibar = None
    
    if not os.path.exists(file_path):
        print(f"Warning: Log file {file_path} does not exist.")
        return None

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            if 'ibar:' in line:
                latest_ibar = line.strip()

    if latest_ibar is None:
        print("Warning: 'ibar:' not found in log file.")
        return None

    parts = latest_ibar.split()
    return (float(parts[1]), float(parts[3]))

# ==============================================================================
# INTEGRATION EXAMPLE FOR YOUR OPTIMIZATION LOOP
# ==============================================================================
if __name__ == "__main__":
    OF_CASE_DIR = "run/PEMFC_capillaryPressure"
    
    # Step 1: Check and prepare the mesh/decomposition (Run once)
    prepare_openfoam_case(OF_CASE_DIR)
    
    # Step 2: In your optimization loop, update the dictionaries here
    # import of_modifier
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    RAW_DATA_CSV = os.path.join(BASE_DIR, "raw_data.csv")
    PARAMS_CSV = os.path.join(BASE_DIR, "params.csv")

    update_all_openfuelcell_files(
            case_dir=OF_CASE_DIR,
            raw_data_csv=RAW_DATA_CSV,
            params_csv=PARAMS_CSV,
            target_zones=["cathode", "mpl"]
        )
    
    # Step 3: Run the solver
    success = run_openfuelcell(OF_CASE_DIR)
    
    if success:
        print("Ready to extract current density for the optimizer.")
        # Step 4: Add your CSV reading logic here to extract the current density 
        # and return it to the optimizer.

    # Example usage:
    result = get_latest_ibar(os.path.join(OF_CASE_DIR, 'log.run'))
    print(result)