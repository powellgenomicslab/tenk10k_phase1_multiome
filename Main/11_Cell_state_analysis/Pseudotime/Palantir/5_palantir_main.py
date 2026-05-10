# Load packages
import scanpy as sc
import palantir
import pandas as pd

# Load object with diffusion maps
print("Reading combined data...")
ad = sc.read_h5ad("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4/CD4_combined_diffusion.h5ad")

# Define terminal states
print("Defining terminal states...")
terminal_states = pd.Series(
    ["Treg", "CTL"],
    index = ["GCCCGAAAGAACGACC-1_S0269_1", "TTTGAGGGTAGACACG-1_S0236_1"],
)

# Run Palantir
print("Running Palantir...")
palantir_result = palantir.core.run_palantir(
    ad,
    early_cell = "TTGCAGAAGTCATCTG-1_S0289_1",  
    knn = 50, 
    num_waypoints = 8000,
    terminal_states = terminal_states,
    n_jobs = 48
)

# Save results to AnnData object
print("Saving AnnData object with Palantir results...")
ad.write("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4/CD4_combined_palantir.h5ad")
print("Finished!")
