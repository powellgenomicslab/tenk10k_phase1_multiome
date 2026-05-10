# Load packages
import palantir
import pandas as pd
import numpy as np

# Load multiscale space
print("Reading multiscale space...")
ms_data = pd.read_csv("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CSVs/DM_multiscale_space.csv", index_col=0)

# Identify terminal states
print("Identifying terminal states...")
terminal_states, excluded_boundaries = palantir.core.identify_terminal_states(
    ms_data, 
    knn = 50, 
    early_cell = "CACAACAAGGCCTCGT-1_S0334_1",  
    num_waypoints = 8000,
    n_jobs = 48
)

# Print terminal states
print("Terminal states:")
print(terminal_states)

print("Excluded boundaries:")
print(excluded_boundaries)

# Save terminal states
terminal_df = pd.DataFrame(terminal_states, columns=["terminal_cell"])
excluded_df = pd.DataFrame(excluded_boundaries, columns=["excluded_boundary_cell"])

terminal_df.to_csv("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CSVs/terminal_states.csv", index=False)
excluded_df.to_csv("/g/data/fy54/od8037/CHIP_Project/Palantir/CD4_All/CSVs/excluded_boundaries.csv", index=False)

print(f"Saved {len(terminal_df)} terminal states and {len(excluded_df)} excluded boundaries.")