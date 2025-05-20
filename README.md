# lca-inter-flows
# Intermediate Flows Analysis

A Python toolkit and Jupyter notebook for analyzing intermediate and elementary flows of extractive resource use in life cycle assessment (LCA) projects, including supply‐risk calculation and resource category breakdowns.

## Key Components

- **`main/intermediate_flows.py`**  
  - Defines `CONFIG` for project settings  
  - `get_filtered_activities()`: select activities by keyword  
  - `compute_lca_flows()`: compute intermediate/elementary flows, mass metrics, supply risk, fossil/metal tags  
  - `summarize_records()`: aggregate metrics per activity & category  
  - `plot_scatter()`: scatter intermediate vs. elementary metrics  
  - `plot_activity_contributions()`: 3×4 bar charts of top resource contributions  

- **`intermediate_flows_notebook.ipynb`**  
  Interactive Jupyter notebook to:  
  1. Choose Brightway2 database & filter keyword  
  2. Import functions without auto‐running the script  
  3. Execute analysis pipeline & generate plots  
  4. Display contribution charts for selected activities  

## Quick Start

```bash
# Clone repository and install dependencies
git clone https://github.com/philoubet/lca-inter-flows.git
cd lca-inter-flows
pip install brightway2 pandas matplotlib jupyter
```

## Run the analysis script (optional)
python main/intermediate_flows.py

## Open the interactive notebook
jupyter notebook intermediate_flows_notebook.ipynb
In the notebook, set database index and filter keyword in the first cells.
Run each cell to produce scatter and contribution charts.
