# Updated LCA Flow Script with Enhanced Plot Labeling and Risky Mass Plotting
import logging
import random
from brightway2 import projects, Database, LCA
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- Configuration ---
CONFIG = {
    "project": "intermediate_flows",
    "databases": ["ecoinvent-3.10-cutoff"],
    "matched_csv": "processed_b_public_with_percentages.csv",
    "risk_csv": "model_csv_geopolrisk_with_colors_corrected.csv",
    
    "filter_keyword": "market for battery,",
    
    "sample_size": 50,
    
    "fossil_resources": [
    'Coal, brown', 'Coal, hard', 'Gas, natural', 'Shale', 'Oil, crude',
    'lignite', 'waste natural gas, sour', 'sweet gas, burned in gas turbine',
    'sour gas, burned in gas turbine', "refinery gas", "petrol, unleaded",
    "naphtha", "liquefied petroleum gas", "light fuel oil", "diesel",
    "petroleum combustion, in drilling tests", "petroleum",
    "natural gas, high pressure", "natural gas, vented", "burnt shale",
    "hard coal", "hard coal, run-of-mine", "waste natural gas, sweet", "sweet gas"
    ],
}

# --- Setup ---
projects.set_current(CONFIG["project"])
matches_df = pd.read_csv(CONFIG["matched_csv"], delimiter=";")
supply_risk_factors_elem = pd.read_csv(CONFIG["risk_csv"], delimiter=";")
supply_risk_factors_int = (
    matches_df[["Activity", "Product", "Geography", "Matched Substance", "Final_Percentage"]]
    .merge(
        supply_risk_factors_elem[["Color","Substance", "Supplyrisk"]],
        left_on="Matched Substance",
        right_on="Substance",
        how="left",
    )
    .drop(columns="Substance")
)
supply_risk_factors_int["Supplyrisk"] *= supply_risk_factors_int["Final_Percentage"] / 100

# Pre-load metadata
activity_meta = {db: {act.key: act for act in Database(db)} for db in CONFIG["databases"]}
bios_meta = {flow.key: {"name": flow["name"], "type": flow["type"]}
             for flow in Database("biosphere3")}

# --- Functions definition ---


def get_filtered_activities(db_name):
    db = Database(db_name)
    acts = [act for act in db if CONFIG["filter_keyword"] in act["name"].lower()]
    return random.sample(acts, min(len(acts), CONFIG["sample_size"]))

# Precompute once
FOSSIL_SET = set(CONFIG['fossil_resources'])
def classify(name):
    return 'fossil' if name in FOSSIL_SET else 'metal'


# Function to get "extractive" intermediate flows, elementary flows, compute supply risk for selected activity
def compute_lca_flows(db_name, activities):
    records = []
    for act in activities:
        func_str = f"{act['name']} | {act.get('reference product','')} | {act['location']}"
        lca = LCA({act.key: 1})
        lca.lci()

        # --- Intermediate flows ---
        rows_i = [
            {
                'Activity':      meta['name'],
                'Product':       meta.get('reference product',''),
                'Geography':     meta['location'],
                'Supply_Amount': amt,
            }
            for idx, key in enumerate(lca.activity_dict)
            if (amt := lca.supply_array[idx]) > 0
            and key in activity_meta[db_name]
            for meta in [activity_meta[db_name][key]]
        ]
        df_i = (
            pd.DataFrame(rows_i)
              .merge(supply_risk_factors_int, on=['Activity','Product','Geography'], how='inner')
              .assign(
                  risky_int_mass=lambda df: df['Supply_Amount'].where(df['Supplyrisk'] > 0, 0),
                  risk_int=lambda df: df['Supplyrisk'] * df['Supply_Amount'],
                  Category=lambda df: df['Product'].map(classify),
              )
        )

        # --- Elementary flows ---
        sums = lca.inventory.sum(axis=1).A1
        rows_e = [
            {'Substance': bm['name'], 'Mass': mass}
            for idx, key in enumerate(lca.biosphere_dict)
            if (mass := sums[idx]) > 0
            and key in bios_meta
            and (bm := bios_meta[key])['type'] == 'natural resource'
        ]
        df_e = (
            pd.DataFrame(rows_e)
              .merge(supply_risk_factors_elem, on='Substance', how='inner')
              .assign(
                  risky_elem_mass=lambda df: df['Mass'].where(df['Supplyrisk'] > 0, 0),
                  risk_elem=lambda df: df['Supplyrisk'] * df['Mass'],
                  Category=lambda df: df['Substance'].map(classify),
              )
        )

        records.append((func_str, df_i, df_e))

    return records


# Function to summarize
def summarize_records(records):
    summary = []
    for func_str, df_i, df_e in records:
        for cat in ['fossil', 'metal']:
            int_mass        = df_i.loc[df_i.Category == cat, 'Supply_Amount'].sum()
            risky_int_mass  = df_i.loc[df_i.Category == cat, 'risky_int_mass'].sum()
            risk_int        = df_i.loc[df_i.Category == cat, 'risk_int'].sum()

            elem_mass       = df_e.loc[df_e.Category == cat, 'Mass'].sum()
            risky_elem_mass = df_e.loc[df_e.Category == cat, 'risky_elem_mass'].sum()
            risk_elem       = df_e.loc[df_e.Category == cat, 'risk_elem'].sum()

            summary.append({
                'Activity':         func_str,
                'Category':         cat,
                'int_mass':         int_mass,
                'risky_int_mass':   risky_int_mass,
                'risk_int':         risk_int,
                'elem_mass':        elem_mass,
                'risky_elem_mass':  risky_elem_mass,
                'risk_elem':        risk_elem,
                'Activity Label':   extract_label(func_str),
            })
    return pd.DataFrame(summary)

# Extracts a concise label from the activity string via commas (for plotting)
def extract_label(activity_name):
    parts = activity_name.split(',')
    if len(parts) >= 4:
        return f"{parts[1].strip()} - {parts[2].strip()}"
    elif len(parts) >= 2:
        return parts[1].strip()
    else:
        return activity_name

# Function to plot results with intermediate VS elementary
def plot_scatter(df, x_column, y_column, xlabel, ylabel, title, offset=0.02, dpi=300):
    # Compute dynamic limits
    x_min, x_max = df[x_column].min(), df[x_column].max()
    y_min, y_max = df[y_column].min(), df[y_column].max()
    x_range = max(x_max - x_min, 1e-5)
    y_range = max(y_max - y_min, 1e-5)
    xlim = (max(x_min - 0.1 * x_range, 0), x_max + 0.1 * x_range)
    ylim = (max(y_min - 0.1 * y_range, 0), y_max + 0.1 * y_range)
    plt.figure(figsize=(7, 5), dpi=dpi)
    ref_line = np.linspace(min(xlim[0], ylim[0]), max(xlim[1], ylim[1]), 100)
    plt.plot(ref_line, ref_line, 'k--', label="y = x (Reference)")
    plt.scatter(df[x_column], df[y_column], marker='x', color='black')
    labels = df.get('Activity Label', df['Activity'])
    for i, label in enumerate(labels):
        plt.annotate(
            label,
            (df[x_column].iloc[i] + offset * x_range, df[y_column].iloc[i]),
            fontsize=10
        )
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.title(title, fontsize=14)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.xticks(np.round(plt.xticks()[0], 2), fontsize=12)
    plt.yticks(np.round(plt.yticks()[0], 2), fontsize=12)
    plt.show()


def plot_activity_contributions(records, top_n=5, max_activities=None, dpi=300):


    def ideal_text_color(hexcolor):
        rgb = np.array(mcolors.to_rgb(hexcolor))
        lum = 0.299*rgb[0] + 0.587*rgb[1] + 0.114*rgb[2]
        return 'white' if lum < 0.6 else 'black'

    ylabels = {
        'mass':       "Mass (kg)",
        'risky_mass': "Risky mass (kg)",
        'risk':       "Supply risk (kg Cu eq)"
    }

    # if max_activities is set, only plot that many
    to_plot = records if max_activities is None else records[:max_activities]

    for func_str, df_i, df_e in to_plot:
        # split title into two lines at first '|'
        title_lines = func_str.split('|', 1)
        big_title   = "\n".join([s.strip() for s in title_lines])

        # taller figure for top space
        fig, axes = plt.subplots(3, 4, figsize=(10, 16), sharey='row', dpi=dpi)
        fig.suptitle(big_title, fontsize=12, y=0.94)

        def plot_stack(ax, data, value_col, color_col, title):
            data = data.copy()
            data[color_col] = data[color_col].fillna('#D3D3D3')
            d = data.sort_values(value_col, ascending=False).head(top_n).copy()
            rest = data[value_col].sum() - d[value_col].sum()
            if rest > 0:
                d = pd.concat([
                    d,
                    pd.DataFrame({
                        data.columns[0]: ["Rest"],
                        value_col:       [rest],
                        color_col:       ["#D3D3D3"]
                    })
                ], ignore_index=True)

            bottom = 0
            max_val = data[value_col].max()
            ax.set_title(title, fontsize=9)
            ax.set_xticks([])

            for _, row in d.iterrows():
                h = row[value_col]
                bar = ax.bar(0, h, bottom=bottom,
                             color=row[color_col], width=0.6)[0]
                y_center = bottom + h/2
                label = str(row[data.columns[0]])
                txt_color = ideal_text_color(row[color_col])

                if h >= max_val * 0.05:
                    ax.text(0, y_center, label,
                            va='center', ha='center',
                            fontsize=7, color=txt_color)
                else:
                    x_center = bar.get_x() + bar.get_width()/2
                    y_start  = bottom + h
                    y_end    = y_start + max_val*0.05
                    ax.plot([x_center, x_center], [y_start, y_end],
                            'k-', linewidth=0.5)
                    ax.text(x_center, y_end + max_val*0.01, label,
                            va='bottom', ha='center',
                            fontsize=7, color='black')
                bottom += h

        panels = [
            ('metal',   'intermediate', df_i),
            ('metal',   'elementary',   df_e),
            ('fossil',  'intermediate', df_i),
            ('fossil',  'elementary',   df_e),
        ]
        colmap = {
            'mass':        {'intermediate':'Supply_Amount',   'elementary':'Mass'},
            'risky_mass':  {'intermediate':'risky_int_mass',  'elementary':'risky_elem_mass'},
            'risk':        {'intermediate':'risk_int',        'elementary':'risk_elem'},
        }
        metrics = ['mass', 'risky_mass', 'risk']

        for i, metric in enumerate(metrics):
            for j, (cat, flow, df) in enumerate(panels):
                ax = axes[i, j]
                sub = df[df.Category == cat]
                group_col = 'Product' if flow=='intermediate' else 'Substance'
                grouped = (
                    sub.groupby(group_col, as_index=False)
                       .agg({colmap[metric][flow]:'sum', 'Color':'first'})
                )
                title = f"{cat.title()} – {flow.title()}"
                plot_stack(
                    ax,
                    grouped.rename(columns={group_col:grouped.columns[0]}),
                    colmap[metric][flow],
                    'Color',
                    title
                )
                if j == 0:
                    ax.set_ylabel(ylabels[metric], fontsize=9)

        plt.subplots_adjust(
            top=0.88,
            bottom=0.06,
            left=0.04,
            right=0.96,
            hspace=0.4,
            wspace=0.2
        )
        plt.show()




# Execution
if __name__ == "__main__":

    activities = get_filtered_activities(CONFIG['databases'][0])
    all_results = {}
    
    for db in CONFIG['databases']:
        logger.info(f"Processing DB: {db}")
        # compute and summarize
        rec        = compute_lca_flows(db, activities)
        summary_df = summarize_records(rec)
    
        # split into fossil/metalframes in one shot
        all_results[db] = {
            cat: summary_df[summary_df['Category'] == cat]
            for cat in ['fossil', 'metal']
        }
    
    # Plot each category
    for db, cats in all_results.items():
        for cat, df in cats.items():
            plot_scatter(
                df,
                'risky_int_mass', 'risky_elem_mass',
                'mass - intermediate (kg/kg)',
                'mass - elementary (kg/kg)',
                f'Mass (with risk) - {cat} – {db}'
            )
            plot_scatter(
                df,
                'risk_int', 'risk_elem',
                'risk - intermediate (kg Cu eq/kg)',
                'risk - elementary (kg Cu eq/kg)',
                f'Supply risk - {cat} – {db}'
            )
            
    plot_activity_contributions(rec, top_n=5, max_activities=5, dpi=300)
