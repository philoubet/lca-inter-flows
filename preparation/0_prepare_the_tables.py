import pandas as pd
import pubchempy as pcp
import periodictable
import chemparse
import re
import os

# 📌 STEP 1: Load "Elementary Exchanges" from ecoinvent database
excel_file = "Database-Overview-for-ecoinvent-v3.10.xlsx"
df_eco = pd.read_excel(excel_file, sheet_name="Elementary Exchanges")

# Filter for natural resources in ground
df_flows = df_eco[
    (df_eco["Compartment"] == "natural resource") & 
    (df_eco["Subcompartment"] == "in ground")
][["Name", "CAS Number", "Formula"]]

print(f"✅ Extracted {len(df_flows)} natural resources from 'Elementary Exchanges'.")

# 📌 STEP 2: Load "B_public.csv"
csv_b_public = "universal_matrix_export_3.10_cut-off/B_public.csv"
df_b_public = pd.read_csv(csv_b_public, delimiter=";")

# Keep relevant columns and filter system process = 2
df_b_public = df_b_public[df_b_public["system process"] == 2][["activityName", "product", "geography"]]

print(f"✅ Extracted {len(df_b_public)} products from B_public.csv (filtered on system process = 2).")

cache_file = "pubchem_cache.csv"

if os.path.exists(cache_file):
    pubchem_cache = pd.read_csv(cache_file, delimiter=";").set_index("Product")["Formula"].to_dict()
    print(f"✅ Loaded {len(pubchem_cache)} cached PubChem formulas.")
else:
    pubchem_cache = {}

# 📌 STEP 3: Load manual matches from "intermediate_flows_with_colors_corrected.csv"
csv_matches = "intermediate_flows_with_colors_corrected.csv"
df_matches = pd.read_csv(csv_matches, delimiter=";")[["Product", "Matched Substance", "Color"]].dropna()

# Create mapping dictionaries
ecoinvent_dict = {row["Name"]: row["Formula"] for _, row in df_flows.iterrows()}
manual_matches = {row["Product"]: row["Matched Substance"] for _, row in df_matches.iterrows()}
color_dict = {row["Product"]: row["Color"] for _, row in df_matches.iterrows()}  # Color mapping

# 📌 STEP 4: Extract chemical component from product names
def extract_chemical_part(product_name):
    """Extract the most likely chemical component from a product name."""
    if product_name in manual_matches:
        return manual_matches[product_name]

    match = re.search(r"(\d+(\.\d+)?)\s?%\s?([\w\s-]+)", product_name)
    if match:
        return match.group(3).strip()

    return product_name  # Default to full product name

df_b_public["Matched Substance"] = df_b_public["product"].apply(extract_chemical_part)

# 📌 STEP 5: Assign color from df_matches
default_color = ""  # Leave empty if no color found
df_b_public["Color"] = df_b_public["product"].map(color_dict).fillna(default_color)

# 📌 STEP 5: Fetch formula (ecoinvent first, then PubChem)
def get_product_formula(chemical_name):
    """Get the chemical formula from ecoinvent or PubChem."""
    if chemical_name in ecoinvent_dict:
        return ecoinvent_dict[chemical_name]

    try:
        compound = pcp.get_compounds(chemical_name, 'name')
        if compound:
            return compound[0].molecular_formula  # Return PubChem formula
    except Exception as e:
        print(f"❌ Error fetching formula for {chemical_name}: {e}")

    return None  # Not found

df_b_public["Formula"] = df_b_public["Matched Substance"].apply(get_product_formula)

# 📌 STEP 6: Extract percentage from product names
def extract_percentage(product_name):
    """Extract percentage factor from the product name (e.g., '95% titanium dioxide')."""
    match = re.search(r"(\d+(\.\d+)?)\s?%", product_name)
    return float(match.group(1)) / 100 if match else 1  # Default to 100%

df_b_public["Product_Percentage"] = df_b_public["product"].apply(extract_percentage)

# 📌 STEP 7: Compute molar mass
def molar_mass(formula):
    """Calculate the molar mass of a compound."""
    try:
        parsed_formula = chemparse.parse_formula(formula)
        total_mass = sum(periodictable.elements.symbol(element).mass * count for element, count in parsed_formula.items())
        return total_mass, parsed_formula
    except Exception as e:
        print(f"❌ Error parsing formula {formula}: {e}")
        return None, None

# 📌 STEP 8: Compute atomic percentage
element_name_to_symbol = {el.name.capitalize(): el.symbol for el in periodictable.elements}

def compute_element_percentage(product, substance):
    """Computes the percentage of an element in a product using the chemical formula."""
    formula = pubchem_cache.get(product)  # Use cached formulas
    if not formula:
        return None  # If no formula found

    total_mass, parsed_formula = molar_mass(formula)
    if total_mass is None or parsed_formula is None:
        return None  # Return None if molar mass calculation failed

    # Convert substance name to atomic symbol
    substance_symbol = element_name_to_symbol.get(substance.capitalize(), None)
    if not substance_symbol or substance_symbol not in parsed_formula:
        return None  # No unnecessary print statements

    element_mass = periodictable.elements.symbol(substance_symbol).mass * parsed_formula[substance_symbol]
    return (element_mass / total_mass) * 100  # Percentage

# Apply atomic percentage calculation
df_b_public["Atomic_Percentage"] = df_b_public.apply(
    lambda row: compute_element_percentage(row["product"], row["Matched Substance"]) if pd.notna(row["Matched Substance"]) else None,
    axis=1
)

# Ensure no NaN values remain in "Atomic_Percentage"
df_b_public["Atomic_Percentage"].fillna(100, inplace=True)  # Just in case

# 📌 STEP 9: Apply percentage factor
df_b_public["Final_Percentage"] = df_b_public["Atomic_Percentage"] * df_b_public["Product_Percentage"]

df_b_public.rename(columns={"activityName": "Activity", "product": "Product", "geography":"Geography"}, inplace=True)

# 📌 STEP 10: Save and display results
df_b_public.to_csv("processed_b_public_with_percentages.csv", sep=";", index=False)

print("🎉 Process completed! You can now review the results.")
