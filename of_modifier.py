import csv
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

# ==============================================================================
# 1. DICTIONARY PARSING & REGEX HELPERS
# ==============================================================================

def get_block_bounds(text, block_name):
    """
    Finds block boundaries robustly by anchoring search to line starts and
    jumping over inline OpenFOAM comments before the opening brace.
    """
    pattern = r'(?m)^[ \t]*' + re.escape(block_name) + r'\b'
    match = re.search(pattern, text)
    if not match:
        return None, None

    start_idx = match.start()
    open_brace_idx = text.find('{', match.end())
    if open_brace_idx == -1:
        return None, None

    brace_count = 0
    for i in range(open_brace_idx, len(text)):
        if text[i] == '{':
            brace_count += 1
        elif text[i] == '}':
            brace_count -= 1
            if brace_count == 0:
                return start_idx, i + 1
    return None, None


def remove_table_entry(text, keyword):
    """Removes existing Function1 table entries from dictionary blocks robustly."""
    pattern = r'\b' + re.escape(keyword) + r'\b'
    while True:
        match = re.search(pattern, text)
        if not match:
            break
        start_idx = match.start()
        
        open_idx = -1
        for i in range(match.end(), len(text)):
            if text[i] in ['{', '(']:
                open_idx = i
                break
        
        if open_idx == -1:
            break
            
        open_char = text[open_idx]
        close_char = '}' if open_char == '{' else ')'
        
        brace_count = 0
        end_idx = -1
        for i in range(open_idx, len(text)):
            if text[i] == open_char:
                brace_count += 1
            elif text[i] == close_char:
                brace_count -= 1
                if brace_count == 0:
                    end_idx = i + 1
                    break
        
        if end_idx != -1:
            for i in range(end_idx, len(text)):
                if text[i] == ';':
                    end_idx = i + 1
                    break
                elif text[i] not in [' ', '\t', '\n', '\r']:
                    break
            text = text[:start_idx] + text[end_idx:]
        else:
            break
    return text


def build_foam_table(table_name, x_vals, y_vals, indent="        "):
    """Formats 1D arrays into strictly unambiguous OpenFOAM Function1 dictionary syntax."""
    lines = [
        f"{indent}{table_name}",
        f"{indent}{{",
        f"{indent}    type            table;",
        f"{indent}    outOfBounds     clamp;",
        f"{indent}    values",
        f"{indent}    ("
    ]
    for x, y in zip(x_vals, y_vals):
        lines.append(f"{indent}        ({x} {y})")
    lines.append(f"{indent}    );")  # <--- FIXED: Closes the values list correctly
    lines.append(f"{indent}}}")
    return "\n".join(lines)


# ==============================================================================
# 2. SPLINE GENERATION & TABULAR DATA
# ==============================================================================

def generate_spline_tables(input_csv, output_csv=None, n_points=100, plot=False):
    """
    Fits monotone PCHIP splines to raw data and samples uniform grid points.
    Returns pandas DataFrame and optionally writes to output_csv.
    """
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Raw data file '{input_csv}' not found.")

    df_raw = pd.read_csv(input_csv)
    required_cols = ['s', 'krL', 'krG', 'pc']
    for col in required_cols:
        if col not in df_raw.columns:
            raise KeyError(f"Missing column '{col}' in {input_csv}")

    df_raw = df_raw.sort_values(by='s').drop_duplicates(subset=['s'])

    s_raw = df_raw['s'].values
    krl_raw = df_raw['krL'].values
    krg_raw = df_raw['krG'].values
    pc_raw = df_raw['pc'].values

    # Fit PCHIP splines
    krl_spline = PchipInterpolator(s_raw, krl_raw)
    krg_spline = PchipInterpolator(s_raw, krg_raw)
    pc_spline = PchipInterpolator(s_raw, pc_raw)
    dpcds_spline = pc_spline.derivative()  

    s_grid = np.linspace(0.0, 1.0, n_points)
    krl_grid = np.clip(krl_spline(s_grid), 0.0, 1.0)
    krg_grid = np.clip(krg_spline(s_grid), 0.0, 1.0)
    pc_grid = pc_spline(s_grid)
    dpcds_grid = np.clip(dpcds_spline(s_grid), 1e-3, 1e5)  

    out_df = pd.DataFrame({
        's': s_grid,
        'krL': krl_grid,
        'krG': krg_grid,
        'pc': pc_grid,
        'dpcds': dpcds_grid
    })

    if output_csv:
        out_df.to_csv(output_csv, index=False, float_format='%.6e')
        print(f"Generated {n_points} points and saved to '{output_csv}'.")

    if plot:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        ax1.plot(s_raw, pc_raw, 'ko', label='Raw $P_c$')
        ax1.plot(s_grid, pc_grid, 'b-', label='PCHIP Fit $P_c$')
        ax1.set_xlabel('Liquid Saturation $S_l$')
        ax1.set_ylabel('Capillary Pressure $P_c$ [Pa]', color='b')
        ax1.grid(True)

        ax1_twin = ax1.twinx()
        ax1_twin.plot(s_grid, dpcds_grid, 'r--', label='Derivative $dP_c/dS_l$')
        ax1_twin.set_ylabel('Derivative [Pa]', color='r')
        ax1.set_title('Capillary Pressure & Derivative')

        ax2.plot(s_raw, krl_raw, 'go', label='Raw $K_{r,l}$')
        ax2.plot(s_raw, krg_raw, 'mo', label='Raw $K_{r,g}$')
        ax2.plot(s_grid, krl_grid, 'g-', label='PCHIP $K_{r,l}$')
        ax2.plot(s_grid, krg_grid, 'm-', label='PCHIP $K_{r,g}$')
        ax2.set_xlabel('Liquid Saturation $S_l$')
        ax2.set_ylabel('Relative Permeability $K_r$')
        ax2.set_title('Relative Permeabilities')
        ax2.set_ylim(-0.05, 1.05)
        ax2.grid(True)
        ax2.legend()

        plt.tight_layout()
        plt.show()

    return out_df


# ==============================================================================
# 3. OPENFOAM DICTIONARY MODIFIERS
# ==============================================================================

def update_porous_zones_tables(file_path, tables_data, target_zones):
    """Inserts generated spline tables into porousZones dictionary."""
    if isinstance(tables_data, str):
        df = pd.read_csv(tables_data)
    else:
        df = tables_data

    s_vals = df['s'].values
    krL_vals = df['krL'].values
    krG_vals = df['krG'].values
    pc_vals = df['pc'].values
    dpcds_vals = df['dpcds'].values

    krL_str = build_foam_table("krLTable", s_vals, krL_vals)
    krG_str = build_foam_table("krGTable", s_vals, krG_vals)
    pc_str = build_foam_table("pcTable", s_vals, pc_vals)
    dpcds_str = build_foam_table("dpcdsTable", s_vals, dpcds_vals)

    new_tables_block = f"{krL_str}\n\n{krG_str}\n\n{pc_str}\n\n{dpcds_str}"

    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read().replace('\r\n', '\n')

    for zone in target_zones:
        start_idx, end_idx = get_block_bounds(content, zone)
        if start_idx is None:
            continue

        zone_block = content[start_idx:end_idx]
        sub_start, sub_end = get_block_bounds(zone_block, "DarcyForchheimerCoeffs")

        if sub_start is not None:
            sub_block = zone_block[sub_start:sub_end]
            for tbl in ["krLTable", "krGTable", "pcTable", "dpcdsTable"]:
                sub_block = remove_table_entry(sub_block, tbl)

            coord_idx = sub_block.find("coordinateSystem")
            cap_idx = sub_block.find("capillaryPressureModel")

            if cap_idx != -1:
                semi_idx = sub_block.find(';', cap_idx)
                insert_pos = (semi_idx + 1) if semi_idx != -1 else cap_idx
                updated_sub_block = (
                    sub_block[:insert_pos] + "\n\n" + new_tables_block + "\n" + sub_block[insert_pos:].lstrip('\n')
                )
            elif coord_idx != -1:
                insert_pos = sub_block.rfind('\n', 0, coord_idx)
                insert_pos = (insert_pos + 1) if insert_pos != -1 else coord_idx
                updated_sub_block = sub_block[:insert_pos] + new_tables_block + "\n\n" + sub_block[insert_pos:]
            else:
                last_brace_pos = sub_block.rfind('}')
                updated_sub_block = sub_block[:last_brace_pos].rstrip() + "\n\n" + new_tables_block + "\n    }"

            updated_zone_block = zone_block[:sub_start] + updated_sub_block + zone_block[sub_end:]
        else:
            for tbl in ["krLTable", "krGTable", "pcTable", "dpcdsTable"]:
                zone_block = remove_table_entry(zone_block, tbl)
            last_brace_pos = zone_block.rfind('}')
            updated_zone_block = zone_block[:last_brace_pos].rstrip() + "\n\n" + new_tables_block + "\n}"

        content = content[:start_idx] + updated_zone_block + content[end_idx:]

    with open(file_path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(content)

    print(f"Successfully updated spline tables in '{file_path}' for: {target_zones}")


def update_porosity_and_d(file_path, update_data):
    """Updates root-level porosity and Darcy coefficients 'd' in porousZones."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read().replace('\r\n', '\n')

    for zone, data in update_data.items():
        start_idx, end_idx = get_block_bounds(content, zone)
        if start_idx is None:
            print(f"Warning: Zone '{zone}' not found in {file_path}. Skipping.")
            continue

        zone_block = content[start_idx:end_idx]

        if 'porosity' in data:
            new_porosity = data['porosity']
            porosity_pattern = r'(?m)^([ \t]*porosity[ \t]+)[^;]+;'
            zone_block = re.sub(porosity_pattern, rf'\g<1>{new_porosity};', zone_block)

        if 'd' in data:
            sub_start, sub_end = get_block_bounds(zone_block, "DarcyForchheimerCoeffs")
            if sub_start is not None:
                sub_block = zone_block[sub_start:sub_end]
                new_d = data['d']
                d_pattern = r'(?m)^([ \t]*d[ \t]+)\([^)]+\)[ \t]*;'
                sub_block = re.sub(d_pattern, rf'\g<1>{new_d};', sub_block)
                zone_block = zone_block[:sub_start] + sub_block + zone_block[sub_end:]

        content = content[:start_idx] + zone_block + content[end_idx:]

    with open(file_path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(content)

    print(f"Successfully updated porosity/d in '{file_path}' for zones: {list(update_data.keys())}")


def update_diffusivity_model(file_path, update_data):
    """Updates porosity and dPore in diffusivityModel.air."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read().replace('\r\n', '\n')

    for zone, data in update_data.items():
        start_idx, end_idx = get_block_bounds(content, zone)
        if start_idx is None:
            print(f"Warning: Zone '{zone}' not found in {file_path}. Skipping.")
            continue

        zone_block = content[start_idx:end_idx]
        sub_start, sub_end = get_block_bounds(zone_block, "porousFSGCoeffs")

        if sub_start is not None:
            sub_block = zone_block[sub_start:sub_end]

            if 'porosity' in data:
                new_porosity = data['porosity']
                porosity_pattern = r'(?m)^([ \t]*porosity[ \t]+)[^;]+;'
                sub_block = re.sub(porosity_pattern, rf'\g<1>{new_porosity};', sub_block)

            if 'dPore' in data:
                new_dpore = data['dPore']
                dpore_pattern = r'(?m)^([ \t]*dPore[ \t]+dPore[ \t]+\[.*?\][ \t]+)[^;]+;'
                sub_block = re.sub(dpore_pattern, rf'\g<1>{new_dpore};', sub_block)

            zone_block = zone_block[:sub_start] + sub_block + zone_block[sub_end:]

        content = content[:start_idx] + zone_block + content[end_idx:]

    with open(file_path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(content)

    print(f"Successfully updated diffusivity in '{file_path}' for zones: {list(update_data.keys())}")


# ==============================================================================
# 4. PARAMETER CSV PARSER & HIGH-LEVEL PIPELINE
# ==============================================================================

def load_parameters_from_csv(csv_path):
    """
    Parses a CSV containing zone parameter values.
    Expected CSV format:
    zone,porosity,dPore,d
    cathode,0.75,20.0e-06,(2e11 2e11 2e11)
    mpl,0.45,50.0e-06,(2e12 2e12 2e12)
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Parameter CSV file '{csv_path}' not found.")

    param_df = pd.read_csv(csv_path)
    update_data = {}

    for _, row in param_df.iterrows():
        zone = str(row['zone']).strip()
        zone_data = {}

        if 'porosity' in row and pd.notna(row['porosity']):
            zone_data['porosity'] = str(row['porosity']).strip()
        if 'dPore' in row and pd.notna(row['dPore']):
            zone_data['dPore'] = str(row['dPore']).strip()
        if 'd' in row and pd.notna(row['d']):
            zone_data['d'] = str(row['d']).strip()

        if zone_data:
            update_data[zone] = zone_data

    return update_data


def update_all_openfuelcell_files(case_dir, raw_data_csv, params_csv=None, param_dict=None, target_zones=["cathode", "mpl"]):
    """
    High-level orchestrator function to update all OpenFOAM files in a single pass.
    Parameters can be passed either via `params_csv` OR directly via `param_dict`.
    """
    porous_zones_path = os.path.join(case_dir, "constant/air/porousZones")
    diffusivity_path = os.path.join(case_dir, "constant/air/diffusivityModel.air")

    # 1. Generate spline tables and update porousZones
    tables_df = generate_spline_tables(raw_data_csv, plot=False)
    update_porous_zones_tables(porous_zones_path, tables_df, target_zones)

    # 2. Resolve parameter dictionary
    if params_csv:
        update_data = load_parameters_from_csv(params_csv)
    elif param_dict:
        update_data = param_dict
    else:
        print("No parameters provided for porosity/dPore/d. Skipping scalar updates.")
        return

    # 3. Apply scalar updates to OpenFOAM dictionaries
    update_porosity_and_d(porous_zones_path, update_data)
    update_diffusivity_model(diffusivity_path, update_data)


# ==============================================================================
# 5. CLI EXECUTION ENTRY POINT
# ==============================================================================

if __name__ == "__main__":
    # Example usage for testing standalone script
    CASE_DIRECTORY = "run/PEMFC_capillaryPressure"
    RAW_DATA_CSV = "raw_data.csv"
    PARAMS_CSV = "params.csv"

    # Ensure demo parameters CSV exists if running directly
    if not os.path.exists(PARAMS_CSV):
        demo_params = pd.DataFrame([
            {"zone": "cathode", "porosity": "0.75", "dPore": "20.0e-06", "d": "(2e11 2e11 2e11)"},
            {"zone": "mpl", "porosity": "0.45", "dPore": "50.0e-06", "d": "(2e12 2e12 2e12)"}
        ])
        demo_params.to_csv(PARAMS_CSV, index=False)
        print(f"Created template parameter file '{PARAMS_CSV}'.")

    update_all_openfuelcell_files(
        case_dir=CASE_DIRECTORY,
        raw_data_csv=RAW_DATA_CSV,
        params_csv=PARAMS_CSV,
        target_zones=["cathode", "mpl"]
    )