import tables             # Library to read Nastran HDF5 binary database files
import pandas as pd       # Library for data manipulation, table creation, and Excel export
import numpy as np        # Library for high-performance vectorized mathematical operations
import os                 # Library to interact with the Operating System
import glob               # Library to search for .h5 files
import re                 # Library for Regular Expressions (text cleaning)
from tqdm import tqdm     # Library to display progress bars in terminal
from xlsxwriter.utility import xl_col_to_name # Library to format Excel columns
import matplotlib.pyplot as plt # Standard library for generating 2D graphs
import seaborn as sns     # Advanced visualization library for better-looking charts
import warnings
import gc                 # Garbage collector for RAM management
warnings.filterwarnings("ignore") # Suppress unnecessary matplotlib warnings

# ==============================================================================
# 1. STRING CLEANER & UTILITIES
# ==============================================================================
def clean_subtitle(text):
    """ Cleans byte strings from Nastran H5 into normal strings """
    if isinstance(text, bytes): 
        text = text.decode('utf-8', errors='ignore')
    text = str(text).strip()
    if (text.startswith("b'") and text.endswith("'")) or (text.startswith('b"') and text.endswith('"')):
        text = text[2:-1]
    text = text.replace("'", "").replace('"', '')
    return text.strip()

def safe_filename(text):
    """ Replaces illegal characters in a string for saving image files """
    return re.sub(r'[^A-Za-z0-9_\-\.]', '_', text)

# ==============================================================================
# 2. COMPOSITE FAILURE PHYSICS ENGINE (ALL CRITERIA)
# ==============================================================================
def analyze_composite_failure(df, props):
    """
    Core mathematical engine to calculate Failure Indices (FI) and Margins of Safety (MS).
    Includes: Max Stress, Max Strain, Tsai-Hill, Tsai-Wu, Hashin, and Interlaminar Shear.
    """
    E1, E2, nu12, G12 = props['E1'], props['E2'], props['nu12'], props['G12']
    
    # SAFEGUARD: Prevent ZeroDivisionError for strengths
    def enforce_limit(val): return abs(val) if abs(val) > 1e-6 else 1e12 
        
    Xt, Xc = enforce_limit(props['Xt']), enforce_limit(props['Xc'])
    Yt, Yc = enforce_limit(props['Yt']), enforce_limit(props['Yc'])
    S = enforce_limit(props['S'])
    S13 = enforce_limit(props.get('S13', S)) # Interlaminar Shear XZ (Default to S if missing)
    S23 = enforce_limit(props.get('S23', S)) # Interlaminar Shear YZ (Default to S if missing)

    # Extract numpy arrays for extreme calculation speed
    s1, s2, t12 = df['Sigma_11'].values, df['Sigma_22'].values, df['Tau_12'].values
    t13, t23 = df.get('Tau_13', np.zeros_like(s1)).values, df.get('Tau_23', np.zeros_like(s1)).values

    # --- A. MAXIMUM STRESS CRITERION ---
    fi_1_maxs = np.where(s1 > 0, s1 / Xt, np.abs(s1) / Xc)
    fi_2_maxs = np.where(s2 > 0, s2 / Yt, np.abs(s2) / Yc)
    fi_12_maxs = np.abs(t12) / S 
    fi_max_stress = np.maximum.reduce([fi_1_maxs, fi_2_maxs, fi_12_maxs])
    ms_max_stress = np.where(fi_max_stress > 0, (1.0 / fi_max_stress) - 1.0, 999.9)

    modes_ms = []
    for f1, f2, fm, sig1, sig2 in np.nditer([fi_1_maxs, fi_2_maxs, fi_max_stress, s1, s2]):
        if fm == 0: modes_ms.append("Safe")
        elif fm == f1 and sig1 > 0: modes_ms.append("Fiber Tension")
        elif fm == f1 and sig1 <= 0: modes_ms.append("Fiber Compression")
        elif fm == f2 and sig2 > 0: modes_ms.append("Matrix Tension")
        elif fm == f2 and sig2 <= 0: modes_ms.append("Matrix Compress")
        else: modes_ms.append("In-Plane Shear")

    # --- B. MAXIMUM STRAIN CRITERION ---
    eps_1T, eps_1C = Xt / E1, Xc / E1
    eps_2T, eps_2C = Yt / E2, Yc / E2
    gamma_12u = S / G12
    nu21 = (nu12 * E2) / E1 if E1 > 1e-6 else 0.0

    eps1_act = (s1 / E1) - (nu21 * s2 / E2)
    eps2_act = -(nu12 * s1 / E1) + (s2 / E2)
    gamma12_act = t12 / G12

    fi_1_maxe = np.where(eps1_act > 0, eps1_act / eps_1T, np.abs(eps1_act) / eps_1C)
    fi_2_maxe = np.where(eps2_act > 0, eps2_act / eps_2T, np.abs(eps2_act) / eps_2C)
    fi_12_maxe = np.abs(gamma12_act) / gamma_12u

    fi_max_strain = np.maximum.reduce([fi_1_maxe, fi_2_maxe, fi_12_maxe])
    ms_max_strain = np.where(fi_max_strain > 0, (1.0 / fi_max_strain) - 1.0, 999.9)

    # --- C. TSAI-HILL CRITERION ---
    X_hill = np.where(s1 > 0, Xt, Xc)
    Y_hill = np.where(s2 > 0, Yt, Yc)
    fi_tsai_hill = (s1 / X_hill)**2 - (s1 * s2) / (X_hill**2) + (s2 / Y_hill)**2 + (t12 / S)**2
    ms_tsai_hill = np.where(fi_tsai_hill > 0, (1.0 / np.sqrt(fi_tsai_hill)) - 1.0, 999.9)

    # --- D. TSAI-WU CRITERION ---
    F1 = (1.0 / Xt) - (1.0 / Xc)
    F2 = (1.0 / Yt) - (1.0 / Yc)
    F11 = 1.0 / (Xt * Xc)
    F22 = 1.0 / (Yt * Yc)
    F66 = 1.0 / (S**2)
    F12 = -0.5 * np.sqrt(F11 * F22) 

    fi_tsai_wu = F1*s1 + F2*s2 + F11*(s1**2) + F22*(s2**2) + F66*(t12**2) + 2*F12*s1*s2

    A_tw = F11*(s1**2) + F22*(s2**2) + F66*(t12**2) + 2*F12*s1*s2
    B_tw = F1*s1 + F2*s2
    
    R_tw = np.zeros_like(s1, dtype=float)
    mask_quad = A_tw > 1e-12 
    R_tw[mask_quad] = (-B_tw[mask_quad] + np.sqrt(B_tw[mask_quad]**2 + 4 * A_tw[mask_quad])) / (2 * A_tw[mask_quad])
    mask_lin = (A_tw <= 1e-12) & (B_tw > 1e-12)
    R_tw[mask_lin] = 1.0 / B_tw[mask_lin]
    R_tw[(A_tw <= 1e-12) & (B_tw <= 1e-12)] = 1000.0
    ms_tsai_wu = R_tw - 1.0

    # --- E. HASHIN CRITERION (1980 2D) ---
    fi_hashin_ft = np.where(s1 > 0, (s1/Xt)**2 + (t12/S)**2, 0)
    fi_hashin_fc = np.where(s1 <= 0, (np.abs(s1)/Xc)**2, 0)
    fi_hashin_mt = np.where(s2 > 0, (s2/Yt)**2 + (t12/S)**2, 0)
    fi_hashin_mc = np.where(s2 <= 0, (np.abs(s2)/Yc)**2 + (t12/S)**2, 0)
    
    fi_hashin = np.maximum.reduce([fi_hashin_ft, fi_hashin_fc, fi_hashin_mt, fi_hashin_mc])
    ms_hashin = np.where(fi_hashin > 0, (1.0 / np.sqrt(fi_hashin)) - 1.0, 999.9)

    hashin_modes = []
    for ft, fc, mt, mc, fm in np.nditer([fi_hashin_ft, fi_hashin_fc, fi_hashin_mt, fi_hashin_mc, fi_hashin]):
        if fm == 0: hashin_modes.append("Safe")
        elif fm == ft: hashin_modes.append("Fiber Tension")
        elif fm == fc: hashin_modes.append("Fiber Compress")
        elif fm == mt: hashin_modes.append("Matrix Tension")
        else: hashin_modes.append("Matrix Compress")

    # --- F. INTERLAMINAR SHEAR (DELAMINATION) ---
    fi_ils = (t13 / S13)**2 + (t23 / S23)**2
    ms_ils = np.where(fi_ils > 0, (1.0 / np.sqrt(fi_ils)) - 1.0, 999.9)

    # Build the final output DataFrame
    out_df = df.copy()
    out_df['Max_Stress_FI'] = fi_max_stress; out_df['Max_Stress_MS'] = ms_max_stress; out_df['Max_Stress_Mode'] = modes_ms
    out_df['Max_Strain_FI'] = fi_max_strain; out_df['Max_Strain_MS'] = ms_max_strain
    out_df['Tsai_Hill_FI'] = fi_tsai_hill; out_df['Tsai_Hill_MS'] = ms_tsai_hill
    out_df['Tsai_Wu_FI'] = fi_tsai_wu; out_df['Tsai_Wu_MS'] = ms_tsai_wu
    out_df['Hashin_FI'] = fi_hashin; out_df['Hashin_MS'] = ms_hashin; out_df['Hashin_Mode'] = hashin_modes
    out_df['Interlaminar_FI'] = fi_ils; out_df['Interlaminar_MS'] = ms_ils

    return out_df

# ==============================================================================
# 3. VISUALIZATION (TSAI-WU FAILURE ENVELOPES)
# ==============================================================================
def create_failure_envelopes(df_res, props, base_name):
    """ Plots the Tsai-Wu Failure Ellipse and scatters the element stresses """
    print("\n[STEP 6] 🎨 Generating Tsai-Wu Failure Envelopes...")
    sns.set_style("darkgrid")
    
    Xt, Xc = abs(props['Xt']), abs(props['Xc'])
    Yt, Yc = abs(props['Yt']), abs(props['Yc'])
    
    # Tsai-Wu coefficients for plotting
    F1 = (1.0 / Xt) - (1.0 / Xc); F2 = (1.0 / Yt) - (1.0 / Yc)
    F11 = 1.0 / (Xt * Xc); F22 = 1.0 / (Yt * Yc)
    F12 = -0.5 * np.sqrt(F11 * F22) 
    
    # Create meshgrid for the ellipse contour
    s1_vals = np.linspace(-Xc*1.3, Xt*1.3, 400)
    s2_vals = np.linspace(-Yc*1.3, Yt*1.3, 400)
    S1, S2 = np.meshgrid(s1_vals, s2_vals)
    Z = F11*(S1**2) + F22*(S2**2) + F1*S1 + F2*S2 + 2*F12*S1*S2
    
    load_cases = df_res['Load_Case'].unique()
    images_dict = {}

    for lc in load_cases:
        safe_lc = safe_filename(str(lc))
        df_lc = df_res[df_res['Load_Case'] == lc]
        
        plt.figure(figsize=(9, 7))
        plt.contour(S1, S2, Z, levels=[1.0], colors='red', linewidths=2.5, linestyles='dashed')
        
        safe_mask = df_lc['Tsai_Wu_FI'] < 1.0
        failed_mask = df_lc['Tsai_Wu_FI'] >= 1.0
        
        plt.scatter(df_lc[safe_mask]['Sigma_11'], df_lc[safe_mask]['Sigma_22'], 
                    color='green', alpha=0.6, s=25, label='Safe Elements (FI < 1)')
        plt.scatter(df_lc[failed_mask]['Sigma_11'], df_lc[failed_mask]['Sigma_22'], 
                    color='red', marker='X', s=45, label='Failed Elements (FI >= 1)')
        
        plt.axhline(0, color='black', linewidth=1)
        plt.axvline(0, color='black', linewidth=1)
        
        plt.title(f"Tsai-Wu Failure Envelope: {lc}\n(Plotted for Tau_12 = 0)", fontsize=13, fontweight='bold')
        plt.xlabel("Sigma 11 (Fiber Direction) [MPa]", fontweight='bold')
        plt.ylabel("Sigma 22 (Matrix Direction) [MPa]", fontweight='bold')
        plt.legend(loc='upper left', frameon=True, shadow=True)
        plt.tight_layout()
        
        img_name = f"Envelope_{safe_lc}.png"
        plt.savefig(img_name, dpi=150)
        plt.close()
        images_dict[lc] = img_name
        
    return images_dict

# ==============================================================================
# 4. H5 DATA MINING (STRESSES, MAT8, PCOMP, SUBTITLES)
# ==============================================================================
def parse_bdf_subtitles(bdf_filepath):
    subcase_dict = {}
    if not bdf_filepath or not os.path.exists(bdf_filepath): return subcase_dict
    try:
        with open(bdf_filepath, 'r', errors='ignore') as f:
            current_subcase = None
            for line in f:
                line = line.strip()
                if line.startswith("SUBCASE"):
                    try: current_subcase = int(line.split()[1])
                    except: pass
                elif line.startswith("SUBTITLE") and current_subcase is not None:
                    if "=" in line:
                        subtitle_text = line.split("=")[1].strip()
                        subcase_dict[current_subcase] = f"Subcase {current_subcase}: {clean_subtitle(subtitle_text)}"
    except Exception: pass
    return subcase_dict

def extract_mat8_h5(input_file, conv_factor):
    mat_dict = {}
    try:
        h5 = tables.open_file(input_file, mode="r")
        if hasattr(h5.root.NASTRAN.INPUT.MATERIAL, 'MAT8'):
            mat8_table = h5.root.NASTRAN.INPUT.MATERIAL.MAT8.read()
            cols = mat8_table.dtype.names
            for row in mat8_table:
                mid = row['ID'] if 'ID' in cols else (row['MID'] if 'MID' in cols else 1)
                def get_val(col_name): return float(row[col_name]) if col_name in cols else 0.0
                mat_dict[mid] = {
                    'E1': get_val('E1') * conv_factor, 'E2': get_val('E2') * conv_factor,
                    'nu12': get_val('NU12'), 'G12': get_val('G12') * conv_factor,
                    'Xt': get_val('XT') * conv_factor, 'Xc': get_val('XC') * conv_factor,
                    'Yt': get_val('YT') * conv_factor, 'Yc': get_val('YC') * conv_factor,
                    'S': get_val('S') * conv_factor,
                    'S13': get_val('S1Z') * conv_factor if 'S1Z' in cols else get_val('S') * conv_factor,
                    'S23': get_val('S2Z') * conv_factor if 'S2Z' in cols else get_val('S') * conv_factor
                }
        h5.close()
    except Exception: pass
    return mat_dict

def extract_composite_h5(input_file, bdf_file=None):
    """ Main extraction function capturing Stresses, Subtitles, and PCOMP details. """
    print(f"    [+] Extracting Data & PCOMP from H5: {os.path.basename(input_file)} ...")
    h5 = tables.open_file(input_file, mode="r") 
    
    # 1. READ LOAD CASES
    bdf_subtitles = parse_bdf_subtitles(bdf_file)
    domain_to_subcase = {} 
    try:
        if hasattr(h5.root.NASTRAN.RESULT, 'DOMAINS'):
            for row in h5.root.NASTRAN.RESULT.DOMAINS.read():
                d_id = row['ID']
                sub_id = row['SUBCASE'] if 'SUBCASE' in row.dtype.names else d_id
                title_str = ""
                for t_col in ['SUBTITLE', 'TITLE', 'LABEL']:
                    if t_col in row.dtype.names:
                        clean_val = clean_subtitle(row[t_col])
                        if clean_val: title_str = clean_val; break
                domain_to_subcase[d_id] = f"Subcase {sub_id}: {title_str}" if title_str else f"Subcase {sub_id}"
    except Exception: pass

    def get_load_case_name(domain_id): 
        real_subcase = domain_to_subcase.get(domain_id, domain_id)
        return bdf_subtitles.get(real_subcase, f"Subcase {real_subcase}")

    # 2. READ PCOMP PLY INFO (ANGLE & THICKNESS)
    ply_metadata = {}
    try:
        eid_to_pid = {}
        if hasattr(h5.root.NASTRAN.INPUT.ELEMENT, 'CQUAD4'):
            for row in h5.root.NASTRAN.INPUT.ELEMENT.CQUAD4.read(): eid_to_pid[row['EID']] = row['PID']
        if hasattr(h5.root.NASTRAN.INPUT.ELEMENT, 'CTRIA3'):
            for row in h5.root.NASTRAN.INPUT.ELEMENT.CTRIA3.read(): eid_to_pid[row['EID']] = row['PID']
                
        pid_to_plies = {}
        if hasattr(h5.root.NASTRAN.INPUT.PROPERTY, 'PCOMP'):
            for row in h5.root.NASTRAN.INPUT.PROPERTY.PCOMP.read():
                pid = row['ID'] if 'ID' in row.dtype.names else row['PID']
                if 'THETA' in row.dtype.names and 'T' in row.dtype.names:
                    pid_to_plies[pid] = {}
                    ply_counter = 1
                    for theta, t in zip(row['THETA'], row['T']):
                        if t > 0: # Only map actual physical layers
                            pid_to_plies[pid][ply_counter] = f"Ply {ply_counter} [{theta}° | t={t}]"
                            ply_counter += 1

        for eid, pid in eid_to_pid.items():
            if pid in pid_to_plies: ply_metadata[eid] = pid_to_plies[pid]
    except Exception: pass

    def get_ply_string(eid, ply_id):
        try: return ply_metadata[eid][int(ply_id)]
        except: return f"Ply {ply_id}"

    # 3. READ COMPOSITE STRESSES (QUAD & TRIA)
    df_list = []
    element_types = ['QUAD4_COMP', 'TRIA3_COMP']
    
    try:
        stress_group = h5.root.NASTRAN.RESULT.ELEMENTAL.STRESS
        for elem_type in element_types:
            if hasattr(stress_group, elem_type):
                data = getattr(stress_group, elem_type).read()
                cols = data.dtype.names
                
                df_temp = pd.DataFrame()
                df_temp['Element ID'] = data['EID']
                df_temp['Element Type'] = elem_type 
                
                if 'PLY' in cols: raw_plies = data['PLY']
                elif 'LAYER' in cols: raw_plies = data['LAYER']
                else: raw_plies = np.ones_like(data['EID'])
                
                # Apply the smart PCOMP mapping
                df_temp['Ply Info'] = [get_ply_string(e, p) for e, p in zip(data['EID'], raw_plies)]
                
                df_temp['Load_Case'] = [get_load_case_name(d_id) for d_id in data['DOMAIN_ID']] if 'DOMAIN_ID' in cols else "Subcase 1"
                
                # In-Plane
                if 'X1' in cols and 'Y1' in cols and 'T1' in cols:
                    df_temp['Sigma_11'] = data['X1']; df_temp['Sigma_22'] = data['Y1']; df_temp['Tau_12'] = data['T1']
                elif 'X' in cols and 'Y' in cols and 'TXY' in cols:
                    df_temp['Sigma_11'] = data['X']; df_temp['Sigma_22'] = data['Y']; df_temp['Tau_12'] = data['TXY']
                elif 'X_NORMAL' in cols and 'Y_NORMAL' in cols and 'XY_SHEAR' in cols:
                    df_temp['Sigma_11'] = data['X_NORMAL']; df_temp['Sigma_22'] = data['Y_NORMAL']; df_temp['Tau_12'] = data['XY_SHEAR']
                else: continue
                
                # Interlaminar
                if 'L1' in cols and 'L2' in cols:
                    df_temp['Tau_13'] = data['L1']; df_temp['Tau_23'] = data['L2']
                elif 'XZ_SHEAR' in cols and 'YZ_SHEAR' in cols:
                    df_temp['Tau_13'] = data['XZ_SHEAR']; df_temp['Tau_23'] = data['YZ_SHEAR']
                else:
                    df_temp['Tau_13'] = 0.0; df_temp['Tau_23'] = 0.0

                df_list.append(df_temp)
                print(f"    -> Extracted {len(df_temp)} records for {elem_type}")

        if not df_list:
            print("  [X] Error: No composite stress data found.")
            h5.close(); return None
            
        df_raw = pd.concat(df_list, ignore_index=True)
        df_raw['Load_Case'] = df_raw['Load_Case'].astype('category')
        df_raw['Element Type'] = df_raw['Element Type'].astype('category')
        
    except Exception as e: 
        print(f"  [X] Fatal extraction error: {e}"); h5.close(); return None
        
    h5.close()
    gc.collect() 
    return df_raw

# ==============================================================================
# 5. 5-SHEET PROFESSIONAL EXCEL EXPORT
# ==============================================================================
def save_composite_report_5sheets(df_res, props, output_excel, base_name, img_dict):
    print(f"\n[STEP 7] 📑 Compiling 5-Sheet Professional Report: {os.path.basename(output_excel)}...")
    writer = pd.ExcelWriter(output_excel, engine='xlsxwriter') 
    wb = writer.book
    
    fmt_title = wb.add_format({'bold':True, 'font_size':16, 'fg_color':'#1F497D', 'font_color':'white', 'border':1, 'align':'center', 'valign':'vcenter'})
    fmt_head  = wb.add_format({'bold':True, 'fg_color':'#2F75B5', 'font_color':'white', 'border':1, 'align':'center', 'valign':'vcenter', 'text_wrap': True})
    fmt_gen   = wb.add_format({'border':1, 'align':'center'}) 
    fmt_num2  = wb.add_format({'num_format':'0.00', 'border':1, 'align':'center'}) 
    fmt_num4  = wb.add_format({'num_format':'0.0000', 'border':1, 'align':'center'}) 
    fmt_red   = wb.add_format({'font_color':'#9C0006', 'bg_color':'#FFC7CE', 'border':1, 'align':'center'}) 
    fmt_grn   = wb.add_format({'font_color':'#006100', 'bg_color':'#C6EFCE', 'border':1, 'align':'center'}) 

    # --- SHEET 1: EXECUTIVE SUMMARY ---
    ws_sum = wb.add_worksheet('1_Executive_Summary')
    ws_sum.merge_range('B2:H3', f"AEROSPACE COMPOSITE ANALYSIS REPORT", fmt_title)
    ws_sum.write('B4', f"File: {base_name}.h5", wb.add_format({'italic':True}))
    
    ws_sum.write('B6', "Applied Material Properties", fmt_head)
    row_m = 7
    for k, v in props.items():
        ws_sum.write(row_m, 1, k, fmt_gen); ws_sum.write(row_m, 2, v, fmt_gen)
        row_m += 1

    headers = ["Load Case", "Most Critical Element", "Element Type", "Critical Ply Info", "Lowest MS (Tsai-Wu)", "Driving Failure Mode"]
    for col, head in enumerate(headers): ws_sum.write(row_m + 1, 1+col, head, fmt_head)
    ws_sum.set_row(row_m + 1, 30) 
    
    row_idx = row_m + 2
    if 'Load_Case' in df_res.columns:
        grouped = df_res.groupby('Load_Case', observed=True)
        for lc_name, group in grouped:
            if group.empty: continue
            worst_row = group.loc[group['Tsai_Wu_MS'].idxmin()]
            
            ws_sum.write(row_idx, 1, str(lc_name), fmt_gen)
            ws_sum.write(row_idx, 2, worst_row['Element ID'], fmt_gen)
            ws_sum.write(row_idx, 3, worst_row['Element Type'], fmt_gen)
            ws_sum.write(row_idx, 4, worst_row['Ply Info'], fmt_gen)
            ws_sum.write(row_idx, 5, worst_row['Tsai_Wu_MS'], fmt_num4)
            ws_sum.write(row_idx, 6, worst_row['Hashin_Mode'], fmt_gen)
            row_idx += 1
            
    ws_sum.set_column('B:B', 30); ws_sum.set_column('C:H', 20)

    # --- HELPER FUNCTION FOR DATA SHEETS ---
    def write_sheet(df, name, specific_cols):
        df_sub = df[specific_cols]
        df_sub.to_excel(writer, sheet_name=name, index=False)
        ws = writer.sheets[name]
        ws.autofilter(0, 0, len(df_sub), len(df_sub.columns)-1) 
        ws.freeze_panes(1, 0) 
        
        last_row = len(df_sub) + 1
        for i, col in enumerate(df_sub.columns):
            ws.write(0, i, col, fmt_head)
            if col == 'Load_Case': ws.set_column(i, i, 25, fmt_gen)
            elif 'Element' in col or 'Ply' in col: ws.set_column(i, i, 16, fmt_gen)
            elif 'Sigma' in col or 'Tau' in col: ws.set_column(i, i, 12, fmt_num2)
            elif 'Mode' in col: ws.set_column(i, i, 20, fmt_gen)
            elif 'FI' in col: ws.set_column(i, i, 15, fmt_num4)
            elif 'MS' in col: 
                ws.set_column(i, i, 15, fmt_num4)
                col_let = xl_col_to_name(i)
                ws.conditional_format(f'{col_let}2:{col_let}{last_row}', {'type':'cell', 'criteria':'<', 'value':0, 'format':fmt_red})
                ws.conditional_format(f'{col_let}2:{col_let}{last_row}', {'type':'cell', 'criteria':'>=', 'value':0, 'format':fmt_grn})
            else: ws.set_column(i, i, 12, fmt_gen)

    # --- SHEET 2: RAW STRESSES ---
    raw_cols = ['Load_Case', 'Element ID', 'Element Type', 'Ply Info', 'Sigma_11', 'Sigma_22', 'Tau_12', 'Tau_13', 'Tau_23']
    write_sheet(df_res, '2_Raw_Stresses', raw_cols)

    # --- SHEET 3: BASIC CRITERIA ---
    basic_cols = ['Load_Case', 'Element ID', 'Ply Info', 'Max_Stress_FI', 'Max_Stress_MS', 'Max_Stress_Mode', 'Max_Strain_FI', 'Max_Strain_MS']
    write_sheet(df_res, '3_Basic_Criteria', basic_cols)

    # --- SHEET 4: ADVANCED CRITERIA ---
    adv_cols = ['Load_Case', 'Element ID', 'Ply Info', 'Tsai_Hill_FI', 'Tsai_Hill_MS', 'Tsai_Wu_FI', 'Tsai_Wu_MS', 'Hashin_FI', 'Hashin_MS', 'Hashin_Mode', 'Interlaminar_FI', 'Interlaminar_MS']
    write_sheet(df_res, '4_Advanced_Criteria', adv_cols)

    # --- SHEET 5: VISUAL GRAPHS ---
    ws_g = wb.add_worksheet('5_Visual_Graphs')
    ws_g.merge_range('B2:M3', "FAILURE ENVELOPE VISUALIZATION", fmt_title)
    ws_g.hide_gridlines(2)
    
    row_img = 5
    for lc, img_file in img_dict.items():
        ws_g.write(f'B{row_img}', f"📈 RESULTS FOR: {lc}", wb.add_format({'bold':True, 'font_size':14, 'bg_color':'#D9E1F2', 'border':1}))
        if os.path.exists(img_file):
            ws_g.insert_image(f'B{row_img+2}', img_file)
        row_img += 38

    writer.close()
    
    # Auto-cleanup PNGs
    for img_file in img_dict.values():
        try: os.remove(img_file)
        except: pass

# ==============================================================================
# 6. SMART MANUAL CALCULATOR (MODE 2)
# ==============================================================================
def manual_composite_input():
    print("\n" + "╔════════════════════════════════════════════════════════════╗")
    print("║        🧮 KALKULATOR MANUAL KOMPOSIT (SINGLE PLY)          ║")
    print("╚════════════════════════════════════════════════════════════╝")
    print("  [STEP 1] INPUT KEKUATAN MATERIAL (STRENGTHS)")
    while True:
        try:
            Xt = float(input("  Xt (Kekuatan Tarik Serat)  : "))
            Xc = float(input("  Xc (Kekuatan Tekan Serat)  : "))
            Yt = float(input("  Yt (Kekuatan Tarik Matriks): "))
            Yc = float(input("  Yc (Kekuatan Tekan Matriks): "))
            S  = float(input("  S  (Kekuatan Geser In-Plane): "))
            S13 = float(input("  S13 (Kekuatan Geser ILS XZ) : "))
            S23 = float(input("  S23 (Kekuatan Geser ILS YZ) : "))
            break
        except ValueError: print("  [X] Harap masukkan angka saja.\n")

    print("\n  [STEP 1B] INPUT KEKAKUAN MATERIAL (STIFFNESS)")
    while True:
        try:
            E1   = float(input("  E1 (Modulus Searah Serat) : "))
            E2   = float(input("  E2 (Modulus Tegak Lurus)  : "))
            nu12 = float(input("  nu12 (Poisson Ratio)      : "))
            G12  = float(input("  G12 (Shear Modulus)       : "))
            break
        except ValueError: print("  [X] Harap masukkan angka saja.\n")

    print("\n  [STEP 2] INPUT TEGANGAN YANG BEKERJA (APPLIED STRESSES)")
    while True:
        try:
            s1_inp  = float(input("  Sigma_11 (Arah Serat)   : "))
            s2_inp  = float(input("  Sigma_22 (Arah Matriks) : "))
            t12_inp = float(input("  Tau_12 (Geser In-Plane) : "))
            t13_inp = float(input("  Tau_13 (Geser ILS XZ)   : "))
            t23_inp = float(input("  Tau_23 (Geser ILS YZ)   : "))
            break
        except ValueError: print("  [X] Harap masukkan angka saja.\n")

    mat_props = {'E1': E1, 'E2': E2, 'nu12': nu12, 'G12': G12, 'Xt': Xt, 'Xc': Xc, 'Yt': Yt, 'Yc': Yc, 'S': S, 'S13': S13, 'S23': S23}
    
    # Create a dummy dataframe to pass into the exact same engine
    dummy_df = pd.DataFrame({'Sigma_11': [s1_inp], 'Sigma_22': [s2_inp], 'Tau_12': [t12_inp], 'Tau_13': [t13_inp], 'Tau_23': [t23_inp]})
    res = analyze_composite_failure(dummy_df, mat_props)
    
    print("\n" + "═"*60)
    print(" 📊 HASIL ANALISIS KOMPOSIT (ALL CRITERIA)")
    print("═"*60)
    
    print("  [MAXIMUM STRESS]")
    print(f"  > Failure Index    : {res['Max_Stress_FI'].iloc[0]:.4f}")
    print(f"  > Margin of Safety : {res['Max_Stress_MS'].iloc[0]:.4f}")
    print(f"  > Mode Kegagalan   : {res['Max_Stress_Mode'].iloc[0]}\n")
        
    print("  [MAXIMUM STRAIN]")
    print(f"  > Failure Index    : {res['Max_Strain_FI'].iloc[0]:.4f}")
    print(f"  > Margin of Safety : {res['Max_Strain_MS'].iloc[0]:.4f}\n")
        
    print("  [TSAI-HILL]")
    print(f"  > Failure Index    : {res['Tsai_Hill_FI'].iloc[0]:.4f}")
    print(f"  > Margin of Safety : {res['Tsai_Hill_MS'].iloc[0]:.4f}\n")
        
    print("  [TSAI-WU]")
    print(f"  > Failure Index    : {res['Tsai_Wu_FI'].iloc[0]:.4f}")
    print(f"  > Margin of Safety : {res['Tsai_Wu_MS'].iloc[0]:.4f}\n")

    print("  [HASHIN (1980)]")
    print(f"  > Failure Index    : {res['Hashin_FI'].iloc[0]:.4f}")
    print(f"  > Margin of Safety : {res['Hashin_MS'].iloc[0]:.4f}")
    print(f"  > Mode Kegagalan   : {res['Hashin_Mode'].iloc[0]}\n")

    print("  [INTERLAMINAR SHEAR (DELAMINATION)]")
    print(f"  > Failure Index    : {res['Interlaminar_FI'].iloc[0]:.4f}")
    print(f"  > Margin of Safety : {res['Interlaminar_MS'].iloc[0]:.4f}\n")
        
    print("════════════════════════════════════════════════════════════\n")
    input("Tekan Enter untuk kembali...")

# ==============================================================================
# 7. MAIN UI / CLI APPLICATION FLOW
# ==============================================================================
def main():
    while True:
        os.system('cls' if os.name == 'nt' else 'clear')
        print("╔═════════════════════════════════════════════════════════════════╗")
        print("║      🚀 AEROSPACE COMPOSITE EXPERT SUITE (MULTI-CRITERIA)       ║")
        print("╚═════════════════════════════════════════════════════════════════╝")
        print("Select the operation mode you want to use:")
        print("  [1] 📁 Analyze .h5 File from Nastran (Full Batch & 5-Sheet Excel)")
        print("  [2] 🧮 Manual Calculator (Single Ply Analysis)")
        print("  [0] Keluar Aplikasi")
        
        mode = input("\n  -> Enter choice (0/1/2): ").strip()
        
        if mode == '0':
            break
        elif mode == '2':
            manual_composite_input()
        elif mode == '1':
            script_dir = os.path.dirname(os.path.abspath(__file__))
            
            print("\n[STEP 1] 📁 SELECT NASTRAN H5 FILE")
            h5_files = glob.glob(os.path.join(script_dir, "*.h5")) 
            if not h5_files:
                print(f"  [X] No .h5 files found in folder:\n{script_dir}"); input(); return

            if len(h5_files) == 1:
                selected_h5 = h5_files[0]
                print(f"  -> Auto-selected the only H5 file: {os.path.basename(selected_h5)}")
            else:
                for idx, file in enumerate(h5_files): print(f"  [{idx + 1}] {os.path.basename(file)}")
                try: choice = int(input("  -> Select H5 file number: ").strip()); selected_h5 = h5_files[choice - 1]
                except: print("  [X] Invalid input. Aborting."); input(); continue

            bdf_files = glob.glob(os.path.join(script_dir, "*.bdf")) + glob.glob(os.path.join(script_dir, "*.dat"))
            selected_bdf = None
            if bdf_files:
                print("\n  📝 Load Case Names (Subtitles) available from:")
                print("  [0] Skip (Use default IDs)")
                for idx, file in enumerate(bdf_files): print(f"  [{idx + 1}] {os.path.basename(file)}")
                try:
                    choice = int(input("  -> Select BDF/DAT file number for Subtitles: ").strip())
                    if choice > 0: selected_bdf = bdf_files[choice - 1]
                except: pass

            print("\n[STEP 2] 🔍 PRE-SCANNING MODEL")
            df_raw = extract_composite_h5(selected_h5, bdf_file=selected_bdf)
            if df_raw is None or df_raw.empty:
                print("  [X] File is empty or composite stresses not found."); input(); continue

            print("\n[STEP 3] ⚖️ UNIT CONVERSION")
            print("  [1] Pascal (Pa) \n  [2] Megapascal (MPa) - Default\n  [3] Gigapascal (GPa)\n  [4] PSI")
            unit_choice = input("  Select unit (1/2/3/4): ").strip()
            to_mpa = {'1': 1e-6, '2': 1.0, '3': 1000.0, '4': 0.00689475}
            conv_factor = to_mpa.get(unit_choice, 1.0)
                
            if conv_factor != 1.0:
                df_raw['Sigma_11'] *= conv_factor; df_raw['Sigma_22'] *= conv_factor; df_raw['Tau_12'] *= conv_factor
                df_raw['Tau_13'] *= conv_factor; df_raw['Tau_23'] *= conv_factor

            print("\n[STEP 4] ⚙️ MATERIAL CONFIGURATION (Auto-Detect from .H5)")
            extracted_mats = extract_mat8_h5(selected_h5, conv_factor)
            mat_selected = None
            
            if extracted_mats:
                print(f"  🕵️‍♂️ Nastran found {len(extracted_mats)} MAT8 materials!")
                for mid, props in extracted_mats.items():
                    print(f"  [ Material ID: {mid} ]")
                    print(f"    E1 = {props['E1']:,.1f} | E2 = {props['E2']:,.1f} | G12 = {props['G12']:,.1f} | nu12 = {props['nu12']}")
                    print(f"    Xt = {props['Xt']:,.1f} | Xc = {props['Xc']:,.1f}")
                    print(f"    Yt = {props['Yt']:,.1f} | Yc = {props['Yc']:,.1f}")
                    print(f"    S  = {props['S']:,.1f} | S13= {props['S13']:,.1f} | S23= {props['S23']:,.1f}")
                    print("-" * 55)
                    
                use_auto = input("\n  -> Gunakan data material otomatis dari H5 ini? (Y/N): ").strip().upper()
                if use_auto == 'Y':
                    if len(extracted_mats) == 1: mat_selected = list(extracted_mats.values())[0]
                    else:
                        try:
                            mid_sel = int(input("  -> Pilih Material ID: "))
                            mat_selected = extracted_mats.get(mid_sel, None)
                        except: pass
                                
            if not mat_selected:
                print("\n  --- INPUT MANUAL PROPERTI MATERIAL ---")
                while True:
                    try:
                        mat_selected = {
                            'E1': float(input("  E1 (Modulus Serat)  : ")), 'E2': float(input("  E2 (Modulus Matriks): ")),
                            'nu12': float(input("  nu12 (Poisson)      : ")), 'G12': float(input("  G12 (Shear Modulus) : ")),
                            'Xt': float(input("  Xt (Tarik Serat)    : ")), 'Xc': float(input("  Xc (Tekan Serat)    : ")),
                            'Yt': float(input("  Yt (Tarik Matriks)  : ")), 'Yc': float(input("  Yc (Tekan Matriks)  : ")),
                            'S': float(input("  S (Geser In-Plane)  : ")),
                            'S13': float(input("  S13 (Geser ILS XZ)  : ")), 'S23': float(input("  S23 (Geser ILS YZ)  : "))
                        }
                        break 
                    except ValueError: print("  [X] Input tidak valid!\n")

            print("\n[STEP 5] 📊 OUTPUT FILTERING")
            try:
                top_n = int(input("  -> Export Top 'N' worst plies per load case? (0 for ALL): ").strip() or 0)
            except: top_n = 0
                
            print("\n" + "═"*65)
            print(" 🚀 INITIATING EXPERT COMPOSITE ENGINE...")
            
            base_name = os.path.splitext(os.path.basename(selected_h5))[0]
            output_excel = os.path.join(script_dir, f"Composite_Report_{base_name}.xlsx")
            try: open(output_excel, 'a').close()
            except: print(f"\n[ERROR] Close the file {os.path.basename(output_excel)} first!"); input(); continue
            
            df_results = analyze_composite_failure(df_raw, mat_selected)
            
            if top_n > 0:
                print(f"  -> Filtering to Top {top_n} worst plies...")
                df_results = df_results.sort_values(['Load_Case', 'Tsai_Wu_FI'], ascending=[True, False])
                df_results = df_results.groupby('Load_Case', observed=True).head(top_n).reset_index(drop=True)
                
            img_dict = create_failure_envelopes(df_results, mat_selected, base_name)
            
            save_composite_report_5sheets(df_results, mat_selected, output_excel, base_name, img_dict)
            
            print("\n" + "═"*65)
            print(f" 🎉 SUCCESS! 5-SHEET EXCEL GENERATED: {os.path.basename(output_excel)}")
            print("═"*65 + "\n")
            
            input("Tekan Enter untuk kembali ke Menu Utama...")

if __name__ == "__main__":
    main()
