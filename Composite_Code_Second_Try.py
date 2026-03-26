# ==============================================================================
# SCRIPT LIMITATIONS & ASSUMPTIONS
# ==============================================================================
# 1. ELEMENT TYPES: This code only reads and calculates failure for 2D Shell 
#    composite elements (CQUAD4 and CTRIA3). 3D Solid elements are not supported.
#
# 2. MATERIAL DEFINITION: Assumes 2D Orthotropic materials (MAT8 property).
#
# 3. STRUCTURE TYPE: Designed for solid laminates. Failure modes unique to 
#    Sandwich Panels (e.g., core shear, face wrinkling) are not evaluated.
#
# 4. ENVIRONMENTAL EFFECTS: Failure is based on mechanical stresses only. 
#    Hygrothermal expansion effects are not explicitly calculated.
#
# 5. DATA SOURCE: Strictly parses MSC Nastran HDF5 (.h5) binary result files.
#
# 6. PUCK PARAMETERS: Uses standard inclination values (p12+=0.3, p12-=0.25, p22-=0.2).
# ==============================================================================

# ==============================================================================
# GOVERNING EQUATIONS
# ==============================================================================
# 1. MAXIMUM STRESS:
#    FI = max( s1/X, s2/Y, t12/S ) | RF = 1 / FI
#
# 2. MAXIMUM STRAIN (Hooke's Law Applied First):
#    eps1 = (s1/E1) - (nu21*s2/E2) | eps2 = -(nu12*s1/E1) + (s2/E2)
#    FI = max( eps1/eps1_u, eps2/eps2_u, gamma12/gamma12_u )
#
# 3. TSAI-HILL (Quadratic Interactive):
#    FI = (s1/X)^2 - (s1*s2)/X^2 + (s2/Y)^2 + (t12/S)^2 | RF = 1 / sqrt(FI)
#
# 4. TSAI-WU (Tensor Polynomial):
#    F1*s1 + F2*s2 + F11*s1^2 + F22*s2^2 + F66*t12^2 + 2*F12*s1*s2 = FI
#    RF calculated via ABC Quadratic Formula: (A)R^2 + (B)R - 1 = 0
#
# 5. HASHIN (1980):
#    Fiber Tension (s1>0): FI = (s1/Xt)^2 + (t12/S)^2
#    Fiber Compress (s1<0): FI = (|s1|/Xc)^2
#    Matrix Tension (s2>0): FI = (s2/Yt)^2 + (t12/S)^2
#    Matrix Compress (s2<0): FI = (|s2|/Yc)^2 + (t12/S)^2
#
# 6. PUCK 2D (Action Plane Theory):
#    Fiber Failure (FF): Same as Max Stress.
#    Inter-Fiber Failure (IFF):
#    Mode A (s2>0): FI = sqrt((t12/S)^2 + (1-p12+*Yt/S)^2*(s2/Yt)^2) + p12+*s2/S
#    Mode B (s2<0, low slope): FI = 1/S * (sqrt(t12^2 + (p12-*s2)^2) + p12-*s2)
#    Mode C (s2<0, high slope): FI = [(t12/(2*(1+p22-)*S))^2 + (s2/Yc)^2] * (Yc/|s2|)
#
# 7. INTERLAMINAR SHEAR:
#    FI = (t13/S13)^2 + (t23/S23)^2 | RF = 1 / sqrt(FI)
#
# 8. LAMINATE ABD MATRIX (CLPT):
#    [A] = sum(Qbar * dt)           -> Extensional Stiffness
#    [B] = 1/2 * sum(Qbar * dz^2)   -> Coupling Stiffness
#    [D] = 1/3 * sum(Qbar * dz^3)   -> Bending Stiffness
# ==============================================================================

# ==============================================================================
# IMPORT REQUIRED LIBRARIES
# ==============================================================================
import tables             # Library to read Nastran HDF5 binary database files
import pandas as pd       # Library for data manipulation, table creation, and Excel export
import numpy as np        # Library for high-performance vectorized mathematical operations
import os                 # Library to interact with the Operating System (paths, files)
import glob               # Library to search for files with specific extensions (like .h5)
import re                 # Library for Regular Expressions (used for text cleaning)
from tqdm import tqdm     # Library to display progress bars in terminal
from xlsxwriter.utility import xl_col_to_name # Library to format Excel columns
import matplotlib.pyplot as plt # Standard library for generating 2D graphs
import seaborn as sns     # Advanced visualization library for better-looking charts
import warnings           # Library to handle system warnings
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
    return text.replace("'", "").replace('"', '').strip()

def safe_filename(text):
    """ Replaces illegal characters in a string for safely saving image files """
    return re.sub(r'[^A-Za-z0-9_\-\.]', '_', text)

# ==============================================================================
# 2. PHYSICS ENGINE (INCLUDING PUCK, HASHIN & RESERVE FACTORS)
# ==============================================================================
def analyze_composite_failure(df, props):
    """
    Core mathematical engine to calculate Failure Indices (FI), Reserve Factors (RF), and Margins of Safety (MS).
    Includes: Max Stress, Max Strain, Tsai-Hill, Tsai-Wu, Hashin, Puck, and Interlaminar Shear.
    """
    E1, E2, nu12, G12 = props['E1'], props['E2'], props['nu12'], props['G12']
    
    # SAFEGUARD: Prevent ZeroDivisionError for strengths
    def enforce_limit(val): return abs(val) if abs(val) > 1e-6 else 1e12 
        
    Xt, Xc = enforce_limit(props['Xt']), enforce_limit(props['Xc'])
    Yt, Yc = enforce_limit(props['Yt']), enforce_limit(props['Yc'])
    S = enforce_limit(props['S'])
    S13 = enforce_limit(props.get('S13', S)) 
    S23 = enforce_limit(props.get('S23', S)) 

    # Extract applied stresses as numpy arrays for extreme calculation speed
    s1, s2, t12 = df['Sigma_11'].values, df['Sigma_22'].values, df['Tau_12'].values
    t13 = df['Tau_13'].values if 'Tau_13' in df.columns else np.zeros_like(s1)
    t23 = df['Tau_23'].values if 'Tau_23' in df.columns else np.zeros_like(s1)

    # Helper function for RF and MS calculation
    def calc_rf_ms(fi_array, is_quadratic=False):
        if is_quadratic:
            rf = np.where(fi_array > 1e-8, 1.0 / np.sqrt(fi_array), 999.9)
        else:
            rf = np.where(fi_array > 1e-8, 1.0 / fi_array, 999.9)
        ms = rf - 1.0
        return rf, ms

    # --- A. MAXIMUM STRESS CRITERION ---
    fi_1_maxs = np.where(s1 > 0, s1 / Xt, np.abs(s1) / Xc)
    fi_2_maxs = np.where(s2 > 0, s2 / Yt, np.abs(s2) / Yc)
    fi_12_maxs = np.abs(t12) / S 
    fi_max_stress = np.maximum.reduce([fi_1_maxs, fi_2_maxs, fi_12_maxs])
    rf_max_stress, ms_max_stress = calc_rf_ms(fi_max_stress, is_quadratic=False)

    modes_ms = []
    for f1, f2, fm, sig1, sig2 in zip(fi_1_maxs, fi_2_maxs, fi_max_stress, s1, s2):
        if fm < 1e-6: modes_ms.append("Safe")
        elif fm == f1 and sig1 > 0: modes_ms.append("Fiber Tension")
        elif fm == f1 and sig1 <= 0: modes_ms.append("Fiber Compression")
        elif fm == f2 and sig2 > 0: modes_ms.append("Matrix Tension")
        elif fm == f2 and sig2 <= 0: modes_ms.append("Matrix Compression")
        else: modes_ms.append("In-Plane Shear")

    # --- B. MAXIMUM STRAIN CRITERION ---
    eps_1T, eps_1C = Xt / E1, Xc / E1; eps_2T, eps_2C = Yt / E2, Yc / E2; gamma_12u = S / G12
    nu21 = (nu12 * E2) / E1 if E1 > 1e-6 else 0.0
    eps1_act = (s1 / E1) - (nu21 * s2 / E2); eps2_act = -(nu12 * s1 / E1) + (s2 / E2); gamma12_act = t12 / G12
    
    fi_1_maxe = np.where(eps1_act > 0, eps1_act / eps_1T, np.abs(eps1_act) / eps_1C)
    fi_2_maxe = np.where(eps2_act > 0, eps2_act / eps_2T, np.abs(eps2_act) / eps_2C)
    fi_12_maxe = np.abs(gamma12_act) / gamma_12u
    fi_max_strain = np.maximum.reduce([fi_1_maxe, fi_2_maxe, fi_12_maxe])
    rf_max_strain, ms_max_strain = calc_rf_ms(fi_max_strain, is_quadratic=False)

    # --- C. TSAI-HILL CRITERION ---
    X_hill = np.where(s1 > 0, Xt, Xc); Y_hill = np.where(s2 > 0, Yt, Yc)
    fi_tsai_hill = (s1 / X_hill)**2 - (s1 * s2) / (X_hill**2) + (s2 / Y_hill)**2 + (t12 / S)**2
    rf_tsai_hill, ms_tsai_hill = calc_rf_ms(fi_tsai_hill, is_quadratic=True)

    # --- D. TSAI-WU CRITERION ---
    F1 = (1.0 / Xt) - (1.0 / Xc); F2 = (1.0 / Yt) - (1.0 / Yc)
    F11 = 1.0 / (Xt * Xc); F22 = 1.0 / (Yt * Yc); F66 = 1.0 / (S**2)
    F12 = -0.5 * np.sqrt(F11 * F22) 
    fi_tsai_wu = F1*s1 + F2*s2 + F11*(s1**2) + F22*(s2**2) + F66*(t12**2) + 2*F12*s1*s2
    
    A_tw = F11*(s1**2) + F22*(s2**2) + F66*(t12**2) + 2*F12*s1*s2; B_tw = F1*s1 + F2*s2
    rf_tsai_wu = np.zeros_like(s1, dtype=float)
    mask_quad = A_tw > 1e-12 
    rf_tsai_wu[mask_quad] = (-B_tw[mask_quad] + np.sqrt(B_tw[mask_quad]**2 + 4 * A_tw[mask_quad])) / (2 * A_tw[mask_quad])
    mask_lin = (A_tw <= 1e-12) & (B_tw > 1e-12)
    rf_tsai_wu[mask_lin] = 1.0 / B_tw[mask_lin]
    rf_tsai_wu[(A_tw <= 1e-12) & (B_tw <= 1e-12)] = 999.9
    ms_tsai_wu = rf_tsai_wu - 1.0

    # --- E. HASHIN CRITERION (1980 2D) ---
    fi_h_ft = np.where(s1 > 0, (s1/Xt)**2 + (t12/S)**2, 0)
    fi_h_fc = np.where(s1 <= 0, (np.abs(s1)/Xc)**2, 0)
    fi_h_mt = np.where(s2 > 0, (s2/Yt)**2 + (t12/S)**2, 0)
    fi_h_mc = np.where(s2 <= 0, (np.abs(s2)/Yc)**2 + (t12/S)**2, 0)
    fi_hashin = np.maximum.reduce([fi_h_ft, fi_h_fc, fi_h_mt, fi_h_mc])
    rf_hashin, ms_hashin = calc_rf_ms(fi_hashin, is_quadratic=True)

    hashin_modes = []
    for ft, fc, mt, mc, fm in zip(fi_h_ft, fi_h_fc, fi_h_mt, fi_h_mc, fi_hashin):
        if fm < 1e-6: hashin_modes.append("Safe")
        elif fm == ft: hashin_modes.append("Fiber Tension")
        elif fm == fc: hashin_modes.append("Fiber Compression")
        elif fm == mt: hashin_modes.append("Matrix Tension")
        else: hashin_modes.append("Matrix Compression")

    # --- F. PUCK 2D CRITERION (Action Plane) ---
    p12_plus = 0.3; p12_minus = 0.25; p22_minus = 0.2 
    
    fi_p_ff_t = np.where(s1 >= 0, s1/Xt, 0)
    fi_p_ff_c = np.where(s1 < 0, np.abs(s1)/Xc, 0)
    
    termA = np.sqrt((t12/S)**2 + (1 - p12_plus*Yt/S)**2 * (s2/Yt)**2) + p12_plus*s2/S
    fi_p_iff_a = np.where(s2 > 0, termA, 0)
    
    R_tt_A = Yc / (2 * (1 + p22_minus))
    slope = np.abs(s2 / (np.abs(t12) + 1e-12))
    slope_limit = R_tt_A / S
    
    fi_p_iff_b = np.where((s2 <= 0) & (slope <= slope_limit), (np.sqrt(t12**2 + (p12_minus*s2)**2) + p12_minus*s2)/S, 0)
    fi_p_iff_c = np.where((s2 <= 0) & (slope > slope_limit), ((t12/(2*(1+p22_minus)*S))**2 + (s2/Yc)**2) * (Yc/(np.abs(s2)+1e-12)), 0)
    
    fi_puck = np.maximum.reduce([fi_p_ff_t, fi_p_ff_c, fi_p_iff_a, fi_p_iff_b, fi_p_iff_c])
    rf_puck, ms_puck = calc_rf_ms(fi_puck, is_quadratic=False) 

    puck_modes = []
    for fft, ffc, ia, ib, ic, fm in zip(fi_p_ff_t, fi_p_ff_c, fi_p_iff_a, fi_p_iff_b, fi_p_iff_c, fi_puck):
        if fm < 1e-6: puck_modes.append("Safe")
        elif fm == fft: puck_modes.append("FF Tension")
        elif fm == ffc: puck_modes.append("FF Compression")
        elif fm == ia: puck_modes.append("IFF Mode A (Tension)")
        elif fm == ib: puck_modes.append("IFF Mode B (Compression)")
        else: puck_modes.append("IFF Mode C (Compression)")

    # --- G. INTERLAMINAR SHEAR (DELAMINATION) ---
    fi_ils = (t13 / S13)**2 + (t23 / S23)**2
    rf_ils, ms_ils = calc_rf_ms(fi_ils, is_quadratic=True)

    # --- COMPILE OUTPUT DATAFRAME ---
    out_df = df.copy()
    out_df['Max_Stress_FI'] = fi_max_stress; out_df['Max_Stress_RF'] = rf_max_stress; out_df['Max_Stress_MS'] = ms_max_stress; out_df['Max_Stress_Mode'] = modes_ms
    out_df['Max_Strain_FI'] = fi_max_strain; out_df['Max_Strain_RF'] = rf_max_strain; out_df['Max_Strain_MS'] = ms_max_strain
    out_df['Tsai_Hill_FI'] = fi_tsai_hill; out_df['Tsai_Hill_RF'] = rf_tsai_hill; out_df['Tsai_Hill_MS'] = ms_tsai_hill
    out_df['Tsai_Wu_FI'] = fi_tsai_wu; out_df['Tsai_Wu_RF'] = rf_tsai_wu; out_df['Tsai_Wu_MS'] = ms_tsai_wu
    out_df['Hashin_FI'] = fi_hashin; out_df['Hashin_RF'] = rf_hashin; out_df['Hashin_MS'] = ms_hashin; out_df['Hashin_Mode'] = hashin_modes
    out_df['Puck_FI'] = fi_puck; out_df['Puck_RF'] = rf_puck; out_df['Puck_MS'] = ms_puck; out_df['Puck_Mode'] = puck_modes
    out_df['Interlaminar_FI'] = fi_ils; out_df['Interlaminar_RF'] = rf_ils; out_df['Interlaminar_MS'] = ms_ils

    return out_df

# ==============================================================================
# 3. ABD LAMINATE MATRIX CALCULATOR (CLASSICAL LAMINATED PLATE THEORY)
# ==============================================================================
def calculate_abd_matrices(pcomp_raw_dict, props):
    """ Builds the A (Extensional), B (Coupling), and D (Bending) matrices for PCOMP """
    if not pcomp_raw_dict: return None
    E1, E2, nu12, G12 = props['E1'], props['E2'], props['nu12'], props['G12']
    if E1 < 1e-6 or E2 < 1e-6: return None
    
    nu21 = nu12 * E2 / E1
    den = 1.0 - nu12 * nu21
    Q11 = E1 / den; Q22 = E2 / den; Q12 = nu12 * E2 / den; Q66 = G12
    
    abd_rows = []
    for pid, plies in pcomp_raw_dict.items():
        total_t = sum([p['t'] for p in plies])
        z_curr = -total_t / 2.0
        A = np.zeros((3,3)); B = np.zeros((3,3)); D = np.zeros((3,3))
        
        for ply in plies:
            t = ply['t']; theta = np.radians(ply['theta'])
            z_prev = z_curr; z_curr += t
            m = np.cos(theta); n = np.sin(theta)
            m2, n2, m4, n4 = m**2, n**2, m**4, n**4
            
            Qbar = np.zeros((3,3))
            Qbar[0,0] = Q11*m4 + 2*(Q12 + 2*Q66)*m2*n2 + Q22*n4
            Qbar[1,1] = Q11*n4 + 2*(Q12 + 2*Q66)*m2*n2 + Q22*m4
            Qbar[0,1] = (Q11 + Q22 - 4*Q66)*m2*n2 + Q12*(m4 + n4); Qbar[1,0] = Qbar[0,1]
            Qbar[2,2] = (Q11 + Q22 - 2*Q12 - 2*Q66)*m2*n2 + Q66*(m4 + n4)
            Qbar[0,2] = (Q11 - Q12 - 2*Q66)*n*m**3 + (Q12 - Q22 + 2*Q66)*n**3*m; Qbar[2,0] = Qbar[0,2]
            Qbar[1,2] = (Q11 - Q12 - 2*Q66)*m*n**3 + (Q12 - Q22 + 2*Q66)*m**3*n; Qbar[2,1] = Qbar[1,2]

            A += Qbar * (z_curr - z_prev)
            B += 0.5 * Qbar * (z_curr**2 - z_prev**2)
            D += (1.0/3.0) * Qbar * (z_curr**3 - z_prev**3)
            
        abd_rows.append({'Property ID': pid, 'Matrix Type': 'A (Extensional)', '11': A[0,0], '12': A[0,1], '16': A[0,2], '22': A[1,1], '26': A[1,2], '66': A[2,2]})
        abd_rows.append({'Property ID': pid, 'Matrix Type': 'B (Coupling)', '11': B[0,0], '12': B[0,1], '16': B[0,2], '22': B[1,1], '26': B[1,2], '66': B[2,2]})
        abd_rows.append({'Property ID': pid, 'Matrix Type': 'D (Bending)', '11': D[0,0], '12': D[0,1], '16': D[0,2], '22': D[1,1], '26': D[1,2], '66': D[2,2]})
        
    return pd.DataFrame(abd_rows)

# ==============================================================================
# 4. VISUALIZATION (TSAI-WU FAILURE ENVELOPES)
# ==============================================================================
def create_failure_envelopes(df_res, props, base_name, unit_str):
    """ Plots the Tsai-Wu Failure Ellipse and scatters the element stresses inside/outside it """
    print("\n[STEP 6] 🎨 Generating Tsai-Wu Failure Envelopes...")
    sns.set_style("darkgrid")
    
    Xt, Xc = abs(props['Xt']), abs(props['Xc'])
    Yt, Yc = abs(props['Yt']), abs(props['Yc'])
    
    F1 = (1.0 / Xt) - (1.0 / Xc); F2 = (1.0 / Yt) - (1.0 / Yc)
    F11 = 1.0 / (Xt * Xc); F22 = 1.0 / (Yt * Yc)
    F12 = -0.5 * np.sqrt(F11 * F22) 
    
    s1_vals = np.linspace(-Xc*1.3, Xt*1.3, 400)
    s2_vals = np.linspace(-Yc*1.3, Yt*1.3, 400)
    S1, S2 = np.meshgrid(s1_vals, s2_vals)
    Z = F11*(S1**2) + F22*(S2**2) + F1*S1 + F2*S2 + 2*F12*S1*S2
    
    load_cases = df_res['Load_Case'].unique()
    images_dict = {}

    col_s1 = 'Sigma_11'
    col_s2 = 'Sigma_22'

    for lc in load_cases:
        safe_lc = safe_filename(str(lc))
        df_lc = df_res[df_res['Load_Case'] == lc]
        
        plt.figure(figsize=(9, 7))
        plt.contour(S1, S2, Z, levels=[1.0], colors='red', linewidths=2.5, linestyles='dashed')
        
        safe_mask = df_lc['Tsai_Wu_FI'] < 1.0
        failed_mask = df_lc['Tsai_Wu_FI'] >= 1.0
        
        plt.scatter(df_lc[safe_mask][col_s1], df_lc[safe_mask][col_s2], color='green', alpha=0.6, s=25, label='Safe Elements (FI < 1)')
        plt.scatter(df_lc[failed_mask][col_s1], df_lc[failed_mask][col_s2], color='red', marker='X', s=45, label='Failed Elements (FI >= 1)')
        
        plt.axhline(0, color='black', linewidth=1)
        plt.axvline(0, color='black', linewidth=1)
        
        plt.title(f"Tsai-Wu Failure Envelope: {lc}\n(Plotted for Tau_12 = 0)", fontsize=13, fontweight='bold')
        plt.xlabel(f"Sigma 11 (Fiber Direction) [{unit_str}]", fontweight='bold')
        plt.ylabel(f"Sigma 22 (Matrix Direction) [{unit_str}]", fontweight='bold')
        plt.legend(loc='upper left', frameon=True, shadow=True)
        plt.tight_layout()
        
        img_name = f"Envelope_{safe_lc}.png"
        plt.savefig(img_name, dpi=150)
        plt.close()
        images_dict[lc] = img_name
        
    return images_dict

# ==============================================================================
# 5. H5 DATA MINING (STRESSES, MAT8, PCOMP, SUBTITLES)
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
    print(f"    [+] Extracting Data & PCOMP from H5: {os.path.basename(input_file)} ...")
    h5 = tables.open_file(input_file, mode="r") 
    
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

    ply_metadata = {}; pcomp_raw = {}
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
                    pid_to_plies[pid] = {}; raw_plies = []
                    ply_counter = 1
                    for theta, t in zip(row['THETA'], row['T']):
                        if t > 0: 
                            pid_to_plies[pid][ply_counter] = f"Ply {ply_counter} [{theta}° | t={t}]"
                            raw_plies.append({'ply': ply_counter, 'theta': theta, 't': t})
                            ply_counter += 1
                    pcomp_raw[pid] = raw_plies

        for eid, pid in eid_to_pid.items():
            if pid in pid_to_plies: ply_metadata[eid] = pid_to_plies[pid]
    except Exception: pass

    def get_ply_string(eid, ply_id):
        try: return ply_metadata[eid][int(ply_id)]
        except: return f"Ply {ply_id}"

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
                
                df_temp['Ply Info'] = [get_ply_string(e, p) for e, p in zip(data['EID'], raw_plies)]
                df_temp['Load_Case'] = [get_load_case_name(d_id) for d_id in data['DOMAIN_ID']] if 'DOMAIN_ID' in cols else "Subcase 1"
                
                if 'X1' in cols and 'Y1' in cols and 'T1' in cols:
                    df_temp['Sigma_11'] = data['X1']; df_temp['Sigma_22'] = data['Y1']; df_temp['Tau_12'] = data['T1']
                elif 'X' in cols and 'Y' in cols and 'TXY' in cols:
                    df_temp['Sigma_11'] = data['X']; df_temp['Sigma_22'] = data['Y']; df_temp['Tau_12'] = data['TXY']
                elif 'X_NORMAL' in cols and 'Y_NORMAL' in cols and 'XY_SHEAR' in cols:
                    df_temp['Sigma_11'] = data['X_NORMAL']; df_temp['Sigma_22'] = data['Y_NORMAL']; df_temp['Tau_12'] = data['XY_SHEAR']
                else: continue
                
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
            h5.close(); return None, None
            
        df_raw = pd.concat(df_list, ignore_index=True)
        df_raw['Load_Case'] = df_raw['Load_Case'].astype('category')
        df_raw['Element Type'] = df_raw['Element Type'].astype('category')
        
    except Exception as e: 
        print(f"  [X] Fatal extraction error: {e}"); h5.close(); return None, None
        
    h5.close()
    gc.collect() 
    return df_raw, pcomp_raw

# ==============================================================================
# 6. EXCEL EXPORT (NO NUMBERING IN SHEET NAMES)
# ==============================================================================
def save_composite_report_6sheets(df_res, df_abd, props, output_excel, base_name, img_dict, unit_str):
    print(f"\n[STEP 7] 📑 Compiling Expert Report: {os.path.basename(output_excel)}...")
    writer = pd.ExcelWriter(output_excel, engine='xlsxwriter') 
    wb = writer.book
    
    fmt_title = wb.add_format({'bold':True, 'font_size':16, 'fg_color':'#1F497D', 'font_color':'white', 'border':1, 'align':'center', 'valign':'vcenter'})
    fmt_head  = wb.add_format({'bold':True, 'fg_color':'#2F75B5', 'font_color':'white', 'border':1, 'align':'center', 'valign':'vcenter', 'text_wrap': True})
    fmt_gen   = wb.add_format({'border':1, 'align':'center'}) 
    fmt_num2  = wb.add_format({'num_format':'0.00', 'border':1, 'align':'center'}) 
    fmt_num4  = wb.add_format({'num_format':'0.0000', 'border':1, 'align':'center'}) 
    fmt_red   = wb.add_format({'font_color':'#9C0006', 'bg_color':'#FFC7CE', 'border':1, 'align':'center'}) 
    fmt_grn   = wb.add_format({'font_color':'#006100', 'bg_color':'#C6EFCE', 'border':1, 'align':'center'}) 

    # Dynamic Column renaming for units (Applied to headers only)
    df_formatted = df_res.rename(columns={
        'Sigma_11': f'Sigma_11 [{unit_str}]',
        'Sigma_22': f'Sigma_22 [{unit_str}]',
        'Tau_12': f'Tau_12 [{unit_str}]',
        'Tau_13': f'Tau_13 [{unit_str}]',
        'Tau_23': f'Tau_23 [{unit_str}]'
    })

    # --- SHEET: EXECUTIVE SUMMARY ---
    ws_sum = wb.add_worksheet('Executive_Summary')
    ws_sum.merge_range('B2:H3', f"AEROSPACE COMPOSITE ANALYSIS REPORT", fmt_title)
    ws_sum.write('B4', f"File: {base_name}.h5", wb.add_format({'italic':True}))
    
    ws_sum.write('B6', "Applied Material Properties", fmt_head)
    row_m = 7
    for k, v in props.items():
        unit_label = "" if k == 'nu12' else f" [{unit_str}]"
        ws_sum.write(row_m, 1, f"{k}{unit_label}", fmt_gen)
        ws_sum.write(row_m, 2, v, fmt_gen)
        row_m += 1

    headers = ["Load Case", "Most Critical Element", "Critical Ply Info", "Lowest MS (Tsai-Wu)", "Tsai-Wu RF", "Driving Failure Mode"]
    for col, head in enumerate(headers): ws_sum.write(row_m + 1, 1+col, head, fmt_head)
    ws_sum.set_row(row_m + 1, 30) 
    
    row_idx = row_m + 2
    if 'Load_Case' in df_formatted.columns:
        grouped = df_formatted.groupby('Load_Case', observed=True)
        for lc_name, group in grouped:
            if group.empty: continue
            worst_row = group.loc[group['Tsai_Wu_MS'].idxmin()]
            
            ws_sum.write(row_idx, 1, str(lc_name), fmt_gen)
            ws_sum.write(row_idx, 2, worst_row['Element ID'], fmt_gen)
            ws_sum.write(row_idx, 3, worst_row['Ply Info'], fmt_gen)
            ws_sum.write(row_idx, 4, worst_row['Tsai_Wu_MS'], fmt_num4)
            ws_sum.write(row_idx, 5, worst_row['Tsai_Wu_RF'], fmt_num4)
            ws_sum.write(row_idx, 6, worst_row['Puck_Mode'], fmt_gen)
            row_idx += 1
            
    ws_sum.set_column('B:B', 30); ws_sum.set_column('C:H', 20)

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
            elif 'Element' in col or 'Ply' in col: ws.set_column(i, i, 18, fmt_gen)
            elif 'Sigma' in col or 'Tau' in col: ws.set_column(i, i, 15, fmt_num2)
            elif 'Mode' in col: ws.set_column(i, i, 22, fmt_gen)
            elif 'FI' in col or 'RF' in col: ws.set_column(i, i, 15, fmt_num4)
            elif 'MS' in col: 
                ws.set_column(i, i, 15, fmt_num4)
                col_let = xl_col_to_name(i)
                ws.conditional_format(f'{col_let}2:{col_let}{last_row}', {'type':'cell', 'criteria':'<', 'value':0, 'format':fmt_red})
                ws.conditional_format(f'{col_let}2:{col_let}{last_row}', {'type':'cell', 'criteria':'>=', 'value':0, 'format':fmt_grn})
            else: ws.set_column(i, i, 12, fmt_gen)

    # --- SHEET: RAW STRESSES ---
    raw_cols = ['Load_Case', 'Element ID', 'Element Type', 'Ply Info', 
                f'Sigma_11 [{unit_str}]', f'Sigma_22 [{unit_str}]', 
                f'Tau_12 [{unit_str}]', f'Tau_13 [{unit_str}]', f'Tau_23 [{unit_str}]']
    write_sheet(df_formatted, 'Raw_Stresses', raw_cols)

    # --- SHEET: BASIC CRITERIA ---
    basic_cols = ['Load_Case', 'Element ID', 'Ply Info', 'Max_Stress_FI', 'Max_Stress_RF', 'Max_Stress_MS', 'Max_Stress_Mode', 'Max_Strain_FI', 'Max_Strain_RF', 'Max_Strain_MS']
    write_sheet(df_formatted, 'Basic_Criteria', basic_cols)

    # --- SHEET: ADVANCED CRITERIA ---
    adv_cols = ['Load_Case', 'Element ID', 'Ply Info', 'Tsai_Hill_FI', 'Tsai_Hill_RF', 'Tsai_Hill_MS', 'Tsai_Wu_FI', 'Tsai_Wu_RF', 'Tsai_Wu_MS', 'Hashin_FI', 'Hashin_RF', 'Hashin_MS', 'Hashin_Mode', 'Puck_FI', 'Puck_RF', 'Puck_MS', 'Puck_Mode', 'Interlaminar_FI', 'Interlaminar_RF', 'Interlaminar_MS']
    write_sheet(df_formatted, 'Advanced_Criteria', adv_cols)

    # --- SHEET: ABD MATRICES ---
    if df_abd is not None and not df_abd.empty:
        df_abd.to_excel(writer, sheet_name='Laminate_ABD', index=False)
        ws_abd = writer.sheets['Laminate_ABD']
        ws_abd.freeze_panes(1, 0)
        for i, c in enumerate(df_abd.columns): 
            ws_abd.write(0, i, c, fmt_head)
            ws_abd.set_column(i, i, 20, fmt_gen if i < 2 else fmt_num2)

    # --- SHEET: VISUAL GRAPHS ---
    ws_g = wb.add_worksheet('Visual_Graphs')
    ws_g.merge_range('B2:M3', "FAILURE ENVELOPE VISUALIZATION", fmt_title)
    ws_g.hide_gridlines(2)
    
    row_img = 5
    for lc, img_file in img_dict.items():
        ws_g.write(f'B{row_img}', f"📈 RESULTS FOR: {lc}", wb.add_format({'bold':True, 'font_size':14, 'bg_color':'#D9E1F2', 'border':1}))
        if os.path.exists(img_file):
            ws_g.insert_image(f'B{row_img+2}', img_file)
        row_img += 38

    writer.close()
    for img_file in img_dict.values():
        try: os.remove(img_file)
        except: pass

# ==============================================================================
# 7. SMART MANUAL CALCULATOR (MENU FOR THEORIES & RF)
# ==============================================================================
def manual_composite_input():
    print("\n" + "╔════════════════════════════════════════════════════════════╗")
    print("║        🧮 MANUAL COMPOSITE CALCULATOR (SINGLE PLY)         ║")
    print("╚════════════════════════════════════════════════════════════╝")
    
    print("\n  [STEP 0] SELECT UNIT")
    print("  [1] Pascal (Pa) \n  [2] Megapascal (MPa) - Default\n  [3] Gigapascal (GPa)\n  [4] PSI")
    u_choice = input("  Select unit (1/2/3/4): ").strip()
    unit_str = {'1': 'Pa', '2': 'MPa', '3': 'GPa', '4': 'PSI'}.get(u_choice, 'MPa')

    print("\n  Select the failure criteria you want to calculate:")
    print("  [1] Maximum Stress")
    print("  [2] Maximum Strain")
    print("  [3] Tsai-Hill")
    print("  [4] Tsai-Wu")
    print("  [5] Hashin (1980)")
    print("  [6] Puck (Action Plane)")
    print("  [7] Calculate ALL (Compare)")
    
    while True:
        choice = input("\n  -> Your Choice (1-7): ").strip()
        if choice in ['1', '2', '3', '4', '5', '6', '7']: break
        print("  [X] Invalid choice.")

    print("\n  [STEP 1] INPUT MATERIAL STRENGTHS")
    while True:
        try:
            Xt = float(input(f"  Xt (Fiber Tensile Strength)      [{unit_str}]: "))
            Xc = float(input(f"  Xc (Fiber Compressive Strength)  [{unit_str}]: "))
            Yt = float(input(f"  Yt (Matrix Tensile Strength)     [{unit_str}]: "))
            Yc = float(input(f"  Yc (Matrix Compressive Strength) [{unit_str}]: "))
            S  = float(input(f"  S  (In-Plane Shear Strength)     [{unit_str}]: "))
            S13 = float(input(f"  S13 (Interlaminar Shear XZ)      [{unit_str}]: "))
            S23 = float(input(f"  S23 (Interlaminar Shear YZ)      [{unit_str}]: "))
            break
        except ValueError: print("  [X] Please enter numbers only.\n")

    E1 = 1e6; E2 = 1e6; nu12 = 0.3; G12 = 1e6 # Safe defaults for bypass
    if choice in ['2', '7']:
        print("\n  [STEP 1B] INPUT MATERIAL STIFFNESS")
        print("  (Required because you selected Max Strain / ALL)")
        while True:
            try:
                E1   = float(input(f"  E1 (Longitudinal Modulus) [{unit_str}]: "))
                E2   = float(input(f"  E2 (Transverse Modulus)   [{unit_str}]: "))
                nu12 = float(input(f"  nu12 (Poisson's Ratio)          : "))
                G12  = float(input(f"  G12 (Shear Modulus)       [{unit_str}]: "))
                break
            except ValueError: print("  [X] Please enter numbers only.\n")

    print("\n  [STEP 2] INPUT APPLIED STRESSES")
    while True:
        try:
            s1_inp  = float(input(f"  Sigma_11 (Fiber Direction)  [{unit_str}]: "))
            s2_inp  = float(input(f"  Sigma_22 (Matrix Direction) [{unit_str}]: "))
            t12_inp = float(input(f"  Tau_12 (In-Plane Shear)     [{unit_str}]: "))
            t13_inp = float(input(f"  Tau_13 (ILS XZ Direction)   [{unit_str}]: "))
            t23_inp = float(input(f"  Tau_23 (ILS YZ Direction)   [{unit_str}]: "))
            break
        except ValueError: print("  [X] Please enter numbers only.\n")

    mat_props = {'E1': E1, 'E2': E2, 'nu12': nu12, 'G12': G12, 'Xt': Xt, 'Xc': Xc, 'Yt': Yt, 'Yc': Yc, 'S': S, 'S13': S13, 'S23': S23}
    
    # Create a dummy dataframe to pass into the exact same engine used by the batch processor
    dummy_df = pd.DataFrame({'Sigma_11': [s1_inp], 'Sigma_22': [s2_inp], 'Tau_12': [t12_inp], 'Tau_13': [t13_inp], 'Tau_23': [t23_inp]})
    res = analyze_composite_failure(dummy_df, mat_props)
    
    print("\n" + "═"*60)
    print(" 📊 COMPOSITE ANALYSIS RESULTS")
    print("═"*60)
    
    if choice in ['1', '7']:
        print("  [MAXIMUM STRESS]")
        print(f"  > Failure Index    : {res['Max_Stress_FI'].iloc[0]:.4f}")
        print(f"  > Reserve Factor   : {res['Max_Stress_RF'].iloc[0]:.4f}")
        print(f"  > Margin of Safety : {res['Max_Stress_MS'].iloc[0]:.4f}")
        print(f"  > Failure Mode     : {res['Max_Stress_Mode'].iloc[0]}\n")
        
    if choice in ['2', '7']:
        print("  [MAXIMUM STRAIN]")
        print(f"  > Failure Index    : {res['Max_Strain_FI'].iloc[0]:.4f}")
        print(f"  > Reserve Factor   : {res['Max_Strain_RF'].iloc[0]:.4f}")
        print(f"  > Margin of Safety : {res['Max_Strain_MS'].iloc[0]:.4f}\n")
        
    if choice in ['3', '7']:
        print("  [TSAI-HILL]")
        print(f"  > Failure Index    : {res['Tsai_Hill_FI'].iloc[0]:.4f}")
        print(f"  > Reserve Factor   : {res['Tsai_Hill_RF'].iloc[0]:.4f}")
        print(f"  > Margin of Safety : {res['Tsai_Hill_MS'].iloc[0]:.4f}\n")
        
    if choice in ['4', '7']:
        print("  [TSAI-WU]")
        print(f"  > Failure Index    : {res['Tsai_Wu_FI'].iloc[0]:.4f}")
        print(f"  > Reserve Factor   : {res['Tsai_Wu_RF'].iloc[0]:.4f}")
        print(f"  > Margin of Safety : {res['Tsai_Wu_MS'].iloc[0]:.4f}\n")

    if choice in ['5', '7']:
        print("  [HASHIN (1980)]")
        print(f"  > Failure Index    : {res['Hashin_FI'].iloc[0]:.4f}")
        print(f"  > Reserve Factor   : {res['Hashin_RF'].iloc[0]:.4f}")
        print(f"  > Margin of Safety : {res['Hashin_MS'].iloc[0]:.4f}")
        print(f"  > Failure Mode     : {res['Hashin_Mode'].iloc[0]}\n")

    if choice in ['6', '7']:
        print("  [PUCK 2D ACTION PLANE]")
        print(f"  > Failure Index    : {res['Puck_FI'].iloc[0]:.4f}")
        print(f"  > Reserve Factor   : {res['Puck_RF'].iloc[0]:.4f}")
        print(f"  > Margin of Safety : {res['Puck_MS'].iloc[0]:.4f}")
        print(f"  > Failure Mode     : {res['Puck_Mode'].iloc[0]}\n")

    if choice == '7':
        print("  [INTERLAMINAR SHEAR (DELAMINATION)]")
        print(f"  > Failure Index    : {res['Interlaminar_FI'].iloc[0]:.4f}")
        print(f"  > Reserve Factor   : {res['Interlaminar_RF'].iloc[0]:.4f}")
        print(f"  > Margin of Safety : {res['Interlaminar_MS'].iloc[0]:.4f}\n")
        
    print("════════════════════════════════════════════════════════════\n")
    input("Press Enter to return...")

# ==============================================================================
# 8. MAIN UI / CLI APPLICATION FLOW
# ==============================================================================
def main():
    while True:
        os.system('cls' if os.name == 'nt' else 'clear')
        print("╔═════════════════════════════════════════════════════════════════╗")
        print("║      🚀 AEROSPACE COMPOSITE GOD-TIER SUITE (ALL FEATURES)       ║")
        print("╚═════════════════════════════════════════════════════════════════╝")
        print("Select the operation mode you want to use:")
        print("  [1] 📁 Analyze .h5 File (Export Excel, ABD Matrix & Graphs)")
        print("  [2] 🧮 Manual Calculator (Select Theory & Reserve Factor)")
        print("  [0] Exit Application")
        
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
                print(f"  [X] No .h5 files found in directory:\n{script_dir}"); input(); continue

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
            df_raw, pcomp_raw = extract_composite_h5(selected_h5, bdf_file=selected_bdf)
            if df_raw is None or df_raw.empty:
                print("  [X] File is empty or composite stresses not found."); input(); continue

            print("\n[STEP 3] ⚖️ UNIT CONVERSION")
            print("  [1] Pascal (Pa) \n  [2] Megapascal (MPa) - Default\n  [3] Gigapascal (GPa)\n  [4] PSI")
            unit_choice = input("  Select unit (1/2/3/4): ").strip()
            unit_map = {'1': (1e-6, 'Pa'), '2': (1.0, 'MPa'), '3': (1000.0, 'GPa'), '4': (0.00689475, 'PSI')}
            conv_factor, unit_str = unit_map.get(unit_choice, (1.0, 'MPa'))
                
            if conv_factor != 1.0:
                for col in ['Sigma_11', 'Sigma_22', 'Tau_12', 'Tau_13', 'Tau_23']:
                    if col in df_raw.columns: df_raw[col] *= conv_factor

            print("\n[STEP 4] ⚙️ MATERIAL CONFIGURATION (Auto-Detect from .H5)")
            extracted_mats = extract_mat8_h5(selected_h5, conv_factor)
            mat_selected = None
            
            if extracted_mats:
                print(f"  🕵️‍♂️ Nastran found {len(extracted_mats)} MAT8 materials!")
                for mid, props in extracted_mats.items():
                    print(f"  [ Material ID: {mid} ]")
                    print(f"    E1 = {props['E1']:,.1f} {unit_str} | E2 = {props['E2']:,.1f} {unit_str} | G12 = {props['G12']:,.1f} {unit_str} | nu12 = {props['nu12']}")
                    print(f"    Xt = {props['Xt']:,.1f} {unit_str} | Xc = {props['Xc']:,.1f} {unit_str}")
                    print(f"    Yt = {props['Yt']:,.1f} {unit_str} | Yc = {props['Yc']:,.1f} {unit_str}")
                    print(f"    S  = {props['S']:,.1f} {unit_str} | S13= {props['S13']:,.1f} {unit_str} | S23= {props['S23']:,.1f} {unit_str}")
                    print("-" * 55)
                    
                use_auto = input("\n  -> Use automatic material data from this H5? (Y/N): ").strip().upper()
                if use_auto == 'Y':
                    if len(extracted_mats) == 1: mat_selected = list(extracted_mats.values())[0]
                    else:
                        try:
                            mid_sel = int(input("  -> Select Material ID: "))
                            mat_selected = extracted_mats.get(mid_sel, None)
                        except: pass
                                
            if not mat_selected:
                print("\n  --- MANUAL MATERIAL PROPERTY INPUT ---")
                while True:
                    try:
                        mat_selected = {
                            'E1': float(input(f"  E1 (Longitudinal Modulus) [{unit_str}]: ")), 
                            'E2': float(input(f"  E2 (Transverse Modulus)   [{unit_str}]: ")),
                            'nu12': float(input("  nu12 (Poisson's Ratio)          : ")), 
                            'G12': float(input(f"  G12 (Shear Modulus)       [{unit_str}]: ")),
                            'Xt': float(input(f"  Xt (Fiber Tensile Str)    [{unit_str}]: ")), 
                            'Xc': float(input(f"  Xc (Fiber Compressive Str)[{unit_str}]: ")),
                            'Yt': float(input(f"  Yt (Matrix Tensile Str)   [{unit_str}]: ")), 
                            'Yc': float(input(f"  Yc (Matrix Compress Str)  [{unit_str}]: ")),
                            'S': float(input(f"  S (In-Plane Shear Str)    [{unit_str}]: ")),
                            'S13': float(input(f"  S13 (ILS XZ Direction)    [{unit_str}]: ")), 
                            'S23': float(input(f"  S23 (ILS YZ Direction)    [{unit_str}]: "))
                        }
                        break 
                    except ValueError: print("  [X] Invalid input! Numbers only.\n")

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
            
            # --- CALCULATE ALL CRITERIA ---
            df_results = analyze_composite_failure(df_raw, mat_selected)
            
            # --- CALCULATE ABD MATRIX ---
            df_abd = calculate_abd_matrices(pcomp_raw, mat_selected) if pcomp_raw else None
            
            if top_n > 0:
                print(f"  -> Filtering to Top {top_n} worst plies...")
                df_results = df_results.sort_values(['Load_Case', 'Tsai_Wu_MS'], ascending=[True, True])
                df_results = df_results.groupby('Load_Case', observed=True).head(top_n).reset_index(drop=True)
                
            img_dict = create_failure_envelopes(df_results, mat_selected, base_name, unit_str)
            
            save_composite_report_6sheets(df_results, df_abd, mat_selected, output_excel, base_name, img_dict, unit_str)
            
            print("\n" + "═"*65)
            print(f" 🎉 SUCCESS! EXCEL REPORT GENERATED: {os.path.basename(output_excel)}")
            print("═"*65 + "\n")
            
            input("Press Enter to return to Main Menu...")

if __name__ == "__main__":
    main()
