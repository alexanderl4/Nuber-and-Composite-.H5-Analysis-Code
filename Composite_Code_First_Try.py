import tables             # Library to read Nastran HDF5 binary database files
import pandas as pd       # Library for data manipulation, table creation, and Excel export
import numpy as np        # Library for high-performance vectorized mathematical operations
import os                 # Library to interact with the Operating System (paths, files)
import glob               # Library to search for files with specific extensions (like .h5)
import re                 # Library for Regular Expressions to clean strings
from tqdm import tqdm     # Library to display progress bars in the terminal
from xlsxwriter.utility import xl_col_to_name # Library to format Excel columns

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

# ==============================================================================
# 2. COMPOSITE FAILURE PHYSICS ENGINE
# ==============================================================================
def analyze_composite_failure(s1, s2, t12, props):
    """
    Core mathematical engine to calculate Failure Indices (FI) and Margins of Safety (MS) 
    for 4 classical composite criteria. Uses numpy arrays for extreme speed.
    """
    E1, E2, nu12, G12 = props['E1'], props['E2'], props['nu12'], props['G12']
    
    # SAFEGUARD: Prevent ZeroDivisionError
    def enforce_limit(val):
        return abs(val) if abs(val) > 1e-6 else 1e12 
        
    Xt, Xc = enforce_limit(props['Xt']), enforce_limit(props['Xc'])
    Yt, Yc = enforce_limit(props['Yt']), enforce_limit(props['Yc'])
    S = enforce_limit(props['S'])

    s1 = np.asarray(s1); s2 = np.asarray(s2); t12 = np.asarray(t12)

    # --- A. MAXIMUM STRESS CRITERION ---
    fi_1_maxs = np.where(s1 > 0, s1 / Xt, np.abs(s1) / Xc)
    fi_2_maxs = np.where(s2 > 0, s2 / Yt, np.abs(s2) / Yc)
    fi_12_maxs = np.abs(t12) / S 
    fi_max_stress = np.maximum.reduce([fi_1_maxs, fi_2_maxs, fi_12_maxs])
    ms_max_stress = np.where(fi_max_stress > 0, (1.0 / fi_max_stress) - 1.0, 999.9)

    modes = []
    for f1, f2, fm, sig1, sig2 in np.nditer([fi_1_maxs, fi_2_maxs, fi_max_stress, s1, s2]):
        if fm == 0: modes.append("Safe (Zero Stress)")
        elif fm == f1 and sig1 > 0: modes.append("Fiber Tension")
        elif fm == f1 and sig1 <= 0: modes.append("Fiber Compression")
        elif fm == f2 and sig2 > 0: modes.append("Matrix Tension")
        elif fm == f2 and sig2 <= 0: modes.append("Matrix Compression")
        else: modes.append("In-Plane Shear")

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

    return {
        'Max_Stress_FI': fi_max_stress, 'Max_Stress_MS': ms_max_stress, 'Critical_Mode': modes,
        'Max_Strain_FI': fi_max_strain, 'Max_Strain_MS': ms_max_strain,
        'Tsai_Hill_FI': fi_tsai_hill, 'Tsai_Hill_MS': ms_tsai_hill,
        'Tsai_Wu_FI': fi_tsai_wu, 'Tsai_Wu_MS': ms_tsai_wu
    }

# ==============================================================================
# 3. H5 DATA MINING (STRESSES & MATERIAL PROPERTIES)
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
                    'S': get_val('S') * conv_factor
                }
        h5.close()
    except Exception: pass
    return mat_dict

def extract_composite_h5(input_file, bdf_file=None):
    print(f"    [+] Extracting Stresses from H5: {os.path.basename(input_file)} ...")
    h5 = tables.open_file(input_file, mode="r") 
    
    bdf_subtitles = parse_bdf_subtitles(bdf_file)
    domain_to_subcase = {} 
    
    try:
        if hasattr(h5.root.NASTRAN.RESULT, 'DOMAINS'):
            domains_table = h5.root.NASTRAN.RESULT.DOMAINS.read()
            cols_dom = domains_table.dtype.names
            for row in domains_table:
                d_id = row['ID']
                sub_id = row['SUBCASE'] if 'SUBCASE' in cols_dom else d_id
                title_str = ""
                for t_col in ['SUBTITLE', 'TITLE', 'LABEL']:
                    if t_col in cols_dom:
                        clean_val = clean_subtitle(row[t_col])
                        if clean_val: title_str = clean_val; break
                domain_to_subcase[d_id] = f"Subcase {sub_id}: {title_str}" if title_str else f"Subcase {sub_id}"
    except Exception: pass

    def get_load_case_name(domain_id): return domain_to_subcase.get(domain_id, f"Subcase {domain_id}")

    try:
        comp_stress_node = h5.root.NASTRAN.RESULT.ELEMENTAL.STRESS.QUAD4_COMP
        data = comp_stress_node.read()
        cols = data.dtype.names
        
        df_raw = pd.DataFrame()
        df_raw['Element ID'] = data['EID']
        if 'PLY' in cols: df_raw['Ply ID'] = data['PLY']
        elif 'LAYER' in cols: df_raw['Ply ID'] = data['LAYER']
        else: df_raw['Ply ID'] = 1 
        
        df_raw['Load_Case'] = [get_load_case_name(d_id) for d_id in data['DOMAIN_ID']] if 'DOMAIN_ID' in cols else "Subcase 1"
        
        if 'X1' in cols and 'Y1' in cols and 'T1' in cols:
            df_raw['Sigma_11'] = data['X1']; df_raw['Sigma_22'] = data['Y1']; df_raw['Tau_12']   = data['T1']
        elif 'X' in cols and 'Y' in cols and 'TXY' in cols:
            df_raw['Sigma_11'] = data['X']; df_raw['Sigma_22'] = data['Y']; df_raw['Tau_12']   = data['TXY']
        elif 'X_NORMAL' in cols and 'Y_NORMAL' in cols and 'XY_SHEAR' in cols:
            df_raw['Sigma_11'] = data['X_NORMAL']; df_raw['Sigma_22'] = data['Y_NORMAL']; df_raw['Tau_12']   = data['XY_SHEAR']
        else:
            print("  [X] Error: Could not find matching stress columns in QUAD4_COMP.")
            h5.close(); return None

        df_raw['Load_Case'] = df_raw['Load_Case'].astype('category')
        
    except Exception as e: 
        print(f"  [X] Fatal extraction error: {e}"); h5.close(); return None
        
    h5.close() 
    return df_raw

# ==============================================================================
# 4. EXCEL EXPORT
# ==============================================================================
def save_composite_report(df_res, output_excel, base_name):
    print(f"\n[STEP 6] 📑 Compiling Composite Report: {os.path.basename(output_excel)}...")
    writer = pd.ExcelWriter(output_excel, engine='xlsxwriter') 
    wb = writer.book
    
    fmt_head = wb.add_format({'bold':True, 'fg_color':'#1F497D', 'font_color':'white', 'border':1, 'align':'center', 'valign':'vcenter'})
    fmt_gen  = wb.add_format({'border':1, 'align':'center'}) 
    fmt_num2 = wb.add_format({'num_format':'0.00', 'border':1, 'align':'center'}) 
    fmt_red  = wb.add_format({'font_color':'#9C0006', 'bg_color':'#FFC7CE', 'border':1, 'align':'center'}) 
    fmt_grn  = wb.add_format({'font_color':'#006100', 'bg_color':'#C6EFCE', 'border':1, 'align':'center'}) 
    
    cols_order = ['Load_Case', 'Element ID', 'Ply ID', 'Sigma_11', 'Sigma_22', 'Tau_12', 
                  'Max_Stress_FI', 'Max_Stress_MS', 'Critical_Mode', 
                  'Max_Strain_FI', 'Max_Strain_MS',
                  'Tsai_Hill_FI', 'Tsai_Hill_MS', 
                  'Tsai_Wu_FI', 'Tsai_Wu_MS']
    df_sorted = df_res[cols_order]

    df_sorted.to_excel(writer, sheet_name='Composite_Analysis', index=False)
    ws = writer.sheets['Composite_Analysis']
    ws.autofilter(0, 0, len(df_sorted), len(df_sorted.columns)-1) 
    ws.freeze_panes(1, 0) 
    
    for i, col in enumerate(cols_order):
        ws.write(0, i, col, fmt_head)
        if col == 'Load_Case': ws.set_column(i, i, 25, fmt_gen)
        elif 'Sigma' in col or 'Tau' in col: ws.set_column(i, i, 12, fmt_num2)
        elif 'Mode' in col: ws.set_column(i, i, 20, fmt_gen)
        elif 'FI' in col or 'MS' in col: ws.set_column(i, i, 15, fmt_num2)
        else: ws.set_column(i, i, 12, fmt_gen)

    last_row = len(df_sorted) + 1
    for ms_col_name in ['Max_Stress_MS', 'Max_Strain_MS', 'Tsai_Hill_MS', 'Tsai_Wu_MS']:
        try:
            col_letter = xl_col_to_name(df_sorted.columns.get_loc(ms_col_name))
            ws.conditional_format(f'{col_letter}2:{col_letter}{last_row}', {'type':'cell', 'criteria':'<', 'value':0, 'format':fmt_red})
            ws.conditional_format(f'{col_letter}2:{col_letter}{last_row}', {'type':'cell', 'criteria':'>=', 'value':0, 'format':fmt_grn})
        except: pass

    writer.close()

# ==============================================================================
# 5. SMART MANUAL CALCULATOR (MODE 2)
# ==============================================================================
def manual_composite_input():
    print("\n" + "╔════════════════════════════════════════════════════════════╗")
    print("║        🧮 KALKULATOR MANUAL KOMPOSIT (SINGLE PLY)          ║")
    print("╚════════════════════════════════════════════════════════════╝")
    print("  Pilih kriteria kegagalan yang ingin Anda hitung:")
    print("  [1] Maximum Stress")
    print("  [2] Maximum Strain")
    print("  [3] Tsai-Hill")
    print("  [4] Tsai-Wu")
    print("  [5] Hitung SEMUANYA (Bandingkan)")
    
    while True:
        choice = input("\n  -> Pilihan Anda (1/2/3/4/5): ").strip()
        if choice in ['1', '2', '3', '4', '5']: break
        print("  [X] Pilihan tidak valid.")

    print("\n  [STEP 1] INPUT KEKUATAN MATERIAL (STRENGTHS)")
    print("  (Satuan bebas, asalkan konsisten. Misal: MPa)")
    while True:
        try:
            Xt = float(input("  Xt (Kekuatan Tarik Serat)  : "))
            Xc = float(input("  Xc (Kekuatan Tekan Serat)  : "))
            Yt = float(input("  Yt (Kekuatan Tarik Matriks): "))
            Yc = float(input("  Yc (Kekuatan Tekan Matriks): "))
            S  = float(input("  S  (Kekuatan Geser / Shear): "))
            break
        except ValueError:
            print("  [X] Harap masukkan angka saja.\n")

    # Inisialisasi Stiffness Default (Akan dipakai jika tidak butuh Max Strain)
    E1 = 1e6; E2 = 1e6; nu12 = 0.3; G12 = 1e6

    if choice in ['2', '5']:
        print("\n  [STEP 1B] INPUT KEKAKUAN MATERIAL (STIFFNESS)")
        print("  (Karena Anda memilih Max Strain / All, data Modulus diperlukan)")
        while True:
            try:
                E1   = float(input("  E1 (Modulus Searah Serat) : "))
                E2   = float(input("  E2 (Modulus Tegak Lurus)  : "))
                nu12 = float(input("  nu12 (Poisson Ratio)      : "))
                G12  = float(input("  G12 (Shear Modulus)       : "))
                break
            except ValueError:
                print("  [X] Harap masukkan angka saja.\n")

    print("\n  [STEP 2] INPUT TEGANGAN YANG BEKERJA (APPLIED STRESSES)")
    while True:
        try:
            s1_inp  = float(input("  Sigma_11 (Arah Serat)   : "))
            s2_inp  = float(input("  Sigma_22 (Arah Matriks) : "))
            t12_inp = float(input("  Tau_12 (Tegangan Geser) : "))
            break
        except ValueError:
            print("  [X] Harap masukkan angka saja.\n")

    # Bundle ke dalam format yang dibaca oleh Engine
    mat_props = {'E1': E1, 'E2': E2, 'nu12': nu12, 'G12': G12, 'Xt': Xt, 'Xc': Xc, 'Yt': Yt, 'Yc': Yc, 'S': S}
    
    # Hitung
    res = analyze_composite_failure([s1_inp], [s2_inp], [t12_inp], mat_props)
    
    print("\n" + "═"*60)
    print(" 📊 HASIL ANALISIS KOMPOSIT")
    print("═"*60)
    print(f"  Tegangan: S1 = {s1_inp}, S2 = {s2_inp}, T12 = {t12_inp}\n")
    
    if choice in ['1', '5']:
        print("  [MAXIMUM STRESS]")
        print(f"  > Failure Index    : {res['Max_Stress_FI'][0]:.4f}")
        print(f"  > Margin of Safety : {res['Max_Stress_MS'][0]:.4f}")
        print(f"  > Mode Kegagalan   : {res['Critical_Mode'][0]}\n")
        
    if choice in ['2', '5']:
        print("  [MAXIMUM STRAIN]")
        print(f"  > Failure Index    : {res['Max_Strain_FI'][0]:.4f}")
        print(f"  > Margin of Safety : {res['Max_Strain_MS'][0]:.4f}\n")
        
    if choice in ['3', '5']:
        print("  [TSAI-HILL]")
        print(f"  > Failure Index    : {res['Tsai_Hill_FI'][0]:.4f}")
        print(f"  > Margin of Safety : {res['Tsai_Hill_MS'][0]:.4f}\n")
        
    if choice in ['4', '5']:
        print("  [TSAI-WU]")
        print(f"  > Failure Index    : {res['Tsai_Wu_FI'][0]:.4f}")
        print(f"  > Margin of Safety : {res['Tsai_Wu_MS'][0]:.4f}\n")
        
    print("════════════════════════════════════════════════════════════\n")
    input("Tekan Enter untuk kembali...")

# ==============================================================================
# 6. MAIN UI / CLI APPLICATION FLOW
# ==============================================================================
def main():
    while True:
        os.system('cls' if os.name == 'nt' else 'clear')
        print("╔═════════════════════════════════════════════════════════════════╗")
        print("║     🚀 AEROSPACE COMPOSITE FAILURE ANALYSIS (ULTIMATE UI)       ║")
        print("╚═════════════════════════════════════════════════════════════════╝")
        print("Select the operation mode you want to use:")
        print("  [1] 📁 Analyze .h5 File from Nastran (Full Batch Processing)")
        print("  [2] 🧮 Manual Calculator (Smart Input Menu)")
        print("  [0] Keluar Aplikasi")
        
        mode = input("\n  -> Enter choice (0/1/2): ").strip()
        
        if mode == '0':
            break
        elif mode == '2':
            manual_composite_input()
        elif mode == '1':
            script_dir = os.path.dirname(os.path.abspath(__file__))
            
            # --- [STEP 1] DATA SELECTION ---
            print("\n[STEP 1] 📁 SELECT NASTRAN H5 FILE")
            h5_files = glob.glob(os.path.join(script_dir, "*.h5")) 
            if not h5_files:
                print(f"  [X] No .h5 files found in folder:\n{script_dir}"); return

            if len(h5_files) == 1:
                selected_h5 = h5_files[0]
                print(f"  -> Auto-selected the only H5 file: {os.path.basename(selected_h5)}")
            else:
                for idx, file in enumerate(h5_files): print(f"  [{idx + 1}] {os.path.basename(file)}")
                try: choice = int(input("  -> Select H5 file number: ").strip()); selected_h5 = h5_files[choice - 1]
                except: print("  [X] Invalid input. Aborting."); return

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

            # --- [STEP 2] PRE-SCAN STRESSES ---
            print("\n[STEP 2] 🔍 PRE-SCANNING MODEL")
            df_raw = extract_composite_h5(selected_h5, bdf_file=selected_bdf)
            if df_raw is None or df_raw.empty:
                print("  [X] File is empty or composite stresses not found."); return
            print(f"  -> Successfully loaded {len(df_raw)} ply stress records.")

            # --- [STEP 3] UNIT CONVERSION ---
            print("\n[STEP 3] ⚖️ UNIT CONVERSION")
            print("  What unit was exported by Nastran? (This applies to Stresses & Materials)")
            print("  [1] Pascal (Pa)")
            print("  [2] Megapascal (MPa) - Default")
            print("  [3] Gigapascal (GPa)")
            print("  [4] PSI (Pounds per Square Inch)")
            
            while True:
                unit_choice = input("  Select unit (1/2/3/4): ").strip()
                if unit_choice == '1': conv_factor = 1e-6; break
                elif unit_choice == '2': conv_factor = 1.0; break
                elif unit_choice == '3': conv_factor = 1000.0; break
                elif unit_choice == '4': conv_factor = 0.00689475729; break
                elif unit_choice == '': conv_factor = 1.0; break
                else: print("  [!] Invalid choice.")
                
            if conv_factor != 1.0:
                df_raw['Sigma_11'] *= conv_factor
                df_raw['Sigma_22'] *= conv_factor
                df_raw['Tau_12'] *= conv_factor

            # --- [STEP 4] MATERIAL CONFIGURATION (SMART H5 EXTRACTION) ---
            print("\n[STEP 4] ⚙️ MATERIAL CONFIGURATION (Auto-Detect from .H5)")
            extracted_mats = extract_mat8_h5(selected_h5, conv_factor)
            mat_selected = None
            
            if extracted_mats:
                print(f"  🕵️‍♂️ Detektif Nastran menemukan {len(extracted_mats)} material MAT8 di dalam file Anda!")
                print("  (Nilai di bawah ini sudah dikalikan dengan konversi unit pilihan Anda)\n")
                
                for mid, props in extracted_mats.items():
                    print(f"  [ Material ID: {mid} ]")
                    print(f"    E1 = {props['E1']:,.1f} | E2 = {props['E2']:,.1f} | G12 = {props['G12']:,.1f} | nu12 = {props['nu12']}")
                    print(f"    Xt = {props['Xt']:,.1f} | Xc = {props['Xc']:,.1f}")
                    print(f"    Yt = {props['Yt']:,.1f} | Yc = {props['Yc']:,.1f}")
                    print(f"    S  = {props['S']:,.1f}")
                    print("-" * 55)
                    
                use_auto = input("\n  -> Gunakan data material dari H5 ini? (Y untuk Ya / N untuk Input Manual): ").strip().upper()
                
                if use_auto == 'Y':
                    if len(extracted_mats) == 1:
                        mat_selected = list(extracted_mats.values())[0]
                        print("  -> [OK] Material otomatis diaplikasikan!")
                    else:
                        while True:
                            try:
                                chosen_mid = int(input("  -> Masukkan Material ID yang ingin digunakan: "))
                                if chosen_mid in extracted_mats:
                                    mat_selected = extracted_mats[chosen_mid]
                                    print(f"  -> [OK] Material ID {chosen_mid} diaplikasikan!")
                                    break
                                else:
                                    print("  [X] ID tidak ditemukan dalam daftar.")
                            except ValueError:
                                print("  [X] Harap masukkan angka.")
                                
            if not mat_selected:
                print("\n  --- INPUT MANUAL PROPERTI MATERIAL ---")
                while True:
                    try:
                        mat_selected = {
                            'E1': float(input("  E1 (Longitudinal Modulus)  : ")),
                            'E2': float(input("  E2 (Transverse Modulus)    : ")),
                            'nu12': float(input("  nu12 (Major Poisson Ratio) : ")),
                            'G12': float(input("  G12 (In-Plane Shear Mod)   : ")),
                            'Xt': float(input("  Xt (Fiber Tensile Str)     : ")),
                            'Xc': float(input("  Xc (Fiber Compressive Str) : ")),
                            'Yt': float(input("  Yt (Matrix Tensile Str)    : ")),
                            'Yc': float(input("  Yc (Matrix Compress Str)   : ")),
                            'S': float(input("  S (In-Plane Shear Str)     : "))
                        }
                        break 
                    except ValueError:
                        print("  [X] Input tidak valid! Harap hanya masukkan angka.\n")

            # --- [STEP 5] OUTPUT FILTERING ---
            print("\n[STEP 5] 📊 OUTPUT FILTERING")
            print("  Do you want to export ALL plies, or only the most critical ones?")
            try:
                top_n_inp = input("  -> Export Top 'N' worst plies per load case? (0 for ALL): ").strip()
                top_n = int(top_n_inp) if top_n_inp else 0
            except: top_n = 0
                
            print("\n" + "═"*65)
            print(" 🚀 INITIATING COMPOSITE ENGINE...")
            
            base_name = os.path.splitext(os.path.basename(selected_h5))[0]
            output_excel = os.path.join(script_dir, f"Composite_Report_{base_name}.xlsx")
            try: open(output_excel, 'a').close()
            except: print(f"\n[ERROR] Close the file {os.path.basename(output_excel)} first!"); return
            
            print("  -> Crunching Failure Criteria Arrays...")
            res_dict = analyze_composite_failure(df_raw['Sigma_11'], df_raw['Sigma_22'], df_raw['Tau_12'], mat_selected)
            
            df_results = df_raw.copy()
            for key, val in res_dict.items():
                df_results[key] = val
                
            if top_n > 0:
                print(f"  -> Filtering data to Top {top_n} worst plies based on Tsai-Wu FI...")
                df_results = df_results.sort_values(['Load_Case', 'Tsai_Wu_FI'], ascending=[True, False])
                df_results = df_results.groupby('Load_Case', observed=True).head(top_n).reset_index(drop=True)
                
            save_composite_report(df_results, output_excel, base_name)
            
            print("\n" + "═"*65)
            print(f" 🎉 SUCCESS! EXCEL GENERATED: {os.path.basename(output_excel)}")
            print("═"*65 + "\n")
            
            input("Tekan Enter untuk kembali ke Menu Utama...")

if __name__ == "__main__":
    main()
