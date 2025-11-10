import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import phase_screen_generator as psg
import concurrent.futures
import threading
from calculating_dne import compute_iri_path

# === Simulation Parameters ===
n_runs = 1000
n = 16384
lx = 30000.0
wave_amp = np.ones(n, dtype=complex)
dx = lx / n

root = None
progress_bar = None
result_text = None
run_button = None
entries = {}
psd_var = None

current_screens = 4
mean_max_dphi = 0.0
dz_list = []
ne_local_list = []
dne_list = []
z_len = 0.0
lambda_wave = None
wave_freq = None  


def run_propagation(dz_list, dne_list, wave_freq):
    grid = psg.init_grid(lx, n)  
    amp = wave_amp.copy()
    max_phase = 0.0

    for i in range(len(dz_list)):
        dz_local = dz_list[i]
        dne_local = dne_list[i]

        phase_screen = psg.generate_phase_screen(dne_local, dz_local, wave_freq, grid)
        current_max = np.max(np.abs(phase_screen))
        if current_max > max_phase:
            max_phase = current_max
        amp = psg.propagate_wave(amp, wave_freq, dz_local, phase_screen, grid)

    intensity = np.abs(amp)**2
    mean_I = np.mean(intensity)
    s4 = np.sqrt(np.mean(intensity**2) - mean_I**2) / mean_I
    return s4, max_phase


def create_gui():
    global root, progress_bar, result_text, run_button, entries, psd_var

    root = tk.Tk()
    root.title("Phase Screen")
    root.geometry("620x900")

    tk.Label(root, text="Parameters", font=("Arial", 16, "bold")).pack(pady=10)

    tk.Label(root, text="PSD Type:").pack(anchor='w', padx=20)
    psd_var = tk.IntVar(value=2)
    frame_psd = tk.Frame(root)
    frame_psd.pack(pady=5)
    for val, text in [(0, "Gaussian"), (1, "von Karman"), (2, "Shkarofsky")]:
        tk.Radiobutton(frame_psd, text=text, variable=psd_var, value=val).pack(side='left', padx=15)

    params = [
        ("Outer Scale [m]:", "5000.0"),
        ("Inner Scale [m]:", "10.0"),
        ("Spectral Index:", "1.5"),
        ("Num Screens:", "4"),
        ("dN/N:", "0.1"),
        ("Frequency [GHz]:", "1"),
    ]
    for label, default in params:
        frame = tk.Frame(root)
        frame.pack(pady=3, padx=20, anchor='w')
        tk.Label(frame, text=label, width=18).pack(side='left')
        entry = tk.Entry(frame, width=12)
        entry.insert(0, default)
        entry.pack(side='left', padx=5)
        entries[label] = entry

    tk.Label(root, text="Source and reciever (lat, lon, alt [km]):", font=("Arial", 12)).pack(pady=(15,5), anchor='w', padx=20)
    traj = [
        ("Src Lat:", "75.0"), ("Src Lon:", "0.0"), ("Src Alt:", "400.0"),
        ("Rcv Lat:", "75.0"), ("Rcv Lon:", "1.0"), ("Rcv Alt:", "200.0")
    ]
    for label, default in traj:
        frame = tk.Frame(root)
        frame.pack(pady=2, padx=20, anchor='w')
        tk.Label(frame, text=label, width=12).pack(side='left')
        entry = tk.Entry(frame, width=10)
        entry.insert(0, default)
        entry.pack(side='left', padx=5)
        entries[label] = entry

    tk.Label(root, text="Date, time, F10.7:", font=("Arial", 12)).pack(pady=(15,5), anchor='w', padx=20)
    time_p = [("Year:", "2020"), ("Month:", "06"), ("Day:", "21"), ("Hour:", "14"), ("F10.7:", "150")]
    for label, default in time_p:
        frame = tk.Frame(root)
        frame.pack(pady=2, padx=20, anchor='w')
        tk.Label(frame, text=label, width=10).pack(side='left')
        entry = tk.Entry(frame, width=8)
        entry.insert(0, default)
        entry.pack(side='left', padx=5)
        entries[label] = entry

    run_button = tk.Button(root, text="RUN SIMULATION", font=("Arial", 12, "bold"), bg="#4CAF50", fg="white", command=start_simulation)
    run_button.pack(pady=20)

    progress_bar = ttk.Progressbar(root, mode='determinate', maximum=n_runs)
    progress_bar.pack(fill='x', padx=20, pady=10)

    result_text = tk.Text(root, height=22, width=70, font=("Courier", 9), state='disabled', bg="#f0f0f0")
    result_text.pack(pady=10, padx=20)

    root.mainloop()

# === Validation ===
def validate_pre(outer, inner):
    L = lx
    cond20 = L > 5 * outer
    cond21 = dx < inner / 3
    errors = []
    if not cond20: errors += [f"cond20: L≤5*L0 ({L}≤{5*outer})", "→ Увеличь L или уменьши L0"]
    if not cond21: errors += [f"cond21: dx≥l0/3 ({dx:.3f}≥{inner/3:.3f})", "→ Увеличь N или уменьши l0"]
    return cond20 and cond21, errors

def validate_post(mean_dphi, dz_list):
    global lambda_wave, z_len
    L = lx
    warnings = validate_cond29(mean_dphi)
    
    max_dz = np.max(dz_list)
    cond27_limit = 2 * L * dx / lambda_wave
    if max_dz >= cond27_limit:
        warnings += [f"cond27: max(dz)={max_dz:.1f}m ≥ 2 L dx/λ ({cond27_limit:.1f}m)", "→ Уменьши dz (увеличь num_screens)"]
    
    return warnings

def validate_cond29(mean_dphi):
    global wave_freq, z_len
    L = lx
    k_wave = 2 * np.pi * wave_freq / psg.C_LIGHT
    mean_max_dphi_dx = mean_dphi / dx
    threshold = (z_len / k_wave) * abs(mean_max_dphi_dx)
    if L <= threshold:
        return [f"cond29: L ({L:.1f}m) ≤ (z/k)*|mean_dphi_dx| ({threshold:.1f}m)", "→ Увеличь L"]
    return []

# === Launch ===
def start_simulation():
    global current_screens, dz_list, ne_local_list, dne_list, z_len, mean_max_dphi, lambda_wave, wave_freq

    try:
        outer = float(entries["Outer Scale [m]:"].get())
        inner = float(entries["Inner Scale [m]:"].get())
        spectral = float(entries["Spectral Index:"].get())
        psd_type = psd_var.get()
        num_screens = int(entries["Num Screens:"].get())
        dN_N = float(entries["dN/N:"].get())
        freq_ghz = float(entries["Frequency [GHz]:"].get())

        src_lat = float(entries["Src Lat:"].get())
        src_lon = float(entries["Src Lon:"].get())
        src_alt = float(entries["Src Alt:"].get()) * 1000
        rcv_lat = float(entries["Rcv Lat:"].get())
        rcv_lon = float(entries["Rcv Lon:"].get())
        rcv_alt = float(entries["Rcv Alt:"].get()) * 1000

        year = int(entries["Year:"].get())
        month = int(entries["Month:"].get())
        day = int(entries["Day:"].get())
        hour = float(entries["Hour:"].get())
        f107 = float(entries["F10.7:"].get())
    

        valid, errors = validate_pre(outer, inner)
        if not valid:
            messagebox.showerror("Error", "\n".join(errors))
            return

        psg.OUTER_SCALE = outer
        psg.INNER_SCALE = inner
        psg.SPECTRAL_INDEX = spectral
        psg.PSD_TYPE = psd_type

        wave_freq = freq_ghz * 1e9  
        lambda_wave = psg.C_LIGHT / wave_freq
        
        z_len, dz_list, ne_local_list, dne_list = compute_iri_path(
            src_lat, src_lon, src_alt, rcv_lat, rcv_lon, rcv_alt,
            year, month, day, hour, num_screens, dN_N, f107
        )
        current_screens = len(dz_list)

        run_button.config(state='disabled', text="RUNNING...")
        progress_bar['value'] = 0
        result_text.config(state='normal')
        result_text.delete(1.0, tk.END)
        result_text.insert(tk.END, f"z_len: {z_len/1000:.1f} km\nScreens: {current_screens} \nStarting...\n")
        result_text.config(state='disabled')

        threading.Thread(target=run_simulation_background, daemon=True).start()

    except Exception as e:
        messagebox.showerror("Error", str(e))

# === Running ===
def run_simulation_background():
    global current_screens, dz_list, ne_local_list, dne_list, mean_max_dphi, z_len, wave_freq

    iter_count = 0
    max_iters = 10  
    while iter_count < max_iters:
        iter_count += 1
        results = []
        completed = 0
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = [executor.submit(run_propagation, dz_list, dne_list, wave_freq) for _ in range(n_runs)]
            for future in concurrent.futures.as_completed(futures):
                s4, max_phase = future.result()
                results.append((s4, max_phase))
                completed += 1
                update_progress(completed)

        results = np.array(results)
        mean_s4 = np.mean(results[:, 0])
        error_s4 = np.std(results[:, 0], ddof=1) / np.sqrt(n_runs)
        mean_max_dphi = np.mean(results[:, 1])

        if mean_max_dphi > 0.25:

            num_screens = current_screens * 2
            msg = f"mean_max_dφ={mean_max_dphi:.3f}>0.25 → screens {current_screens}→{num_screens}\n"
            src_lat = float(entries["Src Lat:"].get())
            src_lon = float(entries["Src Lon:"].get())
            src_alt = float(entries["Src Alt:"].get()) * 1000
            rcv_lat = float(entries["Rcv Lat:"].get())
            rcv_lon = float(entries["Rcv Lon:"].get())
            rcv_alt = float(entries["Rcv Alt:"].get()) * 1000
            year = int(entries["Year:"].get())
            month = int(entries["Month:"].get())
            day = int(entries["Day:"].get())
            hour = float(entries["Hour:"].get())
            dN_N = float(entries["dN/N:"].get())
            f107 = float(entries["F10.7:"].get())

            z_len, dz_list, ne_local_list, dne_list = compute_iri_path(
                src_lat, src_lon, src_alt, rcv_lat, rcv_lon, rcv_alt,
                year, month, day, hour, num_screens, dN_N, f107
            )
            current_screens = len(dz_list)
            
            result_text.config(state='normal')
            result_text.insert(tk.END, msg)
            result_text.config(state='disabled')
            continue
        break

    if iter_count >= max_iters:
        show_result("ERROR: Max iterations reached, simulation not converging.")

    warnings = validate_post(mean_max_dphi, dz_list)
    result_str = f"""SIMULATION COMPLETE

frequency = {wave_freq * 1e-6:.1f} MHz
S4 = {mean_s4:.6f} (δS4 = {error_s4:.6f})
dφ_max (mean) = {mean_max_dphi:.3f} rad

n_runs: {n_runs}
number_of_screens: {current_screens}
z_len: {z_len/1000:.1f} km
max dz: {np.max(dz_list):.1f} m
"""
    if warnings:
        result_str += "\nWARNINGS:\n" + "\n".join(warnings)

    show_result(result_str)

def update_progress(value):
    progress_bar['value'] = value
    root.update_idletasks()

def show_result(text):
    result_text.config(state='normal')
    result_text.delete(1.0, tk.END)
    result_text.insert(tk.END, text)
    result_text.config(state='disabled')
    run_button.config(state='normal', text="RUN SIMULATION")
    progress_bar['value'] = 0
    entries["Num Screens:"].delete(0, tk.END)
    entries["Num Screens:"].insert(0, str(current_screens))

if __name__ == "__main__":
    create_gui()

