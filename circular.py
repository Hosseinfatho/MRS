import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pywt
from scipy.interpolate import interp1d
import matplotlib.cm as cm

# Display settings for better visualization
plt.style.use('default')
plt.rcParams['figure.figsize'] = (18, 18)
plt.rcParams['font.size'] = 10

# === Load and preprocess MRS data ===
spectra = pd.read_csv("./data/spectra.csv", header=None)
chemical_shifts = pd.read_csv("./data/chemical_shifts.csv", header=None).values.flatten()

# === Select spectra to visualize ===
selected_indices = [0, 1, 2, 3]
colors = ['#3288bd', '#66c2a5', '#fdae61', '#f46d43']

# === Align spectra by PCr peak (around 0 ppm) ===
sort_idx = np.argsort(chemical_shifts)[::-1]
ppm = chemical_shifts[sort_idx]
spectra = spectra.iloc[sort_idx]

ref_signal = spectra.iloc[:, selected_indices[0]].values
pc_range = (ppm > -1) & (ppm < 1)
pcr_index = np.argmax(ref_signal[pc_range])
pcr_ppm = ppm[pc_range][pcr_index]
aligned_ppm = ppm - pcr_ppm
ppm_mask = (aligned_ppm >= -20) & (aligned_ppm <= 20)
aligned_ppm = aligned_ppm[ppm_mask]

# === Denoising using wavelet with BayesShrink ===
def denoise_signal(signal, wavelet='db4'):
    level = pywt.dwt_max_level(len(signal), pywt.Wavelet(wavelet).dec_len)
    coeffs = pywt.wavedec(signal, wavelet, level=level)
    sigma = np.median(np.abs(coeffs[-1])) / 0.6745

    def bayes_shrink_threshold(cD, sigma):
        var_signal = np.var(cD)
        return sigma**2 / np.sqrt(max(var_signal - sigma**2, 1e-10))

    denoised_coeffs = [coeffs[0]]
    for i in range(1, len(coeffs)):
        threshold = bayes_shrink_threshold(coeffs[i], sigma)
        denoised_coeffs.append(pywt.threshold(coeffs[i], value=threshold, mode='soft'))

    return pywt.waverec(denoised_coeffs, wavelet)[:len(signal)]

# === Denoise and normalize selected spectra ===
processed = []
for i, idx in enumerate(selected_indices):
    raw = spectra.iloc[:, idx].values
    denoised = denoise_signal(raw)[sort_idx][ppm_mask]
    normed = (denoised - np.min(denoised)) / (np.max(denoised) - np.min(denoised))
    processed.append((normed, denoised, colors[i], f"Spectrum {idx+1}"))

R = 9  # base radius of the circular plot (1.5x larger)

# === Helper functions for polar plotting ===
def ppm_to_angle(ppm_val):
    return np.pi / 2 - 2 * np.pi * ppm_val / 40

def radial_stretch(ppm_val, signal_val, base_radius=R, scale_factor=4, compress_factor=0.4):
    if -20 <= ppm_val <= 20:
        return base_radius + signal_val * scale_factor
    else:
        return base_radius + signal_val * scale_factor * compress_factor

# === Target metabolite peaks to track ===
radar_targets = {
    "ATPβ": -16.15, "NAD+": -8.31, "NADH": -8.13, "ATPα": -7.56,
    "ATPγ": -2.53, "PCr:0""PME": 3.5, "Pi": 4.82, "PDE": 5.5
}

# === Compute CRLB, SNR, and FWHM metrics at each peak ===
def compute_quality_metrics(x, y, signal, ppm_val):
    window_mask = (x > ppm_val - 0.3) & (x < ppm_val + 0.3)
    if np.sum(window_mask) < 5:
        return ('red', 'red', 'red')  # insufficient data

    x_win = x[window_mask]
    y_win = y[window_mask]

    # CRLB (Coefficient of variation)
    crlb_val = np.std(y_win) / max(np.mean(y_win), 1e-6)
    crlb_color = 'green' if crlb_val < 0.2 else 'yellow' if crlb_val < 0.5 else 'red'

    # SNR estimation
    noise_mask = (x > 9)
    noise_std = np.std(signal[noise_mask])
    snr_val = (max(y_win) - np.min(y_win)) / max(noise_std, 1e-6)
    snr_color = 'green' if snr_val > 10 else 'yellow' if snr_val > 5 else 'red'

    # FWHM (Full width at half maximum)
    half_max = (np.max(y_win) + np.min(y_win)) / 2
    interp = interp1d(x_win, y_win - half_max, kind='linear', fill_value="extrapolate")
    try:
        roots = np.where(np.diff(np.sign(interp(x_win))))[0]
        fwhm_val = abs(x_win[roots[-1]] - x_win[roots[0]]) if len(roots) >= 2 else 999
    except:
        fwhm_val = 999
    fwhm_color = 'green' if fwhm_val < 0.1 else 'yellow' if fwhm_val < 0.2 else 'red'

    return (crlb_color, snr_color, fwhm_color)

# === Evaluate quality metrics across spectra ===
quality_metrics = []
for i, (normed, denoised, color, label) in enumerate(processed):
    qm = {}
    for name, ppm_val in radar_targets.items():
        if -20 <= ppm_val <= 20:
            crlb_c, snr_c, fwhm_c = compute_quality_metrics(aligned_ppm, denoised, denoised, ppm_val)
            qm[name] = [crlb_c == 'green', snr_c == 'green', fwhm_c == 'green']
    quality_metrics.append(qm)

# === Set up circular plot ===
fig, ax = plt.subplots(figsize=(18, 18))
ax.set_aspect('equal')
ax.axis('off')
ax.set_xlim(-2.5*R, 2.5*R)
ax.set_ylim(-2.5*R, 2.5*R)
ax.add_patch(Circle((0, 0), R, edgecolor='black', facecolor='none', lw=2))

# === Draw circular bands ===
# === Draw circular bands (5 orbits total inside main circle) ===
n_spectra = len(processed)
total_orbits = n_spectra + 1  # 4 spectrum orbits + 1 extra orbit = 5 total

# Calculate orbit radii to ensure extra orbit is inside main circle with 3 units gap
# Main circle radius = R, we want extra orbit at R-3
extra_orbit_radius = R - 3  # 3 units inside main circle (1.5x larger gap)

# Calculate spacing for spectrum orbits (1, 2, 3, 4)
spectrum_orbit_radii = []
for i in range(1, n_spectra + 1):
    # Distribute spectrum orbits evenly between 0 and extra_orbit_radius
    radius = (i / (n_spectra + 1)) * extra_orbit_radius
    spectrum_orbit_radii.append(radius)

for i in range(1, total_orbits + 1):
    if i <= n_spectra:
        # Spectrum orbits (gray) - evenly spaced
        radius = spectrum_orbit_radii[i-1]
        ax.add_patch(Circle((0, 0), radius=radius, edgecolor='gray', facecolor='none',
                            linestyle='--', linewidth=1.0, alpha=0.6))
    else:
        # Extra orbit (blue) - 2 units inside main circle
        ax.add_patch(Circle((0, 0), radius=extra_orbit_radius, edgecolor='blue', facecolor='none',
                            linestyle='--', linewidth=1.5, alpha=0.7))


# === Plot each denoised spectrum in polar coordinates ===
margin = 0.015  # 1.5x larger
step = 0.75  # 1.5x larger
angles = ppm_to_angle(aligned_ppm)

for i, (norm, denoised, color, label) in enumerate(processed):
    base_r = R + margin + i * step
    radii = [radial_stretch(p, s, base_radius=base_r) for p, s in zip(aligned_ppm, norm)]
    x = [r * np.cos(a) for r, a in zip(radii, angles)]
    y = [r * np.sin(a) for r, a in zip(radii, angles)]
    ax.plot(x, y, color=color, lw=1.5, label=label)

# === Annotate known metabolite peaks ===
metabolites = {
    'ATPβ': -16.15, 'NAD+': -8.31, 'NADH': -8.13, 'ATPα': -7.56, 'ATPγ': -2.53,
    'PCr': 0, 'MP': 2.30, 'GPC': 2.94, 'GPE': 3.50, 'Pi': 4.82, 'DPG2': 5.23,
    'DPG1': 5.71, 'PC': 6.23, 'PE': 6.77
}
for name, ppm_val in metabolites.items():
    if -20 <= ppm_val <= 20:
        angle = ppm_to_angle(ppm_val)
        r = radial_stretch(ppm_val, 1.2, R, 1.5)
        x0, y0 = r * np.cos(angle), r * np.sin(angle)
        x1, y1 = (r - 3.75) * np.cos(angle), (r - 3.75) * np.sin(angle)  # 1.5x larger
        ax.plot([x0, x1], [y0, y1], color='red', lw=1.5)  # 1.5x thicker
        ax.text(x1, y1, name, ha='center', va='center', fontsize=9)  # 1.5x larger font

# === Draw ppm radial ticks ===
ppm_ticks = np.arange(-20, 21, 5)
for ppm_val in ppm_ticks:
    angle = ppm_to_angle(ppm_val)
    r_outer = R + 0.6  # 1.5x larger
    r_inner = R - 0.15  # 1.5x larger
    x_outer, y_outer = r_outer * np.cos(angle), r_outer * np.sin(angle)
    x_inner, y_inner = r_inner * np.cos(angle), r_inner * np.sin(angle)
    ax.plot([x_inner, x_outer], [y_inner, y_outer], color='red', lw=1.2)  # 1.5x thicker
    ax.text(R * np.cos(angle) + 0.3, R * np.sin(angle) + 0.3, f'{ppm_val}', fontsize=15, ha='center', va='center')  # 1.5x larger font

# === Draw biomarker ratios as partial circles based on quality metrics passed ===
factor = 1.0
for i, (norm, denoised, color, label) in enumerate(processed):
    interp_fn = interp1d(aligned_ppm, denoised, kind='linear', fill_value="extrapolate")
    pcr_intensity = interp_fn(0.0)
    for name, ppm_val in radar_targets.items():
        if -20 <= ppm_val <= 20:
            met_intensity = interp_fn(ppm_val)
            ratio = met_intensity / pcr_intensity if pcr_intensity != 0 else 0
            angle = ppm_to_angle(ppm_val)
            # Use spectrum orbit radius for proper alignment
            orbit_radius = spectrum_orbit_radii[i]
            cx = orbit_radius * np.cos(angle)
            cy = orbit_radius * np.sin(angle)
            passes = sum(quality_metrics[i].get(name, [False, False, False]))
            
            # Determine circle fraction based on quality metrics passed
            if passes == 0:
                # No quality metrics passed - draw 1/4 circle
                circle_fraction = 0.25
                start_angle = 0
                end_angle = 90
            elif passes == 1:
                # One quality metric passed - draw 1/2 circle
                circle_fraction = 0.5
                start_angle = 0
                end_angle = 180
            elif passes == 2:
                # Two quality metrics passed - draw 3/4 circle
                circle_fraction = 0.75
                start_angle = 0
                end_angle = 270
            else:
                # All quality metrics passed - draw full circle
                circle_fraction = 1.0
                start_angle = 0
                end_angle = 360
            
            # Create partial circle using Wedge
            from matplotlib.patches import Wedge
            wedge = Wedge((cx, cy), r=ratio * factor, theta1=start_angle, theta2=end_angle,
                         edgecolor=color, facecolor=color, lw=1.2)
            ax.add_patch(wedge)

# === Side legends: spectra + quality + ratio guide ===
# Spectra color legend
spectra_legend_ax = fig.add_axes([0.85, 0.75, 0.12, 0.15])
spectra_legend_ax.axis('off')
for i, col in enumerate(colors[:n_spectra], start=1):
    y = 0.65- i * 0.2  # Increased spacing by 20%
    spectra_legend_ax.add_patch(plt.Circle((0.15, y), +0.05, color=col))
    spectra_legend_ax.text(0.25, y - 0.005, f"Spectrum {i}", fontsize=8, verticalalignment='center')
# spectra_legend_ax.text(0.0, 0.5, "Spectra", fontsize=12, fontweight='bold')

# Quality legend showing circle fractions
quality_ax = fig.add_axes([0.85, 0.55, 0.12, 0.12])
quality_ax.axis('off')
quality_ax.text(0, 0.9, "Quality Metrics", fontsize=8, fontweight='bold')

# Draw example partial circles for quality legend
quality_examples = [
    (0.25, "0 passed", 0, 90),
    (0.5, "1 passed", 0, 180), 
    (0.75, "2 passed", 0, 270),
    (1.0, "3 passed", 0, 360)
]

for i, (fraction, label, start, end) in enumerate(quality_examples):
    y = 0.7 - i * 0.18  # Increased spacing by 20%
    # Draw example wedge
    wedge = Wedge((0.15, y), r=0.08, theta1=start, theta2=end,
                  edgecolor='gray', facecolor='gray', lw=1)
    quality_ax.add_patch(wedge)
    quality_ax.text(0.3, y, label, fontsize=6, verticalalignment='center')

# Ratio guide
ratio_ax = fig.add_axes([0.85, 0.35, 0.12, 0.12])
ratio_ax.axis('off')
for i, r in enumerate([0.25, 0.5, 0.75]):
    x = 0.2 + i * 0.3
    ratio_ax.add_patch(plt.Circle((x, 0.5), 0.1 * r, color='gray', alpha=0.6))
    # ratio_ax.text(x, 0.28, f'{int(r*100)}%PCr', fontsize=10, ha='center', va='center')
ratio_ax.text(0,0.75, "Ratio to PCr", fontsize=6, fontweight='bold')

# === Quadrant-style glyphs for all radar targets / PCr ratios on extra orbit ===
# Extra orbit radius is 3 units inside main circle (1.5x larger gap)
extra_orbit_radius = R - 3

# Create quadrant plots for all radar targets (only those <= +5 ppm)
target_items = list(radar_targets.items())
for target_name, target_ppm in target_items:
    if -20 <= target_ppm <= 5:  # Only for targets <= +5 ppm
        # Calculate ratios for this target
        target_ratios = []
        for norm, denoised, color, label in processed:
            interp_fn = interp1d(aligned_ppm, denoised, kind='linear', fill_value="extrapolate")
            target_intensity = interp_fn(target_ppm)
            pcr = interp_fn(0.0)
            ratio = target_intensity / pcr if pcr != 0 else 0
            target_ratios.append((ratio, color))

        # Normalize ratios for better visibility
        raw_ratios = [r[0] for r in target_ratios]
        max_ratio = max(raw_ratios) if max(raw_ratios) > 0 else 1
        normalized_ratios = [r / max_ratio * 0.6 for r in raw_ratios]  # 1.5x larger size

        # Position on extra orbit at target angle
        target_angle = ppm_to_angle(target_ppm)
        
        # Calculate base position on extra orbit
        base_x = extra_orbit_radius * np.cos(target_angle)
        base_y = extra_orbit_radius * np.sin(target_angle)

        # Create quadrant circle with 4 sections
        n_spectra = len(processed)
        section_angle = 360 / n_spectra  # 90 degrees for 4 spectra
        
        for i, (ratio, color) in enumerate(zip(normalized_ratios, colors)):
            start_angle = i * section_angle
            end_angle = (i + 1) * section_angle
            
            # Create wedge for this spectrum's section
            wedge = Wedge((base_x, base_y), r=ratio, theta1=start_angle, theta2=end_angle,
                         facecolor=color, edgecolor='black', lw=0.5, alpha=0.85)
            ax.add_patch(wedge)
        
        # Add metabolite label in the center (removed)
        # ax.text(base_x, base_y, target_name, fontsize=6, ha='center', va='center', 
        #        fontweight='bold', color='black')

plt.tight_layout()

# For Jupyter notebooks, you can also use:
# plt.show()
# Or simply display the figure without plt.show() in Jupyter
plt.show()

print("Plot completed! Now you have 5 orbits inside the main circle:")
print(f"- {n_spectra} spectrum orbits (gray dashed lines)")
print(f"- 1 extra orbit (blue dashed line)")
print(f"Total: {total_orbits} orbits inside the main circle")
