import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pywt
from scipy.interpolate import interp1d

# --- Load and preprocess data ---
spectra = pd.read_csv("./data/spectra.csv", header=None)
chemical_shifts = pd.read_csv("./data/chemical_shifts.csv", header=None).values.flatten()

selected_indices = [0, 1, 2, 3]  # spectra to plot
colors = ['#3288bd', '#66c2a5', '#fdae61', '#f46d43']

# Sort ppm descending (common MRS style)
sort_idx = np.argsort(chemical_shifts)[::-1]
ppm = chemical_shifts[sort_idx]
spectra = spectra.iloc[sort_idx]

# Align ppm by PCr peak
ref_signal = spectra.iloc[:, selected_indices[0]].values
pc_range = (ppm > -1) & (ppm < 1)
pcr_index = np.argmax(ref_signal[pc_range])
pcr_ppm = ppm[pc_range][pcr_index]
aligned_ppm = ppm - pcr_ppm
ppm_mask = (aligned_ppm >= -18) & (aligned_ppm <= 8)
aligned_ppm = aligned_ppm[ppm_mask]

# --- Denoising ---
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

processed = []
raw_denoised_signals = []
for i, idx in enumerate(selected_indices):
    raw = spectra.iloc[:, idx].values
    denoised = denoise_signal(raw)[sort_idx][ppm_mask]
    normed = (denoised - np.min(denoised)) / (np.max(denoised) - np.min(denoised))
    processed.append((normed, denoised, colors[i], f"Spectrum {idx+1}"))
    raw_denoised_signals.append(denoised)

# --- Metabolite targets ---
radar_targets = {
    "ATPβ": -16.15, "NAD+": -8.31, "NADH": -8.13, "ATPα": -7.56,
    "ATPγ": -2.53,  "PME": 3.5, "Pi": 4.82, "PDE": 5.5
}

metabolites = list(radar_targets.keys())
met_ppms = np.array([radar_targets[m] for m in metabolites])
sorted_idx = np.argsort(met_ppms)[::-1]
metabolites = [metabolites[i] for i in sorted_idx]
met_ppms = met_ppms[sorted_idx]

print(f"Original metabolites: {metabolites}")
print(f"Original met_ppms: {met_ppms}")
print(f"Metabolite centers for transformation: {metabolite_centers}")

ratios = {}
for i, signal in enumerate(raw_denoised_signals):
    spectrum_ratios = {}
    
    # Integrate over ±2 indices (5 points total) around PCr at 0.0 ppm
    pcr_idx = np.argmin(np.abs(aligned_ppm - 0.0))
    pcr_start = max(pcr_idx - 2, 0)
    pcr_end = min(pcr_idx + 3, len(aligned_ppm))  # +3 because slicing is exclusive
    pcr_auc = np.trapz(signal[pcr_start:pcr_end], x=aligned_ppm[pcr_start:pcr_end])

    for met in metabolites:
        met_ppm = radar_targets[met]
        met_idx = np.argmin(np.abs(aligned_ppm - met_ppm))
        start = max(met_idx - 2, 0)
        end = min(met_idx + 3, len(aligned_ppm))
        met_auc = np.trapz(signal[start:end], x=aligned_ppm[start:end])
        
        ratio = (met_auc / pcr_auc if pcr_auc != 0 else 0) * 3  # scaled for visibility
        spectrum_ratios[met] = ratio

    ratios[i] = spectrum_ratios

# --- Non-uniform ppm transformation for better metabolite visualization ---
def transform_ppm(ppm_vals, metabolite_centers, stretch_factor=3.0, compress_factor=0.4):
    """
    Transform ppm values to stretch around specific metabolites and compress elsewhere
    Focus on radar_targets: ATPβ(-16.15), NAD+(-8.31), NADH(-8.13), ATPα(-7.56), 
    ATPγ(-2.53), PME(3.5), Pi(4.82), PDE(5.5)
    """
    transformed = np.copy(ppm_vals).astype(float)
    
    # Create mask for regions near metabolites (stretch regions)
    stretch_mask = np.zeros_like(ppm_vals, dtype=bool)
    
    # Define stretch regions around each metabolite
    for center in metabolite_centers:
        stretch_region = 0.8  # ppm range around each metabolite to stretch
        metabolite_mask = np.abs(ppm_vals - center) <= stretch_region
        stretch_mask |= metabolite_mask
        
        if np.any(metabolite_mask):
            # Stretch around metabolite
            distance = ppm_vals[metabolite_mask] - center
            transformed[metabolite_mask] = center + distance * stretch_factor
    
    # Compress all non-metabolite regions
    compress_mask = ~stretch_mask
    if np.any(compress_mask):
        transformed[compress_mask] = transformed[compress_mask] * compress_factor
    
    # Additional compression for -10 to -15 region
    extra_compression_mask = (ppm_vals >= -15) & (ppm_vals <= -10)
    if np.any(extra_compression_mask):
        # Apply extra compression to this specific region
        transformed[extra_compression_mask] = transformed[extra_compression_mask] * 0.3
    
    # Ensure monotonicity (prevent crossing of values)
    transformed = np.sort(transformed)
    
    return transformed

# Get metabolite centers for transformation
metabolite_centers = list(radar_targets.values())
metabolite_centers = [c for c in metabolite_centers if -18 <= c <= 8]
# Remove duplicates and sort
metabolite_centers = sorted(list(set(metabolite_centers)))

# Transform ppm values
# Use uniform ppm (no transformation)
transformed_ppm = aligned_ppm.copy()
transformed_met_ppms = met_ppms.copy()

# Alternative: Use original ppm for spectra if transformation causes issues
# Uncomment the next line if you want to use original ppm for spectra
# transformed_ppm = aligned_ppm

# --- Plotting ---
fig = plt.figure(figsize=(14, 9))
gs = fig.add_gridspec(2, 1, height_ratios=[2, 1], hspace=0.1)
ax_spectra = fig.add_subplot(gs[0])
ax_table = fig.add_subplot(gs[1])  # Remove sharex to have independent control

# Set limits based on transformed values - positive to negative (left to right)
ppm_min, ppm_max = transformed_ppm.min(), transformed_ppm.max()
# Both plots: positive to negative (left to right)
ax_spectra.set_xlim(ppm_max, ppm_min)  # max (positive) to min (negative)
ax_table.set_xlim(ppm_max, ppm_min)    # max (positive) to min (negative)
# No invert_xaxis() needed - we want positive to negative

# Vertical lines - use the same transformed positions as rectangles
for met in metabolites:
    original_ppm = radar_targets[met]
    # Find the closest ppm value in aligned_ppm
    closest_idx = np.argmin(np.abs(aligned_ppm - original_ppm))
    ppm_val = transformed_ppm[closest_idx]
    ax_spectra.axvline(ppm_val, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax_table.axvline(ppm_val, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

# Test alignment lines removed as requested

# Plot spectra with proper ordering
n_spectra = len(processed)
vertical_spacing = 0.1
for i, (normed, _, color, label) in enumerate(processed):
    y_offset = i * vertical_spacing
    # Ensure proper ordering of data points
    sorted_indices = np.argsort(transformed_ppm)
    ax_spectra.plot(transformed_ppm[sorted_indices], normed[sorted_indices] + y_offset, 
                   color=color, label=label, linewidth=1.2)

yticks = [i * vertical_spacing for i in range(n_spectra)]
yticklabels = [p[3] for p in processed]
yticklabels.append("Summary")
ax_spectra.set_yticks(yticks)
ax_spectra.set_yticklabels(yticklabels[:-1])
ax_spectra.set_ylabel("Normalized Intensity + Offset")
ax_spectra.set_title("Stacked Spectra with Corresponding Ratio Table Below")
ax_spectra.grid(True, linestyle='--', alpha=0.3)
ax_spectra.legend(loc='upper left', fontsize=9)
plt.setp(ax_spectra.get_xticklabels(), visible=False)

# --- Ratio table plot ---
ax_table.set_ylim(-0.5, (n_spectra + 1.5) * rect_spacing)

ax_table.set_yticks(range(n_spectra + 1))
ax_table.set_yticklabels(yticklabels)
ax_table.set_xlabel("Chemical Shift (ppm)")
ax_table.set_ylabel("Spectra")
ax_table.grid(axis='y', linestyle='--', alpha=0.3)
ax_table.xaxis.set_ticks_position('bottom')
ax_table.spines['top'].set_visible(False)

# Rectangles for spectra and 1x1 summary overlay
rect_spacing = 1.1
for row in range(n_spectra + 1):
    for col, met in enumerate(metabolites):
        # Find the exact position of this metabolite in transformed_ppm
        original_ppm = radar_targets[met]
        # Find the closest ppm value in aligned_ppm
        closest_idx = np.argmin(np.abs(aligned_ppm - original_ppm))
        x_center = transformed_ppm[closest_idx]
        y_center = row * rect_spacing

        if row < n_spectra:
            val = ratios[row][met]
            color = processed[row][2]
            quality = abs(val)
            width = max(abs(val), 0.1) * 1.25  # Scale by 1.5
            height = max(quality, 0.05) * 1.25  # Scale by 1.5
        else:
            width, height = 1.0, 1.0  # fixed size for summary cell
            color = 'none'
            continue  # skip base summary box (we'll plot overlaid values next)

        rect = Rectangle(
            (x_center - width / 2, y_center - height / 2),
            width, height,
            facecolor=color,
            edgecolor='white' if row == n_spectra else 'black',
            linewidth=0.8,
            alpha=0.8
        )
        ax_table.add_patch(rect)

        # Arrow for direction (skip for summary row)
        if row < n_spectra:
            arrow_dir = 1 if val >= 0 else -1
            arrow_dx = 0
            arrow_dy = 0.3 * arrow_dir * abs(val)
            ax_table.arrow(
                x_center, y_center,
                arrow_dx, arrow_dy,
                head_width=0.1,
                head_length=0.05,
                fc='black',
                ec='black',
                linewidth=0.8,
                alpha=0.9,
                length_includes_head=True
            )

# --- Overlay summary rectangles in 1x1 cells ---
summary_y = (n_spectra + 0.5) * rect_spacing
for col, met in enumerate(metabolites):
    # Find the exact position of this metabolite in transformed_ppm
    original_ppm = radar_targets[met]
    # Find the closest ppm value in aligned_ppm
    closest_idx = np.argmin(np.abs(aligned_ppm - original_ppm))
    x_center = transformed_ppm[closest_idx]
    y_center = summary_y
    for i, val in enumerate(summary_cells[met]):
        quality = abs(val)
        width = max(abs(val), 0.1) * 1.2  # Scale by 1.2
        height = max(quality, 0.05) * 1.2  # Scale by 1.2
        color = processed[i][2]
        rect = Rectangle(
            (x_center - width / 2, y_center - height / 2),
            width, height,
            facecolor=color,
            edgecolor='black',
            linewidth=0.5,
            alpha=0.4
        )
        ax_table.add_patch(rect)

# --- Add direction legend arrows ---
# Create a separate axes for legends
legend_ax = fig.add_axes([0.02, 0.15, 0.15, 0.25])
legend_ax.axis('off')

# Ratio and Quality legend
legend_ax.text(0.0, 2.8, "Ratio →", fontsize=12, ha='left', va='center', fontweight='bold')
legend_ax.text(0.0, 2.7, "Quality ↓", fontsize=12, ha='left', va='center', fontweight='bold')

# Add arrows for ratio and quality
legend_ax.annotate("", xy=(0.3, 2.8), xytext=(0.1, 0.8),
                   arrowprops=dict(arrowstyle='->', lw=2.0, color='black'))
legend_ax.annotate("", xy=(0.2, 2.9), xytext=(0.2, 0.7),
                   arrowprops=dict(arrowstyle='->', lw=2.0, color='black'))

# Metabolite labels
for i, met in enumerate(metabolites):
    # Find the exact position of this metabolite in transformed_ppm
    original_ppm = radar_targets[met]
    # Find the closest ppm value in aligned_ppm
    closest_idx = np.argmin(np.abs(aligned_ppm - original_ppm))
    x_pos = transformed_ppm[closest_idx]
    ax_table.text(x_pos, -1.4, met, rotation=90, ha='center', va='top', fontsize=8)

ax_table.invert_yaxis()
plt.tight_layout()
plt.show()
