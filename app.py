from flask import Flask, render_template, request, jsonify, send_file
import os
import io
import base64
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import pywt
from scipy.interpolate import interp1d
import matplotlib.cm as cm

app = Flask(__name__)

# Global variables for data
spectra = None
chemical_shifts = None
processed_data = None

def load_data():
    """Load and preprocess MRS data"""
    global spectra, chemical_shifts, processed_data
    
    spectra = pd.read_csv("./data/spectra.csv", header=None)
    chemical_shifts = pd.read_csv("./data/chemical_shifts.csv", header=None).values.flatten()
    
    # Default settings
    selected_indices = [0, 1, 2, 3]
    default_colors = ['#440154', '#31688e', '#35b779', '#fde725']  # Viridis palette
    colors = default_colors
    
    # Align spectra by PCr peak
    sort_idx = np.argsort(chemical_shifts)[::-1]
    ppm = chemical_shifts[sort_idx]
    spectra_sorted = spectra.iloc[sort_idx]
    
    ref_signal = spectra_sorted.iloc[:, selected_indices[0]].values
    pc_range = (ppm > -1) & (ppm < 1)
    pcr_index = np.argmax(ref_signal[pc_range])
    pcr_ppm = ppm[pc_range][pcr_index]
    aligned_ppm = ppm - pcr_ppm
    ppm_mask = (aligned_ppm >= -20) & (aligned_ppm <= 20)
    aligned_ppm = aligned_ppm[ppm_mask]
    
    # Denoising function
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
    
    # Process selected spectra
    processed = []
    for i, idx in enumerate(selected_indices):
        raw = spectra_sorted.iloc[:, idx].values
        denoised = denoise_signal(raw)
        denoised_sorted = denoised[sort_idx]
        denoised_filtered = denoised_sorted[ppm_mask]
        normed = (denoised_filtered - np.min(denoised_filtered)) / (np.max(denoised_filtered) - np.min(denoised_filtered))
        processed.append((normed, denoised_filtered, colors[i], f"Spectrum {idx+1}"))
    
    processed_data = {
        'aligned_ppm': aligned_ppm,
        'processed': processed,
        'spectra_sorted': spectra_sorted,
        'ppm_mask': ppm_mask
    }

def create_circular_plot(selected_metabolites=None, ppm_range=None, saturation_range=0.15, color_palette=None):
    """Create circular plot with selected metabolites and ppm range"""
    if processed_data is None:
        load_data()
    
    aligned_ppm = processed_data['aligned_ppm']
    processed = processed_data['processed']
    
    # Use custom color palette if provided
    if color_palette and len(color_palette) >= 4:
        colors = color_palette
    else:
        colors = ['#440154', '#31688e', '#35b779', '#fde725']  # Default Viridis
    
    # Filter by ppm range if specified
    if ppm_range:
        ppm_min, ppm_max = ppm_range
        ppm_mask = (aligned_ppm >= ppm_min) & (aligned_ppm <= ppm_max)
        plot_ppm = aligned_ppm[ppm_mask]
        plot_processed = []
        for i, (norm, denoised, old_color, label) in enumerate(processed):
            new_color = colors[i % 4]  # Ensure we only use 4 colors
            plot_processed.append((norm[ppm_mask], denoised[ppm_mask], new_color, label))
    else:
        plot_ppm = aligned_ppm
        plot_processed = []
        for i, (norm, denoised, old_color, label) in enumerate(processed):
            new_color = colors[i % 4]  # Ensure we only use 4 colors
            plot_processed.append((norm, denoised, new_color, label))
    
    # Metabolite targets
    radar_targets = {
        "ATPβ": -16.15, "NAD+": -8.31, "NADH": -8.13, "ATPα": -7.56,
        "ATPγ": -2.53, "PCr": 0, "PME": 3.5, "Pi": 4.82, "PDE": 5.5
    }
    
    # Filter metabolites if specified
    if selected_metabolites:
        radar_targets = {k: v for k, v in radar_targets.items() if k in selected_metabolites}
    
    # Create plot
    fig, ax = plt.subplots(figsize=(18, 18))  # 1.5x larger
    ax.set_aspect('equal')
    ax.axis('off')
    
    R = 9  # 1.5x larger
    ax.set_xlim(-2.5*R, 2.5*R)
    ax.set_ylim(-2.5*R, 2.5*R)
    ax.add_patch(Circle((0, 0), R, edgecolor='black', facecolor='none', lw=2))
    
    # Helper functions
    def ppm_to_angle(ppm_val):
        return np.pi / 2 - 2 * np.pi * ppm_val / 40
    
    def radial_stretch(ppm_val, signal_val, base_radius=R, scale_factor=4, compress_factor=0.4):
        if -20 <= ppm_val <= 20:
            return base_radius + signal_val * scale_factor
        else:
            return base_radius + signal_val * scale_factor * compress_factor
    
    # Draw circular bands
    n_spectra = len(plot_processed)
    extra_orbit_radius = R - 3  # 1.5x larger gap
    spectrum_orbit_radii = []
    for i in range(1, n_spectra + 1):
        radius = (i / (n_spectra + 1)) * extra_orbit_radius
        spectrum_orbit_radii.append(radius)
    
    for i in range(1, n_spectra + 1):
        radius = spectrum_orbit_radii[i-1]
        ax.add_patch(Circle((0, 0), radius=radius, edgecolor='gray', facecolor='none',
                            linestyle='--', linewidth=1.0, alpha=0.6))
    
    ax.add_patch(Circle((0, 0), radius=extra_orbit_radius, edgecolor='blue', facecolor='none',
                        linestyle='--', linewidth=1.5, alpha=0.7))
    
    # Plot spectra with saturation around selected metabolites
    margin = 0.015  # 1.5x larger
    step = 0.75  # 1.5x larger
    angles = ppm_to_angle(plot_ppm)
    
    for i, (norm, denoised, color, label) in enumerate(plot_processed):
        base_r = R + margin + i * step
        radii = [radial_stretch(p, s, base_radius=base_r) for p, s in zip(plot_ppm, norm)]
        x = [r * np.cos(a) for r, a in zip(radii, angles)]
        y = [r * np.sin(a) for r, a in zip(radii, angles)]
        
        # Plot base spectrum
        ax.plot(x, y, color=color, lw=1.5, label=label)
        
        # Add saturation around selected metabolites
        for name, ppm_val in radar_targets.items():
            if -20 <= ppm_val <= 20:
                # Define saturation range around metabolite
                sat_mask = (plot_ppm >= ppm_val - saturation_range) & (plot_ppm <= ppm_val + saturation_range)
                
                if np.any(sat_mask):
                    sat_x = [x[j] for j in range(len(x)) if sat_mask[j]]
                    sat_y = [y[j] for j in range(len(y)) if sat_mask[j]]
                    
                    # Plot saturated region with thicker line and higher alpha
                    ax.plot(sat_x, sat_y, color=color, lw=3.0, alpha=0.8, zorder=10)
    
    # Annotate metabolites
    for name, ppm_val in radar_targets.items():
        if -20 <= ppm_val <= 20:
            angle = ppm_to_angle(ppm_val)
            r = radial_stretch(ppm_val, 1.2, R, 1.5)
            x0, y0 = r * np.cos(angle), r * np.sin(angle)
            x1, y1 = (r - 3.75) * np.cos(angle), (r - 3.75) * np.sin(angle)  # 1.5x larger
            ax.plot([x0, x1], [y0, y1], color='red', lw=1.5)  # 1.5x thicker
            ax.text(x1, y1, name, ha='center', va='center', fontsize=9)  # 1.5x larger font
    
    # Draw ppm ticks
    ppm_ticks = np.arange(-20, 21, 5)
    for ppm_val in ppm_ticks:
        angle = ppm_to_angle(ppm_val)
        r_outer = R + 0.6  # 1.5x larger
        r_inner = R - 0.15  # 1.5x larger
        x_outer, y_outer = r_outer * np.cos(angle), r_outer * np.sin(angle)
        x_inner, y_inner = r_inner * np.cos(angle), r_inner * np.sin(angle)
        ax.plot([x_inner, x_outer], [y_inner, y_outer], color='red', lw=1.2)  # 1.5x thicker
        ax.text(R * np.cos(angle) + 0.3, R * np.sin(angle) + 0.3, f'{ppm_val}', 
                fontsize=15, ha='center', va='center')  # 1.5x larger font
    
    # Compute quality metrics for each metabolite
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

    # Evaluate quality metrics across spectra
    quality_metrics = []
    for i, (normed, denoised, color, label) in enumerate(plot_processed):
        qm = {}
        for name, ppm_val in radar_targets.items():
            if -20 <= ppm_val <= 20:
                crlb_c, snr_c, fwhm_c = compute_quality_metrics(plot_ppm, denoised, denoised, ppm_val)
                qm[name] = [crlb_c == 'green', snr_c == 'green', fwhm_c == 'green']
        quality_metrics.append(qm)

    # Add partial circles on spectrum orbits based on quality metrics passed
    factor = 1.0
    for i, (norm, denoised, color, label) in enumerate(plot_processed):
        interp_fn = interp1d(plot_ppm, denoised, kind='linear', fill_value="extrapolate")
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
                    start_angle = 0
                    end_angle = 90
                elif passes == 1:
                    # One quality metric passed - draw 1/2 circle
                    start_angle = 0
                    end_angle = 180
                elif passes == 2:
                    # Two quality metrics passed - draw 3/4 circle
                    start_angle = 0
                    end_angle = 270
                else:
                    # All quality metrics passed - draw full circle
                    start_angle = 0
                    end_angle = 360
                
                # Create partial circle using Wedge
                from matplotlib.patches import Wedge
                wedge = Wedge((cx, cy), r=ratio * factor, theta1=start_angle, theta2=end_angle,
                             edgecolor=color, facecolor=color, lw=1.2)
                ax.add_patch(wedge)
    
    # Add quadrant-style glyphs on extra orbit
    for target_name, target_ppm in radar_targets.items():
        if -20 <= target_ppm <= 5:
            # Calculate ratios for this target metabolite
            target_ratios = []
            for norm, denoised, color, label in plot_processed:
                interp_fn = interp1d(plot_ppm, denoised, kind='linear', fill_value="extrapolate")
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
            base_x = extra_orbit_radius * np.cos(target_angle)
            base_y = extra_orbit_radius * np.sin(target_angle)
            
            # Create quadrant circle with 4 sections
            n_spectra = len(plot_processed)
            section_angle = 360 / n_spectra  # 90 degrees for 4 spectra
            
            for i, (ratio, color) in enumerate(zip(normalized_ratios, [p[2] for p in plot_processed])):
                start_angle = i * section_angle
                end_angle = (i + 1) * section_angle
                
                # Create wedge for this spectrum's section
                wedge = Wedge((base_x, base_y), r=ratio, theta1=start_angle, theta2=end_angle,
                             facecolor=color, edgecolor='black', lw=0.5, alpha=0.85)
                ax.add_patch(wedge)
            
            # Add metabolite label in the center (removed)
            # ax.text(base_x, base_y, target_name, fontsize=6, ha='center', va='center', 
            #        fontweight='bold', color='black')
    
    # === Right side vertical legend ===
    # Create separate legend areas on the right side, stacked vertically
    
    # Spectra legend (top)
    spectra_legend_ax = fig.add_axes([0.85, 0.55, 0.12, 0.15])  # Moved much lower
    spectra_legend_ax.axis('off')
    spectra_legend_ax.text(0.0, 0.75, "Spectra", fontsize=12, fontweight='bold')  # Moved down
    for i, col in enumerate([p[2] for p in plot_processed], start=1):
        y = 0.75 - i * 0.12  # Much reduced spacing
        spectra_legend_ax.add_patch(plt.Circle((0.15, y), 0.05, color=col))
        spectra_legend_ax.text(0.5, y - 0.015, f"Spectrum {i}", fontsize=10, verticalalignment='center')

    # Quality legend (middle)
    quality_ax = fig.add_axes([0.85, 0.42, 0.12, 0.12])  # Moved much lower and closer
    quality_ax.axis('off')
    quality_ax.text(0, 0.9, "Quality Metrics", fontsize=12, fontweight='bold')

    # Draw example partial circles for quality legend
    quality_examples = [
        (0.25, "0 passed", 0, 90),
        (0.5, "1 passed", 0, 180), 
        (0.75, "2 passed", 0, 270),
        (1.0, "3 passed", 0, 360)
    ]

    for i, (fraction, label, start, end) in enumerate(quality_examples):
        y = 0.75 - i * 0.12  # Much reduced spacing
        # Draw example wedge
        wedge = Wedge((0.15, y), r=0.08, theta1=start, theta2=end,
                      edgecolor='gray', facecolor='gray', lw=1)
        quality_ax.add_patch(wedge)
        quality_ax.text(0.3, y, label, fontsize=10, verticalalignment='center')

    # Ratio guide (bottom)
    ratio_ax = fig.add_axes([0.85, 0.29, 0.12, 0.12])  # Moved much lower and closer
    ratio_ax.axis('off')
    ratio_ax.text(0, 0.75, "Ratio to PCr", fontsize=12, fontweight='bold')
    for i, r in enumerate([0.25, 0.5, 0.75]):
        x = 0.2 + i * 0.3
        ratio_ax.add_patch(plt.Circle((x, 0.5), 0.1 * r, color='gray', alpha=0.6))
    
    # Adjust layout to prevent overlap
    plt.subplots_adjust(right=0.82)
    
    # Convert to base64
    img = io.BytesIO()
    plt.savefig(img, format='png', dpi=150, bbox_inches='tight')
    img.seek(0)
    plt.close()
    
    return base64.b64encode(img.getvalue()).decode()

def create_rectangular_plot(selected_metabolites=None, ppm_range=None, saturation_range=0.15, color_palette=None):
    """Create rectangular plot with selected metabolites and ppm range based on rectangle.py"""
    if processed_data is None:
        load_data()
    
    aligned_ppm = processed_data['aligned_ppm']
    processed = processed_data['processed']
    
    # Use custom color palette if provided
    if color_palette and len(color_palette) >= 4:
        colors = color_palette
    else:
        colors = ['#440154', '#31688e', '#35b779', '#fde725']  # Default Viridis
    
    # Filter by ppm range if specified
    if ppm_range:
        ppm_min, ppm_max = ppm_range
        ppm_mask = (aligned_ppm >= ppm_min) & (aligned_ppm <= ppm_max)
        plot_ppm = aligned_ppm[ppm_mask]
        plot_processed = []
        for i, (norm, denoised, old_color, label) in enumerate(processed):
            new_color = colors[i % 4]  # Ensure we only use 4 colors
            plot_processed.append((norm[ppm_mask], denoised[ppm_mask], new_color, label))
    else:
        plot_ppm = aligned_ppm
        plot_processed = []
        for i, (norm, denoised, old_color, label) in enumerate(processed):
            new_color = colors[i % 4]  # Ensure we only use 4 colors
            plot_processed.append((norm, denoised, new_color, label))
    
    # Metabolite targets
    radar_targets = {
        "ATPβ": -16.15, "NAD+": -8.31, "NADH": -8.13, "ATPα": -7.56,
        "ATPγ": -2.53, "PCr": 0, "PME": 3.5, "Pi": 4.82, "PDE": 5.5
    }
    
    # Filter metabolites if specified
    if selected_metabolites:
        radar_targets = {k: v for k, v in radar_targets.items() if k in selected_metabolites}
    
    metabolites = list(radar_targets.keys())
    met_ppms = np.array([radar_targets[m] for m in metabolites])
    sorted_idx = np.argsort(met_ppms)[::-1]
    metabolites = [metabolites[i] for i in sorted_idx]
    met_ppms = met_ppms[sorted_idx]
    
    # Calculate ratios for each spectrum
    ratios = {}
    for i, (norm, denoised, color, label) in enumerate(plot_processed):
        spectrum_ratios = {}
        
        # Integrate over ±2 indices around PCr at 0.0 ppm
        pcr_idx = np.argmin(np.abs(plot_ppm - 0.0))
        pcr_start = max(pcr_idx - 2, 0)
        pcr_end = min(pcr_idx + 3, len(plot_ppm))
        pcr_auc = np.trapz(denoised[pcr_start:pcr_end], x=plot_ppm[pcr_start:pcr_end])
        
        for met in metabolites:
            met_ppm = radar_targets[met]
            met_idx = np.argmin(np.abs(plot_ppm - met_ppm))
            start = max(met_idx - 2, 0)
            end = min(met_idx + 3, len(plot_ppm))
            met_auc = np.trapz(denoised[start:end], x=plot_ppm[start:end])
            
            ratio = (met_auc / pcr_auc if pcr_auc != 0 else 0) * 3  # scaled for visibility
            spectrum_ratios[met] = ratio
        
        ratios[i] = spectrum_ratios
    
    # Create summary cells for overlay
    summary_cells = {}
    for met in metabolites:
        summary_cells[met] = [ratios[i][met] for i in range(len(plot_processed))]
    
    # Create plot
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 1, height_ratios=[2, 1], hspace=0.1)
    ax_spectra = fig.add_subplot(gs[0])
    ax_table = fig.add_subplot(gs[1])
    
    # Set limits - positive to negative (left to right)
    ppm_min, ppm_max = plot_ppm.min(), plot_ppm.max()
    ax_spectra.set_xlim(ppm_max, ppm_min)
    ax_table.set_xlim(ppm_max, ppm_min)
    
    # Vertical lines for metabolites
    for met in metabolites:
        original_ppm = radar_targets[met]
        closest_idx = np.argmin(np.abs(plot_ppm - original_ppm))
        ppm_val = plot_ppm[closest_idx]
        ax_spectra.axvline(ppm_val, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax_table.axvline(ppm_val, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # Plot stacked spectra with saturation around selected metabolites
    n_spectra = len(plot_processed)
    vertical_spacing = 0.1
    for i, (normed, _, color, label) in enumerate(plot_processed):
        y_offset = i * vertical_spacing
        sorted_indices = np.argsort(plot_ppm)
        
        # Plot base spectrum
        ax_spectra.plot(plot_ppm[sorted_indices], normed[sorted_indices] + y_offset, 
                       color=color, label=label, linewidth=1.2)
        
        # Add saturation around selected metabolites
        for name, ppm_val in radar_targets.items():
            if -20 <= ppm_val <= 20:
                # Define saturation range around metabolite
                sat_mask = (plot_ppm >= ppm_val - saturation_range) & (plot_ppm <= ppm_val + saturation_range)
                
                if np.any(sat_mask):
                    sat_ppm = plot_ppm[sat_mask]
                    sat_normed = normed[sat_mask]
                    
                    # Plot saturated region with thicker line and higher alpha
                    ax_spectra.plot(sat_ppm, sat_normed + y_offset, 
                                   color=color, linewidth=2.5, alpha=0.8, zorder=10)
    
    yticks = [i * vertical_spacing for i in range(n_spectra)]
    yticklabels = [p[3] for p in plot_processed]
    yticklabels.append("Summary")
    ax_spectra.set_yticks(yticks)
    ax_spectra.set_yticklabels(yticklabels[:-1])
    ax_spectra.set_ylabel("Normalized Intensity + Offset")
    ax_spectra.set_title("Stacked Spectra with Corresponding Ratio Table Below")
    ax_spectra.grid(True, linestyle='--', alpha=0.3)
    ax_spectra.legend(loc='upper left', fontsize=9)
    plt.setp(ax_spectra.get_xticklabels(), visible=False)
    
    # Ratio table plot
    rect_spacing = 1.1
    ax_table.set_ylim(-0.5, (n_spectra + 1.5) * rect_spacing)
    ax_table.set_yticks(range(n_spectra + 1))
    ax_table.set_yticklabels(yticklabels)
    ax_table.set_xlabel("Chemical Shift (ppm)")
    ax_table.set_ylabel("Spectra")
    ax_table.grid(axis='y', linestyle='--', alpha=0.3)
    ax_table.xaxis.set_ticks_position('bottom')
    ax_table.spines['top'].set_visible(False)
    
    # Rectangles for spectra
    for row in range(n_spectra + 1):
        for col, met in enumerate(metabolites):
            original_ppm = radar_targets[met]
            closest_idx = np.argmin(np.abs(plot_ppm - original_ppm))
            x_center = plot_ppm[closest_idx]
            y_center = row * rect_spacing
            
            if row < n_spectra:
                val = ratios[row][met]
                color = plot_processed[row][2]
                quality = abs(val)
                width = max(abs(val), 0.1) * 1.25
                height = max(quality, 0.05) * 1.25
            else:
                width, height = 1.0, 1.0
                color = 'none'
                continue
            
            rect = Rectangle(
                (x_center - width / 2, y_center - height / 2),
                width, height,
                facecolor=color,
                edgecolor='white' if row == n_spectra else 'black',
                linewidth=0.8,
                alpha=0.8
            )
            ax_table.add_patch(rect)
            
            # Arrow for direction
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
    
    # Overlay summary rectangles
    summary_y = (n_spectra + 0.5) * rect_spacing
    for col, met in enumerate(metabolites):
        original_ppm = radar_targets[met]
        closest_idx = np.argmin(np.abs(plot_ppm - original_ppm))
        x_center = plot_ppm[closest_idx]
        y_center = summary_y
        for i, val in enumerate(summary_cells[met]):
            quality = abs(val)
            width = max(abs(val), 0.1) * 1.2
            height = max(quality, 0.05) * 1.2
            color = plot_processed[i][2]
            rect = Rectangle(
                (x_center - width / 2, y_center - height / 2),
                width, height,
                facecolor=color,
                edgecolor='black',
                linewidth=0.5,
                alpha=0.4
            )
            ax_table.add_patch(rect)
    
    # Add direction legend
    legend_ax = fig.add_axes([0.02, 0.15, 0.15, 0.25])
    legend_ax.axis('off')
    legend_ax.text(0.0, 2.8, "Ratio →", fontsize=12, ha='left', va='center', fontweight='bold')
    legend_ax.text(0.0, 2.7, "Quality ↓", fontsize=12, ha='left', va='center', fontweight='bold')
    legend_ax.annotate("", xy=(0.3, 2.8), xytext=(0.1, 0.8),
                       arrowprops=dict(arrowstyle='->', lw=2.0, color='black'))
    legend_ax.annotate("", xy=(0.2, 2.9), xytext=(0.2, 0.7),
                       arrowprops=dict(arrowstyle='->', lw=2.0, color='black'))
    
    # Metabolite labels
    for i, met in enumerate(metabolites):
        original_ppm = radar_targets[met]
        closest_idx = np.argmin(np.abs(plot_ppm - original_ppm))
        x_pos = plot_ppm[closest_idx]
        ax_table.text(x_pos, -1.4, met, rotation=90, ha='center', va='top', fontsize=8)
    
    ax_table.invert_yaxis()
    plt.tight_layout()
    
    # Convert to base64
    img = io.BytesIO()
    plt.savefig(img, format='png', dpi=150, bbox_inches='tight')
    img.seek(0)
    plt.close()
    
    return base64.b64encode(img.getvalue()).decode()

@app.route('/')
def index():
    """Main page with controls"""
    return render_template('index.html')

@app.route('/generate_plot', methods=['POST'])
def generate_plot():
    """Generate plot based on user selections"""
    data = request.get_json()
    plot_type = data.get('plot_type', 'circular')
    selected_metabolites = data.get('metabolites', [])
    ppm_range = data.get('ppm_range', None)
    saturation_range = data.get('saturation_range', 0.15)
    color_palette = data.get('color_palette', None)
    
    if ppm_range:
        ppm_range = (float(ppm_range[0]), float(ppm_range[1]))
    
    if plot_type == 'circular':
        plot_data = create_circular_plot(selected_metabolites, ppm_range, saturation_range, color_palette)
    else:
        plot_data = create_rectangular_plot(selected_metabolites, ppm_range, saturation_range, color_palette)
    
    return jsonify({'plot_data': plot_data})

@app.route('/get_metabolites')
def get_metabolites():
    """Get available metabolites"""
    metabolites = {
        "ATPβ": -16.15, "NAD+": -8.31, "NADH": -8.13, "ATPα": -7.56,
        "ATPγ": -2.53, "PCr": 0, "PME": 3.5, "Pi": 4.82, "PDE": 5.5
    }
    return jsonify(metabolites)

if __name__ == '__main__':
    # Load data on startup
    load_data()
    app.run(debug=True, host='0.0.0.0', port=5000) 