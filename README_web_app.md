# MRS Visualization Web Application

## Overview
This web application provides an interactive interface for visualizing Magnetic Resonance Spectroscopy (MRS) data using two different visualization approaches: Circular and Rectangular plots.

## Features

### Plot Type Selection
- **Circular Plot**: Displays MRS spectra in a circular format with metabolite annotations and quality metrics
- **Rectangular Plot**: Shows stacked spectra with corresponding ratio tables in a traditional rectangular layout
- **Single Plot Display**: Only one plot type is shown at a time for better focus and performance

### Interactive Controls
- **Plot Type Buttons**: Switch between Circular and Rectangular visualizations
- **PPM Range Selection**: Adjust the chemical shift range to focus on specific regions
- **Metabolite Selection**: Choose which metabolites to highlight in the visualization
- **Real-time Updates**: Plots automatically update when settings change

### Data Processing
- **Wavelet Denoising**: Automatic signal denoising using BayesShrink thresholding
- **Peak Alignment**: Spectra aligned by PCr peak (around 0 ppm)
- **Quality Metrics**: CRLB, SNR, and FWHM calculations for metabolite peaks
- **Ratio Calculations**: Metabolite-to-PCr ratios with visual indicators

## Installation

1. Install required dependencies:
```bash
pip install -r requirements.txt
```

2. Ensure data files are in the `data/` directory:
   - `spectra.csv`: MRS spectral data
   - `chemical_shifts.csv`: Chemical shift values

## Usage

1. Start the application:
```bash
python app.py
```

2. Open your web browser and navigate to `http://localhost:5000`

3. Use the interface:
   - Select plot type (Circular or Rectangular)
   - Adjust PPM range as needed
   - Select metabolites to highlight
   - Click "Generate Plot" to update the visualization

## Plot Types

### Circular Plot
- **Radial Layout**: Spectra displayed in circular orbits
- **Metabolite Annotations**: Red lines pointing to known metabolite positions
- **Quality Indicators**: Partial circles showing quality metrics passed
- **Quadrant Glyphs**: Summary ratios on outer orbit
- **PPM Scale**: Radial ticks showing chemical shift values

### Rectangular Plot
- **Stacked Spectra**: Multiple spectra displayed vertically
- **Ratio Table**: Rectangular cells showing metabolite-to-PCr ratios
- **Color Coding**: Each spectrum has a unique color
- **Direction Indicators**: Arrows showing ratio direction and magnitude
- **Summary Overlay**: Combined ratios in bottom row

## Technical Details

### Backend (Flask)
- **Data Loading**: Automatic loading and preprocessing of MRS data
- **Signal Processing**: Wavelet denoising and normalization
- **Plot Generation**: Dynamic plot creation based on user selections
- **API Endpoints**: RESTful API for plot generation and metabolite data

### Frontend (HTML/JavaScript)
- **Responsive Design**: Works on desktop and mobile devices
- **Interactive Controls**: Real-time plot updates
- **Modern UI**: Clean, professional interface with gradient backgrounds
- **Error Handling**: Graceful error display and recovery

## Data Format

### Input Files
- **spectra.csv**: Matrix of spectral data (rows: chemical shifts, columns: spectra)
- **chemical_shifts.csv**: Vector of chemical shift values in ppm

### Metabolites
The application tracks these key metabolites:
- ATPβ (-16.15 ppm)
- NAD+ (-8.31 ppm)
- NADH (-8.13 ppm)
- ATPα (-7.56 ppm)
- ATPγ (-2.53 ppm)
- PCr (0 ppm)
- PME (3.5 ppm)
- Pi (4.82 ppm)
- PDE (5.5 ppm)

## Performance Notes

- **Single Plot Display**: Improved performance by showing only one plot at a time
- **Caching**: Data is loaded once and cached for faster plot generation
- **Optimized Rendering**: High-quality plots with reasonable file sizes
- **Responsive Updates**: Automatic regeneration when parameters change

## Browser Compatibility

- Chrome/Chromium (recommended)
- Firefox
- Safari
- Edge

## Troubleshooting

1. **Plot not generating**: Check that data files are in the correct location
2. **Slow performance**: Try reducing the PPM range or selecting fewer metabolites
3. **Display issues**: Ensure your browser supports modern CSS and JavaScript features

## Future Enhancements

- Additional plot types (3D visualization, heatmaps)
- Export functionality (PNG, PDF, SVG)
- Advanced filtering and analysis tools
- User preference saving
- Batch processing capabilities 