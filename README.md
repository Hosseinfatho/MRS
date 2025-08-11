# MRS Visualization Tool

An interactive web application for visualizing Magnetic Resonance Spectroscopy (MRS) data with both circular and rectangular plot representations.

## Overview

This tool provides a comprehensive interface for analyzing and visualizing MRS data, featuring:

- **Interactive metabolite selection** with 9 common metabolites
- **Dual visualization modes**: Circular and rectangular plots
- **Advanced customization options**: PPM ranges, color palettes, saturation controls
- **Export functionality**: Download plots in PNG, SVG, and PDF formats
- **Responsive design**: Works on desktop and mobile devices

## Features

### Metabolite Analysis
- **ATPβ** (-16.15 PPM)
- **NAD+** (-8.31 PPM)
- **NADH** (-8.13 PPM)
- **ATPα** (-7.56 PPM)
- **ATPγ** (-2.53 PPM)
- **PCr** (0 PPM)
- **PME** (3.5 PPM)
- **Pi** (4.82 PPM)
- **PDE** (5.5 PPM)

### Visualization Options
- **Circular Plot**: Polar representation of MRS data
- **Rectangular Plot**: Traditional linear MRS spectrum
- **PPM Range Selection**: Full (-20 to 20), Central (-5 to 5), Extended (-10 to 10), Limited (-2 to 2)
- **Color Palettes**: Viridis, Plasma, Inferno, Magma, Cividis
- **Saturation Control**: Adjustable intensity scaling (0.05 to 0.30)

### Export Features
- High-resolution PNG export
- Vector-based SVG export
- Print-ready PDF export
- Customizable dimensions

## Live Demo

Visit the live application: [MRS Visualization Tool](https://[username].github.io/[repository-name]/)

## Installation

### Option 1: GitHub Pages (Recommended)
1. Fork this repository
2. Enable GitHub Pages in repository settings
3. Set source to "Deploy from a branch"
4. Select main branch and root folder
5. Access your site at `https://[username].github.io/[repository-name]/`

### Option 2: Local Development
1. Clone the repository
2. Open `index.html` in a web browser
3. No additional setup required

## Usage

1. **Select Metabolites**: Choose which metabolites to visualize using the checkboxes
2. **Choose Plot Type**: Select between circular or rectangular visualization
3. **Adjust Parameters**: Modify PPM range, saturation, and color palette
4. **Generate Plot**: Click "Generate Plot" to create the visualization
5. **Export**: Download the plot in your preferred format

## Technical Details

### Technologies Used
- **HTML5**: Structure and semantics
- **CSS3**: Styling and responsive design
- **JavaScript**: Interactive functionality
- **Plotly.js**: Advanced plotting library
- **PapaParse**: CSV data processing

### Browser Compatibility
- Chrome 60+
- Firefox 55+
- Safari 12+
- Edge 79+

### Performance
- Optimized for real-time interaction
- Efficient data processing
- Responsive design for all screen sizes

## Development

### Project Structure
```
MRS/
├── index.html              # Main application file
├── README.md              # Project documentation
├── README_gh_pages.md     # GitHub Pages setup guide
├── 404.html              # Custom error page
├── CNAME                 # Custom domain configuration
├── sitemap.xml           # SEO sitemap
├── .gitignore            # Git ignore rules
└── data/                 # Sample data files
    ├── chemical_shifts.csv
    ├── single_spectrum.csv
    └── spectra.csv
```

### Customization
- Modify metabolite definitions in the JavaScript section
- Adjust color palettes and styling in CSS
- Add new plot types by extending the plotting functions
- Integrate with real data sources

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Plotly.js** for the excellent plotting library
- **GitHub Pages** for hosting
- **MRS Community** for scientific guidance

## Support

For questions, issues, or feature requests:
- Create an issue in the GitHub repository
- Contact the development team
- Check the documentation

## Version History

- **v1.0.0**: Initial release with basic functionality
- **v1.1.0**: Added metabolite selection and advanced controls
- **v1.2.0**: Enhanced visualization options and export features
- **v1.3.0**: Improved responsive design and performance

---

**Note**: This tool is designed for educational and research purposes. For clinical applications, please ensure compliance with relevant medical device regulations.
