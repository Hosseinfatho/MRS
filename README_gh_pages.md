# MRS Visualization Tool - GitHub Pages

This project is an interactive visualization tool for Magnetic Resonance Spectroscopy (MRS) data that runs as a static application on GitHub Pages.

## Features

- **Circular Plot**: Display MRS data in polar format
- **Rectangular Plot**: Traditional MRS data visualization
- **Interactive Controls**: Change PPM range, color palette, and other parameters
- **Download Plots**: Export plots in PNG, SVG, and PDF formats
- **English Interface**: Designed for English-speaking users

## How to Use

### View Online
The application is available online at:
```
https://[username].github.io/[repository-name]/
```

### Run Locally
1. Download the `index.html` file
2. Open it in your browser
3. Select your desired settings
4. Click "Generate Plot"

## GitHub Pages Setup

### Step 1: Create New Repository
1. Create a new repository on GitHub
2. Name it `[username].github.io` (for personal pages)
3. Or choose any name (for project pages)

### Step 2: Upload Files
1. Upload the `index.html` file to the repository
2. Upload the `README.md` file

### Step 3: Enable GitHub Pages
1. Go to Settings > Pages
2. Set Source to "Deploy from a branch"
3. Set Branch to "main" and folder to "/ (root)"
4. Click Save

### Step 4: Wait for Deployment
- GitHub Pages typically activates within a few minutes
- The site URL will be displayed in the Pages section

## File Structure

```
repository/
├── index.html          # Main application file
├── README.md           # Project description
└── README_gh_pages.md  # This file
```

## Technologies Used

- **HTML5**: Main structure
- **CSS3**: Styling and design
- **JavaScript**: Interactive logic
- **Plotly.js**: Chart library
- **PapaParse**: CSV data processing

## Technical Features

### Circular Plot
- Uses `scatterpolar` in Plotly
- Displays data in polar format
- Adjustable angular range

### Rectangular Plot
- Uses `scatter` in Plotly
- Traditional MRS data display
- Configurable axes

### Color Palettes
- Viridis (default)
- Plasma
- Inferno
- Magma
- Cividis

## Troubleshooting

### Common Issues

1. **Plot Not Displaying**
   - Ensure internet connection is active
   - Plotly and PapaParse CDNs must be accessible

2. **Font Display Issues**
   - Browser must support system fonts
   - Check system font settings

3. **Download Not Working**
   - Plotly.js must be properly loaded
   - Browser must support file downloads

## Contributing

To contribute to improving this project:

1. Fork the repository
2. Create a new branch
3. Make your changes
4. Submit a Pull Request

## License

This project is released under the MIT License.

## Contact

For questions and suggestions, please create an Issue or contact via email.

---

**Note**: This static version is designed for demonstration purposes. For use with real data, use the original Flask version.
