# Ionospheric Scintillation Simulator

## Project Description

This software simulates the propagation of radio waves through ionosphere using the phase screen method. It calculates the scintillation index (S4), integrating electron density profiles from the International Reference Ionosphere (IRI) model via PyIRI. The simulation supports various Power Spectral Density (PSD) models: Gaussian, von Karman, and Shkarofsky. Key features include adaptive phase screen placement weighted by electron density, validation conditions for simulation accuracy (e.g., grid resolution and thin-screen approximations), and a graphical user interface (GUI) for parameter input and result visualization.



## Dependencies

- Python 3.8+

- Required packages (see requirements.txt):

  - numpy>=1.26.0

  - scipy>=1.14.0

  - matplotlib>=3.9.0

  - pyproj>=3.6.1

  - PyIRI>=1.0.1

## Installation

1. Clone the repository: `git clone https://github.com/Margo-10/phase_screen_calculations.git`

2. Create and activate a virtual environment:  

   - `python -m venv venv`  

   - Linux/Mac: `source venv/bin/activate`  

   - Windows: `venv\Scripts\activate`

3. Install dependencies: `pip install -r requirements.txt`

4. Run the GUI: `python test.py`

## Usage

- Launch the application: `python test.py`

- Input parameters in the GUI:

  - PSD Type: Gaussian (0), von Karman (1), Shkarofsky (2)

  - Turbulence scales: Outer/Inner Scale [m], Spectral Index

  - Simulation: Number of Screens, dN/N (fluctuation level), Frequency [GHz]

  - Path: Source/Receiver coordinates (lat, lon, alt [km])

  - IRI: Date (Year, Month, Day, Hour), F10.7 solar flux

- Click "RUN SIMULATION" to start. The tool performs 1000 runs (by default), adapts screen count if phase deviation >0.25 rad, and displays S4, errors, and warnings.

- Results include S4 mean/error, max phase, path length, and validation warnings.

## Test Calculation

To demonstrate functionality, run with default parameters:

- PSD Type: Shkarofsky (2)

- Outer Scale: 5000 m

- Inner Scale: 10 m

- Spectral Index: 1.5

- Num Screens: 4

- dN/N: 0.1

- Frequency: 1 GHz

- Source: Lat 75.0, Lon 0.0, Alt 400 km

- Receiver: Lat 75.0, Lon 1.0, Alt 200 km

- Date: 2020-06-21, Hour 14, F10.7: 150

Expected output (approximate, due to randomness):

- Path length (z_len): 202.3 km

- S4: ~0.123 (weak scintillation)

- Mean max dφ: 0.213 rad 

- Warnings: None if conditions met

Screenshots of the interface:

- Initial GUI (screenshots/start.png):  

  (Shows parameter inputs and "RUN SIMULATION" button)

- During Run (screenshots/running.png):

  (Progress bar and status messages)

- Results (screenshots/results.png):

  (Displays S4, errors, warnings, etc. in the text box)

For reproducibility, set `np.random.seed(42)` in `phase_screen_generator.py`.
