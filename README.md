# Ionospheric Scintillation Simulator

## Project Description

This software implements a **phase screen model** for simulating the propagation of radio waves through a turbulent ionosphere along the **z-direction**. The transverse spatial dimension is represented by a **1D grid along the x-axis**, and the turbulence is modeled using **one-dimensional Power Spectral Density (PSD)** functions of spatial frequency \( **k<sub>x</sub>** \). The model uses split-step Fourier method for forward propagation between phase screens. It calculates the scintillation index (**S4**), integrating electron density profiles from the International Reference Ionosphere (IRI) model via PyIRI. The simulation supports various Power Spectral Density (PSD) models: Gaussian, von Karman, and Shkarofsky. Key features include adaptive phase screen placement weighted by electron density, validation conditions for simulation accuracy (e.g., grid resolution and thin-screen approximations), and a graphical user interface (GUI) for parameter input and result visualization.



## Dependencies

- Python 3.12.3

- Required packages (see requirements.txt):

  - numpy==1.26.4

  - SciPy==1.16.2

  - matplotlib==3.10.7

  - PyIRI==0.1.3

  - pyproj==3.7.2

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

- Frequency: 1000.0 MHz

- S4: ~0.125 (weak scintillation)

- dφ_max: 0.217 rad (average for all implementations)

- n_runs: 1000 (number of implementations)

- number_of_screens: 16

- z_len: 202.3 km (path length)

- max dz: 20117.0 m (max dz from the entire path for all implementations)
  
- Warnings: None if conditions met

Screenshots of the interface:

- Initial GUI (screenshots/start.png):  

  Shows parameter inputs and "RUN SIMULATION" button

- During Run (screenshots/running.png):

  Progress bar and status messages

- Results (screenshots/results.png):

  Displays S4, errors, warnings, etc. in the text box

For reproducibility, set `np.random.seed(42)` in `phase_screen_generator.py`.
