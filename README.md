# RASPA3 Porosity Analyser for Mercury

*[Version française : LISEZMOI.md](LISEZMOI.md)*

A Python tool for **Mercury** (CCDC) that compares the **measured** and **simulated** porosity of a porous material (MOF, ZIF, COF…). It takes the crystal structure open in Mercury, runs gas-adsorption simulations with **RASPA3**, computes BET surface areas and geometric pore descriptors, and puts everything side by side with your experimental isotherm.

Version 1.4.

## What it does

1. Takes the structure open in Mercury, or any CIF file, and builds a clean P1 model. Disorder is resolved, and solvent molecules not bonded to the framework can be removed.
2. Writes the RASPA3 input files and runs one **GCMC simulation per pressure point** of the isotherm. Several points can run in parallel.
3. Computes the **BET surface area** of the measured and simulated isotherms. The pressure range is chosen automatically and checked against the four **Rouquerol criteria**.
4. Computes the **pore volume** with the Gurvich rule.
5. Computes the **geometric properties**: accessible surface area, pore volume, pore limiting diameter (PLD) and largest cavity diameter (LCD). It uses the CSD Pore Analyser when possible, otherwise a documented internal method.
6. Compares experiment, simulation and geometry in one table with figures, and exports **Excel, CSV, PDF, PNG and an HTML report** shown in Mercury.

**Adsorbates available:** N₂ (77 K), Ar (87 K), CO₂ (273 K).
**Experimental data:** `.xlsx`, `.xls`, `.csv` or a table pasted into the window. BELSORP/BELMaster exports are recognised.

The window has five tabs: *Structure*, *Experimental data*, *RASPA settings*, *Calculation*, *Results*. The interface, messages and reports are in English.

## Requirements

| What | Why | Where to get it |
|---|---|---|
| **Mercury with the CSD Python API** | reads the structure, Pore Analyser | [CCDC](https://www.ccdc.cam.ac.uk/) (requires a CSD licence) |
| **RASPA3** | adsorption simulations (GCMC) | [RASPA3 website](https://iraspa.github.io/RASPA3/index.html) |
| numpy, pandas, matplotlib, openpyxl | calculations, figures, Excel export | `pip` (see below) |
| **PyQt6** (or PyQt5, PySide6, PySide2) | the graphical window | `pip`, in a **separate folder** (see below) |
| scipy, xlrd *(optional)* | faster geometry, reading old `.xls` files | `pip` |

## Installation

### Step 1: install RASPA3

The simplest method on **macOS, Windows and Linux** is conda (Miniconda or Anaconda). conda-forge provides builds for all three systems.

```
conda create -n raspa -c conda-forge raspa3 -y
```

Then find the path of the `raspa3` program. You will need it later:

- macOS: `conda activate raspa` then `which raspa3`
- Windows (Anaconda Prompt): `conda activate raspa` then `where raspa3`

On macOS you can also use the installer (`.pkg`) from the [RASPA3 releases page](https://github.com/iRASPA/raspa3/releases). The program is then usually in `/usr/local/bin/raspa3`, where the tool finds it automatically.

### Step 2: add the Python modules to Mercury's Python

The tool runs with the Python of the CSD Python API. Typical locations:

- macOS: `~/CCDC/ccdc-software/csd-python-api/miniconda/bin/python`
- Windows: look for `python.exe` inside the `csd-python-api` folder of your CCDC installation.

In a terminal, replace `CSD_PYTHON` below with that path (in quotes if it contains spaces):

```
CSD_PYTHON -m pip install --target ~/RASPA_Mercury/qt_libs PyQt6
CSD_PYTHON -m pip install pandas matplotlib openpyxl xlrd scipy
```

On Windows, write the first command with `%USERPROFILE%\RASPA_Mercury\qt_libs` instead of `~/RASPA_Mercury/qt_libs`.

> ⚠️ **Do not install PyQt directly into Mercury's Python.** On macOS it conflicts with the Qt libraries of the CSD Python API, and Mercury scripts stop working. That is why PyQt goes into the separate folder `RASPA_Mercury/qt_libs`.

**Easier:** run the tool once (step 3). If a module is missing, it displays the **exact commands** to copy into the terminal for your computer.

### Step 3: add the script to Mercury

1. Put `raspa3_mercury_porosity.py` in a folder that you keep for your Mercury scripts.
2. In Mercury, open the **CSD Python API** menu → **Options…** → **Add Location**, and choose that folder.
3. The tool now appears in the **CSD Python API** menu (on some Mercury versions: **Calculate → Python Scripts**).

## Usage

1. Open the structure in Mercury and run **CSD Python API → raspa3_mercury_porosity**. The window opens.
2. **Structure** tab: check the P1 model. Disorder and solvent options are on by default.
3. **Experimental data** tab: load your isotherm file or paste a table (P/P0 and adsorbed amount in cm³ STP/g).
4. **RASPA settings** tab: choose the gas, check the cycle numbers, and give the path of `raspa3` the first time. It is remembered afterwards. The **Auto** button picks the smallest supercell allowed by the cut-off.
5. **Calculation** tab: click **Run RASPA Simulation**. A progress bar shows the remaining time. You can also click **Analyse without simulation** for a quick experiment-plus-geometry comparison.
6. **Results** tab: figures and the comparison table. All files are saved in `~/RASPA_Mercury/<structure>_<gas>_<date>/`: `_results.xlsx`, `_report.pdf`, `_report.html`, CSV files and PNG figures.

**Without Mercury (batch mode, for tests or remote machines):**

```
python raspa3_mercury_porosity.py structure.cif --batch --exp measurement.xlsx --gas N2 --raspa /path/to/raspa3
python raspa3_mercury_porosity.py --help      # all options
```

## Modelling choices (to interpret the results)

- Rigid framework; UFF or DREIDING Lennard-Jones parameters per element; Lorentz–Berthelot mixing rules; potential truncated and shifted at 12 Å.
- Adsorbates: TraPPE models for N₂ and CO₂ (3 charged sites); single-site Lennard-Jones model for Ar.
- Framework charges: RASPA3 charge equilibration, charges from the CIF, or no charges.
- The simulation gives **absolute** adsorption, whereas experiments measure **excess** adsorption. At 77 or 87 K and P ≤ P₀, the difference is negligible.

## Limitations

- Simulations can take several hours, depending on the structure size, the number of points and the number of cycles.
- Flexible frameworks and chemisorption are not modelled.
- Developed and used on macOS. Windows should work, but has not been fully tested yet. Feedback is welcome.

## Credits and references

- **RASPA3:** « RASPA3: A Monte Carlo code for computing adsorption and diffusion in nanoporous materials and thermodynamics properties of fluids », *J. Chem. Phys.* **2024**, 161, 114106. <https://iraspa.github.io/RASPA3/>
- **BET consistency criteria:** J. Rouquerol, P. Llewellyn, F. Rouquerol, *Stud. Surf. Sci. Catal.* **2007**, 160, 49.
- **Accessible surface area:** T. Düren *et al.*, *J. Phys. Chem. C* **2007**, 111, 15350. **Probe-occupiable volume:** D. Ongari *et al.*, *Langmuir* **2017**.
- **Force fields:** UFF (Rappé *et al.*, *J. Am. Chem. Soc.* **1992**, 114, 10024); DREIDING (Mayo *et al.*, *J. Phys. Chem.* **1990**, 94, 8897); TraPPE (Potoff & Siepmann, 2001).
- The tool uses the CSD Python API from the Cambridge Crystallographic Data Centre (CCDC).
- Written by Yoann Rousselin with the help of Claude (Anthropic) as a coding assistant.

## Version history

- **1.4** (November 2026): the whole tool is now in English (interface, messages, reports, Excel sheets and columns, file names, code comments). Output files are renamed accordingly, for example `_results.xlsx` and `_report.pdf` instead of `_resultats.xlsx` and `_rapport.pdf`.
- **1.3** (November 2026): first public version.

## Licence

MIT. See [LICENSE](LICENSE). RASPA3 is distributed under its own licence.

## Author

**Yoann Rousselin**, CNRS research engineer (X-ray crystallography), ICMUB, Dijon (France).
Questions or bug reports: please open an *Issue* on this page.
