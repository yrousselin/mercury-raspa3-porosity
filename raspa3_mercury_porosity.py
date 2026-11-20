#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
===============================================================================
 RASPA3 Porosity Analyser for Mercury (CCDC)
===============================================================================
 Automatic comparison of the EXPERIMENTAL and SIMULATED porous properties
 of a porous material (MOF, ZIF, COF...):

   1. reads the structure currently open in Mercury (CSD Python API);
   2. exports it to a temporary CIF, then to a "cleaned" P1 CIF;
   3. prepares the RASPA3 input files (simulation.json, force_field.json,
      <molecule>.json, framework.cif);
   4. runs RASPA3 (GCMC) for each pressure point of the isotherm;
   5. reads the results and computes the BET area (Rouquerol criteria);
   6. computes the geometric properties (accessible surface area, pore volume,
      PLD, LCD) with the CSD Python API Pore Analyser or, failing that, with a
      documented internal implementation;
   7. compares with the experimental data and draws the figures;
   8. exports Excel, CSV, PDF and an HTML report shown in Mercury.

 Available adsorbates: N2 (77 K), Ar (87 K), CO2 (273 K).

-------------------------------------------------------------------------------
 INSTALLATION IN MERCURY
-------------------------------------------------------------------------------
   * Copy this file into a scripts folder, then in Mercury:
       CSD Python API > Options... > Add Location  (add that folder)
     The script then appears in the "CSD Python API" menu
     ("Calculate > Python Scripts" menu depending on the Mercury version).
   * Modules required in the Python environment used by Mercury:
       numpy, pandas, matplotlib, openpyxl, and a Qt binding
       (PyQt6, PyQt5, PySide6 or PySide2 - the first one found is used).
       Optional: scipy (faster geometric calculation), xlrd (.xls files).
     Example (terminal, with the CSD Python API python):
       python -m pip install --target ~/RASPA_Mercury/qt_libs PyQt6
       python -m pip install matplotlib openpyxl xlrd scipy
     WARNING: do not install PyQt directly into Mercury's Python
     (conflict with the Qt libraries of the CSD Python API on macOS); the
     script prints the exact commands if something is missing.
   * RASPA3: install the macOS package (.pkg) or "conda install raspa3".
     The executable path is set in the "RASPA Settings" tab
     (it is remembered between sessions).

-------------------------------------------------------------------------------
 USE OUTSIDE MERCURY (testing, batch processing)
-------------------------------------------------------------------------------
   python raspa3_mercury_porosity.py structure.cif            # GUI
   python raspa3_mercury_porosity.py structure.cif --batch \
          --exp measurement.xls --gas N2 --raspa /usr/local/bin/raspa3
   In this mode the CIF is read by the internal parser (CIFParser class):
   this is the documented fallback when the CSD Python API is not available.

-------------------------------------------------------------------------------
 MODELLING CHOICES (worth knowing when interpreting the results)
-------------------------------------------------------------------------------
   * Rigid framework; UFF or DREIDING Lennard-Jones parameters per ELEMENT,
     Lorentz-Berthelot mixing rules, truncated-and-shifted potential at 12 Å.
   * Adsorbates: TraPPE N2 and CO2 (3 sites, charges), 1-site Lennard-Jones Ar.
   * Framework charges: RASPA3 charge equilibration (Wilmer-Snurr),
     charges from the CIF (_atom_site_charge), or none.
   * "Equilibration cycles" in the GUI = RASPA3 NumberOfInitializationCycles
     (the RASPA3 "Equilibration" phase is only used for CFCMC).
   * The simulation gives ABSOLUTE adsorption; the experiment measures excess
     adsorption. At 77/87 K and P <= P0 the difference is negligible.
   * Crystallographic disorder: partially occupied sites that clash sterically
     are resolved by keeping a single component; guest species not bonded to
     the framework (solvent) can be removed.
===============================================================================
"""

from __future__ import annotations

import argparse
import base64
import datetime as _dt
import glob
import importlib
import io as _io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import traceback
from collections import deque
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import pandas as pd

APP_NAME = "RASPA3 Porosity Analyser"
APP_VERSION = "1.4"

# =============================================================================
# 1. PHYSICAL CONSTANTS AND TABLES
# =============================================================================
N_AVOGADRO = 6.02214076e23          # mol-1
V_MOLAR_STP = 22414.0               # cm3(STP)/mol (0 °C, 1 atm)
KCAL_TO_K = 503.2195                # 1 kcal/mol expressed in kelvin (E/kB)

# Atomic masses (g/mol)
ATOMIC_MASS = {
    "H": 1.008, "D": 2.014, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.95, "Tc": 98.0, "Ru": 101.07, "Rh": 102.91,
    "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91,
    "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24,
    "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.50,
    "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97,
    "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23,
    "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59, "Tl": 204.38,
    "Pb": 207.2, "Bi": 208.98, "Th": 232.04, "U": 238.03,
}

ATOMIC_NUMBER = {s: i + 1 for i, s in enumerate(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I "
    "Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt "
    "Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U".split())}

# Covalent radii (Å, Cordero 2008): bond / disorder detection
COVALENT_RADIUS = {
    "H": 0.31, "D": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84,
    "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66,
    "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Ar": 1.06, "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53,
    "Cr": 1.39, "Mn": 1.50, "Fe": 1.42, "Co": 1.38, "Ni": 1.24, "Cu": 1.32,
    "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20,
    "Kr": 1.16, "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64,
    "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45,
    "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39,
    "Xe": 1.40, "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Pr": 2.03,
    "Nd": 2.01, "Sm": 1.98, "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92,
    "Ho": 1.92, "Er": 1.89, "Tm": 1.90, "Yb": 1.87, "Lu": 1.87, "Hf": 1.75,
    "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36,
    "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "Th": 2.06,
    "U": 1.96,
}

# UFF (Rappé et al., JACS 1992, 114, 10024): (x_i [Å], D_i [kcal/mol]).
# Lennard-Jones conversion: sigma = x_i / 2^(1/6); epsilon/kB = D_i * 503.22
_UFF_RAW = {
    "H": (2.886, 0.044), "D": (2.886, 0.044), "He": (2.362, 0.056),
    "Li": (2.451, 0.025), "Be": (2.745, 0.085), "B": (4.083, 0.180),
    "C": (3.851, 0.105), "N": (3.660, 0.069), "O": (3.500, 0.060),
    "F": (3.364, 0.050), "Ne": (3.243, 0.042), "Na": (2.983, 0.030),
    "Mg": (3.021, 0.111), "Al": (4.499, 0.505), "Si": (4.295, 0.402),
    "P": (4.147, 0.305), "S": (4.035, 0.274), "Cl": (3.947, 0.227),
    "Ar": (3.868, 0.185), "K": (3.812, 0.035), "Ca": (3.399, 0.238),
    "Sc": (3.295, 0.019), "Ti": (3.175, 0.017), "V": (3.144, 0.016),
    "Cr": (3.023, 0.015), "Mn": (2.961, 0.013), "Fe": (2.912, 0.013),
    "Co": (2.872, 0.014), "Ni": (2.834, 0.015), "Cu": (3.495, 0.005),
    "Zn": (2.763, 0.124), "Ga": (4.383, 0.415), "Ge": (4.280, 0.379),
    "As": (4.230, 0.309), "Se": (4.205, 0.291), "Br": (4.189, 0.251),
    "Kr": (4.141, 0.220), "Rb": (4.114, 0.040), "Sr": (3.641, 0.235),
    "Y": (3.345, 0.072), "Zr": (3.124, 0.069), "Nb": (3.165, 0.059),
    "Mo": (3.052, 0.056), "Tc": (2.998, 0.048), "Ru": (2.963, 0.056),
    "Rh": (2.929, 0.053), "Pd": (2.899, 0.048), "Ag": (3.148, 0.036),
    "Cd": (2.848, 0.228), "In": (4.463, 0.599), "Sn": (4.392, 0.567),
    "Sb": (4.420, 0.449), "Te": (4.470, 0.398), "I": (4.500, 0.339),
    "Xe": (4.404, 0.332), "Cs": (4.517, 0.045), "Ba": (3.703, 0.364),
    "La": (3.522, 0.017), "Ce": (3.556, 0.013), "Pr": (3.606, 0.010),
    "Nd": (3.575, 0.010), "Sm": (3.520, 0.008), "Eu": (3.493, 0.008),
    "Gd": (3.368, 0.009), "Tb": (3.451, 0.007), "Dy": (3.428, 0.007),
    "Ho": (3.409, 0.007), "Er": (3.391, 0.007), "Tm": (3.374, 0.006),
    "Yb": (3.355, 0.228), "Lu": (3.640, 0.041), "Hf": (3.141, 0.072),
    "Ta": (3.170, 0.081), "W": (3.069, 0.067), "Re": (2.954, 0.066),
    "Os": (3.120, 0.037), "Ir": (2.840, 0.073), "Pt": (2.754, 0.080),
    "Au": (3.293, 0.039), "Hg": (2.705, 0.385), "Tl": (4.347, 0.680),
    "Pb": (4.297, 0.663), "Bi": (4.370, 0.518), "Th": (3.396, 0.026),
    "U": (3.395, 0.022),
}
# DREIDING (Mayo et al., J. Phys. Chem. 1990, 94, 8897): (R0 [Å], D0 [kcal/mol])
_DREIDING_RAW = {
    "H": (3.195, 0.0152), "B": (4.02, 0.095), "C": (3.8983, 0.0951),
    "N": (3.6621, 0.0774), "O": (3.4046, 0.0957), "F": (3.4720, 0.0725),
    "Al": (4.39, 0.31), "Si": (4.27, 0.31), "P": (4.15, 0.32),
    "S": (4.03, 0.344), "Cl": (3.9503, 0.2833), "Br": (3.95, 0.37),
}


def _lj_from_raw(raw: Tuple[float, float]) -> Tuple[float, float]:
    """(distance of the minimum, well depth kcal/mol) -> (epsilon/kB [K], sigma [Å])."""
    return raw[1] * KCAL_TO_K, raw[0] / 2.0 ** (1.0 / 6.0)


def framework_lj(element: str, forcefield: str) -> Tuple[float, float, str]:
    """LJ parameters (epsilon/kB [K], sigma [Å], source) of a framework element.

    `forcefield` contains "DREIDING" or "UFF". DREIDING only covers a few
    elements, so the others (metals in particular) are taken from UFF:
    this is the usual practice for MOFs.
    """
    if "DREIDING" in forcefield.upper() and element in _DREIDING_RAW:
        eps, sig = _lj_from_raw(_DREIDING_RAW[element])
        return eps, sig, "DREIDING"
    if element in _UFF_RAW:
        eps, sig = _lj_from_raw(_UFF_RAW[element])
        return eps, sig, "UFF"
    raise ValueError("No Lennard-Jones parameters for element '%s'." % element)


FORCE_FIELDS = ["TraPPE (adsorbate) + UFF (framework)",
                "TraPPE (adsorbate) + DREIDING/UFF (framework)"]

# -----------------------------------------------------------------------------
# Adsorbates. Each entry gives: default conditions, critical constants
# (RASPA3 derives the fugacity coefficient from them with Peng-Robinson),
# molecular model (sites, LJ, charges), BET cross-sectional area, liquid molar
# volume (Gurvich rule) and geometric probe diameter.
# -----------------------------------------------------------------------------
ADSORBATES: Dict[str, dict] = {
    "N2": {
        "label": "N2 (77 K)", "temperature": 77.0, "p0": 97190.0,
        "Tc": 126.192, "Pc": 3395800.0, "omega": 0.0372,
        "cross_section": 0.162, "liquid_molar_volume": 34.67,
        "probe_diameter": 3.681, "molar_mass": 28.0134,
        "model": "TraPPE (Potoff & Siepmann 2001)",
        "pseudo_atoms": [
            {"name": "N_n2", "element": "N", "mass": 14.00674, "charge": -0.482,
             "epsilon": 36.0, "sigma": 3.31},
            {"name": "N_com", "element": "N", "mass": 0.0, "charge": 0.964,
             "epsilon": 0.0, "sigma": 1.0},
        ],
        "sites": [["N_n2", [0.0, 0.0, 0.55]], ["N_com", [0.0, 0.0, 0.0]],
                  ["N_n2", [0.0, 0.0, -0.55]]],
    },
    "Ar": {
        "label": "Ar (87 K)", "temperature": 87.0, "p0": 98250.0,
        "Tc": 150.687, "Pc": 4863000.0, "omega": -0.00219,
        "cross_section": 0.142, "liquid_molar_volume": 28.6,
        "probe_diameter": 3.40, "molar_mass": 39.948,
        "model": "Lennard-Jones 1 site (eps/kB = 119.8 K, sigma = 3.405 Å)",
        "pseudo_atoms": [
            {"name": "Ar_ads", "element": "Ar", "mass": 39.948, "charge": 0.0,
             "epsilon": 119.8, "sigma": 3.405},
        ],
        "sites": [["Ar_ads", [0.0, 0.0, 0.0]]],
    },
    "CO2": {
        "label": "CO2 (273 K)", "temperature": 273.15, "p0": 3485100.0,
        "Tc": 304.1282, "Pc": 7377300.0, "omega": 0.22394,
        "cross_section": 0.170, "liquid_molar_volume": 47.4,
        "probe_diameter": 3.30, "molar_mass": 44.0095,
        "model": "TraPPE (Potoff & Siepmann 2001)",
        "pseudo_atoms": [
            {"name": "C_co2", "element": "C", "mass": 12.0107, "charge": 0.70,
             "epsilon": 27.0, "sigma": 2.80},
            {"name": "O_co2", "element": "O", "mass": 15.9994, "charge": -0.35,
             "epsilon": 79.0, "sigma": 3.05},
        ],
        "sites": [["O_co2", [0.0, 0.0, 1.16]], ["C_co2", [0.0, 0.0, 0.0]],
                  ["O_co2", [0.0, 0.0, -1.16]]],
    },
}

CHARGE_MODES = ["Charge equilibration (RASPA3, Wilmer-Snurr)",
                "Charges from the CIF (_atom_site_charge)",
                "No framework charges"]

# Simplified CPK colours for the structure picture
ELEMENT_COLORS = {"H": "#e8e8e8", "C": "#505050", "N": "#3050f8", "O": "#ff2010",
                  "F": "#90e050", "Cl": "#1ff01f", "Br": "#a62929", "S": "#e6c820",
                  "P": "#ff8000", "Zn": "#7d80b0", "Cu": "#c88033", "Zr": "#94e0e0",
                  "Co": "#f090a0", "Fe": "#e06633", "Ni": "#50d050", "Al": "#bfa6a6",
                  "Cr": "#8a99c7", "Mg": "#8aff00", "Si": "#f0c8a0", "Ti": "#bfc2c7"}

# Figure palette (readable in greyscale and by colour-blind readers)
COLOR_EXP = "#1f4e9c"
COLOR_SIM = "#d1495b"
COLOR_GEO = "#6b8e23"


# =============================================================================
# 2. GENERAL UTILITIES
# =============================================================================
def to_float(text) -> float:
    """Converts '0.03150(13)', '1,5', '?' ... to float (NaN if impossible)."""
    if text is None:
        return float("nan")
    if isinstance(text, (int, float, np.floating, np.integer)):
        return float(text)
    s = str(text).strip().replace("−", "-")
    s = re.sub(r"\(.*?\)", "", s)                  # strip the standard uncertainty (13)
    if "," in s and "." not in s:
        s = s.replace(",", ".")                    # decimal comma
    try:
        return float(s)
    except ValueError:
        return float("nan")


def element_from_symbol(symbol: str) -> str:
    """'Zn2+', 'ZN', 'C1A' -> 'Zn', 'C'. Returns '' if not recognised."""
    m = re.match(r"([A-Za-z]{1,2})", str(symbol).strip())
    if not m:
        return ""
    s = m.group(1)
    cand = s[0].upper() + s[1:].lower()
    if cand in ATOMIC_MASS:
        return cand
    cand = s[0].upper()
    return cand if cand in ATOMIC_MASS else ""


def hill_formula(counts: Dict[str, float]) -> str:
    """Empirical formula in Hill notation (C, H, then alphabetical order)."""
    keys = sorted(counts)
    if "C" in counts:
        keys = ["C"] + (["H"] if "H" in counts else []) + \
               [k for k in keys if k not in ("C", "H")]
    parts = []
    for k in keys:
        n = counts[k]
        if n <= 1e-6:
            continue
        txt = "%d" % round(n) if abs(n - round(n)) < 1e-3 else ("%.2f" % n).rstrip("0").rstrip(".")
        parts.append(k + ("" if txt == "1" else txt))
    return " ".join(parts)


def pct_diff(value: float, reference: float) -> float:
    """Relative difference (%) of `value` with respect to `reference`."""
    if reference is None or value is None:
        return float("nan")
    if not np.isfinite(reference) or not np.isfinite(value) or reference == 0:
        return float("nan")
    return 100.0 * (value - reference) / reference


def fmt(value, digits: int = 2, unit: str = "") -> str:
    """Tolerant formatting (None/NaN -> '—')."""
    try:
        if value is None or not np.isfinite(float(value)):
            return "—"
        return ("%." + str(digits) + "f") % float(value) + ((" " + unit) if unit else "")
    except (TypeError, ValueError):
        return str(value)


# =============================================================================
# 3. INTERNAL CIF READER (documented alternative to the CSD Python API)
# =============================================================================
class CIFParser:
    """Minimal but robust CIF parser (CIF1: SHELX, Olex2, CSD, RASPA...).

    It does not replace the CSD Python API: it is used (a) outside Mercury and
    (b) to re-read the original file to obtain occupancies and disorder groups,
    which are needed to build a clean P1 model.
    Only the first `data_` block that contains atoms is used.
    """

    def __init__(self, text: str):
        self.items: Dict[str, str] = {}
        self.loops: List[Tuple[List[str], List[List[str]]]] = []
        self._parse(text)

    # -- tokenisation ----------------------------------------------------------
    @staticmethod
    def _tokenize(text: str) -> List[Tuple[str, str]]:
        """Returns (type, value) pairs; type = 'tok' or 'str'."""
        tokens: List[Tuple[str, str]] = []
        lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.startswith(";"):                       # multi-line text field
                buf = [line[1:]]
                i += 1
                while i < len(lines) and not lines[i].startswith(";"):
                    buf.append(lines[i])
                    i += 1
                tokens.append(("str", "\n".join(buf)))
                rest = lines[i][1:] if i < len(lines) else ""
                i += 1
                line = rest
                if not line.strip():
                    continue
            pos, n = 0, len(line)
            while pos < n:
                ch = line[pos]
                if ch in " \t":
                    pos += 1
                elif ch == "#":
                    break
                elif ch in "'\"":
                    end = pos + 1
                    while True:                            # closing quote = followed by blank
                        end = line.find(ch, end)
                        if end == -1:
                            end = n
                            break
                        if end + 1 >= n or line[end + 1] in " \t":
                            break
                        end += 1
                    tokens.append(("str", line[pos + 1:end]))
                    pos = end + 1
                else:
                    end = pos
                    while end < n and line[end] not in " \t":
                        end += 1
                    tokens.append(("tok", line[pos:end]))
                    pos = end
            i += 1
        return tokens

    def _parse(self, text: str) -> None:
        tokens = self._tokenize(text)
        i, n = 0, len(tokens)
        blocks_seen = 0
        while i < n:
            kind, val = tokens[i]
            low = val.lower() if kind == "tok" else ""
            if kind == "tok" and low.startswith("data_"):
                blocks_seen += 1
                if blocks_seen > 1 and self._has_atoms():
                    break                                   # first useful block read
                i += 1
            elif kind == "tok" and low == "loop_":
                i += 1
                tags = []
                while i < n and tokens[i][0] == "tok" and tokens[i][1].startswith("_"):
                    tags.append(tokens[i][1].lower().replace("_atom_site.", "_atom_site_"))
                    i += 1
                values = []
                while i < n:
                    k, v = tokens[i]
                    if k == "tok" and (v.startswith("_") or v.lower() == "loop_"
                                       or v.lower().startswith("data_")):
                        break
                    values.append(v)
                    i += 1
                if tags:
                    rows = [values[j:j + len(tags)]
                            for j in range(0, len(values) - len(tags) + 1, len(tags))]
                    self.loops.append((tags, rows))
            elif kind == "tok" and val.startswith("_"):
                if i + 1 < n and not (tokens[i + 1][0] == "tok" and tokens[i + 1][1].startswith("_")):
                    self.items[val.lower()] = tokens[i + 1][1]
                    i += 2
                else:
                    i += 1
            else:
                i += 1

    def _has_atoms(self) -> bool:
        return any("_atom_site_fract_x" in tags for tags, _ in self.loops)

    # -- access --------------------------------------------------------------
    def get(self, *names: str, default: Optional[str] = None) -> Optional[str]:
        for name in names:
            if name.lower() in self.items:
                return self.items[name.lower()]
        return default

    def loop_with(self, tag: str) -> Optional[Tuple[List[str], List[List[str]]]]:
        for tags, rows in self.loops:
            if tag in tags:
                return tags, rows
        return None

    def symmetry_operators(self) -> List[str]:
        for tag in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz",
                    "_space_group_symop.operation_xyz"):
            found = self.loop_with(tag)
            if found:
                tags, rows = found
                k = tags.index(tag)
                return [r[k] for r in rows]
            if tag in self.items:
                return [self.items[tag]]
        return []


# =============================================================================
# 4. STRUCTURE MODEL
# =============================================================================
@dataclass
class Atom:
    """Atom of the asymmetric unit (or of the P1 cell after expansion)."""
    label: str
    element: str
    frac: Tuple[float, float, float]
    occupancy: float = 1.0
    disorder_group: str = "."
    charge: Optional[float] = None


def parse_symop(op: str) -> Tuple[np.ndarray, np.ndarray]:
    """'-x+1/2, y, z+0.5' -> (3x3 rotation matrix, translation)."""
    rot = np.zeros((3, 3))
    trans = np.zeros(3)
    parts = op.lower().replace(" ", "").split(",")
    if len(parts) != 3:
        raise ValueError("Unreadable symmetry operator: %r" % op)
    for i, part in enumerate(parts):
        for sign, num, den, axis in re.findall(
                r"([+-]?)(\d*\.?\d*)(?:/(\d+))?\*?([xyz]?)", part):
            if not num and not axis:
                continue
            sgn = -1.0 if sign == "-" else 1.0
            val = float(num) if num else 1.0
            if den:
                val /= float(den)
            if axis:
                rot[i, "xyz".index(axis)] += sgn * val
            else:
                trans[i] += sgn * val
    return rot, trans


class StructureData:
    """Crystal structure, independent of any external library.

    It can be built from a CSD Python API entry (`from_csd_entry`)
    or from a CIF file (`from_cif_file`).
    """

    def __init__(self):
        self.name = "structure"
        self.source_file = ""
        self.cell = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)      # a b c alpha beta gamma
        self.spacegroup = "P 1"
        self.symops: List[str] = ["x,y,z"]
        self.atoms: List[Atom] = []
        self.formula_reported = ""                         # formula from the CIF / the CSD
        self.formula_weight_reported = float("nan")
        self.z_value = float("nan")
        self.density_reported = float("nan")
        self.notes: List[str] = []                         # warnings
        self.csd_cif_path = ""                             # CIF readable by the CSD Python API (Pore Analyser)
        # --- result of the P1 expansion (build_p1) ---
        self.p1_elements: List[str] = []
        self.p1_frac = np.zeros((0, 3))
        self.p1_charges: Optional[np.ndarray] = None
        self.p1_removed: Dict[str, int] = {}
        self.p1_log: List[str] = []

    # -------------------------------------------------------------------- cell
    @property
    def cell_matrix(self) -> np.ndarray:
        """3x3 matrix whose ROWS are the a, b, c vectors (Å)."""
        a, b, c, al, be, ga = self.cell
        al, be, ga = np.radians([al, be, ga])
        cx = c * math.cos(be)
        cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
        cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
        return np.array([[a, 0.0, 0.0],
                         [b * math.cos(ga), b * math.sin(ga), 0.0],
                         [cx, cy, cz]])

    @property
    def volume(self) -> float:
        return float(abs(np.linalg.det(self.cell_matrix)))

    @property
    def perpendicular_widths(self) -> np.ndarray:
        """Distances between opposite faces of the cell (Å)."""
        m = self.cell_matrix
        v = self.volume
        return np.array([v / np.linalg.norm(np.cross(m[1], m[2])),
                         v / np.linalg.norm(np.cross(m[2], m[0])),
                         v / np.linalg.norm(np.cross(m[0], m[1]))])

    def minimal_supercell(self, cutoff: float = 12.0) -> Tuple[int, int, int]:
        """Smallest supercell that satisfies the minimum image convention."""
        return tuple(int(math.ceil(2.0 * cutoff / w - 1e-9)) for w in self.perpendicular_widths)

    # -------------------------------------------------------------- construction
    @classmethod
    def from_cif_file(cls, path: str) -> "StructureData":
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
        st = cls.from_cif_text(text)
        st.source_file = path
        st.name = os.path.splitext(os.path.basename(path))[0]
        return st

    @classmethod
    def from_cif_text(cls, text: str) -> "StructureData":
        cif = CIFParser(text)
        st = cls()
        try:
            st.cell = tuple(to_float(cif.get("_cell_length_" + k)) for k in "abc") + \
                      tuple(to_float(cif.get("_cell_angle_" + k)) for k in ("alpha", "beta", "gamma"))
        except Exception as exc:                            # pragma: no cover
            raise ValueError("Cell parameters missing from the CIF.") from exc
        if not all(np.isfinite(st.cell)):
            raise ValueError("Cell parameters missing or unreadable in the CIF.")
        st.spacegroup = (cif.get("_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m",
                                 "_space_group_name_hall", default="?") or "?").strip()
        ops = cif.symmetry_operators()
        if ops:
            st.symops = ops
        else:
            st.notes.append("No symmetry operators in the CIF: structure treated as P1.")
        st.formula_reported = (cif.get("_chemical_formula_sum", default="") or "").strip()
        st.formula_weight_reported = to_float(cif.get("_chemical_formula_weight"))
        st.z_value = to_float(cif.get("_cell_formula_units_z"))
        st.density_reported = to_float(cif.get("_exptl_crystal_density_diffrn"))
        found = cif.loop_with("_atom_site_fract_x")
        if not found:
            raise ValueError("No _atom_site_fract_x loop in the CIF.")
        tags, rows = found
        col = {t: i for i, t in enumerate(tags)}

        def cell_value(row, tag, default=None):
            return row[col[tag]] if tag in col else default

        for row in rows:
            label = cell_value(row, "_atom_site_label", "X")
            element = element_from_symbol(cell_value(row, "_atom_site_type_symbol", "") or "") \
                or element_from_symbol(label)
            if not element:
                st.notes.append("Atom ignored (unknown element): %s" % label)
                continue
            frac = tuple(to_float(cell_value(row, "_atom_site_fract_" + k)) for k in "xyz")
            if not all(np.isfinite(frac)):
                continue
            occ = to_float(cell_value(row, "_atom_site_occupancy", "1"))
            charge = to_float(cell_value(row, "_atom_site_charge")) \
                if "_atom_site_charge" in col else None
            st.atoms.append(Atom(label, element, frac,
                                 1.0 if not np.isfinite(occ) else occ,
                                 str(cell_value(row, "_atom_site_disorder_group", ".")),
                                 charge if charge is None or np.isfinite(charge) else None))
        if not st.atoms:
            raise ValueError("No usable atoms in the CIF.")
        return st

    @classmethod
    def from_csd_entry(cls, entry, original_file: str = "") -> "StructureData":
        """Builds the structure from a `ccdc.entry.Entry` object.

        Strategy (most faithful first):
          1. if Mercury loaded a CIF file whose path is known, the atoms are
             re-read from that file (occupancies, disorder);
          2. otherwise, atoms of `crystal.asymmetric_unit_molecule`;
          3. otherwise, the CIF exported by `crystal.to_string('cif')` is parsed.
        Metadata (formula, space group...) come from the API.
        """
        crystal = entry.crystal
        st: Optional[StructureData] = None
        if original_file and original_file.lower().endswith(".cif") and os.path.isfile(original_file):
            try:
                st = cls.from_cif_file(original_file)
            except Exception:
                st = None
        if st is None:
            try:
                st = cls()
                lengths, angles = crystal.cell_lengths, crystal.cell_angles
                st.cell = (lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2])
                st.symops = [str(o) for o in crystal.symmetry_operators] or ["x,y,z"]
                try:
                    molecule = crystal.asymmetric_unit_molecule
                except Exception:
                    molecule = crystal.molecule
                for at in molecule.atoms:
                    fc = at.fractional_coordinates
                    if fc is None:
                        continue
                    element = element_from_symbol(at.atomic_symbol)
                    if not element:
                        continue
                    occ = getattr(at, "occupancy", 1.0)
                    st.atoms.append(Atom(str(at.label), element,
                                         (float(fc[0]), float(fc[1]), float(fc[2])),
                                         1.0 if occ is None else float(occ)))
                if not st.atoms:
                    raise ValueError("no atoms with coordinates")
            except Exception:
                st = cls.from_cif_text(crystal.to_string("cif"))
        # Metadata from the API (each access is guarded: the available
        # attributes vary slightly between CSD Python API versions).
        for setter in (
            lambda: setattr(st, "name", str(entry.identifier)),
            lambda: setattr(st, "spacegroup", str(crystal.spacegroup_symbol)),
            lambda: setattr(st, "formula_reported", str(entry.formula or crystal.formula)),
            lambda: setattr(st, "z_value", float(crystal.z_value)),
            lambda: setattr(st, "density_reported", float(crystal.calculated_density)),
            lambda: setattr(st, "formula_weight_reported", float(entry.molecule.molecular_weight)),
        ):
            try:
                setter()
            except Exception:
                pass
        return st

    # ------------------------------------------------------------- descriptors
    def asymmetric_summary(self) -> Dict[str, float]:
        """Unit-cell contents weighted by occupancy (before cleaning)."""
        elements, frac, occ, _ = self._expand(dedupe=True)
        counts: Dict[str, float] = {}
        for el, o in zip(elements, occ):
            counts[el] = counts.get(el, 0.0) + o
        return counts

    @staticmethod
    def mass_of(counts: Dict[str, float]) -> float:
        return float(sum(ATOMIC_MASS[el] * n for el, n in counts.items()))

    def density_from_mass(self, mass: float) -> float:
        """Cell mass (g/mol) -> density (g/cm3)."""
        return mass / N_AVOGADRO / (self.volume * 1e-24)

    # ---------------------------------------------------------------- expansion
    def _expand(self, dedupe: bool = True):
        """Applies all symmetry operators.

        Returns (elements, fractional coordinates in [0,1[, occupancies,
        index of the parent atom). Duplicates created by special positions
        (same parent atom closer than 0.2 Å) are merged.
        """
        ops = [parse_symop(o) for o in self.symops]
        cellm = self.cell_matrix
        elements, fracs, occs, parents = [], [], [], []
        for idx, atom in enumerate(self.atoms):
            f0 = np.array(atom.frac, dtype=float)
            images = np.array([(rot @ f0 + trans) % 1.0 for rot, trans in ops])
            if dedupe:
                keep: List[np.ndarray] = []
                for img in images:
                    if keep:
                        d = np.array(keep) - img
                        d -= np.round(d)
                        if np.min(np.linalg.norm(d @ cellm, axis=1)) < 0.2:
                            continue
                    keep.append(img)
                images = np.array(keep)
            for img in images:
                elements.append(atom.element)
                fracs.append(img)
                occs.append(atom.occupancy)
                parents.append(idx)
        return elements, np.array(fracs).reshape(-1, 3) % 1.0, np.array(occs), np.array(parents)

    def neighbor_pairs(self, frac: np.ndarray, rmax: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Pairs (i<j) of atoms closer than `rmax` (minimum image).

        Block-wise calculation in pure numpy (no dependency), valid for the
        short distances used here (bonds, disorder clashes).
        """
        cellm = self.cell_matrix
        n = len(frac)
        out_i, out_j, out_d = [], [], []
        block = max(1, int(4e6 // max(n, 1)))
        for start in range(0, n, block):
            d = frac[start:start + block, None, :] - frac[None, :, :]
            d -= np.round(d)
            dist = np.linalg.norm(d @ cellm, axis=2)
            ii, jj = np.nonzero(dist < rmax)
            ii = ii + start
            mask = ii < jj
            out_i.append(ii[mask])
            out_j.append(jj[mask])
            out_d.append(dist[ii[mask] - start, jj[mask]])
        return np.concatenate(out_i), np.concatenate(out_j), np.concatenate(out_d)

    def build_p1(self, resolve_disorder: bool = True, remove_guests: bool = True,
                 guest_fraction: float = 0.05) -> None:
        """Builds the P1 model used by RASPA3 and by the geometric calculations.

        resolve_disorder: among partially occupied sites that clash sterically
            (two components of the same disordered group), only one position
            is kept (greedy algorithm, decreasing occupancy).
            Example: ZIF-8 methyl described by 6 half-H -> 3 H.
        remove_guests: removes fragments not bonded to the framework whose
            size is below `guest_fraction` times that of the largest fragment
            (solvent molecules, pseudo-located water oxygens...).
        """
        elements, frac, occ, parents = self._expand(dedupe=True)
        n0 = len(elements)
        keep = np.ones(n0, dtype=bool)
        log: List[str] = ["Raw P1 cell: %d sites" % n0]
        removed: Dict[str, int] = {}
        radii = np.array([COVALENT_RADIUS.get(e, 1.5) for e in elements])

        # --- 1) disorder resolution ---------------------------------------
        partial = occ < 0.999
        if resolve_disorder and partial.any():
            pi, pj, pd_ = self.neighbor_pairs(frac, 2.6)
            conflicts: Dict[int, List[int]] = {}
            for i, j, d in zip(pi, pj, pd_):
                if not (partial[i] or partial[j]):
                    continue
                both_h = elements[i] == "H" and elements[j] == "H"
                limit = 1.2 if both_h else 0.7 * (radii[i] + radii[j])
                if d < limit:
                    conflicts.setdefault(int(i), []).append(int(j))
                    conflicts.setdefault(int(j), []).append(int(i))
            order = sorted(np.nonzero(partial)[0], key=lambda k: (-occ[k], k))
            accepted = ~partial                              # fully occupied sites: always kept
            for k in order:
                if any(accepted[m] for m in conflicts.get(int(k), [])):
                    keep[k] = False
                    removed["disorder:" + self.atoms[parents[k]].label] = \
                        removed.get("disorder:" + self.atoms[parents[k]].label, 0) + 1
                else:
                    accepted[k] = True
            log.append("Disorder: %d partial sites, %d positions discarded"
                       % (int(partial.sum()), int((~keep).sum())))
        elif partial.any():
            log.append("WARNING: %d partially occupied sites kept as they are."
                       % int(partial.sum()))

        # --- 2) removal of guest species -------------------------------
        if remove_guests:
            idx = np.nonzero(keep)[0]
            sub = frac[idx]
            pi, pj, pd_ = self.neighbor_pairs(sub, 3.2)
            parent = list(range(len(idx)))

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x

            sub_r = radii[idx]
            bonded = pd_ < 1.25 * (sub_r[pi] + sub_r[pj])
            for i, j in zip(pi[bonded], pj[bonded]):
                ri, rj = find(int(i)), find(int(j))
                if ri != rj:
                    parent[ri] = rj
            roots = np.array([find(k) for k in range(len(idx))])
            uniq, counts = np.unique(roots, return_counts=True)
            biggest = counts.max() if len(counts) else 0
            small_roots = set(uniq[counts < guest_fraction * biggest].tolist())
            n_guest = 0
            for local, root in enumerate(roots):
                if root in small_roots:
                    k = idx[local]
                    keep[k] = False
                    n_guest += 1
                    key = "guest:" + self.atoms[parents[k]].label
                    removed[key] = removed.get(key, 0) + 1
            log.append("Guest species: %d atoms removed (%d fragments)"
                       % (n_guest, len(small_roots)))

        self.p1_elements = [e for e, k in zip(elements, keep) if k]
        self.p1_frac = frac[keep]
        charges = [self.atoms[p].charge for p in parents[keep]]
        self.p1_charges = np.array(charges, dtype=float) if charges and all(c is not None for c in charges) else None
        self.p1_removed = removed
        log.append("Final P1 model: %d atoms, %s" % (len(self.p1_elements),
                                                         hill_formula(self.p1_counts())))
        self.p1_log = log

    def removed_summary(self) -> str:
        """Readable summary of the discarded atoms, by reason and element."""
        groups: Dict[Tuple[str, str], int] = {}
        for key, count in self.p1_removed.items():
            reason, label = key.split(":", 1)
            groups[(reason, element_from_symbol(label))] = groups.get((reason, element_from_symbol(label)), 0) + count
        return "; ".join("%s: %d %s" % (reason, count, el) for (reason, el), count in sorted(groups.items())) \
            or "none"

    def p1_counts(self) -> Dict[str, float]:
        counts: Dict[str, float] = {}
        for el in self.p1_elements:
            counts[el] = counts.get(el, 0.0) + 1.0
        return counts

    @property
    def p1_mass(self) -> float:
        """Molar mass of the cleaned P1 cell (g per mol of cells)."""
        return self.mass_of(self.p1_counts())

    @property
    def p1_cart(self) -> np.ndarray:
        return self.p1_frac @ self.cell_matrix

    def write_p1_cif(self, path: str, with_charges: bool = False) -> None:
        """Writes the P1 CIF read by RASPA3 (label = type = element symbol)."""
        a, b, c, al, be, ga = self.cell
        lines = ["data_framework", "",
                 "_audit_creation_method '%s %s'" % (APP_NAME, APP_VERSION),
                 "_cell_length_a %.6f" % a, "_cell_length_b %.6f" % b, "_cell_length_c %.6f" % c,
                 "_cell_angle_alpha %.5f" % al, "_cell_angle_beta %.5f" % be,
                 "_cell_angle_gamma %.5f" % ga, "_cell_volume %.4f" % self.volume, "",
                 "_symmetry_cell_setting triclinic",
                 "_symmetry_space_group_name_Hall 'P 1'",
                 "_symmetry_space_group_name_H-M 'P 1'",
                 "_symmetry_Int_Tables_number 1", "",
                 "loop_", "_symmetry_equiv_pos_as_xyz", " 'x,y,z'", "",
                 "loop_", "_atom_site_label", "_atom_site_type_symbol",
                 "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"]
        use_q = with_charges and self.p1_charges is not None
        if use_q:
            lines.append("_atom_site_charge")
        for k, (el, f) in enumerate(zip(self.p1_elements, self.p1_frac)):
            row = "%-6s %-3s %.6f %.6f %.6f" % ("%s%d" % (el, k + 1), el, f[0], f[1], f[2])
            if use_q:
                row += " %.5f" % self.p1_charges[k]
            lines.append(row)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")

    def write_transfer_cif(self, path: str) -> None:
        """Writes the complete structure (symmetry, occupancies, disorder, charges).

        This CIF passes the structure from the parent process (CSD Python API)
        to the GUI process; it is read back with `from_cif_file`.
        """
        def num(value: float) -> str:
            return "%.6g" % value if np.isfinite(value) else "?"

        a, b, c, al, be, ga = self.cell
        lines = ["data_transfer",
                 "_chemical_formula_sum '%s'" % (self.formula_reported.replace("'", " ") or "?"),
                 "_chemical_formula_weight %s" % num(self.formula_weight_reported),
                 "_cell_formula_units_Z %s" % num(self.z_value),
                 "_exptl_crystal_density_diffrn %s" % num(self.density_reported),
                 "_space_group_name_H-M_alt '%s'" % (self.spacegroup.replace("'", " ") or "?"),
                 "_cell_length_a %.6f" % a, "_cell_length_b %.6f" % b, "_cell_length_c %.6f" % c,
                 "_cell_angle_alpha %.5f" % al, "_cell_angle_beta %.5f" % be, "_cell_angle_gamma %.5f" % ga,
                 "loop_", "_space_group_symop_operation_xyz"]
        lines += ["'%s'" % op.replace("'", "") for op in self.symops]
        with_q = bool(self.atoms) and all(at.charge is not None for at in self.atoms)
        lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
                  "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy",
                  "_atom_site_disorder_group"] + (["_atom_site_charge"] if with_q else [])
        for k, at in enumerate(self.atoms):
            label = re.sub(r"[^A-Za-z0-9_+*-]", "", at.label) or "%s%d" % (at.element, k + 1)
            group = re.sub(r"\s", "", str(at.disorder_group)) or "."
            row = "%s %s %.7f %.7f %.7f %.5f %s" % (label, at.element, at.frac[0], at.frac[1], at.frac[2],
                                                    at.occupancy, group)
            lines.append(row + (" %.6f" % at.charge if with_q else ""))
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")


# =============================================================================
# 5. GEOMETRIC PROPERTIES (accessible surface area, pore volume, PLD, LCD)
# =============================================================================
PORE_MARKER = "PORE_ANALYSER_JSON:"       # prefix of the result line of the helper process


def pore_helper_main(cif_path: str) -> int:
    """Helper process: CSD Python API Pore Analyser on a CIF.

    Prints one line `PORE_ANALYSER_JSON:{...}` containing the fields of a
    `GeometricResult` (or {"error": ...}). Never imports Qt.
    """
    try:
        from ccdc import io as ccdc_io
        from ccdc.descriptors import CrystalDescriptors
        crystal = ccdc_io.CrystalReader(cif_path)[0]
        analyser = CrystalDescriptors.PoreAnalyser(crystal)
        mass = float(analyser.system_mass)                   # g per mol of cells
        notes = []
        try:
            probe = float(analyser.settings.n_probe_sigma)
        except Exception:
            probe = 3.314                                    # CSD default value
        surface = float(analyser.network_accessible_surface_area_per_mass)
        if not surface > 0:
            surface = float(analyser.total_surface_area_per_mass)
            notes.append("Total surface area (network-accessible surface area is zero).")
        to_cm3_g = 1e-24 * N_AVOGADRO / mass
        data = {"probe_diameter": probe, "surface_m2_g": surface,
                "surface_m2_cm3": float(analyser.total_surface_area_per_volume),
                "geometric_volume_cm3_g": float(analyser.total_geometric_volume) * to_cm3_g,
                "pore_volume_cm3_g": float(analyser.total_helium_volume) * to_cm3_g,
                "pore_volume_label": "helium volume (Pore Analyser)",
                "void_fraction": float(analyser.total_geometric_volume) / float(analyser.system_volume),
                "pld": float(analyser.pore_limiting_diameter), "lcd": float(analyser.max_pore_diameter),
                "notes": notes}
    except BaseException as exc:
        data = {"error": "%s: %s" % (type(exc).__name__, exc)}
    print(PORE_MARKER + json.dumps(data))
    return 0 if "error" not in data else 1


@dataclass
class GeometricResult:
    """Geometric results, whichever method was used."""
    method: str = ""
    probe_diameter: float = float("nan")      # Å
    surface_m2_g: float = float("nan")        # accessible surface area (m2/g)
    surface_m2_cm3: float = float("nan")      # accessible surface area (m2/cm3)
    pore_volume_cm3_g: float = float("nan")   # pore volume used for the comparison
    pore_volume_label: str = ""
    geometric_volume_cm3_g: float = float("nan")   # point probe
    void_fraction: float = float("nan")
    pld: float = float("nan")                 # pore limiting diameter (Å)
    lcd: float = float("nan")                 # largest cavity diameter (Å)
    notes: List[str] = field(default_factory=list)


class GeometricAnalyzer:
    """Geometric calculation on the cleaned P1 model.

    Two routes:
      * `run_csd()`: CSD Python API Pore Analyser
        (ccdc.descriptors.CrystalDescriptors.PoreAnalyser), run in a helper
        process and used only if the structure was not modified by the
        cleaning step (otherwise the crystal seen by the CSD would still
        contain the solvent);
      * `run_internal()`: equivalent implementation, documented below.

    Internal implementation
    -----------------------
    Radius of a framework atom: sigma_UFF / 2 (as in Pore Analyser, Zeo++
    and RASPA). Probe: sphere of diameter `probe_diameter`.
      * Accessible surface area (Düren et al., J. Phys. Chem. C 2007, 111,
        15350): points are spread over the sphere of radius
        (sigma_i + d_probe)/2 centred on each atom; the fraction of points
        that do not overlap any other atom gives the accessible area.
      * 3D grid of the distance D to the surface of the nearest atom:
          - geometric volume         = fraction of nodes where D > 0;
          - probe-occupiable volume  = dilation, by a sphere of radius
            r_probe, of the region accessible to the probe centre
            (D >= r_probe) - definition of Ongari et al., Langmuir 2017;
          - LCD = 2 max(D);
          - PLD = 2 r_max, r_max being the largest radius for which the
            region D >= r percolates through the periodic cell
            (bisection + breadth-first search tracking periodic images).
        The PLD/LCD accuracy is of the order of the grid spacing.
    """

    MARGIN = 14.0   # Å: thickness of the periodic-image shell

    def __init__(self, structure: StructureData, probe_diameter: float,
                 log: Optional[Callable[[str], None]] = None):
        self.st = structure
        self.probe_d = float(probe_diameter)
        self.log = log or (lambda msg: None)

    # ------------------------------------------------------------------ CSD API
    def run_csd(self) -> GeometricResult:
        """Gets the Pore Analyser values (CSD Python API >= 3.0.15).

        The calculation runs in a separate Python process (internal option
        --pore-helper): the CSD Python API must not be loaded into the Qt GUI
        process (Qt library conflict).
        """
        command = [sys.executable, os.path.abspath(__file__), "--pore-helper", self.st.csd_cif_path]
        proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, timeout=3600)
        data = None
        for line in proc.stdout.splitlines():
            if line.startswith(PORE_MARKER):
                data = json.loads(line[len(PORE_MARKER):])
        if data is None:
            raise RuntimeError("no answer from the Pore Analyser: %s" % (proc.stderr.strip()[-300:] or "?"))
        if "error" in data:
            raise RuntimeError(data["error"])
        return GeometricResult(method="CSD Python API - Pore Analyser", **data)

    # --------------------------------------------------------------- internal tools
    def _atoms_with_images(self) -> Tuple[np.ndarray, np.ndarray]:
        """P1 atoms + periodic images within a MARGIN Å shell."""
        st = self.st
        frac = st.p1_frac
        radii = np.array([framework_lj(e, "UFF")[1] / 2.0 for e in st.p1_elements])
        margin_frac = self.MARGIN / st.perpendicular_widths
        reach = np.ceil(margin_frac).astype(int)
        all_frac, all_rad = [], []
        for i in range(-reach[0], reach[0] + 1):
            for j in range(-reach[1], reach[1] + 1):
                for k in range(-reach[2], reach[2] + 1):
                    shifted = frac + np.array([i, j, k])
                    mask = np.all((shifted > -margin_frac) & (shifted < 1.0 + margin_frac), axis=1)
                    all_frac.append(shifted[mask])
                    all_rad.append(radii[mask])
        return np.vstack(all_frac) @ st.cell_matrix, np.concatenate(all_rad)

    @staticmethod
    def _surface_distance(points: np.ndarray, atoms: np.ndarray, radii: np.ndarray) -> np.ndarray:
        """min_j(|p - a_j| - R_j) for each point (distance to the atomic surface).

        With scipy: k-d tree (fast). Without scipy: block-wise brute force.
        """
        try:
            from scipy.spatial import cKDTree
            tree = cKDTree(atoms)
            k = min(32, len(atoms))
            out = np.empty(len(points))
            for start in range(0, len(points), 200000):
                dist, idx = tree.query(points[start:start + 200000], k=k)
                if k == 1:
                    dist, idx = dist[:, None], idx[:, None]
                out[start:start + 200000] = np.min(dist - radii[idx], axis=1)
            return out
        except ImportError:
            out = np.empty(len(points))
            block = max(1, int(3e6 // max(len(atoms), 1)))
            for start in range(0, len(points), block):
                d = np.linalg.norm(points[start:start + block, None, :] - atoms[None, :, :], axis=2)
                out[start:start + block] = np.min(d - radii[None, :], axis=1)
            return out

    @staticmethod
    def _have_scipy() -> bool:
        try:
            import scipy.spatial  # noqa: F401
            return True
        except ImportError:
            return False

    @staticmethod
    def _percolates(mask: np.ndarray) -> bool:
        """True if region `mask` (periodic grid) is connected to infinity.

        Principle: each node is tagged with the cell image in which it was
        reached; reaching a known node again through a different image proves
        that a channel runs through the crystal (26 neighbours).
          * with scipy: cluster labelling inside the cell (ndimage.label),
            then weighted union-find over the links crossing the cell faces;
          * without scipy: breadth-first search in pure Python (slower).
        """
        shape = mask.shape
        if not mask.any():
            return False
        neighbours = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)
                      if (i, j, k) != (0, 0, 0)]
        try:
            from scipy import ndimage
        except ImportError:
            ndimage = None
        if ndimage is not None:
            labels, count = ndimage.label(mask, structure=np.ones((3, 3, 3), dtype=int))
            idx = np.indices(shape)
            edges = set()
            for d in neighbours:
                if d < (0, 0, 0):
                    continue                                  # each pair only once
                other = np.roll(labels, shift=(-d[0], -d[1], -d[2]), axis=(0, 1, 2))
                crossing = np.zeros(shape, dtype=bool)
                wraps = []
                for axis in range(3):
                    w = (idx[axis] + d[axis]) // shape[axis]   # -1, 0 or +1
                    wraps.append(w)
                    crossing |= w != 0
                sel = crossing & (labels > 0) & (other > 0)
                if sel.any():
                    rows = np.stack([labels[sel], other[sel], wraps[0][sel], wraps[1][sel],
                                     wraps[2][sel]], axis=1)
                    edges.update(map(tuple, np.unique(rows, axis=0).tolist()))
            parent = list(range(count + 1))
            offset = [(0, 0, 0)] * (count + 1)                 # image relative to the root

            def find(x):
                path = []
                while parent[x] != x:
                    path.append(x)
                    x = parent[x]
                total = (0, 0, 0)
                for node in reversed(path):                    # path compression
                    o = offset[node]
                    total = (total[0] + o[0], total[1] + o[1], total[2] + o[2])
                    parent[node], offset[node] = x, total
                return x

            for la, lb, wx, wy, wz in edges:
                ra, rb = find(la), find(lb)
                oa, ob = offset[la], offset[lb]
                # image(lb) = image(la) + w
                want = (oa[0] + wx - ob[0], oa[1] + wy - ob[1], oa[2] + wz - ob[2])
                if ra == rb:
                    if want != (0, 0, 0):
                        return True
                else:
                    parent[rb], offset[rb] = ra, want
            return False
        nodes = set(map(tuple, np.argwhere(mask)))
        image: Dict[Tuple[int, int, int], Tuple[int, int, int]] = {}
        for seed in nodes:
            if seed in image:
                continue
            image[seed] = (0, 0, 0)
            queue = deque([seed])
            while queue:
                cur = queue.popleft()
                ci = image[cur]
                for d in neighbours:
                    raw = (cur[0] + d[0], cur[1] + d[1], cur[2] + d[2])
                    wrapped = (raw[0] % shape[0], raw[1] % shape[1], raw[2] % shape[2])
                    if wrapped not in nodes:
                        continue
                    img = (ci[0] + raw[0] // shape[0], ci[1] + raw[1] // shape[1],
                           ci[2] + raw[2] // shape[2])
                    known = image.get(wrapped)
                    if known is None:
                        image[wrapped] = img
                        queue.append(wrapped)
                    elif known != img:
                        return True
        return False

    # --------------------------------------------------------- internal calculation
    def run_internal(self, samples_per_atom: int = 300, grid_spacing: float = 0.25) -> GeometricResult:
        st = self.st
        if not st.p1_elements:
            raise ValueError("Empty P1 model: build the structure first.")
        res = GeometricResult(method="Internal calculation (Monte Carlo + grid)",
                              probe_diameter=self.probe_d)
        rp = self.probe_d / 2.0
        cellm = st.cell_matrix
        atoms, radii = self._atoms_with_images()
        own_radii = np.array([framework_lj(e, "UFF")[1] / 2.0 for e in st.p1_elements])
        n_atoms = len(st.p1_elements)
        fast = self._have_scipy()
        if not fast:
            res.notes.append("scipy missing: brute-force calculation (reduced sampling).")
            samples_per_atom = int(np.clip(1.5e9 / (n_atoms * len(atoms)), 30, samples_per_atom))
            grid_spacing = max(grid_spacing, (st.volume / max(1.5e9 / len(atoms), 1.0)) ** (1 / 3))

        # --- accessible surface area: Fibonacci sphere around each atom ---
        self.log("Accessible surface area: %d points per atom..." % samples_per_atom)
        k = np.arange(samples_per_atom) + 0.5
        phi = np.arccos(1.0 - 2.0 * k / samples_per_atom)
        theta = math.pi * (1.0 + 5.0 ** 0.5) * k
        unit = np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)
        centers = st.p1_cart
        shell = own_radii + rp
        points = (centers[:, None, :] + shell[:, None, None] * unit[None, :, :]).reshape(-1, 3)
        inv = np.linalg.inv(cellm)
        points = ((points @ inv) % 1.0) @ cellm                    # wrap into the cell
        dist = self._surface_distance(points, atoms, radii).reshape(n_atoms, samples_per_atom)
        accessible = (dist >= rp - 1e-6).mean(axis=1)
        area = float(np.sum(4.0 * math.pi * shell ** 2 * accessible))   # Å2 per cell
        res.surface_m2_g = area * 1e-20 * N_AVOGADRO / st.p1_mass
        res.surface_m2_cm3 = area * 1e-20 / (st.volume * 1e-24)

        # --- distance grid ---
        lengths = np.linalg.norm(cellm, axis=1)
        max_nodes = 3.0e6 if fast else 2.0e5
        spacing = max(grid_spacing, (st.volume / max_nodes) ** (1 / 3))
        shape = tuple(int(max(4, math.ceil(l / spacing))) for l in lengths)
        self.log("Distance grid: %d x %d x %d nodes (spacing ~%.2f Å)..."
                 % (shape[0], shape[1], shape[2], spacing))
        axes = [(np.arange(n) + 0.5) / n for n in shape]
        gf = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
        dgrid = self._surface_distance(gf @ cellm, atoms, radii).reshape(shape)
        to_cm3_g = st.volume * 1e-24 * N_AVOGADRO / st.p1_mass       # fraction -> cm3/g
        res.void_fraction = float((dgrid > 0).mean())
        res.geometric_volume_cm3_g = res.void_fraction * to_cm3_g
        dmax = float(dgrid.max())
        res.lcd = max(0.0, 2.0 * dmax)
        if dmax > self.MARGIN - radii.max() - 0.5:
            res.notes.append("LCD above the measurable limit (%.0f Å)." % res.lcd)

        # --- probe-occupiable volume: dilation of {D >= r_probe} ---
        centers_mask = dgrid >= rp
        if centers_mask.any():
            steps = [cellm[i] / shape[i] for i in range(3)]
            reach = [int(math.ceil(rp / (st.perpendicular_widths[i] / shape[i]))) for i in range(3)]
            occupiable = centers_mask.copy()
            for i in range(-reach[0], reach[0] + 1):
                for j in range(-reach[1], reach[1] + 1):
                    for kk in range(-reach[2], reach[2] + 1):
                        if (i, j, kk) == (0, 0, 0):
                            continue
                        if np.linalg.norm(i * steps[0] + j * steps[1] + kk * steps[2]) <= rp:
                            occupiable |= np.roll(centers_mask, (i, j, kk), axis=(0, 1, 2))
            res.pore_volume_cm3_g = float(occupiable.mean()) * to_cm3_g
        else:
            res.pore_volume_cm3_g = 0.0
        res.pore_volume_label = "probe-occupiable volume (d = %.2f Å)" % self.probe_d

        # --- PLD by bisection on the percolation radius ---
        self.log("Pore limiting diameter (percolation)...")
        lo, hi = 0.0, max(dmax, 0.0)
        for _ in range(12 if hi > 0 else 0):                  # lo stays zero if nothing percolates
            mid = 0.5 * (lo + hi)
            if self._percolates(dgrid >= mid):
                lo = mid
            else:
                hi = mid
        res.pld = 2.0 * lo
        res.notes.append("PLD/LCD accuracy ~ %.2f Å (grid spacing)." % spacing)
        if res.pld < self.probe_d:
            res.notes.append("PLD < probe diameter: the windows are geometrically closed "
                             "for a rigid framework; the surface area and volume reported "
                             "include non-connected cavities.")
        return res

    def run(self, mode: str = "auto") -> GeometricResult:
        """mode = 'auto' | 'csd' | 'internal'."""
        unchanged = not any(k.startswith("guest:") for k in self.st.p1_removed)
        if mode in ("auto", "csd") and self.st.csd_cif_path and (unchanged or mode == "csd"):
            try:
                self.log("Pore Analyser (CSD Python API)...")
                return self.run_csd()
            except Exception as exc:
                self.log("Pore Analyser unavailable (%s) -> internal calculation." % exc)
        return self.run_internal()


# =============================================================================
# 6. EXPERIMENTAL DATA
# =============================================================================
@dataclass
class Isotherm:
    """Isotherm: relative pressure and amount adsorbed in cm3(STP)/g."""
    p_rel: np.ndarray
    v_cm3: np.ndarray
    label: str = ""
    des_p_rel: np.ndarray = field(default_factory=lambda: np.zeros(0))
    des_v_cm3: np.ndarray = field(default_factory=lambda: np.zeros(0))
    meta: Dict[str, object] = field(default_factory=dict)

    def to_frame(self) -> pd.DataFrame:
        frame = pd.DataFrame({"Branch": "adsorption", "P/P0": self.p_rel,
                              "Adsorption_cm3STPg": self.v_cm3})
        if len(self.des_p_rel):
            frame = pd.concat([frame, pd.DataFrame({"Branch": "desorption", "P/P0": self.des_p_rel,
                                                    "Adsorption_cm3STPg": self.des_v_cm3})],
                              ignore_index=True)
        frame["Adsorption_mmol_g"] = frame["Adsorption_cm3STPg"] / (V_MOLAR_STP / 1000.0)
        return frame


class ExperimentalDataLoader:
    """Tolerant reader for experimental isotherms.

    Accepted formats: .xlsx, .xls, .csv/.txt and pasted text. Headers are
    recognised by keywords ("P/P0", "Adsorption", "Va/cm3(STP)g-1",
    "Quantity Adsorbed"...). BELSORP/BELMaster exports (ADS / DES blocks,
    "Isotherm" sheet) and simple two-column tables are supported.
    If the amount is in mmol/g or mol/kg, it is converted to cm3(STP)/g.
    """

    RE_PRESSURE = re.compile(r"p\s*/\s*p\s*[0o°]|relative\s*pressure|pression\s*relative|p_?rel", re.I)
    RE_AMOUNT = re.compile(r"adsor|va\b|va/|quantit|uptake|loading|cm3|cm\^3|cc/g|mmol|mol/kg|volume", re.I)
    RE_EXCLUDE = re.compile(r"kpa|\bpa\b|mmhg|torr|\bbar\b|^no\.?$|p\s*/\s*va|/\s*\(?va", re.I)
    META_KEYS = {
        "p0_kPa": r"saturation vapor pressure", "adsorptive": r"^adsorptive$",
        "temperature_K": r"adsorption temperature", "cross_section_nm2": r"cross-sectional area",
        "sample_mass_g": r"sample mass", "bet_report_m2g": r"^as,\s*bet", "vm_report_cm3g": r"^vm$",
        "c_report": r"^c$", "pore_volume_report_cm3g": r"total pore volume", "comment": r"^comment$",
    }

    # --------------------------------------------------------------------- files
    @classmethod
    def load_file(cls, path: str) -> Isotherm:
        ext = os.path.splitext(path)[1].lower()
        if ext in (".xlsx", ".xlsm", ".xls"):
            try:
                sheets = pd.read_excel(path, sheet_name=None, header=None)
            except ImportError as exc:
                need = "xlrd" if ext == ".xls" else "openpyxl"
                raise RuntimeError("Cannot read %s files: install the '%s' module "
                                   "(python -m pip install %s)." % (ext, need, need)) from exc
            meta: Dict[str, object] = {}
            for frame in sheets.values():
                for key, value in cls._metadata(frame).items():
                    meta.setdefault(key, value)
            ordered = sorted(sheets.items(), key=lambda kv: 0 if "isotherm" in str(kv[0]).lower() else 1)
            errors = []
            for name, frame in ordered:
                try:
                    iso = cls._from_table(frame)
                    iso.meta.update(meta)
                    iso.meta["sheet"] = str(name)
                    iso.label = os.path.basename(path)
                    return iso
                except ValueError as exc:
                    errors.append("%s: %s" % (name, exc))
            raise ValueError("No isotherm found in %s.\n%s" % (path, "\n".join(errors)))
        with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
            iso = cls.parse_text(handle.read())
        iso.label = os.path.basename(path)
        return iso

    # ---------------------------------------------------------------- pasted text
    @classmethod
    def parse_text(cls, text: str) -> Isotherm:
        """Converts a pasted table (tabs, spaces, ';' or ',')."""
        rows = []
        for line in text.replace("\r", "\n").split("\n"):
            line = line.strip()
            if not line:
                continue
            if "\t" in line:
                cells = line.split("\t")
            elif ";" in line:
                cells = line.split(";")
            elif "," in line and not re.search(r"\d,\d+(\s|$)", line.replace(", ", " ")):
                cells = line.split(",")                      # English-style CSV
            elif "," in line and "." in line:
                cells = line.split(",")
            else:
                cells = line.split()                         # spaces (comma = decimal separator)
            rows.append([c.strip() for c in cells])
        if not rows:
            raise ValueError("No data to read.")
        width = max(len(r) for r in rows)
        frame = pd.DataFrame([r + [None] * (width - len(r)) for r in rows])
        iso = cls._from_table(frame)
        iso.label = "pasted data"
        return iso

    # ----------------------------------------------------------------- parser core
    @classmethod
    def _metadata(cls, frame: pd.DataFrame) -> Dict[str, object]:
        meta: Dict[str, object] = {}
        for row in frame.itertuples(index=False):
            cells = list(row)
            for pos, cell in enumerate(cells):
                if not isinstance(cell, str):
                    continue
                for key, pattern in cls.META_KEYS.items():
                    if key in meta or not re.search(pattern, cell.strip(), re.I):
                        continue
                    for other in cells[pos + 1:pos + 4]:
                        if other is None or (isinstance(other, float) and np.isnan(other)):
                            continue
                        if key in ("adsorptive", "comment"):
                            meta[key] = str(other)
                        elif np.isfinite(to_float(other)):
                            meta[key] = to_float(other)
                        break
        return meta

    @classmethod
    def _from_table(cls, frame: pd.DataFrame) -> Isotherm:
        values = frame.values.tolist()
        # 1) look for the header row
        header_row, p_col, v_col, factor = None, None, None, 1.0
        for r, row in enumerate(values):
            texts = [(c, str(cell)) for c, cell in enumerate(row) if isinstance(cell, str)]
            p_cands = [c for c, t in texts if cls.RE_PRESSURE.search(t) and not cls.RE_EXCLUDE.search(t)]
            if not p_cands:
                continue
            v_cands = [(c, t) for c, t in texts
                       if c != p_cands[0] and cls.RE_AMOUNT.search(t) and not cls.RE_EXCLUDE.search(t)
                       and not cls.RE_PRESSURE.search(t)]
            if v_cands:
                header_row, p_col = r, p_cands[0]
                v_col, v_text = v_cands[0]
                if re.search(r"mmol|mol\s*/\s*kg|mol\s*kg", v_text, re.I):
                    factor = V_MOLAR_STP / 1000.0            # mmol/g -> cm3(STP)/g
                break
        numeric = [[to_float(cell) if not isinstance(cell, str) or re.match(r"^\s*[-+0-9.,eE]+\s*$", cell)
                    else float("nan") for cell in row] for row in values]
        # 2) no header: first two numeric columns, or successive pairs
        def headerless(rows):
            data_rows = [[x for x in row if np.isfinite(x)] for row in rows]
            data_rows = [row for row in data_rows if len(row) >= 2]
            if not data_rows:
                raise ValueError("P/P0 and adsorption columns not found")
            if len(data_rows) == 1 and len(data_rows[0]) >= 4 and len(data_rows[0]) % 2 == 0:
                flat = data_rows[0]                          # "0.0001 15.3 0.001 31.3 ..."
                return [(flat[i], flat[i + 1] * factor) for i in range(0, len(flat), 2)]
            return [(row[0], row[1] * factor) for row in data_rows]

        if header_row is None:
            ads, des = headerless(numeric), []
        else:
            ads, des, branch = [], [], "ads"
            for r in range(header_row + 1, len(values)):
                marks = [str(c).strip().upper() for c in values[r] if isinstance(c, str)]
                if any(m.startswith("DES") for m in marks):
                    branch = "des"
                elif any(m.startswith("ADS") for m in marks):
                    branch = "ads"
                p = numeric[r][p_col] if p_col < len(numeric[r]) else float("nan")
                v = numeric[r][v_col] if v_col < len(numeric[r]) else float("nan")
                if np.isfinite(p) and np.isfinite(v):
                    (ads if branch == "ads" else des).append((p, v * factor))
            if len(ads) < 3:                                 # header, then all values on a single line
                ads, des = headerless(numeric[header_row + 1:]), []
        ads = [(p, v) for p, v in ads if p > 0]
        if len(ads) < 3:
            raise ValueError("fewer than 3 usable adsorption points")
        p_arr = np.array([p for p, _ in ads])
        if p_arr.max() > 1.5:
            raise ValueError("the pressure column is not a relative pressure P/P0 "
                             "(maximum value %.3g)" % p_arr.max())
        return Isotherm(p_arr, np.array([v for _, v in ads]),
                        des_p_rel=np.array([p for p, _ in des]), des_v_cm3=np.array([v for _, v in des]))


# =============================================================================
# 7. BET ANALYSIS (Rouquerol criteria) AND PORE VOLUME (Gurvich)
# =============================================================================
@dataclass
class BETResult:
    """Result of a BET fit."""
    valid: bool = False
    message: str = ""
    p_min: float = float("nan")
    p_max: float = float("nan")
    n_points: int = 0
    slope: float = float("nan")
    intercept: float = float("nan")
    r2: float = float("nan")
    vm_cm3_g: float = float("nan")            # monolayer capacity
    c_constant: float = float("nan")
    surface_m2_g: float = float("nan")
    p_at_vm: float = float("nan")             # P/P0 where n = n_m (criterion 3)
    criterion4_error: float = float("nan")    # |1/(sqrt(C)+1) - p(n_m)| / p(n_m)
    rouquerol_ok: Tuple[bool, bool, bool, bool] = (False, False, False, False)
    plot_p: np.ndarray = field(default_factory=lambda: np.zeros(0))       # all points
    plot_y: np.ndarray = field(default_factory=lambda: np.zeros(0))       # p/[v(1-p)]
    in_range: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))

    def summary(self) -> Dict[str, object]:
        return {"P/P0 range min": self.p_min, "P/P0 range max": self.p_max,
                "Number of points": self.n_points, "Slope": self.slope,
                "Intercept": self.intercept, "R2": self.r2,
                "Vm (cm3 STP/g)": self.vm_cm3_g, "C constant": self.c_constant,
                "BET surface area (m2/g)": self.surface_m2_g, "P/P0 at Vm": self.p_at_vm,
                "Criterion 4 deviation (%)": 100.0 * self.criterion4_error,
                "Rouquerol 1 (n(1-p) increasing)": self.rouquerol_ok[0],
                "Rouquerol 2 (C > 0)": self.rouquerol_ok[1],
                "Rouquerol 3 (p(Vm) within range)": self.rouquerol_ok[2],
                "Rouquerol 4 (1/(sqrt(C)+1) ~ p(Vm))": self.rouquerol_ok[3],
                "Remark": self.message}


class BETAnalyzer:
    """Classical BET method with automatic range selection.

    Linearised equation:  p / [v (1 - p)] = 1/(vm C) + (C - 1)/(vm C) * p ,
    with p = P/P0. Every contiguous window of at least `min_points` points is
    fitted, then filtered with the Rouquerol criteria
    (Rouquerol, Llewellyn, Rouquerol, Stud. Surf. Sci. Catal. 2007, 160, 49):
      1. v(1 - p) increases monotonically with p over the range;
      2. C > 0 (positive intercept);
      3. the pressure p(vm) at which the monolayer is reached lies in the range;
      4. 1/(sqrt(C) + 1) deviates from p(vm) by no more than 20 %.
    Among the compliant windows (R2 >= 0.995), the one with the most points is
    kept, then the best R2 in case of a tie. If no window satisfies everything,
    the requirements are relaxed and the result says so.
    For a simulated isotherm, which is noisy by nature, `noise_tolerance`
    allows small decreases of v(1 - p) (of the order of the statistical error).
    """

    def __init__(self, cross_section_nm2: float = 0.162, min_points: int = 3,
                 max_p_rel: float = 0.5, min_r2: float = 0.995, noise_tolerance: float = 0.0):
        self.cross_section = cross_section_nm2
        # Relative decrease of v(1-p) tolerated between two points for criterion 1:
        # 0 for a measurement; ~2 statistical uncertainties for a simulation.
        self.noise_tolerance = noise_tolerance
        self.min_points = min_points
        self.max_p_rel = max_p_rel
        self.min_r2 = min_r2

    def surface_from_vm(self, vm_cm3_g: float) -> float:
        """Vm [cm3 STP/g] -> surface [m2/g]."""
        return vm_cm3_g / V_MOLAR_STP * N_AVOGADRO * self.cross_section * 1e-18

    def analyse(self, p_rel: Sequence[float], v_cm3: Sequence[float]) -> BETResult:
        p = np.asarray(p_rel, dtype=float)
        v = np.asarray(v_cm3, dtype=float)
        ok = np.isfinite(p) & np.isfinite(v) & (p > 0) & (p < 1) & (v > 0)
        p, v = p[ok], v[ok]
        order = np.argsort(p)
        p, v = p[order], v[order]
        result = BETResult(plot_p=p, plot_y=p / (v * (1.0 - p)) if len(p) else np.zeros(0),
                           in_range=np.zeros(len(p), dtype=bool))
        if len(p) < self.min_points:
            result.message = "Not enough points for a BET analysis."
            return result
        y = p / (v * (1.0 - p))
        rouq = v * (1.0 - p)                                  # Rouquerol function
        log_p = np.log(p)
        candidates = []
        n = len(p)
        for i in range(n - self.min_points + 1):
            for j in range(i + self.min_points - 1, n):
                if p[j] > self.max_p_rel:
                    break
                xs, ys = p[i:j + 1], y[i:j + 1]
                slope, intercept = np.polyfit(xs, ys, 1)
                fit = slope * xs + intercept
                ss_tot = float(np.sum((ys - ys.mean()) ** 2))
                r2 = 1.0 - float(np.sum((ys - fit) ** 2)) / ss_tot if ss_tot > 0 else 0.0
                if slope <= 0:
                    continue
                c1 = bool(np.all(np.diff(rouq[i:j + 1]) >= -self.noise_tolerance * rouq[i:j]))
                c2 = bool(intercept > 0)
                vm = 1.0 / (slope + intercept)
                c_const = slope / intercept + 1.0 if intercept != 0 else float("inf")
                # p(vm): interpolation (linear in ln p) of the isotherm at v = vm
                p_vm = float("nan")
                above = np.nonzero(v >= vm)[0]
                if len(above) and above[0] > 0:
                    k = above[0]
                    t = (vm - v[k - 1]) / (v[k] - v[k - 1]) if v[k] != v[k - 1] else 0.0
                    p_vm = float(np.exp(log_p[k - 1] + t * (log_p[k] - log_p[k - 1])))
                elif len(above):
                    p_vm = float(p[0])
                c3 = bool(np.isfinite(p_vm) and p[i] <= p_vm <= p[j])
                err4 = float("nan")
                if c2 and np.isfinite(p_vm) and p_vm > 0:
                    err4 = abs(1.0 / (math.sqrt(c_const) + 1.0) - p_vm) / p_vm
                c4 = bool(np.isfinite(err4) and err4 <= 0.20)
                candidates.append(dict(i=i, j=j, slope=slope, intercept=intercept, r2=r2, vm=vm,
                                       c=c_const, p_vm=p_vm, err4=err4, flags=(c1, c2, c3, c4)))
        if not candidates:
            result.message = "No usable BET window (negative slope everywhere)."
            return result

        def pick(pool):
            return max(pool, key=lambda d: (d["j"] - d["i"] + 1, d["r2"]))

        tiers = [
            ("", [d for d in candidates if all(d["flags"]) and d["r2"] >= self.min_r2]),
            ("Rouquerol criterion 4 not met (deviation > 20 %).",
             [d for d in candidates if all(d["flags"][:3]) and d["r2"] >= self.min_r2]),
            ("Rouquerol criterion 3 not met: p(Vm) outside the fitted range.",
             [d for d in candidates if all(d["flags"][:2]) and d["r2"] >= self.min_r2]),
            ("Rouquerol criteria not met: indicative value only.",
             [d for d in candidates if d["flags"][1]]),
            ("Negative C constant: the BET method does not apply to these data.", candidates),
        ]
        for message, pool in tiers:
            if pool:
                best = pick(pool)
                result.message = message
                break
        i, j = best["i"], best["j"]
        result.valid = True
        result.p_min, result.p_max, result.n_points = float(p[i]), float(p[j]), j - i + 1
        result.slope, result.intercept, result.r2 = float(best["slope"]), float(best["intercept"]), float(best["r2"])
        result.vm_cm3_g, result.c_constant = float(best["vm"]), float(best["c"])
        result.surface_m2_g = self.surface_from_vm(result.vm_cm3_g)
        result.p_at_vm, result.criterion4_error = best["p_vm"], best["err4"]
        result.rouquerol_ok = best["flags"]
        result.in_range[i:j + 1] = True
        return result


def gurvich_pore_volume(p_rel: Sequence[float], v_cm3: Sequence[float], p_target: float,
                        liquid_molar_volume: float) -> Tuple[float, float]:
    """Pore volume (cm3/g) from the Gurvich rule at P/P0 = p_target.

    The adsorbed amount is assumed to be in the liquid state:
        Vp = v_ads / 22414 * Vm_liquid.
    If p_target is outside the measured range, the nearest point is used.
    Returns (volume, P/P0 actually used).
    """
    p = np.asarray(p_rel, dtype=float)
    v = np.asarray(v_cm3, dtype=float)
    ok = np.isfinite(p) & np.isfinite(v)
    p, v = p[ok], v[ok]
    if len(p) == 0:
        return float("nan"), float("nan")
    order = np.argsort(p)
    p, v = p[order], v[order]
    used = float(np.clip(p_target, p[0], p[-1]))
    amount = float(np.interp(used, p, v))
    return amount / V_MOLAR_STP * liquid_molar_volume, used


# =============================================================================
# 8. RASPA3 PREPARATION AND EXECUTION
# =============================================================================
@dataclass
class SimulationSettings:
    """All parameters that can be set from the GUI."""
    gas: str = "N2"
    temperature: float = 77.0                 # K
    p0: float = 97190.0                       # Pa: saturation pressure at T
    n_init_cycles: int = 2000                 # "equilibration" (RASPA3 Initialization)
    n_prod_cycles: int = 5000
    supercell: Tuple[int, int, int] = (1, 1, 1)
    forcefield: str = FORCE_FIELDS[0]
    charge_mode: str = CHARGE_MODES[0]
    cutoff: float = 12.0                      # Å
    p_rel_min: float = 1e-6
    p_rel_max: float = 1.0
    n_points: int = 35
    cross_section: float = 0.162              # nm2
    gurvich_p_rel: float = 0.95
    raspa_path: str = ""
    n_jobs: int = 1                           # simultaneous simulations
    work_dir: str = ""

    def pressure_grid(self) -> np.ndarray:
        """Relative pressures, logarithmically spaced."""
        return np.logspace(math.log10(self.p_rel_min), math.log10(self.p_rel_max), int(self.n_points))


def find_raspa_executable() -> str:
    """Looks for the raspa3 executable in the usual locations (macOS, Linux, conda)."""
    candidates = [shutil.which("raspa3"), "/usr/local/bin/raspa3", "/opt/homebrew/bin/raspa3",
                  "/usr/bin/raspa3", "/opt/local/bin/raspa3",
                  os.path.join(sys.prefix, "bin", "raspa3"),
                  os.path.expanduser("~/miniconda3/bin/raspa3"),
                  os.path.expanduser("~/miniforge3/bin/raspa3"),
                  os.path.expanduser("~/anaconda3/bin/raspa3")]
    for cand in candidates:
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return ""


class RaspaInputBuilder:
    """Writes the RASPA3 input files (version 3 JSON format)."""

    FRAMEWORK_NAME = "framework"

    def __init__(self, structure: StructureData, settings: SimulationSettings):
        self.st = structure
        self.s = settings
        self.gas = ADSORBATES[settings.gas]

    @property
    def adsorbate_charged(self) -> bool:
        return any(abs(pa["charge"]) > 0 for pa in self.gas["pseudo_atoms"])

    def validate(self) -> List[str]:
        """Returns the list of blocking problems (empty if everything is fine)."""
        problems = []
        if not self.st.p1_elements:
            problems.append("The P1 structure is empty.")
        for el in sorted(set(self.st.p1_elements)):
            try:
                framework_lj(el, self.s.forcefield)
            except ValueError as exc:
                problems.append(str(exc))
        if self.s.charge_mode == CHARGE_MODES[1] and self.st.p1_charges is None:
            problems.append("The CIF contains no charges (_atom_site_charge): "
                            "choose another charge mode.")
        if not (0 < self.s.p_rel_min < self.s.p_rel_max <= 1.5):
            problems.append("Invalid relative pressure range.")
        return problems

    def force_field(self) -> dict:
        """force_field.json: framework pseudo-atoms (one per element) + adsorbate."""
        pseudo, inter = [], []
        for el in sorted(set(self.st.p1_elements)):
            eps, sig, source = framework_lj(el, self.s.forcefield)
            pseudo.append({"name": el, "framework": True, "print_to_output": True, "element": el,
                           "print_as": el, "mass": ATOMIC_MASS[el], "charge": 0.0})
            inter.append({"name": el, "type": "lennard-jones",
                          "parameters": [round(eps, 5), round(sig, 5)], "source": source})
        for pa in self.gas["pseudo_atoms"]:
            pseudo.append({"name": pa["name"], "framework": False, "print_to_output": True,
                           "element": pa["element"], "print_as": pa["element"],
                           "mass": pa["mass"], "charge": pa["charge"]})
            inter.append({"name": pa["name"],
                          "type": "lennard-jones" if pa["epsilon"] > 0 else "none",
                          "parameters": [pa["epsilon"], pa["sigma"]], "source": self.gas["model"]})
        return {"MixingRule": "Lorentz-Berthelot", "TruncationMethod": "shifted",
                "TailCorrections": False, "CutOffVDW": self.s.cutoff, "CutOffCoulomb": self.s.cutoff,
                "PseudoAtoms": pseudo, "SelfInteractions": inter}

    def molecule(self) -> dict:
        """<gas>.json: critical constants (Peng-Robinson fugacity) + rigid geometry."""
        return {"CriticalTemperature": self.gas["Tc"], "CriticalPressure": self.gas["Pc"],
                "AcentricFactor": self.gas["omega"], "pseudoAtoms": self.gas["sites"]}

    def simulation(self, pressure_pa: float) -> dict:
        """simulation.json for one pressure point (GCMC, muVT ensemble)."""
        charged = self.adsorbate_charged
        system = {"Type": "Framework", "Name": self.FRAMEWORK_NAME,
                  "NumberOfUnitCells": list(self.s.supercell),
                  "ChargeMethod": "Ewald" if charged else "None",
                  "ExternalTemperature": self.s.temperature,
                  "ExternalPressure": float(pressure_pa)}
        if charged:
            system["UseChargesFrom"] = {CHARGE_MODES[0]: "ChargeEquilibration",
                                        CHARGE_MODES[1]: "CIF_File",
                                        CHARGE_MODES[2]: "PseudoAtoms"}[self.s.charge_mode]
        component = {"Name": self.s.gas, "TranslationProbability": 0.5,
                     "ReinsertionProbability": 0.5, "SwapProbability": 1.0,
                     "CreateNumberOfMolecules": 0}
        if len(self.gas["sites"]) > 1:
            component["RotationProbability"] = 0.5
        # "FugacityCoefficient" deliberately omitted: RASPA3 computes it with
        # the Peng-Robinson equation of state (critical constants above).
        total = self.s.n_init_cycles + self.s.n_prod_cycles
        return {"SimulationType": "MonteCarlo",
                "NumberOfCycles": int(self.s.n_prod_cycles),
                "NumberOfInitializationCycles": int(self.s.n_init_cycles),
                "NumberOfEquilibrationCycles": 0,
                "PrintEvery": int(max(10, total // 50)),
                "Systems": [system], "Components": [component]}

    def write_point(self, directory: str, pressure_pa: float) -> None:
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, "simulation.json"), "w", encoding="utf-8") as handle:
            json.dump(self.simulation(pressure_pa), handle, indent=2)
        with open(os.path.join(directory, "force_field.json"), "w", encoding="utf-8") as handle:
            json.dump(self.force_field(), handle, indent=2)
        with open(os.path.join(directory, self.s.gas + ".json"), "w", encoding="utf-8") as handle:
            json.dump(self.molecule(), handle, indent=2)
        self.st.write_p1_cif(os.path.join(directory, self.FRAMEWORK_NAME + ".cif"),
                             with_charges=self.s.charge_mode == CHARGE_MODES[1])


class RaspaError(RuntimeError):
    """RASPA3 run error (message meant for the user)."""


class RaspaRunner:
    """Runs RASPA3 for each pressure, tracks progress, reads the results.

    One sub-folder per point (`point_01_...`): the simulations are independent
    and can run in parallel (`n_jobs`). RASPA3 standard output and standard
    error are kept in each folder (raspa_stdout.txt / raspa_stderr.txt) for
    troubleshooting.
    """

    RE_LOADING = re.compile(r"Abs\. loading average\s+([-+0-9.eE]+|nan|inf)\s+\+/-\s+"
                            r"([-+0-9.eE]+|nan|inf)\s+\[mol/kg(?:-framework)?\]", re.I)
    RE_UC = re.compile(r"Abs\. loading average\s+([-+0-9.eE]+)\s+\+/-\s+([-+0-9.eE]+)\s+"
                       r"\[molecules/(?:uc|unit cell)\]", re.I)
    RE_CYCLE = re.compile(r"(Initialization: )?Current cycle:\s*(\d+) out of (\d+)")

    def __init__(self, structure: StructureData, settings: SimulationSettings,
                 log: Optional[Callable[[str], None]] = None,
                 progress: Optional[Callable[[float, str], None]] = None):
        self.st = structure
        self.s = settings
        self.log = log or (lambda msg: None)
        self.progress = progress or (lambda fraction, text: None)
        self._cancel = threading.Event()

    def cancel(self) -> None:
        """Stop request (thread-safe): the RASPA3 processes are terminated."""
        self._cancel.set()

    # ------------------------------------------------------------------- helpers
    def _environment(self) -> dict:
        env = dict(os.environ)
        if not env.get("RASPA_DIR"):
            prefix = os.path.dirname(os.path.dirname(os.path.realpath(self.s.raspa_path)))
            for cand in (os.path.join(prefix, "share", "raspa3"), "/usr/local/share/raspa3",
                         "/usr/share/raspa3", "/opt/homebrew/share/raspa3"):
                if os.path.isdir(cand):
                    env["RASPA_DIR"] = cand
                    break
        return env

    @staticmethod
    def _output_text(directory: str) -> str:
        files = sorted(glob.glob(os.path.join(directory, "output", "*.txt")))
        if not files:
            return ""
        try:
            with open(files[-1], "r", encoding="utf-8", errors="replace") as handle:
                return handle.read()
        except OSError:
            return ""

    def _point_fraction(self, directory: str) -> float:
        """Progress (0-1) of a running simulation, read from its output file."""
        matches = self.RE_CYCLE.findall(self._output_text(directory)[-20000:])
        if not matches:
            return 0.0
        init, cur, _total = matches[-1]
        done = int(cur) + (0 if init else self.s.n_init_cycles)
        return min(1.0, done / max(1, self.s.n_init_cycles + self.s.n_prod_cycles))

    def _parse_point(self, directory: str) -> Tuple[float, float, float]:
        """(loading mol/kg, uncertainty mol/kg, molecules per cell) of a finished point."""
        text = self._output_text(directory)
        found = self.RE_LOADING.findall(text)
        if not found:
            raise RaspaError("loading result not found in the output file")
        loading, error = float(found[-1][0]), float(found[-1][1])
        per_uc = self.RE_UC.findall(text)
        return loading, error, float(per_uc[-1][0]) if per_uc else float("nan")

    @staticmethod
    def _tail(path: str, n_chars: int = 1500) -> str:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                return handle.read()[-n_chars:].strip()
        except OSError:
            return ""

    # ------------------------------------------------------------------- execution
    def run(self) -> pd.DataFrame:
        s = self.s
        if not s.raspa_path or not os.path.isfile(s.raspa_path):
            raise RaspaError("RASPA3 executable not found: '%s'.\nSet its path in "
                             "the \"RASPA Settings\" tab." % s.raspa_path)
        if not os.access(s.raspa_path, os.X_OK):
            raise RaspaError("File '%s' is not executable (chmod +x)." % s.raspa_path)
        builder = RaspaInputBuilder(self.st, s)
        problems = builder.validate()
        if problems:
            raise RaspaError("Invalid parameters:\n - " + "\n - ".join(problems))
        os.makedirs(s.work_dir, exist_ok=True)
        grid = s.pressure_grid()
        points = []
        for k, p_rel in enumerate(grid):
            directory = os.path.join(s.work_dir, "point_%02d_%.3e" % (k + 1, p_rel))
            builder.write_point(directory, p_rel * s.p0)
            points.append({"index": k, "p_rel": float(p_rel), "dir": directory,
                           "state": "pending", "loading": float("nan"),
                           "error": float("nan"), "per_uc": float("nan"), "message": ""})
        n_atoms = len(self.st.p1_elements) * int(np.prod(s.supercell))
        self.log("Input files written to: %s" % s.work_dir)
        self.log("%d points, %s at %.2f K, P0 = %.0f Pa, supercell %s (%d atoms), "
                 "%d + %d cycles, %d simultaneous simulation(s)."
                 % (len(points), s.gas, s.temperature, s.p0, "x".join(map(str, s.supercell)),
                    n_atoms, s.n_init_cycles, s.n_prod_cycles, s.n_jobs))
        env = self._environment()
        pending = deque(points)
        running: List[Tuple[dict, subprocess.Popen, object, object]] = []
        finished = 0
        t0 = time.time()
        try:
            while pending or running:
                if self._cancel.is_set():
                    raise RaspaError("Calculation stopped by the user.")
                while pending and len(running) < max(1, s.n_jobs):
                    point = pending.popleft()
                    out = open(os.path.join(point["dir"], "raspa_stdout.txt"), "w")
                    err = open(os.path.join(point["dir"], "raspa_stderr.txt"), "w")
                    try:
                        proc = subprocess.Popen([s.raspa_path], cwd=point["dir"], stdout=out,
                                                stderr=err, env=env)
                    except OSError as exc:
                        out.close()
                        err.close()
                        raise RaspaError("Cannot start RASPA3: %s" % exc) from exc
                    point["state"] = "running"
                    running.append((point, proc, out, err))
                time.sleep(0.4)
                still = []
                for point, proc, out, err in running:
                    code = proc.poll()
                    if code is None:
                        still.append((point, proc, out, err))
                        continue
                    out.close()
                    err.close()
                    finished += 1
                    try:
                        if code != 0:
                            raise RaspaError("exit code %d" % code)
                        point["loading"], point["error"], point["per_uc"] = self._parse_point(point["dir"])
                        point["state"] = "done"
                        self.log("P/P0 = %.3e: %.4f ± %.4f mol/kg  (%.1f cm3 STP/g)"
                                 % (point["p_rel"], point["loading"], point["error"],
                                    point["loading"] * V_MOLAR_STP / 1000.0))
                    except RaspaError as exc:
                        point["state"] = "failed"
                        details = self._tail(os.path.join(point["dir"], "raspa_stderr.txt")) or \
                            self._tail(os.path.join(point["dir"], "raspa_stdout.txt")) or \
                            self._output_text(point["dir"])[-1500:]
                        point["message"] = "%s\n%s" % (exc, details)
                        self.log("RASPA3 ERROR at P/P0 = %.3e: %s" % (point["p_rel"], point["message"]))
                running = still
                partial = sum(self._point_fraction(pt["dir"]) for pt, _, _, _ in running)
                fraction = (finished + partial) / len(points)
                elapsed = time.time() - t0
                eta = elapsed * (1 - fraction) / fraction if fraction > 0.02 else float("nan")
                self.progress(fraction, "%d/%d points done — elapsed %s%s"
                              % (finished, len(points), _dt.timedelta(seconds=int(elapsed)),
                                 "" if not np.isfinite(eta) else ", remaining ~%s" % _dt.timedelta(seconds=int(eta))))
        finally:
            for point, proc, out, err in running:                # clean shutdown if interrupted
                if proc.poll() is None:
                    proc.terminate()
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                out.close()
                err.close()
        failed = [pt for pt in points if pt["state"] != "done"]
        if len(failed) == len(points):
            raise RaspaError("All RASPA3 simulations failed.\nFirst message:\n%s"
                             % failed[0]["message"])
        if failed:
            self.log("WARNING: %d failed point(s), ignored in the analysis." % len(failed))
        frame = pd.DataFrame({
            "P/P0": [pt["p_rel"] for pt in points],
            "Pressure_Pa": [pt["p_rel"] * s.p0 for pt in points],
            "Adsorption_mmol_g": [pt["loading"] for pt in points],            # mol/kg = mmol/g
            "Uncertainty_mmol_g": [pt["error"] for pt in points],
            "Adsorption_cm3STPg": [pt["loading"] * V_MOLAR_STP / 1000.0 for pt in points],
            "Uncertainty_cm3STPg": [pt["error"] * V_MOLAR_STP / 1000.0 for pt in points],
            "Molecules_per_cell": [pt["per_uc"] for pt in points],
            "Status": [pt["state"] for pt in points],
        })
        return frame


# =============================================================================
# 9. FULL ANALYSIS AND EXPERIMENT / SIMULATION COMPARISON
# =============================================================================
@dataclass
class AnalysisResults:
    """Container for all the results of a session."""
    structure: StructureData
    settings: SimulationSettings
    experimental: Optional[Isotherm] = None
    simulated: Optional[pd.DataFrame] = None
    bet_exp: Optional[BETResult] = None
    bet_sim: Optional[BETResult] = None
    geometric: Optional[GeometricResult] = None
    vpore_exp: float = float("nan")
    vpore_exp_p: float = float("nan")
    vpore_sim: float = float("nan")
    vpore_sim_p: float = float("nan")
    created: str = field(default_factory=lambda: _dt.datetime.now().strftime("%Y-%m-%d %H:%M"))

    # -------------------------------------------------------------- calculations
    def analyse(self) -> None:
        """(Re)computes BET and pore volumes from the available isotherms."""
        gas = ADSORBATES[self.settings.gas]
        bet = BETAnalyzer(self.settings.cross_section)
        if self.experimental is not None:
            self.bet_exp = bet.analyse(self.experimental.p_rel, self.experimental.v_cm3)
            self.vpore_exp, self.vpore_exp_p = gurvich_pore_volume(
                self.experimental.p_rel, self.experimental.v_cm3,
                self.settings.gurvich_p_rel, gas["liquid_molar_volume"])
        if self.simulated is not None:
            ok = self.simulated["Status"] == "done"
            p = self.simulated.loc[ok, "P/P0"].values
            v = self.simulated.loc[ok, "Adsorption_cm3STPg"].values
            err = self.simulated.loc[ok, "Uncertainty_cm3STPg"].values
            relative = err[v > 0] / v[v > 0]
            noise = float(np.clip(2.0 * np.nanmedian(relative), 0.0, 0.05)) if len(relative) else 0.0
            self.bet_sim = BETAnalyzer(self.settings.cross_section, noise_tolerance=noise).analyse(p, v)
            self.vpore_sim, self.vpore_sim_p = gurvich_pore_volume(
                p, v, self.settings.gurvich_p_rel, gas["liquid_molar_volume"])

    # -------------------------------------------------------------------- tables
    def comparison_frame(self) -> pd.DataFrame:
        """Comparison table: experiment / simulation / geometry."""
        be, bs, geo = self.bet_exp, self.bet_sim, self.geometric
        s_exp = be.surface_m2_g if be and be.valid else float("nan")
        s_sim = bs.surface_m2_g if bs and bs.valid else float("nan")
        s_geo = geo.surface_m2_g if geo else float("nan")
        rows = [
            ("BET surface area (m2/g)", s_exp, s_sim, float("nan"), pct_diff(s_sim, s_exp)),
            ("Monolayer Vm (cm3 STP/g)", be.vm_cm3_g if be else float("nan"),
             bs.vm_cm3_g if bs else float("nan"), float("nan"),
             pct_diff(bs.vm_cm3_g if bs else float("nan"), be.vm_cm3_g if be else float("nan"))),
            ("C constant", be.c_constant if be else float("nan"),
             bs.c_constant if bs else float("nan"), float("nan"), float("nan")),
            ("Gurvich pore volume (cm3/g)", self.vpore_exp, self.vpore_sim, float("nan"),
             pct_diff(self.vpore_sim, self.vpore_exp)),
            ("P/P0 used for Gurvich", self.vpore_exp_p, self.vpore_sim_p, float("nan"), float("nan")),
            ("Geometric pore volume (cm3/g)", float("nan"), float("nan"),
             geo.pore_volume_cm3_g if geo else float("nan"), float("nan")),
            ("Accessible geometric surface area (m2/g)", float("nan"), float("nan"), s_geo, float("nan")),
            ("Experimental BET / geometric surface area (%)", float("nan"), float("nan"),
             100.0 * s_exp / s_geo if s_geo and np.isfinite(s_geo) and s_geo > 0 else float("nan"), float("nan")),
            ("Simulated BET / geometric surface area (%)", float("nan"), float("nan"),
             100.0 * s_sim / s_geo if s_geo and np.isfinite(s_geo) and s_geo > 0 else float("nan"), float("nan")),
            ("Pore limiting diameter PLD (Å)", float("nan"), float("nan"),
             geo.pld if geo else float("nan"), float("nan")),
            ("Largest cavity diameter LCD (Å)", float("nan"), float("nan"),
             geo.lcd if geo else float("nan"), float("nan")),
        ]
        frame = pd.DataFrame(rows, columns=["Quantity", "Experimental", "RASPA3 simulation",
                                            "Geometric", "Difference vs experiment (%)"])
        # Geometric vs experimental pore-volume difference, placed on the geometric row
        if geo and np.isfinite(self.vpore_exp):
            frame.loc[5, "Difference vs experiment (%)"] = pct_diff(geo.pore_volume_cm3_g, self.vpore_exp)
        return frame

    def bet_frame(self) -> pd.DataFrame:
        empty = BETResult().summary()
        exp = self.bet_exp.summary() if self.bet_exp else empty
        sim = self.bet_sim.summary() if self.bet_sim else empty
        return pd.DataFrame({"Parameter": list(exp.keys()), "Experimental": list(exp.values()),
                             "RASPA3 simulation": [sim[k] for k in exp]})

    def bet_points_frame(self) -> pd.DataFrame:
        frames = []
        for name, bet in (("experimental", self.bet_exp), ("simulation", self.bet_sim)):
            if bet is not None and len(bet.plot_p):
                frames.append(pd.DataFrame({"Series": name, "P/P0": bet.plot_p,
                                            "p/[v(1-p)] (g/cm3 STP)": bet.plot_y,
                                            "In BET range": bet.in_range}))
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def structure_frame(self) -> pd.DataFrame:
        st = self.structure
        raw = st.asymmetric_summary()
        raw_mass = st.mass_of(raw)
        a, b, c, al, be, ga = st.cell
        rows = [
            ("Name", st.name), ("File", st.source_file or "(Mercury structure)"),
            ("Reported formula", st.formula_reported or "—"),
            ("Formula weight (g/mol)", fmt(st.formula_weight_reported, 2)),
            ("Z", fmt(st.z_value, 0)),
            ("Space group", st.spacegroup),
            ("a, b, c (Å)", "%.4f, %.4f, %.4f" % (a, b, c)),
            ("alpha, beta, gamma (°)", "%.3f, %.3f, %.3f" % (al, be, ga)),
            ("Cell volume (Å3)", fmt(st.volume, 1)),
            ("Cell contents (occupancy-weighted)", hill_formula(raw)),
            ("Cell mass (g/mol)", fmt(raw_mass, 2)),
            ("Calculated density, structure as read (g/cm3)", fmt(st.density_from_mass(raw_mass), 4)),
            ("Density reported in the file (g/cm3)", fmt(st.density_reported, 4)),
            ("P1 simulation model", "%d atoms: %s" % (len(st.p1_elements), hill_formula(st.p1_counts()))),
            ("P1 model mass (g/mol)", fmt(st.p1_mass, 2)),
            ("P1 model density (g/cm3)", fmt(st.density_from_mass(st.p1_mass), 4)),
            ("Atoms discarded from the P1 model", st.removed_summary()),
        ]
        return pd.DataFrame(rows, columns=["Property", "Value"])

    def settings_frame(self) -> pd.DataFrame:
        s = self.settings
        gas = ADSORBATES[s.gas]
        rows = [("Adsorbate", gas["label"]), ("Adsorbate model", gas["model"]),
                ("Temperature (K)", s.temperature), ("P0 (Pa)", s.p0),
                ("Force field", s.forcefield), ("Framework charges", s.charge_mode),
                ("Cut-off radius (Å)", s.cutoff), ("Supercell", "x".join(map(str, s.supercell))),
                ("Equilibration cycles (Initialization)", s.n_init_cycles),
                ("Production cycles", s.n_prod_cycles),
                ("P/P0 min", s.p_rel_min), ("P/P0 max", s.p_rel_max), ("Number of points", s.n_points),
                ("BET cross-sectional area (nm2)", s.cross_section),
                ("Gurvich P/P0", s.gurvich_p_rel),
                ("Liquid molar volume (cm3/mol)", gas["liquid_molar_volume"]),
                ("RASPA3 executable", s.raspa_path), ("Working folder", s.work_dir),
                ("Date", self.created), ("Software", "%s %s" % (APP_NAME, APP_VERSION))]
        return pd.DataFrame(rows, columns=["Parameter", "Value"])


# =============================================================================
# 10. FIGURES (matplotlib without pyplot: usable in Qt and off-screen)
# =============================================================================
class FigureFactory:
    """Builds the figures from an `AnalysisResults`."""

    def __init__(self, results: AnalysisResults):
        self.r = results

    @staticmethod
    def _new(width: float, height: float):
        from matplotlib.figure import Figure
        try:
            return Figure(figsize=(width, height), dpi=100, layout="constrained")
        except TypeError:                                       # matplotlib < 3.5
            return Figure(figsize=(width, height), dpi=100, constrained_layout=True)

    @staticmethod
    def _style(ax) -> None:
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.grid(True, color="#d9d9d9", linewidth=0.6)
        ax.set_axisbelow(True)

    # Figure 1 ------------------------------------------------------------------
    def isotherm(self):
        """Experimental and simulated isotherms (linear and logarithmic scales)."""
        r = self.r
        fig = self._new(10.0, 4.4)
        axes = fig.subplots(1, 2)
        for ax, log in zip(axes, (False, True)):
            self._style(ax)
            if r.experimental is not None:
                e = r.experimental
                ax.plot(e.p_rel, e.v_cm3, "o-", color=COLOR_EXP, ms=4.5, lw=1.0,
                        label="Experimental (adsorption)")
                if len(e.des_p_rel):
                    ax.plot(e.des_p_rel, e.des_v_cm3, "o--", color=COLOR_EXP, ms=4.5, lw=0.8,
                            mfc="white", label="Experimental (desorption)")
            if r.simulated is not None:
                sim = r.simulated[r.simulated["Status"] == "done"]
                ax.errorbar(sim["P/P0"], sim["Adsorption_cm3STPg"], yerr=sim["Uncertainty_cm3STPg"],
                            fmt="s-", color=COLOR_SIM, ms=4.5, lw=1.0, capsize=2,
                            label="Simulation RASPA3 (GCMC)")
            ax.set_xlabel("P/P0" + (" (log scale)" if log else ""))
            ax.set_ylabel("Adsorption (cm3 STP/g)")
            if log:
                ax.set_xscale("log")
            else:
                ax.set_xlim(0, 1.02)
            ax.set_ylim(bottom=0)
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            axes[0].legend(handles, labels, frameon=False, fontsize=9, loc="best")
        fig.suptitle("%s adsorption isotherm at %.4g K — %s"
                     % (r.settings.gas, r.settings.temperature, r.structure.name), fontsize=11)
        return fig

    # Figures 2 and 3 ------------------------------------------------------------
    def bet(self, which: str):
        """BET plot ("exp" or "sim") with the selected range and the fitted line."""
        bet = self.r.bet_exp if which == "exp" else self.r.bet_sim
        color = COLOR_EXP if which == "exp" else COLOR_SIM
        title = "Experimental BET" if which == "exp" else "Simulated BET (RASPA3)"
        fig = self._new(5.6, 4.4)
        ax = fig.subplots()
        self._style(ax)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("P/P0")
        ax.set_ylabel("p / [v (1 - p)]  (g / cm3 STP)")
        if bet is None or not len(bet.plot_p):
            ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                    transform=ax.transAxes, color="#666666")
            return fig
        if not bet.valid:
            ax.plot(bet.plot_p, bet.plot_y, "o", color="#999999", ms=4)
            ax.text(0.5, 0.9, bet.message, ha="center", transform=ax.transAxes, fontsize=9)
            return fig
        x_max = min(max(bet.p_max * 1.8, bet.p_max + 0.02), 1.0)
        shown = bet.plot_p <= x_max
        ax.plot(bet.plot_p[shown & ~bet.in_range], bet.plot_y[shown & ~bet.in_range], "o",
                mfc="white", color="#888888", ms=5, label="Outside range")
        ax.plot(bet.plot_p[bet.in_range], bet.plot_y[bet.in_range], "o", color=color, ms=6,
                label="BET range (%d points)" % bet.n_points)
        xs = np.array([0.0, x_max])
        ax.plot(xs, bet.slope * xs + bet.intercept, "-", color="#222222", lw=1.0, label="Linear fit")
        ax.set_xlim(0, x_max)
        ax.set_ylim(bottom=0)
        flags = "".join("✓" if f else "✗" for f in bet.rouquerol_ok)
        text = ("S(BET) = %.1f m2/g\nVm = %.2f cm3 STP/g\nC = %.4g\nR2 = %.5f\n"
                "Range: %.2e – %.2e\nRouquerol 1-4: %s"
                % (bet.surface_m2_g, bet.vm_cm3_g, bet.c_constant, bet.r2, bet.p_min, bet.p_max, flags))
        ax.text(0.97, 0.05, text, transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", edgecolor="#cccccc"))
        ax.legend(frameon=False, fontsize=9, loc="upper left")
        return fig

    # Figure 4 ------------------------------------------------------------------
    def bars(self):
        """Bar chart: geometric surface area, experimental BET, simulated BET."""
        r = self.r
        values = [r.geometric.surface_m2_g if r.geometric else float("nan"),
                  r.bet_exp.surface_m2_g if r.bet_exp and r.bet_exp.valid else float("nan"),
                  r.bet_sim.surface_m2_g if r.bet_sim and r.bet_sim.valid else float("nan")]
        labels = ["Geometric\nsurface area", "Experimental\nBET", "Simulated\nBET"]
        colors = [COLOR_GEO, COLOR_EXP, COLOR_SIM]
        fig = self._new(5.6, 4.4)
        ax = fig.subplots()
        self._style(ax)
        ax.grid(False, axis="x")
        heights = [v if np.isfinite(v) else 0.0 for v in values]
        bars = ax.bar(labels, heights, color=colors, width=0.6)
        top = max(heights) if max(heights) > 0 else 1.0
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.015 * top,
                    "%.0f" % value if np.isfinite(value) else "n.d.", ha="center", va="bottom", fontsize=10)
        ax.set_ylim(0, top * 1.12)
        ax.set_ylabel("Specific surface area (m2/g)")
        ax.set_title("Specific surface areas compared", fontsize=11)
        return fig

    # Structure picture ----------------------------------------------------------
    def structure(self):
        """View of the P1 cell (orthographic projection along c).

        Documented alternative to a screenshot of the Mercury window, which the
        CSD Python API does not provide to scripts.
        """
        st = self.r.structure
        fig = self._new(4.6, 4.6)
        ax = fig.subplots()
        ax.set_aspect("equal")
        ax.axis("off")
        if not st.p1_elements:
            return fig
        cellm = st.cell_matrix
        c_dir = cellm[2] / np.linalg.norm(cellm[2])
        ex = cellm[0] - np.dot(cellm[0], c_dir) * c_dir
        ex /= np.linalg.norm(ex)
        ey = np.cross(c_dir, ex)
        cart = st.p1_cart
        x, y, depth = cart @ ex, cart @ ey, cart @ c_dir
        elements = np.array(st.p1_elements)
        if len(elements) <= 3000:                               # bonds (minimum image)
            radii = np.array([COVALENT_RADIUS.get(e, 1.5) for e in elements])
            pi, pj, pd_ = st.neighbor_pairs(st.p1_frac, 2.9)
            bonded = pd_ < 1.2 * (radii[pi] + radii[pj])
            delta = st.p1_frac[pj[bonded]] - st.p1_frac[pi[bonded]]
            inside = np.all(np.abs(delta) < 0.5, axis=1)        # bonds that do not cross the cell
            segments = [((x[i], y[i]), (x[j], y[j]))
                        for i, j in zip(pi[bonded][inside], pj[bonded][inside])]
            from matplotlib.collections import LineCollection
            ax.add_collection(LineCollection(segments, colors="#9a9a9a", linewidths=0.8, zorder=1))
        order = np.argsort(depth)
        sizes = np.array([18.0 + 60.0 * COVALENT_RADIUS.get(e, 1.0) ** 2 for e in elements])
        colors = [ELEMENT_COLORS.get(e, "#c080c0") for e in elements]
        ax.scatter(x[order], y[order], s=sizes[order] * 0.5, c=[colors[k] for k in order],
                   edgecolors="#333333", linewidths=0.3, zorder=2)
        corners = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0]]) @ cellm
        ax.plot(corners @ ex, corners @ ey, "-", color="#222222", lw=0.8, zorder=3)
        ax.set_title("%s — P1 cell viewed along c" % st.name, fontsize=10)
        ax.autoscale_view()
        return fig

    def all_figures(self) -> Dict[str, object]:
        return {"isotherm": self.isotherm(), "bet_exp": self.bet("exp"),
                "bet_sim": self.bet("sim"), "surfaces": self.bars()}


# =============================================================================
# 11. EXPORTS: EXCEL, CSV, PDF, HTML
# =============================================================================
class ReportBuilder:
    """Writes all the result files into a folder."""

    def __init__(self, results: AnalysisResults, out_dir: str,
                 log: Optional[Callable[[str], None]] = None):
        self.r = results
        self.out_dir = out_dir
        self.log = log or (lambda msg: None)
        self.figures = FigureFactory(results)
        os.makedirs(out_dir, exist_ok=True)
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", results.structure.name)[:60] or "structure"
        self.stem = "%s_%s" % (safe, results.settings.gas)

    def _path(self, suffix: str) -> str:
        return os.path.join(self.out_dir, self.stem + suffix)

    def _exp_frame(self) -> pd.DataFrame:
        if self.r.experimental is None:
            return pd.DataFrame(columns=["Branch", "P/P0", "Adsorption_cm3STPg", "Adsorption_mmol_g"])
        return self.r.experimental.to_frame()

    def _sim_frame(self) -> pd.DataFrame:
        if self.r.simulated is None:
            return pd.DataFrame(columns=["P/P0", "Adsorption_mmol_g", "Adsorption_cm3STPg"])
        return self.r.simulated

    # ---------------------------------------------------------------------- Excel
    def write_excel(self) -> Optional[str]:
        """Workbook: 1 experimental isotherm, 2 simulated isotherm, 3 BET, 4 comparison
        (+ Structure and Parameters sheets for traceability)."""
        path = self._path("_results.xlsx")
        engine = None
        for cand in ("openpyxl", "xlsxwriter"):
            try:
                importlib.import_module(cand)
                engine = cand
                break
            except ImportError:
                continue
        if engine is None:
            self.log("Excel export impossible: install 'openpyxl' (the CSV files are still written).")
            return None
        with pd.ExcelWriter(path, engine=engine) as writer:
            self._exp_frame().to_excel(writer, sheet_name="Isotherm_exp", index=False)
            self._sim_frame().to_excel(writer, sheet_name="Isotherm_sim", index=False)
            self.r.bet_frame().to_excel(writer, sheet_name="BET", index=False)
            points = self.r.bet_points_frame()
            if len(points):
                points.to_excel(writer, sheet_name="BET", index=False, startrow=len(self.r.bet_frame()) + 3)
            self.r.comparison_frame().to_excel(writer, sheet_name="Comparison", index=False)
            self.r.structure_frame().to_excel(writer, sheet_name="Structure", index=False)
            self.r.settings_frame().to_excel(writer, sheet_name="Parameters", index=False)
            if engine == "openpyxl":                            # column widths
                for sheet in writer.sheets.values():
                    for column in sheet.columns:
                        width = max(len(str(c.value)) if c.value is not None else 0 for c in column)
                        sheet.column_dimensions[column[0].column_letter].width = min(max(12, width + 2), 60)
        return path

    # ------------------------------------------------------------------------ CSV
    def write_csv(self) -> List[str]:
        paths = []
        for suffix, frame in (("_isotherm_exp.csv", self._exp_frame()),
                              ("_isotherm_sim.csv", self._sim_frame()),
                              ("_bet.csv", self.r.bet_frame()),
                              ("_comparison.csv", self.r.comparison_frame())):
            path = self._path(suffix)
            frame.to_csv(path, index=False, encoding="utf-8-sig")
            paths.append(path)
        return paths

    # ------------------------------------------------------------------------ PDF
    @staticmethod
    def _table(ax, frame: pd.DataFrame, col_widths: Sequence[float], font: float = 8.0) -> None:
        ax.axis("off")
        cells = [[c if isinstance(c, str) else fmt(c, 4 if abs(c) < 10 else 1) if isinstance(c, (float, np.floating))
                  else str(c) for c in row] for row in frame.values.tolist()]
        table = ax.table(cellText=cells, colLabels=list(frame.columns), colWidths=list(col_widths),
                         cellLoc="left", loc="upper center")
        table.auto_set_font_size(False)
        table.set_fontsize(font)
        table.scale(1.0, 1.35)
        for (row, _col), cell in table.get_celld().items():
            cell.set_edgecolor("#bbbbbb")
            if row == 0:
                cell.set_facecolor("#e9eef5")
                cell.set_text_props(weight="bold")

    def write_pdf(self) -> str:
        """PDF report: structure, crystallographic parameters, isotherms, BET, comparison."""
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        r = self.r
        path = self._path("_report.pdf")
        a4 = (8.27, 11.69)
        with PdfPages(path) as pdf:
            # --- page 1: structure ---
            fig = Figure(figsize=a4)
            FigureCanvasAgg(fig)
            fig.text(0.5, 0.955, "Porosity report: experiment and simulation", ha="center",
                     fontsize=15, weight="bold")
            fig.text(0.5, 0.93, "%s — %s — %s" % (r.structure.name, ADSORBATES[r.settings.gas]["label"],
                                                   r.created), ha="center", fontsize=10, color="#444444")
            ax_img = fig.add_axes([0.2, 0.52, 0.6, 0.39])
            self._draw_structure(ax_img)
            ax_tab = fig.add_axes([0.06, 0.05, 0.88, 0.45])
            self._table(ax_tab, r.structure_frame(), (0.42, 0.58))
            pdf.savefig(fig)
            # --- page 2: isotherms + parameters ---
            fig = Figure(figsize=a4)
            FigureCanvasAgg(fig)
            fig.text(0.06, 0.955, "Adsorption isotherms", fontsize=13, weight="bold")
            self._embed(fig, self.figures.isotherm(), [0.04, 0.58, 0.92, 0.35])
            ax_tab = fig.add_axes([0.06, 0.04, 0.88, 0.5])
            self._table(ax_tab, r.settings_frame().astype(str), (0.42, 0.58))
            pdf.savefig(fig)
            # --- page 3: BET ---
            fig = Figure(figsize=a4)
            FigureCanvasAgg(fig)
            fig.text(0.06, 0.955, "BET analysis (Rouquerol criteria)", fontsize=13, weight="bold")
            self._embed(fig, self.figures.bet("exp"), [0.03, 0.6, 0.47, 0.33])
            self._embed(fig, self.figures.bet("sim"), [0.5, 0.6, 0.47, 0.33])
            ax_tab = fig.add_axes([0.06, 0.05, 0.88, 0.52])
            self._table(ax_tab, r.bet_frame(), (0.46, 0.27, 0.27))
            pdf.savefig(fig)
            # --- page 4: comparison ---
            fig = Figure(figsize=a4)
            FigureCanvasAgg(fig)
            fig.text(0.06, 0.955, "Comparison: experiment / simulation / geometry", fontsize=13, weight="bold")
            self._embed(fig, self.figures.bars(), [0.2, 0.58, 0.6, 0.35])
            ax_tab = fig.add_axes([0.04, 0.2, 0.92, 0.35])
            self._table(ax_tab, r.comparison_frame(), (0.36, 0.15, 0.17, 0.14, 0.18), font=7.5)
            notes = self._notes()
            fig.text(0.06, 0.18, "\n".join(notes), fontsize=8, va="top", color="#333333", wrap=True)
            pdf.savefig(fig)
        return path

    def _notes(self) -> List[str]:
        r = self.r
        notes = []
        if r.geometric:
            notes.append("Geometry: %s; probe %.2f Å; %s." % (r.geometric.method, r.geometric.probe_diameter,
                                                                 r.geometric.pore_volume_label))
            notes.extend(self._wrap(n) for n in r.geometric.notes)
        for name, bet in (("Experimental BET", r.bet_exp), ("Simulated BET", r.bet_sim)):
            if bet is not None and bet.message:
                notes.append("%s: %s" % (name, bet.message))
        notes.append("Simulation: absolute adsorption, rigid framework; experiment: excess adsorption.")
        return notes

    @staticmethod
    def _wrap(text: str, width: int = 120) -> str:
        import textwrap
        return "\n".join(textwrap.wrap(text, width))

    def _draw_structure(self, ax) -> None:
        self._embed(ax.figure, self.figures.structure(), ax.get_position().bounds)
        ax.axis("off")

    @staticmethod
    def _png_bytes(fig, dpi: int = 200) -> bytes:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        FigureCanvasAgg(fig)
        buffer = _io.BytesIO()
        fig.savefig(buffer, format="png", dpi=dpi)
        return buffer.getvalue()

    def _embed(self, page, fig, rect) -> None:
        """Inserts a figure (rendered as PNG) into a page of the PDF report."""
        import matplotlib.image as mpimg
        image = mpimg.imread(_io.BytesIO(self._png_bytes(fig)), format="png")
        ax = page.add_axes(rect)
        ax.imshow(image)
        ax.axis("off")

    # ----------------------------------------------------------------------- PNG
    def write_png(self) -> List[str]:
        paths = []
        figures = dict(self.figures.all_figures(), structure=self.figures.structure())
        for name, fig in figures.items():
            path = self._path("_%s.png" % name)
            with open(path, "wb") as handle:
                handle.write(self._png_bytes(fig, dpi=160))
            paths.append(path)
        return paths

    # ---------------------------------------------------------------------- HTML
    def html(self) -> str:
        """Self-contained HTML report (embedded images) shown in Mercury."""
        r = self.r

        def table(frame: pd.DataFrame) -> str:
            return frame.to_html(index=False, border=0, na_rep="—", float_format=lambda v: "%.4g" % v,
                                 classes="t")

        def image(fig) -> str:
            return '<img src="data:image/png;base64,%s" style="max-width:100%%"/>' % \
                base64.b64encode(self._png_bytes(fig, dpi=110)).decode("ascii")

        parts = [
            "<html><head><meta charset='utf-8'><style>",
            "body{font-family:Helvetica,Arial,sans-serif;margin:18px;color:#222}",
            "h1{font-size:20px}h2{font-size:15px;border-bottom:1px solid #ccc;padding-bottom:3px;margin-top:22px}",
            "table.t{border-collapse:collapse;font-size:12px}table.t td,table.t th{border:1px solid #ccc;"
            "padding:3px 8px;text-align:left}table.t th{background:#e9eef5}",
            "</style></head><body>",
            "<h1>%s — %s</h1><p>%s, %s</p>" % (APP_NAME, r.structure.name,
                                               ADSORBATES[r.settings.gas]["label"], r.created),
            "<h2>Comparison</h2>", table(r.comparison_frame()), image(self.figures.bars()),
            "<h2>Isotherms</h2>", image(self.figures.isotherm()),
            "<h2>BET analysis</h2>", image(self.figures.bet("exp")), image(self.figures.bet("sim")),
            table(r.bet_frame().astype(str)),
            "<h2>Structure</h2>", image(self.figures.structure()), table(r.structure_frame()),
            "<h2>Parameters</h2>", table(r.settings_frame().astype(str)),
            "<p>" + "<br/>".join(n.replace("\n", " ") for n in self._notes()) + "</p>",
            "<p>Files: %s</p>" % self.out_dir, "</body></html>"]
        return "\n".join(parts)

    def write_all(self) -> Dict[str, object]:
        """Writes Excel, CSV, PDF, PNG and HTML; returns the created paths."""
        written: Dict[str, object] = {"csv": self.write_csv()}
        written["excel"] = self.write_excel()
        written["pdf"] = self.write_pdf()
        written["png"] = self.write_png()
        html_path = self._path("_report.html")
        with open(html_path, "w", encoding="utf-8") as handle:
            handle.write(self.html())
        written["html"] = html_path
        self.log("Results exported to: %s" % self.out_dir)
        return written


# =============================================================================
# 12. MERCURY LINK (CSD Python API)
# =============================================================================
class MercuryBridge:
    """Access to the structure open in Mercury.

    Mercury runs the script in a separate Python process and passes it, as
    arguments, the identifier of the current structure and the output files.
    Two generations of interface exist:
      * `ccdc.utilities.ApplicationInterface` (recent versions);
      * `mercury_interface.MercuryInterface` (older versions).
    If neither can be used (script run outside Mercury), `available` stays
    False and the program falls back to reading a CIF file directly.
    """

    def __init__(self):
        self.interface = None
        self.available = False
        self.error = ""
        saved_stderr = sys.stderr
        try:
            sys.stderr = _io.StringIO()                         # keep the API's argparse quiet
            try:
                from ccdc.utilities import ApplicationInterface
                interface = ApplicationInterface(parse_commandline=False)
                interface.parse_commandline()
            except ImportError as first:                        # older interface, or broken ccdc
                self.error = "%s: %s" % (type(first).__name__, first)
                from mercury_interface import MercuryInterface
                interface = MercuryInterface()
                self.error = ""
            self.interface = interface
            self.available = True
        except BaseException as exc:                            # including argparse's SystemExit
            second = "%s: %s" % (type(exc).__name__, exc)
            self.error = second if not self.error or self.error == second else self.error + " | " + second
        finally:
            sys.stderr = saved_stderr

    def current_entry(self):
        entry = getattr(self.interface, "current_entry")
        return entry() if callable(entry) else entry

    def original_cif(self) -> str:
        """Path of the CIF file loaded in Mercury, if it is among the options."""
        values = []
        for name in ("database_path", "input_file", "options"):
            value = getattr(self.interface, name, None)
            values.extend(value.values() if isinstance(value, dict) else [value])
        values.extend(sys.argv[1:])
        for value in values:
            if isinstance(value, str) and value.lower().endswith(".cif") and os.path.isfile(value):
                return value
        return ""

    def export_temp_cif(self, entry, directory: str) -> str:
        """Exports the current structure to a temporary CIF (CSD Python API)."""
        from ccdc import io as ccdc_io
        path = os.path.join(directory, "%s_mercury_export.cif" % re.sub(r"\W+", "_", str(entry.identifier)))
        with ccdc_io.CrystalWriter(path) as writer:
            writer.write(entry.crystal)
        return path

    def output_html_file(self) -> str:
        return str(getattr(self.interface, "output_html_file", "") or "") if self.available else ""

    def output_directory(self) -> str:
        if not self.available:
            return ""
        value = getattr(self.interface, "output_directory_path", "") or ""
        if not value and isinstance(getattr(self.interface, "options", None), dict):
            value = self.interface.options.get("working_directory_path", "")
        return str(value)


def ccdc_installed() -> bool:
    """True if the CSD Python API is installed (without loading it)."""
    import importlib.util
    try:
        return importlib.util.find_spec("ccdc") is not None
    except (ImportError, ValueError):
        return False


class ChildBridge:
    """Replaces `MercuryBridge` in the GUI process (which has no CSD API)."""

    def __init__(self, html_out: str = "", from_mercury: bool = False):
        self._html_out = html_out
        self.available = from_mercury

    def output_html_file(self) -> str:
        return self._html_out


def load_structure(cif_path: str = "") -> Tuple[StructureData, MercuryBridge, str]:
    """Loads the structure: from Mercury first, otherwise from the given CIF file.

    Returns (structure, Mercury bridge, path of the exported temporary CIF).
    """
    bridge = MercuryBridge()
    temp_cif = ""
    if bridge.available:
        try:
            entry = bridge.current_entry()
            temp_dir = tempfile.mkdtemp(prefix="raspa_mercury_")
            try:
                temp_cif = bridge.export_temp_cif(entry, temp_dir)
            except Exception as exc:
                bridge.error = "CIF export through the API failed: %s" % exc
            original = bridge.original_cif()
            structure = StructureData.from_csd_entry(entry, original)
            structure.csd_cif_path = temp_cif or original
            if temp_cif:
                structure.notes.append("Temporary CIF exported by the CSD Python API: %s" % temp_cif)
            if bridge.error:
                structure.notes.append(bridge.error)
            return structure, bridge, temp_cif
        except Exception as exc:                                # no usable current structure
            bridge.error = "Mercury structure unreadable: %s" % exc
            if not cif_path:
                raise RuntimeError(bridge.error) from exc
    if cif_path:
        structure = StructureData.from_cif_file(cif_path)
        structure.csd_cif_path = cif_path if ccdc_installed() else ""
        return structure, bridge, temp_cif
    raise RuntimeError("No structure: run the script from Mercury or give a CIF file.\n"
                       "(Mercury interface: %s)" % (bridge.error or "unavailable"))


def default_work_dir(structure: StructureData, gas: str, base: str = "") -> str:
    base = base or os.path.join(os.path.expanduser("~"), "RASPA_Mercury")
    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", structure.name)[:40] or "structure"
    return os.path.join(base, "%s_%s_%s" % (safe, gas, stamp))


# =============================================================================
# 13. Qt GRAPHICAL INTERFACE (PyQt6 / PyQt5 / PySide6 / PySide2)
# =============================================================================
# IMPORTANT (macOS especially): the CSD Python API ships its own Qt libraries
# (ccdc/_lib/Qt*.framework). Loading them in the same process as PyQt/PySide
# causes a symbol clash ("Symbol not found: QWidgetPrivate...").
# The script therefore runs as TWO processes:
#   * parent process (started by Mercury): CSD Python API, never Qt;
#   * child process (--gui-child option): Qt GUI, never the CSD Python API.
# The structure is passed from one to the other through a transfer CIF file.
# On macOS, the mere PRESENCE of PyQt in the site-packages of Mercury's Python
# can be enough to break the CSD Python API. PyQt is therefore installed
# separately, in the private folder QT_PRIVATE_DIR (pip install --target),
# which only the GUI process adds to its search path.
QT_CANDIDATES = ("PyQt6", "PyQt5", "PySide6", "PySide2")
QT_PRIVATE_DIR = os.environ.get("RASPA_POROSITY_QT_DIR") or \
    os.path.join(os.path.expanduser("~"), "RASPA_Mercury", "qt_libs")


def qt_in_private_dir() -> str:
    """Qt binding present in the private folder ('' if none)."""
    for name in QT_CANDIDATES:
        if os.path.isdir(os.path.join(QT_PRIVATE_DIR, name)):
            return name
    return ""


def qt_in_site_packages() -> str:
    """Qt binding installed in the current Python, WITHOUT loading it ('' if none)."""
    import importlib.util
    for name in QT_CANDIDATES:
        try:
            if importlib.util.find_spec(name) is not None:
                return name
        except (ImportError, ValueError):
            continue
    return ""


def detect_qt_binding() -> str:
    """Qt binding usable by the GUI process ('' if none)."""
    return qt_in_private_dir() or qt_in_site_packages()


def qt_fix_instructions() -> str:
    """Exact Terminal commands to install Qt without disturbing the CSD Python API."""
    python = sys.executable
    lines = ["Copy these commands into the Terminal (one at a time), then run the script again:", ""]
    if qt_in_site_packages():
        lines += ["1) remove PyQt from Mercury's Python (it disturbs the CSD Python API):",
                  '"%s" -m pip uninstall -y PyQt6 PyQt6-Qt6 PyQt6-sip PyQt5 PyQt5-Qt5 PyQt5-sip PySide6 PySide2'
                  % python, ""]
    lines += ["%s install PyQt6 into a separate folder:" % ("2)" if qt_in_site_packages() else "1)"),
              '"%s" -m pip install --target "%s" PyQt6' % (python, QT_PRIVATE_DIR), "",
              "%s additional modules (safe for Mercury):" % ("3)" if qt_in_site_packages() else "2)"),
              '"%s" -m pip install matplotlib openpyxl xlrd scipy' % python, "",
              "(%s %s, Python %s)" % (APP_NAME, APP_VERSION, sys.version.split()[0])]
    return "\n".join(lines)


_LOAD_QT = "--gui-child" in sys.argv[1:] or os.environ.get("RASPA_POROSITY_QT") == "1"
if _LOAD_QT and os.path.isdir(QT_PRIVATE_DIR):
    sys.path.insert(0, QT_PRIVATE_DIR)                          # PyQt from the private folder first
QT_BINDING = ""
QtWidgets = QtCore = QtGui = None
for _binding in (QT_CANDIDATES if _LOAD_QT else ()):
    try:
        QtWidgets = importlib.import_module(_binding + ".QtWidgets")
        QtCore = importlib.import_module(_binding + ".QtCore")
        QtGui = importlib.import_module(_binding + ".QtGui")
        QT_BINDING = _binding
        break
    except ImportError:
        continue

if QT_BINDING:
    os.environ.setdefault("QT_API", QT_BINDING.lower())         # same binding for matplotlib
    import matplotlib
    try:
        matplotlib.use("QtAgg")
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
    except Exception:                                           # matplotlib < 3.5
        matplotlib.use("Qt5Agg")
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
    Signal = getattr(QtCore, "pyqtSignal", None) or getattr(QtCore, "Signal")

    class Worker(QtCore.QThread):
        """Runs a function in a separate thread so the GUI does not freeze.

        The function receives two callbacks: log(text) and progress(fraction, text).
        """
        log = Signal(str)
        progress = Signal(float, str)
        done = Signal(object)
        failed = Signal(str)

        def __init__(self, function: Callable, parent=None):
            super().__init__(parent)
            self._function = function

        def run(self) -> None:                                  # called by QThread.start()
            try:
                self.done.emit(self._function(self.log.emit, self.progress.emit))
            except RaspaError as exc:
                self.failed.emit(str(exc))
            except Exception:
                self.failed.emit(traceback.format_exc())

    class FigureTab(QtWidgets.QWidget):
        """Tab holding a matplotlib figure and its toolbar."""

        def __init__(self, parent=None):
            super().__init__(parent)
            self._layout = QtWidgets.QVBoxLayout(self)
            self._canvas = None
            self._toolbar = None

        def set_figure(self, figure) -> None:
            for widget in (self._toolbar, self._canvas):
                if widget is not None:
                    self._layout.removeWidget(widget)
                    widget.setParent(None)
                    widget.deleteLater()
            self._canvas = FigureCanvasQTAgg(figure)
            self._toolbar = NavigationToolbar2QT(self._canvas, self)
            self._layout.addWidget(self._toolbar)
            self._layout.addWidget(self._canvas)
            self._canvas.draw_idle()

    class MainWindow(QtWidgets.QMainWindow):
        """Main window: Structure, Experimental data, RASPA settings,
        Calculation, Results."""

        def __init__(self, structure: StructureData, bridge: Optional[ChildBridge] = None,
                     cli: Optional[argparse.Namespace] = None):
            super().__init__()
            self.structure = structure
            self.bridge = bridge
            self.experimental: Optional[Isotherm] = None
            self.geometric: Optional[GeometricResult] = None
            self.results: Optional[AnalysisResults] = None
            self.worker: Optional[Worker] = None
            self.runner: Optional[RaspaRunner] = None
            self.last_html = ""
            self.qsettings = QtCore.QSettings("RaspaMercury", "PorosityAnalyser")
            self.setWindowTitle("%s %s — %s" % (APP_NAME, APP_VERSION, structure.name))
            self.resize(1080, 760)
            self.tabs = QtWidgets.QTabWidget()
            self.setCentralWidget(self.tabs)
            self._build_structure_tab()
            self._build_experimental_tab()
            self._build_parameters_tab()
            self._build_run_tab()
            self._build_results_tab()
            self.statusBar().showMessage("Qt: %s — %s" % (
                QT_BINDING, "structure supplied by Mercury" if bridge and bridge.available
                else "standalone mode (internal CIF reader)"))
            self._rebuild_model()
            self._on_gas_changed()
            if cli is not None:
                self._apply_cli(cli)

        # ------------------------------------------------------------------- helpers
        def _error(self, text: str) -> None:
            QtWidgets.QMessageBox.critical(self, APP_NAME, text)

        @staticmethod
        def _html_table(frame: pd.DataFrame) -> str:
            html = frame.to_html(index=False, border=0, na_rep="—",
                                 float_format=lambda v: "%.4g" % v)
            return ("<style>table{border-collapse:collapse}td,th{border:1px solid #bbb;"
                    "padding:3px 8px;text-align:left}th{background:#e9eef5}</style>" + html)

        def _apply_cli(self, cli: argparse.Namespace) -> None:
            if getattr(cli, "gas", None):
                self.gas_combo.setCurrentText(cli.gas)
            if getattr(cli, "raspa", None):
                self.raspa_edit.setText(cli.raspa)
            if getattr(cli, "exp", None):
                self._load_experimental_file(cli.exp)

        # ============================== Tab 1: Structure ==============================
        def _build_structure_tab(self) -> None:
            page = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout(page)
            self.structure_view = QtWidgets.QTextBrowser()
            layout.addWidget(self.structure_view, 3)
            box = QtWidgets.QGroupBox("Preparation of the simulation model (P1 cell)")
            form = QtWidgets.QVBoxLayout(box)
            self.disorder_check = QtWidgets.QCheckBox(
                "Resolve disorder (a single component for partially occupied sites)")
            self.disorder_check.setChecked(True)
            self.guest_check = QtWidgets.QCheckBox(
                "Remove guest species not bonded to the framework (solvent, pseudo-located water)")
            self.guest_check.setChecked(True)
            for check in (self.disorder_check, self.guest_check):
                check.stateChanged.connect(lambda _state: self._rebuild_model())
                form.addWidget(check)
            layout.addWidget(box)
            row = QtWidgets.QHBoxLayout()
            open_button = QtWidgets.QPushButton("Open another CIF file…")
            open_button.clicked.connect(self._open_cif)
            self.geometry_button = QtWidgets.QPushButton("Compute geometric properties")
            self.geometry_button.clicked.connect(self._run_geometry)
            self.geometry_mode = QtWidgets.QComboBox()
            self.geometry_mode.addItems(["Auto (CSD Pore Analyser if possible)", "Pore Analyser (CSD Python API)",
                                         "Internal calculation"])
            row.addWidget(open_button)
            row.addStretch(1)
            row.addWidget(QtWidgets.QLabel("Method:"))
            row.addWidget(self.geometry_mode)
            row.addWidget(self.geometry_button)
            layout.addLayout(row)
            self.geometry_view = QtWidgets.QTextBrowser()
            layout.addWidget(self.geometry_view, 2)
            self.tabs.addTab(page, "1. Structure")

        def _open_cif(self) -> None:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "CIF file", "", "CIF (*.cif)")
            if not path:
                return
            try:
                self.structure = StructureData.from_cif_file(path)
                self.structure.csd_cif_path = path if ccdc_installed() else ""
            except Exception as exc:
                self._error("Cannot read the CIF:\n%s" % exc)
                return
            self.setWindowTitle("%s %s — %s" % (APP_NAME, APP_VERSION, self.structure.name))
            self._rebuild_model()

        def _rebuild_model(self) -> None:
            """Rebuilds the P1 cell according to the options and refreshes the display."""
            try:
                self.structure.build_p1(self.disorder_check.isChecked(), self.guest_check.isChecked())
            except Exception as exc:
                self._error("Cannot build the P1 model:\n%s" % exc)
                return
            self.geometric = None
            self.geometry_view.setHtml("<i>Geometric properties not computed yet "
                                       "(they will also be computed during the analysis).</i>")
            dummy = AnalysisResults(self.structure, SimulationSettings())
            notes = self.structure.notes + self.structure.p1_log
            self.structure_view.setHtml(self._html_table(dummy.structure_frame()) +
                                        "<p>" + "<br/>".join(notes) + "</p>")
            if hasattr(self, "cell_spins"):
                self._auto_supercell()

        def _probe_diameter(self) -> float:
            return ADSORBATES[self.gas_combo.currentText()]["probe_diameter"]

        def _geometry_mode(self) -> str:
            return ("auto", "csd", "internal")[self.geometry_mode.currentIndex()]

        def _show_geometry(self, geo: GeometricResult) -> None:
            self.geometric = geo
            frame = pd.DataFrame([
                ("Method", geo.method), ("Probe diameter (Å)", fmt(geo.probe_diameter, 3)),
                ("Accessible surface area (m2/g)", fmt(geo.surface_m2_g, 1)),
                ("Accessible surface area (m2/cm3)", fmt(geo.surface_m2_cm3, 1)),
                ("Pore volume (cm3/g) — " + geo.pore_volume_label, fmt(geo.pore_volume_cm3_g, 4)),
                ("Geometric volume, point probe (cm3/g)", fmt(geo.geometric_volume_cm3_g, 4)),
                ("Void fraction", fmt(geo.void_fraction, 4)),
                ("Pore limiting diameter PLD (Å)", fmt(geo.pld, 2)),
                ("Largest cavity diameter LCD (Å)", fmt(geo.lcd, 2)),
            ], columns=["Geometric property", "Value"])
            self.geometry_view.setHtml(self._html_table(frame) + "<p>" + "<br/>".join(geo.notes) + "</p>")

        def _run_geometry(self) -> None:
            if self.worker is not None:
                return
            analyser = GeometricAnalyzer(self.structure, self._probe_diameter())
            mode = self._geometry_mode()

            def job(log, progress):
                analyser.log = log
                return analyser.run(mode)

            self.geometry_button.setEnabled(False)
            self.geometry_view.setHtml("<i>Calculation in progress…</i>")
            self._start_worker(job, self._geometry_done)

        def _geometry_done(self, geo: GeometricResult) -> None:
            self.geometry_button.setEnabled(True)
            self._show_geometry(geo)

        # =========================== Tab 2: Experimental data ===========================
        def _build_experimental_tab(self) -> None:
            page = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout(page)
            box_a = QtWidgets.QGroupBox("Option A — Import a file (xlsx, xls, csv)")
            row = QtWidgets.QHBoxLayout(box_a)
            self.exp_path = QtWidgets.QLineEdit()
            self.exp_path.setPlaceholderText("Expected columns: P/P0 and Adsorption_cm3STPg "
                                             "(BELSORP exports are recognised)")
            browse = QtWidgets.QPushButton("Browse…")
            browse.clicked.connect(self._browse_experimental)
            row.addWidget(self.exp_path)
            row.addWidget(browse)
            layout.addWidget(box_a)
            box_b = QtWidgets.QGroupBox("Option B — Paste a table (P/P0, then adsorption in cm3 STP/g)")
            col = QtWidgets.QVBoxLayout(box_b)
            self.paste_edit = QtWidgets.QPlainTextEdit()
            self.paste_edit.setPlaceholderText("P/P0 Adsorption\n0.0001 15.3\n0.0010 31.3\n0.0100 36.6\n"
                                               "0.0500 37.7\n0.1000 38.3\n0.5000 40.0\n0.9990 56.7")
            convert = QtWidgets.QPushButton("Convert the pasted table")
            convert.clicked.connect(self._convert_pasted)
            col.addWidget(self.paste_edit)
            col.addWidget(convert)
            layout.addWidget(box_b, 2)
            self.exp_view = QtWidgets.QTextBrowser()
            self.exp_view.setHtml("<i>No experimental data loaded.</i>")
            layout.addWidget(self.exp_view, 3)
            self.tabs.addTab(page, "2. Experimental data")

        def _browse_experimental(self) -> None:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, "Experimental isotherm", "", "Spreadsheets (*.xlsx *.xls *.csv *.txt);;All files (*)")
            if path:
                self._load_experimental_file(path)

        def _load_experimental_file(self, path: str) -> None:
            try:
                iso = ExperimentalDataLoader.load_file(path)
            except Exception as exc:
                self._error("Cannot read the experimental data:\n%s" % exc)
                return
            self.exp_path.setText(path)
            self._set_experimental(iso)

        def _convert_pasted(self) -> None:
            try:
                iso = ExperimentalDataLoader.parse_text(self.paste_edit.toPlainText())
            except Exception as exc:
                self._error("Unreadable table:\n%s" % exc)
                return
            self._set_experimental(iso)

        def _set_experimental(self, iso: Isotherm) -> None:
            self.experimental = iso
            bet = BETAnalyzer(self._float(self.cross_edit, 0.162)).analyse(iso.p_rel, iso.v_cm3)
            lines = ["<b>%s</b>: %d adsorption points, %d desorption points; P/P0 from %.3g to %.3g."
                     % (iso.label, len(iso.p_rel), len(iso.des_p_rel), iso.p_rel.min(), iso.p_rel.max())]
            if bet.valid:
                lines.append("BET (Rouquerol): <b>%.1f m2/g</b>, Vm = %.2f cm3 STP/g, C = %.4g, range "
                             "%.2e – %.2e. %s" % (bet.surface_m2_g, bet.vm_cm3_g, bet.c_constant,
                                                  bet.p_min, bet.p_max, bet.message))
            meta = iso.meta
            if meta:
                lines.append("File metadata: " + ", ".join("%s = %s" % kv for kv in meta.items()))
            gas = str(meta.get("adsorptive", "")).upper()
            if gas in ADSORBATES and gas != self.gas_combo.currentText():
                lines.append("<span style='color:#b00'>The file states the adsorbate is %s: remember to "
                             "select it in the Settings tab.</span>" % gas)
            self.exp_view.setHtml("<p>" + "</p><p>".join(lines) + "</p>" +
                                  self._html_table(iso.to_frame()))

        # ============================ Tab 3: RASPA settings ============================
        @staticmethod
        def _float(edit, default: float) -> float:
            value = to_float(edit.text())
            return value if np.isfinite(value) else default

        def _build_parameters_tab(self) -> None:
            page = QtWidgets.QWidget()
            outer = QtWidgets.QHBoxLayout(page)
            left, right = QtWidgets.QFormLayout(), QtWidgets.QFormLayout()
            box_left = QtWidgets.QGroupBox("Adsorbate and simulation")
            box_left.setLayout(left)
            box_right = QtWidgets.QGroupBox("Isotherm, analysis and execution")
            box_right.setLayout(right)
            outer.addWidget(box_left)
            outer.addWidget(box_right)

            self.gas_combo = QtWidgets.QComboBox()
            self.gas_combo.addItems(list(ADSORBATES))
            self.gas_combo.currentIndexChanged.connect(lambda _i: self._on_gas_changed())
            self.gas_info = QtWidgets.QLabel()
            self.gas_info.setWordWrap(True)
            self.temp_spin = QtWidgets.QDoubleSpinBox()
            self.temp_spin.setRange(1.0, 2000.0)
            self.temp_spin.setDecimals(2)
            self.temp_spin.setSuffix(" K")
            self.p0_edit = QtWidgets.QLineEdit()
            self.init_spin = QtWidgets.QSpinBox()
            self.init_spin.setRange(0, 100000000)
            self.init_spin.setValue(2000)
            self.prod_spin = QtWidgets.QSpinBox()
            self.prod_spin.setRange(1, 100000000)
            self.prod_spin.setValue(5000)
            self.cell_spins = []
            cell_row = QtWidgets.QHBoxLayout()
            for _ in range(3):
                spin = QtWidgets.QSpinBox()
                spin.setRange(1, 20)
                self.cell_spins.append(spin)
                cell_row.addWidget(spin)
            auto = QtWidgets.QPushButton("Auto")
            auto.setToolTip("Smallest supercell compatible with the cut-off radius (minimum image)")
            auto.clicked.connect(self._auto_supercell)
            cell_row.addWidget(auto)
            self.ff_combo = QtWidgets.QComboBox()
            self.ff_combo.addItems(FORCE_FIELDS)
            self.charge_combo = QtWidgets.QComboBox()
            self.charge_combo.addItems(CHARGE_MODES)
            self.cutoff_spin = QtWidgets.QDoubleSpinBox()
            self.cutoff_spin.setRange(6.0, 30.0)
            self.cutoff_spin.setValue(12.0)
            self.cutoff_spin.setSuffix(" Å")
            self.cutoff_spin.valueChanged.connect(lambda _v: self._auto_supercell())
            left.addRow("Adsorbate:", self.gas_combo)
            left.addRow("", self.gas_info)
            left.addRow("Temperature:", self.temp_spin)
            left.addRow("Saturation pressure P0 (Pa):", self.p0_edit)
            left.addRow("Equilibration cycles:", self.init_spin)
            left.addRow("Production cycles:", self.prod_spin)
            left.addRow("Supercell (a, b, c):", cell_row)
            left.addRow("Force field:", self.ff_combo)
            left.addRow("Framework charges:", self.charge_combo)
            left.addRow("Cut-off radius:", self.cutoff_spin)

            self.pmin_edit = QtWidgets.QLineEdit("1e-6")
            self.pmax_edit = QtWidgets.QLineEdit("1.0")
            self.npoints_spin = QtWidgets.QSpinBox()
            self.npoints_spin.setRange(3, 200)
            self.npoints_spin.setValue(35)
            self.cross_edit = QtWidgets.QLineEdit("0.162")
            self.gurvich_edit = QtWidgets.QLineEdit("0.95")
            self.jobs_spin = QtWidgets.QSpinBox()
            self.jobs_spin.setRange(1, 256)
            self.jobs_spin.setValue(max(1, (os.cpu_count() or 2) - 1))
            self.raspa_edit = QtWidgets.QLineEdit(
                str(self.qsettings.value("raspa_path", "") or "") or find_raspa_executable())
            self.raspa_edit.setPlaceholderText("/usr/local/bin/raspa3")
            raspa_row = QtWidgets.QHBoxLayout()
            raspa_row.addWidget(self.raspa_edit)
            raspa_browse = QtWidgets.QPushButton("Browse…")
            raspa_browse.clicked.connect(self._browse_raspa)
            raspa_row.addWidget(raspa_browse)
            self.out_edit = QtWidgets.QLineEdit(
                str(self.qsettings.value("output_base", "") or "") or
                os.path.join(os.path.expanduser("~"), "RASPA_Mercury"))
            out_row = QtWidgets.QHBoxLayout()
            out_row.addWidget(self.out_edit)
            out_browse = QtWidgets.QPushButton("Browse…")
            out_browse.clicked.connect(self._browse_output)
            out_row.addWidget(out_browse)
            right.addRow("Minimum P/P0:", self.pmin_edit)
            right.addRow("Maximum P/P0:", self.pmax_edit)
            right.addRow("Number of points (log scale):", self.npoints_spin)
            right.addRow("BET cross-sectional area (nm2):", self.cross_edit)
            right.addRow("P/P0 for the pore volume (Gurvich):", self.gurvich_edit)
            right.addRow("Simultaneous simulations:", self.jobs_spin)
            right.addRow("RASPA3 executable:", raspa_row)
            right.addRow("Results folder:", out_row)
            self.tabs.addTab(page, "3. RASPA settings")

        def _browse_raspa(self) -> None:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "raspa3 executable", "/usr/local/bin")
            if path:
                self.raspa_edit.setText(path)

        def _browse_output(self) -> None:
            path = QtWidgets.QFileDialog.getExistingDirectory(self, "Results folder", self.out_edit.text())
            if path:
                self.out_edit.setText(path)

        def _on_gas_changed(self) -> None:
            """Updates T, P0, BET cross-section and charge mode for the selected adsorbate."""
            gas = ADSORBATES[self.gas_combo.currentText()]
            self.temp_spin.setValue(gas["temperature"])
            self.p0_edit.setText("%.0f" % gas["p0"])
            self.cross_edit.setText("%.3f" % gas["cross_section"])
            charged = any(abs(pa["charge"]) > 0 for pa in gas["pseudo_atoms"])
            self.charge_combo.setEnabled(charged)
            self.charge_combo.setCurrentIndex(0 if charged else 2)
            self.gas_info.setText("Model: %s. Geometric probe: %.2f Å." % (gas["model"], gas["probe_diameter"]))
            if self.geometric is not None:                      # the probe depends on the adsorbate
                self.geometric = None
                self.geometry_view.setHtml("<i>Adsorbate changed: geometric properties must be recomputed.</i>")

        def _auto_supercell(self) -> None:
            for spin, value in zip(self.cell_spins, self.structure.minimal_supercell(self.cutoff_spin.value())):
                spin.setValue(value)

        def _collect_settings(self) -> SimulationSettings:
            gas = self.gas_combo.currentText()
            return SimulationSettings(
                gas=gas, temperature=self.temp_spin.value(),
                p0=self._float(self.p0_edit, ADSORBATES[gas]["p0"]),
                n_init_cycles=self.init_spin.value(), n_prod_cycles=self.prod_spin.value(),
                supercell=tuple(spin.value() for spin in self.cell_spins),
                forcefield=self.ff_combo.currentText(), charge_mode=self.charge_combo.currentText(),
                cutoff=self.cutoff_spin.value(),
                p_rel_min=self._float(self.pmin_edit, 1e-6), p_rel_max=self._float(self.pmax_edit, 1.0),
                n_points=self.npoints_spin.value(),
                cross_section=self._float(self.cross_edit, ADSORBATES[gas]["cross_section"]),
                gurvich_p_rel=self._float(self.gurvich_edit, 0.95),
                raspa_path=self.raspa_edit.text().strip(), n_jobs=self.jobs_spin.value(),
                work_dir=default_work_dir(self.structure, gas, self.out_edit.text().strip()))

        # ============================== Tab 4: Calculation ==============================
        def _build_run_tab(self) -> None:
            page = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout(page)
            row = QtWidgets.QHBoxLayout()
            self.run_button = QtWidgets.QPushButton("Run RASPA Simulation")
            self.run_button.setMinimumHeight(38)
            self.run_button.clicked.connect(lambda: self._start_analysis(True))
            self.quick_button = QtWidgets.QPushButton("Analyse without simulation (experiment + geometry)")
            self.quick_button.clicked.connect(lambda: self._start_analysis(False))
            self.stop_button = QtWidgets.QPushButton("Stop")
            self.stop_button.setEnabled(False)
            self.stop_button.clicked.connect(self._stop)
            row.addWidget(self.run_button, 2)
            row.addWidget(self.quick_button, 2)
            row.addWidget(self.stop_button, 1)
            layout.addLayout(row)
            self.progress_bar = QtWidgets.QProgressBar()
            self.progress_bar.setRange(0, 1000)
            self.progress_label = QtWidgets.QLabel("Waiting.")
            layout.addWidget(self.progress_bar)
            layout.addWidget(self.progress_label)
            self.log_edit = QtWidgets.QPlainTextEdit()
            self.log_edit.setReadOnly(True)
            layout.addWidget(self.log_edit)
            self.tabs.addTab(page, "4. Calculation")

        def _log(self, text: str) -> None:
            self.log_edit.appendPlainText("[%s] %s" % (_dt.datetime.now().strftime("%H:%M:%S"), text))

        def _progress(self, fraction: float, text: str) -> None:
            self.progress_bar.setValue(int(1000 * max(0.0, min(1.0, fraction))))
            self.progress_label.setText(text)

        def _start_worker(self, job: Callable, on_done: Callable) -> None:
            self.worker = Worker(job, self)
            self.worker.setStackSize(64 * 1024 * 1024)
            self.worker.log.connect(self._log)
            self.worker.progress.connect(self._progress)
            self.worker.done.connect(lambda value: (self._worker_finished(), on_done(value)))
            self.worker.failed.connect(self._worker_failed)
            self.worker.start()

        def _worker_finished(self) -> None:
            self.worker = None
            self.runner = None
            for button in (self.run_button, self.quick_button, self.geometry_button):
                button.setEnabled(True)
            self.stop_button.setEnabled(False)

        def _worker_failed(self, message: str) -> None:
            self._worker_finished()
            self._log("ERROR: " + message)
            self.progress_label.setText("Failed — see the log.")
            self._error(message[-3000:])

        def _stop(self) -> None:
            if self.runner is not None:
                self._log("Stop requested…")
                self.runner.cancel()

        def _start_analysis(self, with_simulation: bool) -> None:
            """Starts (in a separate thread) geometry, RASPA3 simulation, analysis and exports."""
            if self.worker is not None:
                return
            settings = self._collect_settings()
            self.qsettings.setValue("raspa_path", settings.raspa_path)
            self.qsettings.setValue("output_base", self.out_edit.text().strip())
            if with_simulation:
                problems = RaspaInputBuilder(self.structure, settings).validate()
                if not os.path.isfile(settings.raspa_path):
                    problems.append("RASPA3 executable not found: '%s'." % settings.raspa_path)
                if problems:
                    self._error("Cannot start the simulation:\n - " + "\n - ".join(problems))
                    return
            elif self.experimental is None:
                self._error("Load experimental data first (tab 2).")
                return
            structure, experimental = self.structure, self.experimental
            geometric, geometry_mode = self.geometric, self._geometry_mode()
            runner = RaspaRunner(structure, settings) if with_simulation else None
            self.runner = runner

            def job(log, progress):
                results = AnalysisResults(structure, settings, experimental=experimental)
                if runner is not None:
                    runner.log, runner.progress = log, progress
                    results.simulated = runner.run()
                progress(1.0, "Analysing the results…")
                if geometric is not None:
                    results.geometric = geometric
                else:
                    analyser = GeometricAnalyzer(structure, ADSORBATES[settings.gas]["probe_diameter"], log)
                    results.geometric = analyser.run(geometry_mode)
                results.analyse()
                written = ReportBuilder(results, settings.work_dir, log).write_all()
                return results, written

            for button in (self.run_button, self.quick_button, self.geometry_button):
                button.setEnabled(False)
            self.stop_button.setEnabled(with_simulation)
            self.progress_bar.setValue(0)
            self._log("=== New calculation: %s ===" % ("RASPA3 simulation" if with_simulation else "without simulation"))
            self._start_worker(job, self._analysis_done)

        def _analysis_done(self, payload) -> None:
            results, written = payload
            self.results = results
            self._show_geometry(results.geometric)
            self.last_html = written.get("html", "")
            self._publish_to_mercury()
            self._progress(1.0, "Done. Results: %s" % results.settings.work_dir)
            self._show_results()
            self.tabs.setCurrentIndex(4)

        def _publish_to_mercury(self) -> None:
            """Copies the HTML report into the file Mercury displays when the script ends."""
            target = self.bridge.output_html_file() if self.bridge else ""
            if target and self.last_html and os.path.isfile(self.last_html):
                try:
                    shutil.copyfile(self.last_html, target)
                    self._log("HTML report passed to Mercury: %s" % target)
                except OSError as exc:
                    self._log("HTML report not passed to Mercury: %s" % exc)

        # ================================ Tab 5: Results ================================
        def _build_results_tab(self) -> None:
            page = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout(page)
            self.result_tabs = QtWidgets.QTabWidget()
            self.figure_tabs = {}
            for key, title in (("isotherm", "Figure 1 — Isotherms"), ("bet_exp", "Figure 2 — Experimental BET"),
                               ("bet_sim", "Figure 3 — Simulated BET"), ("surfaces", "Figure 4 — Surface areas")):
                tab = FigureTab()
                self.figure_tabs[key] = tab
                self.result_tabs.addTab(tab, title)
            self.comparison_view = QtWidgets.QTextBrowser()
            self.result_tabs.addTab(self.comparison_view, "Comparison table")
            layout.addWidget(self.result_tabs)
            row = QtWidgets.QHBoxLayout()
            export = QtWidgets.QPushButton("Export again (Excel, CSV, PDF)")
            export.clicked.connect(self._export_again)
            folder = QtWidgets.QPushButton("Open the results folder")
            folder.clicked.connect(self._open_folder)
            row.addStretch(1)
            row.addWidget(export)
            row.addWidget(folder)
            layout.addLayout(row)
            self.tabs.addTab(page, "5. Results")

        def _show_results(self) -> None:
            factory = FigureFactory(self.results)
            for key, figure in factory.all_figures().items():
                self.figure_tabs[key].set_figure(figure)
            notes = ReportBuilder(self.results, self.results.settings.work_dir)._notes()
            self.comparison_view.setHtml(
                self._html_table(self.results.comparison_frame()) + "<br/>" +
                self._html_table(self.results.bet_frame().astype(str)) +
                "<p>" + "<br/>".join(n.replace("\n", " ") for n in notes) + "</p>")

        def _export_again(self) -> None:
            if self.results is None:
                self._error("No results to export.")
                return
            try:
                self.results.settings.cross_section = self._float(self.cross_edit, self.results.settings.cross_section)
                self.results.settings.gurvich_p_rel = self._float(self.gurvich_edit, 0.95)
                self.results.analyse()                          # takes a changed cross-section / P/P0 into account
                written = ReportBuilder(self.results, self.results.settings.work_dir, self._log).write_all()
                self.last_html = written.get("html", "")
                self._publish_to_mercury()
                self._show_results()
            except Exception as exc:
                self._error("Export failed:\n%s" % exc)

        def _open_folder(self) -> None:
            folder = self.results.settings.work_dir if self.results else self.out_edit.text()
            if os.path.isdir(folder):
                QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(folder))

        def closeEvent(self, event) -> None:                    # noqa: N802 (name imposed by Qt)
            """On closing: stops RASPA3 and leaves a report for Mercury."""
            if self.runner is not None:
                self.runner.cancel()
            if self.worker is not None:
                self.worker.wait(8000)
            target = self.bridge.output_html_file() if self.bridge else ""
            if target and not self.last_html:
                try:
                    with open(target, "w", encoding="utf-8") as handle:
                        handle.write("<html><body><h2>%s</h2><p>Session closed without a calculation.</p></body></html>"
                                     % APP_NAME)
                except OSError:
                    pass
            event.accept()


# =============================================================================
# 14. BATCH MODE (no GUI) AND ENTRY POINT
# =============================================================================
def run_batch(structure: StructureData, args: argparse.Namespace) -> AnalysisResults:
    """Complete pipeline without GUI (tests, scripts, remote calculation)."""
    gas = args.gas
    structure.build_p1(resolve_disorder=True, remove_guests=not args.keep_guests)
    for line in structure.p1_log:
        print(line)
    settings = SimulationSettings(
        gas=gas, temperature=args.temperature or ADSORBATES[gas]["temperature"],
        p0=args.p0 or ADSORBATES[gas]["p0"], n_init_cycles=args.init, n_prod_cycles=args.prod,
        supercell=tuple(args.cells) if args.cells else structure.minimal_supercell(12.0),
        charge_mode=CHARGE_MODES[args.charges], n_points=args.points,
        p_rel_min=args.pmin, p_rel_max=args.pmax, cross_section=ADSORBATES[gas]["cross_section"],
        raspa_path=args.raspa or find_raspa_executable(), n_jobs=args.jobs,
        work_dir=args.out or default_work_dir(structure, gas))
    results = AnalysisResults(structure, settings)
    if args.exp:
        results.experimental = ExperimentalDataLoader.load_file(args.exp)
    if not args.no_sim:
        runner = RaspaRunner(structure, settings, log=print,
                             progress=lambda f, t: print("  %5.1f %%  %s" % (100 * f, t), end="\r"))
        results.simulated = runner.run()
        print()
    results.geometric = GeometricAnalyzer(structure, ADSORBATES[gas]["probe_diameter"], print).run(
        "internal" if args.internal_geometry else "auto")
    results.analyse()
    ReportBuilder(results, settings.work_dir, print).write_all()
    print(results.comparison_frame().to_string(index=False))
    return results


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=APP_NAME, add_help=True)
    parser.add_argument("cif", nargs="?", default="", help="CIF file (outside Mercury)")
    parser.add_argument("--batch", action="store_true", help="no graphical interface")
    parser.add_argument("--exp", default="", help="experimental isotherm (xlsx, xls, csv)")
    parser.add_argument("--gas", default="N2", choices=list(ADSORBATES))
    parser.add_argument("--raspa", default="", help="path of the raspa3 executable")
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--p0", type=float, default=0.0, help="saturation pressure (Pa)")
    parser.add_argument("--init", type=int, default=2000, help="equilibration cycles")
    parser.add_argument("--prod", type=int, default=5000, help="production cycles")
    parser.add_argument("--points", type=int, default=35)
    parser.add_argument("--pmin", type=float, default=1e-6)
    parser.add_argument("--pmax", type=float, default=1.0)
    parser.add_argument("--cells", type=int, nargs=3, default=None, help="supercell a b c")
    parser.add_argument("--charges", type=int, default=0, choices=[0, 1, 2],
                        help="0 charge equilibration, 1 charges from the CIF, 2 none")
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--out", default="", help="results folder")
    parser.add_argument("--no-sim", action="store_true", help="analysis without RASPA3 simulation")
    parser.add_argument("--keep-guests", action="store_true", help="keep guest species")
    parser.add_argument("--internal-geometry", action="store_true",
                        help="force the internal geometric calculation")
    # --- internal options (communication between the script's processes) ---
    parser.add_argument("--gui-child", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--child-meta", default="", help=argparse.SUPPRESS)
    parser.add_argument("--pore-helper", default="", help=argparse.SUPPRESS)
    return parser


def run_gui_child(args: argparse.Namespace) -> int:
    """GUI process: Qt only; the structure comes from the transfer CIF."""
    meta = {}
    if args.child_meta and os.path.isfile(args.child_meta):
        with open(args.child_meta, "r", encoding="utf-8") as handle:
            meta = json.load(handle)
    if not QT_BINDING:
        print("No usable Qt binding in the GUI process.", file=sys.stderr)
        return 1
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    try:
        structure = StructureData.from_cif_file(args.cif)
    except Exception as exc:
        QtWidgets.QMessageBox.critical(None, APP_NAME, "Unreadable structure:\n%s" % exc)
        return 1
    structure.name = meta.get("name") or structure.name
    structure.source_file = meta.get("source_file", "")
    structure.csd_cif_path = meta.get("csd_cif_path", "")
    structure.notes = list(meta.get("notes", [])) + structure.notes
    bridge = ChildBridge(meta.get("html_out", ""), bool(meta.get("from_mercury")))
    window = MainWindow(structure, bridge, args)
    window.show()
    window.raise_()
    window.activateWindow()
    close_ms = os.environ.get("RASPA_POROSITY_TEST_CLOSE_MS")       # automated tests only
    if close_ms:
        QtCore.QTimer.singleShot(int(close_ms), window.close)
    run_loop = getattr(app, "exec", None) or getattr(app, "exec_")
    run_loop()
    return 0


def launch_gui_child(structure: StructureData, bridge: MercuryBridge, args: argparse.Namespace) -> int:
    """Parent process: writes the transfer CIF and starts the GUI in another process."""
    temp_dir = tempfile.mkdtemp(prefix="raspa_gui_")
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", structure.name)[:60] or "structure"
    transfer = os.path.join(temp_dir, safe + ".cif")
    structure.write_transfer_cif(transfer)
    meta_path = os.path.join(temp_dir, "meta.json")
    with open(meta_path, "w", encoding="utf-8") as handle:
        json.dump({"name": structure.name, "source_file": structure.source_file,
                   "csd_cif_path": structure.csd_cif_path, "notes": structure.notes,
                   "html_out": bridge.output_html_file(), "from_mercury": bridge.available}, handle)
    command = [sys.executable, os.path.abspath(__file__), transfer, "--gui-child", "--child-meta", meta_path,
               "--gas", args.gas]
    for option, value in (("--exp", args.exp), ("--raspa", args.raspa)):
        if value:
            command += [option, value]
    # The GUI must not inherit the Qt paths of Mercury or of the CSD API:
    # PyQt/PySide must load their own libraries and plugins.
    env = {k: v for k, v in os.environ.items()
           if not (k.startswith("DYLD_") or k in ("QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH",
                                                  "QML2_IMPORT_PATH", "QML_IMPORT_PATH"))}
    # The GUI runs detached from Mercury (own session); its output goes to a log file.
    log = open(os.path.join(temp_dir, "gui_child.log"), "w")
    subprocess.Popen(command, env=env, stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                     start_new_session=True)
    target = bridge.output_html_file()
    if target:
        with open(target, "w", encoding="utf-8") as handle:
            handle.write("<html><body><h2>%s</h2><p>The interface was opened in a separate "
                         "window; the results will be shown there.</p></body></html>" % APP_NAME)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    # parse_known_args: Mercury adds its own arguments, which must be tolerated
    args, _unknown = build_argument_parser().parse_known_args(argv)
    cif_path = args.cif if args.cif and os.path.isfile(args.cif) and args.cif.lower().endswith(".cif") else ""
    if args.pore_helper:                                        # helper process (CSD API only)
        return pore_helper_main(args.pore_helper)
    if args.gui_child:                                          # GUI process (Qt only)
        return run_gui_child(args)
    # ---- parent process: CSD Python API / Mercury, never Qt ----
    try:
        structure, bridge, _temp_cif = load_structure(cif_path)
    except Exception as exc:
        text = str(exc)
        if "Symbol not found" in text or "dlopen" in text:      # Qt / CSD Python API conflict
            print("ERROR (%s %s): the CSD Python API no longer loads, because a Qt library "
                  "installed in Mercury's Python conflicts with its own.\n\n%s\n\n"
                  "Technical detail: %s" % (APP_NAME, APP_VERSION, qt_fix_instructions(), text[:400]),
                  file=sys.stderr)
        else:
            print("ERROR (%s %s): %s" % (APP_NAME, APP_VERSION, text), file=sys.stderr)
        return 1
    if args.batch:
        try:
            run_batch(structure, args)
        except (RaspaError, ValueError, RuntimeError) as exc:
            print("\nERROR: %s" % exc, file=sys.stderr)
            return 1
        return 0
    if not detect_qt_binding():
        message = "No Qt binding is available for the graphical interface.\n\n" + qt_fix_instructions()
        print(message, file=sys.stderr)
        target = bridge.output_html_file()
        if target:
            with open(target, "w", encoding="utf-8") as handle:
                handle.write("<html><body><h2>%s</h2><pre>%s</pre></body></html>" % (APP_NAME, message))
        return 1
    return launch_gui_child(structure, bridge, args)


if __name__ == "__main__":
    _code = main()
    sys.stdout.flush()
    sys.stderr.flush()
    if "ccdc" in sys.modules:
        # On macOS, the Qt libraries of the CSD Python API can crash
        # ("QMutex: destroying locked mutex", segmentation fault) during the
        # normal interpreter shutdown, although all the work is done.
        # Exit immediately, without the clean-up phase, to avoid this false crash.
        os._exit(int(_code or 0))
    sys.exit(_code)
