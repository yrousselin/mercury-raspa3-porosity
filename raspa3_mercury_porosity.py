#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
===============================================================================
 RASPA3 Porosity Analyser pour Mercury (CCDC)
===============================================================================
 Comparaison automatique des propriétés poreuses EXPÉRIMENTALES et SIMULÉES
 d'un matériau poreux (MOF, ZIF, COF...) :

   1. récupère la structure ouverte dans Mercury (CSD Python API) ;
   2. l'exporte vers un CIF temporaire, puis vers un CIF P1 « nettoyé » ;
   3. prépare les fichiers d'entrée RASPA3 (simulation.json, force_field.json,
      <molécule>.json, framework.cif) ;
   4. lance RASPA3 (GCMC) pour chaque point de pression de l'isotherme ;
   5. lit les résultats, calcule la surface BET (critères de Rouquerol) ;
   6. calcule les propriétés géométriques (surface accessible, volume poreux,
      PLD, LCD) via le Pore Analyser de la CSD Python API ou, à défaut, via une
      implémentation interne documentée ;
   7. compare aux données expérimentales, trace les figures ;
   8. exporte Excel, CSV, PDF et un rapport HTML affiché dans Mercury.

 Adsorbats disponibles : N2 (77 K), Ar (87 K), CO2 (273 K).

-------------------------------------------------------------------------------
 INSTALLATION DANS MERCURY
-------------------------------------------------------------------------------
   * Copier ce fichier dans un dossier de scripts, puis dans Mercury :
       CSD Python API > Options... > Add Location  (ajouter ce dossier)
     Le script apparaît alors dans le menu « CSD Python API »
     (menu « Calculate > Python Scripts » selon la version de Mercury).
   * Modules requis dans l'environnement Python utilisé par Mercury :
       numpy, pandas, matplotlib, openpyxl, et une liaison Qt
       (PyQt6, PyQt5, PySide6 ou PySide2 — la première trouvée est utilisée).
       Optionnels : scipy (calcul géométrique plus rapide), xlrd (fichiers .xls).
     Exemple (terminal, avec le python de la CSD Python API) :
       python -m pip install --target ~/RASPA_Mercury/qt_libs PyQt6
       python -m pip install matplotlib openpyxl xlrd scipy
     ATTENTION : ne pas installer PyQt directement dans le Python de Mercury
     (conflit avec les bibliothèques Qt de la CSD Python API sur macOS) ; le
     script affiche les commandes exactes si quelque chose manque.
   * RASPA3 : installer le paquet macOS (.pkg) ou « conda install raspa3 ».
     Le chemin de l'exécutable se règle dans l'onglet « Paramètres RASPA »
     (il est mémorisé d'une session à l'autre).

-------------------------------------------------------------------------------
 UTILISATION HORS MERCURY (test, traitement par lot)
-------------------------------------------------------------------------------
   python raspa3_mercury_porosity.py structure.cif            # interface
   python raspa3_mercury_porosity.py structure.cif --batch \
          --exp mesure.xls --gas N2 --raspa /usr/local/bin/raspa3
   Dans ce mode, le CIF est lu par l'analyseur interne (classe CIFParser) :
   c'est l'alternative documentée quand la CSD Python API n'est pas disponible.

-------------------------------------------------------------------------------
 CHOIX DE MODÉLISATION (à connaître pour interpréter les résultats)
-------------------------------------------------------------------------------
   * Charpente rigide ; Lennard-Jones UFF ou DREIDING par ÉLÉMENT, règles de
     mélange de Lorentz-Berthelot, potentiel tronqué-décalé à 12 Å.
   * Adsorbats : N2 et CO2 TraPPE (3 sites, charges), Ar Lennard-Jones 1 site.
   * Charges de la charpente : équilibration de charges de RASPA3
     (Wilmer-Snurr), charges du CIF (_atom_site_charge) ou aucune.
   * « Cycles d'équilibration » de l'interface = NumberOfInitializationCycles
     de RASPA3 (la phase « Equilibration » de RASPA3 ne sert qu'au CFCMC).
   * La simulation donne l'adsorption ABSOLUE ; l'expérience mesure l'excès.
     À 77/87 K et P <= P0 la différence est négligeable.
   * Désordre cristallographique : les sites d'occupation partielle en conflit
     stérique sont résolus en gardant une seule composante ; les espèces
     invitées non liées à la charpente (solvant) peuvent être supprimées.
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

import numpy as np
import pandas as pd

APP_NAME = "RASPA3 Porosity Analyser"
APP_VERSION = "1.3"

# =============================================================================
# 1. CONSTANTES PHYSIQUES ET TABLES
# =============================================================================
N_AVOGADRO = 6.02214076e23          # mol-1
V_MOLAR_STP = 22414.0               # cm3(STP)/mol (0 °C, 1 atm)
KCAL_TO_K = 503.2195                # 1 kcal/mol exprimé en kelvins (E/kB)

# Masses atomiques (g/mol)
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

# Rayons covalents (Å, Cordero 2008) : détection des liaisons / du désordre
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

# UFF (Rappé et al., JACS 1992, 114, 10024) : (x_i [Å], D_i [kcal/mol]).
# Conversion Lennard-Jones : sigma = x_i / 2^(1/6) ; epsilon/kB = D_i * 503.22
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
# DREIDING (Mayo et al., J. Phys. Chem. 1990, 94, 8897) : (R0 [Å], D0 [kcal/mol])
_DREIDING_RAW = {
    "H": (3.195, 0.0152), "B": (4.02, 0.095), "C": (3.8983, 0.0951),
    "N": (3.6621, 0.0774), "O": (3.4046, 0.0957), "F": (3.4720, 0.0725),
    "Al": (4.39, 0.31), "Si": (4.27, 0.31), "P": (4.15, 0.32),
    "S": (4.03, 0.344), "Cl": (3.9503, 0.2833), "Br": (3.95, 0.37),
}


def _lj_from_raw(raw: Tuple[float, float]) -> Tuple[float, float]:
    """(distance du minimum, profondeur kcal/mol) -> (epsilon/kB [K], sigma [Å])."""
    return raw[1] * KCAL_TO_K, raw[0] / 2.0 ** (1.0 / 6.0)


def framework_lj(element: str, forcefield: str) -> Tuple[float, float, str]:
    """Paramètres LJ (epsilon/kB [K], sigma [Å], source) d'un élément de charpente.

    `forcefield` contient "DREIDING" ou "UFF". DREIDING ne couvrant que
    quelques éléments, les autres (métaux notamment) sont pris dans UFF :
    c'est la pratique usuelle pour les MOFs.
    """
    if "DREIDING" in forcefield.upper() and element in _DREIDING_RAW:
        eps, sig = _lj_from_raw(_DREIDING_RAW[element])
        return eps, sig, "DREIDING"
    if element in _UFF_RAW:
        eps, sig = _lj_from_raw(_UFF_RAW[element])
        return eps, sig, "UFF"
    raise ValueError("Aucun paramètre Lennard-Jones pour l'élément '%s'." % element)


FORCE_FIELDS = ["TraPPE (adsorbat) + UFF (charpente)",
                "TraPPE (adsorbat) + DREIDING/UFF (charpente)"]

# -----------------------------------------------------------------------------
# Adsorbats. Chaque entrée décrit : conditions par défaut, constantes critiques
# (RASPA3 en déduit le coefficient de fugacité par Peng-Robinson), modèle
# moléculaire (sites, LJ, charges), section d'encombrement BET, volume molaire
# du liquide (règle de Gurvich) et diamètre de sonde géométrique.
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

CHARGE_MODES = ["Équilibration de charges (RASPA3, Wilmer-Snurr)",
                "Charges du CIF (_atom_site_charge)",
                "Aucune charge sur la charpente"]

# Couleurs CPK simplifiées pour l'image de la structure
ELEMENT_COLORS = {"H": "#e8e8e8", "C": "#505050", "N": "#3050f8", "O": "#ff2010",
                  "F": "#90e050", "Cl": "#1ff01f", "Br": "#a62929", "S": "#e6c820",
                  "P": "#ff8000", "Zn": "#7d80b0", "Cu": "#c88033", "Zr": "#94e0e0",
                  "Co": "#f090a0", "Fe": "#e06633", "Ni": "#50d050", "Al": "#bfa6a6",
                  "Cr": "#8a99c7", "Mg": "#8aff00", "Si": "#f0c8a0", "Ti": "#bfc2c7"}

# Palette des figures (lisible en niveaux de gris et pour les daltoniens)
COLOR_EXP = "#1f4e9c"
COLOR_SIM = "#d1495b"
COLOR_GEO = "#6b8e23"


# =============================================================================
# 2. OUTILS GÉNÉRAUX
# =============================================================================
def to_float(text) -> float:
    """Convertit '0.03150(13)', '1,5', '?' ... en float (NaN si impossible)."""
    if text is None:
        return float("nan")
    if isinstance(text, (int, float, np.floating, np.integer)):
        return float(text)
    s = str(text).strip().replace("−", "-")
    s = re.sub(r"\(.*?\)", "", s)                  # retire l'écart-type (13)
    if "," in s and "." not in s:
        s = s.replace(",", ".")                    # virgule décimale
    try:
        return float(s)
    except ValueError:
        return float("nan")


def element_from_symbol(symbol: str) -> str:
    """'Zn2+', 'ZN', 'C1A' -> 'Zn', 'C'. Renvoie '' si non reconnu."""
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
    """Formule brute en notation de Hill (C, H puis ordre alphabétique)."""
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
    """Écart relatif (%) de `value` par rapport à `reference`."""
    if reference is None or value is None:
        return float("nan")
    if not np.isfinite(reference) or not np.isfinite(value) or reference == 0:
        return float("nan")
    return 100.0 * (value - reference) / reference


def fmt(value, digits: int = 2, unit: str = "") -> str:
    """Formatage tolérant (None/NaN -> '—')."""
    try:
        if value is None or not np.isfinite(float(value)):
            return "—"
        return ("%." + str(digits) + "f") % float(value) + ((" " + unit) if unit else "")
    except (TypeError, ValueError):
        return str(value)


# =============================================================================
# 3. LECTURE CIF INTERNE (alternative documentée à la CSD Python API)
# =============================================================================
class CIFParser:
    """Analyseur CIF minimal mais robuste (CIF1 : SHELX, Olex2, CSD, RASPA...).

    Il ne remplace pas la CSD Python API : il sert (a) hors Mercury, (b) à
    relire le fichier d'origine pour disposer des occupations et groupes de
    désordre, informations indispensables pour construire un modèle P1 propre.
    Seul le premier bloc `data_` contenant des atomes est exploité.
    """

    def __init__(self, text: str):
        self.items: Dict[str, str] = {}
        self.loops: List[Tuple[List[str], List[List[str]]]] = []
        self._parse(text)

    # -- découpage en jetons ---------------------------------------------------
    @staticmethod
    def _tokenize(text: str) -> List[Tuple[str, str]]:
        """Renvoie des couples (type, valeur) ; type = 'tok' ou 'str'."""
        tokens: List[Tuple[str, str]] = []
        lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.startswith(";"):                       # champ multi-lignes
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
                    while True:                            # guillemet fermant = suivi d'un blanc
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
                    break                                   # premier bloc utile lu
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

    # -- accès ----------------------------------------------------------------
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
# 4. MODÈLE DE STRUCTURE
# =============================================================================
@dataclass
class Atom:
    """Atome de l'unité asymétrique (ou de la maille P1 après expansion)."""
    label: str
    element: str
    frac: Tuple[float, float, float]
    occupancy: float = 1.0
    disorder_group: str = "."
    charge: Optional[float] = None


def parse_symop(op: str) -> Tuple[np.ndarray, np.ndarray]:
    """'-x+1/2, y, z+0.5' -> (matrice de rotation 3x3, translation)."""
    rot = np.zeros((3, 3))
    trans = np.zeros(3)
    parts = op.lower().replace(" ", "").split(",")
    if len(parts) != 3:
        raise ValueError("Opérateur de symétrie illisible : %r" % op)
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
    """Structure cristalline indépendante de toute bibliothèque externe.

    Elle peut être construite depuis une entrée CSD Python API
    (`from_csd_entry`) ou depuis un fichier CIF (`from_cif_file`).
    """

    def __init__(self):
        self.name = "structure"
        self.source_file = ""
        self.cell = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)      # a b c alpha beta gamma
        self.spacegroup = "P 1"
        self.symops: List[str] = ["x,y,z"]
        self.atoms: List[Atom] = []
        self.formula_reported = ""                         # formule du CIF / de la CSD
        self.formula_weight_reported = float("nan")
        self.z_value = float("nan")
        self.density_reported = float("nan")
        self.notes: List[str] = []                         # avertissements
        self.csd_cif_path = ""                             # CIF lisible par la CSD Python API (Pore Analyser)
        # --- résultat de l'expansion P1 (build_p1) ---
        self.p1_elements: List[str] = []
        self.p1_frac = np.zeros((0, 3))
        self.p1_charges: Optional[np.ndarray] = None
        self.p1_removed: Dict[str, int] = {}
        self.p1_log: List[str] = []

    # ------------------------------------------------------------------ maille
    @property
    def cell_matrix(self) -> np.ndarray:
        """Matrice 3x3 dont les LIGNES sont les vecteurs a, b, c (Å)."""
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
        """Distances entre faces opposées de la maille (Å)."""
        m = self.cell_matrix
        v = self.volume
        return np.array([v / np.linalg.norm(np.cross(m[1], m[2])),
                         v / np.linalg.norm(np.cross(m[2], m[0])),
                         v / np.linalg.norm(np.cross(m[0], m[1]))])

    def minimal_supercell(self, cutoff: float = 12.0) -> Tuple[int, int, int]:
        """Plus petite supercellule respectant la convention d'image minimale."""
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
            raise ValueError("Paramètres de maille absents du CIF.") from exc
        if not all(np.isfinite(st.cell)):
            raise ValueError("Paramètres de maille absents ou illisibles dans le CIF.")
        st.spacegroup = (cif.get("_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m",
                                 "_space_group_name_hall", default="?") or "?").strip()
        ops = cif.symmetry_operators()
        if ops:
            st.symops = ops
        else:
            st.notes.append("Aucun opérateur de symétrie dans le CIF : structure traitée en P1.")
        st.formula_reported = (cif.get("_chemical_formula_sum", default="") or "").strip()
        st.formula_weight_reported = to_float(cif.get("_chemical_formula_weight"))
        st.z_value = to_float(cif.get("_cell_formula_units_z"))
        st.density_reported = to_float(cif.get("_exptl_crystal_density_diffrn"))
        found = cif.loop_with("_atom_site_fract_x")
        if not found:
            raise ValueError("Aucune boucle _atom_site_fract_x dans le CIF.")
        tags, rows = found
        col = {t: i for i, t in enumerate(tags)}

        def cell_value(row, tag, default=None):
            return row[col[tag]] if tag in col else default

        for row in rows:
            label = cell_value(row, "_atom_site_label", "X")
            element = element_from_symbol(cell_value(row, "_atom_site_type_symbol", "") or "") \
                or element_from_symbol(label)
            if not element:
                st.notes.append("Atome ignoré (élément inconnu) : %s" % label)
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
            raise ValueError("Aucun atome exploitable dans le CIF.")
        return st

    @classmethod
    def from_csd_entry(cls, entry, original_file: str = "") -> "StructureData":
        """Construit la structure depuis un objet `ccdc.entry.Entry`.

        Stratégie (du plus fidèle au moins fidèle) :
          1. si Mercury a chargé un fichier CIF dont le chemin est connu, les
             atomes sont relus dans ce fichier (occupations, désordre) ;
          2. sinon, atomes de `crystal.asymmetric_unit_molecule` ;
          3. sinon, relecture du CIF exporté par `crystal.to_string('cif')`.
        Les métadonnées (formule, groupe d'espace...) viennent de l'API.
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
                    raise ValueError("aucun atome avec coordonnées")
            except Exception:
                st = cls.from_cif_text(crystal.to_string("cif"))
        # Métadonnées issues de l'API (chaque accès est protégé : les attributs
        # disponibles varient légèrement selon la version de la CSD Python API).
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

    # ------------------------------------------------------------ descripteurs
    def asymmetric_summary(self) -> Dict[str, float]:
        """Contenu de la maille pondéré par les occupations (avant nettoyage)."""
        elements, frac, occ, _ = self._expand(dedupe=True)
        counts: Dict[str, float] = {}
        for el, o in zip(elements, occ):
            counts[el] = counts.get(el, 0.0) + o
        return counts

    @staticmethod
    def mass_of(counts: Dict[str, float]) -> float:
        return float(sum(ATOMIC_MASS[el] * n for el, n in counts.items()))

    def density_from_mass(self, mass: float) -> float:
        """Masse de la maille (g/mol) -> masse volumique (g/cm3)."""
        return mass / N_AVOGADRO / (self.volume * 1e-24)

    # ---------------------------------------------------------------- expansion
    def _expand(self, dedupe: bool = True):
        """Applique tous les opérateurs de symétrie.

        Renvoie (éléments, coordonnées fractionnaires dans [0,1[, occupations,
        index de l'atome parent). Les doublons créés par les positions
        spéciales (même atome parent à moins de 0,2 Å) sont fusionnés.
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
        """Paires (i<j) d'atomes distants de moins de `rmax` (image minimale).

        Calcul par blocs en numpy pur (aucune dépendance), valable pour les
        courtes distances utilisées ici (liaisons, conflits de désordre).
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
        """Construit le modèle P1 utilisé par RASPA3 et les calculs géométriques.

        resolve_disorder : parmi les sites d'occupation partielle en conflit
            stérique (deux composantes d'un même groupement désordonné), une
            seule position est conservée (algorithme glouton, occupation
            décroissante). Exemple : méthyle de ZIF-8 décrit par 6 demi-H -> 3 H.
        remove_guests : supprime les fragments non liés à la charpente dont la
            taille est inférieure à `guest_fraction` fois celle du plus grand
            fragment (molécules de solvant, oxygènes d'eau pseudo-localisés...).
        """
        elements, frac, occ, parents = self._expand(dedupe=True)
        n0 = len(elements)
        keep = np.ones(n0, dtype=bool)
        log: List[str] = ["Maille P1 brute : %d sites" % n0]
        removed: Dict[str, int] = {}
        radii = np.array([COVALENT_RADIUS.get(e, 1.5) for e in elements])

        # --- 1) résolution du désordre ------------------------------------
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
            accepted = ~partial                              # sites pleins : toujours gardés
            for k in order:
                if any(accepted[m] for m in conflicts.get(int(k), [])):
                    keep[k] = False
                    removed["désordre:" + self.atoms[parents[k]].label] = \
                        removed.get("désordre:" + self.atoms[parents[k]].label, 0) + 1
                else:
                    accepted[k] = True
            log.append("Désordre : %d sites partiels, %d positions écartées"
                       % (int(partial.sum()), int((~keep).sum())))
        elif partial.any():
            log.append("ATTENTION : %d sites d'occupation partielle conservés tels quels."
                       % int(partial.sum()))

        # --- 2) suppression des espèces invitées --------------------------
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
                    key = "invité:" + self.atoms[parents[k]].label
                    removed[key] = removed.get(key, 0) + 1
            log.append("Espèces invitées : %d atomes retirés (%d fragments)"
                       % (n_guest, len(small_roots)))

        self.p1_elements = [e for e, k in zip(elements, keep) if k]
        self.p1_frac = frac[keep]
        charges = [self.atoms[p].charge for p in parents[keep]]
        self.p1_charges = np.array(charges, dtype=float) if charges and all(c is not None for c in charges) else None
        self.p1_removed = removed
        log.append("Modèle P1 final : %d atomes, %s" % (len(self.p1_elements),
                                                         hill_formula(self.p1_counts())))
        self.p1_log = log

    def removed_summary(self) -> str:
        """Résumé lisible des atomes écartés, par motif et par élément."""
        groups: Dict[Tuple[str, str], int] = {}
        for key, count in self.p1_removed.items():
            reason, label = key.split(":", 1)
            groups[(reason, element_from_symbol(label))] = groups.get((reason, element_from_symbol(label)), 0) + count
        return " ; ".join("%s : %d %s" % (reason, count, el) for (reason, el), count in sorted(groups.items())) \
            or "aucun"

    def p1_counts(self) -> Dict[str, float]:
        counts: Dict[str, float] = {}
        for el in self.p1_elements:
            counts[el] = counts.get(el, 0.0) + 1.0
        return counts

    @property
    def p1_mass(self) -> float:
        """Masse molaire de la maille P1 nettoyée (g/mol de maille)."""
        return self.mass_of(self.p1_counts())

    @property
    def p1_cart(self) -> np.ndarray:
        return self.p1_frac @ self.cell_matrix

    def write_p1_cif(self, path: str, with_charges: bool = False) -> None:
        """Écrit le CIF P1 lu par RASPA3 (étiquette = type = symbole de l'élément)."""
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
        """Écrit la structure complète (symétrie, occupations, désordre, charges).

        Ce CIF sert à transmettre la structure du processus parent (CSD Python
        API) au processus d'interface ; il est relu par `from_cif_file`.
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
# 5. PROPRIÉTÉS GÉOMÉTRIQUES (surface accessible, volume poreux, PLD, LCD)
# =============================================================================
PORE_MARKER = "PORE_ANALYSER_JSON:"       # préfixe de la ligne de résultat du processus auxiliaire


def pore_helper_main(cif_path: str) -> int:
    """Processus auxiliaire : Pore Analyser de la CSD Python API sur un CIF.

    Imprime une ligne `PORE_ANALYSER_JSON:{...}` contenant les champs d'un
    `GeometricResult` (ou {"error": ...}). N'importe jamais Qt.
    """
    try:
        from ccdc import io as ccdc_io
        from ccdc.descriptors import CrystalDescriptors
        crystal = ccdc_io.CrystalReader(cif_path)[0]
        analyser = CrystalDescriptors.PoreAnalyser(crystal)
        mass = float(analyser.system_mass)                   # g/mol de maille
        notes = []
        try:
            probe = float(analyser.settings.n_probe_sigma)
        except Exception:
            probe = 3.314                                    # valeur par défaut CSD
        surface = float(analyser.network_accessible_surface_area_per_mass)
        if not surface > 0:
            surface = float(analyser.total_surface_area_per_mass)
            notes.append("Surface totale (surface accessible par le réseau nulle).")
        to_cm3_g = 1e-24 * N_AVOGADRO / mass
        data = {"probe_diameter": probe, "surface_m2_g": surface,
                "surface_m2_cm3": float(analyser.total_surface_area_per_volume),
                "geometric_volume_cm3_g": float(analyser.total_geometric_volume) * to_cm3_g,
                "pore_volume_cm3_g": float(analyser.total_helium_volume) * to_cm3_g,
                "pore_volume_label": "volume hélium (Pore Analyser)",
                "void_fraction": float(analyser.total_geometric_volume) / float(analyser.system_volume),
                "pld": float(analyser.pore_limiting_diameter), "lcd": float(analyser.max_pore_diameter),
                "notes": notes}
    except BaseException as exc:
        data = {"error": "%s: %s" % (type(exc).__name__, exc)}
    print(PORE_MARKER + json.dumps(data))
    return 0 if "error" not in data else 1


@dataclass
class GeometricResult:
    """Résultats géométriques, quelle que soit la méthode employée."""
    method: str = ""
    probe_diameter: float = float("nan")      # Å
    surface_m2_g: float = float("nan")        # surface accessible (m2/g)
    surface_m2_cm3: float = float("nan")      # surface accessible (m2/cm3)
    pore_volume_cm3_g: float = float("nan")   # volume poreux retenu pour la comparaison
    pore_volume_label: str = ""
    geometric_volume_cm3_g: float = float("nan")   # sonde ponctuelle
    void_fraction: float = float("nan")
    pld: float = float("nan")                 # diamètre de pore limitant (Å)
    lcd: float = float("nan")                 # diamètre maximal de cavité (Å)
    notes: List[str] = field(default_factory=list)


class GeometricAnalyzer:
    """Calcul géométrique sur le modèle P1 nettoyé.

    Deux voies :
      * `run_csd()` : Pore Analyser de la CSD Python API
        (ccdc.descriptors.CrystalDescriptors.PoreAnalyser), exécuté dans un
        processus auxiliaire et utilisé seulement si la structure n'a pas été
        modifiée par le nettoyage (sinon le cristal vu par la CSD contiendrait
        encore le solvant) ;
      * `run_internal()` : implémentation équivalente, documentée ci-dessous.

    Implémentation interne
    ----------------------
    Rayon d'un atome de charpente : sigma_UFF / 2 (comme Pore Analyser, Zeo++
    et RASPA). Sonde : sphère de diamètre `probe_diameter`.
      * Surface accessible (Düren et al., J. Phys. Chem. C 2007, 111, 15350) :
        des points sont répartis sur la sphère de rayon (sigma_i + d_sonde)/2
        centrée sur chaque atome ; la fraction de points n'entrant en
        collision avec aucun autre atome donne l'aire accessible.
      * Grille 3D de la distance D à la surface de l'atome le plus proche :
          - volume géométrique      = fraction des noeuds où D > 0 ;
          - volume occupable sonde  = dilatation, par une sphère de rayon
            r_sonde, de la région accessible au centre de la sonde
            (D >= r_sonde) — définition d'Ongari et al., Langmuir 2017 ;
          - LCD = 2 max(D) ;
          - PLD = 2 r_max, r_max étant le plus grand rayon pour lequel la
            région D >= r percole à travers la maille périodique
            (dichotomie + parcours en largeur avec suivi des images).
        La précision sur PLD/LCD est de l'ordre du pas de grille.
    """

    MARGIN = 14.0   # Å : épaisseur de la couche d'images périodiques

    def __init__(self, structure: StructureData, probe_diameter: float,
                 log: Optional[Callable[[str], None]] = None):
        self.st = structure
        self.probe_d = float(probe_diameter)
        self.log = log or (lambda msg: None)

    # ------------------------------------------------------------------ CSD API
    def run_csd(self) -> GeometricResult:
        """Récupère les valeurs du Pore Analyser (CSD Python API >= 3.0.15).

        Le calcul est fait dans un processus Python séparé (option interne
        --pore-helper) : la CSD Python API ne doit pas être chargée dans le
        processus de l'interface Qt (conflit de bibliothèques Qt).
        """
        command = [sys.executable, os.path.abspath(__file__), "--pore-helper", self.st.csd_cif_path]
        proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, timeout=3600)
        data = None
        for line in proc.stdout.splitlines():
            if line.startswith(PORE_MARKER):
                data = json.loads(line[len(PORE_MARKER):])
        if data is None:
            raise RuntimeError("pas de réponse du Pore Analyser : %s" % (proc.stderr.strip()[-300:] or "?"))
        if "error" in data:
            raise RuntimeError(data["error"])
        return GeometricResult(method="CSD Python API - Pore Analyser", **data)

    # ------------------------------------------------------------- outils internes
    def _atoms_with_images(self) -> Tuple[np.ndarray, np.ndarray]:
        """Atomes P1 + images périodiques dans une couche de MARGIN Å."""
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
        """min_j(|p - a_j| - R_j) pour chaque point (distance à la surface atomique).

        Avec scipy : arbre k-d (rapide). Sans scipy : force brute par blocs.
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
        """Vrai si la région `mask` (grille périodique) est connexe à l'infini.

        Principe : chaque noeud est repéré par l'image de la maille dans
        laquelle il est atteint ; rejoindre un noeud déjà connu avec une image
        différente prouve l'existence d'un canal traversant (26 voisins).
          * avec scipy : étiquetage des amas dans la maille (ndimage.label),
            puis union-find pondéré sur les liaisons traversant les faces ;
          * sans scipy : parcours en largeur en Python pur (plus lent).
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
                    continue                                  # chaque paire une seule fois
                other = np.roll(labels, shift=(-d[0], -d[1], -d[2]), axis=(0, 1, 2))
                crossing = np.zeros(shape, dtype=bool)
                wraps = []
                for axis in range(3):
                    w = (idx[axis] + d[axis]) // shape[axis]   # -1, 0 ou +1
                    wraps.append(w)
                    crossing |= w != 0
                sel = crossing & (labels > 0) & (other > 0)
                if sel.any():
                    rows = np.stack([labels[sel], other[sel], wraps[0][sel], wraps[1][sel],
                                     wraps[2][sel]], axis=1)
                    edges.update(map(tuple, np.unique(rows, axis=0).tolist()))
            parent = list(range(count + 1))
            offset = [(0, 0, 0)] * (count + 1)                 # image relative à la racine

            def find(x):
                path = []
                while parent[x] != x:
                    path.append(x)
                    x = parent[x]
                total = (0, 0, 0)
                for node in reversed(path):                    # compression de chemin
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

    # --------------------------------------------------------------- calcul interne
    def run_internal(self, samples_per_atom: int = 300, grid_spacing: float = 0.25) -> GeometricResult:
        st = self.st
        if not st.p1_elements:
            raise ValueError("Modèle P1 vide : construire la structure d'abord.")
        res = GeometricResult(method="Calcul interne (Monte-Carlo + grille)",
                              probe_diameter=self.probe_d)
        rp = self.probe_d / 2.0
        cellm = st.cell_matrix
        atoms, radii = self._atoms_with_images()
        own_radii = np.array([framework_lj(e, "UFF")[1] / 2.0 for e in st.p1_elements])
        n_atoms = len(st.p1_elements)
        fast = self._have_scipy()
        if not fast:
            res.notes.append("scipy absent : calcul par force brute (échantillonnage réduit).")
            samples_per_atom = int(np.clip(1.5e9 / (n_atoms * len(atoms)), 30, samples_per_atom))
            grid_spacing = max(grid_spacing, (st.volume / max(1.5e9 / len(atoms), 1.0)) ** (1 / 3))

        # --- surface accessible : sphère de Fibonacci autour de chaque atome ---
        self.log("Surface accessible : %d points par atome..." % samples_per_atom)
        k = np.arange(samples_per_atom) + 0.5
        phi = np.arccos(1.0 - 2.0 * k / samples_per_atom)
        theta = math.pi * (1.0 + 5.0 ** 0.5) * k
        unit = np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)
        centers = st.p1_cart
        shell = own_radii + rp
        points = (centers[:, None, :] + shell[:, None, None] * unit[None, :, :]).reshape(-1, 3)
        inv = np.linalg.inv(cellm)
        points = ((points @ inv) % 1.0) @ cellm                    # repli dans la maille
        dist = self._surface_distance(points, atoms, radii).reshape(n_atoms, samples_per_atom)
        accessible = (dist >= rp - 1e-6).mean(axis=1)
        area = float(np.sum(4.0 * math.pi * shell ** 2 * accessible))   # Å2 par maille
        res.surface_m2_g = area * 1e-20 * N_AVOGADRO / st.p1_mass
        res.surface_m2_cm3 = area * 1e-20 / (st.volume * 1e-24)

        # --- grille de distance ---
        lengths = np.linalg.norm(cellm, axis=1)
        max_nodes = 3.0e6 if fast else 2.0e5
        spacing = max(grid_spacing, (st.volume / max_nodes) ** (1 / 3))
        shape = tuple(int(max(4, math.ceil(l / spacing))) for l in lengths)
        self.log("Grille de distance : %d x %d x %d noeuds (pas ~%.2f Å)..."
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
            res.notes.append("LCD supérieur à la limite mesurable (%.0f Å)." % res.lcd)

        # --- volume occupable par la sonde : dilatation de {D >= r_sonde} ---
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
        res.pore_volume_label = "volume occupable par la sonde (d = %.2f Å)" % self.probe_d

        # --- PLD par dichotomie sur le rayon de percolation ---
        self.log("Diamètre de pore limitant (percolation)...")
        lo, hi = 0.0, max(dmax, 0.0)
        for _ in range(12 if hi > 0 else 0):                  # lo reste nul si rien ne percole
            mid = 0.5 * (lo + hi)
            if self._percolates(dgrid >= mid):
                lo = mid
            else:
                hi = mid
        res.pld = 2.0 * lo
        res.notes.append("Précision PLD/LCD ~ %.2f Å (pas de grille)." % spacing)
        if res.pld < self.probe_d:
            res.notes.append("PLD < diamètre de la sonde : les fenêtres sont géométriquement "
                             "fermées pour une charpente rigide ; la surface et le volume "
                             "indiqués incluent des cavités non connectées.")
        return res

    def run(self, mode: str = "auto") -> GeometricResult:
        """mode = 'auto' | 'csd' | 'interne'."""
        unchanged = not any(k.startswith("invité:") for k in self.st.p1_removed)
        if mode in ("auto", "csd") and self.st.csd_cif_path and (unchanged or mode == "csd"):
            try:
                self.log("Pore Analyser (CSD Python API)...")
                return self.run_csd()
            except Exception as exc:
                self.log("Pore Analyser indisponible (%s) -> calcul interne." % exc)
        return self.run_internal()


# =============================================================================
# 6. DONNÉES EXPÉRIMENTALES
# =============================================================================
@dataclass
class Isotherm:
    """Isotherme : pression relative et quantité adsorbée en cm3(STP)/g."""
    p_rel: np.ndarray
    v_cm3: np.ndarray
    label: str = ""
    des_p_rel: np.ndarray = field(default_factory=lambda: np.zeros(0))
    des_v_cm3: np.ndarray = field(default_factory=lambda: np.zeros(0))
    meta: Dict[str, object] = field(default_factory=dict)

    def to_frame(self) -> pd.DataFrame:
        frame = pd.DataFrame({"Branche": "adsorption", "P/P0": self.p_rel,
                              "Adsorption_cm3STPg": self.v_cm3})
        if len(self.des_p_rel):
            frame = pd.concat([frame, pd.DataFrame({"Branche": "désorption", "P/P0": self.des_p_rel,
                                                    "Adsorption_cm3STPg": self.des_v_cm3})],
                              ignore_index=True)
        frame["Adsorption_mmol_g"] = frame["Adsorption_cm3STPg"] / (V_MOLAR_STP / 1000.0)
        return frame


class ExperimentalDataLoader:
    """Lecture tolérante des isothermes expérimentales.

    Formats acceptés : .xlsx, .xls, .csv/.txt et texte collé. Les en-têtes
    sont reconnus par mots-clés (« P/P0 », « Adsorption », « Va/cm3(STP)g-1 »,
    « Quantity Adsorbed »...). Les exports BELSORP/BELMaster (blocs ADS / DES,
    feuille « Isotherm ») et les tableaux simples à deux colonnes sont gérés.
    Si la quantité est en mmol/g ou mol/kg, elle est convertie en cm3(STP)/g.
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

    # ------------------------------------------------------------------ fichiers
    @classmethod
    def load_file(cls, path: str) -> Isotherm:
        ext = os.path.splitext(path)[1].lower()
        if ext in (".xlsx", ".xlsm", ".xls"):
            try:
                sheets = pd.read_excel(path, sheet_name=None, header=None)
            except ImportError as exc:
                need = "xlrd" if ext == ".xls" else "openpyxl"
                raise RuntimeError("Lecture %s impossible : installer le module '%s' "
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
                    iso.meta["feuille"] = str(name)
                    iso.label = os.path.basename(path)
                    return iso
                except ValueError as exc:
                    errors.append("%s : %s" % (name, exc))
            raise ValueError("Aucune isotherme trouvée dans %s.\n%s" % (path, "\n".join(errors)))
        with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
            iso = cls.parse_text(handle.read())
        iso.label = os.path.basename(path)
        return iso

    # --------------------------------------------------------------- texte collé
    @classmethod
    def parse_text(cls, text: str) -> Isotherm:
        """Convertit un tableau collé (tabulations, espaces, ';' ou ',')."""
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
                cells = line.split(",")                      # CSV anglo-saxon
            elif "," in line and "." in line:
                cells = line.split(",")
            else:
                cells = line.split()                         # espaces (virgule = décimale)
            rows.append([c.strip() for c in cells])
        if not rows:
            raise ValueError("Aucune donnée à lire.")
        width = max(len(r) for r in rows)
        frame = pd.DataFrame([r + [None] * (width - len(r)) for r in rows])
        iso = cls._from_table(frame)
        iso.label = "données collées"
        return iso

    # ------------------------------------------------------------- coeur du parseur
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
        # 1) recherche de la ligne d'en-tête
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
        # 2) sans en-tête : deux premières colonnes numériques, ou paires successives
        def headerless(rows):
            data_rows = [[x for x in row if np.isfinite(x)] for row in rows]
            data_rows = [row for row in data_rows if len(row) >= 2]
            if not data_rows:
                raise ValueError("colonnes P/P0 et adsorption introuvables")
            if len(data_rows) == 1 and len(data_rows[0]) >= 4 and len(data_rows[0]) % 2 == 0:
                flat = data_rows[0]                          # « 0.0001 15.3 0.001 31.3 ... »
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
            if len(ads) < 3:                                 # en-tête puis valeurs sur une seule ligne
                ads, des = headerless(numeric[header_row + 1:]), []
        ads = [(p, v) for p, v in ads if p > 0]
        if len(ads) < 3:
            raise ValueError("moins de 3 points d'adsorption exploitables")
        p_arr = np.array([p for p, _ in ads])
        if p_arr.max() > 1.5:
            raise ValueError("la colonne de pression n'est pas une pression relative P/P0 "
                             "(valeur maximale %.3g)" % p_arr.max())
        return Isotherm(p_arr, np.array([v for _, v in ads]),
                        des_p_rel=np.array([p for p, _ in des]), des_v_cm3=np.array([v for _, v in des]))


# =============================================================================
# 7. ANALYSE BET (critères de Rouquerol) ET VOLUME POREUX (Gurvich)
# =============================================================================
@dataclass
class BETResult:
    """Résultat d'un ajustement BET."""
    valid: bool = False
    message: str = ""
    p_min: float = float("nan")
    p_max: float = float("nan")
    n_points: int = 0
    slope: float = float("nan")
    intercept: float = float("nan")
    r2: float = float("nan")
    vm_cm3_g: float = float("nan")            # capacité de monocouche
    c_constant: float = float("nan")
    surface_m2_g: float = float("nan")
    p_at_vm: float = float("nan")             # P/P0 où n = n_m (critère 3)
    criterion4_error: float = float("nan")    # |1/(sqrt(C)+1) - p(n_m)| / p(n_m)
    rouquerol_ok: Tuple[bool, bool, bool, bool] = (False, False, False, False)
    plot_p: np.ndarray = field(default_factory=lambda: np.zeros(0))       # tous les points
    plot_y: np.ndarray = field(default_factory=lambda: np.zeros(0))       # p/[v(1-p)]
    in_range: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))

    def summary(self) -> Dict[str, object]:
        return {"Plage P/P0 min": self.p_min, "Plage P/P0 max": self.p_max,
                "Nombre de points": self.n_points, "Pente": self.slope,
                "Ordonnée à l'origine": self.intercept, "R2": self.r2,
                "Vm (cm3 STP/g)": self.vm_cm3_g, "Constante C": self.c_constant,
                "Surface BET (m2/g)": self.surface_m2_g, "P/P0 à Vm": self.p_at_vm,
                "Écart critère 4 (%)": 100.0 * self.criterion4_error,
                "Rouquerol 1 (n(1-p) croissant)": self.rouquerol_ok[0],
                "Rouquerol 2 (C > 0)": self.rouquerol_ok[1],
                "Rouquerol 3 (p(Vm) dans la plage)": self.rouquerol_ok[2],
                "Rouquerol 4 (1/(sqrt(C)+1) ~ p(Vm))": self.rouquerol_ok[3],
                "Remarque": self.message}


class BETAnalyzer:
    """Méthode BET classique avec sélection automatique de la plage.

    Équation linéarisée :  p / [v (1 - p)] = 1/(vm C) + (C - 1)/(vm C) * p ,
    avec p = P/P0. Toutes les fenêtres contiguës d'au moins `min_points`
    points sont ajustées, puis filtrées par les critères de Rouquerol
    (Rouquerol, Llewellyn, Rouquerol, Stud. Surf. Sci. Catal. 2007, 160, 49) :
      1. v(1 - p) croît de façon monotone avec p sur la plage ;
      2. C > 0 (ordonnée à l'origine positive) ;
      3. la pression p(vm) où la monocouche est atteinte est dans la plage ;
      4. 1/(sqrt(C) + 1) ne s'écarte pas de p(vm) de plus de 20 %.
    Parmi les fenêtres conformes (R2 >= 0,995), on retient celle qui contient
    le plus de points puis, à égalité, le meilleur R2. Si aucune fenêtre ne
    satisfait tout, les exigences sont relâchées et le résultat le signale.
    Pour une isotherme simulée, bruitée par nature, `noise_tolerance` autorise
    de petites baisses de v(1 - p) (de l'ordre de l'incertitude statistique).
    """

    def __init__(self, cross_section_nm2: float = 0.162, min_points: int = 3,
                 max_p_rel: float = 0.5, min_r2: float = 0.995, noise_tolerance: float = 0.0):
        self.cross_section = cross_section_nm2
        # Baisse relative de v(1-p) tolérée entre deux points pour le critère 1 :
        # 0 pour une mesure ; ~2 incertitudes statistiques pour une simulation.
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
            result.message = "Pas assez de points pour une analyse BET."
            return result
        y = p / (v * (1.0 - p))
        rouq = v * (1.0 - p)                                  # fonction de Rouquerol
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
                # p(vm) : interpolation (linéaire en ln p) de l'isotherme à v = vm
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
            result.message = "Aucune fenêtre BET exploitable (pente négative partout)."
            return result

        def pick(pool):
            return max(pool, key=lambda d: (d["j"] - d["i"] + 1, d["r2"]))

        tiers = [
            ("", [d for d in candidates if all(d["flags"]) and d["r2"] >= self.min_r2]),
            ("Critère 4 de Rouquerol non satisfait (écart > 20 %).",
             [d for d in candidates if all(d["flags"][:3]) and d["r2"] >= self.min_r2]),
            ("Critère 3 de Rouquerol non satisfait : p(Vm) hors de la plage ajustée.",
             [d for d in candidates if all(d["flags"][:2]) and d["r2"] >= self.min_r2]),
            ("Critères de Rouquerol non satisfaits : valeur indicative seulement.",
             [d for d in candidates if d["flags"][1]]),
            ("Constante C négative : la méthode BET ne s'applique pas à ces données.", candidates),
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
    """Volume poreux (cm3/g) par la règle de Gurvich à P/P0 = p_target.

    La quantité adsorbée est supposée à l'état liquide :
        Vp = v_ads / 22414 * Vm_liquide.
    Si p_target est hors de la plage mesurée, le point le plus proche est
    utilisé. Renvoie (volume, P/P0 effectivement utilisé).
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
# 8. PRÉPARATION ET EXÉCUTION DE RASPA3
# =============================================================================
@dataclass
class SimulationSettings:
    """Tous les paramètres réglables depuis l'interface."""
    gas: str = "N2"
    temperature: float = 77.0                 # K
    p0: float = 97190.0                       # Pa : pression de saturation à T
    n_init_cycles: int = 2000                 # « équilibration » (Initialization RASPA3)
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
    n_jobs: int = 1                           # simulations simultanées
    work_dir: str = ""

    def pressure_grid(self) -> np.ndarray:
        """Pressions relatives, espacées logarithmiquement."""
        return np.logspace(math.log10(self.p_rel_min), math.log10(self.p_rel_max), int(self.n_points))


def find_raspa_executable() -> str:
    """Cherche l'exécutable raspa3 aux emplacements usuels (macOS, Linux, conda)."""
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
    """Écrit les fichiers d'entrée RASPA3 (format JSON de la version 3)."""

    FRAMEWORK_NAME = "framework"

    def __init__(self, structure: StructureData, settings: SimulationSettings):
        self.st = structure
        self.s = settings
        self.gas = ADSORBATES[settings.gas]

    @property
    def adsorbate_charged(self) -> bool:
        return any(abs(pa["charge"]) > 0 for pa in self.gas["pseudo_atoms"])

    def validate(self) -> List[str]:
        """Renvoie la liste des problèmes bloquants (vide si tout est correct)."""
        problems = []
        if not self.st.p1_elements:
            problems.append("La structure P1 est vide.")
        for el in sorted(set(self.st.p1_elements)):
            try:
                framework_lj(el, self.s.forcefield)
            except ValueError as exc:
                problems.append(str(exc))
        if self.s.charge_mode == CHARGE_MODES[1] and self.st.p1_charges is None:
            problems.append("Le CIF ne contient pas de charges (_atom_site_charge) : "
                            "choisir un autre mode de charges.")
        if not (0 < self.s.p_rel_min < self.s.p_rel_max <= 1.5):
            problems.append("Plage de pressions relatives invalide.")
        return problems

    def force_field(self) -> dict:
        """force_field.json : pseudo-atomes de charpente (par élément) + adsorbat."""
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
        """<gaz>.json : constantes critiques (fugacité Peng-Robinson) + géométrie rigide."""
        return {"CriticalTemperature": self.gas["Tc"], "CriticalPressure": self.gas["Pc"],
                "AcentricFactor": self.gas["omega"], "pseudoAtoms": self.gas["sites"]}

    def simulation(self, pressure_pa: float) -> dict:
        """simulation.json pour un point de pression (GCMC, ensemble muVT)."""
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
        # « FugacityCoefficient » volontairement absent : RASPA3 le calcule
        # avec l'équation d'état de Peng-Robinson (constantes critiques ci-dessus).
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
    """Erreur d'exécution de RASPA3 (message destiné à l'utilisateur)."""


class RaspaRunner:
    """Lance RASPA3 pour chaque pression, suit la progression, lit les résultats.

    Un sous-dossier par point (`point_01_...`) : les simulations sont
    indépendantes et peuvent tourner en parallèle (`n_jobs`). Les sorties
    standard et d'erreur de RASPA3 sont conservées dans chaque dossier
    (raspa_stdout.txt / raspa_stderr.txt) pour le diagnostic.
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
        """Demande d'arrêt (thread-safe) : les processus RASPA3 sont terminés."""
        self._cancel.set()

    # ---------------------------------------------------------------- utilitaires
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
        """Avancement (0-1) d'une simulation en cours, lu dans son fichier de sortie."""
        matches = self.RE_CYCLE.findall(self._output_text(directory)[-20000:])
        if not matches:
            return 0.0
        init, cur, _total = matches[-1]
        done = int(cur) + (0 if init else self.s.n_init_cycles)
        return min(1.0, done / max(1, self.s.n_init_cycles + self.s.n_prod_cycles))

    def _parse_point(self, directory: str) -> Tuple[float, float, float]:
        """(charge mol/kg, incertitude mol/kg, molécules par maille) d'un point terminé."""
        text = self._output_text(directory)
        found = self.RE_LOADING.findall(text)
        if not found:
            raise RaspaError("résultat de charge introuvable dans le fichier de sortie")
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

    # ------------------------------------------------------------------ exécution
    def run(self) -> pd.DataFrame:
        s = self.s
        if not s.raspa_path or not os.path.isfile(s.raspa_path):
            raise RaspaError("Exécutable RASPA3 introuvable : '%s'.\nIndiquer son chemin dans "
                             "l'onglet « Paramètres RASPA »." % s.raspa_path)
        if not os.access(s.raspa_path, os.X_OK):
            raise RaspaError("Le fichier '%s' n'est pas exécutable (chmod +x)." % s.raspa_path)
        builder = RaspaInputBuilder(self.st, s)
        problems = builder.validate()
        if problems:
            raise RaspaError("Paramètres invalides :\n - " + "\n - ".join(problems))
        os.makedirs(s.work_dir, exist_ok=True)
        grid = s.pressure_grid()
        points = []
        for k, p_rel in enumerate(grid):
            directory = os.path.join(s.work_dir, "point_%02d_%.3e" % (k + 1, p_rel))
            builder.write_point(directory, p_rel * s.p0)
            points.append({"index": k, "p_rel": float(p_rel), "dir": directory,
                           "state": "attente", "loading": float("nan"),
                           "error": float("nan"), "per_uc": float("nan"), "message": ""})
        n_atoms = len(self.st.p1_elements) * int(np.prod(s.supercell))
        self.log("Fichiers d'entrée écrits dans : %s" % s.work_dir)
        self.log("%d points, %s à %.2f K, P0 = %.0f Pa, supercellule %s (%d atomes), "
                 "%d + %d cycles, %d simulation(s) simultanée(s)."
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
                    raise RaspaError("Calcul interrompu par l'utilisateur.")
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
                        raise RaspaError("Impossible de lancer RASPA3 : %s" % exc) from exc
                    point["state"] = "en cours"
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
                            raise RaspaError("code de retour %d" % code)
                        point["loading"], point["error"], point["per_uc"] = self._parse_point(point["dir"])
                        point["state"] = "terminé"
                        self.log("P/P0 = %.3e : %.4f ± %.4f mol/kg  (%.1f cm3 STP/g)"
                                 % (point["p_rel"], point["loading"], point["error"],
                                    point["loading"] * V_MOLAR_STP / 1000.0))
                    except RaspaError as exc:
                        point["state"] = "échec"
                        details = self._tail(os.path.join(point["dir"], "raspa_stderr.txt")) or \
                            self._tail(os.path.join(point["dir"], "raspa_stdout.txt")) or \
                            self._output_text(point["dir"])[-1500:]
                        point["message"] = "%s\n%s" % (exc, details)
                        self.log("ERREUR RASPA3 à P/P0 = %.3e : %s" % (point["p_rel"], point["message"]))
                running = still
                partial = sum(self._point_fraction(pt["dir"]) for pt, _, _, _ in running)
                fraction = (finished + partial) / len(points)
                elapsed = time.time() - t0
                eta = elapsed * (1 - fraction) / fraction if fraction > 0.02 else float("nan")
                self.progress(fraction, "%d/%d points terminés — écoulé %s%s"
                              % (finished, len(points), _dt.timedelta(seconds=int(elapsed)),
                                 "" if not np.isfinite(eta) else ", reste ~%s" % _dt.timedelta(seconds=int(eta))))
        finally:
            for point, proc, out, err in running:                # arrêt propre si interruption
                if proc.poll() is None:
                    proc.terminate()
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                out.close()
                err.close()
        failed = [pt for pt in points if pt["state"] != "terminé"]
        if len(failed) == len(points):
            raise RaspaError("Toutes les simulations RASPA3 ont échoué.\nPremier message :\n%s"
                             % failed[0]["message"])
        if failed:
            self.log("ATTENTION : %d point(s) en échec, ignorés dans l'analyse." % len(failed))
        frame = pd.DataFrame({
            "P/P0": [pt["p_rel"] for pt in points],
            "Pression_Pa": [pt["p_rel"] * s.p0 for pt in points],
            "Adsorption_mmol_g": [pt["loading"] for pt in points],            # mol/kg = mmol/g
            "Incertitude_mmol_g": [pt["error"] for pt in points],
            "Adsorption_cm3STPg": [pt["loading"] * V_MOLAR_STP / 1000.0 for pt in points],
            "Incertitude_cm3STPg": [pt["error"] * V_MOLAR_STP / 1000.0 for pt in points],
            "Molecules_par_maille": [pt["per_uc"] for pt in points],
            "Statut": [pt["state"] for pt in points],
        })
        return frame


# =============================================================================
# 9. ANALYSE COMPLÈTE ET COMPARAISON EXPÉRIENCE / SIMULATION
# =============================================================================
@dataclass
class AnalysisResults:
    """Conteneur de tous les résultats d'une session."""
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

    # ------------------------------------------------------------------ calculs
    def analyse(self) -> None:
        """(Re)calcule BET et volumes poreux à partir des isothermes disponibles."""
        gas = ADSORBATES[self.settings.gas]
        bet = BETAnalyzer(self.settings.cross_section)
        if self.experimental is not None:
            self.bet_exp = bet.analyse(self.experimental.p_rel, self.experimental.v_cm3)
            self.vpore_exp, self.vpore_exp_p = gurvich_pore_volume(
                self.experimental.p_rel, self.experimental.v_cm3,
                self.settings.gurvich_p_rel, gas["liquid_molar_volume"])
        if self.simulated is not None:
            ok = self.simulated["Statut"] == "terminé"
            p = self.simulated.loc[ok, "P/P0"].values
            v = self.simulated.loc[ok, "Adsorption_cm3STPg"].values
            err = self.simulated.loc[ok, "Incertitude_cm3STPg"].values
            relative = err[v > 0] / v[v > 0]
            noise = float(np.clip(2.0 * np.nanmedian(relative), 0.0, 0.05)) if len(relative) else 0.0
            self.bet_sim = BETAnalyzer(self.settings.cross_section, noise_tolerance=noise).analyse(p, v)
            self.vpore_sim, self.vpore_sim_p = gurvich_pore_volume(
                p, v, self.settings.gurvich_p_rel, gas["liquid_molar_volume"])

    # ------------------------------------------------------------------ tableaux
    def comparison_frame(self) -> pd.DataFrame:
        """Tableau comparatif expérience / simulation / géométrie."""
        be, bs, geo = self.bet_exp, self.bet_sim, self.geometric
        s_exp = be.surface_m2_g if be and be.valid else float("nan")
        s_sim = bs.surface_m2_g if bs and bs.valid else float("nan")
        s_geo = geo.surface_m2_g if geo else float("nan")
        rows = [
            ("Surface BET (m2/g)", s_exp, s_sim, float("nan"), pct_diff(s_sim, s_exp)),
            ("Vm monocouche (cm3 STP/g)", be.vm_cm3_g if be else float("nan"),
             bs.vm_cm3_g if bs else float("nan"), float("nan"),
             pct_diff(bs.vm_cm3_g if bs else float("nan"), be.vm_cm3_g if be else float("nan"))),
            ("Constante C", be.c_constant if be else float("nan"),
             bs.c_constant if bs else float("nan"), float("nan"), float("nan")),
            ("Volume poreux Gurvich (cm3/g)", self.vpore_exp, self.vpore_sim, float("nan"),
             pct_diff(self.vpore_sim, self.vpore_exp)),
            ("P/P0 utilisé pour Gurvich", self.vpore_exp_p, self.vpore_sim_p, float("nan"), float("nan")),
            ("Volume poreux géométrique (cm3/g)", float("nan"), float("nan"),
             geo.pore_volume_cm3_g if geo else float("nan"), float("nan")),
            ("Surface géométrique accessible (m2/g)", float("nan"), float("nan"), s_geo, float("nan")),
            ("BET expérimentale / surface géométrique (%)", float("nan"), float("nan"),
             100.0 * s_exp / s_geo if s_geo and np.isfinite(s_geo) and s_geo > 0 else float("nan"), float("nan")),
            ("BET simulée / surface géométrique (%)", float("nan"), float("nan"),
             100.0 * s_sim / s_geo if s_geo and np.isfinite(s_geo) and s_geo > 0 else float("nan"), float("nan")),
            ("Diamètre de pore limitant PLD (Å)", float("nan"), float("nan"),
             geo.pld if geo else float("nan"), float("nan")),
            ("Diamètre maximal de cavité LCD (Å)", float("nan"), float("nan"),
             geo.lcd if geo else float("nan"), float("nan")),
        ]
        frame = pd.DataFrame(rows, columns=["Grandeur", "Expérimental", "Simulation RASPA3",
                                            "Géométrique", "Écart / expérience (%)"])
        # Écart volume géométrique / volume expérimental, placé sur la ligne géométrique
        if geo and np.isfinite(self.vpore_exp):
            frame.loc[5, "Écart / expérience (%)"] = pct_diff(geo.pore_volume_cm3_g, self.vpore_exp)
        return frame

    def bet_frame(self) -> pd.DataFrame:
        empty = BETResult().summary()
        exp = self.bet_exp.summary() if self.bet_exp else empty
        sim = self.bet_sim.summary() if self.bet_sim else empty
        return pd.DataFrame({"Paramètre": list(exp.keys()), "Expérimental": list(exp.values()),
                             "Simulation RASPA3": [sim[k] for k in exp]})

    def bet_points_frame(self) -> pd.DataFrame:
        frames = []
        for name, bet in (("expérimental", self.bet_exp), ("simulation", self.bet_sim)):
            if bet is not None and len(bet.plot_p):
                frames.append(pd.DataFrame({"Série": name, "P/P0": bet.plot_p,
                                            "p/[v(1-p)] (g/cm3 STP)": bet.plot_y,
                                            "Dans la plage BET": bet.in_range}))
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def structure_frame(self) -> pd.DataFrame:
        st = self.structure
        raw = st.asymmetric_summary()
        raw_mass = st.mass_of(raw)
        a, b, c, al, be, ga = st.cell
        rows = [
            ("Nom", st.name), ("Fichier", st.source_file or "(structure Mercury)"),
            ("Formule brute déclarée", st.formula_reported or "—"),
            ("Masse molaire de la formule (g/mol)", fmt(st.formula_weight_reported, 2)),
            ("Z", fmt(st.z_value, 0)),
            ("Groupe d'espace", st.spacegroup),
            ("a, b, c (Å)", "%.4f, %.4f, %.4f" % (a, b, c)),
            ("alpha, beta, gamma (°)", "%.3f, %.3f, %.3f" % (al, be, ga)),
            ("Volume de maille (Å3)", fmt(st.volume, 1)),
            ("Contenu de la maille (occupations)", hill_formula(raw)),
            ("Masse de la maille (g/mol)", fmt(raw_mass, 2)),
            ("Densité calculée, structure lue (g/cm3)", fmt(st.density_from_mass(raw_mass), 4)),
            ("Densité déclarée dans le fichier (g/cm3)", fmt(st.density_reported, 4)),
            ("Modèle de simulation P1", "%d atomes : %s" % (len(st.p1_elements), hill_formula(st.p1_counts()))),
            ("Masse du modèle P1 (g/mol)", fmt(st.p1_mass, 2)),
            ("Densité du modèle P1 (g/cm3)", fmt(st.density_from_mass(st.p1_mass), 4)),
            ("Atomes écartés du modèle P1", st.removed_summary()),
        ]
        return pd.DataFrame(rows, columns=["Propriété", "Valeur"])

    def settings_frame(self) -> pd.DataFrame:
        s = self.settings
        gas = ADSORBATES[s.gas]
        rows = [("Adsorbat", gas["label"]), ("Modèle d'adsorbat", gas["model"]),
                ("Température (K)", s.temperature), ("P0 (Pa)", s.p0),
                ("Champ de force", s.forcefield), ("Charges de la charpente", s.charge_mode),
                ("Rayon de coupure (Å)", s.cutoff), ("Supercellule", "x".join(map(str, s.supercell))),
                ("Cycles d'équilibration (Initialization)", s.n_init_cycles),
                ("Cycles de production", s.n_prod_cycles),
                ("P/P0 min", s.p_rel_min), ("P/P0 max", s.p_rel_max), ("Nombre de points", s.n_points),
                ("Section d'encombrement BET (nm2)", s.cross_section),
                ("P/P0 Gurvich", s.gurvich_p_rel),
                ("Volume molaire du liquide (cm3/mol)", gas["liquid_molar_volume"]),
                ("Exécutable RASPA3", s.raspa_path), ("Dossier de travail", s.work_dir),
                ("Date", self.created), ("Logiciel", "%s %s" % (APP_NAME, APP_VERSION))]
        return pd.DataFrame(rows, columns=["Paramètre", "Valeur"])


# =============================================================================
# 10. FIGURES (matplotlib, sans pyplot : utilisables dans Qt et hors écran)
# =============================================================================
class FigureFactory:
    """Construit les figures à partir d'un `AnalysisResults`."""

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
        """Isothermes expérimentale et simulée (échelles linéaire et logarithmique)."""
        r = self.r
        fig = self._new(10.0, 4.4)
        axes = fig.subplots(1, 2)
        for ax, log in zip(axes, (False, True)):
            self._style(ax)
            if r.experimental is not None:
                e = r.experimental
                ax.plot(e.p_rel, e.v_cm3, "o-", color=COLOR_EXP, ms=4.5, lw=1.0,
                        label="Expérimental (adsorption)")
                if len(e.des_p_rel):
                    ax.plot(e.des_p_rel, e.des_v_cm3, "o--", color=COLOR_EXP, ms=4.5, lw=0.8,
                            mfc="white", label="Expérimental (désorption)")
            if r.simulated is not None:
                sim = r.simulated[r.simulated["Statut"] == "terminé"]
                ax.errorbar(sim["P/P0"], sim["Adsorption_cm3STPg"], yerr=sim["Incertitude_cm3STPg"],
                            fmt="s-", color=COLOR_SIM, ms=4.5, lw=1.0, capsize=2,
                            label="Simulation RASPA3 (GCMC)")
            ax.set_xlabel("P/P0" + (" (échelle log)" if log else ""))
            ax.set_ylabel("Adsorption (cm3 STP/g)")
            if log:
                ax.set_xscale("log")
            else:
                ax.set_xlim(0, 1.02)
            ax.set_ylim(bottom=0)
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            axes[0].legend(handles, labels, frameon=False, fontsize=9, loc="best")
        fig.suptitle("Isotherme d'adsorption %s à %.4g K — %s"
                     % (r.settings.gas, r.settings.temperature, r.structure.name), fontsize=11)
        return fig

    # Figures 2 et 3 -------------------------------------------------------------
    def bet(self, which: str):
        """Tracé BET (« exp » ou « sim ») avec la plage retenue et la droite ajustée."""
        bet = self.r.bet_exp if which == "exp" else self.r.bet_sim
        color = COLOR_EXP if which == "exp" else COLOR_SIM
        title = "BET expérimental" if which == "exp" else "BET simulé (RASPA3)"
        fig = self._new(5.6, 4.4)
        ax = fig.subplots()
        self._style(ax)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("P/P0")
        ax.set_ylabel("p / [v (1 - p)]  (g / cm3 STP)")
        if bet is None or not len(bet.plot_p):
            ax.text(0.5, 0.5, "Données non disponibles", ha="center", va="center",
                    transform=ax.transAxes, color="#666666")
            return fig
        if not bet.valid:
            ax.plot(bet.plot_p, bet.plot_y, "o", color="#999999", ms=4)
            ax.text(0.5, 0.9, bet.message, ha="center", transform=ax.transAxes, fontsize=9)
            return fig
        x_max = min(max(bet.p_max * 1.8, bet.p_max + 0.02), 1.0)
        shown = bet.plot_p <= x_max
        ax.plot(bet.plot_p[shown & ~bet.in_range], bet.plot_y[shown & ~bet.in_range], "o",
                mfc="white", color="#888888", ms=5, label="Hors plage")
        ax.plot(bet.plot_p[bet.in_range], bet.plot_y[bet.in_range], "o", color=color, ms=6,
                label="Plage BET (%d points)" % bet.n_points)
        xs = np.array([0.0, x_max])
        ax.plot(xs, bet.slope * xs + bet.intercept, "-", color="#222222", lw=1.0, label="Ajustement linéaire")
        ax.set_xlim(0, x_max)
        ax.set_ylim(bottom=0)
        flags = "".join("✓" if f else "✗" for f in bet.rouquerol_ok)
        text = ("S(BET) = %.1f m2/g\nVm = %.2f cm3 STP/g\nC = %.4g\nR2 = %.5f\n"
                "Plage : %.2e – %.2e\nRouquerol 1-4 : %s"
                % (bet.surface_m2_g, bet.vm_cm3_g, bet.c_constant, bet.r2, bet.p_min, bet.p_max, flags))
        ax.text(0.97, 0.05, text, transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", edgecolor="#cccccc"))
        ax.legend(frameon=False, fontsize=9, loc="upper left")
        return fig

    # Figure 4 ------------------------------------------------------------------
    def bars(self):
        """Histogramme : surface géométrique, BET expérimentale, BET simulée."""
        r = self.r
        values = [r.geometric.surface_m2_g if r.geometric else float("nan"),
                  r.bet_exp.surface_m2_g if r.bet_exp and r.bet_exp.valid else float("nan"),
                  r.bet_sim.surface_m2_g if r.bet_sim and r.bet_sim.valid else float("nan")]
        labels = ["Surface\ngéométrique", "BET\nexpérimentale", "BET\nsimulée"]
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
        ax.set_ylabel("Surface spécifique (m2/g)")
        ax.set_title("Surfaces spécifiques comparées", fontsize=11)
        return fig

    # Image de la structure -------------------------------------------------------
    def structure(self):
        """Vue de la maille P1 (projection orthographique le long de c).

        Alternative documentée à une capture de la fenêtre Mercury, que la CSD
        Python API ne permet pas d'obtenir depuis un script.
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
        if len(elements) <= 3000:                               # liaisons (image minimale)
            radii = np.array([COVALENT_RADIUS.get(e, 1.5) for e in elements])
            pi, pj, pd_ = st.neighbor_pairs(st.p1_frac, 2.9)
            bonded = pd_ < 1.2 * (radii[pi] + radii[pj])
            delta = st.p1_frac[pj[bonded]] - st.p1_frac[pi[bonded]]
            inside = np.all(np.abs(delta) < 0.5, axis=1)        # liaisons ne traversant pas la maille
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
        ax.set_title("%s — maille P1 vue selon c" % st.name, fontsize=10)
        ax.autoscale_view()
        return fig

    def all_figures(self) -> Dict[str, object]:
        return {"isotherme": self.isotherm(), "bet_exp": self.bet("exp"),
                "bet_sim": self.bet("sim"), "surfaces": self.bars()}


# =============================================================================
# 11. EXPORTS : EXCEL, CSV, PDF, HTML
# =============================================================================
class ReportBuilder:
    """Écrit tous les fichiers de résultats dans un dossier."""

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
            return pd.DataFrame(columns=["Branche", "P/P0", "Adsorption_cm3STPg", "Adsorption_mmol_g"])
        return self.r.experimental.to_frame()

    def _sim_frame(self) -> pd.DataFrame:
        if self.r.simulated is None:
            return pd.DataFrame(columns=["P/P0", "Adsorption_mmol_g", "Adsorption_cm3STPg"])
        return self.r.simulated

    # ---------------------------------------------------------------------- Excel
    def write_excel(self) -> Optional[str]:
        """Classeur : 1 isotherme exp., 2 isotherme simulée, 3 BET, 4 comparaison
        (+ feuilles Structure et Paramètres pour la traçabilité)."""
        path = self._path("_resultats.xlsx")
        engine = None
        for cand in ("openpyxl", "xlsxwriter"):
            try:
                importlib.import_module(cand)
                engine = cand
                break
            except ImportError:
                continue
        if engine is None:
            self.log("Export Excel impossible : installer 'openpyxl' (les CSV restent écrits).")
            return None
        with pd.ExcelWriter(path, engine=engine) as writer:
            self._exp_frame().to_excel(writer, sheet_name="Isotherme_exp", index=False)
            self._sim_frame().to_excel(writer, sheet_name="Isotherme_sim", index=False)
            self.r.bet_frame().to_excel(writer, sheet_name="BET", index=False)
            points = self.r.bet_points_frame()
            if len(points):
                points.to_excel(writer, sheet_name="BET", index=False, startrow=len(self.r.bet_frame()) + 3)
            self.r.comparison_frame().to_excel(writer, sheet_name="Comparaison", index=False)
            self.r.structure_frame().to_excel(writer, sheet_name="Structure", index=False)
            self.r.settings_frame().to_excel(writer, sheet_name="Parametres", index=False)
            if engine == "openpyxl":                            # largeur des colonnes
                for sheet in writer.sheets.values():
                    for column in sheet.columns:
                        width = max(len(str(c.value)) if c.value is not None else 0 for c in column)
                        sheet.column_dimensions[column[0].column_letter].width = min(max(12, width + 2), 60)
        return path

    # ------------------------------------------------------------------------ CSV
    def write_csv(self) -> List[str]:
        paths = []
        for suffix, frame in (("_isotherme_exp.csv", self._exp_frame()),
                              ("_isotherme_sim.csv", self._sim_frame()),
                              ("_bet.csv", self.r.bet_frame()),
                              ("_comparaison.csv", self.r.comparison_frame())):
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
        """Rapport PDF : structure, paramètres cristallographiques, isothermes, BET, comparaison."""
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        r = self.r
        path = self._path("_rapport.pdf")
        a4 = (8.27, 11.69)
        with PdfPages(path) as pdf:
            # --- page 1 : structure ---
            fig = Figure(figsize=a4)
            FigureCanvasAgg(fig)
            fig.text(0.5, 0.955, "Rapport de porosité : expérience et simulation", ha="center",
                     fontsize=15, weight="bold")
            fig.text(0.5, 0.93, "%s — %s — %s" % (r.structure.name, ADSORBATES[r.settings.gas]["label"],
                                                   r.created), ha="center", fontsize=10, color="#444444")
            ax_img = fig.add_axes([0.2, 0.52, 0.6, 0.39])
            self._draw_structure(ax_img)
            ax_tab = fig.add_axes([0.06, 0.05, 0.88, 0.45])
            self._table(ax_tab, r.structure_frame(), (0.42, 0.58))
            pdf.savefig(fig)
            # --- page 2 : isothermes + paramètres ---
            fig = Figure(figsize=a4)
            FigureCanvasAgg(fig)
            fig.text(0.06, 0.955, "Isothermes d'adsorption", fontsize=13, weight="bold")
            self._embed(fig, self.figures.isotherm(), [0.04, 0.58, 0.92, 0.35])
            ax_tab = fig.add_axes([0.06, 0.04, 0.88, 0.5])
            self._table(ax_tab, r.settings_frame().astype(str), (0.42, 0.58))
            pdf.savefig(fig)
            # --- page 3 : BET ---
            fig = Figure(figsize=a4)
            FigureCanvasAgg(fig)
            fig.text(0.06, 0.955, "Analyse BET (critères de Rouquerol)", fontsize=13, weight="bold")
            self._embed(fig, self.figures.bet("exp"), [0.03, 0.6, 0.47, 0.33])
            self._embed(fig, self.figures.bet("sim"), [0.5, 0.6, 0.47, 0.33])
            ax_tab = fig.add_axes([0.06, 0.05, 0.88, 0.52])
            self._table(ax_tab, r.bet_frame(), (0.46, 0.27, 0.27))
            pdf.savefig(fig)
            # --- page 4 : comparaison ---
            fig = Figure(figsize=a4)
            FigureCanvasAgg(fig)
            fig.text(0.06, 0.955, "Comparaison expérience / simulation / géométrie", fontsize=13, weight="bold")
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
            notes.append("Géométrie : %s ; sonde %.2f Å ; %s." % (r.geometric.method, r.geometric.probe_diameter,
                                                                 r.geometric.pore_volume_label))
            notes.extend(self._wrap(n) for n in r.geometric.notes)
        for name, bet in (("BET expérimental", r.bet_exp), ("BET simulé", r.bet_sim)):
            if bet is not None and bet.message:
                notes.append("%s : %s" % (name, bet.message))
        notes.append("Simulation : adsorption absolue, charpente rigide ; expérience : adsorption d'excès.")
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
        """Insère une figure (rendue en PNG) dans une page du rapport PDF."""
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
        """Rapport HTML autonome (images intégrées) affiché dans Mercury."""
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
            "<h2>Comparaison</h2>", table(r.comparison_frame()), image(self.figures.bars()),
            "<h2>Isothermes</h2>", image(self.figures.isotherm()),
            "<h2>Analyse BET</h2>", image(self.figures.bet("exp")), image(self.figures.bet("sim")),
            table(r.bet_frame().astype(str)),
            "<h2>Structure</h2>", image(self.figures.structure()), table(r.structure_frame()),
            "<h2>Paramètres</h2>", table(r.settings_frame().astype(str)),
            "<p>" + "<br/>".join(n.replace("\n", " ") for n in self._notes()) + "</p>",
            "<p>Fichiers : %s</p>" % self.out_dir, "</body></html>"]
        return "\n".join(parts)

    def write_all(self) -> Dict[str, object]:
        """Écrit Excel, CSV, PDF, PNG et HTML ; renvoie les chemins créés."""
        written: Dict[str, object] = {"csv": self.write_csv()}
        written["excel"] = self.write_excel()
        written["pdf"] = self.write_pdf()
        written["png"] = self.write_png()
        html_path = self._path("_rapport.html")
        with open(html_path, "w", encoding="utf-8") as handle:
            handle.write(self.html())
        written["html"] = html_path
        self.log("Résultats exportés dans : %s" % self.out_dir)
        return written


# =============================================================================
# 12. LIAISON AVEC MERCURY (CSD Python API)
# =============================================================================
class MercuryBridge:
    """Accès à la structure ouverte dans Mercury.

    Mercury lance le script dans un processus Python séparé et lui transmet,
    en arguments, l'identifiant de la structure courante et les fichiers de
    sortie. Deux générations d'interface coexistent :
      * `ccdc.utilities.ApplicationInterface` (versions récentes) ;
      * `mercury_interface.MercuryInterface` (anciennes versions).
    Si aucune n'est utilisable (script lancé hors Mercury), `available` reste
    faux et le programme bascule sur la lecture directe d'un fichier CIF.
    """

    def __init__(self):
        self.interface = None
        self.available = False
        self.error = ""
        saved_stderr = sys.stderr
        try:
            sys.stderr = _io.StringIO()                         # argparse de l'API reste silencieux
            try:
                from ccdc.utilities import ApplicationInterface
                interface = ApplicationInterface(parse_commandline=False)
                interface.parse_commandline()
            except ImportError as first:                        # ancienne interface, ou ccdc cassé
                self.error = "%s: %s" % (type(first).__name__, first)
                from mercury_interface import MercuryInterface
                interface = MercuryInterface()
                self.error = ""
            self.interface = interface
            self.available = True
        except BaseException as exc:                            # y compris SystemExit d'argparse
            second = "%s: %s" % (type(exc).__name__, exc)
            self.error = second if not self.error or self.error == second else self.error + " | " + second
        finally:
            sys.stderr = saved_stderr

    def current_entry(self):
        entry = getattr(self.interface, "current_entry")
        return entry() if callable(entry) else entry

    def original_cif(self) -> str:
        """Chemin du fichier CIF chargé dans Mercury, s'il figure parmi les options."""
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
        """Exporte la structure courante vers un CIF temporaire (CSD Python API)."""
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
    """Vrai si la CSD Python API est installée (sans la charger)."""
    import importlib.util
    try:
        return importlib.util.find_spec("ccdc") is not None
    except (ImportError, ValueError):
        return False


class ChildBridge:
    """Remplace `MercuryBridge` dans le processus d'interface (qui n'a pas la CSD API)."""

    def __init__(self, html_out: str = "", from_mercury: bool = False):
        self._html_out = html_out
        self.available = from_mercury

    def output_html_file(self) -> str:
        return self._html_out


def load_structure(cif_path: str = "") -> Tuple[StructureData, MercuryBridge, str]:
    """Charge la structure : Mercury en priorité, sinon le fichier CIF indiqué.

    Renvoie (structure, pont Mercury, chemin du CIF temporaire exporté).
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
                bridge.error = "Export CIF via l'API impossible : %s" % exc
            original = bridge.original_cif()
            structure = StructureData.from_csd_entry(entry, original)
            structure.csd_cif_path = temp_cif or original
            if temp_cif:
                structure.notes.append("CIF temporaire exporté par la CSD Python API : %s" % temp_cif)
            if bridge.error:
                structure.notes.append(bridge.error)
            return structure, bridge, temp_cif
        except Exception as exc:                                # pas de structure courante exploitable
            bridge.error = "Structure Mercury illisible : %s" % exc
            if not cif_path:
                raise RuntimeError(bridge.error) from exc
    if cif_path:
        structure = StructureData.from_cif_file(cif_path)
        structure.csd_cif_path = cif_path if ccdc_installed() else ""
        return structure, bridge, temp_cif
    raise RuntimeError("Aucune structure : lancer le script depuis Mercury ou indiquer un fichier CIF.\n"
                       "(interface Mercury : %s)" % (bridge.error or "indisponible"))


def default_work_dir(structure: StructureData, gas: str, base: str = "") -> str:
    base = base or os.path.join(os.path.expanduser("~"), "RASPA_Mercury")
    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", structure.name)[:40] or "structure"
    return os.path.join(base, "%s_%s_%s" % (safe, gas, stamp))


# =============================================================================
# 13. INTERFACE GRAPHIQUE Qt (PyQt6 / PyQt5 / PySide6 / PySide2)
# =============================================================================
# IMPORTANT (macOS surtout) : la CSD Python API embarque ses propres bibliothèques
# Qt (ccdc/_lib/Qt*.framework). Les charger dans le même processus que PyQt/PySide
# provoque un conflit de symboles (« Symbol not found: QWidgetPrivate... »).
# Le script fonctionne donc en DEUX processus :
#   * processus parent (lancé par Mercury) : CSD Python API, jamais Qt ;
#   * processus enfant (option --gui-child) : interface Qt, jamais la CSD Python API.
# La structure passe de l'un à l'autre par un fichier CIF de transfert.
# Sur macOS, la seule PRÉSENCE de PyQt dans le site-packages du Python de Mercury
# peut suffire à casser la CSD Python API. PyQt est donc installé à part, dans le
# dossier privé QT_PRIVATE_DIR (pip install --target), que seul le processus
# d'interface ajoute à son chemin de recherche.
QT_CANDIDATES = ("PyQt6", "PyQt5", "PySide6", "PySide2")
QT_PRIVATE_DIR = os.environ.get("RASPA_POROSITY_QT_DIR") or \
    os.path.join(os.path.expanduser("~"), "RASPA_Mercury", "qt_libs")


def qt_in_private_dir() -> str:
    """Liaison Qt présente dans le dossier privé ('' si aucune)."""
    for name in QT_CANDIDATES:
        if os.path.isdir(os.path.join(QT_PRIVATE_DIR, name)):
            return name
    return ""


def qt_in_site_packages() -> str:
    """Liaison Qt installée dans le Python courant, SANS la charger ('' si aucune)."""
    import importlib.util
    for name in QT_CANDIDATES:
        try:
            if importlib.util.find_spec(name) is not None:
                return name
        except (ImportError, ValueError):
            continue
    return ""


def detect_qt_binding() -> str:
    """Liaison Qt utilisable par le processus d'interface ('' si aucune)."""
    return qt_in_private_dir() or qt_in_site_packages()


def qt_fix_instructions() -> str:
    """Commandes Terminal exactes pour installer Qt sans perturber la CSD Python API."""
    python = sys.executable
    lines = ["Copier ces commandes dans le Terminal (une par une), puis relancer le script :", ""]
    if qt_in_site_packages():
        lines += ["1) retirer PyQt du Python de Mercury (il perturbe la CSD Python API) :",
                  '"%s" -m pip uninstall -y PyQt6 PyQt6-Qt6 PyQt6-sip PyQt5 PyQt5-Qt5 PyQt5-sip PySide6 PySide2'
                  % python, ""]
    lines += ["%s installer PyQt6 dans un dossier séparé :" % ("2)" if qt_in_site_packages() else "1)"),
              '"%s" -m pip install --target "%s" PyQt6' % (python, QT_PRIVATE_DIR), "",
              "%s modules complémentaires (sans risque pour Mercury) :" % ("3)" if qt_in_site_packages() else "2)"),
              '"%s" -m pip install matplotlib openpyxl xlrd scipy' % python, "",
              "(%s %s, Python %s)" % (APP_NAME, APP_VERSION, sys.version.split()[0])]
    return "\n".join(lines)


_LOAD_QT = "--gui-child" in sys.argv[1:] or os.environ.get("RASPA_POROSITY_QT") == "1"
if _LOAD_QT and os.path.isdir(QT_PRIVATE_DIR):
    sys.path.insert(0, QT_PRIVATE_DIR)                          # PyQt du dossier privé en priorité
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
    os.environ.setdefault("QT_API", QT_BINDING.lower())         # même liaison pour matplotlib
    import matplotlib
    try:
        matplotlib.use("QtAgg")
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
    except Exception:                                           # matplotlib < 3.5
        matplotlib.use("Qt5Agg")
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
    Signal = getattr(QtCore, "pyqtSignal", None) or getattr(QtCore, "Signal")

    class Worker(QtCore.QThread):
        """Exécute une fonction dans un fil séparé pour ne pas figer l'interface.

        La fonction reçoit deux rappels : log(texte) et progress(fraction, texte).
        """
        log = Signal(str)
        progress = Signal(float, str)
        done = Signal(object)
        failed = Signal(str)

        def __init__(self, function: Callable, parent=None):
            super().__init__(parent)
            self._function = function

        def run(self) -> None:                                  # appelé par QThread.start()
            try:
                self.done.emit(self._function(self.log.emit, self.progress.emit))
            except RaspaError as exc:
                self.failed.emit(str(exc))
            except Exception:
                self.failed.emit(traceback.format_exc())

    class FigureTab(QtWidgets.QWidget):
        """Onglet contenant une figure matplotlib et sa barre d'outils."""

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
        """Fenêtre principale : Structure, Données expérimentales, Paramètres RASPA,
        Calcul, Résultats."""

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
            self.statusBar().showMessage("Qt : %s — %s" % (
                QT_BINDING, "structure fournie par Mercury" if bridge and bridge.available
                else "mode autonome (lecture CIF interne)"))
            self._rebuild_model()
            self._on_gas_changed()
            if cli is not None:
                self._apply_cli(cli)

        # ---------------------------------------------------------------- utilitaires
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

        # ============================ Onglet 1 : Structure ============================
        def _build_structure_tab(self) -> None:
            page = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout(page)
            self.structure_view = QtWidgets.QTextBrowser()
            layout.addWidget(self.structure_view, 3)
            box = QtWidgets.QGroupBox("Préparation du modèle de simulation (maille P1)")
            form = QtWidgets.QVBoxLayout(box)
            self.disorder_check = QtWidgets.QCheckBox(
                "Résoudre le désordre (une seule composante pour les sites d'occupation partielle)")
            self.disorder_check.setChecked(True)
            self.guest_check = QtWidgets.QCheckBox(
                "Supprimer les espèces invitées non liées à la charpente (solvant, eau pseudo-localisée)")
            self.guest_check.setChecked(True)
            for check in (self.disorder_check, self.guest_check):
                check.stateChanged.connect(lambda _state: self._rebuild_model())
                form.addWidget(check)
            layout.addWidget(box)
            row = QtWidgets.QHBoxLayout()
            open_button = QtWidgets.QPushButton("Ouvrir un autre fichier CIF…")
            open_button.clicked.connect(self._open_cif)
            self.geometry_button = QtWidgets.QPushButton("Calculer les propriétés géométriques")
            self.geometry_button.clicked.connect(self._run_geometry)
            self.geometry_mode = QtWidgets.QComboBox()
            self.geometry_mode.addItems(["Auto (Pore Analyser CSD si possible)", "Pore Analyser (CSD Python API)",
                                         "Calcul interne"])
            row.addWidget(open_button)
            row.addStretch(1)
            row.addWidget(QtWidgets.QLabel("Méthode :"))
            row.addWidget(self.geometry_mode)
            row.addWidget(self.geometry_button)
            layout.addLayout(row)
            self.geometry_view = QtWidgets.QTextBrowser()
            layout.addWidget(self.geometry_view, 2)
            self.tabs.addTab(page, "1. Structure")

        def _open_cif(self) -> None:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Fichier CIF", "", "CIF (*.cif)")
            if not path:
                return
            try:
                self.structure = StructureData.from_cif_file(path)
                self.structure.csd_cif_path = path if ccdc_installed() else ""
            except Exception as exc:
                self._error("Lecture du CIF impossible :\n%s" % exc)
                return
            self.setWindowTitle("%s %s — %s" % (APP_NAME, APP_VERSION, self.structure.name))
            self._rebuild_model()

        def _rebuild_model(self) -> None:
            """Reconstruit la maille P1 selon les options, met à jour l'affichage."""
            try:
                self.structure.build_p1(self.disorder_check.isChecked(), self.guest_check.isChecked())
            except Exception as exc:
                self._error("Construction du modèle P1 impossible :\n%s" % exc)
                return
            self.geometric = None
            self.geometry_view.setHtml("<i>Propriétés géométriques non calculées "
                                       "(elles le seront aussi lors de l'analyse).</i>")
            dummy = AnalysisResults(self.structure, SimulationSettings())
            notes = self.structure.notes + self.structure.p1_log
            self.structure_view.setHtml(self._html_table(dummy.structure_frame()) +
                                        "<p>" + "<br/>".join(notes) + "</p>")
            if hasattr(self, "cell_spins"):
                self._auto_supercell()

        def _probe_diameter(self) -> float:
            return ADSORBATES[self.gas_combo.currentText()]["probe_diameter"]

        def _geometry_mode(self) -> str:
            return ("auto", "csd", "interne")[self.geometry_mode.currentIndex()]

        def _show_geometry(self, geo: GeometricResult) -> None:
            self.geometric = geo
            frame = pd.DataFrame([
                ("Méthode", geo.method), ("Diamètre de sonde (Å)", fmt(geo.probe_diameter, 3)),
                ("Surface accessible (m2/g)", fmt(geo.surface_m2_g, 1)),
                ("Surface accessible (m2/cm3)", fmt(geo.surface_m2_cm3, 1)),
                ("Volume poreux (cm3/g) — " + geo.pore_volume_label, fmt(geo.pore_volume_cm3_g, 4)),
                ("Volume géométrique, sonde ponctuelle (cm3/g)", fmt(geo.geometric_volume_cm3_g, 4)),
                ("Fraction de vide", fmt(geo.void_fraction, 4)),
                ("Diamètre de pore limitant PLD (Å)", fmt(geo.pld, 2)),
                ("Diamètre maximal de cavité LCD (Å)", fmt(geo.lcd, 2)),
            ], columns=["Propriété géométrique", "Valeur"])
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
            self.geometry_view.setHtml("<i>Calcul en cours…</i>")
            self._start_worker(job, self._geometry_done)

        def _geometry_done(self, geo: GeometricResult) -> None:
            self.geometry_button.setEnabled(True)
            self._show_geometry(geo)

        # ====================== Onglet 2 : Données expérimentales ======================
        def _build_experimental_tab(self) -> None:
            page = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout(page)
            box_a = QtWidgets.QGroupBox("Option A — Importer un fichier (xlsx, xls, csv)")
            row = QtWidgets.QHBoxLayout(box_a)
            self.exp_path = QtWidgets.QLineEdit()
            self.exp_path.setPlaceholderText("Colonnes attendues : P/P0 et Adsorption_cm3STPg "
                                             "(les exports BELSORP sont reconnus)")
            browse = QtWidgets.QPushButton("Parcourir…")
            browse.clicked.connect(self._browse_experimental)
            row.addWidget(self.exp_path)
            row.addWidget(browse)
            layout.addWidget(box_a)
            box_b = QtWidgets.QGroupBox("Option B — Coller un tableau (P/P0 puis adsorption en cm3 STP/g)")
            col = QtWidgets.QVBoxLayout(box_b)
            self.paste_edit = QtWidgets.QPlainTextEdit()
            self.paste_edit.setPlaceholderText("P/P0 Adsorption\n0.0001 15.3\n0.0010 31.3\n0.0100 36.6\n"
                                               "0.0500 37.7\n0.1000 38.3\n0.5000 40.0\n0.9990 56.7")
            convert = QtWidgets.QPushButton("Convertir le tableau collé")
            convert.clicked.connect(self._convert_pasted)
            col.addWidget(self.paste_edit)
            col.addWidget(convert)
            layout.addWidget(box_b, 2)
            self.exp_view = QtWidgets.QTextBrowser()
            self.exp_view.setHtml("<i>Aucune donnée expérimentale chargée.</i>")
            layout.addWidget(self.exp_view, 3)
            self.tabs.addTab(page, "2. Données expérimentales")

        def _browse_experimental(self) -> None:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, "Isotherme expérimentale", "", "Tableurs (*.xlsx *.xls *.csv *.txt);;Tous (*)")
            if path:
                self._load_experimental_file(path)

        def _load_experimental_file(self, path: str) -> None:
            try:
                iso = ExperimentalDataLoader.load_file(path)
            except Exception as exc:
                self._error("Lecture des données expérimentales impossible :\n%s" % exc)
                return
            self.exp_path.setText(path)
            self._set_experimental(iso)

        def _convert_pasted(self) -> None:
            try:
                iso = ExperimentalDataLoader.parse_text(self.paste_edit.toPlainText())
            except Exception as exc:
                self._error("Tableau illisible :\n%s" % exc)
                return
            self._set_experimental(iso)

        def _set_experimental(self, iso: Isotherm) -> None:
            self.experimental = iso
            bet = BETAnalyzer(self._float(self.cross_edit, 0.162)).analyse(iso.p_rel, iso.v_cm3)
            lines = ["<b>%s</b> : %d points d'adsorption, %d de désorption ; P/P0 de %.3g à %.3g."
                     % (iso.label, len(iso.p_rel), len(iso.des_p_rel), iso.p_rel.min(), iso.p_rel.max())]
            if bet.valid:
                lines.append("BET (Rouquerol) : <b>%.1f m2/g</b>, Vm = %.2f cm3 STP/g, C = %.4g, plage "
                             "%.2e – %.2e. %s" % (bet.surface_m2_g, bet.vm_cm3_g, bet.c_constant,
                                                  bet.p_min, bet.p_max, bet.message))
            meta = iso.meta
            if meta:
                lines.append("Métadonnées du fichier : " + ", ".join("%s = %s" % kv for kv in meta.items()))
            gas = str(meta.get("adsorptive", "")).upper()
            if gas in ADSORBATES and gas != self.gas_combo.currentText():
                lines.append("<span style='color:#b00'>Le fichier indique l'adsorbat %s : pensez à le "
                             "sélectionner dans l'onglet Paramètres.</span>" % gas)
            self.exp_view.setHtml("<p>" + "</p><p>".join(lines) + "</p>" +
                                  self._html_table(iso.to_frame()))

        # ======================== Onglet 3 : Paramètres RASPA ========================
        @staticmethod
        def _float(edit, default: float) -> float:
            value = to_float(edit.text())
            return value if np.isfinite(value) else default

        def _build_parameters_tab(self) -> None:
            page = QtWidgets.QWidget()
            outer = QtWidgets.QHBoxLayout(page)
            left, right = QtWidgets.QFormLayout(), QtWidgets.QFormLayout()
            box_left = QtWidgets.QGroupBox("Adsorbat et simulation")
            box_left.setLayout(left)
            box_right = QtWidgets.QGroupBox("Isotherme, analyse et exécution")
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
            auto.setToolTip("Plus petite supercellule compatible avec le rayon de coupure (image minimale)")
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
            left.addRow("Adsorbat :", self.gas_combo)
            left.addRow("", self.gas_info)
            left.addRow("Température :", self.temp_spin)
            left.addRow("Pression de saturation P0 (Pa) :", self.p0_edit)
            left.addRow("Cycles d'équilibration :", self.init_spin)
            left.addRow("Cycles de production :", self.prod_spin)
            left.addRow("Supercellule (a, b, c) :", cell_row)
            left.addRow("Champ de force :", self.ff_combo)
            left.addRow("Charges de la charpente :", self.charge_combo)
            left.addRow("Rayon de coupure :", self.cutoff_spin)

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
            raspa_browse = QtWidgets.QPushButton("Parcourir…")
            raspa_browse.clicked.connect(self._browse_raspa)
            raspa_row.addWidget(raspa_browse)
            self.out_edit = QtWidgets.QLineEdit(
                str(self.qsettings.value("output_base", "") or "") or
                os.path.join(os.path.expanduser("~"), "RASPA_Mercury"))
            out_row = QtWidgets.QHBoxLayout()
            out_row.addWidget(self.out_edit)
            out_browse = QtWidgets.QPushButton("Parcourir…")
            out_browse.clicked.connect(self._browse_output)
            out_row.addWidget(out_browse)
            right.addRow("P/P0 minimal :", self.pmin_edit)
            right.addRow("P/P0 maximal :", self.pmax_edit)
            right.addRow("Nombre de points (échelle log) :", self.npoints_spin)
            right.addRow("Section d'encombrement BET (nm2) :", self.cross_edit)
            right.addRow("P/P0 pour le volume poreux (Gurvich) :", self.gurvich_edit)
            right.addRow("Simulations simultanées :", self.jobs_spin)
            right.addRow("Exécutable RASPA3 :", raspa_row)
            right.addRow("Dossier des résultats :", out_row)
            self.tabs.addTab(page, "3. Paramètres RASPA")

        def _browse_raspa(self) -> None:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Exécutable raspa3", "/usr/local/bin")
            if path:
                self.raspa_edit.setText(path)

        def _browse_output(self) -> None:
            path = QtWidgets.QFileDialog.getExistingDirectory(self, "Dossier des résultats", self.out_edit.text())
            if path:
                self.out_edit.setText(path)

        def _on_gas_changed(self) -> None:
            """Met à jour T, P0, section BET et mode de charges selon l'adsorbat."""
            gas = ADSORBATES[self.gas_combo.currentText()]
            self.temp_spin.setValue(gas["temperature"])
            self.p0_edit.setText("%.0f" % gas["p0"])
            self.cross_edit.setText("%.3f" % gas["cross_section"])
            charged = any(abs(pa["charge"]) > 0 for pa in gas["pseudo_atoms"])
            self.charge_combo.setEnabled(charged)
            self.charge_combo.setCurrentIndex(0 if charged else 2)
            self.gas_info.setText("Modèle : %s. Sonde géométrique : %.2f Å." % (gas["model"], gas["probe_diameter"]))
            if self.geometric is not None:                      # la sonde dépend de l'adsorbat
                self.geometric = None
                self.geometry_view.setHtml("<i>Adsorbat modifié : propriétés géométriques à recalculer.</i>")

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

        # ============================== Onglet 4 : Calcul ==============================
        def _build_run_tab(self) -> None:
            page = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout(page)
            row = QtWidgets.QHBoxLayout()
            self.run_button = QtWidgets.QPushButton("Run RASPA Simulation")
            self.run_button.setMinimumHeight(38)
            self.run_button.clicked.connect(lambda: self._start_analysis(True))
            self.quick_button = QtWidgets.QPushButton("Analyser sans simulation (expérience + géométrie)")
            self.quick_button.clicked.connect(lambda: self._start_analysis(False))
            self.stop_button = QtWidgets.QPushButton("Arrêter")
            self.stop_button.setEnabled(False)
            self.stop_button.clicked.connect(self._stop)
            row.addWidget(self.run_button, 2)
            row.addWidget(self.quick_button, 2)
            row.addWidget(self.stop_button, 1)
            layout.addLayout(row)
            self.progress_bar = QtWidgets.QProgressBar()
            self.progress_bar.setRange(0, 1000)
            self.progress_label = QtWidgets.QLabel("En attente.")
            layout.addWidget(self.progress_bar)
            layout.addWidget(self.progress_label)
            self.log_edit = QtWidgets.QPlainTextEdit()
            self.log_edit.setReadOnly(True)
            layout.addWidget(self.log_edit)
            self.tabs.addTab(page, "4. Calcul")

        def _log(self, text: str) -> None:
            self.log_edit.appendPlainText("[%s] %s" % (_dt.datetime.now().strftime("%H:%M:%S"), text))

        def _progress(self, fraction: float, text: str) -> None:
            self.progress_bar.setValue(int(1000 * max(0.0, min(1.0, fraction))))
            self.progress_label.setText(text)

        def _start_worker(self, job: Callable, on_done: Callable) -> None:
            self.worker = Worker(job, self)
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
            self._log("ERREUR : " + message)
            self.progress_label.setText("Échec — voir le journal.")
            self._error(message[-3000:])

        def _stop(self) -> None:
            if self.runner is not None:
                self._log("Arrêt demandé…")
                self.runner.cancel()

        def _start_analysis(self, with_simulation: bool) -> None:
            """Lance (dans un fil séparé) géométrie, simulation RASPA3, analyse et exports."""
            if self.worker is not None:
                return
            settings = self._collect_settings()
            self.qsettings.setValue("raspa_path", settings.raspa_path)
            self.qsettings.setValue("output_base", self.out_edit.text().strip())
            if with_simulation:
                problems = RaspaInputBuilder(self.structure, settings).validate()
                if not os.path.isfile(settings.raspa_path):
                    problems.append("Exécutable RASPA3 introuvable : '%s'." % settings.raspa_path)
                if problems:
                    self._error("Impossible de lancer la simulation :\n - " + "\n - ".join(problems))
                    return
            elif self.experimental is None:
                self._error("Charger d'abord des données expérimentales (onglet 2).")
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
                progress(1.0, "Analyse des résultats…")
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
            self._log("=== Nouveau calcul : %s ===" % ("simulation RASPA3" if with_simulation else "sans simulation"))
            self._start_worker(job, self._analysis_done)

        def _analysis_done(self, payload) -> None:
            results, written = payload
            self.results = results
            self._show_geometry(results.geometric)
            self.last_html = written.get("html", "")
            self._publish_to_mercury()
            self._progress(1.0, "Terminé. Résultats : %s" % results.settings.work_dir)
            self._show_results()
            self.tabs.setCurrentIndex(4)

        def _publish_to_mercury(self) -> None:
            """Copie le rapport HTML dans le fichier que Mercury affiche à la fin du script."""
            target = self.bridge.output_html_file() if self.bridge else ""
            if target and self.last_html and os.path.isfile(self.last_html):
                try:
                    shutil.copyfile(self.last_html, target)
                    self._log("Rapport HTML transmis à Mercury : %s" % target)
                except OSError as exc:
                    self._log("Rapport HTML non transmis à Mercury : %s" % exc)

        # ============================= Onglet 5 : Résultats =============================
        def _build_results_tab(self) -> None:
            page = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout(page)
            self.result_tabs = QtWidgets.QTabWidget()
            self.figure_tabs = {}
            for key, title in (("isotherme", "Figure 1 — Isothermes"), ("bet_exp", "Figure 2 — BET expérimental"),
                               ("bet_sim", "Figure 3 — BET simulé"), ("surfaces", "Figure 4 — Surfaces")):
                tab = FigureTab()
                self.figure_tabs[key] = tab
                self.result_tabs.addTab(tab, title)
            self.comparison_view = QtWidgets.QTextBrowser()
            self.result_tabs.addTab(self.comparison_view, "Tableau comparatif")
            layout.addWidget(self.result_tabs)
            row = QtWidgets.QHBoxLayout()
            export = QtWidgets.QPushButton("Réexporter (Excel, CSV, PDF)")
            export.clicked.connect(self._export_again)
            folder = QtWidgets.QPushButton("Ouvrir le dossier des résultats")
            folder.clicked.connect(self._open_folder)
            row.addStretch(1)
            row.addWidget(export)
            row.addWidget(folder)
            layout.addLayout(row)
            self.tabs.addTab(page, "5. Résultats")

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
                self._error("Aucun résultat à exporter.")
                return
            try:
                self.results.settings.cross_section = self._float(self.cross_edit, self.results.settings.cross_section)
                self.results.settings.gurvich_p_rel = self._float(self.gurvich_edit, 0.95)
                self.results.analyse()                          # prend en compte section / P/P0 modifiés
                written = ReportBuilder(self.results, self.results.settings.work_dir, self._log).write_all()
                self.last_html = written.get("html", "")
                self._publish_to_mercury()
                self._show_results()
            except Exception as exc:
                self._error("Export impossible :\n%s" % exc)

        def _open_folder(self) -> None:
            folder = self.results.settings.work_dir if self.results else self.out_edit.text()
            if os.path.isdir(folder):
                QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(folder))

        def closeEvent(self, event) -> None:                    # noqa: N802 (nom imposé par Qt)
            """À la fermeture : arrête RASPA3 et laisse un rapport à Mercury."""
            if self.runner is not None:
                self.runner.cancel()
            if self.worker is not None:
                self.worker.wait(8000)
            target = self.bridge.output_html_file() if self.bridge else ""
            if target and not self.last_html:
                try:
                    with open(target, "w", encoding="utf-8") as handle:
                        handle.write("<html><body><h2>%s</h2><p>Session fermée sans calcul.</p></body></html>"
                                     % APP_NAME)
                except OSError:
                    pass
            event.accept()


# =============================================================================
# 14. MODE PAR LOT (sans interface) ET POINT D'ENTRÉE
# =============================================================================
def run_batch(structure: StructureData, args: argparse.Namespace) -> AnalysisResults:
    """Chaîne complète sans interface graphique (tests, scripts, calcul distant)."""
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
        "interne" if args.internal_geometry else "auto")
    results.analyse()
    ReportBuilder(results, settings.work_dir, print).write_all()
    print(results.comparison_frame().to_string(index=False))
    return results


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=APP_NAME, add_help=True)
    parser.add_argument("cif", nargs="?", default="", help="fichier CIF (hors Mercury)")
    parser.add_argument("--batch", action="store_true", help="sans interface graphique")
    parser.add_argument("--exp", default="", help="isotherme expérimentale (xlsx, xls, csv)")
    parser.add_argument("--gas", default="N2", choices=list(ADSORBATES))
    parser.add_argument("--raspa", default="", help="chemin de l'exécutable raspa3")
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--p0", type=float, default=0.0, help="pression de saturation (Pa)")
    parser.add_argument("--init", type=int, default=2000, help="cycles d'équilibration")
    parser.add_argument("--prod", type=int, default=5000, help="cycles de production")
    parser.add_argument("--points", type=int, default=35)
    parser.add_argument("--pmin", type=float, default=1e-6)
    parser.add_argument("--pmax", type=float, default=1.0)
    parser.add_argument("--cells", type=int, nargs=3, default=None, help="supercellule a b c")
    parser.add_argument("--charges", type=int, default=0, choices=[0, 1, 2],
                        help="0 équilibration de charges, 1 charges du CIF, 2 aucune")
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--out", default="", help="dossier des résultats")
    parser.add_argument("--no-sim", action="store_true", help="analyse sans simulation RASPA3")
    parser.add_argument("--keep-guests", action="store_true", help="conserver les espèces invitées")
    parser.add_argument("--internal-geometry", action="store_true",
                        help="forcer le calcul géométrique interne")
    # --- options internes (communication entre les processus du script) ---
    parser.add_argument("--gui-child", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--child-meta", default="", help=argparse.SUPPRESS)
    parser.add_argument("--pore-helper", default="", help=argparse.SUPPRESS)
    return parser


def run_gui_child(args: argparse.Namespace) -> int:
    """Processus d'interface : Qt uniquement, la structure vient du CIF de transfert."""
    meta = {}
    if args.child_meta and os.path.isfile(args.child_meta):
        with open(args.child_meta, "r", encoding="utf-8") as handle:
            meta = json.load(handle)
    if not QT_BINDING:
        print("Aucune liaison Qt utilisable dans le processus d'interface.", file=sys.stderr)
        return 1
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    try:
        structure = StructureData.from_cif_file(args.cif)
    except Exception as exc:
        QtWidgets.QMessageBox.critical(None, APP_NAME, "Structure illisible :\n%s" % exc)
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
    close_ms = os.environ.get("RASPA_POROSITY_TEST_CLOSE_MS")       # tests automatiques seulement
    if close_ms:
        QtCore.QTimer.singleShot(int(close_ms), window.close)
    run_loop = getattr(app, "exec", None) or getattr(app, "exec_")
    run_loop()
    return 0


def launch_gui_child(structure: StructureData, bridge: MercuryBridge, args: argparse.Namespace) -> int:
    """Processus parent : écrit le CIF de transfert et lance l'interface dans un autre processus."""
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
    # L'interface ne doit pas hériter des chemins Qt de Mercury ou de la CSD API :
    # PyQt/PySide doivent charger leurs propres bibliothèques et greffons.
    env = {k: v for k, v in os.environ.items()
           if not (k.startswith("DYLD_") or k in ("QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH",
                                                  "QML2_IMPORT_PATH", "QML_IMPORT_PATH"))}
    code = subprocess.call(command, env=env)
    if code != 0:
        print("L'interface graphique s'est terminée avec le code %d." % code, file=sys.stderr)
    return code


def main(argv: Optional[Sequence[str]] = None) -> int:
    # parse_known_args : Mercury ajoute ses propres arguments, qu'il faut tolérer
    args, _unknown = build_argument_parser().parse_known_args(argv)
    cif_path = args.cif if args.cif and os.path.isfile(args.cif) and args.cif.lower().endswith(".cif") else ""
    if args.pore_helper:                                        # processus auxiliaire (CSD API seule)
        return pore_helper_main(args.pore_helper)
    if args.gui_child:                                          # processus d'interface (Qt seul)
        return run_gui_child(args)
    # ---- processus parent : CSD Python API / Mercury, jamais Qt ----
    try:
        structure, bridge, _temp_cif = load_structure(cif_path)
    except Exception as exc:
        text = str(exc)
        if "Symbol not found" in text or "dlopen" in text:      # conflit Qt / CSD Python API
            print("ERREUR (%s %s) : la CSD Python API ne se charge plus, car une bibliothèque Qt "
                  "installée dans le Python de Mercury entre en conflit avec la sienne.\n\n%s\n\n"
                  "Détail technique : %s" % (APP_NAME, APP_VERSION, qt_fix_instructions(), text[:400]),
                  file=sys.stderr)
        else:
            print("ERREUR (%s %s) : %s" % (APP_NAME, APP_VERSION, text), file=sys.stderr)
        return 1
    if args.batch:
        try:
            run_batch(structure, args)
        except (RaspaError, ValueError, RuntimeError) as exc:
            print("\nERREUR : %s" % exc, file=sys.stderr)
            return 1
        return 0
    if not detect_qt_binding():
        message = "Aucune liaison Qt n'est disponible pour l'interface graphique.\n\n" + qt_fix_instructions()
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
        # Sur macOS, les bibliothèques Qt de la CSD Python API peuvent planter
        # (« QMutex: destroying locked mutex », segmentation fault) pendant la
        # fermeture normale de l'interpréteur, alors que tout le travail est fait.
        # Sortie immédiate, sans phase de nettoyage, pour éviter ce faux plantage.
        os._exit(int(_code or 0))
    sys.exit(_code)
