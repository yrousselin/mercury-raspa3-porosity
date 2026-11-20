# RASPA3 Porosity Analyser pour Mercury

*[English version: README.md](README.md)*

Outil Python pour **Mercury** (CCDC) qui compare la porosité **mesurée** et **simulée** d'un matériau poreux (MOF, ZIF, COF…). Il prend la structure ouverte dans Mercury, lance des simulations d'adsorption de gaz avec **RASPA3**, calcule les surfaces BET et les descripteurs géométriques des pores, puis met le tout en regard de votre isotherme expérimentale.

Version 1.4.

## Ce que fait l'outil

1. Récupère la structure ouverte dans Mercury, ou n'importe quel fichier CIF, et construit un modèle P1 propre. Le désordre est résolu, et les molécules de solvant non liées à la charpente peuvent être supprimées.
2. Écrit les fichiers d'entrée de RASPA3 et lance une **simulation GCMC par point de pression** de l'isotherme. Plusieurs points peuvent tourner en parallèle.
3. Calcule la **surface BET** des isothermes mesurée et simulée. La plage de pression est choisie automatiquement et vérifiée par les quatre **critères de Rouquerol**.
4. Calcule le **volume poreux** par la règle de Gurvich.
5. Calcule les **propriétés géométriques** : surface accessible, volume poreux, diamètre de pore limitant (PLD) et diamètre de la plus grande cavité (LCD). Il utilise le Pore Analyser de la CSD quand c'est possible, sinon une méthode interne documentée.
6. Compare expérience, simulation et géométrie dans un tableau accompagné de figures, et exporte **Excel, CSV, PDF, PNG et un rapport HTML** affiché dans Mercury.

**Gaz disponibles :** N₂ (77 K), Ar (87 K), CO₂ (273 K).
**Données expérimentales :** `.xlsx`, `.xls`, `.csv`, ou tableau collé dans la fenêtre. Les exports BELSORP/BELMaster sont reconnus.

La fenêtre comporte cinq onglets : *Structure*, *Experimental data* (données expérimentales), *RASPA settings* (paramètres RASPA), *Calculation* (calcul), *Results* (résultats). **L'interface, les messages et les rapports sont en anglais.**

## Prérequis

| Quoi | Pourquoi | Où le trouver |
|---|---|---|
| **Mercury avec la CSD Python API** | lecture de la structure, Pore Analyser | [CCDC](https://www.ccdc.cam.ac.uk/) (licence CSD nécessaire) |
| **RASPA3** | simulations d'adsorption (GCMC) | [site de RASPA3](https://iraspa.github.io/RASPA3/index.html) |
| numpy, pandas, matplotlib, openpyxl | calculs, figures, export Excel | `pip` (voir plus bas) |
| **PyQt6** (ou PyQt5, PySide6, PySide2) | la fenêtre graphique | `pip`, dans un **dossier séparé** (voir plus bas) |
| scipy, xlrd *(facultatifs)* | calcul géométrique plus rapide, lecture des anciens `.xls` | `pip` |

## Installation

### Étape 1 : installer RASPA3

La méthode la plus simple sur **macOS, Windows et Linux** passe par conda (Miniconda ou Anaconda). conda-forge fournit des versions pour les trois systèmes.

```
conda create -n raspa -c conda-forge raspa3 -y
```

Repérez ensuite le chemin du programme `raspa3`. Vous en aurez besoin plus tard :

- macOS : `conda activate raspa` puis `which raspa3`
- Windows (Anaconda Prompt) : `conda activate raspa` puis `where raspa3`

Sur macOS, vous pouvez aussi utiliser l'installeur (`.pkg`) de la [page des versions de RASPA3](https://github.com/iRASPA/raspa3/releases). Le programme se trouve alors en général dans `/usr/local/bin/raspa3`, où l'outil le trouve tout seul.

### Étape 2 : ajouter les modules au Python de Mercury

L'outil fonctionne avec le Python de la CSD Python API. Emplacements habituels :

- macOS : `~/CCDC/ccdc-software/csd-python-api/miniconda/bin/python`
- Windows : cherchez `python.exe` dans le dossier `csd-python-api` de votre installation CCDC.

Dans un terminal, remplacez `CSD_PYTHON` ci-dessous par ce chemin (entre guillemets s'il contient des espaces) :

```
CSD_PYTHON -m pip install --target ~/RASPA_Mercury/qt_libs PyQt6
CSD_PYTHON -m pip install pandas matplotlib openpyxl xlrd scipy
```

Sous Windows, écrivez la première commande avec `%USERPROFILE%\RASPA_Mercury\qt_libs` au lieu de `~/RASPA_Mercury/qt_libs`.

> ⚠️ **N'installez pas PyQt directement dans le Python de Mercury.** Sur macOS, il entre en conflit avec les bibliothèques Qt de la CSD Python API, et les scripts Mercury ne fonctionnent plus. C'est pour cela que PyQt va dans le dossier séparé `RASPA_Mercury/qt_libs`.

**Plus simple :** lancez l'outil une fois (étape 3). S'il manque un module, il affiche les **commandes exactes** à copier dans le terminal pour votre ordinateur.

### Étape 3 : ajouter le script à Mercury

1. Placez `raspa3_mercury_porosity.py` dans un dossier réservé à vos scripts Mercury.
2. Dans Mercury, ouvrez le menu **CSD Python API** → **Options…** → **Add Location**, puis choisissez ce dossier.
3. L'outil apparaît alors dans le menu **CSD Python API** (selon la version de Mercury : **Calculate → Python Scripts**).

## Utilisation

1. Ouvrez la structure dans Mercury et lancez **CSD Python API → raspa3_mercury_porosity**. La fenêtre s'ouvre.
2. Onglet **Structure** : vérifiez le modèle P1. Les options de désordre et de solvant sont activées par défaut.
3. Onglet **Experimental data** : chargez votre fichier d'isotherme ou collez un tableau (P/P0 et quantité adsorbée en cm³ STP/g).
4. Onglet **RASPA settings** : choisissez le gaz, vérifiez le nombre de cycles et indiquez le chemin de `raspa3` la première fois. Il est mémorisé ensuite. Le bouton **Auto** choisit la plus petite supercellule permise par le rayon de coupure.
5. Onglet **Calculation** : cliquez sur **Run RASPA Simulation**. Une barre de progression indique le temps restant. Vous pouvez aussi cliquer sur **Analyse without simulation** pour une comparaison rapide entre expérience et géométrie.
6. Onglet **Results** : figures et tableau comparatif. Tous les fichiers sont enregistrés dans `~/RASPA_Mercury/<structure>_<gaz>_<date>/` : `_results.xlsx`, `_report.pdf`, `_report.html`, fichiers CSV et figures PNG.

**Sans Mercury (mode par lot, pour des tests ou une machine distante) :**

```
python raspa3_mercury_porosity.py structure.cif --batch --exp mesure.xlsx --gas N2 --raspa /chemin/vers/raspa3
python raspa3_mercury_porosity.py --help      # toutes les options
```

## Choix de modélisation (pour interpréter les résultats)

- Charpente rigide ; paramètres de Lennard-Jones UFF ou DREIDING par élément ; règles de mélange de Lorentz–Berthelot ; potentiel tronqué et décalé à 12 Å.
- Adsorbats : modèles TraPPE pour N₂ et CO₂ (3 sites chargés) ; modèle de Lennard-Jones à un site pour Ar.
- Charges de la charpente : équilibration de charges de RASPA3, charges lues dans le CIF, ou aucune charge.
- La simulation donne l'adsorption **absolue**, alors que l'expérience mesure l'adsorption **d'excès**. À 77 ou 87 K et pour P ≤ P₀, la différence est négligeable.

## Limites

- Les simulations peuvent durer plusieurs heures, selon la taille de la structure, le nombre de points et le nombre de cycles.
- Les charpentes flexibles et la chimisorption ne sont pas modélisées.
- Développé et utilisé sur macOS. Windows devrait fonctionner, mais n'a pas encore été testé complètement. Vos retours sont les bienvenus.

## Crédits et références

- **RASPA3 :** « RASPA3: A Monte Carlo code for computing adsorption and diffusion in nanoporous materials and thermodynamics properties of fluids », *J. Chem. Phys.* **2024**, 161, 114106. <https://iraspa.github.io/RASPA3/>
- **Critères de cohérence BET :** J. Rouquerol, P. Llewellyn, F. Rouquerol, *Stud. Surf. Sci. Catal.* **2007**, 160, 49.
- **Surface accessible :** T. Düren *et al.*, *J. Phys. Chem. C* **2007**, 111, 15350. **Volume occupable par la sonde :** D. Ongari *et al.*, *Langmuir* **2017**.
- **Champs de force :** UFF (Rappé *et al.*, *J. Am. Chem. Soc.* **1992**, 114, 10024) ; DREIDING (Mayo *et al.*, *J. Phys. Chem.* **1990**, 94, 8897) ; TraPPE (Potoff & Siepmann, 2001).
- L'outil utilise la CSD Python API du Cambridge Crystallographic Data Centre (CCDC).
- Écrit par Yoann Rousselin avec l'aide de Claude (Anthropic) comme assistant de programmation.

## Historique des versions

- **1.4** (novembre 2026) : tout l'outil passe en anglais (interface, messages, rapports, feuilles et colonnes Excel, noms de fichiers, commentaires du code). Les fichiers produits changent de nom en conséquence, par exemple `_results.xlsx` et `_report.pdf` au lieu de `_resultats.xlsx` et `_rapport.pdf`.
- **1.3** (novembre 2026) : première version publique.

## Licence

MIT. Voir [LICENSE](LICENSE). RASPA3 est distribué sous sa propre licence.

## Auteur

**Yoann Rousselin**, ingénieur de recherche CNRS (cristallographie), ICMUB, Dijon.
Questions ou signalement d'un problème : ouvrez une *Issue* sur cette page.
