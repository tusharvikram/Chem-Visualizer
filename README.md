# Chem101 Visualizer

Interactive 3D molecular geometry viewer built for introductory chemistry (VSEPR theory). Select from 21 molecules, rotate them in 3D, and see their VSEPR class, geometry, hybridization, polarity, bond angles, lone pairs, and net dipole rendered in real time.

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## Features

- **21 pre-built molecules** covering linear, bent, trigonal planar, trigonal pyramidal, tetrahedral, trigonal bipyramidal, seesaw, T-shaped, square planar, square pyramidal, and octahedral geometries
- **VSEPR worked out, not just stated** — electron domains, AXₙEₘ class, both geometries, and hybridization, listed in the order the problem is solved
- **Net dipole arrow** computed by summing bond dipoles, the same method taught in class; symmetric molecules show that the dipoles cancel
- **3D interactive viewer** powered by [3Dmol.js](https://3dmol.csb.pitt.edu/) with Ball & Stick and Space Fill modes
- **Lone pair visualization** rendered as pink dummy atoms in chemically accurate positions
- **Hardcoded geometry overrides** for tricky molecules (PCl₅, SF₄, ClF₃, XeF₂, XeF₄, BF₃, SF₆, BrF₅, SO₂) so bond angles are textbook-perfect
- **Info panel** showing molecular formula, average molecular weight, polarity, electron domains, VSEPR class, both geometries, hybridization, and key bond angles
- **Per-molecule legend** listing only the elements actually present, in Jmol colours
- **Pause/resume rotation** so a fixed orientation can be held while explaining a geometry
- **Single-file HTML output** — no web server required, and no CDN at runtime

## Requirements

- Python 3.8+
- [RDKit](https://www.rdkit.org/) (the only third-party dependency; everything else is from the standard library)

Install via conda (recommended):

```bash
conda install -c conda-forge rdkit
```

## Usage

```bash
python chem101_visualizer.py
```

This writes the same page to two places:

- `out/chem_gallery_final.html` — scratch copy, gitignored
- `public/index.html` — the deployable site

Open either in any modern browser, or serve the folder:

```bash
python -m http.server 4321 --directory public
```

3Dmol.js is vendored in `public/vendor/`, so the viewer works with no internet
connection. Only the Inter webfont is remote, and it falls back to a system font.

> **Note on bond angles:** geometries use *idealized* VSEPR angles. Real molecules
> deviate where lone pairs compress bonds — ClF₃ is 87.5° in practice, not 90°.

## Deployment

The site is static, so RDKit never runs on the server — the Python script is a
build step and `public/` is committed. Deploy with:

```bash
vercel --prod
```

After editing molecules, re-run the build so `public/index.html` is regenerated,
then redeploy.

## Molecule Library

| Molecule | VSEPR Class | Molecular Geometry | Electron Geometry | Hybridization | Polarity |
|---|---|---|---|---|---|
| H₂O | AX₂E₂ | Bent (V-shaped) | Tetrahedral | sp³ | Polar |
| NH₃ | AX₃E | Trigonal Pyramidal | Tetrahedral | sp³ | Polar |
| CH₄ | AX₄ | Tetrahedral | Tetrahedral | sp³ | Non-Polar |
| CO₂ | AX₂ | Linear | Linear | sp | Non-Polar |
| C₂H₄ | AX₃ | Trigonal Planar | Trigonal Planar | sp² | Non-Polar |
| C₂H₂ | AX₂ | Linear | Linear | sp | Non-Polar |
| CH₂O | AX₃ | Trigonal Planar | Trigonal Planar | sp² | Polar |
| CHCl₃ | AX₄ | Tetrahedral | Tetrahedral | sp³ | Polar |
| PCl₅ | AX₅ | Trigonal Bipyramidal | Trigonal Bipyramidal | sp³d | Non-Polar |
| SF₄ | AX₄E | Seesaw | Trigonal Bipyramidal | sp³d | Polar |
| ClF₃ | AX₃E₂ | T-Shaped | Trigonal Bipyramidal | sp³d | Polar |
| XeF₂ | AX₂E₃ | Linear | Trigonal Bipyramidal | sp³d | Non-Polar |
| XeF₄ | AX₄E₂ | Square Planar | Octahedral | sp³d² | Non-Polar |
| BF₃ | AX₃ | Trigonal Planar | Trigonal Planar | sp² | Non-Polar |
| C₆H₆ | AX₃ | Trigonal Planar (Ring) | Trigonal Planar | sp² | Non-Polar |
| C₂H₅OH | AX₂E₂ | Tetrahedral (C) / Bent (O) | Tetrahedral | sp³ | Polar |
| C₃H₆O | AX₃ | Trigonal Planar (at C=O) | Trigonal Planar | sp² | Polar |
| H₂O₂ | AX₂E₂ | Open Book (Non-planar) | Tetrahedral | sp³ | Polar |
| SF₆ | AX₆ | Octahedral | Octahedral | sp³d² | Non-Polar |
| BrF₅ | AX₅E | Square Pyramidal | Octahedral | sp³d² | Polar |
| SO₂ | AX₂E | Bent (V-shaped) | Trigonal Planar | sp² | Polar |

## How It Works

1. Each molecule is built from a SMILES string using RDKit
2. 3D coordinates are generated via `EmbedMolecule` + UFF optimization, or overridden with hardcoded coordinates for geometries RDKit handles poorly
3. Lone pairs are calculated from bond vectors and added as dummy atoms
4. AXE notation, steric number and hybridization are derived from one `(bonds, lone pairs)` pair per molecule in `VSEPR_CENTERS`, so those three fields cannot contradict each other
5. The net dipole is computed by summing one vector per bond, pointing at the more electronegative atom and scaled by the electronegativity difference
6. The build compares that computed dipole against each molecule's stated polarity and reports any disagreement — a wrong coordinate set usually breaks a cancellation, so this doubles as a check on the geometry
7. Everything is serialized into a single self-contained HTML file with embedded JSON data and 3Dmol.js for rendering

## License

Released under the MIT License.
