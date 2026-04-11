# Chem101 Visualizer

Interactive 3D molecular geometry viewer built for introductory chemistry (VSEPR theory). Select from 18 molecules, rotate them in 3D, and see their geometry, polarity, bond angles, and lone pairs rendered in real time.

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## Features

- **18 pre-built molecules** covering linear, bent, trigonal planar, tetrahedral, trigonal bipyramidal, seesaw, T-shaped, square planar, and octahedral geometries
- **3D interactive viewer** powered by [3Dmol.js](https://3dmol.csb.pitt.edu/) with Ball & Stick and Space Fill modes
- **Lone pair visualization** rendered as pink dummy atoms in chemically accurate positions
- **Hardcoded geometry overrides** for tricky molecules (PCl₅, SF₄, ClF₃, XeF₂, XeF₄, BF₃) so bond angles are textbook-perfect
- **Info panel** showing molecular formula, weight, polarity, molecular geometry, electron geometry, and key bond angles
- **Single-file HTML output** with no server required

## Requirements

```
rdkit
```

Install via conda (recommended):

```bash
conda install -c conda-forge rdkit
```

## Usage

```bash
python chem_visualizer.py
```

This generates `out/chem_gallery_final.html`. Open it in any modern browser.

## Molecule Library

| Molecule | Molecular Geometry | Electron Geometry | Polarity |
|---|---|---|---|
| H₂O | Bent | Tetrahedral | Polar |
| NH₃ | Trigonal Pyramidal | Tetrahedral | Polar |
| CH₄ | Tetrahedral | Tetrahedral | Non-Polar |
| CO₂ | Linear | Linear | Non-Polar |
| C₂H₄ | Trigonal Planar | Trigonal Planar | Non-Polar |
| C₂H₂ | Linear | Linear | Non-Polar |
| CH₂O | Trigonal Planar | Trigonal Planar | Polar |
| CHCl₃ | Tetrahedral | Tetrahedral | Polar |
| PCl₅ | Trigonal Bipyramidal | Trigonal Bipyramidal | Non-Polar |
| SF₄ | Seesaw | Trigonal Bipyramidal | Polar |
| ClF₃ | T-Shaped | Trigonal Bipyramidal | Polar |
| XeF₂ | Linear | Trigonal Bipyramidal | Non-Polar |
| XeF₄ | Square Planar | Octahedral | Non-Polar |
| BF₃ | Trigonal Planar | Trigonal Planar | Non-Polar |
| C₆H₆ | Trigonal Planar (Ring) | Trigonal Planar | Non-Polar |
| C₂H₅OH | Tetrahedral / Bent | Tetrahedral | Polar |
| C₃H₆O | Trigonal Planar (at C=O) | Trigonal Planar | Polar |
| H₂O₂ | Open Book | Tetrahedral | Polar |

## How It Works

1. Each molecule is built from a SMILES string using RDKit
2. 3D coordinates are generated via `EmbedMolecule` + UFF optimization, or overridden with hardcoded coordinates for geometries RDKit handles poorly
3. Lone pairs are calculated from bond vectors and added as dummy atoms
4. Everything is serialized into a single self-contained HTML file with embedded JSON data and 3Dmol.js for rendering

## License

MIT
