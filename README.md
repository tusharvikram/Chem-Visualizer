# Chem101 Visualizer

Interactive 3D molecular geometry viewer built for introductory chemistry (VSEPR theory). Select from 18 molecules, rotate them in 3D, and see their geometry, polarity, bond angles, and lone pairs rendered in real time.

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## Features

- **18 pre-built molecules** covering linear, bent, trigonal planar, tetrahedral, trigonal bipyramidal, seesaw, T-shaped, square planar, and octahedral geometries
- **3D interactive viewer** powered by [3Dmol.js](https://3dmol.csb.pitt.edu/) with Ball & Stick and Space Fill modes
- **Lone pair visualization** rendered as pink dummy atoms in chemically accurate positions
- **Hardcoded geometry overrides** for tricky molecules (PCl₅, SF₄, ClF₃, XeF₂, XeF₄, BF₃) so bond angles are textbook-perfect
- **Info panel** showing molecular formula, average molecular weight, polarity, molecular geometry, electron geometry, and key bond angles
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

Released under the MIT License.
