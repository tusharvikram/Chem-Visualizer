"""
Chem 101 Visualizer - Ultimate Edition (Perfect Geometry)
---------------------------------------------------------
Updates:
- UI UPDATE: Dropdown names now show "Formula (Name)" format (e.g., H₂O (Water)).
- GEOMETRY FIX: Updated MANUAL_COORDS keys to match new naming convention.
- FORMULA FIX: Uses the formatted formula from the name instead of raw RDKit hill notation.
- POLARITY FIX: Hardcoded accurate polarity for all 18 molecules.
- PCl5 FIX: Corrected coordinate axis to prevent flattening.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Geometry import Point3D
import math
import json
import os

# ============================================================
# 1. DATABASE MOLECOLE
# ============================================================

MOLECULES = [
    # --- BASICS ---
    {
        "name": "H₂O (Water)", 
        "smiles": "O", 
        "desc": "Bent (V-shaped)", 
        "geo_el": "Tetrahedral", 
        "pol": "Polar",
        "fact": "Universal solvent essential for life. 💧"
    },
    {
        "name": "NH₃ (Ammonia)", 
        "smiles": "N", 
        "desc": "Trigonal Pyramidal", 
        "geo_el": "Tetrahedral", 
        "pol": "Polar",
        "fact": "Common in fertilizers and cleaners. 🌿"
    },
    {
        "name": "CH₄ (Methane)", 
        "smiles": "C", 
        "desc": "Tetrahedral", 
        "geo_el": "Tetrahedral", 
        "pol": "Non-Polar",
        "fact": "Main component of natural gas. 🔥"
    },
    
    # --- DOUBLE/TRIPLE BONDS ---
    {
        "name": "CO₂ (Carbon Dioxide)", 
        "smiles": "O=C=O", 
        "desc": "Linear", 
        "geo_el": "Linear", 
        "pol": "Non-Polar",
        "fact": "Product of respiration. 🌬️"
    },
    {
        "name": "C₂H₄ (Ethene)", 
        "smiles": "C=C", 
        "desc": "Trigonal Planar", 
        "geo_el": "Trigonal Planar", 
        "pol": "Non-Polar",
        "fact": "Plant hormone that ripens fruit. 🍌"
    },
    {
        "name": "C₂H₂ (Ethyne)", 
        "smiles": "C#C", 
        "desc": "Linear", 
        "geo_el": "Linear", 
        "pol": "Non-Polar",
        "fact": "Used in welding torches (Acetylene). 💥"
    },

    # --- FUNCTIONAL GROUPS & EXCEPTIONS ---
    {
        "name": "CH₂O (Formaldehyde)",
        "smiles": "C=O",
        "desc": "Trigonal Planar", 
        "geo_el": "Trigonal Planar", 
        "pol": "Polar",
        "fact": "Simplest aldehyde, used as a preservative. 🧪"
    },
    {
        "name": "CHCl₃ (Chloroform)",
        "smiles": "ClC(Cl)Cl",
        "desc": "Tetrahedral", 
        "geo_el": "Tetrahedral", 
        "pol": "Polar",
        "fact": "Historically used as an anesthetic. 😴"
    },
    {
        "name": "PCl₅ (Phosphorus Pentachloride)",
        "smiles": "ClP(Cl)(Cl)(Cl)Cl",
        "desc": "Trigonal Bipyramidal", 
        "geo_el": "Trigonal Bipyramidal", 
        "pol": "Non-Polar",
        "fact": "Example of Expanded Octet (10 electrons on P). 📐"
    },
    {
        "name": "SF₄ (Sulfur Tetrafluoride)",
        "smiles": "FS(F)(F)F",
        "desc": "Seesaw", 
        "geo_el": "Trigonal Bipyramidal", 
        "pol": "Polar",
        "fact": "Classic VSEPR exception (1 lone pair). ⚖️"
    },
    {
        "name": "ClF₃ (Chlorine Trifluoride)",
        "smiles": "Cl(F)(F)F",
        "desc": "T-Shaped", 
        "geo_el": "Trigonal Bipyramidal", 
        "pol": "Polar",
        "fact": "Extremely reactive interhalogen. 🧨"
    },
    {
        "name": "XeF₂ (Xenon Difluoride)",
        "smiles": "F[Xe]F",
        "desc": "Linear", 
        "geo_el": "Trigonal Bipyramidal", 
        "pol": "Non-Polar",
        "fact": "Noble gases CAN form bonds! (3 Lone Pairs). 🎈"
    },
    {
        "name": "XeF₄ (Xenon Tetrafluoride)", 
        "smiles": "F[Xe](F)(F)F", 
        "desc": "Square Planar", 
        "geo_el": "Octahedral", 
        "pol": "Non-Polar",
        "fact": "Noble gas compound with 2 lone pairs. 🧊"
    },
    {
        "name": "BF₃ (Boron Trifluoride)", 
        "smiles": "FB(F)F", 
        "desc": "Trigonal Planar", 
        "geo_el": "Trigonal Planar", 
        "pol": "Non-Polar",
        "fact": "Incomplete octet (only 6 valence electrons). ⚠️"
    },
    {
        "name": "C₆H₆ (Benzene)",
        "smiles": "c1ccccc1",
        "desc": "Trigonal Planar (Ring)", 
        "geo_el": "Trigonal Planar", 
        "pol": "Non-Polar",
        "fact": "Aromatic ring with delocalized electrons. 🍩"
    },
    {
        "name": "C₂H₅OH (Ethanol)",
        "smiles": "CCO",
        "desc": "Tetrahedral (C) / Bent (O)", 
        "geo_el": "Tetrahedral", 
        "pol": "Polar",
        "fact": "The type of alcohol found in beverages. 🍷"
    },
    {
        "name": "C₃H₆O (Acetone)",
        "smiles": "CC(=O)C",
        "desc": "Trigonal Planar (at C=O)", 
        "geo_el": "Trigonal Planar", 
        "pol": "Polar",
        "fact": "Common solvent and nail polish remover. 💅"
    },
    {
        "name": "H₂O₂ (Hydrogen Peroxide)",
        "smiles": "OO",
        "desc": "Open Book (Non-planar)",
        "geo_el": "Tetrahedral",
        "pol": "Polar",
        "fact": "Reactive oxidizer with a twisted structure. 📘"
    },

    # --- OCTAHEDRAL FAMILY & THE SECOND 'BENT' ---
    {
        "name": "SF₆ (Sulfur Hexafluoride)",
        "smiles": "FS(F)(F)(F)(F)F",
        "desc": "Octahedral",
        "geo_el": "Octahedral",
        "pol": "Non-Polar",
        "fact": "Six bonds, no lone pairs — the textbook octahedron. 🎲"
    },
    {
        "name": "BrF₅ (Bromine Pentafluoride)",
        "smiles": "FBr(F)(F)(F)F",
        "desc": "Square Pyramidal",
        "geo_el": "Octahedral",
        "pol": "Polar",
        "fact": "One lone pair turns the octahedron into a pyramid. 🔺"
    },
    {
        "name": "SO₂ (Sulfur Dioxide)",
        "smiles": "O=S=O",
        "desc": "Bent (V-shaped)",
        "geo_el": "Trigonal Planar",
        "pol": "Polar",
        "fact": "Also bent — but from trigonal planar, so ~119° not 104.5°. Compare with H₂O. 🌋"
    }
]

# ============================================================
# 1b. VSEPR CLASSIFICATION
# ============================================================

# The atom each molecule's `desc` describes: (bonds, lone pairs, label).
# AXE notation, steric number and hybridisation are all derived from these two
# integers, so a typo can make one molecule wrong but can never leave the three
# fields contradicting each other.
VSEPR_CENTERS = {
    "H₂O (Water)":                      (2, 2, "O"),
    "NH₃ (Ammonia)":                    (3, 1, "N"),
    "CH₄ (Methane)":                    (4, 0, "C"),
    "CO₂ (Carbon Dioxide)":             (2, 0, "C"),
    "C₂H₄ (Ethene)":                    (3, 0, "C"),
    "C₂H₂ (Ethyne)":                    (2, 0, "C"),
    "CH₂O (Formaldehyde)":              (3, 0, "C"),
    "CHCl₃ (Chloroform)":               (4, 0, "C"),
    "PCl₅ (Phosphorus Pentachloride)":  (5, 0, "P"),
    "SF₄ (Sulfur Tetrafluoride)":       (4, 1, "S"),
    "ClF₃ (Chlorine Trifluoride)":      (3, 2, "Cl"),
    "XeF₂ (Xenon Difluoride)":          (2, 3, "Xe"),
    "XeF₄ (Xenon Tetrafluoride)":       (4, 2, "Xe"),
    "BF₃ (Boron Trifluoride)":          (3, 0, "B"),
    "C₆H₆ (Benzene)":                   (3, 0, "each C"),
    "C₂H₅OH (Ethanol)":                 (2, 2, "O"),
    "C₃H₆O (Acetone)":                  (3, 0, "the C=O carbon"),
    "H₂O₂ (Hydrogen Peroxide)":         (2, 2, "each O"),
    "SF₆ (Sulfur Hexafluoride)":        (6, 0, "S"),
    "BrF₅ (Bromine Pentafluoride)":     (5, 1, "Br"),
    "SO₂ (Sulfur Dioxide)":             (2, 1, "S"),
}

SUBSCRIPTS = "₀₁₂₃₄₅₆₇₈₉"

# Steric number -> hybridisation, as taught in introductory chemistry. sp³d and
# sp³d² are the textbook labels for expanded octets; modern computational work
# finds d-orbital participation to be negligible, but students are examined on
# these, so they are what the panel reports.
HYBRIDIZATION = {2: "sp", 3: "sp²", 4: "sp³", 5: "sp³d", 6: "sp³d²"}

def vsepr_facts(name):
    if name not in VSEPR_CENTERS: return {}
    x, e, center = VSEPR_CENTERS[name]
    axe = "AX" + (SUBSCRIPTS[x] if x > 1 else "")
    if e: axe += "E" + (SUBSCRIPTS[e] if e > 1 else "")
    return {
        "axe": axe,
        "steric": x + e,
        "hybrid": HYBRIDIZATION.get(x + e, "—"),
        "center": center,
        "lone_pairs": e,
    }

# ============================================================
# 2. MATH HELPERS (Vector Logic)
# ============================================================

def v_sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def v_add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def v_scale(v, s): return (v[0]*s, v[1]*s, v[2]*s)
def v_dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def v_cross(a, b): return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])
def v_len(v): return math.sqrt(v_dot(v, v))
def v_norm(v):
    l = v_len(v)
    return (0,0,0) if l==0 else (v[0]/l, v[1]/l, v[2]/l)

def angle_deg(pos_a, pos_b, pos_c) -> float:
    ba = v_sub(pos_a, pos_b)
    bc = v_sub(pos_c, pos_b)
    mba, mbc = v_len(ba), v_len(bc)
    if mba == 0 or mbc == 0: return 0.0
    val = v_dot(ba, bc) / (mba * mbc)
    return math.degrees(math.acos(max(-1.0, min(1.0, val))))

# ============================================================
# 3. MANUAL COORDINATES (Force Perfect Geometry)
# ============================================================
# Hardcoding coordinates for tricky molecules to ensure they look perfect.
# Central atom is always at (0,0,0).
# Units in Angstroms (approx).

MANUAL_COORDS = {
    # Trigonal Bipyramidal (Axial 180, Equatorial 120)
    "PCl₅ (Phosphorus Pentachloride)": [
        (0.0, 0.0, 0.0), # P
        (0.0, 2.14, 0.0), # Cl Axial (Y-axis)
        (0.0, -2.14, 0.0), # Cl Axial (Y-axis)
        (2.02, 0.0, 0.0), # Cl Eq (X-axis)
        (-1.01, 0.0, 1.75), # Cl Eq (XZ plane - Fixed)
        (-1.01, 0.0, -1.75) # Cl Eq (XZ plane - Fixed)
    ],
    # Octahedral (90 degrees everywhere)
    "SF₆ (Sulfur Hexafluoride)": [
        (0.0, 0.0, 0.0), # S
        (1.6, 0.0, 0.0), (-1.6, 0.0, 0.0),
        (0.0, 1.6, 0.0), (0.0, -1.6, 0.0),
        (0.0, 0.0, 1.6), (0.0, 0.0, -1.6)
    ],
    # Square Planar (Atoms in XY plane)
    "XeF₄ (Xenon Tetrafluoride)": [
        (0.0, 0.0, 0.0), # Xe
        (1.95, 0.0, 0.0), (-1.95, 0.0, 0.0),
        (0.0, 1.95, 0.0), (0.0, -1.95, 0.0)
    ],
    # Linear (3 Lone pairs, atoms on Z axis)
    "XeF₂ (Xenon Difluoride)": [
        (0.0, 0.0, 0.0), # Xe
        (0.0, 0.0, 2.0), # F
        (0.0, 0.0, -2.0) # F
    ],
    # Seesaw (TBP with 1 missing equatorial)
    "SF₄ (Sulfur Tetrafluoride)": [
        (0.0, 0.0, 0.0), # S
        (0.0, 0.0, 1.65), # F Axial
        (0.0, 0.0, -1.65), # F Axial
        (1.55, 0.0, 0.0), # F Eq
        (-0.77, 1.34, 0.0) # F Eq (120 deg from other Eq)
        # Missing Eq is at (-0.77, -1.34, 0) -> This is where LP goes
    ],
    # T-Shaped (TBP with 2 missing equatorial)
    "ClF₃ (Chlorine Trifluoride)": [
        (0.0, 0.0, 0.0), # Cl
        (0.0, 0.0, 1.7), # F Axial
        (0.0, 0.0, -1.7), # F Axial
        (1.7, 0.0, 0.0)  # F Eq
        # Missing Eqs are at +/- 120 -> LPs
    ],
    # Trigonal Planar (Perfect Flatness)
    "BF₃ (Boron Trifluoride)": [
        (0.0, 0.0, 0.0), # B
        (1.3, 0.0, 0.0),
        (-0.65, 1.12, 0.0),
        (-0.65, -1.12, 0.0)
    ],
    # Octahedral (6 bonds, no lone pairs) — all angles exactly 90/180
    "SF₆ (Sulfur Hexafluoride)": [
        (0.0, 0.0, 0.0), # S
        (1.6, 0.0, 0.0), (-1.6, 0.0, 0.0),
        (0.0, 1.6, 0.0), (0.0, -1.6, 0.0),
        (0.0, 0.0, 1.6), (0.0, 0.0, -1.6)
    ],
    # Square Pyramidal (octahedron with one vertex replaced by a lone pair).
    # Apex on +Z, square base on the XY plane, so the LP sits on -Z.
    "BrF₅ (Bromine Pentafluoride)": [
        (0.0, 0.0, 0.0),  # Br
        (0.0, 0.0, 1.75), # F apical
        (1.75, 0.0, 0.0), (-1.75, 0.0, 0.0),
        (0.0, 1.75, 0.0), (0.0, -1.75, 0.0)
    ],
    # Bent from TRIGONAL PLANAR (~119°), not from tetrahedral like water.
    # Bisector on +Y so the single lone pair points along -Y.
    "SO₂ (Sulfur Dioxide)": [
        (0.0, 0.0, 0.0),      # S
        (1.229, 0.726, 0.0),  # O
        (-1.229, 0.726, 0.0)  # O
    ]
}

def apply_manual_geometry(mol, name):
    if name not in MANUAL_COORDS: return mol
    
    conf = mol.GetConformer()
    coords = MANUAL_COORDS[name]
    
    atoms = list(mol.GetAtoms())
    # Identify the Central Atom (Highest Degree)
    sorted_atoms = sorted(atoms, key=lambda x: x.GetDegree(), reverse=True)
    central_atom = sorted_atoms[0]
    central_idx = central_atom.GetIdx()

    # Assign Center Coordinate (0,0,0)
    conf.SetAtomPosition(central_idx, Point3D(*coords[0]))
    
    # Assign Neighbors to the rest of the coordinates
    neighbors = central_atom.GetNeighbors()
    neighbor_coords = coords[1:]
    
    for i, neighbor in enumerate(neighbors):
        if i < len(neighbor_coords):
            n_idx = neighbor.GetIdx()
            conf.SetAtomPosition(n_idx, Point3D(*neighbor_coords[i]))
            
    return mol

# ============================================================
# 4. CORE LOGIC
# ============================================================

def get_molecule_data(entry):
    mol = Chem.MolFromSmiles(entry['smiles'])
    strict = mol is not None

    if not strict:
        # Cl is capped at valence 1 in RDKit's model, so ClF₃ is rejected outright
        # (P and S allow hypervalence, which is why PCl₅ and SF₄ parse fine).
        # Re-parse without the valence check and sanitize everything else.
        mol = Chem.MolFromSmiles(entry['smiles'], sanitize=False)
        if not mol: return None
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
        )
        if entry['name'] not in MANUAL_COORDS:
            # AddHs and EmbedMolecule both re-run the valence check, so there is no
            # way to generate coordinates for these without hardcoded ones.
            return None

    if strict:
        mol = Chem.AddHs(mol)

    # 1. Generate Basic 3D
    if strict:
        AllChem.EmbedMolecule(mol)
    else:
        # These molecules carry no hydrogens and always have manual coordinates,
        # so seed an empty conformer for apply_manual_geometry to fill in.
        mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()), assignId=True)

    # 2. OVERRIDE with Manual Coordinates if available (The "Fix")
    if entry['name'] in MANUAL_COORDS:
        mol = apply_manual_geometry(mol, entry['name'])
    else:
        # Only optimize if NOT manual. Manual coords are already perfect.
        try: AllChem.UFFOptimizeMolecule(mol)
        except: pass

    # --- Analysis ---
    analysis = analyze_properties(mol)
    
    # OVERRIDE analysis with corrected data
    # Parse formula from name e.g. "H₂O (Water)" -> "H₂O"
    display_formula = entry['name'].split(' ')[0]
    
    analysis.update({
        'formula': display_formula,  # Use pretty formula
        'polarity': entry['pol'],    # Use corrected polarity
        'desc': entry['desc'],
        'geo_el': entry['geo_el'],
        'fact': entry['fact']
    })
    analysis.update(vsepr_facts(entry['name']))
    analysis['dipole'] = compute_dipole(mol)

    # --- Add Visual Lone Pairs ---
    mol_viz = add_lone_pairs(mol)
    analysis['has_lp'] = any(a.GetAtomicNum() == 0 for a in mol_viz.GetAtoms())

    return {
        "molblock": Chem.MolToMolBlock(mol_viz),
        "analysis": analysis
    }

def add_lone_pairs(mol):
    rwmol = Chem.RWMol(mol)
    conf = rwmol.GetConformer()
    indices = [a.GetIdx() for a in mol.GetAtoms()]

    for idx in indices:
        atom = mol.GetAtomWithIdx(idx)
        sym = atom.GetSymbol()
        hyb = str(atom.GetHybridization())
        deg = atom.GetDegree()
        
        p0_pt = conf.GetAtomPosition(idx)
        p0 = (p0_pt.x, p0_pt.y, p0_pt.z)
        
        vecs = []
        for n in atom.GetNeighbors():
            pn = conf.GetAtomPosition(n.GetIdx())
            vecs.append(v_norm(v_sub((pn.x, pn.y, pn.z), p0)))
        
        lps = []

        # --- N sp3 (Ammonia) ---
        if sym == 'N' and deg == 3:
            sum_v = (0,0,0)
            for v in vecs: sum_v = v_add(sum_v, v)
            lps.append(v_add(p0, v_scale(v_norm(v_scale(sum_v, -1)), 1.0)))
        
        # --- O sp3 (Water) ---
        elif sym == 'O' and deg == 2 and "SP3" in hyb:
            v1, v2 = vecs[0], vecs[1]
            bisector = v_scale(v_norm(v_add(v1, v2)), -1)
            normal = v_norm(v_cross(v1, v2))
            sin_a, cos_a = 0.816, 0.577
            lps.append(v_add(p0, v_scale(v_add(v_scale(bisector, cos_a), v_scale(normal, sin_a)), 0.7)))
            lps.append(v_add(p0, v_scale(v_add(v_scale(bisector, cos_a), v_scale(normal, -sin_a)), 0.7)))

        # --- O sp2 (Formaldehyde) ---
        elif sym == 'O' and deg == 1 and "SP2" in hyb:
            v_bond = vecs[0]
            neighbor = atom.GetNeighbors()[0]
            n_neighbors = neighbor.GetNeighbors()
            ref_vec = None
            for nn in n_neighbors:
                if nn.GetIdx() != idx:
                    pn_grand = conf.GetAtomPosition(nn.GetIdx())
                    p_neighbor = conf.GetAtomPosition(neighbor.GetIdx())
                    v_grand = v_sub((pn_grand.x, pn_grand.y, pn_grand.z), (p_neighbor.x, p_neighbor.y, p_neighbor.z))
                    ref_vec = v_norm(v_grand)
                    break
            if ref_vec:
                normal = v_norm(v_cross(v_bond, ref_vec))
                in_plane_perp = v_norm(v_cross(normal, v_bond))
                dir1 = v_add(v_scale(v_bond, -0.5), v_scale(in_plane_perp, 0.866))
                dir2 = v_add(v_scale(v_bond, -0.5), v_scale(in_plane_perp, -0.866))
                lps.append(v_add(p0, v_scale(dir1, 0.7)))
                lps.append(v_add(p0, v_scale(dir2, 0.7)))

        # --- S sp3d (Seesaw: SF4) ---
        # 4 bonds. 1 LP. Sum of vectors points roughly away from LP.
        elif sym == 'S' and deg == 4:
            sum_v = (0,0,0)
            for v in vecs: sum_v = v_add(sum_v, v)
            # The sum points towards the bonds, LP is opposite
            lps.append(v_add(p0, v_scale(v_norm(v_scale(sum_v, -1)), 1.2)))

        # --- S sp2 (Bent: SO2) ---
        # 2 bonds, 1 LP. The LP completes the trigonal plane, so it sits opposite
        # the bisector of the two bonds.
        elif sym == 'S' and deg == 2:
            sum_v = (0,0,0)
            for v in vecs: sum_v = v_add(sum_v, v)
            lps.append(v_add(p0, v_scale(v_norm(v_scale(sum_v, -1)), 1.1)))

        # --- Br sp3d2 (Square Pyramidal: BrF5) ---
        # 5 bonds, 1 LP. The four basal bonds cancel, so the sum points at the
        # apex and the LP occupies the vacant sixth octahedral site opposite it.
        elif sym == 'Br' and deg == 5:
            sum_v = (0,0,0)
            for v in vecs: sum_v = v_add(sum_v, v)
            lps.append(v_add(p0, v_scale(v_norm(v_scale(sum_v, -1)), 1.2)))

        # --- Cl sp3d (T-Shaped: ClF3) ---
        # 3 bonds. 2 LPs.
        # Structure is T. Bonds at 0, 90, -90 (approx). 
        # Sum points along the "stem" of the T.
        # LPs are in the plane perpendicular to the T-stem? No, TBP arrangement.
        # If Manual Coords are used: Axial (Z), Eq (X). LPs are at 120 and 240 in Eq plane.
        elif sym == 'Cl' and deg == 3:
            # Assume T-shaped planar-ish structure.
            # Find the "odd one out" (Equatorial bond) vs the nearly linear ones (Axial).
            # But with manual coords, we know Eq is +X. Axial are +/- Z.
            # LPs should be rotated 120 from Eq.
            # Let's find the sum vector (points along X).
            sum_v = (0,0,0)
            for v in vecs: sum_v = v_add(sum_v, v)
            base = v_norm(sum_v) # Points +X (the axial bonds cancel out)

            # Rotate around the AXIAL axis, not the normal of the T-plane. Using the
            # T-plane normal puts both lone pairs in the same plane as the axial
            # bonds; they belong in the equatorial plane, 120 deg either side of the
            # equatorial bond. The axial axis is the most antiparallel bond pair.
            axis = None
            most_opposed = 0.0
            for i in range(len(vecs)):
                for j in range(i + 1, len(vecs)):
                    d = v_dot(vecs[i], vecs[j])
                    if d < most_opposed:
                        most_opposed = d
                        axis = v_norm(v_sub(vecs[i], vecs[j]))
            if axis is None:
                axis = v_norm(v_cross(vecs[0], vecs[1]))

            sin_120, cos_120 = 0.866, -0.5
            perp = v_norm(v_cross(axis, base))
            
            lp1 = v_add(v_scale(base, cos_120), v_scale(perp, sin_120))
            lp2 = v_add(v_scale(base, cos_120), v_scale(perp, -sin_120))
            
            lps.append(v_add(p0, v_scale(lp1, 1.0)))
            lps.append(v_add(p0, v_scale(lp2, 1.0)))

        # --- Xe sp3d (Linear: XeF2) ---
        # 2 bonds. 3 LPs equatorial.
        elif sym == 'Xe' and deg == 2:
            # Bonds are axial (opposite).
            # Find plane perpendicular to bonds.
            v_axis = vecs[0]
            # Create an arbitrary perp vector
            if abs(v_axis[2]) < 0.9: arbitrary = (0,0,1)
            else: arbitrary = (0,1,0)
            perp1 = v_norm(v_cross(v_axis, arbitrary))
            perp2 = v_norm(v_cross(v_axis, perp1))
            
            # 3 LPs at 0, 120, 240 in that plane
            lps.append(v_add(p0, v_scale(perp1, 1.3))) # 0 deg
            
            # Rotate perp1 120 deg
            sin_120, cos_120 = 0.866, -0.5
            lp2_dir = v_add(v_scale(perp1, cos_120), v_scale(perp2, sin_120))
            lp3_dir = v_add(v_scale(perp1, cos_120), v_scale(perp2, -sin_120))
            
            lps.append(v_add(p0, v_scale(lp2_dir, 1.3)))
            lps.append(v_add(p0, v_scale(lp3_dir, 1.3)))

        # --- Xe sp3d2 (Square Planar: XeF4) ---
        elif sym == 'Xe' and deg == 4:
            # Square Planar: LPs are axial (perpendicular to plane).
            cross = v_cross(vecs[0], vecs[1])
            if v_len(cross) > 0.1: normal = v_norm(cross)
            else: normal = v_norm(v_cross(vecs[0], vecs[2]))
            
            if normal:
                lps.append(v_add(p0, v_scale(normal, 1.2)))
                lps.append(v_add(p0, v_scale(normal, -1.2)))

        # Add Dummy Atoms
        for pos in lps:
            new_idx = rwmol.AddAtom(Chem.Atom(0)) # 0 = Dummy
            conf.SetAtomPosition(new_idx, Point3D(*pos))
            
    return rwmol

# Pauling electronegativities, for the bond-dipole sum below.
ELECTRONEGATIVITY = {
    'H': 2.20, 'B': 2.04, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
    'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'Xe': 2.60,
}

# Below this the bond dipoles are treated as cancelling. Molecules built by UFF
# rather than from manual coordinates are not perfectly symmetric, so an exact
# zero never occurs; the real gap is wide (~0.05 for cancelling molecules
# against ~0.9+ for polar ones), so the exact cut-off is not delicate.
DIPOLE_CANCEL_THRESHOLD = 0.25

def compute_dipole(mol):
    """Sum one vector per bond, pointing at the more electronegative atom and
    scaled by the electronegativity difference. This is the method taught in
    class, so the arrow it produces is one students can reproduce by hand."""
    conf = mol.GetConformer()
    total = (0.0, 0.0, 0.0)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        en_a = ELECTRONEGATIVITY.get(a.GetSymbol())
        en_b = ELECTRONEGATIVITY.get(b.GetSymbol())
        if en_a is None or en_b is None: continue
        delta = en_b - en_a
        if abs(delta) < 1e-9: continue
        pa, pb = conf.GetAtomPosition(a.GetIdx()), conf.GetAtomPosition(b.GetIdx())
        direction = v_norm(v_sub((pb.x, pb.y, pb.z), (pa.x, pa.y, pa.z)))
        # A negative delta means A is the more electronegative end; flipping the
        # sign points the vector back towards it.
        total = v_add(total, v_scale(direction, delta))

    magnitude = v_len(total)
    if magnitude < DIPOLE_CANCEL_THRESHOLD:
        return {"cancels": True, "magnitude": round(magnitude, 3)}

    # Draw the arrow through the molecule's centre, pointing at the negative end.
    positions = [conf.GetAtomPosition(a.GetIdx()) for a in mol.GetAtoms()]
    centroid = (
        sum(p.x for p in positions) / len(positions),
        sum(p.y for p in positions) / len(positions),
        sum(p.z for p in positions) / len(positions),
    )
    unit = v_norm(total)

    # Span the molecule rather than using a fixed length, and overshoot the far
    # end so the head clears the atoms. A head buried inside the spheres reads as
    # though the arrow points the other way, which inverts the thing being taught.
    projections = [v_dot(v_sub((p.x, p.y, p.z), centroid), unit) for p in positions]
    tail = min(projections) - 0.5
    head = max(projections) + 1.5

    return {
        "cancels": False,
        "magnitude": round(magnitude, 3),
        "start": [round(c, 3) for c in v_add(centroid, v_scale(unit, tail))],
        "end":   [round(c, 3) for c in v_add(centroid, v_scale(unit, head))],
    }

def analyze_properties(mol):
    atoms = list(mol.GetAtoms())
    heavy_atoms = [a for a in atoms if a.GetAtomicNum() > 1]
    
    # Formula is now overridden in get_molecule_data, so this is just placeholder
    # But we still need mass and angles
    
    # Average molecular weight, NOT CalcExactMolWt (monoisotopic). Students add up
    # periodic-table masses, which are averages: CHCl₃ is 119.38 g/mol, not 117.91.
    mw = round(Descriptors.MolWt(mol), 2)

    # Elements actually present, so the legend can be built per molecule instead of
    # being hardcoded to H/C/O (the library also contains F, Cl, S, P, Xe, B, N).
    elements = []
    for a in atoms:
        sym = a.GetSymbol()
        if a.GetAtomicNum() == 0 or sym in elements: continue
        elements.append(sym)

    # Angles (First 4 found)
    conf = mol.GetConformer()
    angles = []
    import itertools
    for a in heavy_atoms:
        neigh = [n.GetIdx() for n in a.GetNeighbors()]
        if len(neigh) >= 2:
            for i, k in itertools.combinations(neigh, 2):
                p_i = conf.GetAtomPosition(i); p_j = conf.GetAtomPosition(a.GetIdx()); p_k = conf.GetAtomPosition(k)
                ang = angle_deg((p_i.x,p_i.y,p_i.z), (p_j.x,p_j.y,p_j.z), (p_k.x,p_k.y,p_k.z))
                angles.append(f"{mol.GetAtomWithIdx(i).GetSymbol()}-{a.GetSymbol()}-{mol.GetAtomWithIdx(k).GetSymbol()}: {round(ang,1)}°")
    
    # Polarity is overridden too
    
    return {"formula": "", "mw": mw, "polarity": "", "angles": angles[:4], "elements": elements}

# ============================================================
# 5. PREMIUM HTML TEMPLATE
# ============================================================

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Chem 101: Ultimate Visualizer</title>
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;800&display=swap" rel="stylesheet">
  <!-- Served locally, not from a CDN: a blocked or down 3Dmol host would otherwise
       leave a blank viewer with no error, which is fatal mid-lecture. -->
  <script src="vendor/3Dmol-min.js"></script>

  <style>
    :root {
      --primary: #6366f1;
      --accent: #ec4899;
      --bg-grad: linear-gradient(135deg, #a5b4fc 0%, #c084fc 100%);
      --glass: rgba(255, 255, 255, 0.7);
      --glass-border: rgba(255, 255, 255, 0.5);
      --text-main: #1f2937;
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: 'Inter', system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
      background: var(--bg-grad);
      min-height: 100vh;
      display: flex; align-items: center; justify-content: center;
      padding: 20px; color: var(--text-main);
    }
    .dashboard {
      display: grid; grid-template-columns: 350px 1fr; gap: 20px;
      width: 100%; max-width: 1100px;
      /* Grow with the window instead of a fixed 650px: the controls and the
         angle note otherwise push the fact box out of the sidebar. */
      height: min(760px, calc(100vh - 40px)); min-height: 560px;
      background: var(--glass); backdrop-filter: blur(20px); -webkit-backdrop-filter: blur(20px);
      border: 1px solid var(--glass-border); border-radius: 24px;
      box-shadow: 0 25px 50px -12px rgba(0, 0, 0, 0.25); overflow: hidden;
    }
    .sidebar {
      padding: 24px 26px; border-right: 1px solid rgba(255,255,255,0.3);
      display: flex; flex-direction: column; gap: 16px; overflow-y: auto;
    }
    h1 { font-weight: 800; font-size: 1.8rem; letter-spacing: -0.02em; line-height: 1.1; }
    h1 span { background: -webkit-linear-gradient(45deg, var(--primary), var(--accent)); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
    .subtitle { font-size: 0.9rem; color: #6b7280; font-weight: 500; margin-top: 5px; }
    label { font-size: 0.85rem; font-weight: 600; color: #4b5563; text-transform: uppercase; letter-spacing: 0.05em; }
    select {
      width: 100%; padding: 12px; border-radius: 12px; border: 1px solid #d1d5db;
      background: rgba(255,255,255,0.5); font-family: inherit; font-size: 1rem; cursor: pointer; transition: all 0.2s;
    }
    select:hover { border-color: var(--primary); background: #fff; }
    .toggles { display: flex; gap: 10px; }
    .btn {
      flex: 1; padding: 10px; border: none; border-radius: 10px; background: #fff;
      font-weight: 600; color: var(--primary); cursor: pointer;
      box-shadow: 0 2px 5px rgba(0,0,0,0.05); transition: all 0.2s;
    }
    .btn:hover { transform: translateY(-2px); box-shadow: 0 4px 12px rgba(99, 102, 241, 0.2); }
    .btn.active { background: var(--primary); color: #fff; }
    .info-card {
      background: rgba(255,255,255,0.4); border-radius: 16px; padding: 15px; border: 1px solid rgba(255,255,255,0.6);
    }
    .info-row { display: flex; justify-content: space-between; margin-bottom: 8px; font-size: 0.9rem; }
    .info-label { color: #6b7280; }
    .info-val { font-weight: 600; }
    .badge {
      display: inline-block; padding: 4px 8px; border-radius: 6px; font-size: 0.75rem; font-weight: 700; text-transform: uppercase;
    }
    .badge-polar { background: #dbeafe; color: #1e40af; }
    .badge-nonpolar { background: #f3f4f6; color: #374151; }
    .fact-box {
      margin-top: auto; background: linear-gradient(135deg, rgba(255,255,255,0.8), rgba(255,255,255,0.4));
      border-left: 4px solid var(--accent); padding: 12px; border-radius: 0 10px 10px 0;
      font-size: 0.85rem; font-style: italic; line-height: 1.4;
    }
    .viewer-area {
      position: relative; background: #fff; margin: 10px; border-radius: 20px;
      overflow: hidden; box-shadow: inset 0 0 40px rgba(0,0,0,0.02);
    }
    #viewer { position: relative; width: 100%; height: 100%; outline: none; }
    .legend {
      position: absolute; bottom: 20px; right: 20px; background: rgba(255,255,255,0.9);
      padding: 10px 15px; border-radius: 12px; font-size: 0.8rem;
      box-shadow: 0 10px 20px rgba(0,0,0,0.1); display: flex; flex-direction: column; gap: 5px;
    }
    .dot {
      display: inline-block; width: 10px; height: 10px; border-radius: 50%;
      margin-right: 6px; border: 1px solid rgba(0,0,0,0.25); vertical-align: middle;
    }
    .angle-note { margin-top: 8px; font-size: 0.7rem; color: #6b7280; line-height: 1.35; font-style: italic; }
    .dipole-note { font-size: 0.72rem; color: #6d28d9; line-height: 1.35; }
    .dipole-note:not(:empty) { margin-top: 6px; }
    .load-error {
      position: absolute; inset: 0; display: flex; flex-direction: column; gap: 8px;
      align-items: center; justify-content: center; text-align: center; padding: 30px;
      background: #fff; color: #7f1d1d; font-size: 0.9rem; line-height: 1.5;
    }
    /* Author styles beat the UA [hidden] rule, so display:none must be restated. */
    .load-error[hidden] { display: none; }
    @media (max-width: 800px) {
      .dashboard { grid-template-columns: 1fr; height: auto; }
      .viewer-area { height: 400px; }
    }
  </style>
</head>
<body>

  <div class="dashboard">
    <div class="sidebar">
      <div>
        <h1>Chem<span>101</span></h1>
        <div class="subtitle">Ultimate VSEPR Edition</div>
      </div>
      <div>
        <label>Select Molecule</label>
        <select id="mol-select" style="margin-top: 8px;"></select>
      </div>
      <div class="toggles">
        <button class="btn active" id="btn-ball" onclick="setStyle('ball')">Ball & Stick</button>
        <button class="btn" id="btn-space" onclick="setStyle('space')">Space Fill</button>
      </div>
      <div class="toggles">
        <button class="btn active" id="btn-spin" onclick="toggleSpin()">Pause Rotation</button>
        <button class="btn" id="btn-dipole" onclick="toggleDipole()">Show Dipole</button>
      </div>
      <div class="info-card">
        <div class="info-row"><span class="info-label">Formula</span><span class="info-val" id="val-formula">--</span></div>
        <div class="info-row"><span class="info-label">Mol. Weight</span><span class="info-val" id="val-mw">--</span></div>
        <div class="info-row"><span class="info-label">Polarity</span><span id="val-polarity">--</span></div>
        <div class="dipole-note" id="dipole-note"></div>
      </div>
      <div class="info-card">
        <!-- Ordered as the problem is worked: count domains, classify, then read
             off the two geometries and the hybridisation. -->
        <div class="info-row"><span class="info-label">Electron Domains</span><span class="info-val" id="val-steric">--</span></div>
        <div class="info-row"><span class="info-label">VSEPR Class</span><span class="info-val" id="val-axe">--</span></div>
        <div class="info-row"><span class="info-label">Electron Geo.</span><span class="info-val" id="val-geoel">--</span></div>
        <div class="info-row"><span class="info-label">Mol. Geometry</span><span class="info-val" style="color:var(--primary)" id="val-desc">--</span></div>
        <div class="info-row"><span class="info-label">Hybridization</span><span class="info-val" id="val-hybrid">--</span></div>
        <div style="margin-top:8px; font-size:0.8rem; color:#666;" id="val-angles"></div>
        <div class="angle-note">Idealized VSEPR angles. Measured values differ where lone pairs compress bonds (ClF₃ is 87.5° in practice).</div>
      </div>
      <div class="fact-box" id="val-fact">Select a molecule to learn more.</div>
    </div>

    <div class="viewer-area">
      <div id="viewer"></div>
      <div class="legend" id="legend"></div>
      <div class="load-error" id="load-error" hidden>
        <strong>3D viewer failed to load.</strong>
        <span>vendor/3Dmol-min.js could not be read, so molecules cannot be drawn.
        The molecule data on the left is still correct.</span>
      </div>
    </div>
  </div>

  <script>
    const DATA = __DATA__;
    let viewer = null;
    let currStyle = 'ball';
    let spinning = true;
    let showDipole = false;
    let currIdx = 0;

    // Jmol/CPK colours, matching the colorscheme the viewer renders atoms with.
    const ELEMENTS = {
      H:  ["Hydrogen",   "#ffffff"],
      B:  ["Boron",      "#ffb5b5"],
      C:  ["Carbon",     "#909090"],
      N:  ["Nitrogen",   "#3050f8"],
      O:  ["Oxygen",     "#ff0d0d"],
      F:  ["Fluorine",   "#90e050"],
      P:  ["Phosphorus", "#ff8000"],
      S:  ["Sulfur",     "#ffff30"],
      Cl: ["Chlorine",   "#1ff01f"],
      Br: ["Bromine",    "#a62929"],
      Xe: ["Xenon",      "#429eb0"]
    };

    window.onload = function() { waitForLib(0); };

    // Poll rather than checking once: the library is ~500KB and can still be
    // settling on a slow disk or connection when load fires. Only give up after
    // a few seconds, and even then keep the data panel usable.
    function waitForLib(attempt) {
      if (typeof $3Dmol !== 'undefined') { initApp(); return; }
      if (attempt > 60) { initFallback(); return; }
      setTimeout(function() { waitForLib(attempt + 1); }, 100);
    }

    function populateSelect(onPick) {
      const sel = document.getElementById('mol-select');
      DATA.forEach((m, i) => {
        let opt = document.createElement('option');
        opt.value = i; opt.textContent = m.name; sel.appendChild(opt);
      });
      sel.addEventListener('change', () => onPick(sel.value));
    }

    function initApp() {
      viewer = $3Dmol.createViewer(document.getElementById('viewer'), { backgroundColor: 'white' });
      populateSelect(loadMol);
      loadMol(0);
    }

    function initFallback() {
      document.getElementById('load-error').hidden = false;
      populateSelect(i => updateUI(DATA[i]));
      updateUI(DATA[0]);
    }

    function loadMol(idx) {
      if(!viewer) return;
      currIdx = idx;
      const m = DATA[idx];
      viewer.clear();
      viewer.addModel(m.molblock, "sdf");
      applyStyle();
      viewer.spin(spinning ? 'y' : false, 0.5);
      viewer.setHoverable({}, true, function(atom, viewer, event, container) {
        if(!atom.label) {
          let txt = atom.elem === '*' ? 'Lone Pair' : atom.elem;
          atom.label = viewer.addLabel(txt, {
            position: atom, backgroundColor: 'rgba(0,0,0,0.7)', fontColor:'white', fontSize: 12, borderRadius: 4, borderThickness: 0
          });
        }
      }, function(atom, viewer) {
        if(atom.label) { viewer.removeLabel(atom.label); atom.label = null; }
      });
      updateUI(m);
      renderDipole();
    }

    function toggleDipole() {
      showDipole = !showDipole;
      const btn = document.getElementById('btn-dipole');
      btn.textContent = showDipole ? 'Hide Dipole' : 'Show Dipole';
      btn.className = showDipole ? 'btn active' : 'btn';
      renderDipole();
    }

    function renderDipole() {
      if(!viewer) return;
      const note = document.getElementById('dipole-note');
      // Shapes are separate from styles, so clear them rather than restyling.
      viewer.removeAllShapes();
      const d = DATA[currIdx].analysis.dipole;
      if(!showDipole || !d) {
        note.textContent = '';
      } else if(d.cancels) {
        note.textContent = 'Bond dipoles cancel — no net dipole.';
      } else {
        note.textContent = 'Arrow points to the negative end (bond-dipole sum ' + d.magnitude + ').';
        viewer.addArrow({
          start: {x: d.start[0], y: d.start[1], z: d.start[2]},
          end:   {x: d.end[0],   y: d.end[1],   z: d.end[2]},
          radius: 0.11, radiusRatio: 3.0, mid: 0.78, color: '#7c3aed'
        });
      }
      viewer.render();
    }

    function applyStyle() {
      if(currStyle === 'ball') {
        // Jmol, not 'spectrum': spectrum rainbows bonds by atom index, which means
        // nothing chemically and fights the element colours on the spheres.
        viewer.setStyle({not:{elem:'*'}}, { stick: {radius: 0.15, colorscheme:'Jmol'}, sphere: {scale: 0.28, colorscheme:'Jmol'} });
      } else {
        viewer.setStyle({not:{elem:'*'}}, { sphere: {scale: 0.9, colorscheme:'Jmol'} });
      }
      viewer.addStyle({elem:'*'}, { sphere: { radius: 0.4, color: "#ec4899", alpha: 0.6 } });
      viewer.zoomTo(); viewer.zoom(1.4); viewer.render();
    }

    function setStyle(s) {
      currStyle = s;
      document.getElementById('btn-ball').className = s==='ball' ? 'btn active' : 'btn';
      document.getElementById('btn-space').className = s==='space' ? 'btn active' : 'btn';
      applyStyle();
    }

    function toggleSpin() {
      spinning = !spinning;
      const btn = document.getElementById('btn-spin');
      btn.textContent = spinning ? 'Pause Rotation' : 'Resume Rotation';
      btn.className = spinning ? 'btn active' : 'btn';
      if (viewer) viewer.spin(spinning ? 'y' : false, 0.5);
    }

    function renderLegend(a) {
      const box = document.getElementById('legend');
      box.innerHTML = '';
      (a.elements || []).forEach(sym => {
        const [label, color] = ELEMENTS[sym] || [sym, '#bbbbbb'];
        const row = document.createElement('div');
        row.innerHTML = '<span class="dot" style="background:' + color + '"></span>' + label;
        box.appendChild(row);
      });
      if (a.has_lp) {
        const row = document.createElement('div');
        row.innerHTML = '<span class="dot" style="background:#ec4899"></span>Lone Pair';
        box.appendChild(row);
      }
    }

    function updateUI(m) {
      const a = m.analysis;
      renderLegend(a);
      document.getElementById('val-formula').textContent = a.formula;
      document.getElementById('val-mw').textContent = a.mw + " g/mol";
      const polEl = document.getElementById('val-polarity');
      polEl.innerHTML = `<span class="badge ${a.polarity === 'Polar' ? 'badge-polar' : 'badge-nonpolar'}">${a.polarity}</span>`;
      document.getElementById('val-desc').textContent = a.desc;
      document.getElementById('val-geoel').textContent = a.geo_el;
      document.getElementById('val-steric').textContent = a.steric || '--';
      document.getElementById('val-hybrid').textContent = a.hybrid || '--';
      // Name the atom only when it is not the single obvious central one, so
      // "each C" and "the C=O carbon" are not mistaken for the whole molecule.
      const multiCentre = a.center && a.center.indexOf(' ') !== -1;
      document.getElementById('val-axe').textContent =
        (a.axe || '--') + (multiCentre ? ' (at ' + a.center + ')' : '');
      document.getElementById('val-fact').textContent = a.fact;
      document.getElementById('val-angles').innerHTML = a.angles.length > 0 ? a.angles.join("<br>") : "No central angles detected.";
    }
  </script>
</body>
</html>
"""

# ============================================================
# 6. BUILDER
# ============================================================

def check_dipole_agrees(entry, analysis):
    """The geometry and the stated polarity are independent claims, so compare
    them. A wrong coordinate set usually breaks a cancellation, which makes this
    a cheap check on the geometry too."""
    dip = analysis.get('dipole') or {}
    computed = "Non-Polar" if dip.get('cancels') else "Polar"
    if computed != entry['pol']:
        return (f"    [Polarity mismatch] {entry['name']}: listed {entry['pol']}, "
                f"but the bond dipoles sum to {dip.get('magnitude')} "
                f"(threshold {DIPOLE_CANCEL_THRESHOLD}) -> {computed}")
    return None

def build_final():
    print("Building Ultimate Chem 101 Gallery...")
    final_data = []
    warnings = []

    for entry in MOLECULES:
        print(f"  - Processing: {entry['name']}")
        try:
            data = get_molecule_data(entry)
            if data:
                warn = check_dipole_agrees(entry, data['analysis'])
                if warn:
                    print(warn)
                    warnings.append(warn)
                final_item = {
                    "name": entry['name'],
                    "molblock": data['molblock'],
                    "analysis": data['analysis']
                }
                final_data.append(final_item)
            else:
                print(f"    [Error] Could not generate {entry['name']}")
        except Exception as e:
            print(f"    [Error] Failed on {entry['name']}: {e}")

    json_str = json.dumps(final_data, indent=0).replace("</", "<\\/")
    html = HTML_TEMPLATE.replace("__DATA__", json_str)

    # public/index.html is what Vercel serves, so write it here rather than copying
    # by hand — a forgotten copy step means deploying a stale page.
    out_paths = ["out/chem_gallery_final.html", "public/index.html"]
    for out_path in out_paths:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(html)

    print("-" * 30)
    print(f"Built {len(final_data)}/{len(MOLECULES)} molecules.")
    if warnings:
        print(f"{len(warnings)} polarity mismatch(es) — geometry and stated polarity disagree:")
        for w in warnings: print(w)
    else:
        print("Polarity cross-check: all molecules agree with their geometry.")
    for out_path in out_paths:
        print(f"SUCCESS! ✨ {os.path.abspath(out_path)}")

if __name__ == "__main__":
    build_final()