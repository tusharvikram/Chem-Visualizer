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
from rdkit.Chem import AllChem
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
    }
]

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
    if not mol: return None
    mol = Chem.AddHs(mol)
    
    # 1. Generate Basic 3D
    AllChem.EmbedMolecule(mol)
    
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

    # --- Add Visual Lone Pairs ---
    mol_viz = add_lone_pairs(mol)
    
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
            base = v_norm(sum_v) # Points +X
            
            # Find a normal to the T-plane (Y axis approx)
            # Cross product of two non-collinear bonds
            normal = v_norm(v_cross(vecs[0], vecs[1])) 
            if v_len(normal) < 0.1: normal = v_norm(v_cross(vecs[0], vecs[2]))

            # Rotate base vector by +/- 120 deg around normal
            # Or simply: LP1 is rotation of Base, LP2 is rotation.
            # Actually LPs are 120 deg from the Bond in the equatorial plane.
            # Ideally: cos(120)*Base + sin(120)*Cross(Normal, Base)
            sin_120, cos_120 = 0.866, -0.5
            perp = v_norm(v_cross(normal, base))
            
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

def analyze_properties(mol):
    atoms = list(mol.GetAtoms())
    heavy_atoms = [a for a in atoms if a.GetAtomicNum() > 1]
    
    # Formula is now overridden in get_molecule_data, so this is just placeholder
    # But we still need mass and angles
    
    mw = round(Chem.rdMolDescriptors.CalcExactMolWt(mol), 2)
    
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
    
    return {"formula": "", "mw": mw, "polarity": "", "angles": angles[:4]}

# ============================================================
# 5. PREMIUM HTML TEMPLATE
# ============================================================

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Chem 101: Ultimate Visualizer</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;800&display=swap" rel="stylesheet">
  <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
  
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
      font-family: 'Inter', sans-serif;
      background: var(--bg-grad);
      min-height: 100vh;
      display: flex; align-items: center; justify-content: center;
      padding: 20px; color: var(--text-main);
    }
    .dashboard {
      display: grid; grid-template-columns: 350px 1fr; gap: 20px;
      width: 100%; max-width: 1100px; height: 650px;
      background: var(--glass); backdrop-filter: blur(20px); -webkit-backdrop-filter: blur(20px);
      border: 1px solid var(--glass-border); border-radius: 24px;
      box-shadow: 0 25px 50px -12px rgba(0, 0, 0, 0.25); overflow: hidden;
    }
    .sidebar {
      padding: 30px; border-right: 1px solid rgba(255,255,255,0.3);
      display: flex; flex-direction: column; gap: 20px; overflow-y: auto;
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
    .dot { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 6px; }
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
      <div class="info-card">
        <div class="info-row"><span class="info-label">Formula</span><span class="info-val" id="val-formula">--</span></div>
        <div class="info-row"><span class="info-label">Mol. Weight</span><span class="info-val" id="val-mw">--</span></div>
        <div class="info-row"><span class="info-label">Polarity</span><span id="val-polarity">--</span></div>
      </div>
      <div class="info-card">
        <div class="info-row"><span class="info-label">Mol. Geometry</span><span class="info-val" style="color:var(--primary)" id="val-desc">--</span></div>
        <div class="info-row"><span class="info-label">Electron Geo.</span><span class="info-val" id="val-geoel">--</span></div>
        <div style="margin-top:8px; font-size:0.8rem; color:#666;" id="val-angles"></div>
      </div>
      <div class="fact-box" id="val-fact">Select a molecule to learn more.</div>
    </div>

    <div class="viewer-area">
      <div id="viewer"></div>
      <div class="legend">
        <div><span class="dot" style="background:#ec4899;"></span>Lone Pair</div>
        <div><span class="dot" style="background:#ddd;"></span>Hydrogen</div>
        <div><span class="dot" style="background:#333;"></span>Carbon</div>
        <div><span class="dot" style="background:#f00;"></span>Oxygen</div>
      </div>
    </div>
  </div>

  <script>
    const DATA = __DATA__;
    let viewer = null;
    let currStyle = 'ball';
    
    window.onload = function() {
      viewer = $3Dmol.createViewer(document.getElementById('viewer'), { backgroundColor: 'white' });
      const sel = document.getElementById('mol-select');
      DATA.forEach((m, i) => {
        let opt = document.createElement('option');
        opt.value = i; opt.textContent = m.name; sel.appendChild(opt);
      });
      sel.addEventListener('change', () => loadMol(sel.value));
      loadMol(0);
    };

    function loadMol(idx) {
      if(!viewer) return;
      const m = DATA[idx];
      viewer.clear();
      viewer.addModel(m.molblock, "sdf");
      applyStyle();
      viewer.spin('y', 0.5); 
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
    }

    function applyStyle() {
      if(currStyle === 'ball') {
        viewer.setStyle({not:{elem:'*'}}, { stick: {radius: 0.15, color:'spectrum'}, sphere: {scale: 0.28, colorscheme:'Jmol'} });
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

    function updateUI(m) {
      const a = m.analysis;
      document.getElementById('val-formula').textContent = a.formula;
      document.getElementById('val-mw').textContent = a.mw + " g/mol";
      const polEl = document.getElementById('val-polarity');
      polEl.innerHTML = `<span class="badge ${a.polarity === 'Polar' ? 'badge-polar' : 'badge-nonpolar'}">${a.polarity}</span>`;
      document.getElementById('val-desc').textContent = a.desc;
      document.getElementById('val-geoel').textContent = a.geo_el;
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

def build_final():
    print("Building Ultimate Chem 101 Gallery...")
    final_data = []

    for entry in MOLECULES:
        print(f"  - Processing: {entry['name']}")
        try:
            data = get_molecule_data(entry)
            if data:
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

    out_path = "out/chem_gallery_final.html"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    json_str = json.dumps(final_data, indent=0).replace("</", "<\\/")
    html = HTML_TEMPLATE.replace("__DATA__", json_str)
    
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    
    print("-" * 30)
    print(f"SUCCESS! ✨ Open: {os.path.abspath(out_path)}")

if __name__ == "__main__":
    build_final()