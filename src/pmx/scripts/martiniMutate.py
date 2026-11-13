#!/usr/bin/env python

"""Mutate Martini CG protein: map beads and generate hybrid structure."""

import sys
import copy
import argparse
import numpy as np
from pmx.model import Model
from pmx.atom import Atom
from pmx.library import _one_letter
from pmx.forcefield import Topology


# =============
# Input Options
# =============
def parse_options():

    parser = argparse.ArgumentParser(
        description="""
Mutate Martini coarse-grained protein: map beads and generate hybrid structure.

BB beads are mapped via sequence alignment or structural superposition.
SC beads are mapped based on bead-name matching rules.
A hybrid PDB is generated with dummy beads (D-prefixed) at mutation sites.

""",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-i1",
        metavar="protein1.pdb",
        dest="i1",
        type=str,
        help='Input CG structure 1. Default is "protein1.pdb"',
        default="protein1.pdb",
    )
    parser.add_argument(
        "-i2",
        metavar="protein2.pdb",
        dest="i2",
        type=str,
        help='Input CG structure 2. Default is "protein2.pdb"',
        default="protein2.pdb",
    )
    parser.add_argument(
        "-itp1",
        metavar="topol1.itp",
        dest="itp1",
        type=str,
        help="Input topology file for structure 1",
        default=None,
    )
    parser.add_argument(
        "-itp2",
        metavar="topol2.itp",
        dest="itp2",
        type=str,
        help="Input topology file for structure 2",
        default=None,
    )
    parser.add_argument(
        "-ff",
        metavar="martini.itp",
        dest="ff",
        type=str,
        help="Martini force field ITP (atomtypes and nonbond_params)",
        default=None,
    )
    parser.add_argument(
        "-o",
        metavar="hybrid.pdb",
        dest="o",
        type=str,
        help='Output hybrid PDB structure. Default is "hybrid.pdb"',
        default="hybrid.pdb",
    )
    parser.add_argument(
        "-o1",
        metavar="pairs.dat",
        dest="o1",
        type=str,
        help='Output pairs file. Default is "pairs.dat"',
        default="pairs.dat",
    )
    parser.add_argument(
        "-opdb1",
        metavar="mapped1.pdb",
        dest="opdb1",
        type=str,
        help="Optional: output PDB with mapped atoms from structure 1",
        default=None,
    )
    parser.add_argument(
        "-opdb2",
        metavar="mapped2.pdb",
        dest="opdb2",
        type=str,
        help="Optional: output PDB with mapped atoms from structure 2",
        default=None,
    )
    parser.add_argument(
        "-log",
        metavar="mapping.log",
        dest="log",
        type=str,
        help='Output log file. Default is "mapping.log"',
        default="mapping.log",
    )
    parser.add_argument(
        "-ffitp",
        metavar="ffhybrid.itp",
        dest="ffitp",
        type=str,
        help='Output FF itp with dummy atomtypes. Default is "ffhybrid.itp"',
        default="ffhybrid.itp",
    )
    parser.add_argument(
        "--method",
        dest="method",
        type=str,
        choices=["alignment", "sequence"],
        help="Method for BB mapping: alignment (structural) or sequence (residue order). "
        'Default is "sequence"',
        default="sequence",
    )

    args, unknown = parser.parse_known_args()
    return args


# ============================
# Martini FF Parsing
# ============================
def parse_martini_ff(filename):
    """Parse martini FF ITP and return dict of {atomtype: mass}."""
    atomtype_masses = {}
    in_atomtypes = False
    with open(filename) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("[ atomtypes ]"):
                in_atomtypes = True
                continue
            if stripped.startswith("[") and in_atomtypes:
                break
            if in_atomtypes and stripped and not stripped.startswith(";"):
                entr = stripped.split()
                atomtype_masses[entr[0]] = float(entr[1])
    return atomtype_masses


def assign_masses(topol, atomtype_masses):
    """Assign masses to topology atoms from atomtype dict."""
    for atom in topol.atoms:
        atype = atom.atomtype.strip()
        if atype in atomtype_masses:
            atom.m = atomtype_masses[atype]


def parse_ifdef_flexible_pairs(filename):
    """Parse raw .itp file to find atom-index pairs inside #ifdef FLEXIBLE
    and #ifndef FLEXIBLE blocks.

    Returns two sets of (idx1, idx2) tuples (1-based atom indices):
    - ifdef_pairs: pairs inside #ifdef FLEXIBLE (bonds)
    - ifndef_pairs: pairs inside #ifndef FLEXIBLE (constraints)
    """
    ifdef_pairs = set()
    ifndef_pairs = set()
    in_ifdef = False
    in_ifndef = False
    with open(filename) as f:
        for line in f:
            stripped = line.strip()
            # Track ifdef/ifndef FLEXIBLE state
            if stripped == "#ifdef FLEXIBLE":
                in_ifdef = True
                continue
            if stripped == "#ifndef FLEXIBLE":
                in_ifndef = True
                continue
            if stripped == "#endif":
                if in_ifdef:
                    in_ifdef = False
                if in_ifndef:
                    in_ifndef = False
                continue
            # Parse data lines inside the blocks
            if in_ifdef or in_ifndef:
                if not stripped or stripped.startswith(";"):
                    continue
                # Skip section headers inside blocks
                if stripped.startswith("["):
                    continue
                entr = stripped.split()
                if len(entr) >= 2:
                    try:
                        idx1 = int(entr[0])
                        idx2 = int(entr[1])
                        if in_ifdef:
                            ifdef_pairs.add((idx1, idx2))
                        else:
                            ifndef_pairs.add((idx1, idx2))
                    except ValueError:
                        continue
    return ifdef_pairs, ifndef_pairs


def is_perturbed(atom):
    """Check if atom has a B-state assigned."""
    return atom.atomtypeB is not None and atom.atomtypeB != ""


def _zero_dihedral_params(ftype, param_str):
    """Construct zeroed dihedral params: k=0, but preserve angle and multiplicity.

    For type 3 (Ryckaert-Bellemans): all 6 coefficients set to 0.
    For type 2 (improper): angle preserved, k=0.
    For types 1, 4, 9 (periodic): angle preserved, k=0, multiplicity preserved.
    """
    parts = param_str.split()
    if ftype == 3:
        # Ryckaert-Bellemans: 6 coefficients, all zeroed
        return "0 0 0 0 0 0"
    elif ftype == 2:
        # Improper: [angle, k]
        angle = parts[0]
        return "%s 0" % angle
    else:
        # Periodic proper (types 1, 4, 9): [angle, k, multiplicity]
        angle = parts[0]
        mult = parts[2]
        return "%s 0 %s" % (angle, mult)


def _atom_id_key(atom):
    """Sort key function: return atom.id for sorting atom objects by ID."""
    return atom.id


def build_excl_map(exclusions):
    """Build dict: atom -> set of atoms it excludes.

    Example: if exclusions = [[atom7, atom8, atom9, atom10], [atom8, atom9, atom10]]
    then emap[atom7] = {atom8, atom9, atom10}
         emap[atom8] = {atom7, atom9, atom10}
         etc.
    """
    emap = {}
    for excl in exclusions:
        for a in excl:
            if a not in emap:
                emap[a] = set()
            emap[a].update(excl)
            emap[a].discard(a)
    return emap


# ====================
# BB Mapping Functions
# ====================
def get_bb_atoms(model):
    """Extract BB beads from a Model, grouped by residue."""
    bb_atoms = []
    for chain in model.chains:
        for residue in chain.residues:
            for atom in residue.atoms:
                if atom.name.strip() == "BB":
                    bb_atoms.append(atom)
                    break
    return bb_atoms


def get_residue_beads(residue):
    """Get all beads in a residue as a dict keyed by bead name."""
    beads = {}
    for atom in residue.atoms:
        beads[atom.name.strip()] = atom
    return beads


def map_bb_by_sequence(model1, model2):
    """Map BB beads by residue sequence order (1-to-1 positional)."""
    bb1 = get_bb_atoms(model1)
    bb2 = get_bb_atoms(model2)

    if len(bb1) != len(bb2):
        sys.exit(
            f"Error: different number of residues: {len(bb1)} vs {len(bb2)}. "
            "Sequence-based mapping requires equal residue counts."
        )

    pairs = []
    for a1, a2 in zip(bb1, bb2):
        pairs.append((a1, a2))
    return pairs


def kabsch_rotation(P, Q):
    """Compute optimal rotation matrix (Kabsch algorithm) to align P onto Q."""
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, d])
    R = Vt.T @ sign_matrix @ U.T
    return R


def map_bb_by_alignment(model1, model2, cutoff=0.5):
    """Map BB beads by structural superposition (Kabsch) + nearest-neighbor.

    Parameters
    ----------
    cutoff : float
        Distance cutoff (nm) for accepting a BB-BB pair.
    """
    bb1 = get_bb_atoms(model1)
    bb2 = get_bb_atoms(model2)

    coords1 = np.array([[a.x[0], a.x[1], a.x[2]] for a in bb1])
    coords2 = np.array([[a.x[0], a.x[1], a.x[2]] for a in bb2])

    # Center both sets
    centroid1 = coords1.mean(axis=0)
    centroid2 = coords2.mean(axis=0)
    P = coords1 - centroid1
    Q = coords2 - centroid2

    # Kabsch rotation
    R = kabsch_rotation(P, Q)

    # Rotate coords2 onto coords1 frame
    coords2_aligned = (Q @ R.T) + centroid1

    # Nearest-neighbor matching
    pairs = []
    for i, c1 in enumerate(coords1):
        dists = np.linalg.norm(coords2_aligned - c1, axis=1)
        j = np.argmin(dists)
        # Convert from Angstrom to nm for cutoff comparison
        if dists[j] / 10.0 < cutoff:
            pairs.append((bb1[i], bb2[j]))

    return pairs


# ====================
# SC Mapping Functions
# ====================
def map_sidechains(bb_pairs, model1, model2):
    """Map sidechain beads based on BB mapping.

    Rules (initial):
    - If both residues have SC1, map SC1<->SC1
    - Higher SC beads: rules TBD (placeholder)
    """
    sc_pairs = []
    for a1_bb, a2_bb in bb_pairs:
        res1 = a1_bb.molecule
        res2 = a2_bb.molecule
        beads1 = get_residue_beads(res1)
        beads2 = get_residue_beads(res2)

        # Map SC1 if both have it
        if "SC1" in beads1 and "SC1" in beads2:
            sc_pairs.append((beads1["SC1"], beads2["SC1"]))

        # Higher SC beads — placeholder for user-defined rules
        for sc_name in ["SC2", "SC3", "SC4"]:
            if sc_name in beads1 and sc_name in beads2:
                sc_pairs.append((beads1[sc_name], beads2[sc_name]))

    return sc_pairs


# ======================
# Hybrid PDB Generation
# ======================
def write_hybrid_pdb(bb_pairs, outfile):
    """Generate hybrid PDB with dummy beads at mutation sites.

    For non-mutated residues: write stateA beads as-is.
    For mutation sites:
      - Mapped beads: write with stateA coords/names
      - StateA-only beads: write as-is (become dummies in stateB)
      - StateB-only beads: write with stateB coords, D-prefixed names (dummies in stateA)
    """
    atom_id = 0
    mutations = []
    lines = []

    for a1_bb, a2_bb in bb_pairs:
        res1 = a1_bb.molecule
        res2 = a2_bb.molecule
        beads1 = get_residue_beads(res1)
        beads2 = get_residue_beads(res2)
        resname1 = res1.resname.strip()
        resname2 = res2.resname.strip()
        chain_id = a1_bb.chain_id if hasattr(a1_bb, "chain_id") and a1_bb.chain_id else "A"
        resid = a1_bb.resnr

        is_mutation = (resname1 != resname2)
        if is_mutation:
            mutations.append(f"{resname1}{resid}{resname2}")
            letter1 = _one_letter.get(resname1, "X")
            letter2 = _one_letter.get(resname2, "X")
            hybrid_resname = f"{letter1}2{letter2}"
        else:
            hybrid_resname = resname1

        if not is_mutation:
            for atom in res1.atoms:
                atom_id += 1
                lines.append(
                    _pdb_line(atom_id, atom.name.strip(), resname1, chain_id, resid, atom.x)
                )
        else:
            # Write all stateA beads (mapped + stateA-only)
            for atom in res1.atoms:
                atom_id += 1
                lines.append(
                    _pdb_line(atom_id, atom.name.strip(), hybrid_resname, chain_id, resid, atom.x)
                )
            # Write all stateB SC beads with D-prefix — all SC go via dummies
            stateb_sc = [name for name in sorted(beads2.keys()) if name != "BB"]
            for bname in stateb_sc:
                atom_id += 1
                dname = "D" + bname
                lines.append(
                    _pdb_line(atom_id, dname, hybrid_resname, chain_id, resid, beads2[bname].x)
                )

    # Write PDB
    with open(outfile, "w") as f:
        for line in lines:
            f.write(line + "\n")
        f.write(f"TER\nEND\n")

    return mutations


def _pdb_line(atom_id, atom_name, resname, chain_id, resid, coords):
    """Format a single PDB ATOM line."""
    # PDB format: cols 1-6 record, 7-11 serial, 13-16 name, 17 altloc,
    # 18-20 resname, 22 chain, 23-26 resseq, 31-38 x, 39-46 y, 47-54 z
    if len(atom_name) <= 3:
        name_field = f" {atom_name:<3s}"
    else:
        name_field = f"{atom_name:<4s}"
    return (
        f"ATOM  {atom_id:>5d} {name_field} {resname:<3s} {chain_id}{resid:>4d}    "
        f"{coords[0]:>8.3f}{coords[1]:>8.3f}{coords[2]:>8.3f}  1.00  0.00"
    )


# ============================
# Hybrid Topology Construction
# ============================
def generate_excl_pairs(topol):
    """Convert exclusion groups into a list of (atom_object, atom_object) pairs.

    Each exclusion group [a, b, c, d] produces pairs: (a,b), (a,c), (a,d).
    The first atom in the group excludes all others.
    """
    excl_pairs = []
    for excl in topol.exclusions:
        if len(excl) < 2:
            continue
        # First atom in the group excludes all the rest
        ref_atom = excl[0]
        for other in excl[1:]:
            excl_pairs.append((ref_atom, other))
    return excl_pairs


def build_hybrid_topology(topol1, topol2, all_pairs):
    """Build hybrid topology from topol1 (stateA) and topol2 (stateB).

    Creates a deep copy of topol1 as the hybrid base.
    topol1 and topol2 remain untouched.
    Returns the hybrid topology object.
    """
    # Create hybrid topology as a deep copy of topol1
    hybrid = copy.deepcopy(topol1)

    # Build mapping dicts from the atom pairs.
    # pair_map_1to2[original_topol1_id] = original_topol2_id
    # pair_map_2to1[original_topol2_id] = original_topol1_id
    # mapped_ids_1: set of topol1 atom IDs that have a partner in topol2
    # mapped_ids_2: set of topol2 atom IDs that have a partner in topol1
    pair_map_1to2 = {}
    pair_map_2to1 = {}
    mapped_ids_1 = set()
    mapped_ids_2 = set()
    for a1, a2 in all_pairs:
        pair_map_1to2[a1.id] = a2.id
        pair_map_2to1[a2.id] = a1.id
        mapped_ids_1.add(a1.id)
        mapped_ids_2.add(a2.id)

    # Build lookup: topol2_by_id[original_topol2_id] = topol2 atom object
    topol2_by_id = {}
    for atom in topol2.atoms:
        topol2_by_id[atom.id] = atom

    # Build lookup: topol1_by_id[original_topol1_id] = topol1 atom object
    topol1_by_id = {}
    for atom in topol1.atoms:
        topol1_by_id[atom.id] = atom

    # --- Assign B-states to hybrid atoms ---
    for atom in hybrid.atoms:
        if atom.id in mapped_ids_1:
            # This atom has a partner in topol2 — check if type/charge differs
            partner_id = pair_map_1to2[atom.id]
            partner = topol2_by_id[partner_id]
            if atom.atomtype != partner.atomtype or atom.q != partner.q:
                # Different type or charge: assign B-state from topol2 partner
                atom.atomtypeB = partner.atomtype
                atom.qB = partner.q
                atom.mB = partner.m
            # else: identical in both states, no B-state needed
        else:
            # This atom has no partner in topol2: it becomes a dummy in stateB
            atom.atomtypeB = "DUM"
            atom.qB = 0.0
            atom.mB = atom.m

    # --- Insert stateB-only atoms as DUM in stateA ---
    for atom_id in sorted(topol2_by_id.keys()):
        if atom_id in mapped_ids_2:
            # Already mapped, skip
            continue
        # This topol2 atom has no partner in topol1: create a DUM-in-A atom
        atom2 = topol2_by_id[atom_id]
        new_atom = copy.deepcopy(atom2)
        # In stateA this atom is a dummy
        new_atom.atomtype = "DUM"
        new_atom.q = 0.0
        # In stateB it has its real topol2 parameters
        new_atom.atomtypeB = atom2.atomtype
        new_atom.qB = atom2.q
        new_atom.mB = atom2.m
        # Insert after the last atom of the same residue in hybrid
        insert_idx = None
        for i, a in enumerate(hybrid.atoms):
            if a.resnr == atom2.resnr:
                insert_idx = i
        if insert_idx is not None:
            hybrid.atoms.insert(insert_idx + 1, new_atom)
        else:
            hybrid.atoms.append(new_atom)

    # --- Build atom object mappings before renumbering ---
    # hybrid_to_topol1[hybrid_atom] = topol1_atom (None if DUM-in-A)
    # hybrid_to_topol2[hybrid_atom] = topol2_atom (None if DUM-in-B)
    # topol1_to_hybrid[topol1_atom] = hybrid_atom
    # topol2_to_hybrid[topol2_atom] = hybrid_atom
    hybrid_to_topol1 = {}
    hybrid_to_topol2 = {}
    topol1_to_hybrid = {}
    topol2_to_hybrid = {}

    for atom in hybrid.atoms:
        if atom.atomtype == "DUM":
            # DUM-in-A: came from topol2, no topol1 counterpart
            hybrid_to_topol1[atom] = None
            t2_atom = topol2_by_id[atom.id]
            hybrid_to_topol2[atom] = t2_atom
            topol2_to_hybrid[t2_atom] = atom
        elif atom.id in mapped_ids_1:
            # Mapped atom: has counterparts in both topologies
            t1_atom = topol1_by_id[atom.id]
            t2_atom = topol2_by_id[pair_map_1to2[atom.id]]
            hybrid_to_topol1[atom] = t1_atom
            hybrid_to_topol2[atom] = t2_atom
            topol1_to_hybrid[t1_atom] = atom
            topol2_to_hybrid[t2_atom] = atom
        else:
            # DUM-in-B: came from topol1, no topol2 counterpart
            t1_atom = topol1_by_id[atom.id]
            hybrid_to_topol1[atom] = t1_atom
            hybrid_to_topol2[atom] = None
            topol1_to_hybrid[t1_atom] = atom

    # Renumber hybrid atoms sequentially
    for i, atom in enumerate(hybrid.atoms):
        atom.id = i + 1
        atom.cgnr = i + 1

    # --- Bond B-states ---
    # Index topol2 bonds by atom ID pair (one ordering only)
    # topol2_bonds[(topol2_id1, topol2_id2)] = bond entry from topol2
    topol2_bonds = {}
    for b in topol2.bonds:
        topol2_bonds[(b[0].id, b[1].id)] = b

    # Pass 1: hybrid bonds involving perturbed atoms — assign B-state from topol2
    for i, b in enumerate(hybrid.bonds):
        if len(b) != 4:
            continue
        if not (is_perturbed(b[0]) or is_perturbed(b[1])):
            continue
        ftype = b[2]
        lA = b[3][0]
        kA = b[3][1]
        # Get topol2 atom objects for both atoms in this bond
        t2_atom1 = hybrid_to_topol2[b[0]]
        t2_atom2 = hybrid_to_topol2[b[1]]
        if t2_atom1 is not None and t2_atom2 is not None:
            # Both atoms have topol2 counterparts: look up bond in topol2
            bB = topol2_bonds.get((t2_atom1.id, t2_atom2.id))
            if bB is None:
                bB = topol2_bonds.get((t2_atom2.id, t2_atom1.id))
            if bB:
                # Found partner bond: B-state params from topol2
                lB = bB[3][0]
                kB = bB[3][1]
                hybrid.bonds[i] = [b[0], b[1], ftype, [ftype, lA, kA], [ftype, lB, kB]]
                continue
        # Atom becomes DUM in B: keep A params for both states
        hybrid.bonds[i] = [b[0], b[1], ftype, [ftype, lA, kA], [ftype, lA, kA]]

    # Pass 2: topol2 bonds where either atom is DUM in stateA — add to hybrid
    for b in topol2.bonds:
        # b[0], b[1] are topol2 atom objects — look up corresponding hybrid atoms
        a1_h = topol2_to_hybrid.get(b[0])
        a2_h = topol2_to_hybrid.get(b[1])
        if a1_h is None or a2_h is None:
            continue
        if not (a1_h.atomtype.startswith("DUM") or a2_h.atomtype.startswith("DUM")):
            continue
        # At least one atom is DUM in stateA: this bond only exists in stateB
        ftype = b[2]
        lB = b[3][0]
        kB = b[3][1]
        # B params for both states (keeps force active on dummy atoms)
        hybrid.bonds.append([a1_h, a2_h, ftype, [ftype, lB, kB], [ftype, lB, kB]])

    # --- Angle B-states ---
    # Index topol2 angles by atom ID triple (one ordering only)
    # topol2_angles[(topol2_id1, topol2_id2, topol2_id3)] = angle entry from topol2
    topol2_angles = {}
    for ang in topol2.angles:
        topol2_angles[(ang[0].id, ang[1].id, ang[2].id)] = ang

    # Pass 1: hybrid angles involving perturbed atoms — assign B-state from topol2
    for i, ang in enumerate(hybrid.angles):
        if len(ang) != 5:
            continue
        if not (is_perturbed(ang[0]) or is_perturbed(ang[1]) or is_perturbed(ang[2])):
            continue
        ftype = ang[3]
        thetaA = ang[4][0]
        kA = ang[4][1]
        # Get topol2 atom objects for all three atoms
        t2_atom1 = hybrid_to_topol2[ang[0]]
        t2_atom2 = hybrid_to_topol2[ang[1]]
        t2_atom3 = hybrid_to_topol2[ang[2]]
        if t2_atom1 is not None and t2_atom2 is not None and t2_atom3 is not None:
            # All atoms have topol2 counterparts: look up angle in topol2
            angB = topol2_angles.get((t2_atom1.id, t2_atom2.id, t2_atom3.id))
            if angB is None:
                angB = topol2_angles.get((t2_atom3.id, t2_atom2.id, t2_atom1.id))
            if angB:
                # Found partner angle: B-state params from topol2
                thetaB = angB[4][0]
                kB = angB[4][1]
                hybrid.angles[i] = [ang[0], ang[1], ang[2], ftype,
                                    [ftype, thetaA, kA], [ftype, thetaB, kB]]
                continue
        # Atom becomes DUM in B: keep A params for both states
        hybrid.angles[i] = [ang[0], ang[1], ang[2], ftype,
                            [ftype, thetaA, kA], [ftype, thetaA, kA]]

    # Pass 2: topol2 angles where any atom is DUM in stateA — add to hybrid
    for ang in topol2.angles:
        # ang[0], ang[1], ang[2] are topol2 atom objects
        a1_h = topol2_to_hybrid.get(ang[0])
        a2_h = topol2_to_hybrid.get(ang[1])
        a3_h = topol2_to_hybrid.get(ang[2])
        if a1_h is None or a2_h is None or a3_h is None:
            continue
        if not (a1_h.atomtype.startswith("DUM") or a2_h.atomtype.startswith("DUM")
                or a3_h.atomtype.startswith("DUM")):
            continue
        # At least one atom is DUM in stateA: this angle only exists in stateB
        ftype = ang[3]
        thetaB = ang[4][0]
        kB = ang[4][1]
        # B params for both states (keeps force active on dummy atoms)
        hybrid.angles.append([a1_h, a2_h, a3_h, ftype,
                              [ftype, thetaB, kB], [ftype, thetaB, kB]])

    # --- Constraint B-states ---
    # Index topol2 constraints by atom ID pair (one ordering only)
    # topol2_constraints[(topol2_id1, topol2_id2)] = constraint entry from topol2
    topol2_constraints = {}
    for c in topol2.constraints:
        topol2_constraints[(c[0].id, c[1].id)] = c

    # Pass 1: hybrid constraints involving perturbed atoms — assign B-state from topol2
    for i, c in enumerate(hybrid.constraints):
        if len(c) != 4:
            continue
        if not (is_perturbed(c[0]) or is_perturbed(c[1])):
            continue
        distA = c[3]
        # Get topol2 atom objects for both atoms
        t2_atom1 = hybrid_to_topol2[c[0]]
        t2_atom2 = hybrid_to_topol2[c[1]]
        if t2_atom1 is not None and t2_atom2 is not None:
            # Both atoms have topol2 counterparts: look up constraint in topol2
            cB = topol2_constraints.get((t2_atom1.id, t2_atom2.id))
            if cB is None:
                cB = topol2_constraints.get((t2_atom2.id, t2_atom1.id))
            if cB:
                # Found partner constraint: B-state distance from topol2
                distB = cB[3]
                hybrid.constraints[i] = [c[0], c[1], c[2], distA, distB]
                continue
        # Atom becomes DUM in B: keep A distance for both states
        hybrid.constraints[i] = [c[0], c[1], c[2], distA, distA]

    # Pass 2: topol2 constraints where either atom is DUM in stateA — add to hybrid
    for c in topol2.constraints:
        # c[0], c[1] are topol2 atom objects
        a1_h = topol2_to_hybrid.get(c[0])
        a2_h = topol2_to_hybrid.get(c[1])
        if a1_h is None or a2_h is None:
            continue
        if not (a1_h.atomtype.startswith("DUM") or a2_h.atomtype.startswith("DUM")):
            continue
        # At least one atom is DUM in stateA: this constraint only exists in stateB
        distB = c[3]
        # B distance for both states (keeps constraint active on dummy atoms)
        hybrid.constraints.append([a1_h, a2_h, c[2], distB, distB])

    # --- Dihedral B-states ---
    # Index topol2 dihedrals by atom ID quad (one ordering only)
    # Multiple dihedrals can exist for same quad (different multiplicities), so store as list
    # topol2_dihedrals[(topol2_id1, topol2_id2, topol2_id3, topol2_id4)] = [dih1, dih2, ...]
    topol2_dihedrals = {}
    for d in topol2.dihedrals:
        key = (d[0].id, d[1].id, d[2].id, d[3].id)
        if key not in topol2_dihedrals:
            topol2_dihedrals[key] = []
        topol2_dihedrals[key].append(d)

    # Pass 1: hybrid dihedrals involving perturbed atoms
    # Rebuild entire list to ensure paired fade-out/fade-in entries are adjacent
    new_dihedrals = []
    for d in hybrid.dihedrals:
        if len(d) != 6:
            new_dihedrals.append(d)
            continue
        if not (is_perturbed(d[0]) or is_perturbed(d[1])
                or is_perturbed(d[2]) or is_perturbed(d[3])):
            # Unperturbed dihedral: keep as-is
            new_dihedrals.append(d)
            continue
        ftype = d[4]
        paramA = d[5]
        # Check if any atom becomes DUM in B
        has_dum_in_B = (d[0].atomtypeB == "DUM" or d[1].atomtypeB == "DUM"
                        or d[2].atomtypeB == "DUM" or d[3].atomtypeB == "DUM")
        if has_dum_in_B:
            # Atom disappears in B: keep A params for both states
            new_dihedrals.append([d[0], d[1], d[2], d[3], ftype, paramA, paramA])
        else:
            # All atoms mapped: fade out A params (k=0, preserve angle and multiplicity)
            zeroA = _zero_dihedral_params(ftype, paramA)
            new_dihedrals.append([d[0], d[1], d[2], d[3], ftype, paramA, zeroA])
            # Find partner dihedral in topol2 and fade in B params
            t2_a1 = hybrid_to_topol2[d[0]]
            t2_a2 = hybrid_to_topol2[d[1]]
            t2_a3 = hybrid_to_topol2[d[2]]
            t2_a4 = hybrid_to_topol2[d[3]]
            if t2_a1 is not None and t2_a2 is not None and t2_a3 is not None and t2_a4 is not None:
                key_B = (t2_a1.id, t2_a2.id, t2_a3.id, t2_a4.id)
                key_B_alt = (t2_a4.id, t2_a3.id, t2_a2.id, t2_a1.id)
                dih_list = topol2_dihedrals.get(key_B)
                if dih_list is None:
                    dih_list = topol2_dihedrals.get(key_B_alt)
                if dih_list:
                    # Find matching dihedral by function type
                    for dB in dih_list:
                        if dB[4] == ftype:
                            paramB = dB[5]
                            # Fade in B (zero -> paramB), adjacent to fade-out above
                            zeroB = _zero_dihedral_params(ftype, paramB)
                            new_dihedrals.append([d[0], d[1], d[2], d[3], ftype, zeroB, paramB])
                            break

    # Pass 2: topol2 dihedrals where any atom is DUM in stateA
    for d in topol2.dihedrals:
        # d[0]..d[3] are topol2 atom objects
        a1_h = topol2_to_hybrid.get(d[0])
        a2_h = topol2_to_hybrid.get(d[1])
        a3_h = topol2_to_hybrid.get(d[2])
        a4_h = topol2_to_hybrid.get(d[3])
        if a1_h is None or a2_h is None or a3_h is None or a4_h is None:
            continue
        if not (a1_h.atomtype.startswith("DUM") or a2_h.atomtype.startswith("DUM")
                or a3_h.atomtype.startswith("DUM") or a4_h.atomtype.startswith("DUM")):
            continue
        # At least one atom is DUM in stateA: this dihedral only exists in stateB
        ftype = d[4]
        paramB = d[5]
        # B params for both states (keeps force active on dummy atoms)
        new_dihedrals.append([a1_h, a2_h, a3_h, a4_h, ftype, paramB, paramB])

    hybrid.dihedrals = new_dihedrals

    # --- Exclusions ---
    # Step 1: generate exclusion pairs for stateA (topol1 atom objects)
    exclListA = generate_excl_pairs(topol1)

    # Step 2: generate exclusion pairs for stateB (topol2 atom objects)
    exclListB = generate_excl_pairs(topol2)

    # Build set of topol2 exclusion pairs for fast lookup in step 3
    exclSetB = set()
    for (a1, a2) in exclListB:
        exclSetB.add((a1.id, a2.id))

    # Build set of topol1 exclusion pairs for fast lookup in step 4
    exclSetA = set()
    for (a1, a2) in exclListA:
        exclSetA.add((a1.id, a2.id))

    # Step 3: go over exclListA
    # Collect valid exclusions as pairs of hybrid atom objects
    final_exclusions = []

    for (a1_t1, a2_t1) in exclListA:
        # Find corresponding hybrid atoms
        a1_h = topol1_to_hybrid.get(a1_t1)
        a2_h = topol1_to_hybrid.get(a2_t1)
        if a1_h is None or a2_h is None:
            sys.exit(f"Error: cannot find hybrid atom for exclusion in topol1 "
                     f"(atoms {a1_t1.id}, {a2_t1.id})")
        if a1_h.atomtypeB == "DUM" or a2_h.atomtypeB == "DUM":
            # At least one atom becomes dummy in B: keep this exclusion
            final_exclusions.append((a1_h, a2_h))
        else:
            # Neither is dummy in B: check if exclusion exists in stateB
            t2_a1 = hybrid_to_topol2[a1_h]
            t2_a2 = hybrid_to_topol2[a2_h]
            if (t2_a1.id, t2_a2.id) in exclSetB or (t2_a2.id, t2_a1.id) in exclSetB:
                # Exists in both states: keep it
                final_exclusions.append((a1_h, a2_h))
            else:
                sys.exit(f"Exclusion for atom {a1_h.id} - {a2_h.id} in state A "
                         f"disappears in state B")

    # Step 4: go over exclListB
    for (a1_t2, a2_t2) in exclListB:
        # Find corresponding hybrid atoms
        a1_h = topol2_to_hybrid.get(a1_t2)
        a2_h = topol2_to_hybrid.get(a2_t2)
        if a1_h is None or a2_h is None:
            sys.exit(f"Error: cannot find hybrid atom for exclusion in topol2 "
                     f"(atoms {a1_t2.id}, {a2_t2.id})")
        if a1_h.atomtype == "DUM" or a2_h.atomtype == "DUM":
            # At least one atom is dummy in A: keep this exclusion
            final_exclusions.append((a1_h, a2_h))
        else:
            # Neither is dummy in A: check if exclusion exists in stateA
            t1_a1 = hybrid_to_topol1[a1_h]
            t1_a2 = hybrid_to_topol1[a2_h]
            if (t1_a1.id, t1_a2.id) in exclSetA or (t1_a2.id, t1_a1.id) in exclSetA:
                # Already handled in step 3: skip to avoid duplicates
                pass
            else:
                sys.exit(f"Exclusion for atom {a1_h.id} - {a2_h.id} in state B "
                         f"disappears in state A")

    # Convert final_exclusions to exclusion groups for topology output
    # Group by first atom: all pairs (ref, other) with same ref form one exclusion line
    excl_groups = {}
    for (a1_h, a2_h) in final_exclusions:
        if a1_h not in excl_groups:
            excl_groups[a1_h] = []
        excl_groups[a1_h].append(a2_h)

    # Write to hybrid.exclusions as list of atom object lists
    # Sort reference atoms by their hybrid ID for deterministic output
    sorted_refs = list(excl_groups.keys())
    sorted_refs.sort(key=_atom_id_key)

    hybrid.exclusions = []
    for ref_atom in sorted_refs:
        group = [ref_atom]
        for other_atom in excl_groups[ref_atom]:
            group.append(other_atom)
        hybrid.exclusions.append(group)

    # Update flag so write_exclusions is called during output
    if len(hybrid.exclusions) > 0:
        hybrid.has_exclusions = True
    else:
        hybrid.has_exclusions = False

    return hybrid


# =====================================
# FLEXIBLE bonds/constraints B-states
# =====================================
def build_ifdef_FLEXIBLE_bonds_hybrid(ifdef_bonds_1, ifdef_bonds_2, hybrid, all_pairs):
    """Build B-state parameters for FLEXIBLE bonds.

    Returns list of hybrid flexible bond entries with B-state info.
    Each entry: [hybrid_atom1, hybrid_atom2, ftype, [lA, kA], [lB, kB]]
    """
    # Build mapping from topol1 atom id -> topol2 atom id
    pair_map_1to2 = {}
    for a1, a2 in all_pairs:
        pair_map_1to2[a1.id] = a2.id

    # Build lookup: hybrid atoms by their original topol1 id (before renumbering)
    # We need to find hybrid atoms corresponding to topol1/topol2 atoms.
    # hybrid.atoms have been renumbered, but we can match by resnr + name.
    hybrid_by_resnr_name = {}
    for atom in hybrid.atoms:
        key = (atom.resnr, atom.name.strip())
        hybrid_by_resnr_name[key] = atom

    # Index topol2 flexible bonds by atom pair (resnr, name)
    ifdef_bonds_2_by_pair = {}
    for b in ifdef_bonds_2:
        key = ((b[0].resnr, b[0].name.strip()), (b[1].resnr, b[1].name.strip()))
        ifdef_bonds_2_by_pair[key] = b

    # Build hybrid flexible bonds with B-states
    hybrid_ifdef_FLEXIBLE_bonds = []
    for b in ifdef_bonds_1:
        # Find hybrid atoms by resnr + name
        key1 = (b[0].resnr, b[0].name.strip())
        key2 = (b[1].resnr, b[1].name.strip())
        h_atom1 = hybrid_by_resnr_name.get(key1)
        h_atom2 = hybrid_by_resnr_name.get(key2)
        if h_atom1 is None or h_atom2 is None:
            sys.exit(f"Error: cannot find hybrid atom for flexible bond "
                     f"({key1}, {key2})")
        ftype = b[2]
        lA = b[3][0]
        kA = b[3][1]

        # Check if both atoms are mapped (have partners in topol2)
        if b[0].id in pair_map_1to2 and b[1].id in pair_map_1to2:
            # Both atoms mapped: look for matching flexible bond in topol2
            t2_key1 = (b[0].resnr, b[0].name.strip())
            t2_key2 = (b[1].resnr, b[1].name.strip())
            # For mapped atoms, find partner names in topol2
            # Use the pair mapping to find topol2 residue/name
            bB = ifdef_bonds_2_by_pair.get((t2_key1, t2_key2))
            if bB is None:
                bB = ifdef_bonds_2_by_pair.get((t2_key2, t2_key1))
            if bB is None:
                # Bond in ifdef stateA but not in ifdef stateB: error
                sys.exit(f"Error: flexible bond ({key1}, {key2}) exists in "
                         f"state A but has no matching flexible bond in state B. "
                         f"ifdef blocks must be consistent between topologies.")
            lB = bB[3][0]
            kB = bB[3][1]
            hybrid_ifdef_FLEXIBLE_bonds.append([h_atom1, h_atom2, ftype, [lA, kA], [lB, kB]])
        else:
            # At least one atom becomes DUM: keep A params for both states
            hybrid_ifdef_FLEXIBLE_bonds.append([h_atom1, h_atom2, ftype, [lA, kA], [lA, kA]])

    # Check for topol2-only flexible bonds (DUM in stateA)
    # These are bonds in ifdef_bonds_2 where at least one atom has no partner in topol1
    mapped_ids_2 = set()
    for a1, a2 in all_pairs:
        mapped_ids_2.add(a2.id)
    for b in ifdef_bonds_2:
        if b[0].id not in mapped_ids_2 or b[1].id not in mapped_ids_2:
            # At least one atom is DUM in stateA
            key1 = (b[0].resnr, b[0].name.strip())
            key2 = (b[1].resnr, b[1].name.strip())
            # Find hybrid atoms — DUM-in-A atoms have D-prefixed names
            h_atom1 = hybrid_by_resnr_name.get(key1)
            if h_atom1 is None:
                # Try with D-prefix
                dkey1 = (b[0].resnr, "D" + b[0].name.strip())
                h_atom1 = hybrid_by_resnr_name.get(dkey1)
            h_atom2 = hybrid_by_resnr_name.get(key2)
            if h_atom2 is None:
                dkey2 = (b[1].resnr, "D" + b[1].name.strip())
                h_atom2 = hybrid_by_resnr_name.get(dkey2)
            if h_atom1 is None or h_atom2 is None:
                continue
            ftype = b[2]
            lB = b[3][0]
            kB = b[3][1]
            # B params for both states (keeps force active on dummy atoms)
            hybrid_ifdef_FLEXIBLE_bonds.append([h_atom1, h_atom2, ftype, [lB, kB], [lB, kB]])

    return hybrid_ifdef_FLEXIBLE_bonds


def build_ifndef_FLEXIBLE_constraints_hybrid(ifndef_constr_1, ifndef_constr_2, hybrid, all_pairs):
    """Build B-state parameters for FLEXIBLE constraints.

    Returns list of hybrid flexible constraint entries with B-state info.
    Each entry: [hybrid_atom1, hybrid_atom2, ftype, distA, distB]
    """
    # Build mapping from topol1 atom id -> topol2 atom id
    pair_map_1to2 = {}
    for a1, a2 in all_pairs:
        pair_map_1to2[a1.id] = a2.id

    # Build lookup: hybrid atoms by resnr + name
    hybrid_by_resnr_name = {}
    for atom in hybrid.atoms:
        key = (atom.resnr, atom.name.strip())
        hybrid_by_resnr_name[key] = atom

    # Index topol2 flexible constraints by atom pair (resnr, name)
    ifndef_constr_2_by_pair = {}
    for c in ifndef_constr_2:
        key = ((c[0].resnr, c[0].name.strip()), (c[1].resnr, c[1].name.strip()))
        ifndef_constr_2_by_pair[key] = c

    # Build hybrid flexible constraints with B-states
    hybrid_ifndef_FLEXIBLE_constraints = []
    for c in ifndef_constr_1:
        key1 = (c[0].resnr, c[0].name.strip())
        key2 = (c[1].resnr, c[1].name.strip())
        h_atom1 = hybrid_by_resnr_name.get(key1)
        h_atom2 = hybrid_by_resnr_name.get(key2)
        if h_atom1 is None or h_atom2 is None:
            sys.exit(f"Error: cannot find hybrid atom for flexible constraint "
                     f"({key1}, {key2})")
        ftype = c[2]
        distA = c[3]

        if c[0].id in pair_map_1to2 and c[1].id in pair_map_1to2:
            # Both atoms mapped: look for matching flexible constraint in topol2
            cB = ifndef_constr_2_by_pair.get((key1, key2))
            if cB is None:
                cB = ifndef_constr_2_by_pair.get((key2, key1))
            if cB is None:
                sys.exit(f"Error: flexible constraint ({key1}, {key2}) exists in "
                         f"state A but has no matching flexible constraint in state B. "
                         f"ifdef blocks must be consistent between topologies.")
            distB = cB[3]
            hybrid_ifndef_FLEXIBLE_constraints.append([h_atom1, h_atom2, ftype, distA, distB])
        else:
            # At least one atom becomes DUM: keep A distance for both states
            hybrid_ifndef_FLEXIBLE_constraints.append([h_atom1, h_atom2, ftype, distA, distA])

    # Check for topol2-only flexible constraints (DUM in stateA)
    mapped_ids_2 = set()
    for a1, a2 in all_pairs:
        mapped_ids_2.add(a2.id)
    for c in ifndef_constr_2:
        if c[0].id not in mapped_ids_2 or c[1].id not in mapped_ids_2:
            key1 = (c[0].resnr, c[0].name.strip())
            key2 = (c[1].resnr, c[1].name.strip())
            h_atom1 = hybrid_by_resnr_name.get(key1)
            if h_atom1 is None:
                dkey1 = (c[0].resnr, "D" + c[0].name.strip())
                h_atom1 = hybrid_by_resnr_name.get(dkey1)
            h_atom2 = hybrid_by_resnr_name.get(key2)
            if h_atom2 is None:
                dkey2 = (c[1].resnr, "D" + c[1].name.strip())
                h_atom2 = hybrid_by_resnr_name.get(dkey2)
            if h_atom1 is None or h_atom2 is None:
                continue
            ftype = c[2]
            distB = c[3]
            hybrid_ifndef_FLEXIBLE_constraints.append([h_atom1, h_atom2, ftype, distB, distB])

    return hybrid_ifndef_FLEXIBLE_constraints


def write_ffhybrid_itp(outfile):
    """Write ffhybrid.itp containing the DUM atomtype definition.

    Format matches martini_v3.0.0.itp: name mass charge A sigma epsilon
    """
    with open(outfile, "w") as fp:
        print("; Force field extension for hybrid topology dummy atoms", file=fp)
        print("", file=fp)
        print("[ atomtypes ]", file=fp)
        print("DUM 72.0 0.000 A 0.0 0.0", file=fp)


# ============================
# Custom Hybrid ITP Writer
# ============================
def write_hybrid_itp(hybrid, ifdef_FLEXIBLE_bonds, ifndef_FLEXIBLE_constraints, outfile):
    """Write hybrid topology with #ifdef FLEXIBLE / #ifndef FLEXIBLE structure.

    Uses hybrid object's write methods for atoms, angles, dihedrals, exclusions.
    Handles bonds and constraints manually to insert ifdef blocks.
    """
    with open(outfile, "w") as fp:
        # Write moleculetype and atoms sections
        hybrid.write_moleculetype(fp)
        hybrid.write_atoms(fp, charges="AB", atomtypes="AB",
                           dummy_qA="off", dummy_qB="off")

        # Write bonds section: regular bonds + #ifdef FLEXIBLE block
        hybrid.write_bonds(fp, state="AB")
        if ifdef_FLEXIBLE_bonds:
            print("\n#ifdef FLEXIBLE", file=fp)
            print("; Sidechain bonds (flexible alternative to constraints)", file=fp)
            for b in ifdef_FLEXIBLE_bonds:
                lA = b[3][0]
                kA = b[3][1]
                lB = b[4][0]
                kB = b[4][1]
                print("%6d %6d %6d %14.6f %14.6f %14.6f %14.6f"
                      % (b[0].id, b[1].id, b[2], lA, kA, lB, kB), file=fp)
            print("#endif", file=fp)

        # Write constraints section: #ifndef FLEXIBLE block
        # In Martini, ALL constraints are inside #ifndef FLEXIBLE
        if hybrid.have_constraints or ifndef_FLEXIBLE_constraints:
            print("\n[ constraints ]", file=fp)
            # Regular constraints (if any exist outside ifndef block)
            if hybrid.have_constraints:
                for p in hybrid.constraints:
                    if len(p) == 3:
                        print("%6d %6d %6d" % (p[0].id, p[1].id, p[2]), file=fp)
                    elif len(p) == 4:
                        print("%6d %6d %6d %8s" % (p[0].id, p[1].id, p[2], p[3]),
                              file=fp)
                    elif len(p) == 5:
                        print("%6d %6d %6d %8s %8s"
                              % (p[0].id, p[1].id, p[2], p[3], p[4]), file=fp)
            if ifndef_FLEXIBLE_constraints:
                print("#ifndef FLEXIBLE", file=fp)
                for c in ifndef_FLEXIBLE_constraints:
                    if len(c) == 5:
                        print("%6d %6d %6d %8s %8s"
                              % (c[0].id, c[1].id, c[2], c[3], c[4]), file=fp)
                    else:
                        print("%6d %6d %6d %8s"
                              % (c[0].id, c[1].id, c[2], c[3]), file=fp)
                print("#endif", file=fp)

        # Write remaining sections
        hybrid.write_pairs(fp)
        hybrid.write_angles(fp, state="AB")
        hybrid.write_dihedrals(fp, state="AB")
        if hybrid.has_exclusions:
            hybrid.write_exclusions(fp)


# ====
# Main
# ====
def main(args):
    # Read structures
    model1 = Model(args.i1)
    model2 = Model(args.i2)

    # Read Martini FF atomtypes
    atomtype_masses = {}
    if args.ff:
        atomtype_masses = parse_martini_ff(args.ff)

    # Read topologies
    if args.itp1:
        topol1 = Topology(args.itp1, is_itp=True, assign_types=False)
        if atomtype_masses:
            assign_masses(topol1, atomtype_masses)
        topol1.write("out_itp1.itp", stateBonded="A", stateTypes="A", stateQ="A")
    if args.itp2:
        topol2 = Topology(args.itp2, is_itp=True, assign_types=False)
        if atomtype_masses:
            assign_masses(topol2, atomtype_masses)
        topol2.write("out_itp2.itp", stateBonded="A", stateTypes="A", stateQ="A")

    # Separate FLEXIBLE bonds/constraints from regular ones
    ifdef_FLEXIBLE_bonds_1 = []
    ifdef_FLEXIBLE_bonds_2 = []
    ifndef_FLEXIBLE_constraints_1 = []
    ifndef_FLEXIBLE_constraints_2 = []
    if args.itp1 and args.itp2:
        # Parse raw itp files to identify ifdef/ifndef FLEXIBLE atom pairs
        ifdef_pairs_1, ifndef_pairs_1 = parse_ifdef_flexible_pairs(args.itp1)
        ifdef_pairs_2, ifndef_pairs_2 = parse_ifdef_flexible_pairs(args.itp2)

        # Split topol1 bonds into regular vs ifdef_FLEXIBLE
        regular_bonds_1 = []
        for b in topol1.bonds:
            pair = (b[0].id, b[1].id)
            if pair in ifdef_pairs_1 or (pair[1], pair[0]) in ifdef_pairs_1:
                ifdef_FLEXIBLE_bonds_1.append(b)
            else:
                regular_bonds_1.append(b)
        topol1.bonds = regular_bonds_1

        # Split topol2 bonds into regular vs ifdef_FLEXIBLE
        regular_bonds_2 = []
        for b in topol2.bonds:
            pair = (b[0].id, b[1].id)
            if pair in ifdef_pairs_2 or (pair[1], pair[0]) in ifdef_pairs_2:
                ifdef_FLEXIBLE_bonds_2.append(b)
            else:
                regular_bonds_2.append(b)
        topol2.bonds = regular_bonds_2

        # Split topol1 constraints into regular vs ifndef_FLEXIBLE
        regular_constraints_1 = []
        for c in topol1.constraints:
            pair = (c[0].id, c[1].id)
            if pair in ifndef_pairs_1 or (pair[1], pair[0]) in ifndef_pairs_1:
                ifndef_FLEXIBLE_constraints_1.append(c)
            else:
                regular_constraints_1.append(c)
        topol1.constraints = regular_constraints_1

        # Split topol2 constraints into regular vs ifndef_FLEXIBLE
        regular_constraints_2 = []
        for c in topol2.constraints:
            pair = (c[0].id, c[1].id)
            if pair in ifndef_pairs_2 or (pair[1], pair[0]) in ifndef_pairs_2:
                ifndef_FLEXIBLE_constraints_2.append(c)
            else:
                regular_constraints_2.append(c)
        topol2.constraints = regular_constraints_2

    # BB mapping
    if args.method == "sequence":
        bb_pairs = map_bb_by_sequence(model1, model2)
    else:
        bb_pairs = map_bb_by_alignment(model1, model2)

    # SC mapping
    sc_pairs = map_sidechains(bb_pairs, model1, model2)

    # Combine all pairs
    all_pairs = bb_pairs + sc_pairs

    # Identify mutated residues (where resname differs between strA and strB)
    mutated_resids = set()
    for a1_bb, a2_bb in bb_pairs:
        res1 = a1_bb.molecule
        res2 = a2_bb.molecule
        if res1.resname.strip() != res2.resname.strip():
            mutated_resids.add(a1_bb.resnr)

    # Remove SC pairs at mutated residues — those SC beads go via dummies
    filtered_pairs = []
    for a1, a2 in all_pairs:
        if a1.name.strip() != "BB" and a1.resnr in mutated_resids:
            continue
        filtered_pairs.append((a1, a2))
    all_pairs = filtered_pairs

    # Write pairs file
    with open(args.o1, "w") as f:
        f.write("# Martini CG bead mapping\n")
        f.write("# atom_idx1  atom_idx2  name1  resname1  resid1  name2  resname2  resid2\n")
        for a1, a2 in all_pairs:
            f.write(
                f"{a1.id:>5d} {a2.id:>5d}  "
                f"{a1.name.strip():<4s} {a1.resname.strip():<5s} {a1.resnr:>4d}  "
                f"{a2.name.strip():<4s} {a2.resname.strip():<5s} {a2.resnr:>4d}\n"
            )

    # Build hybrid topology
    if args.itp1 and args.itp2:
        hybrid = build_hybrid_topology(topol1, topol2, all_pairs)

        # Build B-states for FLEXIBLE bonds and constraints
        hybrid_ifdef_FLEXIBLE_bonds = build_ifdef_FLEXIBLE_bonds_hybrid(
            ifdef_FLEXIBLE_bonds_1, ifdef_FLEXIBLE_bonds_2, hybrid, all_pairs)
        hybrid_ifndef_FLEXIBLE_constraints = build_ifndef_FLEXIBLE_constraints_hybrid(
            ifndef_FLEXIBLE_constraints_1, ifndef_FLEXIBLE_constraints_2, hybrid, all_pairs)

        # Write hybrid topology with ifdef structure
        write_hybrid_itp(hybrid, hybrid_ifdef_FLEXIBLE_bonds, hybrid_ifndef_FLEXIBLE_constraints,
                         "hybrid.itp")

        # Write ffhybrid.itp with DUM atomtype
        write_ffhybrid_itp(args.ffitp)

    # Write hybrid PDB
    mutations = write_hybrid_pdb(bb_pairs, args.o)

    # Log
    with open(args.log, "w") as f:
        f.write(f"Martini Mutate Log\n")
        f.write(f"==================\n")
        f.write(f"Structure 1: {args.i1}\n")
        f.write(f"Structure 2: {args.i2}\n")
        f.write(f"Method: {args.method}\n")
        f.write(f"BB pairs: {len(bb_pairs)}\n")
        f.write(f"SC pairs: {len(sc_pairs)}\n")
        f.write(f"Total pairs: {len(all_pairs)}\n")
        f.write(f"Mutations: {', '.join(mutations) if mutations else 'none'}\n")
        f.write(f"Hybrid PDB: {args.o}\n")

    print(f"Mapped {len(bb_pairs)} BB pairs and {len(sc_pairs)} SC pairs.")
    print(f"Total: {len(all_pairs)} pairs written to {args.o1}")
    if mutations:
        print(f"Mutations detected: {', '.join(mutations)}")
    print(f"Hybrid structure written to {args.o}")


def entry_point():
    args = parse_options()
    main(args)


if __name__ == "__main__":
    entry_point()
