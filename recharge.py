#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import yaml
import shutil
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Optional

# Biopython
from Bio.PDB import PDBParser, MMCIFParser, MMCIFIO, Select, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley

# Analysis
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ANSI Colors
class Colors:
    RED = "\033[91m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    CYAN = "\033[96m"
    MAGENTA = "\033[95m"
    RESET = "\033[0m"

AA_COLORS = {
    'D': Colors.RED, 'E': Colors.RED,  # Negative
    'K': Colors.BLUE, 'R': Colors.BLUE, # Positive
    'H': Colors.CYAN,                   # Histidine (Polar/Positive)
    'N': Colors.GREEN, 'Q': Colors.GREEN, 'S': Colors.GREEN, 'T': Colors.GREEN, 'Y': Colors.GREEN, # Polar
    'A': Colors.YELLOW, 'I': Colors.YELLOW, 'L': Colors.YELLOW, 'M': Colors.YELLOW, 'F': Colors.YELLOW, 'W': Colors.YELLOW, 'V': Colors.YELLOW, 'G': Colors.YELLOW, 'P': Colors.YELLOW, # Non-polar
    'C': Colors.MAGENTA # Cysteine
}

def color_aa(aa):
    return f"{AA_COLORS.get(aa, '')}{aa}{Colors.RESET}"

def print_alignment(seq1, seq2, width=60):
    """Prints a simple alignment of two sequences with colored mutations."""
    for i in range(0, len(seq1), width):
        chunk1 = seq1[i:i+width]
        chunk2 = seq2[i:i+width]
        
        # Original sequence with colored positions that will be changed
        res1 = []
        for a, b in zip(chunk1, chunk2):
            if a != b:
                res1.append(color_aa(a))
            else:
                res1.append(a)
        print(f"Orig: {''.join(res1)}")
        
        # New sequence with colored mutations
        res2 = []
        for a, b in zip(chunk1, chunk2):
            if a != b:
                res2.append(color_aa(b))
            else:
                res2.append(b)
        print(f"New : {''.join(res2)}")
        print()

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
@dataclass
class Config:
    out_folder: str
    ligandmpnn_run_py: str
    ligandmpnn_conda_env: str
    
    input_structure: Optional[str] = None
    input_dir: Optional[str] = None
    
    redesign_chain: str = "A"
    avoid_residues: List[str] = field(default_factory=list)
    avoid_residues_distance: float = 0.0
    avoid_interface_residues_with_chain: List[str] = field(default_factory=list)
    avoid_interface_residues_with_chain_distance: float = 0.0
    
    avoid_salt_bridge_residues: bool = False
    avoid_H_bonding_residues: bool = False
    
    desired_chain_charge: float = 0.0
    tolerance: float = 1.0
    initial_bias: float = 0.0

    # SASA-based selection parameters
    sasa_fraction: float = 0.25        # fraction of max residue SASA to consider exposed
    sasa_threshold: Optional[float] = None  # absolute SASA (Å^2) threshold; set to a number to use an absolute cutoff
    
    seed: int = 42
    batch_size: int = 20 
    contact_radius: float = 6.0
    min_contacts_for_buried: int = 11
    max_evals: int = 10

def load_config(path):
    with open(path) as f:
        data = yaml.safe_load(f)
    return Config(**{k: v for k, v in data.items() if k in Config.__dataclass_fields__})

class RedesignedColorMapper(Select):
    def __init__(self, mutated_res_ids, redesign_chain):
        self.mutated_res_ids = mutated_res_ids
        self.redesign_chain = redesign_chain

    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        if chain_id != self.redesign_chain:
            val = 40.0 # NON_BINDER color (orange)
        else:
            if residue.id[1] in self.mutated_res_ids:
                val = 60.0 # REDESIGNABLE color (yellow) -> REDESIGNED
            else:
                val = 100.0 # BURIED color (blue) -> IDENTICAL

        for atom in residue:
            atom.set_bfactor(val)
        return 1

# ------------------------------------------------------------------
# Geometry & Selection
# ------------------------------------------------------------------

def find_salt_bridges(model, redesign_chain, cutoff=4.0):
    """Identifies residues on the redesign_chain involved in salt bridges."""
    salt_bridges = set()
    
    neg_res = {'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2']}
    pos_res = {'ARG': ['NH1', 'NH2', 'NE'], 'LYS': ['NZ'], 'HIS': ['ND1', 'NE2']}
    
    pos_atoms = []
    for chain in model:
        for res in chain:
            rname = res.get_resname()
            if rname in pos_res:
                for atom_id in pos_res[rname]:
                    if atom_id in res:
                        pos_atoms.append(res[atom_id])
                        
    if not pos_atoms:
        return salt_bridges
        
    ns = NeighborSearch(pos_atoms)
    
    for chain in model:
        for res in chain:
            rname = res.get_resname()
            if rname in neg_res:
                for atom_id in neg_res[rname]:
                    if atom_id in res:
                        close_pos_atoms = ns.search(res[atom_id].get_coord(), cutoff, level='A')
                        if close_pos_atoms:
                            # We found a salt bridge! Process both sides.
                            if chain.id == redesign_chain:
                                salt_bridges.add(f"{chain.id}{res.id[1]}")
                            for pa in close_pos_atoms:
                                if pa.get_parent().get_parent().id == redesign_chain:
                                    salt_bridges.add(f"{redesign_chain}{pa.get_parent().id[1]}")
                                    
    return salt_bridges

def find_h_bonds(model, redesign_chain, cutoff=3.5):
    """Identifies residues on the redesign_chain whose side-chains form hydrogen bonds."""
    h_bonds = set()
    
    sc_polar = {
        'ARG': ['NH1', 'NH2', 'NE'], 'LYS': ['NZ'], 'HIS': ['ND1', 'NE2'],
        'ASN': ['ND2', 'OD1'], 'GLN': ['NE2', 'OE1'],
        'SER': ['OG'], 'THR': ['OG1'], 'TYR': ['OH'], 'TRP': ['NE1'],
        'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2']
    }
    
    neg_res = {'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2']}
    pos_res = {'ARG': ['NH1', 'NH2', 'NE'], 'LYS': ['NZ'], 'HIS': ['ND1', 'NE2']}
    
    def is_salt_bridge(res1, atom1, res2, atom2):
        r1, a1 = res1.get_resname(), atom1.get_name()
        r2, a2 = res2.get_resname(), atom2.get_name()
        
        is_neg1 = r1 in neg_res and a1 in neg_res.get(r1, [])
        is_pos1 = r1 in pos_res and a1 in pos_res.get(r1, [])
        is_neg2 = r2 in neg_res and a2 in neg_res.get(r2, [])
        is_pos2 = r2 in pos_res and a2 in pos_res.get(r2, [])
        
        return (is_neg1 and is_pos2) or (is_pos1 and is_neg2)

    all_polar_atoms = []
    for chain in model:
        for res in chain:
            for atom in res:
                if atom.element in ['N', 'O']:
                    all_polar_atoms.append(atom)
                    
    if not all_polar_atoms:
        return h_bonds
        
    ns = NeighborSearch(all_polar_atoms)
    
    for chain in model:
        for res in chain:
            rname = res.get_resname()
            if rname in sc_polar:
                for atom_id in sc_polar[rname]:
                    if atom_id in res:
                        close_atoms = ns.search(res[atom_id].get_coord(), cutoff, level='A')
                        for pa in close_atoms:
                            other_res = pa.get_parent()
                            if other_res != res:
                                # Ignore salt bridges
                                if is_salt_bridge(res, res[atom_id], other_res, pa):
                                    continue
                                    
                                # Current residue has a side-chain atom involved in an H-bond
                                if chain.id == redesign_chain:
                                    h_bonds.add(f"{chain.id}{res.id[1]}")
                                    
                                # Check if the other residue is ALSO using a sidechain atom for this H-bond
                                # If it's a backbone atom like 'N' or 'O', we don't protect it
                                other_rname = other_res.get_resname()
                                if other_rname in sc_polar and pa.get_name() in sc_polar[other_rname]:
                                    if other_res.get_parent().id == redesign_chain:
                                        h_bonds.add(f"{redesign_chain}{other_res.id[1]}")
                                    
    return h_bonds

def get_residue_key(res, chain_id):
    return f"{chain_id}{res.id[1]}"

def select_residues(structure, cfg: Config):
    model = next(structure.get_models())
    atoms = list(model.get_atoms())
    ns = NeighborSearch(atoms)
    
    target_chain = model[cfg.redesign_chain]
    surface_candidates = []
    
    def _center_coord(res):
        # Prefer CB for side-chain-centered neighborhood; fallback to CA (e.g., glycine)
        if 'CB' in res:
            return res['CB'].get_coord()
        if 'CA' in res:
            return res['CA'].get_coord()
        return None

    for res in target_chain:
        center = _center_coord(res)
        if center is None:
            continue
        neighbors = ns.search(center, cfg.contact_radius, level='R')
        neighbors = [n for n in neighbors if n != res]
        if len(neighbors) < cfg.min_contacts_for_buried:
            surface_candidates.append(get_residue_key(res, cfg.redesign_chain))

    excluded = set(cfg.avoid_residues)

    if getattr(cfg, 'avoid_salt_bridge_residues', False):
        excluded.update(find_salt_bridges(model, cfg.redesign_chain))
        
    if getattr(cfg, 'avoid_H_bonding_residues', False):
        excluded.update(find_h_bonds(model, cfg.redesign_chain))

    if cfg.avoid_residues_distance > 0:
        for seed_id in cfg.avoid_residues:
            try:
                r_obj = next(r for r in target_chain if get_residue_key(r, cfg.redesign_chain) == seed_id)
                # Use CB coordinate when available, otherwise CA (glycine)
                coord = None
                if 'CB' in r_obj:
                    coord = r_obj['CB'].get_coord()
                elif 'CA' in r_obj:
                    coord = r_obj['CA'].get_coord()
                if coord is None:
                    continue
                close = ns.search(coord, cfg.avoid_residues_distance, level='R')
                for f in close:
                    if f.get_parent().id == cfg.redesign_chain:
                        excluded.add(get_residue_key(f, cfg.redesign_chain))
            except StopIteration:
                pass

    interface_chains = cfg.avoid_interface_residues_with_chain
    if isinstance(interface_chains, str):
        interface_chains = [interface_chains]
    
    if interface_chains:
        dist = max(0.1, cfg.avoid_interface_residues_with_chain_distance)
        for target_id in interface_chains:
            if target_id in model:
                other_atoms = list(model[target_id].get_atoms())
                for atom in other_atoms:
                    for r in ns.search(atom.get_coord(), dist, level='R'):
                        if r.get_parent().id == cfg.redesign_chain:
                            excluded.add(get_residue_key(r, cfg.redesign_chain))

    return sorted(list(set(surface_candidates) - excluded), key=lambda x: int(x[1:]))


def select_residues_sasa(structure, cfg: Config):
    """Select residues by solvent exposure (SASA) using Biopython's Shrake-Rupley.

    Selection rules:
      - Compute per-residue SASA (sum of atom.sasa)
      - Select residues with SASA >= max_sasa * cfg.sasa_fraction OR SASA >= cfg.sasa_threshold
      - Apply `cfg.avoid_residues` and `cfg.avoid_residues_distance` exclusions and interface exclusions

    Returns a sorted list of residue keys (e.g., ['A10', 'A11', ...]).
    """
    model = next(structure.get_models())

    # Compute SASA at atom level
    sr = ShrakeRupley()
    sr.compute(structure)

    target_chain = model[cfg.redesign_chain]

    sasa_by_res = {}
    for res in target_chain:
        atoms = list(res.get_atoms())
        if not atoms:
            continue
        total = 0.0
        for a in atoms:
            total += getattr(a, 'sasa', 0.0)
        sasa_by_res[res] = total

    if not sasa_by_res:
        return []

    # Determine selection set based on fraction and absolute threshold
    sorted_res = sorted(sasa_by_res.items(), key=lambda kv: kv[1], reverse=True)
    max_sasa = sorted_res[0][1]

    selected = set()
    for res, val in sorted_res:
        meets_fraction = (max_sasa > 0 and val >= max_sasa * cfg.sasa_fraction)
        meets_abs = (cfg.sasa_threshold is not None and val >= cfg.sasa_threshold)
        if meets_fraction or meets_abs:
            selected.add(get_residue_key(res, cfg.redesign_chain))

    # Apply exclusions (avoid_residues and distance-based and interface)
    excluded = set(cfg.avoid_residues)
    
    if getattr(cfg, 'avoid_salt_bridge_residues', False):
        excluded.update(find_salt_bridges(model, cfg.redesign_chain))
        
    if getattr(cfg, 'avoid_H_bonding_residues', False):
        excluded.update(find_h_bonds(model, cfg.redesign_chain))

    if cfg.avoid_residues_distance > 0:
        for seed_id in cfg.avoid_residues:
            try:
                r_obj = next(r for r in target_chain if get_residue_key(r, cfg.redesign_chain) == seed_id)
                coord = None
                if 'CB' in r_obj:
                    coord = r_obj['CB'].get_coord()
                elif 'CA' in r_obj:
                    coord = r_obj['CA'].get_coord()
                if coord is None:
                    continue
                ns = NeighborSearch(list(model.get_atoms()))
                close = ns.search(coord, cfg.avoid_residues_distance, level='R')
                for f in close:
                    if f.get_parent().id == cfg.redesign_chain:
                        excluded.add(get_residue_key(f, cfg.redesign_chain))
            except StopIteration:
                pass

    interface_chains = cfg.avoid_interface_residues_with_chain
    if isinstance(interface_chains, str):
        interface_chains = [interface_chains]

    if interface_chains:
        dist = max(0.1, cfg.avoid_interface_residues_with_chain_distance)
        for target_id in interface_chains:
            if target_id in model:
                other_atoms = list(model[target_id].get_atoms())
                ns = NeighborSearch(list(model.get_atoms()))
                for atom in other_atoms:
                    for r in ns.search(atom.get_coord(), dist, level='R'):
                        if r.get_parent().id == cfg.redesign_chain:
                            excluded.add(get_residue_key(r, cfg.redesign_chain))

    final = sorted(list(selected - excluded), key=lambda x: int(x[1:]))
    return final

# ------------------------------------------------------------------
# LigandMPNN Execution & Analysis
# ------------------------------------------------------------------

def run_ligandmpnn(cfg: Config, out_dir, residues, bias_val, extra_args=None):
    bias_str = f"D:{-bias_val:.3f},E:{-bias_val:.3f},K:{bias_val:.3f},R:{bias_val:.3f}"
    res_str = " ".join(residues)
    
    cmd = [
        "python", cfg.ligandmpnn_run_py,
        "--pdb_path", f'"{Path(cfg.input_structure).resolve()}"',
        "--out_folder", f'"{Path(out_dir).resolve()}"',
        "--redesigned_residues", f'"{res_str}"',
        "--bias_AA", f'"{bias_str}"',
        "--batch_size", str(cfg.batch_size),
        "--number_of_batches", "1", 
        "--model_type", "soluble_mpnn",
        "--seed", str(cfg.seed),
        "--verbose", "0"
    ]
    if extra_args: cmd.extend(extra_args)

    # Build the command string; only activate conda if an env was provided
    cmd_str = " ".join(cmd)
    if cfg.ligandmpnn_conda_env:
        full_cmd = f'source $(conda info --base)/etc/profile.d/conda.sh && conda activate {cfg.ligandmpnn_conda_env} && {cmd_str}'
    else:
        full_cmd = cmd_str

    # Run in the directory where the LigandMPNN runner lives so relative model paths resolve
    cwd_dir = Path(cfg.ligandmpnn_run_py).parent
    subprocess.run(full_cmd, shell=True, check=True, executable="/bin/bash", cwd=str(cwd_dir))

def get_native_sequence(step_dir):
    p = Path(step_dir)
    seq_dir = p / "seqs"
    if not seq_dir.exists(): return None
    
    fa_file = next(seq_dir.glob("*.fa"), None)
    if not fa_file: return None
    
    lines = fa_file.read_text().splitlines()
    if len(lines) >= 2:
        # Native is always index 1 (the first sequence)
        return lines[1].split(":")[0]
    return None

def process_iteration_results(step_dir):
    results = []
    p = Path(step_dir)
    seq_dir = p / "seqs"
    if not seq_dir.exists(): return []

    for f in seq_dir.glob("*.fa"):
        lines = f.read_text().splitlines()
        # skip index 0 as it is the native sequence in LigandMPNN outputs
        for i in range(3, len(lines), 2): 
            header, seq = lines[i-1], lines[i]
            s = seq.split(":")[0]
            charge = (s.count("K") + s.count("R")) - (s.count("D") + s.count("E"))
            
            name = header.split(",")[0].replace(">", "").strip()
            pdb_path = p / "backbones" / f"{name}.pdb"
            
            results.append({"name": name, "path": str(pdb_path.absolute()), "charge": charge, "sequence": s})
    
    if results:
        pd.DataFrame(results).to_csv(p / "design_stats.csv", index=False)
    return results

def update_plot(history, out_folder):
    data = []
    for s in history:
        for c in s['charges']:
            data.append({"Bias": s['bias'], "Net Charge": c})
    
    if not data: return
    df = pd.DataFrame(data)
    plt.figure(figsize=(8, 5))
    sns.scatterplot(data=df, x="Bias", y="Net Charge", alpha=0.6)
    plt.axhline(y=history[0].get('desired', 0), color='r', linestyle='--', label='Target')
    plt.title(f"Charge Optimization Path")
    plt.savefig(Path(out_folder) / "charge_plot.png")
    plt.close()

# ------------------------------------------------------------------
# Main Logic
# ------------------------------------------------------------------

def run_recharge_pipeline(struct_path: Path, cfg: Config, selection_dir: Path):
    print(f"\n>>> Processing: {struct_path.name}")
    
    protein_out = Path(cfg.out_folder) / struct_path.stem
    protein_out.mkdir(parents=True, exist_ok=True)

    parser = MMCIFParser(QUIET=True) if struct_path.suffix == '.cif' else PDBParser(QUIET=True)
    struct = parser.get_structure(struct_path.stem, str(struct_path))
    
    # Update current path for the subprocess call
    cfg.input_structure = str(struct_path)

    redesign_res = select_residues(struct, cfg)
    print(f"Selected {len(redesign_res)} residues for redesign.")
    
    bias, step_size, direction, history = cfg.initial_bias, 0.2, 0, []
    best_overall_results = []

    for i in range(cfg.max_evals):
        step_dir = protein_out / f"step_{i}_bias_{bias:.3f}"
        step_dir.mkdir(exist_ok=True)
        
        run_ligandmpnn(cfg, str(step_dir), redesign_res, bias)
        current_results = process_iteration_results(step_dir)
        
        if not current_results:
            continue
            
        best_overall_results = current_results # Store the most recent pool
        charges = [r['charge'] for r in current_results]
        mean_q = sum(charges)/len(charges)
        
        history.append({'bias': bias, 'mean_charge': mean_q, 'charges': charges, 'desired': cfg.desired_chain_charge})
        update_plot(history, protein_out)

        print(f"Iteration {i} | Bias: {bias:.3f} | Mean Charge: {mean_q:.2f}")

        if abs(mean_q - cfg.desired_chain_charge) <= cfg.tolerance:
            print("Tolerance reached.")
            break
        
        new_dir = 1.0 if mean_q < cfg.desired_chain_charge else -1.0
        if i > 0 and new_dir != direction: 
            step_size *= 0.5
        direction = new_dir
        bias += direction * step_size

    with open(protein_out / "summary.json", "w") as f:
        json.dump(history, f, indent=2)

    # Selection: Find sequence closest to desired charge in the final pool
    if best_overall_results:
        best_pick = sorted(best_overall_results, key=lambda x: abs(x['charge'] - cfg.desired_chain_charge))[0]
        
        out_fasta = selection_dir / f"{struct_path.stem}_recharged.fasta"
        out_pdb = selection_dir / f"{struct_path.stem}_recharged.pdb"
        
        # Write Fasta
        with open(out_fasta, "w") as f:
            f.write(f">{struct_path.stem}_charge_{best_pick['charge']}\n{best_pick['sequence']}\n")
        
        # Get native for alignment
        # Usually step_0 is run first, so we might find native there or any step.
        native_seq = get_native_sequence(protein_out / f"step_0_bias_{cfg.initial_bias:.3f}")
        if native_seq:
            print(f"\nAlignment for {struct_path.name}:")
            new_seq = best_pick['sequence']
            print_alignment(native_seq, new_seq)
            
            # Map sequence mutations to structure residue IDs
            # LigandMPNN sequence represents residues with CA atoms in the redesign chain
            chain = struct[0][cfg.redesign_chain]
            ca_residues = [res for res in chain if 'CA' in res]
            
            mutated_res_ids = set()
            for i, (orig_aa, new_aa) in enumerate(zip(native_seq, new_seq)):
                if orig_aa != new_aa and i < len(ca_residues):
                    mutated_res_ids.add(ca_residues[i].id[1])
            
            check_cif_path = selection_dir / f"{struct_path.stem}_check_redesigned.cif"
            io = MMCIFIO()
            io.set_structure(struct)
            io.save(str(check_cif_path), RedesignedColorMapper(mutated_res_ids, cfg.redesign_chain))
            print(f"Saved redesigned structural check map to: {check_cif_path}")

        # Copy PDB
        if Path(best_pick['path']).exists():
            shutil.copy(best_pick['path'], out_pdb)
            
        print(f"Saved best design (Charge {best_pick['charge']}) to selection folder.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    
    cfg = load_config(args.config)
    out_root = Path(cfg.out_folder)
    out_root.mkdir(parents=True, exist_ok=True)
    
    selection_dir = out_root / "selection"
    selection_dir.mkdir(exist_ok=True)

    input_files = []
    if cfg.input_dir:
        p = Path(cfg.input_dir)
        input_files = list(p.glob("*.pdb")) + list(p.glob("*.cif"))
    elif cfg.input_structure:
        input_files = [Path(cfg.input_structure)]

    for struct_path in input_files:
        try:
            run_recharge_pipeline(struct_path, cfg, selection_dir)
        except Exception as e:
            print(f"Error processing {struct_path.name}: {e}")

if __name__ == "__main__":
    main()