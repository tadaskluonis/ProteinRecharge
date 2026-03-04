#!/usr/bin/env python3
import argparse
import yaml
import os
from pathlib import Path
from Bio.PDB import PDBParser, MMCIFParser, MMCIFIO, Select, NeighborSearch
from recharge import load_config, get_input_files, get_protected_residues, select_residues_sasa, select_residues, estimate_charge_range, Config

# For quick hyperparameter testing, show default SASA params



class ChargeColorMapper(Select):
    def __init__(self, redesign_res, protected_res, redesign_chain):
        self.redesign_res = set(redesign_res)
        self.protected_res = set(protected_res)
        self.redesign_chain = redesign_chain

    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        res_key = f"{chain_id}{residue.id[1]}"
        
        # B-factor Legend for Mol* Coloring:
        # 100.0: BURRIED (Buried residues) -> PLDDT = blue
        # 80.0:  PROTECTED (Interface / Avoided residues) -> PLDDT = blue
        # 60.0:  REDESIGNABLE (Selected as exposed) -> PLDDT = yellow
        # 40.0:  Non-Binder (other chains) -> PLDDT = orange
        
        # Default: Non-binder (other chains)
        if chain_id != self.redesign_chain:
            val = 40.0
        else:
            # Determine explicit membership to ensure protected overrides buried
            is_protected = res_key in self.protected_res
            is_redesignable = res_key in self.redesign_res

            if is_protected:
                # Protected residues (interface or user-specified avoid) take precedence
                val = 80.0
            elif is_redesignable:
                # Selected as exposed
                val = 60.0
            else:
                # Neither protected nor redesignable -> buried/other
                val = 100.0
            
        for atom in residue:
            atom.set_bfactor(val)
        return 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()

    # 1. Load Config
    cfg = load_config(args.config)

    # 2. Setup Check Directory
    check_dir = Path(cfg.out_folder) / "check"
    check_dir.mkdir(parents=True, exist_ok=True)

    # 3. Identify Files
    all_files = get_input_files(cfg)

    if not all_files:
        print(f"No files found.")
        return

    print(f"Processing {len(all_files)} files for selection check...")

    # 4. Loop through all structures
    for in_file in all_files:
        try:
            p = MMCIFParser(QUIET=True) if in_file.suffix == '.cif' else PDBParser(QUIET=True)
            struct = p.get_structure(in_file.stem, str(in_file))
            
            # Run selection logic based on config
            if getattr(cfg, 'selection_method', 'contact') == 'sasa':
                print(f"Using SASA selection: fraction={cfg.sasa_fraction}, threshold={cfg.sasa_threshold}")  # threshold=None means only fraction is used
                redesignable = select_residues_sasa(struct, cfg)
            else:
                print(f"Using contact-based selection: radius={cfg.contact_radius}, min_contacts={cfg.min_contacts_for_buried}")
                redesignable = select_residues(struct, cfg)

            # Build protected set using shared utility
            protected = get_protected_residues(struct, cfg)

            # Calculate excluded (buried or otherwise not redesignable)
            all_res = [f"{cfg.redesign_chain}{r.id[1]}" for r in struct[0][cfg.redesign_chain] if 'CA' in r]
            buried_or_other = [r for r in all_res if r not in redesignable and r not in protected]

            # Counts per category
            redesignable_count = len(redesignable)
            protected_count = len(protected)
            buried_count = len(buried_or_other)
            # Count residues on other chains (non-binder)
            non_binder_count = 0
            for ch in struct[0]:
                if ch.id == cfg.redesign_chain:
                    continue
                non_binder_count += sum(1 for r in ch if 'CA' in r)
                
            # Estimated charge range
            current_charge, min_possible_charge, max_possible_charge = estimate_charge_range(struct, cfg, redesignable)

            plddt_colors = {
                "very_low_plddt": "0xef821e",
                "low_plddt": "0xf6ed12",
                "confident_plddt": "0x10cff1",
                "very_high_plddt": "0x106dff",
            }
            mapping = {
                "REDESIGNABLE": plddt_colors["low_plddt"],
                "PROTECTED": plddt_colors["confident_plddt"],
                "BURIED": plddt_colors["very_high_plddt"],
                "NON_BINDER": plddt_colors["very_low_plddt"],
            }

            def hex_to_rgb(h):
                h = h.lower().replace("0x", "")
                return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)

            def color_text(text, hexcode):
                r, g, b = hex_to_rgb(hexcode)
                return f"\x1b[38;2;{r};{g};{b}m{text}\x1b[0m"

            parts = []
            for cat, cnt in [("REDESIGNABLE", redesignable_count), ("PROTECTED", protected_count), ("BURIED", buried_count), ("NON_BINDER", non_binder_count)]:
                parts.append(f"{color_text(cat, mapping[cat])}: {cnt}")

            print("Counts -> " + ", ".join(parts))
            print(f"Current charge: {current_charge} | Possible range: [{min_possible_charge}, {max_possible_charge}]")
            print(f"Target charge: {cfg.desired_chain_charge}")
            
            if cfg.desired_chain_charge < min_possible_charge or cfg.desired_chain_charge > max_possible_charge:
                print(f"⚠️ WARNING: Desired charge {cfg.desired_chain_charge} is outside the possible range based on the number of redesignable residues.")
                closest = max(min_possible_charge, min(max_possible_charge, cfg.desired_chain_charge))
                print(f"⚠️ Falling back to closest available value: {closest}")
            else:
                print(f"✅ Target charge {cfg.desired_chain_charge} is achievable.")

            # Save check file
            io = MMCIFIO()
            io.set_structure(struct)
            out_path = check_dir / f"{in_file.stem}_check.cif"
            io.save(str(out_path), ChargeColorMapper(redesignable, protected, cfg.redesign_chain))

        except Exception as e:
            print(f"Error processing {in_file.name}: {e}")

    print(f"\nAll check files saved to: {check_dir}")
    print(f"Open {out_path} in Nano Viewer (VScode extention) and color by PLDDT to visualise the selections.")

if __name__ == "__main__":
    main()