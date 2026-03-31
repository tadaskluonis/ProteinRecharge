# ProteinRecharge

<p align="center">
  <img src="images/ProteinRecharge_transparent.png" alt="ProteinRecharge logo" width="500"/>
</p>

ProteinRecharge is a small utility for recharging protein surfaces after binder design pipelines produce inappropriately charged binders. Requires LigandMPNN for sequence redesign.

## Quick start ✅

1. Clone this repository and install LigandMPNN:

```bash
git clone https://github.com/profdocpizza/ProteinRecharge.git
cd ProteinRecharge

```

2. Follow [LigandMPNN](https://github.com/dauparas/LigandMPNN/) instructions to install it and download weights 

3. Edit `examples/binder_recharge.yaml` and set the following fields for your environment:

```yaml
ligandmpnn_run_py: "/path/to/LigandMPNN/run.py"   # path to the LigandMPNN runner
ligandmpnn_conda_env: "ligandmpnn_env"          # or "" to use the current env
omit_AA: "C"                                    # default; override to allow cysteine
```

4. Check your configuration by running:

```bash
python check.py examples/binder_recharge.yaml
```
It will use biopython's SASA to identify burried residues. Also, it will account for user-specified residues as well as any residues in the interface with user-specified chains. To inspect the selections, open generated `.cif` file in Nano Viewer (VScode extention) and color by PLDDT: 

<p align="center">
  <img src="images/colored_transparent.png" alt="ProteinRecharge logo" width="400"/>
</p>

<table>
  <tr>
    <th>Category</th>
    <th>Description</th>
    <th>Color</th>
  </tr>

  <tr>
    <td><b>REDESIGNABLE</b></td>
    <td>Surface-exposed residues selected for redesign</td>
    <td>
      <img src="https://placehold.co/15x15/f6ed12/f6ed12.png" alt="redesignable">
    </td>
  </tr>

  <tr>
    <td><b>PROTECTED</b></td>
    <td>Interface or user-protected residues (excluded from redesign)</td>
    <td>
      <img src="https://placehold.co/15x15/10cff1/10cff1.png" alt="protected">
    </td>
  </tr>

  <tr>
    <td><b>BURIED</b></td>
    <td>Non-exposed residues on the redesign chain</td>
    <td>
      <img src="https://placehold.co/15x15/106dff/106dff.png" alt="buried">
    </td>
  </tr>

  <tr>
    <td><b>NON_BINDER</b></td>
    <td>Residues on non-redesign chains</td>
    <td><img src="https://placehold.co/15x15/ef821e/ef821e.png" alt="non_binder"></td>
  </tr>
</table>




4. Run the pipeline:

```bash
python recharge.py examples/binder_recharge.yaml
```

## Usage

Basic usage:

```bash
# Check configuration and selected residues
python check.py examples/binder_recharge.yaml

# Run the redesign pipeline
python recharge.py examples/binder_recharge.yaml
```

Where to look for outputs:

- Generated design inputs and logs will be placed alongside the provided YAML (see examples/ for sample outputs).
- A `.cif` file visualizing selected residues will be produced by the check step.

### Full example configuration

Below is a full example YAML configuration covering commonly used options. Copy this to `examples/binder_recharge_full.yaml` and edit paths for your environment.

```yaml
# Full example config for ProteinRecharge
# Edit the paths and settings to match your environment

# Path to LigandMPNN runner (run.py from LigandMPNN repo)
ligandmpnn_run_py: "/path/to/LigandMPNN/run.py"
# Conda environment name to activate for LigandMPNN, or empty string to use current env
ligandmpnn_conda_env: "ligandmpnn_env"

# Input structures: specify individual files and/or directories
inputs:
  - input_pdb_path: "inputs/my_binder.pdb"

# Chain on which redesign should occur (single letter or chain id)
redesign_chain: "A"

# Residues to avoid (explicit list of residue numbers on redesign chain)
avoid_residues: []
# Interface residues with specified chains to avoid (e.g. ["B", "C"])
avoid_interface_residues_with_chain: []

# Selection method for redesignable residues: "sasa" or "contact"
selection_method: "sasa"
# If using sasa, fraction threshold (0.0 means all residues are considered exposed)
sasa_fraction: 0.15

# Omit amino acids from redesign (string, e.g. "C" to omit cysteine), default is "CX"
omit_AA: "C"

# Charge target optimization (optional)
desired_chain_charge: 0.0     # desired net charge for redesign chain
charge_tolerance: 1.0         # allowed deviation from desired_chain_charge

# Output directory
out_folder: "outputs/my_binder_recharge"
```

### Small config examples

Below are two compact config examples showing common variants.

1) Full-chain redesign on target chain (make every residue on the redesign chain redesignable):

```yaml
# full-chain redesign: set sasa_fraction to 0.0 and disable explicit avoids
inputs:
  - input_pdb_path: "inputs/my_binder.pdb"
redesign_chain: "A"
ligandmpnn_run_py: "/path/to/LigandMPNN/run.py"
ligandmpnn_conda_env: "ligandmpnn_env"
selection_method: "sasa"
sasa_fraction: 0.0   # 0.0 -> treat all residues on redesign_chain as exposed
avoid_residues: []
avoid_interface_residues_with_chain: []
omit_AA: ""          # empty string allows cysteines and all residues
out_folder: "outputs/full_chain_recharge"
```

Notes: keep avoid lists empty to allow redesign of the whole chain; set omit_AA to empty string to allow all amino acids.

2) Running the pipeline over many structures in a directory

Use the `inputs` list in the config to specify both individual files and directories. The pipeline will process all `.pdb` and `.cif` files found in listed directories.

```yaml
# Specify multiple inputs: single files and/or directories
inputs:
  - input_pdb_path: "inputs/binder_1.pdb"
  - input_pdb_path: "inputs/binder_2.pdb"
  - input_pdb_dir: "inputs/batch_designs/"  # Processes all *.pdb and *.cif in this directory

# Common settings (applied to all input structures)
ligandmpnn_run_py: "/path/to/LigandMPNN/run.py"
ligandmpnn_conda_env: "ligandmpnn_env"
redesign_chain: "A"
selection_method: "sasa"
sasa_fraction: 0.15
out_folder: "outputs/batch_recharge"
```

Run:

```bash
python recharge.py examples/batch_config.yaml
```

The pipeline will:
- Collect all structures from specified paths and directories.
- Apply the same redesign parameters to each structure.
- Generate outputs in `out_folder` for each input file.


## Chain/context behavior (important)

- ProteinRecharge redesigns only `redesign_chain`, while keeping other protein chains as fixed context in LigandMPNN (`--chains_to_design`).
- This preserves multi-chain structural context during design.
- For multi-chain FASTA outputs (`A_seq:B_seq:...`), ProteinRecharge now maps and analyzes the configured `redesign_chain` segment for:
  - sequence alignment/mutation reporting
  - charge calculations and charge-target optimization

## Full-chain redesign on target chain

If you want **all residues** on the redesign chain to be redesignable (not only exposed residues), set:

```yaml
selection_method: "sasa"
sasa_fraction: 0.0
```

Keep `avoid_residues: []` and related avoid settings disabled if you do not want any explicit exclusions.
