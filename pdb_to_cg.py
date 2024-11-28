import argparse
from Bio.PDB import PDBParser, PDBIO, Select

COARSE_GRAIN_ATOMS = {
    "backbone": ["P", "C4'"],
    "purine": ["N9", "C2", "C6"],
    "pyrimidine": ["N1", "C2", "C4"]
}

ALLOWED_RESIDUES = ["A", "G", "C", "U", "T"]

class CoarseGrainSelect(Select):
    def accept_atom(self, atom):
        residue = atom.get_parent()
        residue_name = residue.get_resname()
        atom_name = atom.get_name()

        if residue_name not in ALLOWED_RESIDUES:
            return False

        if atom_name in COARSE_GRAIN_ATOMS["backbone"]:
            return True

        if residue_name in ["A", "G"] and atom_name in COARSE_GRAIN_ATOMS["purine"]:
            return True

        if residue_name in ["C", "U", "T"] and atom_name in COARSE_GRAIN_ATOMS["pyrimidine"]:
            return True
        return False


def parse_arguments():
    parser = argparse.ArgumentParser(description="Conversion of full-atomic structure to coarse-grained representation.")
    parser.add_argument("--input", required=True, help="Path to the PDB input file.")
    parser.add_argument("--output", required=True, help="Path to the PDB output file.")
    return parser.parse_args()


def convert_to_coarse_grain(input_file, output_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", input_file)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, select=CoarseGrainSelect())
    print(f"The coarse-grained structure has been saved to the file: {output_file}")


if __name__ == "__main__":
    args = parse_arguments()
    convert_to_coarse_grain(args.input, args.output)
