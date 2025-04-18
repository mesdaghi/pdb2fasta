from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO import FastaIO
import sys

# Mapping from 3-letter to 1-letter amino acid codes
three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'MSE': 'M',  # Handle selenomethionine
}

def pdb_to_fasta(pdb_file, output_fasta):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    model = structure[0]
    with open(output_fasta, "w") as out_handle:
        fasta_writer = FastaIO.FastaWriter(out_handle, wrap=None)  # No line wrapping
        fasta_writer.write_header()
        for chain in model:
            sequence = ''
            seen_residues = set()
            for residue in chain:
                if residue.id[0] != ' ':  # Skip heteroatoms and waters
                    continue
                res_id = residue.id[1]
                if res_id in seen_residues:
                    continue
                seen_residues.add(res_id)
                resname = residue.get_resname()
                one_letter = three_to_one.get(resname, 'X')
                sequence += one_letter
            record = SeqRecord(Seq(sequence), id=chain.id, description=f"Chain {chain.id}")
            fasta_writer.write_record(record)
        fasta_writer.write_footer()

# Example usage:
# pdb_to_fasta("example.pdb", "output.fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pdb_to_fasta.py input.pdb output.fasta")
    else:
        pdb_to_fasta(sys.argv[1], sys.argv[2])

