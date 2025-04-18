
# PDB to FASTA Converter

This Python script converts a Protein Data Bank (PDB) file into a FASTA format file by extracting the amino acid sequence from the atomic coordinates. Each chain in the PDB is written as a separate FASTA entry, and the sequence is written on a **single line**.

##  Features

- Converts 3-letter residue codes from PDB to 1-letter FASTA format.
- Skips heteroatoms and water molecules.
- Supports multiple chains.
- Outputs sequences with no line wrapping.
- Handles common modified residues like MSE (selenomethionine).

##  Requirements

- Python 3
- [Biopython](https://biopython.org/)

Install Biopython using pip:

```bash
pip install biopython
