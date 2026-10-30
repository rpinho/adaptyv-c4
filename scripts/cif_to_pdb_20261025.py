#!/usr/bin/env python3
"""Convert an mmCIF complex to PDB for ProteinMPNN (heavy atoms, chain-preserving)."""
import sys
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO, Select
import warnings
warnings.filterwarnings('ignore')


class Heavy(Select):
    def accept_atom(self, atom):
        return atom.element.strip() != 'H' and atom.get_altloc() in (' ', 'A')


def main(src, dst):
    structure = MMCIFParser(QUIET=True).get_structure('x', str(src))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(dst), Heavy())
    print(dst)


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
