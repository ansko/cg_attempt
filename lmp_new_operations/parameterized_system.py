from .structure import Structure
from .forcefield import Forcefield


class ParameterizedSystem:
    def __init__(self, structure, ff):
        self.structure = structure
        self.ff = ff


    def write_json(self, json_fname):
        import json
        json_data = {
            'structure': {
                'atoms_number': self.structure.atoms_number,
                'atoms': self.structure.atoms,
                'atom_types': self.ff.atom_types,
                'bonds_number': self.structure.bonds_number,
                'bonds': self.structure.bonds,
                'bond_types': self.ff.bond_types,
                'angles_number': self.structure.angles_number,
                'angles': self.structure.angles,
                'angle_types': self.ff.angle_types,
                'dihedrals_number': self.structure.dihedrals_number,
                'dihedrals': self.structure.dihedrals,
                'dihedral_types': self.ff.dihedral_types,
                'impropers_number': self.structure.impropers_number,
                'impropers': self.structure.impropers,
                'improper_types': self.ff.improper_types,
            },
            'cell': {
                'xlo': self.structure.xlo, 'xhi': self.structure.xhi,
                'ylo': self.structure.ylo, 'yhi': self.structure.yhi,
                'zlo': self.structure.zlo, 'zhi': self.structure.zhi,
            },
            'masses': self.ff.masses,
            'forcefield': {
                'pair_type': self.ff.pair_type,
                'pair_coeffs': self.ff.pair_coeffs,
                'bonds_type': self.ff.bonds_type,
                'bond_coeffs': self.ff.bond_coeffs,
                'angles_type': self.ff.angles_type,
                'angle_coeffs': self.ff.angle_coeffs,
                'dihedrals_type': self.ff.dihedrals_type,
                'dihedral_coeffs': self.ff.dihedral_coeffs,
                'impropers_type': self.ff.impropers_type,
                'improper_coeffs': self.ff.improper_coeffs,
            }
        }
        json.dump(json_data, open(json_fname, 'w'), indent=4)
