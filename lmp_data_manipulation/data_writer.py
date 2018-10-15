class DataWriter:
    def __init__(self, structure, ff=None):
        self.structure = structure
        self.ff = ff

    def write_data(self, fname):
        with open(fname, 'w') as f:
            f.write('LAMMPS data file by DataWriter.py\n\n')

            f.write('{0} atoms\n'.format(self.structure.atoms_number))
            if self.structure.bonds_number:
                f.write('{0} bonds\n'.format(self.structure.bonds_number))
                if self.structure.angles_number:
                    f.write('{0} angles\n'.format(self.structure.angles_number))
                    if self.structure.dihedrals_number:
                        f.write('{0} dihedrals\n'.format(
                            self.structure.dihedrals_number))
                        if self.structure.impropers_number:
                            f.write('{0} impropers\n'.format(
                                self.structure.impropers_number))

            f.write('\n')

            if self.ff is not None:
                f.write('{0} atom types\n'.format(self.ff.atom_types))
                if self.ff.bond_types:
                    f.write('{0} bond typess\n'.format(self.ff.bond_types))
                    if self.ff.angle_types:
                        f.write('{0} angle typess\n'.format(self.ff.angle_types))
                        if self.ff.dihedral_types:
                            f.write('{0} dihedral types\n'.format(
                                self.ff.dihedral_types))
                            if self.ff.improper_types:
                                f.write('{0} improper types\n'.format(
                                    self.ff.dihedral_types))

            f.write('\n{0} {1} xlo xhi\n'.format(
                self.structure.xlo, self.structure.xhi))
            f.write('{0} {1} ylo yhi\n'.format(
                self.structure.ylo, self.structure.yhi))
            f.write('{0} {1} zlo zhi\n'.format(
                self.structure.zlo, self.structure.zhi))

            f.write('\nMasses\n\n')

            for idx, mass in enumerate(self.ff.masses):
                f.write('{0} {1}\n'.format(idx+1, mass))

            f.write('\nPair Coeffs\n\n')

            for idx, pair_coeff in enumerate(self.ff.pair_coeffs):
                f.write('{0} {1}\n'.format(idx+1,
                    ' '.join([str(coeff) for coeff in pair_coeff])))

            if self.ff.bond_coeffs:
                f.write('\nBond Coeffs\n\n')
                for idx, bond_coeff in enumerate(self.ff.bond_coeffs):
                    f.write('{0} {1}\n'.format(idx+1,
                        ' '.join([str(coeff) for coeff in bond_coeff])))

            if self.ff.angle_coeffs:
                f.write('\nAngle Coeffs\n\n')
                for idx, angle_coeff in enumerate(self.ff.angle_coeffs):
                    f.write('{0} {1}\n'.format(idx+1,
                        ' '.join([str(coeff) for coeff in angle_coeff])))

            if self.ff.dihedral_coeffs:
                f.write('\nDihedral Coeffs\n\n')
                for idx, dihedral_coeff in enumerate(self.ff.dihedral_coeffs):
                    f.write('{0} {1}\n'.format(idx+1,
                        ' '.join([str(coeff) for coeff in dihedral_coeff])))

            if self.ff.improper_coeffs:
                f.write('\nImproper Coeffs\n\n')
                for idx, improper_coeff in enumerate(self.ff.improper_coeffs):
                    f.write('{0} {1}\n'.format(idx+1,
                        ' '.join([str(coeff) for coeff in improper_coeff])))

            f.write('\nAtoms\n\n')
            for atom in self.structure.atoms:
                f.write('{0}\n'.format(' '.join([
                    str(atom['id']), str(atom['molecule_tag']),
                    str(atom['type']), str(atom['charge']),
                    str(atom['x']), str(atom['y']), str(atom['z']),
                    str(atom['nx']), str(atom['ny']), str(atom['nz'])])))
            if self.structure.bonds:
                f.write('\nBonds\n\n')
                for bond in self.structure.bonds:
                    f.write('{0}\n'.format(' '.join([
                        str(bond['id']), str(bond['type']),
                        str(bond['atom_1']), str(bond['atom_2'])])))
            if self.structure.angles:
                f.write('\nAngles\n\n')
                for angle in self.structure.angles:
                    f.write('{0}\n'.format(' '.join([
                        str(angle['id']), str(angle['type']),
                        str(angle['atom_1']), str(angle['atom_2']),
                        str(angle['atom_3'])])))
            if self.structure.dihedrals:
                f.write('\nDihedrals\n\n')
                for dihedral in self.structure.dihedrals:
                    f.write('{0}\n'.format(' '.join([
                        str(dihedral['id']), str(dihedral['type']),
                        str(dihedral['atom_1']), str(dihedral['atom_2']),
                        str(dihedral['atom_3']), str(dihedral['atom_4'])])))
            if self.structure.impropers:
                f.write('\nImpropers\n\n')
                for improper in self.structure.impropers:
                    f.write('{0}\n'.format(' '.join([
                        str(improper['id']), str(improper['type']),
                        str(improper['atom_1']), str(improper['atom_2']),
                        str(improper['atom_3']), str(improper['atom_4'])])))

#            #TODO
#            if self.structure.velocities:
#                for vel in self.structure.velocities:
#                    pass
