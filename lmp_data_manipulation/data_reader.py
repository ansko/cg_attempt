
# incomplete, look https://lammps.sandia.gov/doc/2001/data_format.html

from .structure import Structure
from .forcefield import Forcefield


class DataReader:
    def __init__(self):
        self.structure = Structure()
        self.ff = Forcefield()


    def read_data(self, data_fname):
        # Get lines, remove empty lines and line breakes
        lines = open(data_fname).readlines()
        lines = [line for line in lines if len(line) > 1]
        lines = [line[:-1] for line in lines if line.endswith('\n')]

        self.LAMMPS_description = lines[0]

        for idx, line in enumerate(lines[1:]):
            if line.endswith('atoms'):
                self.structure .atoms_number = int(line.split()[0])
            elif line.endswith('bonds'):
                self.structure.bonds_number = int(line.split()[0])
            elif line.endswith('angles'):
                self.structure.angles_number = int(line.split()[0])
            elif line.endswith('dihedrals'):
                self.structure.dihedrals_number = int(line.split()[0])
            elif line.endswith('impropers'):
                self.structure.impropers_number = int(line.split()[0])

            elif line.endswith('atom types'):
                self.ff.atom_types = int(line.split()[0])
            elif line.endswith('bond types'):
                self.ff.bond_types = int(line.split()[0])
            elif line.endswith('angle types'):
                self.ff.angle_types = int(line.split()[0])
            elif line.endswith('dihedral types'):
                self.ff.dihedral_types = int(line.split()[0])
            elif line.endswith('improper types'):
                self.ff.improper_types = int(line.split()[0])

            elif line.endswith('xlo xhi'):
                self.structure.xlo = float(line.split()[0])
                self.structure.xhi = float(line.split()[1])
            elif line.endswith('ylo yhi'):
                self.structure.ylo = float(line.split()[0])
                self.structure.yhi = float(line.split()[1])
            elif line.endswith('zlo zhi'):
                self.structure.zlo = float(line.split()[0])
                self.structure.zhi = float(line.split()[1])

            # [ idx+2 : ...] +1 because iterating lines[1:] and
            #                +1 because data starts from the next line
            elif line.startswith('Masses'):
                masses_lines = lines[idx + 2 : idx + 2 + self.ff.atom_types]
            elif line.startswith('Pair Coeffs'):
                try:
                    self.ff.pair_type = line.split()[3]
                except IndexError:
                    pass
                pair_coeffs_lines = lines[idx + 2 : idx + 2 + self.ff.atom_types]
            elif line.startswith('Bond Coeffs'):
                try:
                    self.ff.bonds_type = line.split()[3]
                except IndexError:
                    pass
                bond_coeffs_lines = lines[idx + 2 : idx + 2 + self.ff.bond_types]
            elif line.startswith('Angle Coeffs'):
                try:
                    self.ff.angles_type = line.split()[3]
                except IndexError:
                    pass
                angle_coeffs_lines = lines[idx + 2 : idx + 2 + self.ff.angle_types]
            elif line.startswith('Dihedral Coeffs'):
                try:
                    self.ff.dihedrals_type = line.split()[3]
                except IndexError:
                    pass
                dihedral_coeffs_lines = lines[idx + 2 :
                                              idx + 2 + self.ff.dihedral_types]
            elif line.startswith('Improper Coeffs'):
                try:
                    self.fff.impropers_type = line.split()[3]
                except IndexError:
                    pass
                improper_coeffs_lines = lines[idx + 2 :
                    idx + 2 + self.ff.improper_types]
            elif line.startswith('Atoms'):
                atoms_lines = lines[idx + 2:
                    idx + 2 + self.structure.atoms_number]
            elif line.startswith('Velocities'):
                velocities_lines = lines[idx + 2 :
                    idx + 2 + self.structure.atoms_number]
            elif line.startswith('Bonds'):
                bonds_lines = lines[idx + 2 :
                    idx + 2 + self.structure.bonds_number]
            elif line.startswith('Angles'):
                angles_lines = lines[idx + 2 :
                    idx + 2 + self.structure.angles_number]
            elif line.startswith('Dihedrals'):
                dihedrals_lines = lines[idx + 2 :
                    idx + 2 + self.structure.dihedrals_number]
            elif line.startswith('Impropers'):
                impropers_lines = lines[idx + 2 :
                    idx + 2 + self.structure.impropers_number]

        for line in atoms_lines:
            ls = line.split() # id molecule-tag atom-type q x y z nx ny nz
            self.structure.atoms.append({
                'id': int(ls[0]),
                'molecule_tag': int(ls[1]),
                'type': int(ls[2]),
                'charge': float(ls[3]),
                'x': float(ls[4]), 'y': float(ls[5]), 'z': float(ls[6]),
                'nx': int(ls[7]), 'ny': int(ls[8]), 'nz': int(ls[9]),
            })
        try:  # Bonds are optional
            for line in bonds_lines:
                ls = line.split() # id bond-type atom-1 atom-2
                self.structure.bonds.append({
                    'id': int(ls[0]),
                    'type': int(ls[1]),
                    'atom_1': int(ls[2]),
                    'atom_2': int(ls[3])
                })
        except UnboundLocalError:
            pass
        try:  # Angles are optional
            for line in angles_lines:
                ls = line.split() # 1 angle-type atom-1 atom-2 atom-3
                self.structure.angles.append({
                    'id': int(ls[0]),
                    'type': int(ls[1]),
                    'atom_1': int(ls[2]),
                    'atom_2': int(ls[3]),
                    'atom_3': int(ls[4])
                })
        except UnboundLocalError:
            pass
        try:  # Dihedrals are optional
            for line in dihedrals_lines:
                ls = line.split() # 1 dihedral-type atom-1 atom-2 atom-3 atom-4
                self.structure.dihedrals.append({
                    'id': int(ls[0]),
                    'type': int(ls[1]),
                    'atom_1': int(ls[2]),
                    'atom_2': int(ls[3]),
                    'atom_3': int(ls[4]),
                    'atom_4': int(ls[5])
                })
        except UnboundLocalError:
            pass
        try:  # Impropers are optional
            for line in impropers_lines:
                ls = line.split() # 1 improper-type atom-1 atom-2 atom-3 atom-4
                self.structure.impropers.append({
                    'id': int(ls[0]),
                    'type': int(ls[1]),
                    'atom_1': int(ls[2]),
                    'atom_2': int(ls[3]),
                    'atom_3': int(ls[4]),
                    'atom_4': int(ls[5])
                })
        except UnboundLocalError:
            pass

        for line in masses_lines:
            self.ff.masses.append(float(line.split()[1]))
        try:  # Pair Coeffs are optional (?)
            for line in pair_coeffs_lines:
                if '#' in line:
                    line = line[:line.index('#')]
                ls = line.split()
                self.ff.pair_coeffs.append(
                    [float(coeff) for coeff in ls[1:]])
        except UnboundLocalError:
            pass
        try:  # Bond Coeffs are optional (?)
            for line in bond_coeffs_lines:
                if '#' in line:
                    line = line[:line.index('#')]
                ls = line.split()
                self.ff.bond_coeffs.append(
                    [float(coeff) for coeff in ls[1:]])
        except UnboundLocalError:
            pass
        try:  # Angle Coeffs are optional (?)
            for line in angle_coeffs_lines:
                if '#' in line:
                    line = line[:line.index('#')]
                ls = line.split()
                self.ff.angle_coeffs.append(
                    [float(coeff) for coeff in ls[1:]])
        except UnboundLocalError:
            pass
        try:  # Dihedral Coeffs are optional (?)
            for line in dihedral_coeffs_lines:
                if '#' in line:
                    line = line[:line.index('#')]
                ls = line.split()
                self.ff.dihedral_coeffs.append(
                    [float(coeff) for coeff in ls[1:]])
        except UnboundLocalError:
            pass
        try:  # Improper Coeffs are optional (?)
            for line in improper_coeffs_lines:
                if '#' in line:
                    line = line[:line.index('#')]
                ls = line.split()
                self.ff.improper_coeffs.append(
                    [float(coeff) for coeff in ls[1:]])
        except UnboundLocalError:
            pass


    def print_summary(self):
        print('Structure')
        print('    {0} atoms'.format(self.structure.atoms_number))
        if self.structure.bonds_number:
            print('    {0} bonds'.format(self.structure.bonds_number))
        if self.structure.angles_number:
            print('    {0} angles'.format(self.structure.angles_number))
        if self.structure.dihedrals_number:
            print('    {0} dihedrals'.format(self.structure.dihedrals_number))
        if self.structure.impropers_number:
            print('    {0} impropers'.format(self.structure.impropers_number))

        print('\nIn cell with dimentions:')
        print('    {0} x {1} x {2}'.format(
            self.structure.xhi - self.structure.xlo,
            self.structure.yhi - self.structure.ylo,
            self.structure.zhi - self.structure.zlo))

        print('\nWith masses')
        for idx, mass in enumerate(self.ff.masses):
            print('    {0} {1}'.format(idx+1, mass))

        print('\nParameterized as:')
        if self.ff.pair_coeffs:
            print('    Pair Coeffs: {0}'.format(self.ff.pair_type))
            for idx, coeffs in enumerate(self.ff.pair_coeffs):
                print(' '*7, idx+1, ' '.join(str(coeff) for coeff in coeffs))
        if self.ff.bond_coeffs:
            print('    Bond Coeffs: {0}'.format(self.ff.bonds_type))
            for idx, coeffs in enumerate(self.ff.bond_coeffs):
                print(' '*7, idx+1, ' '.join(str(coeff) for coeff in coeffs))
        if self.ff.angle_coeffs:
            print('    Angle Coeffs: {0}'.format(self.ff.angles_type))
            for idx, coeffs in enumerate(self.ff.angle_coeffs):
                print(' '*7, idx+1, ' '.join(str(coeff) for coeff in coeffs))
        if self.ff.dihedral_coeffs:
            print('    Dihedral Coeffs: {0}'.format(self.ff.dihedrals_type))
            for idx, coeffs in enumerate(self.ff.dihedral_coeffs):
                print(' '*7, idx+1, ' '.join(str(coeff) for coeff in coeffs))
        if self.ff.improper_coeffs:
            print('    Improper Coeffs: {0}'.format(self.ff.impropers_type))
            for idx, coeffs in enumerate(self.ff.improper_coeffs):
                print(' '*7, idx+1, ' '.join(str(coeff) for coeff in coeffs))
