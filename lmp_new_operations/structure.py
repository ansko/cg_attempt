class Structure:
    def __init__(self):
        self.atoms = []
        self.atoms_number = 0
        self.bonds = []
        self.bonds_number = 0 
        self.angles = []
        self.angles_number = 0
        self.dihedrals = []
        self.dihedrals_number = 0
        self.impropers = []
        self.impropers_number = 0

        self.xlo = 0
        self.xhi = 0
        self.ylo = 0
        self.yhi = 0
        self.zlo = 0
        self.zhi = 0


    def replicate(self, **kwargs):
        x_times = 1 if not 'x' in kwargs.keys() else kwargs['x']
        y_times = 1 if not 'y' in kwargs.keys() else kwargs['y']
        z_times = 1 if not 'z' in kwargs.keys() else kwargs['z']
        MULT = x_times * y_times * z_times
        lx = self.xhi - self.xlo
        ly = self.yhi - self.ylo
        lz = self.zhi - self.zlo

        max_molecule_tag = max([atom['molecule_tag'] for atom in self.atoms])

        for idx_x in range(x_times):
            for idx_y in range(y_times):
                for idx_z in range(z_times):
                    if idx_x == idx_y == idx_z == 0:
                        continue
                    sc_id = idx_z*(x_times*y_times) + idx_y*x_times + idx_x

                    AN = len(self.atoms)

                    for idx in range(AN):
                        at = self.atoms[idx]
                        new_atom = {
                            'id': sc_id * AN + at['id'],
                            'molecule_tag': (max_molecule_tag * sc_id +
                                             at['molecule_tag']),
                            'type': at['type'],
                            'charge': at['charge'],
                            'x': at['x'] + lx * idx_x,
                            'y': at['y'] + ly * idx_y,
                            'z': at['z'] + lz * idx_z,
                            'nx': at['nx'],
                            'ny': at['ny'],
                            'nz': at['nz']
                        }
                        self.atoms.append(new_atom)

                    for idx in range(len(self.bonds)):
                        bond = self.bonds[idx]
                        new_bond = {
                            'id': sc_id * len(self.bonds),
                            'type': bond['type'],
                            'atom_1': bond['atom_1'] + sc_id * AN,
                            'atom_2': bond['atom_2'] + sc_id * AN,
                        }
                        self.bonds.append(new_bond)

                    for idx in range(len(self.angles)):
                        angle = self.angles[idx]
                        new_angle = {
                            'id': sc_id * len(self.angles),
                            'type': angle['type'],
                            'atom_1': angle['atom_1'] + sc_id * AN,
                            'atom_2': angle['atom_2'] + sc_id * AN,
                            'atom_3': angle['atom_3'] + sc_id * AN,
                        }
                        self.angles.append(new_bond)

                    for idx in range(len(self.dihedrals)):
                        dihedral = self.dihedrals[idx]
                        new_dihedral = {
                            'id': sc_id * len(self.dihedrals),
                            'type': dihedral['type'],
                            'atom_1': dihedral['atom_1'] + sc_id * AN,
                            'atom_2': dihedral['atom_2'] + sc_id * AN,
                            'atom_3': dihedral['atom_3'] + sc_id * AN,
                            'atom_4': dihedral['atom_4'] + sc_id * AN,
                        }
                        self.dihedrals.append(new_bond)

                    for idx in range(len(self.impropers)):
                        improper = self.bonds[idx]
                        new_improper = {
                            'id': sc_id * len(self.impropers),
                            'type': improper['type'],
                            'atom_1': improper['atom_1'] + sc_id * AN,
                            'atom_2': improper['atom_2'] + sc_id * AN,
                            'atom_3': improper['atom_3'] + sc_id * AN,
                            'atom_4': improper['atom_4'] + sc_id * AN,
                        }
                        self.impropers.append(new_improper)

        self.atoms_number *= MULT
        self.bonds_number *= MULT
        self.angles_number *= MULT
        self.dihedrals_number *= MULT
        self.impropers_number *= MULT

        self.xlo *= x_times
        self.xhi *= x_times
        self.ylo *= y_times
        self.yhi *= y_times
        self.zlo *= z_times
        self.zhi *= z_times
