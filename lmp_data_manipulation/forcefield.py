class Forcefield:
    """
    LAMMPS parameterization as it is in .data file

    """

    atom_types = 0
    bond_types = 0
    angle_types = 0
    dihedral_types = 0
    improper_types = 0
    masses = []

    pair_coeffs = []
    bond_coeffs = []
    angle_coeffs = []
    dihedral_coeffs = []
    improper_coeffs = []
    pair_type = 'unknown'
    bonds_type = 'unknown'
    angles_type = 'unknown'
    dihedrals_type = 'unknown'
    impropers_type = 'unknown'
