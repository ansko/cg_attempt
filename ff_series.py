import pprint
pprint=pprint.PrettyPrinter(indent=4).pprint

from lmp_new_operations.structure import Structure
from lmp_new_operations.forcefield import Forcefield
from lmp_new_operations.data_reader import DataReader
from lmp_new_operations.data_writer import DataWriter
from lmp_new_operations.parameterized_system import ParameterizedSystem

from lmp_new_operations.get_info import neighbors_types
from lmp_new_operations.process_info import process_neighbors_oxygens

from lmp_new_operations.get_info import neighbors_types


def make_ff_series(input_fname='pyr_111.data', data_folder='data_adjusting'):
    """
    Read file `input_fname`, coarse-grain it and produce a lot of
    .data files with coordinates of coarse-grained structure
    but different forcefield parameters. Filenames are:
        `data_folder/cg_{0}_{1}_{2}_{3}.data` where
            (0, 1, 2, 3) = (eps_1, sig_1, eps_2, sig_2)
    """

    lmp_reader = DataReader()
    lmp_reader.read_data(input_fname)

    structure, ff = lmp_reader.structure, lmp_reader.ff

    distances = neighbors_types(structure)
    neighbors_counts = {
        1: 6, # ao   6 bonds
        2: 4, # st   4 bonds
        3: 3, # ob   3 bonds
        4: 3, # oh   3 bonds
        5: 1, # ho   1 bond
        6: 6, # mgo  6 bonds
        7: 3, # obos 3 bonds
        8: 3, # ohs  3 bonds
        9: 6, # Na   no bonds
    }
    labels = {
        1: 'ao  ', 2: 'st  ', 3: 'ob  ', 4: 'oh  ', 5: 'ho  ', 6: 'mgo ',
        7: 'obos', 8: 'ohs ', 9: 'na  '
    }
    oxygen_types = { 3: 'ob', 4: 'oh', 7: 'obos', 8: 'ohs' }

    def closest_oxygens(structure, atom_idx, oxygen_types, cutoff=5, epsilon=0.03):
        lx = structure.xhi - structure.xlo
        ly = structure.yhi - structure.ylo
        lz = structure.zhi - structure.zlo
        #for idx_1, atom_1 in enumerate(structure.atoms):
        idx_1 = atom_idx
        atom_1 = structure.atoms[atom_idx]
        if True:
            closest_oxygens = {}  # { distance: oxygen_type }
            for idx_2, atom_2 in enumerate(structure.atoms):
                if atom_2['type'] not in oxygen_types:
                    continue
                dx = abs(atom_2['x'] - atom_1['x'])
                dy = abs(atom_2['y'] - atom_1['y'])
                dz = abs(atom_2['z'] - atom_1['z'])
                dx = min(dx, lx - dx)
                dy = min(dy, ly - dy)
                dz = min(dz, lz - dz)
                if dx > cutoff or dy > cutoff or dz > cutoff:
                    continue
                d = (dx**2 + dy**2 + dz**2)**0.5
                if d > cutoff:
                    continue

                closest_oxygens[d] = [atom_2['type'], atom_2['id']]

        return closest_oxygens

    """
        T       si with 4 neighbors: [ob x4]
        Tns     si with 4 neighbors: [ob x3, obos]
        O       al with 6 neighbors: [ob x4, oh x2]
        Ons     al with 6 neighbors: [ob x2, oh x2, obos x2]
        Onso    al with 6 neighbors: [ob x2, ohs x2, obos x2] ---
        Onso    al with 6 neighbors: [ob x4, ohs x2]          +++
        Os      mgo
        Na      na
    """
    bead_types = ['T', 'Tns', 'O', 'Ons', 'Onso', 'Os', 'Na']
    bead_charges = {
        'T': 0.175,
        'Tns': 0.131066667,
        'O': -0.35,
        'Ons': 0.437866667,
        'Onso': 0.568666667,
        'Os': 0.871533333,
        'Na': 1.0,
    }
    bead_masses = {
        'T': 57.4177,
        'Tns': 57.4177,
        'O': 65.3214,
        'Ons': 65.3214,
        'Onso': 65.3214,
        'Os': 62.6449,
        'Na': 22.99
    }
    hysto = {
        'T': 0, 'Tns': 0, 'O': 0, 'Ons': 0, 'Onso': 0, 'Os': 0, 'Na': 0,
        'unknown': 0
    }

    cg_structure = Structure()
    cg_atoms = []
    cg_id = 1
    for atom_idx, atom in enumerate(structure.atoms):
        if atom['type'] not in [1, 2,  6, 9]:
            continue

        co = closest_oxygens(structure, atom_idx, oxygen_types)

        min_d = sorted(co.keys())
        min_d = min_d[:neighbors_counts[structure.atoms[atom_idx]['type']]]
        co = {k: co[k] for k in min_d}

        me = labels[structure.atoms[atom_idx]['type']]

        names = [labels[co[d][0]] for d in sorted(co.keys())]

        bead_type = 'unknown'

        if me == 'st  ' and names.count('ob  ') == 4:
            bead_type = 'T'
        elif (me == 'st  ' and
              names.count('ob  ') == 3 and
              names.count('obos') == 1):
                  bead_type = 'Tns'
        elif (me == 'ao  ' and
              names.count('ob  ') == 4 and
              names.count('oh  ') == 2):
                  bead_type = 'O'
        elif (me == 'ao  ' and
              names.count('ob  ') == 2 and
              names.count('oh  ') == 2 and
              names.count('obos') == 2):
                  bead_type = 'Ons'
        elif (me == 'ao  ' and
              names.count('ob  ') == 4 and
              names.count('ohs ') == 2):
                  bead_type = 'Onso'
        elif me == 'mgo ':
            bead_type = 'Os'
        elif me == 'na  ':
            bead_type = 'Na'

        hysto[bead_type] += 1

        cg_atoms.append({
            'id': cg_id,
            'molecule_tag': 1,
            'type': bead_types.index(bead_type) + 1,
            'charge': bead_charges[bead_type],
            'x': atom['x'],
            'y': atom['y'],
            'z': atom['z'],
            'nx': 0, 'ny': 0, 'nz': 0
        })
        cg_id += 1

    clean_bead_types = [t for t in bead_types if hysto[t] != 0]
    clean_bead_charges = {
        k: bead_charges[k] for k in bead_types if hysto[k] != 0
    }
    clean_bead_masses = {
        k: bead_masses[k] for k in bead_types if hysto[k] != 0
    }

    new_cg_atoms = []
    for atom in cg_atoms:
        new_type = clean_bead_types.index(bead_types[atom['type'] - 1]) + 1
        new_cg_atoms.append({
            'id': atom['id'],
            'molecule_tag': 1,
            'type': new_type,
            'charge': clean_bead_charges[clean_bead_types[new_type-1]],
            'x': atom['x'],
            'y': atom['y'],
            'z': atom['z'],
            'nx': 0, 'ny': 0, 'nz': 0
        })

    cg_structure.atoms = new_cg_atoms
    cg_structure.atoms_number = len(new_cg_atoms)
    cg_structure.bonds = []
    cg_structure.bonds_number = 0 
    cg_structure.angles = []
    cg_structure.angles_number = 0
    cg_structure.dihedrals = []
    cg_structure.dihedrals_number = 0
    cg_structure.impropers = []
    cg_structure.impropers_number = 0

    cg_structure.xlo = structure.xlo
    cg_structure.xhi = structure.xhi
    cg_structure.ylo = structure.ylo
    cg_structure.yhi = structure.yhi
    cg_structure.zlo = structure.zlo
    cg_structure.zhi = structure.zhi

    cg_ff = Forcefield()
    cg_ff.atom_types = len(clean_bead_types)
    cg_ff.bond_types = 0
    cg_ff.angle_types = 0
    cg_ff.dihedral_types = 0
    cg_ff.improper_types = 0
    cg_ff.masses = [bead_masses[k] for k in clean_bead_types]
    cg_ff.bond_coeffs = []
    cg_ff.angle_coeffs = []
    cg_ff.dihedral_coeffs = []
    cg_ff.improper_coeffs = []
    cg_ff.bonds_type = 'unknown'
    cg_ff.angles_type = 'unknown'
    cg_ff.dihedrals_type = 'unknown'
    cg_ff.impropers_type = 'unknown'

    # Some unknown values #
    #cg_ff.pair_coeffs = [
        #[0., 0, 2.6309744427952944], # si
        #[0, 2.6309744427952944],
        #[0., 2.662315859650639], # al
        #[0, 2.662315859650639],
        #[0, 2.662315859650639],
        #[0, 2.662315859650639],
        #[0.1300999871, 2.3500126639], # Na from ClayFF
    #]
    # rmin = 2.9531689621096175 in all T...-s (si)
    # rmin = 2.988348513069986  in all O...-s (ao)

    """
    There are just 4 parameters to adjust:
        eps_T, eps_0, sig_T, sig_0

    In ClayFF these are approximate parameters values:
        species    energy    r_min
        st         1.8e-6    3.7
        ao         1.3e-6    4.8
        mgo        9.0e-7    5.9
        na         0.13      2.6

    For pure MMT energies are close; r-s are also close.
    Potential energy values are:
        311 MMT -29100 Kcal/mole
        311 Pyr -28900  Kcal/mole
        111 Pyr -9600  Kcal/mole

        species    energy    r_min
        T          1.8e-6    3.7
        O          1.3e-6    4.8

    Rough approximation:
        eps, sig of T, 0 are equal and equal to
        [10,  4]    ->  PotEng = -465
        [25,  4]    ->           -1100
        [100, 4]    ->           -7100
        [150, 4]    ->           -10600

    """

    eps_s1 = [0.01, 0.1, 1, 10, 100, 1000]
    eps_s2 = [0.01, 0.1, 1, 10, 100, 1000]
    sig_s1 = [6]#[r for r in range(1, 11, 1)]
    sig_s2 = [6]#[r for r in range(1, 11, 1)]



    cg_ff.pair_type = 'lj/cut/coul/long'

    all_N = len(eps_s1)*len(eps_s2) * len(sig_s1)*len(sig_s2)
    idx = 0
    for eps_T in eps_s1:
        for eps_O in eps_s2:
            for sig_T in sig_s1:
                for sig_O in sig_s2:
                    print(idx+1,  all_N)
                    idx += 1

                    cg_ff.pair_coeffs = [
                        [eps_T, sig_T],
                        [eps_O, sig_O]
                    ]
                    out_data_fname = '{0}/cg_{1}_{2}_{3}_{4}.data'.format(
                        data_folder,
                        cg_ff.pair_coeffs[0][0], cg_ff.pair_coeffs[0][1],
                        cg_ff.pair_coeffs[1][0], cg_ff.pair_coeffs[1][1])

                    lmp_writer = DataWriter(cg_structure, cg_ff)
                    lmp_writer.write_data(out_data_fname)

    print('done')


if __name__ == '__main__':
    # clean previous results
    import subprocess
    subprocess.call(['rm', '-rf', 'data_adjusting'])
    subprocess.call(['mkdir', 'data_adjusting'])

    # produce new data
    make_ff_series()
