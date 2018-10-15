import os
import pprint
pprint=pprint.PrettyPrinter(indent=4).pprint
import sys


from lmp_utils.data_reader import DataReader


### Some simmple utilities to compare floats #################################

def ae(a, b, epsilon=1e-5):
    """
    ae stands for 'almost equal'
    return true if floats ar ealmost equal

    """
    return abs(a-b) < epsilon*(abs(a) + abs(b))

def aecount(list_of_data, value):
    """
    acount stands for 'approximate-equals count'
    Returns how much time value almost equals to
    element of list_of_data.

    """
    matches = 0
    for v in list_of_data:
        if ae(v, value):
            matches += 1
    return matches

##############################################################################


def main():
    fnames = os.listdir('data')
    #process_single_fname('data/mmt_311_10.data')
    #process_single_fname('mmt_311.data')

    #average_cell('mmt_311.data')

    neighbors_of_al4()


def label(charge):
        if ae(charge, 1.575):
            return 'al'
        elif ae(charge, 1.36):
            return 'mgo'
        elif ae(charge, 2.1):
            return 'si'
        elif ae(charge, -1.05):
            return 'ob'
        elif ae(charge, -0.95):
            return 'oh'
        elif ae(charge, -1.1818):
            return 'obos'
        elif ae(charge, -1.0808):
            return 'ohs'
        elif ae(charge, 0.425):
            return 'ho'
        else:
            return str(charge)#'UNKNOWN'


def process_single_fname(full_fname):
    """
    Trying to extract bead types from atomic coordinates
    using closest neighbors types:
        T       si with 4 neighbors: [ob x4]
        Tns     si with 4 neighbors: [ob x3, obos]
        O       al with 6 neighbors: [ob x4, oh x2]
        Ons     al with 6 neighbors: [ob x2, oh x2, obos x2]
        Onso    al with 6 neighbors: [ob x2, ohs x2, obos x2]
        Os      mgo
        Na      na

    ClayFF charges:
        -0.95    oh
        -1.05    ob
        -1.0808  ohs
        -1.1818  obos    [57, 64, 111, 116]
    """
    lmp_reader = DataReader()
    lmp_reader.read_data(full_fname)

    lx = lmp_reader.xhi - lmp_reader.xlo
    ly = lmp_reader.yhi - lmp_reader.ylo
    lz = lmp_reader.zhi - lmp_reader.zlo

    print(lx, ly, lz)

    atoms = [[] for _ in range(lmp_reader.atoms_number)]

    for atom in lmp_reader.atoms:
        atoms[atom['id'] - 1] = [
            atom['id'],
            atom['charge'],
            atom['x'], atom['y'], atom['z'],
            atom['type']
        ]

    chosen_atoms = []
    for atom in atoms:
        if atom[5] not in [1, 2, 6, 9]: # ao, st, mgo, na
            continue
        chosen_atoms.append(atom)

    hysto = {
        'T':    0,
        'Tns':  0,
        'O':    0,
        'Ons':  0,
        'Onso': 0,
        'Os':   0,
        'Na':   0
    }

    for atom in chosen_atoms:
        neighbors = {} # 6 closest neighbors: { distance: atom }
        for atom_1 in atoms:
            if atom[0] == atom_1[0]:
                continue
            if (not ae(atom_1[1], -0.95) and  # ignore neighboring not-oxygens
                not ae(atom_1[1], -1.05) and
                not ae(atom_1[1], -1.0808) and
                not ae(atom_1[1], -1.1818)):
                    continue

            dx = abs(atom_1[2] - atom[2])
            dy = abs(atom_1[3] - atom[3])
            dz = abs(atom_1[4] - atom[4])
            dx = min(dx, lx-dx)
            dy = min(dy, ly-dy)
            dz = min(dz, lz-dz)

            #if dx > 5 or dy > 5 or dz > 5:
            #    continue
            d = (dx**2 + dy**2 + dz**2)**0.5
            #if len(neighbors.keys()) < 6 or d < max(neighbors.keys()):
            neighbors[d] = atom_1

        #print('neighbors')
        #for k in sorted(neighbors.keys()):
        #    print(k, neighbors[k])

        if ae(atom[1], 1):
            #print('Na', atom)
            hysto['Na'] += 1
        elif ae(atom[1], 1.36):
            #print('Os', atom)
            hysto['Os'] += 1
        elif ae(atom[1], 2.1):  # si
            ds = sorted(neighbors.keys())
            neighbors_charges = [
                neighbors[ds[0]][1], neighbors[ds[1]][1],
                neighbors[ds[2]][1], neighbors[ds[3]][1]
            ]
            if aecount(neighbors_charges, -1.05) == 4:
                #print('T', atom)
                hysto['T'] += 1
            elif (aecount(neighbors_charges, -1.05) == 3 and  # ob
                  aecount(neighbors_charges, -1.1818) == 1):  # obos
                      hysto['Tns'] += 1
                      #print('Tns', atom)
            else:
                print('unknown si')
                print('si', atom)
                for k in sorted(neighbors.keys())[:6]:
                    print(k, neighbors[k])
                sys.exit()
        elif ae(atom[1], 1.575): # al
            ds = sorted(neighbors.keys())
            neighbors_charges = [
                neighbors[ds[0]][1], neighbors[ds[1]][1],
                neighbors[ds[2]][1], neighbors[ds[3]][1],
                neighbors[ds[4]][1], neighbors[ds[5]][1]
            ]
            if (aecount(neighbors_charges, -1.05) == 4 and  # ob
                aecount(neighbors_charges, -0.95) == 2):    # oh
                    hysto['O'] += 1
                    #print('O', atom)
            elif (aecount(neighbors_charges, -1.05) == 2 and  # ob
                  aecount(neighbors_charges, -0.95) == 2 and  # oh
                  aecount(neighbors_charges, -1.1818)):       # obos
                      hysto['Ons'] += 1
                      #print('Ons', atom)
            elif (aecount(neighbors_charges, -1.05) == 2 and    # ob
                  aecount(neighbors_charges, -1.0808) == 2 and  # ohs
                  aecount(neighbors_charges, -1.1818)):         # obos
                      hysto['Onso'] += 1
                      #print('Onso', atom)
            else:
                print('unknown al')
                print('al', atom)
                for k in sorted(neighbors.keys())[:6]:
                    print(k, neighbors[k])
                #sys.exit()
        else:
            print('Completely unknown atom', atom)
            for k in sorted(neighbors.keys())[:6]:
                print(k, neighbors[k])
            sys.exit()

    print('***')
    pprint(hysto)


def average_cell(full_fname):
    lmp_reader = DataReader()
    lmp_reader.read_data(full_fname)

    lx = lmp_reader.xhi - lmp_reader.xlo
    ly = lmp_reader.yhi - lmp_reader.ylo
    lz = lmp_reader.zhi - lmp_reader.zlo

    atoms = lmp_reader.atoms

    for idx, atom in enumerate(atoms):
        if ae(atom['charge'], 1.575):
            neighbors_count = 6
        elif ae(atom['charge'], 2.1):
            neighbors_count = 4
        else:
            continue
        neighbors = {}

        for idx_1, atom_1 in enumerate(atoms):
            if idx == idx_1:
                continue
            if atom_1['type'] == 5:  # H
                continue

            dx = abs(atom_1['x'] - atom['x'])
            dy = abs(atom_1['y'] - atom['y'])
            dz = abs(atom_1['z'] - atom['z'])
            dx = min(dx, lx-dx)
            dy = min(dy, ly-dy)
            dz = min(dz, lz-dz)
            if dx > 5 or dy > 5 or dz > 5:
                continue
            d = (dx**2 + dy**2 + dz**2)**0.5
            neighbors[d] = atom_1

        ks = sorted(neighbors.keys())[:neighbors_count]
        closest_neighbors = [label(neighbors[k]['charge']) for k in ks]

        bead = 'UNKNOWN'
        """
        T       si with 4 neighbors: [ob x4]
        Tns     si with 4 neighbors: [ob x3, obos]
        O       al with 6 neighbors: [ob x4, oh x2]
        Ons     al with 6 neighbors: [ob x2, oh x2, obos x2]
        Onso    al with 6 neighbors: [ob x2, ohs x2, obos x2]
        """
        if closest_neighbors.count('ob') == 4:
            bead = 'T'
        elif (closest_neighbors.count('ob') == 3 and
              closest_neighbors.count('obos') == 1):
                bead = 'Tns'
        elif (closest_neighbors.count('ob') == 4 and
              closest_neighbors.count('oh') == 2):
                bead = 'O'
        elif (closest_neighbors.count('ob') == 2 and
              closest_neighbors.count('oh') == 2 and
              closest_neighbors.count('obos') == 2):
                bead = 'Ons'
        elif (closest_neighbors.count('ob') == 2 and
              closest_neighbors.count('ohs') == 2 and
              closest_neighbors.count('obos') == 2):
                bead = 'Onso'

        print(atom['id'], label(atom['charge']), closest_neighbors, bead)


def neighbors_of_al4():
    lmp_reader = DataReader()
    lmp_reader.read_data('mmt_311.data')

    lx = lmp_reader.xhi - lmp_reader.xlo
    ly = lmp_reader.yhi - lmp_reader.ylo
    lz = lmp_reader.zhi - lmp_reader.zlo

    atoms = lmp_reader.atoms
    al = atoms[3]
    ds = {}

    charge = 0.

    for atom in atoms:
        dx = atom['x'] - al['x']
        dy = atom['y'] - al['y']
        dz = atom['z'] - al['z']
        dx = min(dx, lx-dx)
        dy = min(dy, ly-dy)
        dz = min(dz, lz-dz)
        d = (dx**2 + dy**2 + dz**2)**0.5
        ds[d] = label(atom['charge'])
        charge += atom['charge']

    for idx, k in enumerate(sorted(ds.keys())):
        print(idx, k, ds[k])

    print('charge', charge)


if __name__ == '__main__':
    main()
