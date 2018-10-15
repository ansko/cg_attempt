import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint

from lmp_new_operations.data_reader import DataReader
from lmp_new_operations.utils import f_list


def ave(some_list):
    return sum(some_list) / len(some_list)

def get_distances(data_fname, interesting_types):
    r = DataReader()
    r.read_data(data_fname)

    lx = r.structure.xhi - r.structure.xlo
    ly = r.structure.yhi - r.structure.ylo
    lz = r.structure.zhi - r.structure.zlo
    distances = {
        k1 : {
            k2: [] for k2 in interesting_types
        } for k1 in interesting_types
    }
    for idx_1, atom_1 in enumerate(r.structure.atoms):
        for idx_2, atom_2 in enumerate(r.structure.atoms):
            if idx_1 == idx_2:
                continue
            if (atom_1['type'] not in interesting_types or
                atom_2['type'] not in interesting_types):
                    continue
            dx = abs(atom_2['x'] - atom_1['x'])
            dy = abs(atom_2['y'] - atom_1['y'])
            dz = abs(atom_2['z'] - atom_1['z'])
            dx = min(dx, lx-dx)
            dy = min(dy, ly-dy)
            dz = min(dz, lz-dz)
            d = (dx**2 + dy**2 + dz**2)**0.5
            distances[atom_1['type']][atom_2['type']].append(d)

    distances = {
        k1: {
            k2: distances[k1][k2] for k2 in distances[k1].keys()
                if distances[k1][k2]
        } for k1 in distances.keys()
    }

    for k1 in interesting_types:
        for k2 in interesting_types:
            if len(distances[k1][k2]) == 1:
                distances[k1][k2] = distances[k1][k2][0]
            d = distances[k1][k2][0]
            count = 0
            ave = d
            for d2 in distances[k1][k2]:
                if abs(d2 - d) < 0.01 * abs(d2 + d):
                    ave += d2
                    count += 1
                else:
                    break
            distances[k1][k2] = ave / count
    pprint(distances)


if __name__ == '__main__':
    get_distances('pyr_111.100000.data', [1, 2])
