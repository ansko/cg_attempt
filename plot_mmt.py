import math
import os
import random

from lmp_new_operations.data_reader import DataReader


def main(fname, axis='xz'):
    """
    Make .svg image of MMT structure

    """

    lmp_reader = DataReader()
    lmp_reader.read_data(fname)

    structure, ff = lmp_reader.structure, lmp_reader.ff

    colors = {
        1: 'black',    # ao
        2: 'yellow',   # st
        3: 'grey',     # ob
        4: 'darkgrey', # oh
        5: 'red',      # ho
        6: 'violet',   # mgo
        7: 'orange',   # obos
        8: 'blue',     # ohs
        9: 'green'     # na
    }
    labels = {
        1: 'ao',
        2: 'st',
        3: 'ob',
        4: 'oh',
        5: 'ho',
        6: 'mgo',
        7: 'obos',
        8: 'ohs',
        9: 'na',
    }

    r = 0.01
    lens = { 'x': 15.5, 'y': 7.5, 'z': 10 }

    left_axis = ''.join(set('xyz') - set(axis))


    MULTIP = 30
    r *= MULTIP
    mag = 4*r
    fractions = 4  # label displacement

    svg_fname = fname[:-5] + '_{0}_'.format(axis) + '.svg'

    with open(svg_fname, 'w') as f:
        min_x = structure.atoms[0][axis[0]]
        max_x = structure.atoms[0][axis[0]]
        min_y = structure.atoms[0][axis[1]]
        max_y = structure.atoms[0][axis[1]]

        for atom in structure.atoms:
            min_x = min(min_x, atom[axis[0]])
            max_x = max(max_x, atom[axis[0]])
            min_y = min(min_x, atom[axis[1]])
            max_y = max(max_x, atom[axis[1]])

        f.write('<?xml version="1.0" encoding="UTF-8" ?>\n')
        f.write('<svg xmlns="http://www.w3.org/2000/svg" version="1.1">\n')

        f.write('<rect x="0" y="0" width="{0}" height="{1}" opacity="0."/>\n'.format(
            (max_x - min_x) * MULTIP, (max_y - min_y) * MULTIP))

        for atom in structure.atoms:
            dr = 0
            if atom['type'] == 8:
                dr = 10

            x = (atom[axis[0]] - min_x) * MULTIP
            y = (atom[axis[1]] - min_y) * MULTIP

            c = colors[atom['type']]
            f.write('<circle cx="{0}" cy="{1}" r="{2}" fill="{3}" />\n'.format(
                 r + x, r + y, r * MULTIP + dr, c))

            phase = 2*math.pi * atom[left_axis] / lens[left_axis]
            dx = MULTIP*r * math.sin(phase)
            dy = MULTIP*r * math.cos(phase)
            atom_str = '_'.join([labels[atom['type']], str(atom['id'])])

            f.write('<text x="{0}" y="{1}">{2}</text>\n'.format(
                x + dx, y + dy, atom_str))

        f.write('</svg>\n')

    return



if __name__ == '__main__':
    #for axis in ['xy', 'xz', 'yz']:
    datafiles = [
        fname for fname in os.listdir('tmp') if fname.endswith('.data')
    ]
    fname = 'mmt_321.data'
    #for fname in datafiles:
    for axis in ['xy', 'xz', 'yz']:
        #fname = 'tmp/{0}'.format(fname)
        main(fname, axis)
