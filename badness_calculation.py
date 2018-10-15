import os
import pprint
pprint=pprint.PrettyPrinter(indent=4).pprint
import shutil
import subprocess

from lmp_new_operations.utils import ae # ae stands for almost equal


not_all_results = 'not_all_results_tmp'
estimates = 'estimates'


defaults = {
    'lx': 5.18, 'ly': 8.98, 'lz': 15.,

    'd11': 6.040216556807043,
    'd12': 6.044135283498653,
    'd22': 6.121346385550226,
}


def call_gnuplot(task=None):
    subprocess.call('./plotter_fits.gnu')
    subprocess.call('./plotter_fits_cell.gnu')
    subprocess.call('./plotter_fits_bonds.gnu')
#    subprocess.call('./plotter_slice.gnu')


def get_badnesses_detailed():
    """
    Get list of entries:
    {
        badness: VALLUE,
        badness_cell: VALLUE,
        badness_bonds: VALUE
        others ....
    }
    Max value of badness is 1, min is 0.

    """

    shutil.copyfile('tmp/not_all_results', './{0}'.format(not_all_results))

    all_badnesses = []
    collected_data = []

    with open(not_all_results) as f:
        lines = f.readlines()[1:-1]
        for line in lines:#range(1, 10):
            if '*' in line:
                continue

            structure, lx, ly, lz, e, d11, d12, d22 = line.split()
            lx = float(lx)
            ly = float(ly)
            lz = float(lz)
            e = float(e)
            d11 = float(d11)
            d12 = float(d12)
            d22 = float(d22)

            badness = 0
            badness_cell = 0
            badness_bonds = 0

            eps = abs(lx - defaults['lx']) / (defaults['lx'] + lx) * 2
            badness += eps**2
            badness_cell += eps**2
            eps = abs(ly - defaults['ly']) / (defaults['ly'] + ly) * 2
            badness += eps**2
            badness_cell += eps**2
            eps = abs(lz - defaults['lz']) / (defaults['lz'] + lz) * 2
            badness += eps**2
            badness_cell += eps**2

            eps = abs(d11 - defaults['d11']) / (d11 + defaults['d11']) * 2
            badness += eps**2
            badness_bonds += eps**2
            eps = abs(d11 - defaults['d12']) / (d11 + defaults['d12']) * 2
            badness += eps**2
            badness_bonds += eps**2
            eps = abs(d11 - defaults['d22']) / (d11 + defaults['d22']) * 2
            badness += eps**2
            badness_bonds += eps**2

            badness /= 6
            badness_cell /= 3
            badness_bonds /= 3
            badness = badness**0.5
            badness_cell = badness_cell**0.5
            badness_bonds = badness_bonds**0.5
            all_badnesses.append(badness)

            with open(estimates, 'a') as fo:
                fo.write('{0} {1}\n'.format(structure, badness))
            collected_data.append({
                'structure': structure,
                'lx': lx, 'ly': ly, 'lz': lz,
                'd11': d11, 'd12': d12, 'd22': d22,
                'epot': e,
                'badness': badness,
                'badness_cell': badness_cell,
                'badness_bonds': badness_bonds,
            })
        #print('best is:', min(all_badnesses),
        #      'in', lines[all_badnesses.index(min(all_badnesses))].split()[0],
        #      'with idx', all_badnesses.index(min(all_badnesses)))
    return collected_data


def explore_all(data):
    idxs = { 'e1': 1, 's1': 2, 'e2': 3, 's2': 4 }
    all_params = ['e1', 'e2', 's1', 's2']
    for param in all_params:
        vals = {}
        vals_cell = {}
        vals_bonds = {}
        for entry in data:
            param_val = float(entry['structure'].split('_')[idxs[param]])
            badness = float(entry['badness'])
            badness_cell = float(entry['badness_cell'])
            badness_bonds = float(entry['badness_bonds'])
            if param_val in vals.keys():
                vals[param_val].append(badness)
                vals_cell[param_val].append(badness_cell)
                vals_bonds[param_val].append(badness_bonds)
            else:
                vals[param_val] = [badness]
                vals_cell[param_val] = [badness_cell]
                vals_bonds[param_val] = [badness_bonds]
        with open('ave_fit_{0}'.format(param), 'w') as f:
            f.write('{0}\tbadness\tbadness_cell\tbadness_bonds\n'.format(param))
            for k in sorted(vals.keys()):
                f.write('{0}\t{1}\t{2}\t{3}\n'.format(k,
                    sum(vals[k]) / len(vals[k]),
                    sum(vals_cell[k]) / len(vals_cell[k]),
                    sum(vals_bonds[k]) / len(vals_bonds[k])))


def explore_minima_sections(data):
    idxs = { 'e1': 1, 's1': 2, 'e2': 3, 's2': 4 }
    minimal_badness = min(entry['badness'] for entry in data)
    min_idx = [i for i, v in enumerate(data)
        if data[i]['badness'] == minimal_badness][0]
    all_params = ['e1', 'e2', 's1', 's2']
    param_min_vals = {}
    for param in all_params:
        param_min_vals[param] = float(
            data[min_idx]['structure'].split('_')[idxs[param]])
    for param in all_params:
        vals = {}
        other_param_vals = {
            k: param_min_vals[k] for k in set(all_params) - set([param])
        }
        for entry in data:
            is_in_slice = True
            for k in other_param_vals.keys():
                entry_param = float(entry['structure'].split('_')[idxs[k]])
                min_param = float(param_min_vals[k])
                if not ae(entry_param, float(min_param)):
                    is_in_slice = False
                    break
            if not is_in_slice:
                continue

            param_val = float(entry['structure'].split('_')[idxs[param]])
            badness = float(entry['badness'])
            if param_val in vals.keys():
                vals[param_val].append(badness)
            else:
                vals[param_val] = [badness]

        with open('slice_{0}'.format(param), 'w') as f:
            f.write('{0}\tbadness\n'.format(param))
            for k in sorted(vals.keys()):
                f.write('{0}\t{1}\n'.format(k, sum(vals[k]) / len(vals[k])))


if __name__ == '__main__':
    #check_all()
    data = get_badnesses_detailed()
    explore_all(data)
    #explore_minima_sections(data)
    call_gnuplot()
    cwd = os.getcwd()
    for fname in os.listdir():
        if fname.startswith('ave_fit_') or fname.startswith('slice_'):
            os.remove(fname)
