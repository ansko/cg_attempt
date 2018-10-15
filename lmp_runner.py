data_files_dir = 'data_adjusting'
tmp_dir = 'tmp'

all_results = 'tmp/all_results'
not_all_results = 'tmp/not_all_results'


def make_lmp_in(fname, lmp_fname):
    content = """
units           real
atom_style      full
pair_style      lj/cut/coul/long 10.0
bond_style      harmonic
angle_style     harmonic
dihedral_style  harmonic
improper_style  cvff
kspace_style    ewald             0.00001
neighbor        1.0               nsq
neigh_modify    once no every 1 delay 0 check yes page 100000000 one 10000000
thermo          1
timestep        0.1
read_data       data_adjusting/{0}
thermo_style    custom step lx ly lz pe
#dump            d0 all image 100 tmp/lmp_img.*.jpg type type
fix             ensemble all npt iso 0.0 0.0 10000 temp 298.15 298.15 100
run             100000
write_data      tmp/{0}.100000.data
    """
    with open(lmp_fname, 'w') as f:
        f.write(content.format(fname))


def get_params_log(out_fname):
    f = open(out_fname)
    lines = f.readlines()[-100038:-38]
    ave_lxs = []
    ave_lys = []
    ave_lzs = []
    ave_Eps = []

    for idx in range(5):
        idx_lines = lines[idx*20000:idx*20000+20000-1]
        idx_lines = [line.split() for line in idx_lines]
        ave_lxs.append(
            sum([float(line[1]) for line in idx_lines]) / len(idx_lines))
        ave_lys.append(
            sum([float(line[2]) for line in idx_lines]) / len(idx_lines))
        ave_lzs.append(
            sum([float(line[3]) for line in idx_lines]) / len(idx_lines))
        ave_Eps.append(
            sum([float(line[4]) for line in idx_lines]) / len(idx_lines))


    #with open(all_results, 'a') as f:
    #    f.write('\t'.join(str(num) for num in
    #        [*ave_lxs, *ave_lys, *ave_lzs, *ave_Eps]
    #    ))
    with open(not_all_results, 'a') as f:
        f.write('\t'.join(str(num) for num in
            [ave_lxs[-1], ave_lys[-1], ave_lzs[-1], ave_Eps[-1]]
        ))


def get_params_data(data_fname):
    f = open(data_fname)
    lines = f.readlines()
    xlo = float(lines[5].split()[0])
    xhi = float(lines[5].split()[1])
    ylo = float(lines[6].split()[0])
    yhi = float(lines[6].split()[1])
    zlo = float(lines[7].split()[0])
    zhi = float(lines[7].split()[1])
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo
    distances = {
        1: { 1: [], 2: [] },
        2: { 1: [], 2: [] }
    }
    cutoff = 10000
    atoms = lines[21:33]

    if not (lx < 100 and ly < 100 and lz < 100):
        with open(not_all_results, 'a') as f:
            f.write('*\t*\t*\n')
        return

    for atom_1 in atoms:
        for atom_2 in atoms:
            if atom_1 is atom_2:
                continue
            dx = abs(float(atom_1.split()[4]) - float(atom_2.split()[4]))
            dy = abs(float(atom_1.split()[5]) - float(atom_2.split()[5]))
            dz = abs(float(atom_1.split()[6]) - float(atom_2.split()[6]))
            dx = min(dx, lx-dx)
            dy = min(dy, ly-dy)
            dz = min(dz, lz-dz)
            if dx > cutoff or dy > cutoff or dz > cutoff:
                continue
            d = (dx**2 + dy**2 + dz**2)**0.5
            distances[int(atom_1.split()[2])][int(atom_2.split()[2])].append(d)

    for idx_1 in [1, 2]:
        for idx_2 in [1, 2]:
            distances[idx_1][idx_2].sort()
            if len(distances[idx_1][idx_2]) == 0:
                print(distances, data_fname)
            if len(distances[idx_1][idx_2]) == 1:
                distances[idx_1][idx_2] = distances[idx_1][idx_2][0]
                continue
            d = distances[idx_1][idx_2][0]
            ave_d = 0
            count = 0
            for d2 in distances[idx_1][idx_2]:
                if abs(d2 - d) > 0.01 * (d2 + d):
                    break
                ave_d += d2
                count += 1
            distances[idx_1][idx_2] = ave_d / count
    #with open(all_results, 'a') as f:
    #    f.write('\t'.join(str(num) for num in 
    #        ['', distances[1][1], distances[1][2], distances[2][2]]
    #    ))
    #    f.write('\n')
    with open(not_all_results, 'a') as f:
        f.write('\t'.join(str(num) for num in
            ['', distances[1][1], distances[1][2], distances[2][2], '\n']
        ))


def main():
    import os
    import subprocess
    import time

    fnames = os.listdir(data_files_dir)

    print('start', len(fnames))

    with open(all_results, 'w') as f:
        f.write('\t'.join(['structure',
            'lx0', 'lx1', 'lx2', 'lx3', 'lx4',
            'ly0', 'ly1', 'ly2', 'ly3', 'ly4',
            'lz0', 'lz1', 'lz2', 'lz3', 'lz4',
            'Ep0', 'Ep1', 'Ep2', 'Ep3', 'Epr4',
            'd11', 'd12', 'd22\n'
        ]))
    with open(not_all_results, 'w') as f:
        f.write('\t'.join([
            'structure',
            'lx4', 'ly4', 'lz4', 'Ep4',
            'd11', 'd12', 'd22\n'
        ]))

    for idx in range(len(fnames)):
        it_start = time.time()
        print('{0} of {1}'.format(idx+1, len(fnames)), end=' ')

        with open(all_results, 'a') as f:
            f.write('{0}\t'.format(fnames[idx][:-5]))

        with open(not_all_results, 'a') as f:
            f.write('{0}\t'.format(fnames[idx][:-5]))

        lmp_fname = 'tmp/' + fnames[idx][:-5]
        lmp_out_fname = 'tmp/' + fnames[idx][:-5] + '_out'
        make_lmp_in(fnames[idx], lmp_fname)
        code = subprocess.call(['lammps-daily', '-in', lmp_fname],
            stdout=open(lmp_out_fname, 'w'))
        if code != 0:
            with open(not_all_results, 'a') as f:
                f.write('* '*7 + '\n')
            continue
        get_params_log('tmp/' + fnames[idx][:-5] + '_out')
        get_params_data('tmp/{0}.100000.data'.format(fnames[idx]))

        print('took', time.time() - it_start)


if __name__ == '__main__':
    # clean previous results
    import subprocess
    subprocess.call(['rm', '-rf', 'tmp'])
    subprocess.call(['mkdir', 'tmp'])

    main()
