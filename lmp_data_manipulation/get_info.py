def neighbors_types(structure,
        cutoff=10):
    """
    For every type of atom get its closest neighbors.
    Input parameters:
        cutoff=10 - maximal distance to neighbor

    """

    AN = structure.atoms_number
    atoms = structure.atoms
    lx = structure.xhi - structure.xlo
    ly = structure.yhi - structure.ylo
    lz = structure.zhi - structure.zlo

    # simplify atoms for better performance
    atoms = [
        [atom['id'], atom['type'], atom['x'], atom['y'], atom['z']]
            for atom in atoms
    ]

    results = {}

    for idx_1, atom_1 in enumerate(atoms):
        for idx_2, atom_2 in enumerate(atoms):
            if idx_1 == idx_2:
                continue
            dx = abs(atom_2[2] - atom_1[2])
            dy = abs(atom_2[3] - atom_1[3])
            dz = abs(atom_2[4] - atom_1[4])
            dx = min(lx - dx, dx)
            dy = min(ly - dy, dy)
            dz = min(lz - dz, dz)
            if dx > cutoff or dy > cutoff or dz > cutoff:
                continue
            d = (dx**2 + dy**2 + dz**2)**0.5
            if d > cutoff:
                continue

            if atom_1[1] in results.keys():
                if atom_2[1] in results[atom_1[1]].keys():
                    results[atom_1[1]][atom_2[1]].append(d)
                else:
                    results[atom_1[1]][atom_2[1]] = [d]
            else:
                results[atom_1[1]] = { atom_2[1]: [d] }

    return results
