import pprint
pprint=pprint.PrettyPrinter(indent=4).pprint

import sys


def process_neighbors_types(distances, neighbors_counts, epsilon=0.03):
    """
    Process the dict distances:
    { type_1: { type_2: [distance, ... ], ... }, ... }
    by averaging close ones (differing less than epsilon*distnace).
    Input parameters:
        neighbors_count - number of neighbors by atom type
    """
    ave_results = {}

    for k_1 in distances.keys():
        for k_2 in distances[k_1].keys():
            values = distances[k_1][k_2]

            groups = {} # not the same values

            for value in values:
                flag_appended = False
                for ave in groups.keys():
                    if abs(ave-value) < epsilon * (ave+value):
                        flag_appended = True
                        groups[ave].append(value)
                if not flag_appended:
                    groups[value] = [value]

            for ave in groups.keys():
                groups[ave] = sum(groups[ave]) / len(groups[ave])

            num = neighbors_counts[k_1]
            #for ave in sorted(groups.keys())[num-1:]:
            #    del groups[ave]
            groups = { k: groups[k] for k in sorted(groups.keys())[:num] }
            print(len(groups))

            if k_1 in ave_results.keys():
                ave_results[k_1][k_2] = groups
            else:
                ave_results[k_1] = { k_2: groups }

    return ave_results


def process_neighbors_oxygens(distances, neighbors_counts, epsilon=0.03):
    """
    Process the dict distances:
    { type_1: { type_2: [distance, ... ], ... }, ... }
    by averaging close ones (differing less than epsilon*distnace).
    Input parameters:
        neighbors_count - number of neighbors by atom type
    """

    return 0
