import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

int_filename = sys.argv[1]

LINKAGE = "single"
THRESHOLD = 0.7071  # np.sqrt(2)/2


def read_int_matrix(filename):
    """Read the interface matrix"""
    matrices = {}
    ligands = {}
    lines = open(filename, "r").read().split("\n")
    for ln in lines:
        if not ln:
            logging.warning("empty line, are we at the end of file?")
            continue
        key_values = ln.split(",")
        # extracting receptor
        rec_id = key_values[0].split()[0].lstrip('"')
        key_values[0] = key_values[0].split("{")[1]
        print("rec_id = ", rec_id)
        # check empty
        if len(key_values) == 1:
            logging.warning("Empty interface matrix. Discarding...")
            continue
        len_kv = len(key_values)
        print(f"{len_kv} string elements")
        # creating uncondensed distance matrix
        nstructs = int(np.sqrt(len_kv))
        npairs = nstructs * (nstructs - 1) // 2
        print(f"{nstructs} structures, {npairs} npairs")
        matrices[rec_id] = np.zeros(npairs)
        ligands[rec_id] = []
        entries_list = []
        knt = 0
        for el in key_values:
            # print(el)
            key = el.split(":")[0].lstrip(" '").rstrip("'")
            value = el.split(":")[1].rstrip("").lstrip().rstrip("}\n")
            # print(key)
            first_lig = key.split("_")[0]
            second_lig = key.split("_")[1]
            reverse_lig = second_lig + "_" + first_lig
            if (first_lig == second_lig) or (reverse_lig in entries_list):
                continue
            else:
                matrices[rec_id][knt] = float(value)
                entries_list.append(key)
                knt += 1
            if first_lig not in ligands[rec_id]:
                ligands[rec_id].append(first_lig)
            if second_lig not in ligands[rec_id]:
                ligands[rec_id].append(second_lig)
        print(f"receptor {rec_id} matrix {matrices[rec_id]}")
    return matrices, ligands


def cluster_distance_matrix(matrix, ligands, key):
    """Apply hierarchical clustering."""
    print(f"creating dendrogram for {key}")
    Z = linkage(matrix, LINKAGE)
    dendrogram_figure_filename = key + "_" + LINKAGE + ".png"
    plt.figure()
    dn = dendrogram(Z, color_threshold=THRESHOLD, labels=ligands)
    plt.savefig(dendrogram_figure_filename)
    plt.close()
    # clustering
    clusters = fcluster(Z, t=THRESHOLD, criterion="distance")
    print(f"clusters = {clusters}")
    return clusters


def write_clusters(clusters, ligands, cl_filename):
    """writes clusters to file."""
    print(f"Writing clusters to file {cl_filename}")
    with open(cl_filename, "w") as wfile:
        for key in clusters.keys():
            cl_list = clusters[key]
            lig_list = ligands[key]
            wfile.write(f"Receptor {key}\n")
            cl_dict = {}
            for cl in range(len(cl_list)):
                if clusters[key][cl] not in cl_dict.keys():
                    cl_dict[cl_list[cl]] = [lig_list[cl]]
                else:
                    cl_dict[cl_list[cl]].append(lig_list[cl])
            # writing
            unique_cls = np.unique(clusters[key])
            for cl_id in unique_cls:
                cl_string = " ".join(cl_dict[cl_id])
                wfile.write(f"Cluster {cl_id} -> " + cl_string + os.linesep)


def main():
    if os.path.exists(int_filename) is True:
        matrices, ligands = read_int_matrix(int_filename)
    else:
        raise Exception(f"input path {int_filename} does not exist!")
    # clustering
    clusters = {}
    for el in matrices.keys():
        matrix_filename = el + "_matrix.txt"
        np.savetxt(matrix_filename, matrices[el], fmt="%.6lf")
        clusters[el] = cluster_distance_matrix(matrices[el], ligands[el], el)
    # write clusters
    cl_filename = "clusters_" + LINKAGE + "_thr-" + str(THRESHOLD) + ".out"
    write_clusters(clusters, ligands, cl_filename)


main()
