# This script calculates the sizes of the the molecules in a given benchmark. It then compares the sizes of positive
# and negatives ligand to check if they are within a similitude threshold set by user. If the comparison is successful
# the pair of ligands will be copied in a new benchmark under the name that the user has decided for it.
#
# As output files, it generates an histogram for the distribution of sizes of all the molecules and another one
# comparing the sizes of positive ligands vs negative ones. It also creates two files stroing the sizes of the molecules:
# sizes_dic.txt which is a python dictionary that can be used in other scripts and sizes_table.txt that contains the
# same information but with a much better visualization.
#
# Finally, groups_table.txt is a file generated that contains the information of similar size positive-negative ligands for
# the same receptor. It will be needed later when the analysis to know the pairs of similar size.
#
# Similitude threshold can accept a percentage or a decimal value: 90 = 90% = 0.9
#
# This script does not support interface clustered benchmarks as input, so if wanted to filter a clustered benchmark, it is
# advisable to first filter by size and then proceed with the interface clustering.
#
# Execution: python size_pairs_benchmark.py location/of/the/benchmark name_of_new_benchmark similitude_threshold

# e.g. $ python ./scripts/size_pairs_benchmark.py ./benchmarks/BM-230222/ NEW_BENCHMARK 90
#
# MSc. Jesús López Rivera, Bonvin Lab :)


import argparse
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from brokenaxes import brokenaxes
from tabulate import tabulate


def aminoacid_count(pdb_file):
    """Counting the number of residues in a pdb file"""
    pdb = open(pdb_file, "r").read()
    pdb_lines = pdb.splitlines()
    num_residues_pdb = 0
    for atom in range(
        0, len(pdb_lines) - 1
    ):  # Iterate every atom and the following one
        atom1 = pdb_lines[atom]
        atom2 = pdb_lines[atom + 1]
        atom1_features = atom1.split()
        atom2_features = atom2.split()
        if (
            len(atom1_features) > 1 and len(atom2_features) > 1
        ):  # To avoid last lines with END and empty spaces
            if atom1_features[0] == "ATOM" and (
                atom2_features[0] == "ATOM" or atom2_features[0] == "TER"
            ):
                if int(atom1_features[1]) < int(
                    atom2_features[1]
                ):  # To deal with more than one ensemble and discontinuities
                    if "CA" in str(atom1):
                        num_residues_pdb += 1

    return num_residues_pdb


def compare_size(num_residues_pos_lig, num_residues_neg_lig, similitude):
    """Compare if the size of two molecules is similar given a similitude percentage"""

    if similitude == 0:  # No need to check
        return True

    else:
        if float(num_residues_pos_lig / num_residues_neg_lig) >= similitude and (
            num_residues_pos_lig / num_residues_neg_lig
        ) <= float(1 / similitude):
            return True  # They have similar size
        else:
            return False


def size_benchmark_histogram(directory, sizes_list, benchmark_name):
    """Create a histogram of the sizes of the ligands in a benchmark with broken axes if needed"""

    # Check if it is needed to habe broken axes
    braxes = [100]  # Minimum point to break axes
    xlims = []
    for s in range(0, len(sizes_list) - 1):
        size1 = sorted(sizes_list)[s]  # Smallest size
        size2 = sorted(sizes_list)[(s + 1)]  # Following smallest size
        if (
            size1 + 500
        ) < size2:  # If the difference is bigger than 500, introduce a break in the axis
            braxes += [size1, size2]
    braxes += [max(sizes_list) + 100]  # Maximum point to break axes

    # Take only the even breaks of the axis and create the limits
    for break_axis in range(0, len(braxes) - 1):
        if break_axis % 2 == 0:
            xlim = (braxes[break_axis] - 100, braxes[break_axis + 1] + 100)
            xlims += [xlim]

    # Plot the histogram
    fig = plt.figure(figsize=(8, 5))
    bax = brokenaxes(xlims=tuple(xlims), wspace=0.05)
    bax.hist(
        sizes_list,
        color="blue",
        edgecolor="black",
        bins=np.arange(min(sizes_list), max(sizes_list) + 50, 50),
    )
    bax.set_xlabel("Size (Residues)", labelpad=25)
    bax.set_ylabel("Proteins")
    bax.set_title(f"{benchmark_name} Distribution of sizes")
    plt.savefig(f"{str(directory)}/size_distribution.png")
    plt.close()

    return


def pos_neg_benchmark_histogram(directory, pos_table, neg_table, benchmark_name):
    """Create a histogram of the sizes of the positive-negative ligands in a benchmark"""

    # Check if it is needed to habe broken axes
    all_sizes = pos_table + neg_table  # To obtain the breaks in the axes
    braxes = [100]  # Minimum point to break axes
    xlims = []
    for s in range(0, len(all_sizes) - 1):
        size1 = sorted(all_sizes)[s]  # Smallest size
        size2 = sorted(all_sizes)[(s + 1)]  # Following smallest size
        if (
            size1 + 500
        ) < size2:  # If the difference is bigger than 500, introduce a break in the axis
            braxes += [size1, size2]
    braxes += [max(all_sizes) + 100]  # Maximum point to break axes

    # Take only the even breaks of the axis and create the limits
    for break_axis in range(0, len(braxes) - 1):
        if break_axis % 2 == 0:
            xlim = (braxes[break_axis] - 100, braxes[break_axis + 1] + 100)
            xlims += [xlim]

    # Plot the histogram
    fig = plt.figure(figsize=(8, 5))
    bax = brokenaxes(xlims=tuple(xlims), wspace=0.05)
    bax.hist(
        pos_table,
        bins=np.arange(min(all_sizes), max(all_sizes) + 50, 50),
        alpha=0.4,
        label="Positive ligands",
    )
    bax.hist(
        neg_table,
        bins=np.arange(min(all_sizes), max(all_sizes) + 50, 50),
        alpha=0.4,
        label="Negative ligands",
    )
    bax.set_title(f"{benchmark_name} Distribution of sizes of ligands")
    bax.set_xlabel("Size (Residues)", labelpad=25)
    bax.set_ylabel("Ligands")
    bax.legend(loc="upper right")

    plt.savefig(f"{str(directory)}/pos_neg_size_distribution.png")
    plt.close()
    return


def benchmark_table(directory, output_table, filename):
    """Create a text file with the output table of a benchmark"""

    table_address = f"{str(directory)}/{filename}"
    size_table = open(table_address, "w")
    size_table.write(tabulate(output_table))
    size_table.close()
    return


def write_dic(directory, dictionary, filename):
    """Write a text file for a dictionary"""

    address = f"{str(directory)}/{filename}"
    groups_table = open(address, "w")
    groups_table.write(str(dictionary))
    groups_table.close()
    return


if __name__ == "__main__":

    # Different libraries to store information.
    partners_library = {}
    structures_library = {}

    # Process the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("og_benchmark", help="Introduce benchmark location")
    parser.add_argument("new_benchmark", help="Introduce size-filtered benchmark name")
    parser.add_argument(
        "similitude", help="Introduce percentage of similitude of size. E.g: 0.9 (90%)"
    )
    args = parser.parse_args()
    og_bm_directory = Path(args.og_benchmark)  # Original benchmark path
    og_bm_name = og_bm_directory.stem
    size_bm_name = str(args.new_benchmark)  # New benchmark name
    og_bm_parent_folder = og_bm_directory.parent
    if str(og_bm_parent_folder) == str(
        og_bm_directory
    ):  # The new benchmark is created in the og benchmark parent folder
        og_bm_parent_folder = Path("./..")
    size_bm = f"{og_bm_parent_folder}/{size_bm_name}"  # New benchmark path

    # Process the similitude. If the user gives a percentage value, e.g. 90, transform it to decimal = 0.9
    similitude = float(args.similitude)
    if similitude > 1:
        similitude = similitude / 100

    """ Loop each folder in the benchmark to access the pdb files in them """
    for folder in og_bm_directory.iterdir():
        folder_name = folder.stem
        molecules = folder_name.split("_")  # [Receptor, Ligand, "positive/negative"]
        if len(molecules) == 3 and folder.is_dir():  # Avoid files and other folders
            receptor = molecules[0]
            ligand = molecules[1]
            interaction = molecules[2]
            molecule_pair = f"{receptor}_{ligand}"

            """ Each positive ligand needs to be compared with each negative ligand for each receptor,
            so it is needed a library in which it is stored every positive and negative ligand for every receptor"""
            if receptor not in partners_library.keys():
                partners_library[receptor] = {"positive": [], "negative": []}
            if ligand not in partners_library[receptor][interaction]:
                partners_library[receptor][interaction] += [ligand]

            # partners_library = {
            #               receptor1: {
            #                   "positive":[pos_ligand1, pos_ligand2, ...], "negative": [neg_ligand1, ...],
            #                       },
            #               receptor2: {...},
            #                ...}

            """ Create a library of matching the Uniprot_id to the PDB structure path """
            for file in folder.iterdir():
                if str(file).endswith(".pdb"):  # Avoid other folders or files
                    file_name = file.name
                    structures = file_name.split(
                        "_"
                    )  # ["12ab", "r", "b.pdb"] or ["45de", "positive"/"negative", "l", "b.pdb"]
                    if structures[2] == "l":  # It's a ligand
                        structures_library[ligand] = file
                    elif structures[1] == "r":  # It's a receptor
                        structures_library[receptor] = file
                    # structures_library = {
                    #           protein_id1: pdb_structure_route1, protein_id2: pdb_structure_route2,
                    #                   ...}

    """Create a library of similar size positive and negative ligands and a table with the sizes of all,
    one for the original benchmark and one for the filtered one"""

    size_benchmark = {}

    output_table_og_bm = [
        ["Receptor_ID", "Ligand_ID", "Interaction", "Size"]
    ]  # For better visualization of results
    output_table_size_bm = [
        ["Receptor_ID", "Ligand_ID", "Interaction", "Size"]
    ]  # For better visualization of results

    # Lists of molecules sizes for the histogram distribution of sizes
    sizes_og_bm = []
    sizes_size_bm = []

    # Lists of positive/negative sizes for the size distribution of positive vs negative ligands in the original and filtered benchmark
    pos_list_og_bm = []
    neg_list_og_bm = []
    pos_list_size_bm = []
    neg_list_size_bm = []

    pos_sizes_og_bm = []
    neg_sizes_og_bm = []
    pos_sizes_size_bm = []
    neg_sizes_size_bm = []

    pair_ligands = {}  # Dictionary to create groups_table.txt
    og_sizes_dic = {}  # Store already counted proteins.
    new_sizes_dic = {}  # Dictionary of sizes for the new filtered benchmark

    # Count number of pairs to know the number of analysis pairs
    og_bm_n_pairs = 0
    size_bm_n_pairs = 0

    """ Iterate every receptor """
    for receptor in partners_library.keys():
        if receptor not in og_sizes_dic.keys():  # Not processed before
            # Calculate size of receptor
            receptor_pdb = structures_library[receptor]  # PDB structure of the receptor
            num_residues_rec = aminoacid_count(receptor_pdb)  # Total number of residues
            og_sizes_dic[receptor] = num_residues_rec
        else:  # Already processed
            num_residues_rec = og_sizes_dic[receptor]
        # Write original benchmark output table
        table_line_og_bm = [receptor, "-", "-", num_residues_rec]
        output_table_og_bm += [table_line_og_bm]

        """ Iterate every positive ligand for the receptor """
        for positive in partners_library[receptor]["positive"]:
            if positive not in og_sizes_dic.keys():  # Not processed before
                # Calculate size of positive ligand
                positive_ligand = structures_library[
                    positive
                ]  # PDB structure of the postive
                num_residues_pos = aminoacid_count(
                    positive_ligand
                )  # Total number of residues
                og_sizes_dic[positive] = num_residues_pos
            else:  # Already processed
                num_residues_pos = og_sizes_dic[positive]
            # Write original benchmark output tables
            table_line_og_bm = [receptor, positive, "positive", num_residues_pos]
            output_table_og_bm += [table_line_og_bm]
            if (
                positive not in pos_list_og_bm
            ):  # Add to list of positives for the size distribution
                pos_list_og_bm += [positive]

            # Iterate every negative ligand for each positive ligand
            for negative in partners_library[receptor]["negative"]:
                og_bm_n_pairs += 1  # A positive-negative analysis pair
                if negative not in og_sizes_dic.keys():  # Not processed before
                    # Calculate size of negative ligand
                    negative_ligand = structures_library[
                        negative
                    ]  # PDB structure of the negative
                    num_residues_neg = aminoacid_count(
                        negative_ligand
                    )  # Total number of resides
                    og_sizes_dic[negative] = num_residues_neg
                else:  # Already processed
                    num_residues_neg = og_sizes_dic[negative]
                # Write original benchmark output table
                table_line_og_bm = [receptor, negative, "negative", num_residues_neg]
                # Negative ligands are repeated for every positive, so they are not rewritten if they are already
                if table_line_og_bm not in output_table_og_bm:
                    output_table_og_bm += [table_line_og_bm]

                # output_table_og_bm = [
                #               ["Receptor_ID", "Ligand_ID", "Interaction", "Size"],
                #               [receptor1, "     -   ", "    -   ", size1]
                #               [receptor1, pos_ligand1, "positive", size2],
                #               [receptor1, neg_ligand1, "negative", size3],
                #               [receptor2, neg_ligand2, "negative", size4],
                #               ...]

                if (
                    negative not in neg_list_og_bm
                ):  # Add to list of negatives for the size distribution
                    neg_list_og_bm += [negative]

                """ Check if sizes of ligands are similar """
                if compare_size(
                    num_residues_pos, num_residues_neg, similitude
                ):  # They are!
                    size_bm_n_pairs += (
                        1  # It's an analysis pair in the size-filtered benchmark
                    )

                    # Write size_benchmark dictionary
                    positive_pair = f"{positive}_{negative}"
                    if receptor not in size_benchmark.keys():  # Not processed before
                        size_benchmark[receptor] = []
                    size_benchmark[receptor] += [positive_pair]
                    # size_benchmark = {
                    #           receptor1:[positive1-negative1, positive1-negative2, ...],
                    #           receptor2:[positive1-negative1],
                    #           ... }

                    # Write filtered benchmark output table
                    table_line_size_bm = [receptor, "-", "-", num_residues_rec]
                    if (
                        table_line_size_bm not in output_table_size_bm
                    ):  # Avoid repeating the receptor line for every pair
                        output_table_size_bm += [table_line_size_bm]

                    table_line_size_bm = [
                        receptor,
                        positive,
                        "positive",
                        num_residues_pos,
                    ]
                    if (
                        table_line_size_bm not in output_table_size_bm
                    ):  # Avoid repeating the positive line
                        output_table_size_bm += [table_line_size_bm]
                    if (
                        positive not in pos_list_size_bm
                    ):  # Add to list of positives for the size distribution
                        pos_list_size_bm += [positive]

                    table_line_size_bm = [
                        receptor,
                        negative,
                        "negative",
                        num_residues_neg,
                    ]
                    if (
                        table_line_size_bm not in output_table_size_bm
                    ):  # Avoid repeating the negative line
                        output_table_size_bm += [table_line_size_bm]
                    if (
                        negative not in neg_list_size_bm
                    ):  # Add to list of negatives for the size distribution
                        neg_list_size_bm += [negative]

                    # Write new_sizes_dic for later save it in a file
                    if receptor not in new_sizes_dic.keys():
                        new_sizes_dic[receptor] = num_residues_rec
                    if positive not in new_sizes_dic.keys():
                        new_sizes_dic[positive] = num_residues_pos
                    if negative not in new_sizes_dic.keys():
                        new_sizes_dic[negative] = num_residues_neg

                    # output_table_size_bm = [
                    #               ["Receptor_ID", "Ligand_ID", "Interaction", "Size"],
                    #               [receptor1, pos_ligand1, positive, size1],
                    #               [receptor1, neg_ligand1, negative, size2],
                    #                   ...]

                    """Knowing the pairs will be useful in the analysis, so it is useful to create a 
                    dictionary to store this information -> groups_table.txt"""
                    rec_lig_positive = f"{receptor}_{positive}_positive"
                    rec_lig_negative = f"{receptor}_{negative}_negative"
                    if rec_lig_positive not in pair_ligands.keys():
                        pair_ligands[rec_lig_positive] = []
                    pair_ligands[rec_lig_positive] += [rec_lig_negative]
                    # pair_ligands_pos = {
                    #           receptor1_positive1:[receptor1_negative1, receptor1-negative2],
                    #           receptor1_positive2:[receptor1_negative3],
                    #               ...}

    """ Create the size filtered benchmark from the data in size_benchmark """
    for receptor in size_benchmark.keys():
        for ligand_pair in size_benchmark[receptor]:
            ligands = ligand_pair.split("_")
            pos_folder_name = f"{receptor}_{ligands[0]}_positive"
            og_folder_pos = f"{str(og_bm_directory)}/{pos_folder_name}"  # Original folder of the benchmark
            size_folder_pos = f"{size_bm}/{pos_folder_name}"  # Destination folder in the filtered benchmark
            shutil.copytree(
                og_folder_pos, size_folder_pos, dirs_exist_ok=True
            )  # Avoid repeated folders
            neg_folder_name = f"{receptor}_{ligands[1]}_negative"
            og_folder_neg = f"{str(og_bm_directory)}/{neg_folder_name}"  # Original folder of the benchmark
            size_folder_neg = f"{size_bm}/{neg_folder_name}"  # Destination folder in the filtered benchmark
            shutil.copytree(
                og_folder_neg, size_folder_neg, dirs_exist_ok=True
            )  # Avoid repeated folders

    """ Create the distribution of sizes histograms of both benchmarks """
    for (
        protein
    ) in og_sizes_dic.keys():  # Add the size of every protein in the original dataset
        sizes_og_bm += [og_sizes_dic[protein]]
    for (
        protein
    ) in new_sizes_dic.keys():  # Add the size oof every protein in the filtered dataset
        sizes_size_bm += [new_sizes_dic[protein]]

    size_benchmark_histogram(og_bm_directory, sizes_og_bm, og_bm_name)
    size_benchmark_histogram(size_bm, sizes_size_bm, size_bm_name)

    """ Create the positive vs negative distribution of sizes of both benchmarks """
    for (
        positive
    ) in pos_list_og_bm:  # Add the size of every positive in the original dataset
        pos_sizes_og_bm += [og_sizes_dic[positive]]
    for (
        negative
    ) in neg_list_og_bm:  # Add the size of every negative in the original dataset
        neg_sizes_og_bm += [og_sizes_dic[negative]]

    for (
        positive
    ) in pos_list_size_bm:  # Add the size of every positive in the filtered dataset
        pos_sizes_size_bm += [new_sizes_dic[positive]]
    for (
        negative
    ) in neg_list_size_bm:  # Add the size of every negative in the filtered dataset
        neg_sizes_size_bm += [new_sizes_dic[negative]]

    pos_neg_benchmark_histogram(
        og_bm_directory, pos_sizes_og_bm, neg_sizes_og_bm, og_bm_name
    )
    pos_neg_benchmark_histogram(
        size_bm, pos_sizes_size_bm, neg_sizes_size_bm, size_bm_name
    )

    """ Create the table for both benchmarks """
    benchmark_table(og_bm_directory, output_table_og_bm, "size_table.txt")
    benchmark_table(size_bm, output_table_size_bm, "size_table.txt")

    """ Generate the dictionary files for the analysis """
    write_dic(size_bm, pair_ligands, "groups_table.txt")
    write_dic(og_bm_directory, og_sizes_dic, "sizes_dic.txt")
    write_dic(size_bm, new_sizes_dic, "sizes_dic.txt")

    """ Finally, print the number of positive-negative analysis pairs """
    print(f"There are {og_bm_n_pairs} pairs in the original benchmark")
    print(f"There are {size_bm_n_pairs} pairs in the size-filtered benchmark.")
