# This script obtains the ranked solutions from Haddock runs stored in /structures/it*/file.list files,
# from an analysis pair receptor+positive_ligand - receptor+negative_ligand, combines them and sort the solutions
# according to their score and creates a combined file. If size_pairs_benchmark.py has been run
# before, the pairs are stored in the file groups_table.txt. If not, it can determine the pairs by default.
#
# For the scores, it sets different scenarios: Haddock score (HS), Haddock score by interface size (HSIS)
# and Haddock score by average interface size for different number of solutions (HSIAS_n)
#
# After that, it takes the previously generated files and calculates the Right Answer Rate (RAR) for different
# thresholds (1, 10, 50, 100, 200). The RAR is defined as the number of positive solutions divided by the number
# of total solutions under a certain threshold. E.g: under threshold 10, there are 6 positive and 4 negative
# solutions, so the RAR = 60%.
#
# Finally, it calculates the average RAR for each threshold for the whole benchmark as well as its standard
# deviation and plots a bar graph with that information. All the files and graphs will be stored in a new
# folder under the benchmark location with the name analysis/.
#
# The mode -full creates extra plots for each RAR threshold. These plots have the receptor in the x-axis sorted
# by size from smallest to biggest. On the y-axis there is the RAR for each of the analysis pairs. Each of the
# dot sizes represent the average size of both the ligands (it is useful if they have siilar size, but no point
# if they don't). Finally, the color of the dots represent the average energy (HADDOCK score) of the top threshold
# solutions and the thickness of the border, the variation of energy among solutions.
#
# Execution: python interaction_pair_ana.py [-full] location/of/the/benchmark [location/of/groups_table.txt] [location/of/sizes_dic.txt]
# e.g. $ python ./scripts/interaction_pairs_ana.py [-full] ./benchmarks/OutHM-220622/BM-AO-220622 [groups_table] [sizes_dic]
#
# MSc. Jesús López Rivera, Bonvin Lab :)


import argparse
import ast
import logging
import os
import shutil
from pathlib import Path
from statistics import stdev

import matplotlib.pyplot as plt

FORMAT = " %(asctime)s L%(lineno)d %(levelname)s - %(message)s"
logging.basicConfig(format=FORMAT, level="INFO")


def get_ligands_pair(bm_directory):
    """Define pairs of positive-negative pairs for every receptor in the dataset"""

    partners_library = {}
    ligands_pair = {}
    clusters_library = {}
    # Proceding in default mode.
    for folder in bm_directory.iterdir():
        folder_name = folder.stem
        group = folder_name.split("-")[0]
        molecules = group.split("_")  # [Receptor, Ligand, "positive"]
        if len(molecules) == 3 and folder.is_dir():  # Avoid files and other folders
            if group not in clusters_library.keys():
                clusters_library[group] = {}
            receptor = molecules[0]
            ligand = molecules[1]
            interaction = molecules[2]

            it0_file_list = f"{str(folder)}/run-interface_all/structures/it0/file.list"
            if Path(it0_file_list).exists():
                clusters_library[group][folder_name] = True
                """ Each positive ligand needs to be compared with each negative ligand for each receptor,
                so it is needed a library in which it is stored every positive and negative ligand for every receptor"""
                if receptor not in partners_library.keys():
                    partners_library[receptor] = {"positive": [], "negative": []}
                if ligand not in partners_library[receptor][interaction]:
                    partners_library[receptor][interaction] += [ligand]
            else:
                clusters_library[group][folder_name] = False

    # Iterate the previously generated dictionary to create the groups
    del_groups = []
    del_receptors = []
    for group in clusters_library.keys():
        remove_group = False
        for cluster in clusters_library[group].keys():
            if not clusters_library[group][cluster]:
                remove_group = True
                break
        if remove_group:
            del_groups += [group]

    for group in del_groups:
        molecules = group.split("_")
        receptor = molecules[0]
        ligand = molecules[1]
        interaction = molecules[2]
        try:
            partners_library[receptor][interaction].remove(ligand)
        except:
            pass

    for receptor in partners_library.keys():
        if (
            not len(partners_library[receptor]["positive"]) == 0
            and not len(partners_library[receptor]["negative"]) == 0
        ):
            for positive in partners_library[receptor]["positive"]:
                rec_lig_positive = f"{receptor}_{positive}_positive"
                if rec_lig_positive not in ligands_pair.keys():
                    ligands_pair[rec_lig_positive] = []
                for negative in partners_library[receptor]["negative"]:
                    rec_lig_negative = f"{receptor}_{negative}_negative"
                    if rec_lig_negative not in ligands_pair[rec_lig_positive]:
                        ligands_pair[rec_lig_positive] += [rec_lig_negative]
        # ligands_pair: {
        #   'Receptor1_Positive1_positive': ['Receptor1_Negative1_negative', 'Receptor1_Negative2_negative'],
        #   'Receptor1_Positive2_positive': [...],
        #  ...}
        else:
            del_receptors += [receptor]

    for receptor in del_receptors:
        try:
            del partners_library[receptor]
        except:
            pass

    return ligands_pair, partners_library


def aminoacid_count(pdb_file):
    """Counting the number of residues in a pdb file"""

    pdb = open(pdb_file, "r").read()
    pdb_lines = pdb.splitlines()
    num_residues_pdb = 0
    for atom in range(0, len(pdb_lines) - 1):
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


def get_benchmark_sizes(bm_directory):
    """Calculate the size of the different molecules in the benchmark"""

    sizes_dic = {}
    for folder in bm_directory.iterdir():
        folder_name = str(folder.stem)
        molecules = folder_name.split("_")  # [Receptor, Ligand, "positive(-cluster-n)"]
        if len(molecules) == 3 and folder.is_dir():  # Avoid files and other folders
            receptor = molecules[0]
            ligand = molecules[1]

            for file in folder.iterdir():
                file_name = file.name
                if str(file).endswith(".pdb"):
                    file_name = file.name
                    structures = file_name.split(
                        "_"
                    )  # ["12ab", "r", "b.pdb"] or ["45de", "positive"/"negative", "l", "b.pdb"]
                    if structures[2] == "l":
                        if ligand not in sizes_dic.keys():
                            sizes_dic[ligand] = aminoacid_count(file)
                    elif structures[1] == "r":
                        if receptor not in sizes_dic.keys():
                            sizes_dic[receptor] = aminoacid_count(file)

    return sizes_dic


def get_BSA(file_list, pdb_solution):
    "Get the buried surface area for a certain pdb solution"

    bsa = False
    file_list_folder = Path(file_list).parent
    pdb_solution_name = (
        str(pdb_solution).split(":")[-1].strip('"')
    )  # Receptor1_Ligand1_positive/negative-(cluster-n)_94.pdb
    pdb_solution_file = f"{file_list_folder}/{pdb_solution_name}"

    read_pdb_solution = open(pdb_solution_file, "r").read()
    pdb_solution_lines = read_pdb_solution.splitlines()
    for line in pdb_solution_lines:
        if "buried surface area" in line:
            bsa = float(line.split()[-1])

    return bsa


def get_av_bsa(file_list, approach):
    """Calculate the average buried surface area for a given number of solutions"""

    bsas = []
    score_results = open(file_list, "r").read()
    results_lines = score_results.splitlines()
    for line in results_lines[0 : int(approach)]:
        results_list = str(
            line
        ).split()  # ['"PREVIT:Receptor1_Ligand1_positive/negative-(cluster-n)_94.pdb"', '{', 'score', '}']
        bsa = get_BSA(file_list, results_list[0])
        bsas += [bsa]
    av_bsa = float(sum(bsas) / len(bsas))

    return av_bsa


def results_library(it_path, results, scores, approach):
    """Read file.list files and store the results"""

    if "HSAIS" in approach:
        # Need to calculate the average BSA por the given number of solutions
        approach = int(approach.split("_")[-1])
        av_bsa = get_av_bsa(it_path, approach)

    score_results = open(it_path, "r").read()
    results_lines = score_results.splitlines()

    for line in results_lines:  # Each line is a solution with its HADDOCK score
        # Calculate the score according to the approach
        results_list = str(line).split()
        if approach == "HS":  # HADDOCK score
            score = float(results_list[2])
            line = f"{results_list[0]} | HS: {str(score)} |"
        elif approach == "HSIS":  # HADDOCK score by Interface Size
            bsa = get_BSA(it_path, results_list[0])
            score = float(results_list[2]) / bsa
            line = f"{results_list[0]} | HSIS: {str(score)} | HS: {results_list[2]} |"
        else:  # HADDOCK score by Average Interface Size
            score = float(results_list[2]) / av_bsa
            line = f"{results_list[0]} | HSAIS: {str(score)} | HS: {results_list[2]} |"

        # Store the scores
        results[score] = line
        scores += [score]

    return results, scores


def sort_results(directory, step, results, scores, group):
    """Sort the results obtained from the negative and positive file.list and stored in results and scores."""

    sorted_results = ""  # Text that will go inside the file
    n = 1  # This is optional, but it is useful for easier visualisation of the ranking.
    for score in sorted(scores):
        # Add the line for each score to the text
        sorted_results += f"{str(n)} {results[score]} \n"
        n += 1

    # Store the file in benchmark/analysis/scenario/scoring_approach/step/file.txt
    itn_ana_file_address = f"{str(directory)}/{step}_sorted_results_{group}.txt"
    itn_ana_file = open(itn_ana_file_address, "w")
    itn_ana_file.write(sorted_results)
    itn_ana_file.close()
    return


def right_answer_rate(itn_file, group_ID, output_table):
    """Calculate the right answer rate of a sorted results file under different thresholds"""
    output_dic = {}
    itn_results = open(itn_file, "r").read()
    itn_lines = itn_results.splitlines()
    thresholds = [1, 10, 50, 100, 200]
    rates_line = [group_ID]
    for threshold in thresholds:
        solution = 0
        pos_count = 0
        while solution < threshold:
            score_features = itn_lines[solution].split(
                "_"
            )  # ['1 "PREVIT:Receptor1', 'Ligand1', 'positive/negative', '94.pdb" { score }']
            if score_features[2].split("-")[0] == "positive":
                pos_count += 1
            solution += 1
        rate_threshold = float(pos_count) / float(solution) * 100
        output_dic[str(threshold)] = rate_threshold
        rates_line += [rate_threshold]
    output_table += [rates_line]
    return output_dic, output_table


def group_energies(itn_file):
    """Calculate the average energy of a sorted results file under different thresholds"""

    output_dic = {}
    itn_results = open(itn_file, "r").read()
    itn_lines = itn_results.splitlines()
    thresholds = [1, 10, 50, 100, 200]
    for threshold in thresholds:
        solution = 0  # Count number of solutions until threshold
        energies_list = []
        while solution < threshold:
            score_features = itn_lines[solution].split()
            energies_list += [float(score_features[4])]  # HADDOCK score
            solution += 1
        energy_average = float(sum(energies_list)) / float(len(energies_list))
        if len(energies_list) == 1:
            energy_sd = 0  # Avoid crashing with threshold 1
        else:
            # Calculate the standard deviation
            energy_sd = float(stdev(energies_list)) / float(energy_average) * 1.2
            # 1.2 is a correcting factor for better visualization of the std in the plot
        output_dic[str(threshold)] = {"average": energy_average, "std": energy_sd}
    return output_dic


def average_sd(Tn):
    """Calculate the average and standard deviation of a given set of measures"""
    logging.info(f" There are {len(Tn)} analysis pairs")
    average_Tn = float(sum(Tn)) / float(len(Tn))
    sd_Tn = stdev(Tn)
    return [average_Tn, sd_Tn]


def average_rar(table, scenario, step, approach, folder):
    """Plot a bar graph of the average right answer rates for different thresholds"""
    T1 = []
    T10 = []
    T50 = []
    T100 = []
    T200 = []
    final_table = [["Success_rate", "SD"]]

    # Store the RAR values for each threshold separatedly
    for group in table[1:]:
        T1 += [group[1]]
        T10 += [group[2]]
        T50 += [group[3]]
        T100 += [group[4]]
        T200 += [group[5]]

    # Calculate the average and standard deviation for each threshold and store the results in final table
    final_table += [average_sd(T1)]
    final_table += [average_sd(T10)]
    final_table += [average_sd(T50)]
    final_table += [average_sd(T100)]
    final_table += [average_sd(T200)]

    # Plot a bar graph with the previous data
    labels = ["T1", "T10", "T50", "T100", "T200"]
    rar_means = []
    rar_std = []
    for rar in final_table[1:]:
        rar_means += [rar[0]]
        rar_std += [rar[1]]
    width = 0.85
    fig, ax = plt.subplots()
    plt.ylim(top=110)
    ax.bar(
        labels,
        rar_means,
        width,
        yerr=rar_std,
        align="center",
        alpha=0.5,
        color="yellowgreen",
        edgecolor="black",
        capsize=10,
    )
    ax.set_ylabel("Right Answer Rate (%)")
    ax.set_title(f"Average RAR {scenario}_{step}_{approach}")
    # The place to save it is: bm_directory/analysis/scenario/step/scoring_approach/
    plt.savefig(f"{folder}/average_rar_{scenario}_{step}_{approach}.png")
    plt.close()
    return


def write_dic(directory, dictionary, filename):
    """Write a text file for a dictionary"""
    address = f"{str(directory)}/{filename}"
    groups_table = open(address, "w")
    groups_table.write(str(dictionary))
    groups_table.close()
    return


if __name__ == "__main__":

    # Dictionary to store the paths to the file.list files for each Receptor_Ligand_positive/negative.
    molecule_to_results = {}

    scenarios = ["interface_all"]  # Different run setups of HADDOCK
    steps = ["it0"]  # Different HADDOCK steps
    rar_thresholds = ["1", "10", "50", "100", "200"]
    approaches = ["HS", "HSIS", "HSAIS_10", "HSAIS_50", "HSAIS_100", "HSAIS_200"]

    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", help="Introduce benchmark location")
    parser.add_argument(
        "groups_table",
        nargs="?",
        help="Introduce groups_table.txt file if previously filtered by size",
    )
    parser.add_argument(
        "sizes_dic",
        nargs="?",
        help="Introduce sizes_dic.txt file if previously filtered by size",
    )
    parser.add_argument(
        "-full",
        action="store_true",
        help="Generate a detailed plot for every threshold showing the individual groups",
    )
    args = parser.parse_args()
    bm_directory = Path(args.benchmark)
    groups_table = args.groups_table
    sizes_dic_address = args.sizes_dic
    full = args.full  # Activates full analysis

    """groups_table.txt is a file generated after size_benchmark.py which contains the groups 
    receptor-positive_ligand-negative_ligand. It is important to have it if the dataset has been filtered by similar
    size of positive-negative ligands, since the groups will be different."""
    if groups_table:
        if Path(groups_table).exists():
            logging.info("<<< groups_table.txt found >>>")
            read_groups_table = open(groups_table, "r").read()
            ligands_pair = ast.literal_eval(read_groups_table)
        # ligands_pair = {Receptor1_Ligand1_positive: [Receptor1_Ligand2_negative, Receptor1_Ligand3_negative],
        #               Receptor 2_Ligand4_positive: [Receptor2_Ligand5_negative], ...}
        else:
            logging.warning(
                "<<< groups_table.txt not found, proceding to deafult mode >>>"
            )
            ligands_pair, partners_library = get_ligands_pair(bm_directory)
    else:
        logging.info("<<< groups_table.txt not provided, proceding to defautl mode >>>")
        ligands_pair, partners_library = get_ligands_pair(bm_directory)
        logging.info(ligands_pair)

    # Store all different receptors
    receptors = []
    for ligand_pair in ligands_pair.keys():
        rec = ligand_pair.split("_")[0]
        if rec not in receptors:
            receptors += [rec]

    if full:
        """sizes_dic.txt is a file generated after size_benchmark.py which contains the sizes of the different
        molecules in the benchmark, if it has been done, the script can be faster avoiding that step."""
        if sizes_dic_address:
            if Path(sizes_dic_address).exists():
                logging.info("<<< sizes_dic.txt found >>>")
                read_sizes_dic = open(sizes_dic_address, "r").read()
                sizes_dic = ast.literal_eval(read_sizes_dic)
                # sizes_dic = {Protein1: 156, Protein2: 235, Protein3: 201, ...}
            else:
                logging.info("<<< sizes_dic.txt not found, calculating sizes... >>>")
                sizes_dic = get_benchmark_sizes(bm_directory)
        else:
            logging.info("<<< No sizes_dic.txt providing, calculating sizes... >>>")
            sizes_dic = get_benchmark_sizes(bm_directory)

        energies_dic = {}  # To store information on the energies of the solutions later

    """The output information of this script is organised by pos/neg pairs, but the benchmark is structured
    individually, so to better handling it, it creates a new folder to store all the output"""
    analysis_folder = f"{str(bm_directory)}/analysis"
    if Path(analysis_folder).exists():
        shutil.rmtree(analysis_folder)
    os.mkdir(analysis_folder)

    if full:
        """Sort receptors by size"""
        sorted_receptors = []
        size_receptor_dic = {}
        for receptor in partners_library.keys():
            if (
                receptor not in receptors
            ):  # Allows you to check only those that appear in groups_table.txt in case of size-filter
                continue

            receptor_size = sizes_dic[receptor]
            if (
                receptor_size not in size_receptor_dic.keys()
            ):  # Check if it has been processed before
                size_receptor_dic[receptor_size] = []
            size_receptor_dic[receptor_size] += [
                receptor
            ]  # {74:[Receptor1], 124:[Receptor2], ...}

        sorted_receptor_sizes = sorted(size_receptor_dic.keys())  # [124, 74, ...]
        for rec_size in sorted_receptor_sizes:
            for receptor_id in size_receptor_dic[rec_size]:
                if receptor_id not in sorted_receptors:
                    sorted_receptors += [receptor_id]  # [Receptor2, Receptor1, ...]
        """Receptors are sorted by size for an additional visualization in the full analysis mode."""

    for scenario in scenarios:
        logging.info(f"<<<< Processing scenario {scenario} >>>> ")
        # Make a folder for each HADDOCK setup scenario
        ana_scenario_folder = f"{analysis_folder}/{scenario}"
        if not Path(ana_scenario_folder).exists():
            os.mkdir(ana_scenario_folder)

        """ Get the path to file.list for each Receptor_Ligand_positive/negative """
        for folder in bm_directory.iterdir():
            folder_name = folder.stem
            group_name = folder_name.split("-")[0]
            molecules = group_name.split("_")
            if folder.is_dir() and len(molecules) == 3:  # Skip files and other folders
                file_list_addresses = {}
                file_list_addresses[
                    "it0"
                ] = f"{str(folder)}/run-{scenario}/structures/it0/file.list"
                file_list_addresses[
                    "it1"
                ] = f"{str(folder)}/run-{scenario}/structures/it1/file.list"
                file_list_addresses[
                    "itw"
                ] = f"{str(folder)}/run-{scenario}/structures/it1/water/file.list"

                for step in steps:
                    if Path(file_list_addresses[step]).exists:
                        logging.info(f"{file_list_addresses[step]} exists")
                        # Store the paths in a dictionary
                        if group_name not in molecule_to_results.keys():
                            molecule_to_results[group_name] = {
                                "it0": {},
                                "it1": {},
                                "itw": {},
                            }
                        molecule_to_results[group_name][step][
                            folder_name
                        ] = file_list_addresses[step]

                # molecule_to_results = {
                #               receptor1_ligand1_positive: {
                #                           "it0":{receptor1_ligand1_positive(-cluster-n):it0_results_address, ...}
                #                           "it1":{...}
                #                           "itw":{...}
                #                         },
                #               receptor1_ligand1_negative: {...},
                #                   ...}

        for step in steps:
            logging.info(f"<<< Processing step {step} >>>")
            # Make a folder for every HADDOCK step
            ana_step_folder = f"{ana_scenario_folder}/{step}"
            if not Path(ana_step_folder).exists():
                os.mkdir(ana_step_folder)

            for approach in approaches:
                logging.info(f"<< Processing scoring approach {approach} >>")
                # Make a folder for every scoring approach
                ana_approach_folder = f"{ana_step_folder}/{approach}"
                if not Path(ana_approach_folder).exists():
                    os.mkdir(ana_approach_folder)

                for positive in ligands_pair.keys():
                    for negative in ligands_pair[positive]:

                        full_group = f"{positive}_{negative}"  # Receptor1_Ligand1_positive_Receptor1_Ligand2_negative
                        full_group_names = full_group.split("_")
                        receptor_ID = full_group_names[0]
                        positive_ID = full_group_names[1]
                        negative_ID = full_group_names[4]
                        rec_posl_negl = f"{receptor_ID}-{positive_ID}-{negative_ID}"  # For the name of the file

                        logging.info(f"< Processing pair {rec_posl_negl} >")

                        # Variables to store the information
                        results = {}
                        scores = []

                        """The information in file.list looks like:
                        " solution1.pdb { score }
                          solution2.pdb { score }
                        ... "
                        The function results_library extracts the score into the list scores and the information
                        linked to every score into the dictionary results to later sort the scores. """
                        # Get the scores of the positive ligand run
                        for cluster in molecule_to_results[positive][step].keys():
                            results, scores = results_library(
                                molecule_to_results[positive][step][cluster],
                                results,
                                scores,
                                approach,
                            )
                        # Get the scores of the negative ligand run
                        for cluster in molecule_to_results[negative][step].keys():
                            results, scores = results_library(
                                molecule_to_results[negative][step][cluster],
                                results,
                                scores,
                                approach,
                            )

                        # results = {-10.000: '"PREVIT:Receptor1_Ligand1_positive_123.pdb"  { -10.000 }',
                        # -9.000: '"PREVIT:Receptor1_Ligand2_negative_123.pdb"  { -9.800 }', ...}
                        # scores = [-10.000, -9.800, ...]

                        """ Sort results for each positive-negative pair and store them in a file """
                        sort_results(
                            ana_approach_folder, step, results, scores, rec_posl_negl
                        )

                # List to store the right answer rates (RAR) for different Thresholds for each positive/negative pair.
                output_table = [["Group_IDs", "T1", "T10", "T50", "T100", "T200"]]
                rar_dic = {}

                """ Calculate the RAR for different thresholds for each of the files generated before """
                # RAR = n positives / total of solutions(threshold).
                logging.info(f"< Calulating RARs for {scenario}-{step}-{approach} >")
                for file in Path(ana_approach_folder).iterdir():
                    file_name = file.stem
                    file_features = file_name.split(
                        "_"
                    )  # ["itn", "sorted", "results", "Receptor-PosLigand-NegLigand"]
                    if len(file_features) == 4 and str(file).endswith(".txt"):
                        group_name = file_features[-1]
                        group_molecules = group_name.split("-")
                        receptor = group_molecules[0]
                        positive = group_molecules[1]
                        negative = group_molecules[2]
                        pos_neg = f"{positive}_{negative}"

                        if receptor not in rar_dic.keys():
                            rar_dic[receptor] = {}
                        # Store the RARs in the dictionary rar_dic and in output_table. rar_dic is useful for the -full analysis
                        # while output_table is for the basic analysis
                        rar_dic[receptor][pos_neg], output_table = right_answer_rate(
                            file, group_name, output_table
                        )

                        # rar_dic = {
                        #          receptor1:{
                        #                  group1:{
                        #                       '1': 100
                        #                       '10': 70
                        #                       '50': 66
                        #                       '100': 63
                        #                       '1000': 60.8
                        #                           },
                        #                   group2:{...}, ...
                        #                    }
                        #          receptor2:{...}, ...
                        #            }

                        if full:
                            """Get the average energies for the analysis groups for each threshold"""
                            if group_name not in energies_dic.keys():
                                energies_dic[group_name] = group_energies(file)

                            # energies_dic = {
                            #           group1:{
                            #               '1': {
                            #                   'average':-267.5
                            #                   'std': 5.67
                            #               '10': {...},
                            #               ...},
                            #           group2:{...},
                            #              ...}

                # Write a file with the information in rar_dic
                rar_dic_name = f"{step}_{approach}_rar_dic.txt"
                write_dic(ana_approach_folder, rar_dic, rar_dic_name)

                """ Calculate the average and sd for each threshold and each itn, and plot a bar graph """
                average_rar(output_table, scenario, step, approach, ana_approach_folder)

                if full:
                    """Plot a detailed graph containing the information of the individual RARs"""
                    logging.info(f"< Plotting full analysis graphs >")
                    for (
                        rar_threshold
                    ) in rar_thresholds:  # Each threshold has its own plot
                        rar_values = []
                        receptors = []
                        sizes = []
                        energies_av = []  # average
                        energies_sd = []  # standard deviation

                        for (
                            receptor
                        ) in (
                            sorted_receptors
                        ):  # Loop every receptor to get the RAR of the individual groups
                            receptors += [receptor]
                            receptor_rars = []
                            receptor_sizes = []
                            receptor_energies_av = []
                            receptor_energies_sd = []

                            for group in rar_dic[
                                receptor
                            ]:  # Loop every analysis pair (group) to get their sizes and RARs
                                receptor_rars += [
                                    rar_dic[receptor][group][rar_threshold]
                                ]  # Analysis pair RAR
                                positive = group.split("_")[0]
                                positive_size = sizes_dic[
                                    positive
                                ]  # Size of the positive
                                negative = group.split("_")[1]
                                negative_size = sizes_dic[
                                    negative
                                ]  # Size of the negative
                                average_pos_neg_size = (
                                    float(positive_size + negative_size) / 2
                                )  # Average size of the ligands
                                receptor_sizes += [average_pos_neg_size]
                                whole_group = f"{receptor}-{positive}-{negative}"
                                receptor_energies_av += [
                                    energies_dic[whole_group][rar_threshold]["average"]
                                ]  # Average energy of the group.
                                receptor_energies_sd += [
                                    energies_dic[whole_group][rar_threshold]["std"]
                                ]  # Std of the energy of the group.

                            # Tuples are needed for the generation of the plot.
                            # Transfor lists into tuples for each receptor
                            receptor_rars = tuple(receptor_rars)
                            receptor_sizes = tuple(receptor_sizes)
                            receptor_energies_av = tuple(receptor_energies_av)
                            receptor_energies_sd = tuple(receptor_energies_sd)
                            # Add the tuples to the general lists with the information of all receptors
                            rar_values += [receptor_rars]
                            sizes += [receptor_sizes]
                            energies_av += [receptor_energies_av]
                            energies_sd += [receptor_energies_sd]

                        # Generate the plot
                        plt.figure(
                            figsize=(24, 8)
                        )  # This size allows to have the labels and titles well separated, however it can
                        # be better to adjust it to the number of receptors to analyse.
                        for xe, ye, se, ce, le in zip(
                            receptors, rar_values, sizes, energies_av, energies_sd
                        ):
                            plt.scatter(
                                [xe] * len(ye),
                                ye,
                                s=se,
                                c=ce,
                                cmap="viridis_r",
                                edgecolors="black",
                                linewidth=le,
                                alpha=0.35,
                            )
                        plt.title(
                            f"RAR distribution {scenario}_{step}_{approach}_{rar_threshold}",
                            fontsize=20,
                        )
                        plt.ylabel("Right Answer Rate (%)", fontsize=15)
                        plt.ylim(-5, 105)
                        plt.yticks([0, 20, 40, 60, 80, 100])
                        plt.xticks(fontsize=12, rotation=65)
                        plt.colorbar()
                        plt.tight_layout()
                        # The place to save the file is bm_directory/analysis/scenario/step/scoring_approach/
                        plt.savefig(
                            f"{ana_approach_folder}/rar_{scenario}_{step}_{approach}_{rar_threshold}.png"
                        )
                        plt.close()
