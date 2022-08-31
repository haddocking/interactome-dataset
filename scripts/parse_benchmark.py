# This script identifies different mistakes and errors in a dataset:
#    - Sth wrong with the structure: two types of failure, a) Not a protein and b) Not correct structure
#    - No positive/negative ligands only positives/negatives: the receptor only has a type of ligands in the benchmark
#    - No interface API: there is no interface information for the molecule in the Uniprot API
#    - Positive = negative: the same ligand appears as positive and negative for the same receptor
#
# It is recommended to be used in the mode 'eval' in which it will generate a parse_benchmark.txt file in which
# it can be seen the analysis and evaluation of every molecule in the dataset. It is a way to check that the
# benchmark creation was successful.
#
# Additionally, it contains the mode 'clean' in which it will create a new benchmark avoiding all the mistakes.
# Note that the preference way must be first fixing the dataset generation, but it can be used as a complement
# in some cases.
#
# Execution: python parse_bm_linux.py path/to/benchmark mode
# E.g: $ python parse_bm_linux.py ./benchmarks/BM-230222 eval

import argparse
import shutil
import time
from pathlib import Path
from urllib.request import urlopen

from tabulate import tabulate

INTERFACE_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues"


def aminoacid_count(pdb_file):
    """Count the number of residues in a pdb file"""

    pdb_name = Path(pdb_file).stem.split("_")[0]
    chain = pdb_name[4:]
    pdb = open(pdb_file, "r").read()
    pdb_lines = pdb.splitlines()
    residues = []
    num_residues_pdb = 0

    for atom in pdb_lines:
        atom_features = atom.split()
        if atom_features[0] == "ATOM" and "CA" in atom_features:
            if len(atom_features[4]) == len(chain):
                if atom_features[4] == chain:
                    res = int(atom_features[5])
                    if res not in residues:
                        residues += [res]
            # When number of residue > 1000
            else:
                in_chain = atom_features[4][: len(chain)]
                res = int(atom_features[4][len(chain) :])
                if (
                    res > 5000
                ):  # There is a problem with the numbering, better to discard
                    return num_residues_pdb
                else:
                    if in_chain == chain:
                        if res not in residues:
                            residues += [res]

    num_residues_pdb = len(residues)
    return num_residues_pdb


def discontinuities_check(pdb_file):
    "Check if there are discontinuities in the pdb file"

    pdb_name = Path(pdb_file).stem.split("_")[0]
    chain = pdb_name[4:]
    pdb = open(pdb_file, "r").read()
    pdb_lines = pdb.splitlines()
    discontinuities = False

    # Check discontinuities
    for atom in range(0, len(pdb_lines) - 1):
        atom1_features = pdb_lines[atom].split()
        atom2_features = pdb_lines[atom + 1].split()
        if atom1_features[0] == "ATOM" and atom2_features[0] == "ATOM":
            if len(atom1_features[4]) == len(chain):
                res_1 = int(atom1_features[5])
            # When number of residue > 1000
            else:
                res_1 = int(atom1_features[4][len(chain) :])
            if len(atom2_features[4]) == len(chain):
                res_2 = int(atom2_features[5])
            # When number of residue > 1000
            else:
                res_2 = int(atom2_features[4][len(chain) :])

            if res_1 != res_2:
                if res_1 != res_2 - 1:
                    discontinuities = True
        elif atom1_features[0] == "TER" and atom2_features[0] == "ATOM":
            discontinuities = True

    return discontinuities


def evaluate(prot_id, prot_type, interaction_type, pdb_file):
    """Helper function to check the lenght of the pdb files and if exists an interface API for the uniprot ID"""
    if prot_type == "receptor":
        interaction_type = " "

    if (
        prot_id in parse_dict.keys()
    ):  # Avoid analyzing already done receptors as ligands.
        num_residues = parse_dict[prot_id][0]
        evaluation = parse_dict[prot_id][1]
        parse_table_line = [
            prot_id,
            prot_type,
            interaction_type,
            num_residues,
            evaluation,
        ]
        return parse_table_line

    else:
        num_residues = aminoacid_count(pdb_file)
        parse_dict[prot_id] = [num_residues]
        if num_residues == 0:  # 2nd mistake
            evaluation = "Sth wrong with structure"
            parse_dict[prot_id] += [evaluation]
            parse_table_line = [
                prot_id,
                prot_type,
                interaction_type,
                num_residues,
                evaluation,
            ]
            return parse_table_line

        else:
            discontinuities = discontinuities_check(pdb_file)
            if discontinuities:
                evaluation = "Discontinuities"
                parse_dict[prot_id] += [evaluation]
                parse_table_line = [
                    prot_id,
                    prot_type,
                    interaction_type,
                    num_residues,
                    evaluation,
                ]
                return parse_table_line

            else:  # Only if the receptor has passed the first two mistakes, it goes for the 3rd one.
                if str(prot_id[-1]).islower():
                    real_prot_id = prot_id[:-1]
                else:
                    real_prot_id = prot_id
                url = f"{INTERFACE_URL}/{real_prot_id}"
                try:
                    response = urlopen(url)
                    evaluation = " "
                    parse_dict[prot_id] += [evaluation]
                    parse_table_line = [
                        prot_id,
                        prot_type,
                        interaction_type,
                        num_residues,
                        evaluation,
                    ]
                    return parse_table_line
                except:
                    # Check 3 times because the server sometimes fail
                    time.sleep(0.1)
                    try:
                        response = urlopen(url)
                        evaluation = " "
                        parse_dict[prot_id] += [evaluation]
                        parse_table_line = [
                            prot_id,
                            prot_type,
                            interaction_type,
                            num_residues,
                            evaluation,
                        ]
                        return parse_table_line
                    except:
                        time.sleep(0.1)
                        try:
                            response = urlopen(url)
                            evaluation = " "
                            parse_dict[prot_id] += [evaluation]
                            parse_table_line = [
                                prot_id,
                                prot_type,
                                interaction_type,
                                num_residues,
                                evaluation,
                            ]
                            return parse_table_line
                        except:
                            evaluation = "No interface API"  # 3rd mistake
                            parse_dict[prot_id] += [evaluation]
                            # We need to add to output table here, otherwise the variable evaluation won't be stored when we continue
                            parse_table_line = [
                                prot_id,
                                prot_type,
                                interaction_type,
                                num_residues,
                                evaluation,
                            ]
                            return parse_table_line


def create_file(directory, output_table, mistakes):
    """Creates a text file with the results of the parsing"""
    table_address = f"{str(directory)}/parse-table.txt"
    parse_file = open(table_address, "w")
    parse_text = f"This is a txt file which contains mistakes in the parsing of the benchmark.\n\
There are different types of mistakes: \n\
    - Sth wrong with the structure: two types of failure, a) Not a protein and b) Not correct structure \n\
    - No positive/negative ligands only positives/negatives: the receptor only has a type of ligands in the benchmark \n\
    - No interface API: there is no interface information for the molecule in the Uniprot API \n\
    - Positive = negative: the same ligand appears as positive and negative for the same receptor \n\
The total number of mistakes is: {mistakes}\n"
    parse_file.write(parse_text)
    parse_file.write(tabulate(output_table))
    parse_file.close()
    return


if __name__ == "__main__":

    # Process the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", help="Introduce benchmark location")
    parser.add_argument(
        "mode",
        help="Introduce 'eval' for report of mistakes or 'clean' for generating a new fixed benchmark",
    )
    args = parser.parse_args()
    bm_directory = Path(args.benchmark)
    mode = str(args.mode)

    if mode != "eval" and mode != "clean":
        error_msg = "Incorrect use of the script, please introduce 'eval' or 'clean' as second argument. \n\
    E.g: python parse_clean_bm_linux.py BM-123456 eval"
        print(error_msg)
        exit()

    elif mode == "clean":  # Ask for the new benchmark name.
        clean_bm_name = str(input("Please introduce the name of the new benchmark: "))
        og_bm_parent_folder = bm_directory.parent
        if str(og_bm_parent_folder) == str(
            bm_directory
        ):  # The new benchmark is created in the og benchmark parent folder
            og_bm_parent_folder = Path("./..")
        clean_bm = f"{og_bm_parent_folder}/{clean_bm_name}"

    # Different libraries to store information.
    receptor_list = {}
    ligand_list = {}
    pair_list = []
    group_list = []

    partners_library = {}
    parse_dict = {}

    # Output table of the parsing.
    output_parse_table = [["Protein_ID", "Type", "Interaction", "Size", "Evaluation"]]

    # This step will allow us to detect receptors that contains only positive or only negative ligands.
    for folder in bm_directory.iterdir():
        group = folder.name  # group = "Receptor_Ligand_positive/negative"
        molecule_names = group.split("_")
        if (
            folder.is_dir() and len(molecule_names) == 3
        ):  # This allow us to skip files and other folders
            molecule_names = group.split("_")
            receptor = molecule_names[0]
            ligand = molecule_names[1]
            interaction = molecule_names[2]  # positive/negative

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

    # We loop each folder again to look for the different types of mistakes.
    for folder in bm_directory.iterdir():
        group = folder.name  # group = "Receptor_Ligand_positive/negative"
        molecule_names = group.split("_")
        if (
            folder.is_dir() and len(molecule_names) == 3
        ):  # This allow us to skip files and other folders.
            receptor = molecule_names[0]
            ligand = molecule_names[1]
            interaction = molecule_names[2]  # positive/negative
            pair = f"{receptor}_{ligand}"

            for file in folder.iterdir():
                file_name = file.name
                num_residues = (
                    "Not counted"  # We avoid to count if we first find another problem.
                )

                # Receptors. We need to go from top to bottom in the ranking:
                # 1. Only positive/negative ligands: it is the fastest so it is better to set it first.
                # 2. Structure: if sth is wrong here it is better not to continue.
                # 3. Interface.
                if file_name.endswith("r_b.pdb"):
                    if (
                        receptor not in receptor_list.keys()
                    ):  # Some receptors are also ligands for other receptors, so we use different lists.
                        receptor_list[receptor] = False
                        # Address the first mistake, which is the fastest to check and only applies to receptors.
                        if (
                            len(partners_library[receptor]["positive"]) == 0
                            and len(partners_library[receptor]["negative"]) != 0
                        ):
                            evaluation = "No positive ligands, only negatives"
                            parse_table_line = [
                                receptor,
                                "receptor",
                                " ",
                                num_residues,
                                evaluation,
                            ]
                            output_parse_table += [parse_table_line]
                        elif (
                            len(partners_library[receptor]["positive"]) != 0
                            and len(partners_library[receptor]["negative"]) == 0
                        ):
                            evaluation = "No negative ligands, only positives"
                            parse_table_line = [
                                receptor,
                                "receptor",
                                " ",
                                num_residues,
                                evaluation,
                            ]
                            output_parse_table += [parse_table_line]
                        # We don't add the evaluation to parse_dict because it only applies to receptors

                        else:
                            parse_table_line = evaluate(
                                receptor, "receptor", interaction, str(file)
                            )
                            output_parse_table += [parse_table_line]
                            if parse_table_line[4] == " ":
                                receptor_list[receptor] = True
                            continue

                # Ligands. We need to go from top to bottom in the ranking:
                # 1. Same ligand as positive/negative.
                # 2. Structure: if sth is wrong here it is better not to continue.
                # 3. Interface.
                elif file_name.endswith("l_b.pdb"):
                    # 1st mistake
                    if (
                        ligand in ligand_list.keys()
                        and pair in pair_list
                        and group not in group_list
                    ):
                        evaluation = "Positive = negative"
                        parse_table_line = [
                            ligand,
                            "ligand",
                            interaction,
                            num_residues,
                            evaluation,
                        ]
                        output_parse_table += [parse_table_line]
                        ligand_list[ligand] = False
                        # Don't add the evaluation to parse_dict because it only applies to ligands

                    elif ligand not in ligand_list.keys():
                        ligand_list[
                            ligand
                        ] = False  # ligand_list = ['ligand1', 'ligand2', 'ligand3', ...]
                        pair_list += [
                            pair
                        ]  # pair_list = ['receptor1_ligand1', 'receptor1_ligand2', 'receptor2_ligand3', ...]
                        group_list += [
                            group
                        ]  # group_list = ['receptor1_ligand1_positive', 'receptor1_ligand1_negative', ...]
                        parse_table_line = evaluate(
                            ligand, "ligand", interaction, str(file)
                        )
                        output_parse_table += [parse_table_line]
                        if (
                            parse_table_line[4] == " "
                            or parse_table_line[4] == "Discontinuities"
                        ):
                            ligand_list[ligand] = True
                        continue

            if mode == "clean":
                if receptor_list[receptor] and ligand_list[ligand]:
                    new_bm_folder = f"{str(clean_bm)}/{group}"
                    shutil.copytree(folder, new_bm_folder, dirs_exist_ok=True)

    n_mistakes = -1
    for line in output_parse_table:
        if (line)[4] != " ":
            n_mistakes += 1
    print(("Assessment finished with {} mistakes found.").format(n_mistakes))

    # Creating the output file.
    create_file(bm_directory, output_parse_table, n_mistakes)
