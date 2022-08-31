# This script creates a folder organized benchmark from the data obtained from create-db.py and
# stored in interacting_db.csv.
#
# The generated benchmark has the following structure:
# Benchmark/
#   |_ ReceptorID_LigandID_positive/
#           |_ pdbIDchain_r_b.pdb               (receptor pdb file)
#           |_ pdbIDchain_l_b_positive.pdb      (positive ligand pdb file)
#   |_ ReceptorID_LigandID_negative/
#           |_ pdbIDchain_r_b.pdb               (receptor pdb file)
#           |_ pdbIDchain_l_b_negative.pdb      (negative ligand pdb file)
#   |_ ...
#
# It can be run with different filters:
#    - outHD: Exclude homodimers
#    - outIg: Exclude Immunoglobulins
#    - inIg: Take only Immunoglobulins
#    - outCA: Exclude molecules whose pdb structure is only Alpha Carbons
#    - outMS: Exclude multiple structures for the same molecule. It takes the one with more number of partners and if tie, the one with more residues.
#    - AO: Exclude homodimers, Immunoglobulins and AlphaCarbons

# The script also stores information of Igs from previous runs, to make it faster to run,
# if wanted to reprocess from the beginning, use the flag -abini (ab initio).
#
# Execution: python generate_benchmark.py [-outHD] [-outIg] [-inIg] [-outCA] [-outMS] [-AO] [-abini] location/of/interacting_db.csv location/of/create_benchmark_folder location/of/the/benchmark [location/of/interactome_db.txt]
# e.g: $ python generate_benchmark.py -outHD -outIg ./create_benchmark/interacting_db.csv ./create_benchmark ./BMAO-230222 .create_benchmark/interactome_db.txt
#
#    MSc. Jesús López Rivera, Bonvin Lab :)

import argparse
import ast
import json
import logging
import os
import shutil
import time
from pathlib import Path

import requests

FORMAT = " %(asctime)s L%(lineno)d %(levelname)s - %(message)s"
logging.basicConfig(format=FORMAT, level="INFO")

PDB_MOLECULES_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/"
INTERFACE_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues"
PROTEINS_URL = "https://www.ebi.ac.uk/proteins/api/proteins"
QUERY_DB = {}

Ig_words = [
    "Ig",
    "IG",
    "Immunoglobulin",
    "immunoglobulin",
    "IMMUNOGLOBULIN",
    "Antibody",
    "antibody",
    "ANTIBODY",
    "SAB",
]
# Doesn't include "Ab" or "ab" due to some exceptions like "Rab" proteins or "probable" in name.


def make_request(url, data):
    """Helper function to make the requests."""
    # The requests might be done multiple times,
    # use QUERY_DB dictionary to store the results
    if url not in QUERY_DB:
        # first time this query is done
        response = requests.get(url)

        if response.status_code == 404:
            # not found == empty or failure of seraver, retry 2 times
            time.sleep(0.1)
            response = requests.get(url)
            if response.status_code == 404:
                time.sleep(0.1)
                response = requests.get(url)
                if response.status_code == 404:
                    data = None
                else:
                    data = response.json()
            else:
                data = response.json()
        else:
            data = response.json()

        QUERY_DB[url] = data
    else:
        # Done before
        data = QUERY_DB[url]

    # Wait here not to get banned out of the API
    time.sleep(0.1)
    return data


def check_if_HD(protein_id):
    url = f"{INTERFACE_URL}/{protein_id}"
    homodimer = False

    # The request might fail
    try:
        interface_data = make_request(url, None)
    except Exception as e:
        logging.debug(f"Could not make InterfaceResidues request for {protein_id}, {e}")
        pass

    if interface_data and len(interface_data) != 0:
        if len(interface_data[protein_id]) == 0:  # Check if the entry is empty.
            logging.debug(f"InterfaceResidues for {protein_id} was empty")
        else:
            # Iterate every interacting partner
            for partner in interface_data[protein_id]["data"]:
                partner_id = partner["accession"]
                if partner_id == protein_id:
                    homodimer = True
                    break
    return homodimer


def check_if_Ig(pdb_id, uniprot_ID, Ig_words):
    "Use the PDB API to find the real chain ID."

    url = f"{PDB_MOLECULES_URL}/{pdb_id}"
    url_name = f"{PROTEINS_URL}/{uniprot_ID}"
    Ig = False

    try:
        pdb_data = make_request(url, None)
    except Exception as e:
        logging.info(f"Could not make pdb/molecules API request for {pdb_id}, {e}")
        return None

    if pdb_data:
        try:
            protein_data = make_request(url_name, None)
        except Exception as e:
            logging.info(
                f"Could not make InterfaceResidues request for {uniprot_ID}, {e}"
            )
            return False

        if protein_data:
            if len(protein_data) == 0:  # Check if the entry is empty.
                logging.info(f"Proteins API for {uniprot_ID} was empty")
                return False
            else:
                all_names = []
                if "recommendedName" in protein_data["protein"].keys():
                    try:
                        main_name = protein_data["protein"]["recommendedName"][
                            "fullName"
                        ]["value"]
                        all_names += [main_name]
                    except:
                        pass
                if "alternativeName" in protein_data["protein"].keys():
                    try:
                        for alt_name in protein_data["protein"]["alternativeName"]:
                            all_names += [alt_name["fullName"]["value"]]
                    except:
                        pass

                if len(all_names) == 0:
                    logging.info(
                        f"Protein API for {uniprot_ID} did NOT provide a valid name"
                    )
                    return None

        if len(pdb_data[pdb_id]) == 0:  # Check if the entry is empty.
            logging.info(f"PDB_Moleucles API for {pdb_id} was empty")
            return None
        else:
            for entity in pdb_data[pdb_id]:
                if entity["molecule_name"] in all_names:
                    for ig_word in Ig_words:
                        if "synonym" in entity.keys():
                            if ig_word in str(entity["synonym"]) or ig_word in str(
                                entity["molecule_name"]
                            ):
                                Ig = True
                                break
                        else:
                            if ig_word in str(entity["molecule_name"]):
                                Ig = True
                                break
    return Ig


def check_if_CA(pdb_file):
    "Check if the structure is an only alpha carbons"

    CA = False
    pdb = open(pdb_file, "r").read()
    pdb_lines = pdb.splitlines()
    atom_types = []

    for atom in pdb_lines:
        atom_features = atom.split()
        if atom_features[0] == "ATOM":
            atom_type = atom_features[2]
            if atom_type not in atom_types:
                atom_types += [atom_type]

    if len(atom_types) == 1:
        CA = True

    return CA


def sort_by_positives(receptor):
    rec_pos_dic = {}
    n_pos = []
    sorted_pos_pdbs = []
    for receptor_pdb in receptor.keys():
        n_pos_rec_pdb = len(receptor[receptor_pdb]["positives"])
        if n_pos_rec_pdb not in rec_pos_dic.keys():
            n_pos += [n_pos_rec_pdb]
            rec_pos_dic[n_pos_rec_pdb] = []
        rec_pos_dic[n_pos_rec_pdb] += [receptor_pdb]
    sorted_n_pos = sorted(n_pos, reverse=True)
    for s_n_pos in sorted_n_pos:
        rec_pdbs = rec_pos_dic[s_n_pos]
        if len(rec_pdbs) > 1:
            logging.info("Sorting receptor by structure")
            rec_pdbs = sort_by_structure(receptor, rec_pdbs)
        sorted_pos_pdbs += rec_pdbs

    return sorted_pos_pdbs


def sort_by_structure(protein, pdbs):
    if pdbs == None:
        pdbs = protein.keys()

    protein_len_dic = {}
    lengths = []
    sorted_len_pdbs = []
    for pdb in pdbs:
        length = len(protein[pdb]["residues"])
        if length not in protein_len_dic.keys():
            lengths += [length]
            protein_len_dic[length] = []
        protein_len_dic[length] += [pdb]
    sorted_lengths = sorted(lengths, reverse=True)
    for s_length in sorted_lengths:
        sorted_len_pdbs += protein_len_dic[s_length]
    return sorted_len_pdbs


def write_datafile(data, address):
    """Helper function to write a file"""
    groups_table = open(address, "w")
    groups_table.write(str(data))
    groups_table.close()
    return


def write_dic(uniprot_id, pdb_file, uniprot_pdb_dic, pdb_uniprot_dic):
    """Write two dictionaries, assigning the correct pdbID to the correct UniprotID and vicerversa"""
    # Need to be careful because one single UniprotID may have different pdbID files. To distinguish
    # them a lowercase letter will be added to the UniprotID, from b to z.
    for letter in [
        "b",
        "c",
        "d",
        "e",
        "f",
        "g",
        "h",
        "i",
        "j",
        "k",
        "l",
        "m",
        "n",
        "o",
        "p",
        "q",
        "r",
        "s",
        "t",
        "u",
        "v",
        "w",
        "x",
        "y",
        "z",
    ]:
        if uniprot_id not in uniprot_pdb_dic:
            uniprot_pdb_dic[uniprot_id] = pdb_file
            pdb_uniprot_dic[pdb_file] = uniprot_id
            break
        else:
            if pdb_file != uniprot_pdb_dic[uniprot_id]:
                if str(uniprot_id[-1]).islower():
                    uniprot_id = uniprot_id[:-1]
                uniprot_id = str(uniprot_id) + letter

    return uniprot_pdb_dic, pdb_uniprot_dic


def mkdir_cpfiles(new_folder, receptor_pdb_file, ligand_pdb_file, interaction):
    """Make a directory in the new benchmark and copy the file from the storage folder"""
    try:
        new_folder.mkdir(
            exist_ok=False
        )  # Just to make sure there is nothing strange in the data
    except:
        new_folder_name = Path(new_folder).stem
        logging.warning(f"The folder {new_folder_name} is a duplicate")
    receptor_pdb_name = (
        f"{(Path(receptor_pdb_file).stem)[:4]}{(Path(receptor_pdb_file).stem)[5]}"
    )
    receptor_pdb_file_address = f"{create_benchmark_directory}/{receptor_pdb_file}"
    new_receptor_pdb = f"{new_folder}/{receptor_pdb_name}_r_b.pdb"
    shutil.copy(receptor_pdb_file_address, new_receptor_pdb)

    ligand_pdb_name = (
        f"{(Path(ligand_pdb_file).stem)[:4]}{(Path(ligand_pdb_file).stem)[5]}"
    )
    ligand_pdb_file_address = f"{create_benchmark_directory}/{ligand_pdb_file}"
    new_ligand_pdb = f"{new_folder}/{ligand_pdb_name}_{interaction}_l_b.pdb"
    shutil.copy(ligand_pdb_file_address, new_ligand_pdb)
    return


if __name__ == "__main__":

    # Process the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", help="Introduce 'interacting_db.csv' file")
    parser.add_argument(
        "create_benchmark", help="Introduce location of the create_benchmark folder"
    )
    parser.add_argument(
        "benchmark",
        help="Introduce the new benchmark path, e.g: /path/to/new/benchmark/BM-123456",
    )
    parser.add_argument(
        "-outIg", action="store_true", help="Exclude Igs from the dataset"
    )
    parser.add_argument(
        "-inIg", action="store_true", help="Include only Igs in the benchmark"
    )
    parser.add_argument(
        "-outHD", action="store_true", help="Exclude homodimers from the dataset"
    )
    parser.add_argument(
        "-outCA",
        action="store_true",
        help="Exclude structures that only have Alpha Carbons (CA)",
    )
    parser.add_argument(
        "-outMS",
        action="store_true",
        help="Exclude multiple structures for a same pair receptor-ligand",
    )
    parser.add_argument(
        "-AO", action="store_true", help="All Out, exclude Igs, HDs, CAs and MS"
    )
    parser.add_argument(
        "-abini",
        action="store_true",
        help="ab initio: start from the beginning not reusing anything",
    )
    parser.add_argument(
        "interactome_db",
        nargs="?",
        help="Introduce interactome_db.txt file if possible",
    )
    args = parser.parse_args()
    bm_directory = Path(args.benchmark)
    if not bm_directory.exists():
        bm_directory.mkdir()

    csv_file_address = Path(args.csv_file)
    create_benchmark_directory = Path(
        args.create_benchmark
    )  # Where create-db.py stores everything.

    outIg = args.outIg
    inIg = args.inIg
    outHD = args.outHD
    outCA = args.outCA
    outMS = args.outMS
    AO = args.AO
    abini = args.abini

    if AO:
        outIg = outHD = outMS = outCA = True
        logging.info(f"Modes outHD, outIg, outCA and outHM selected")

    if outMS or outCA:
        interactome_db_file = args.interactome_db
        if interactome_db_file is None:
            interactome_db_file = f"{csv_file_address.parent}/interactome_db.txt"

        if Path(interactome_db_file).exists():
            logging.info(">>>>>> interactome_db.txt file found, reusing it. >>>>>>")
            interactome_read = open(interactome_db_file, "r").read()
            interactome_db = ast.literal_eval(interactome_read)
            default_MS = False
        else:
            logging.warning(
                "<<<<<< No interactome_db.txt file found, proceding to default mode. Please consider providing one. <<<<<<"
            )
            default_MS = True

    # Different libraries to store information
    groups = []
    uniprot_pdb_dic = {}
    pdb_uniprot_dic = {}
    dataset = {}
    MS = []

    # Open the csv file
    csv_file = open(csv_file_address, "r").read()
    csv_lines = csv_file.splitlines()

    if outHD:
        """Exclude homodimers"""

        HDs = {}
        HDs_address = f"{create_benchmark_directory}/homodimers.txt"
        not_HD_lines = ["Not Homodimers"]
        n_HDs = 0

        """Storing data on homodimers can be useful for future analysis, so the molecule IDs of those 
        is stored later in a file called homodimers.txt that can also be reused to run faster"""
        if Path(HDs_address).exists():
            if not abini:
                # Reuse old processed data
                logging.info("homodimers.txt found, reusing it.")
                read_HDs = open(HDs_address, "r").read()
                HDs = ast.literal_eval(read_HDs)
            else:
                # Discard old data
                os.remove(HDs_address)

        # Iterate the csv file
        for line in csv_lines[1:]:
            line_data = line.split(",")

            receptor_id = line_data[0]
            ligand_id = line_data[2]
            interaction = line_data[4]

            if receptor_id not in HDs.keys():
                receptor_is_HD = check_if_HD(receptor_id)
                HDs[receptor_id] = receptor_is_HD
            else:
                receptor_is_HD = HDs[receptor_id]

            if ligand_id not in HDs.keys():
                ligand_is_HD = check_if_HD(ligand_id)
                HDs[ligand_id] = ligand_is_HD
            else:
                ligand_is_HD = HDs[ligand_id]

            if receptor_is_HD or ligand_is_HD:
                n_HDs += 1
                if receptor_is_HD:
                    logging.info(f"{receptor_id} is a Homodimer.")
                elif ligand_is_HD:
                    logging.info(f"{ligand_id} is a Homodimer.")

            elif not receptor_is_HD and not ligand_is_HD:
                not_HD_lines += [line]

        write_datafile(HDs, HDs_address)  # Write a text file with the data in Igs.

        csv_lines = not_HD_lines
        print(f"The number of homodimers in the dataset is {n_HDs}")

    if outIg or inIg:

        Igs_address = f"{create_benchmark_directory}/immunoglobulins.txt"
        Igs_dict = {}

        """Storing data on immunoglobulins can be useful for future analysis, so the molecule IDs of those 
        is stored later in a file called immunoglobulins.txt that can also be reused to run faster"""
        if Path(Igs_address).exists():
            if not abini:
                # Reuse old processed data
                logging.info("immunoglobulins.txt found, reusing it.")
                read_Igs = open(Igs_address, "r").read()
                Igs_dict = ast.literal_eval(read_Igs)
            else:
                # Discard old data
                os.remove(Igs_address)

        not_Ig_lines = ["Not Igs"]
        Ig_lines = ["Igs"]

        # Iterate the csv file
        for line in csv_lines[1:]:
            line_data = line.split(",")

            receptor_id = line_data[0]
            receptor_pdb = line_data[1]
            receptor_pdb_name = str(Path(receptor_pdb).stem)
            ligand_id = line_data[2]
            ligand_pdb = line_data[3]
            ligand_pdb_name = str(Path(ligand_pdb).stem)
            interaction = line_data[4]

            if receptor_id not in Igs_dict:
                receptor_is_Ig = check_if_Ig(
                    receptor_pdb_name[:4], receptor_id, Ig_words
                )  # check_if_Ig(pdbID, chain, Ig_words)
                Igs_dict[receptor_id] = receptor_is_Ig
            else:
                receptor_is_Ig = Igs_dict[receptor_id]

            if ligand_id not in Igs_dict:
                ligand_is_Ig = check_if_Ig(
                    ligand_pdb_name[:4], ligand_id, Ig_words
                )  # check_if_Ig(pdbID, chain, Ig_words)
                Igs_dict[ligand_id] = ligand_is_Ig
            else:
                ligand_is_Ig = Igs_dict[ligand_id]

            if receptor_is_Ig or ligand_is_Ig:
                if receptor_is_Ig:
                    logging.info(f"{receptor_id} is an Immunoglobulin.")
                elif ligand_is_Ig:
                    logging.info(f"{ligand_id} is an Immunoglobulin.")

                Ig_lines += [line]

            elif not receptor_is_Ig and not ligand_is_Ig:
                not_Ig_lines += [line]

        write_datafile(Igs_dict, Igs_address)  # Write a text file with the data in Igs.

        if outIg:
            csv_lines = not_Ig_lines
        elif inIg:
            csv_lines = Ig_lines

        n_Igs = len(Ig_lines) - 1
        print(f"The number of Igs interactions in the dataset is {n_Igs}")

    if outCA:
        """Exclude only CA"""

        not_CA_lines = ["Not CA"]
        CAs = {}
        CAs_address = f"{create_benchmark_directory}/onlyCA.txt"
        n_CAs = 0

        """Storing data on alpha carbons can be useful for future analysis, so the molecule IDs of those 
        is stored later in a file called onlyCA.txt that can also be reused to run faster"""
        if Path(CAs_address).exists():
            if not abini:
                # Reuse old processed data
                logging.info("onlyCA.txt found, reusing it.")
                read_CAs = open(CAs_address, "r").read()
                CAs_dict = ast.literal_eval(read_CAs)
            else:
                # Discard old data
                os.remove(CAs_address)

        # Iterate the csv file
        for line in csv_lines[1:]:
            line_data = line.split(",")

            receptor_id = line_data[0]
            receptor_pdb = line_data[1]
            receptor_pdb_name = Path(receptor_pdb).stem
            ligand_id = line_data[2]
            ligand_pdb = line_data[3]
            ligand_pdb_name = Path(ligand_pdb).stem
            interaction = line_data[4]

            if receptor_pdb not in CAs.keys():
                receptor_pdb_path = f"{create_benchmark_directory}/{receptor_pdb}"
                receptor_is_CA = check_if_CA(receptor_pdb_path)
                CAs[receptor_pdb] = receptor_is_CA
            else:
                receptor_is_CA = CAs[receptor_pdb]

            if ligand_pdb not in CAs.keys():
                ligand_pdb_path = f"{create_benchmark_directory}/{ligand_pdb}"
                ligand_is_CA = check_if_CA(ligand_pdb_path)
                CAs[ligand_pdb] = ligand_is_CA
            else:
                ligand_is_CA = CAs[ligand_pdb]

            if receptor_is_CA or ligand_is_CA:
                n_CAs += 1
                if receptor_is_CA:
                    logging.info(
                        f"Receptor {receptor_id} - {receptor_pdb_name} is only alpha carbons."
                    )
                    if receptor_id in interactome_db.keys():
                        if receptor_pdb_name in interactome_db[receptor_id].keys():
                            del interactome_db[receptor_id][receptor_pdb_name]
                elif ligand_is_CA:
                    logging.info(
                        f"Ligand {ligand_id} - {ligand_pdb_name} is only alpha carbons."
                    )
                    if receptor_id in interactome_db.keys():
                        if receptor_pdb_name in interactome_db[receptor_id].keys():
                            if (
                                ligand_id
                                in interactome_db[receptor_id][receptor_pdb_name][
                                    f"{interaction}s"
                                ].keys()
                            ):
                                if (
                                    ligand_pdb_name
                                    in interactome_db[receptor_id][receptor_pdb_name][
                                        f"{interaction}s"
                                    ][ligand_id].keys()
                                ):
                                    del interactome_db[receptor_id][receptor_pdb_name][
                                        f"{interaction}s"
                                    ][ligand_id][ligand_pdb_name]

            elif not receptor_is_CA and not ligand_is_CA:
                not_CA_lines += [line]

        write_datafile(CAs, CAs_address)  # Write a text file with the data in Igs.

        csv_lines = not_CA_lines
        print(f"The number of only CA structures in the dataset is {n_CAs}")

    if outMS and not default_MS:
        outMS_lines = ["outMS"]
        # First selection: higher number of positive ligands
        # Second selection higher lenght of structure
        for receptor in interactome_db.keys():
            all_positives = []
            best_negatives = []
            for receptor_pdb in interactome_db[receptor].keys():
                for positive_partner in interactome_db[receptor][receptor_pdb][
                    "positives"
                ].keys():
                    if positive_partner not in all_positives:
                        all_positives += [positive_partner]

            logging.info(f"Sorting receptor {receptor} by positives and/or structure")
            sort_by_pos_rec_pdbs = sort_by_positives(interactome_db[receptor])

            while len(all_positives) != 0:
                for rec_pdb in sort_by_pos_rec_pdbs:
                    best_ligands = []
                    useful_structure = False
                    if len(best_negatives) == 0:
                        for negative in interactome_db[receptor][rec_pdb][
                            "negatives"
                        ].keys():
                            # logging.info("Sorting negatives by structure")
                            sort_by_str_neg_pdbs = sort_by_structure(
                                interactome_db[receptor][rec_pdb]["negatives"][
                                    negative
                                ],
                                None,
                            )
                            best_negatives += [sort_by_str_neg_pdbs[0]]
                    best_ligands += best_negatives
                    for positive in interactome_db[receptor][rec_pdb][
                        "positives"
                    ].keys():
                        if positive in all_positives:
                            all_positives.remove(positive)
                            useful_structure = True
                            # logging.info("Sorting positives by structure")
                            sort_by_str_pos_pdbs = sort_by_structure(
                                interactome_db[receptor][rec_pdb]["positives"][
                                    positive
                                ],
                                None,
                            )
                            best_ligands += [sort_by_str_pos_pdbs[0]]

                    if useful_structure:
                        for line in csv_lines[1:]:
                            line_data = line.split(",")
                            receptor_id = line_data[0]
                            receptor_pdb = Path(line_data[1]).stem
                            ligand_id = line_data[2]
                            ligand_pdb = Path(line_data[3]).stem
                            interaction = line_data[4]
                            if receptor_id == receptor and receptor_pdb == rec_pdb:
                                if ligand_pdb in best_ligands:
                                    outMS_lines += [line]

        csv_lines = outMS_lines

    for line in csv_lines[1:]:
        line_data = line.split(",")

        receptor_id = line_data[0]
        receptor_pdb = line_data[1]
        ligand_id = line_data[2]
        ligand_pdb = line_data[3]
        interaction = line_data[4]

        # Avoid multiple structures for the same pair receptor - ligand
        if outMS and default_MS:
            pair = f"{receptor_id}-{ligand_id}-{interaction}"
            if pair in MS:
                continue
            else:
                MS += [pair]

        # Process UniprotID and pdbID of the receptor
        uniprot_pdb_dic, pdb_uniprot_dic = write_dic(
            receptor_id, receptor_pdb, uniprot_pdb_dic, pdb_uniprot_dic
        )
        real_receptor_id = pdb_uniprot_dic[receptor_pdb]
        # Write a dictionary containing all positive and all negative ligands for each receptor
        if real_receptor_id not in dataset.keys():
            dataset[real_receptor_id] = {"positive": [], "negative": []}
        # Process UniprotID and pdbID of the ligand
        uniprot_pdb_dic, pdb_uniprot_dic = write_dic(
            ligand_id, ligand_pdb, uniprot_pdb_dic, pdb_uniprot_dic
        )
        real_ligand_id = pdb_uniprot_dic[ligand_pdb]
        # Add the ligand to the corresponding interaction of its receptor.
        if real_ligand_id not in dataset[real_receptor_id][interaction]:
            dataset[real_receptor_id][interaction] += [real_ligand_id]
        else:
            logging.warning(
                f"{real_receptor_id}_{real_ligand_id}_{interaction} is a multiple structure"
            )

    # Generate the benchmark
    for receptor in dataset.keys():
        if (
            len(dataset[receptor]["positive"]) > 0
            and len(dataset[receptor]["negative"]) > 0
        ):
            receptor_pdb_file = uniprot_pdb_dic[receptor]
            for interaction in dataset[receptor].keys():
                for ligand in dataset[receptor][interaction]:
                    ligand_pdb_file = uniprot_pdb_dic[ligand]
                    new_folder = Path(
                        f"{str(bm_directory)}/{receptor}_{ligand}_{interaction}"
                    )
                    if not new_folder.exists():
                        mkdir_cpfiles(
                            new_folder, receptor_pdb_file, ligand_pdb_file, interaction
                        )
                    else:
                        logging.warning(f"{new_folder.name} is a duplicate")

        # Check for mistakes
        elif (
            len(dataset[receptor]["positive"]) > 0
            and len(dataset[receptor]["negative"]) == 0
        ):
            logging.warning(f"{receptor} has positive ligands, but no negative ones")
        elif (
            len(dataset[receptor]["positive"]) == 0
            and len(dataset[receptor]["negative"]) > 0
        ):
            logging.warning(f"{receptor} has negative ligands, but no positive ones")

    print("Dataset creation finished")
