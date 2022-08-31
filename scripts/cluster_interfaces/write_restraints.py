# This scripts creates a clustered benchmark from the clustered interfaces obtained after running cluster_interfaces.py and
# cluster_distance_matrix.py. It gets the interfaces from an interface_dictionary.txt file and the different clustered
# ligands from the output file of cluster_distance_matrix.py

# It creates a restraints .tbl file combining the active residues of the interface at the receptors and the passive residues
# from the execution of calc_accessibility.py at the ligands PDBs.

# Execution: python write_restraints.py location/of/the/benchmark location/of/interface_dictionary location/of/clusters_file name_of_new_benchmark
# e.g. $ python ../scripts/write_restraints.py ./OutHM-220622/BM-AO-220622 ./OutHM-220622/BM-AO-220622/interface_analysis/interface_dictionary.txt
# ./OutHM-220622/BM-AO-220622/cluster_interfaces/clusters_single_thr-0.7071.out NESTOR

# MSc. Jesús López Rivera, Bonvin Lab :)

import argparse
import ast
import logging
import os
import shutil
import subprocess
from pathlib import Path

FORMAT = " %(asctime)s L%(lineno)d %(levelname)s - %(message)s"
logging.basicConfig(format=FORMAT, level="INFO")

CALC_ACCESSIBILITY = (
    "/trinity/login/jlopez/repos/haddock-CSB-tools/calc-accessibility.py"
)


def get_accessibility(pdb_file):
    """Helper function to run calc_accessibility.py and get solvent accessible residues of a pdb"""

    cmd = f"python {CALC_ACCESSIBILITY} {pdb_file}"
    p = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )

    out, _ = p.communicate()
    try:
        accessible_residues = (
            out.decode("utf-8").splitlines()[-1].split()[-1].split(",")
        )
    except Exception as e:
        logging.debug(f"{str(pdb_file)}, {e}")
        accessible_residues = None
    return accessible_residues


def contacts2file(file, protein_id, cluster):
    """Write contacts .txt files from a contact dictionary and below a certain appearance threshold"""

    with open(file, "w") as fh:
        if str(file).endswith("interface.txt"):
            for residue in interface_residues[protein_id][cluster]:
                fh.write(f"{str(residue)} ")
        elif str(file).endswith("accessible.txt"):
            for residue in accessibility_dic[protein_id]:
                fh.write(f"{str(residue)} ")
    fh.close()
    read_file = open(file, "r").read()
    if len(read_file) == 0:
        os.remove(file)
    return


def load_file(residue_file):
    """Read a space-separated file and return a list of integers"""

    output_list = []
    with open(residue_file, "r") as fh:
        for line in fh.readlines():
            # all the residues are in one single line
            residue_str_list = line.split()  # remove spaces
            residue_int_list = map(int, residue_str_list)  # convert to integers
            # add it to the residue_list
            for res in residue_int_list:
                output_list.append(res)
    return output_list


def generate_restraint(interface_r, interface_l, tbl_file, segid_r="A", segid_l="B"):
    """Write a .tbl restraints file from the two given .txt files"""

    tbl_string = ""
    for resi_r in interface_r:
        tbl_string += (
            "assign (resi {} and segid {})".format(resi_r, segid_r) + os.linesep
        )
        tbl_string += "(" + os.linesep
        c = 0
        for resi_l in interface_l:
            tbl_string += (
                "       (resi {} and segid {})".format(resi_l, segid_l) + os.linesep
            )
            c += 1
            if c != len(interface_l):
                tbl_string += "        or" + os.linesep
        tbl_string += ") 2.0 2.0 0.0" + os.linesep

    with open(tbl_file, "w") as out_fh:
        out_fh.write(tbl_string)
    out_fh.close()

    return


if __name__ == "__main__":

    # Process the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", help="Introduce benchmark location")
    parser.add_argument(
        "interface_dictionary", help="Please introduce interface_dictionary.txt file"
    )
    parser.add_argument("clusters_file", help="Please introduce clusters file")
    parser.add_argument(
        "new_benchmark_name", help="Please introduce name for the new benchmark"
    )

    args = parser.parse_args()
    bm_directory = Path(args.benchmark)
    interface_dictionary_file = args.interface_dictionary
    clusters_file = args.clusters_file
    new_bm_name = args.new_benchmark_name

    # Check if interface_dictionary.txt file exists and process it
    if Path(interface_dictionary_file).exists():
        logging.info("<<<<<< interface_dictionary.txt file found, reusing it. >>>>>>")
        read_interface_dictionary = open(interface_dictionary_file, "r").read()
        interface_dictionary = ast.literal_eval(read_interface_dictionary)
    else:
        logging.info("<<<<<< interface_dictionary.txt file not found. >>>>>>")
        exit()

    # Check if clusters file exists and process it
    if Path(clusters_file).exists():
        logging.info("<<<<<< clusters file found, reusing it. >>>>>>")
        clusters = open(clusters_file, "r").read()
    else:
        logging.info("<<<<<< clusters file not found. >>>>>>")
        exit()

    new_bm = f"{bm_directory.parent}/{new_bm_name}"  # Path to new clustered benchmark
    os.mkdir(new_bm)

    # Dictionaries to store information
    clusters_dic = {}
    interface_residues = {}
    accessibility_dic = {}

    """ Process clusters file and get the name of the ligands for each cluster """
    for line in clusters.splitlines():
        words = line.split()
        if words[0] == "Receptor":
            receptor = words[1]
            if receptor not in clusters_dic.keys():
                clusters_dic[receptor] = []
        elif words[0] == "Cluster":
            rec_cluster = words[3:]
            clusters_dic[receptor] += [rec_cluster]

    # clusters_dic = {
    #   Receptor1: [[Ligand1,Ligand2], [Ligand3,Ligand 4], ...],
    #   Receptor2: [...],
    #   ...}

    """ If receptor only has 1 ligand, it has not been clustered, so it is needed to add it to clusters_dic """
    for receptor in interface_dictionary.keys():
        if receptor not in clusters_dic.keys():
            clusters_dic[receptor] = [interface_dictionary[receptor].keys()]

    """ Get the interface residues for each ligand in the cluster """
    # Iterate every receptor
    for receptor in clusters_dic.keys():
        interface_residues[receptor] = []
        # Iterate every cluster for the receptor
        for cluster in clusters_dic[receptor]:
            interface_cluster = []
            # Iterate every ligand in the cluster
            for ligand in cluster:
                # Iterate every interface residue for the receptor-ligand interaction
                for residue in interface_dictionary[receptor][ligand]:
                    if residue not in interface_cluster:  # Avoid repeated ones
                        interface_cluster += [
                            residue
                        ]  # Add residue to the list of residues
            interface_cluster = sorted(
                interface_cluster
            )  # Sort them for a better visualization
            interface_residues[receptor] += [interface_cluster]  # Cluster of residues

    # Process every folder in the benchmark
    for folder in bm_directory.iterdir():
        folder_name = folder.stem
        molecules = folder_name.split("_")  # [Receptor, Ligand, "positive"/"negative"]
        if len(molecules) == 3 and folder.is_dir():  # Avoid files and other folders

            logging.info(f"Writing restraints files for folder {folder_name}")
            receptor = molecules[0]
            ligand = molecules[1]

            if (
                receptor in interface_residues.keys()
            ):  # The receptor is processed before

                """Find the PDB of the ligand and get the solvent accessible ligands of it"""
                for file in folder.iterdir():
                    if str(file).endswith("l_b.pdb"):
                        if (
                            ligand not in accessibility_dic.keys()
                        ):  # Not processed before
                            # Get the accessible residues
                            accessible_residues = get_accessibility(file)
                            accessibility_dic[ligand] = accessible_residues
                        else:
                            accessible_residues = accessibility_dic[ligand]

                if (
                    accessible_residues is None
                ):  # Something did not work with get accessibility
                    logging.warning(f"No accessible residues for {ligand}, continue")
                    continue
            else:
                # Something did not work with the receptor
                logging.warning(f"Receptor {receptor} not clustered, continue")
                continue

            # Store the information
            interface_folder = f"{str(folder)}/interface"
            if Path(interface_folder).exists():  # In case of already done
                shutil.rmtree(interface_folder)
            os.mkdir(interface_folder)

            """ Make a file for the accessible residues """
            contacts_file_l = f"{interface_folder}/{ligand}_l_accessible.txt"
            contacts2file(contacts_file_l, ligand, cluster=None)

            """ Make a .txt file for every cluster """
            n_clusters = len(clusters_dic[receptor])
            for cluster in range(0, n_clusters):
                str_cluster = f"cluster-{str(cluster+1)}"  # cluster-1 / cluster-2 / ...

                contacts_file_r = (
                    f"{interface_folder}/{receptor}_{str_cluster}_r_interface.txt"
                )
                contacts2file(contacts_file_r, receptor, cluster)

                """ Make a .tbl restraints file for every cluster """
                tbl_file = f"{interface_folder}/{str_cluster}_interface_restraints.tbl"
                # Check if both .txt files exist before doing the .tbl file
                if Path(contacts_file_r).exists() and Path(contacts_file_l).exists():
                    interface_r = load_file(contacts_file_r)
                    interface_l = load_file(contacts_file_l)
                    # Write the .tbl file
                    generate_restraint(interface_r, interface_l, tbl_file)

                else:  # One of them is missing
                    logging.warning(
                        f"One of the restraints file is missing for the folder {folder_name}"
                    )
                    continue

                # Make new folder for the cluster
                new_folder = f"{new_bm}/{folder_name}-{str_cluster}"
                os.mkdir(new_folder)

                """ Copy everything from the original benchmark to the clustered one """
                new_interface_folder = f"{str(new_folder)}/interface"
                os.mkdir(new_interface_folder)
                for item in folder.iterdir():
                    if str(item).endswith(".pdb"):  # PDB files
                        old_pdb_name = item.name
                        new_pdb_name = Path(f"{str(new_folder)}/{old_pdb_name}")
                        shutil.copyfile(item, new_pdb_name)
                    elif str(item).endswith("interface"):  # interface folder
                        for inter_file in item.iterdir():
                            old_inter_file_name = inter_file.name
                            if str_cluster in str(
                                inter_file
                            ):  # files that have cluster in name (.txt and .tbl)
                                new_name = old_inter_file_name.split(f"{str_cluster}_")
                                new_inter_filename = f"{new_name[0]}{new_name[1]}"
                                new_inter_file = Path(
                                    f"{str(new_interface_folder)}/{new_inter_filename}"
                                )
                                shutil.copyfile(inter_file, new_inter_file)
                            elif str(inter_file).endswith(
                                "accessible.txt"
                            ):  # acccessible residues file
                                new_inter_file = Path(
                                    f"{str(new_interface_folder)}/{old_inter_file_name}"
                                )
                                shutil.copyfile(inter_file, new_inter_file)
