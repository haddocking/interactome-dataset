# This is a script to identify possible interfaces in the receptor molecules of a benchmark and cluster
# them according to the relative distances from the observed interfaces with different ligands.

# It gets the information from https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues and can accept
# different filters:
#   - Exclusion of Immunoglobulins: -outIg
#   - Exclusion of Homodimers: -outHD
#   - Only 1:1 PPI: -outHM (HomoMultimers)
#   - Only 2 protein molecules in a PDB: -t2m
#   - Check partners coverage: -cpc. Both protein in a PPI, need to have a coverage higher than the established cutoff
#   - Check biological interaction: -cbi. Performed by ProdigyCrystal
#   - All Out: -AO. Exclude Ig, HD and HM.

# It gets the information of the PDBs in a benchmark to check if the interfaces are present in the strcuture. For that it
# is crucial that those PDBs are correctly renumbered. For doing that, it is strongly recommended to use the Open Software
# PDBrenum.

# For doing the clustering, it first calculates the residue-residue distances in the PDB structure from the euclidean distance
# of their Alpha Carbons. To calculate the D-values for the clustering it does 4 different approaches, of those, the preferred one
# is the sine one.

# The optional files prodigy_results_pdb.txt and interactome_db.txt stores information on biological tests and lists of residues
# for the receptors in a benchmark, respectively.

# Execution: python cluster_interfaces.py location/of/the/benchmark [prodigy_results_pdb] [interactome_db] [-outIg] [-outHD] [-outHM] [-AO] [-cpc] [-cbi]
# e.g. $ python ./scripts/cluster_interfaces.py -AO -cbi ./benchmarks/BM-230222/ ./create_benchmark/prodigy_results_pdb.txt
#
# MSc. Jesús López Rivera, Bonvin Lab :)

import argparse
import ast
import json
import logging
import math
import multiprocessing
import os
import shutil
import subprocess
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import scipy.spatial.distance

FORMAT = " %(asctime)s L%(lineno)d %(levelname)s - %(message)s"
logging.basicConfig(format=FORMAT, level="INFO")

INTERFACE_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues"
ALLPDB_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot"
PDB_MOLECULES_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules"
NEGATOME_URL = (
    "http://mips.helmholtz-muenchen.de/proj/ppi/negatome/combined_stringent.txt"
)
PROTEINS_URL = "https://www.ebi.ac.uk/proteins/api/proteins"
QUERY_DB = {}
PDB_DB = []

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

SIGMA = 1.9
COV_CUTOFF = 0.8

PRODIGY = "/trinity/login/jlopez/software/prodigy-cryst/interface_classifier.py"  # Location of Prodigy Crystal


def make_request(url, data):
    """Helper function to make the requests."""
    # The requests might be done multiple times,
    #  use QUERY_DB dictionary to store the results
    #  this way we do not send the same query multiple times
    if url not in QUERY_DB:
        # First time this query is done
        response = requests.get(url)
        if response.status_code == 404:
            # not found == empty or failure of server, retry 2 times
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


def get_pdb_info(pdb_file):
    """Obtain the residues list of a given PDB"""

    pdb_name = Path(pdb_file).stem.split("_")[0]
    chain = pdb_name[4:]
    pdb = open(pdb_file, "r").read()
    pdb_lines = pdb.splitlines()
    residues = []

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
                ):  # There is a problem with the renumbering, better to discard
                    return None
                else:
                    if in_chain == chain:
                        if res not in residues:
                            residues += [res]

    """
    # Check discontinuities (Not used anymore)
    discontinuities = False
    for atom in range(0, len(pdb_lines)-1):
        atom1_features = pdb_lines[atom].split()
        atom2_features = pdb_lines[atom+1].split()
        if atom1_features[0] == "ATOM" and atom2_features[0] == "ATOM":
            if len(atom1_features[4]) == len(chain):
                res_1 = int(atom1_features[5])
            # When number of residue > 1000
            else:
                res_1 = int(atom1_features[4][len(chain):])
            if len(atom2_features[4]) == len(chain):
                res_2 = int(atom2_features[5])
            # When number of residue > 1000
            else:
                res_2 = int(atom2_features[4][len(chain):])

            if res_1 != res_2:
                if res_1 != res_2-1:
                    discontinuities = True
        elif atom1_features[0] == "TER" and atom2_features[0] == "ATOM":
            discontinuities = True
    """

    if len(residues) == 0:
        residues = None

    return residues


def check_if_Ig(uniprot_ID, Ig_words):
    "Use the Uniprot API to find if the given Uniprot ID is an Immunoglobulin."

    url = f"{PROTEINS_URL}/{uniprot_ID}"
    Ig = False

    # Access the API
    try:
        protein_data = make_request(url, None)
    except Exception as e:
        logging.info(f"Could not make InterfaceResidues request for {uniprot_ID}, {e}")
        return None

    if protein_data:
        if len(protein_data) == 0:  # Check if the entry is empty.
            logging.info(f"Proteins API for {uniprot_ID} was empty")
            return False
        else:
            all_names = []  # Store all possible names for the given Uniprot ID
            if "recommendedName" in protein_data["protein"].keys():  # Main name
                try:
                    main_name = protein_data["protein"]["recommendedName"]["fullName"][
                        "value"
                    ]
                    all_names += [main_name]
                except:
                    pass
            if (
                "alternativeName" in protein_data["protein"].keys()
            ):  # Other posible names
                try:
                    for alt_name in protein_data["protein"]["alternativeName"]:
                        all_names += [alt_name["fullName"]["value"]]
                except:
                    pass

            if (
                len(all_names) == 0
            ):  # Something is wrong with this UniprotID API or it's empty
                logging.info(
                    f"Protein API for {uniprot_ID} did NOT provide a valid name"
                )
                return None

            else:
                # Check if any of the Immunoglobulin words is in any of the names that the Uniprot ID has
                for Ig_word in Ig_words:
                    if Ig_word in all_names:
                        Ig = True
                        break  # No need to check more

    return Ig


def check_multimers_receptor(ligand_id, receptor_id, rec_int_pdbs, outHD):
    """Check if the receptor protein interacts always as a multimer in a given list of PDBs"""

    confirmed_pdbs = []  # PDBs to accept

    # Get info from the API
    url = f"{INTERFACE_URL}/{ligand_id}"
    # The request might fail
    try:
        interface_data = make_request(url, None)
    except Exception as e:
        logging.debug(f"Could not make InterfaceResidues request for {ligand_id}, {e}")
        return False, None

    if interface_data and len(interface_data) != 0:  # Check if the entry is empty
        if len(interface_data[ligand_id]) == 0:  # Check if the entry is empty
            logging.debug(f"InterfaceResidues for {ligand_id} was empty")
        else:
            # Iterate every interacting partner
            for partner in interface_data[ligand_id]["data"]:
                partner_id = partner["accession"]
                if partner_id == ligand_id:  # Ligand is a Homodimer
                    if outHD:  # Exclude homodimers
                        logging.info(f"{ligand_id} discarded due to being an homodimer")
                        return False, None

                elif partner_id == receptor_id:  # Receptor-partner interaction
                    for residue in partner["residues"]:  # Iterate every residue
                        # Check if any of the PDBs for the residue is in the preselected PDBs
                        for lig_int_pdb in residue["interactingPDBEntries"]:
                            lig_int_pdb_id = lig_int_pdb["pdbId"]
                            if (
                                lig_int_pdb_id in rec_int_pdbs
                            ):  # PDB is in the preselected ones
                                # Count the number of molecules for the receptor protein in the interacting pdb
                                if "chainIds" in lig_int_pdb.keys():
                                    chain_Ids = lig_int_pdb["chainIds"].split(",")
                                    if (
                                        len(chain_Ids) == 1
                                    ):  # There is only one partner molecule
                                        if lig_int_pdb_id not in confirmed_pdbs:
                                            confirmed_pdbs += [
                                                lig_int_pdb_id
                                            ]  # Accept PDB_ID
                                    rec_int_pdbs.remove(
                                        lig_int_pdb_id
                                    )  # Avoid reprocessing same data

    if len(confirmed_pdbs) == 0:
        confirmed_pdbs = None

    return confirmed_pdbs


def check_only_two(candidate_pdbs):
    """Check if there are only 2 protein molecules in the PDBs from a given list"""

    confirmed_pdbs = []

    for pdb in candidate_pdbs:
        # Get data from the PDB Molecules API
        pdb_url = f"{PDB_MOLECULES_URL}/{pdb}"
        # The request might fail.
        try:
            pdb_data = make_request(pdb_url, None)
        except Exception as e:
            logging.debug(f"Could not make PDB request for {pdb}, {e}")
            return None

        if pdb_data:
            n_proteins = 0
            # Iterate every molecule in the pdb
            for entity in pdb_data[pdb]:
                # Check if entity is a protein, to deal with cofactors and other atoms
                if "molecule_type" in entity.keys():
                    if "polypeptide" in entity["molecule_type"]:  # It is a protein
                        if "in_struct_asyms" in entity.keys():
                            n_proteins += len(
                                entity["in_struct_asyms"]
                            )  # Every chain is a different protein
                            if (
                                n_proteins > 2
                            ):  # Avoid wasting time in processing everything
                                break

            if n_proteins == 2:  # There are only 2 proteins, accept PDB
                confirmed_pdbs += [pdb]

    if len(confirmed_pdbs) == 0:
        confirmed_pdbs = None

    return confirmed_pdbs


def check_partners_coverage(receptor_id, rec_len, ligand_id, int_pdbs):
    "Check if the coverage of the interacting partners surpass the established threshold"

    # Get the Uniprot length for the ligand
    lig_interface_url = f"{INTERFACE_URL}/{ligand_id}"
    try:
        lig_interface_data = make_request(lig_interface_url, None)
    except Exception as e:
        logging.debug(f"Could not make INTERFACE request for {ligand_id}, {e}")
        return None

    if lig_interface_data:
        if (
            lig_interface_data and len(lig_interface_data) != 0
        ):  # Check if the entry is empty
            if len(lig_interface_data[ligand_id]) == 0:  # Check if the entry is empty.
                logging.debug(f"InterfaceResidues for {ligand_id} was empty")
            else:  # There is information
                lig_len = len(
                    lig_interface_data[ligand_id]["sequence"]
                )  # Lenght of the receptor according to Uniprot

    confirmed_pdbs = []

    # Information from the PDBs
    rec_url = f"{ALLPDB_URL}/{receptor_id}"
    lig_url = f"{ALLPDB_URL}/{ligand_id}"

    # Check receptor coverage
    try:
        rec_allpdb_data = make_request(rec_url, None)
    except Exception as e:
        logging.debug(f"Could not make AllPDB request for {receptor_id}, {e}")
        return None

    if rec_allpdb_data:
        rec_pdbs = []

        # Iterate every pdb for the receptor to check if it is one of the preselected ones
        for rec_pdb_structure in rec_allpdb_data[receptor_id]["mappings"]:
            pdb_id = rec_pdb_structure["entry_id"]
            if pdb_id in int_pdbs:  # It is
                if (
                    len(rec_pdb_structure["segments"]) == 1
                ):  # *Jesus: This can be improved I guess
                    rec_pdb_len = len(rec_pdb_structure["segments"][0]["pdb_sequence"])
                    rec_pdb_coverage = float(rec_pdb_len / rec_len)
                    if rec_pdb_coverage >= COV_CUTOFF:  # Successful receptor coverage
                        rec_pdbs += [pdb_id]
                int_pdbs.remove(pdb_id)
            if len(int_pdbs) == 0:
                break

    if len(rec_pdbs) > 0:
        # Check ligand coverage
        try:
            lig_allpdb_data = make_request(lig_url, None)
        except Exception as e:
            logging.debug(f"Could not make AllPDB request for {ligand_id}, {e}")
            return None

        if lig_allpdb_data:
            # Iterate every pdb for the ligand to check if it is one of the receptor-selected ones
            for lig_pdb_structure in lig_allpdb_data[ligand_id]["mappings"]:
                pdb_id = lig_pdb_structure["entry_id"]
                if pdb_id in rec_pdbs:
                    if len(lig_pdb_structure["segments"]) == 1:
                        lig_pdb_len = len(
                            lig_pdb_structure["segments"][0]["pdb_sequence"]
                        )
                        lig_pdb_coverage = float(lig_pdb_len / lig_len)
                        if lig_pdb_coverage >= COV_CUTOFF:  # Successful ligand coverage
                            confirmed_pdbs += [pdb_id]  # It's a valid PDB

                    rec_pdbs.remove(pdb_id)
                if len(rec_pdbs) == 0:
                    break

    if len(confirmed_pdbs) == 0:
        confirmed_pdbs = None

    return confirmed_pdbs


def get_interface(partner, confirmed_int_pdbs):
    "Get the receptor interface residues for the interaction with a partner"

    residue_list = {}
    for resi in partner["residues"]:
        for int_pdb in resi["interactingPDBEntries"]:
            int_pdb_id = int_pdb["pdbId"]
            if int_pdb_id in confirmed_int_pdbs:
                if int_pdb_id not in residue_list.keys():
                    residue_list[int_pdb_id] = []

                interface_start = int(resi["startIndex"])
                interface_end = int(resi["endIndex"])

                if interface_start != interface_end:
                    # Interface has a few residues.
                    interface = list(range(interface_start - 1, interface_end))
                else:
                    # Interface is only one residue.
                    interface = [interface_start]

                    # Keep track of the residues
                    for residue in interface:
                        if residue not in residue_list[int_pdb_id]:
                            residue_list[int_pdb_id] += [residue]

    if len(residue_list) == 0:
        residue_list = None

    return residue_list


def get_interface_data(
    receptor_id,
    real_receptor_id,
    receptor_pdb,
    rec_residues,
    prodigy_results_pdb,
    data_folder,
    outHD,
    outIg,
    Ig_words,
    outHM,
    t2m,
    cpc,
    cbi,
):
    """Obtain interface data of a given uniprot ID"""

    # Inputs:
    #   receptor id: uniprot_id of the receptor (+lowecase letter in case of more than one structure for that receptor)
    #   real_receptor_id: uniprot_id of the receptor
    #   receptor_pdb: path/to/pdb file of the receptor structure
    #   rec_residues: list of residues of the receptor structure
    #   prodigy_results_pdb: dictionary that stores information of previous ProdigyCrystal runs that tell if the interactions are
    #                        biological or not
    #   data_folder: place where ProdigyCrystal is going to download the pdbs to proceed with the evaluation
    #   outHD: flag to exclude homodimers
    #   outIg: flag to exclude immunoglobulins
    #   Ig_words: list of words used to find out if the protein is an immunoglobulin
    #   outHM: flag to exclude homomultimers (i.e. only 1:1 PPI)
    #   t2m: flag to consider only PDBs in which there are only 2 protein molecules
    #   cpc: flag to "check partners coverage", i.e. check that in the PDB where the interaction has been seen, the coverage
    #        of both receptor and ligand is greater than the established threshold.
    #   cbi: flah to "check biological interaction", i.e. run ProdigyCrystal to know if the interaction is biological or not

    receptor_interface_data = {}  # Dictionary to store information

    # Get info from the API
    url = f"{INTERFACE_URL}/{real_receptor_id}"
    # The request might fail
    try:
        interface_data = make_request(url, None)
    except Exception as e:
        logging.debug(
            f"Could not make InterfaceResidues request for {real_receptor_id}, {e}"
        )
        pass

    if interface_data and len(interface_data) != 0:  # Check if the entry is empty
        if len(interface_data[real_receptor_id]) == 0:  # Check if the entry is empty.
            logging.debug(f"InterfaceResidues for {receptor_id} was empty")
        else:  # There is information

            rec_len = len(
                interface_data[real_receptor_id]["sequence"]
            )  # Lenght of the receptor according to Uniprot
            # Iterate every interacting partner
            for partner in interface_data[real_receptor_id]["data"]:
                partner_id = partner["accession"]  # Uniprot_id of the partner
                logging.info(
                    f"Evaluating interaction between {receptor_id} and {partner_id}"
                )

                # Check if the interaction is in the Negatome (negative interaction)
                logging.info(
                    f"Checking if {receptor_id} or {partner_id} are positive-negative partner duplicates"
                )
                if real_receptor_id in NEGATOME_DIC.keys():
                    if (
                        partner_id in NEGATOME_DIC[real_receptor_id]
                    ):  # Receptor:Negative partner interaction
                        logging.warning(
                            f"{partner_id} discarded due to positive-negative partner duplicate"
                        )
                        continue  # Avoid misleading information
                if partner_id in NEGATOME_DIC.keys():
                    if (
                        real_receptor_id in NEGATOME_DIC[partner_id]
                    ):  # Negative partner: Receptor interaction
                        logging.warning(
                            f"{partner_id} discarded due to positive-negative partner duplicate"
                        )
                        continue  # Avoid misleading information

                if outHD:  # Exclude Homodimers
                    if (
                        partner_id == real_receptor_id
                    ):  # The receptor itself is a homodimer
                        logging.warning(
                            f"{receptor_id} discarded due to being a homodimer"
                        )
                        receptor_interface_data = None
                        break  # Discard the entire receptor. *Jesus: This can be handled in a more thorough manner

                if outIg:  # Exclude Immunoglobulins
                    logging.info(f"Checking if {partner_id} is an Immunoglobulin")
                    isIg = check_if_Ig(partner_id, Ig_words)
                    if isIg:  # Interacting partner is an Immunoglobulin
                        logging.warning(
                            f"{partner_id} is an Immunoglobulin, discarding it"
                        )
                        continue

                if (
                    outHM
                ):  # Check if receptor and partner interact as 1 vs 1 in any of the PDB structures
                    # First step: check if ligand is always a multimer
                    logging.info(f"Checking ligand multimers for {partner_id}")
                    excluded_int_pdbs = []  # Multimer pdbs
                    included_int_pdbs = []  # Not multimer pdbs
                    for residue in partner["residues"]:
                        for int_pdb in residue[
                            "interactingPDBEntries"
                        ]:  # Every PDB in which the residue was seen in the interface
                            int_pdb_id = int_pdb[
                                "pdbId"
                            ]  # PDB ID of the structure where the interaction was seen
                            if (
                                int_pdb_id not in excluded_int_pdbs
                                and int_pdb_id not in included_int_pdbs
                            ):  # Avoid repeated pdbs
                                # Count the number of molecules for the ligand protein in the interacting pdb
                                if "chainIds" in int_pdb.keys():
                                    chain_Ids = int_pdb["chainIds"].split(",")
                                    if (
                                        len(chain_Ids) == 1
                                    ):  # There is only one partner molecule
                                        included_int_pdbs += [
                                            int_pdb_id
                                        ]  # Accept PDB_ID
                                    else:  # There is more than one partner molecule
                                        excluded_int_pdbs += [
                                            int_pdb_id
                                        ]  # Discard PDB_ID

                    if (
                        len(included_int_pdbs) > 0
                    ):  # There are PDBs in which there is only 1 molecule of the partner
                        # Second step: check if receptor is always a multimer in the included interacting PDBs.
                        logging.info(f"Checking receptor multimers for {receptor_id}")
                        confirmed_int_pdbs = check_multimers_receptor(
                            partner_id, real_receptor_id, included_int_pdbs, outHD
                        )
                        if confirmed_int_pdbs != None:  # There are 1:1 PPI interactions

                            # (Optional step): check if only two protein molecules in the pdb structure
                            if t2m:
                                logging.info(
                                    f"Checking only 2 proteins in the pdb files."
                                )
                                confirmed_int_pdbs = check_only_two(confirmed_int_pdbs)
                                if confirmed_int_pdbs == None:
                                    logging.warning(
                                        f"{receptor_id}_{partner_id} discarded due to no pdb with only 2 proteins"
                                    )
                                    continue
                        else:
                            logging.warning(
                                f"{receptor_id}_{partner_id} discarded due to receptor only interacts as multimer"
                            )
                            continue
                    else:
                        logging.warning(
                            f"{receptor_id}_{partner_id} discarded due to ligand only interacts as multimer"
                        )
                        continue
                else:
                    confirmed_int_pdbs = partner["allPDBEntries"]

                if (
                    cpc
                ):  # Check if coverage of the proteins in the interacting PDBs is over the established threshold
                    logging.info(
                        f"Checking partners coverage for {receptor_id}_{partner_id}"
                    )
                    confirmed_int_pdbs = check_partners_coverage(
                        real_receptor_id, rec_len, partner_id, confirmed_int_pdbs
                    )
                    if confirmed_int_pdbs == None:
                        logging.warning(
                            f"{receptor_id}_{partner_id} discarded due to coverage under cutoff"
                        )
                        continue

                if (
                    cbi
                ):  # Check: if the interactions for the confirmed pdbs are biological by Prodigy Crystal
                    logging.info(
                        f"Checking biological interaction for {receptor_id}-{partner_id}"
                    )
                    confirmed_int_pdbs, prodigy_results_pdb = validate_interaction(
                        real_receptor_id,
                        partner_id,
                        confirmed_int_pdbs,
                        prodigy_results_pdb,
                        data_folder,
                    )
                    if confirmed_int_pdbs == None:
                        logging.warning(
                            f"{partner_id} discarded due to no biological interaction"
                        )
                        continue

                # Get interface residues of the single PPI interaction
                logging.info(
                    f"Checking matching interfaces in the structure of {receptor_id}"
                )
                interface = get_interface(partner, confirmed_int_pdbs)
                if interface != None:  # There are valid PDBs and valid residues

                    # Check the interface of each PDB separatedly
                    for int_pdb in interface.keys():
                        inter_check = True
                        # Check if interface residues are present in the receptor structure
                        for residue in interface[int_pdb]:
                            if residue not in rec_residues:
                                inter_check = False
                                break

                        if inter_check:  # The interface is valid
                            if partner_id not in receptor_interface_data:
                                logging.info("It's a match")
                                receptor_interface_data[partner_id] = []

                            # Add residues to partner interface
                            for residue in interface[int_pdb]:
                                if residue not in receptor_interface_data[partner_id]:
                                    receptor_interface_data[partner_id] += [residue]

                    if (
                        partner_id not in receptor_interface_data.keys()
                    ):  # No valid interfaces in receptor structure
                        logging.warning(
                            f"{partner_id} discarded due to interface with {receptor_id}_{receptor_pdb} not in receptor structure"
                        )
                else:
                    logging.warning(
                        f"{partner_id} discarded due to interface with {receptor_id}_{receptor_pdb} not valid"
                    )

    if (
        receptor_interface_data != None
    ):  # In outHD it can be None if receptor is a homodimer
        if len(receptor_interface_data) == 0:
            receptor_interface_data = None

    return receptor_interface_data, prodigy_results_pdb


def validate_interaction(
    prot_a, prot_b, pre_sel_pdbs, prodigy_results_pdb, data_folder
):
    """Check if the interaction between two proteins is biological."""

    biological_pdbs = []

    # Interface data Protein A
    url_a = f"{INTERFACE_URL}/{prot_a}"
    try:
        interface_data_a = make_request(url_a, None)
    except Exception as e:
        logging.debug(f"Could not make InterfaceResidues request for {prot_a}, {e}")
        return None

    # Interface data Protein B
    url_b = f"{INTERFACE_URL}/{prot_b}"
    try:
        interface_data_b = make_request(url_b, None)
    except Exception as e:
        logging.debug(f"Could not make InterfaceResidues request for {prot_b}, {e}")
        return None

    # Check if has been processed totally or partially before
    protein_pair = f"{prot_a}-{prot_b}"
    if (
        protein_pair not in prodigy_results_pdb.keys()
    ):  # The pair has not been processed before
        prodigy_results_pdb[protein_pair] = {}
    else:
        # The pair has been processed before
        not_redo_pdbs = []
        for pre_sel_pdb in pre_sel_pdbs:
            if (
                pre_sel_pdb in prodigy_results_pdb[protein_pair].keys()
            ):  # The pdb for the pair has been processed before
                not_redo_pdbs += [pre_sel_pdb]
        # Avoid redoing already done pdbs.
        for not_redo_pdb in not_redo_pdbs:
            if not_redo_pdb in pre_sel_pdbs:
                pre_sel_pdbs.remove(not_redo_pdb)

    # Obtain a dictionary with pdb_id : chain_id for both proteins. *Jesus: I used a previous script for this
    pdb_chain_a = dict(
        [
            (k["pdbId"], k["chainIds"])
            for e in interface_data_a[prot_a]["data"]
            if e["accession"] == prot_b
            for j in e["residues"]
            for k in j["interactingPDBEntries"]
            if k["pdbId"] in pre_sel_pdbs
        ]
    )

    pdb_chain_b = dict(
        [
            (k["pdbId"], k["chainIds"])
            for e in interface_data_b[prot_b]["data"]
            if e["accession"] == prot_a
            for j in e["residues"]
            for k in j["interactingPDBEntries"]
            if k["pdbId"] in pre_sel_pdbs
        ]
    )

    # Check that the pdbs are the same for both proteins
    checked_pdb_chain_a = {}
    checked_pdb_chain_b = {}

    for pdb in pdb_chain_a.keys():
        if pdb in pdb_chain_b.keys():
            checked_pdb_chain_a[pdb] = pdb_chain_a[pdb]
            checked_pdb_chain_b[pdb] = pdb_chain_b[pdb]
    assert checked_pdb_chain_a.keys() == checked_pdb_chain_b.keys()

    # Run Prodigy
    for pdb in checked_pdb_chain_a:
        renum_pdb_f = Path(f"{data_folder}/{pdb}_renum.pdb")
        if (
            renum_pdb_f.exists()
        ):  # Check if the pdb has already been downloaded and renumbered
            pdb_f = renum_pdb_f
        else:  # Check if the pdb has already been downloaded
            pdb_f = Path(f"{data_folder}/{pdb}.pdb")

        if (
            not pdb_f.exists()
        ):  # Download the pdb. *Jesus: this can be done by PDBrenum too
            cmd = f'pdb_fetch {pdb} | grep "ATOM" | pdb_tidy'
            p = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )

            out, err = p.communicate()

            with open(pdb_f, "w") as fh:
                fh.write(out.decode("utf-8"))

        chain_a_list = checked_pdb_chain_a[pdb].split(",")
        chain_b_list = checked_pdb_chain_b[pdb].split(",")
        job_list = []
        pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
        for chain_a in chain_a_list:
            for chain_b in chain_b_list:
                if chain_a != chain_b:
                    # Prodigy wont work if chains are the same..!
                    job_list.append(
                        (
                            pdb_f,
                            chain_a,
                            chain_b,
                            pool.apply_async(
                                run_prodigy, args=(pdb_f, chain_a, chain_b)
                            ),
                        )
                    )

        for element in job_list:
            pdb_f, chain_a, chain_b, proc = element
            interface_class = proc.get()
            PDB_DB.append((prot_a, prot_b, pdb_f, chain_a, chain_b, interface_class))
            logging.info(
                f"===> {prot_a}_{prot_b} = {pdb} {chain_a}-{chain_b} "
                f"= {interface_class}"
            )
            if interface_class == "BIO":
                # If any chain-chain interaction for this pdb is biological, the evaluation is biological
                if pdb not in prodigy_results_pdb[protein_pair].keys():
                    prodigy_results_pdb[protein_pair][pdb] = "biological"

        pool.close()

    for pdb in pre_sel_pdbs:
        if (
            pdb not in prodigy_results_pdb[protein_pair].keys()
        ):  # Those pdbs haven't been evaluated as biological before, so they are not
            prodigy_results_pdb[protein_pair][pdb] = "NOT biological"

    # Confirm biological pdbs
    for processsed_pdb in prodigy_results_pdb[protein_pair].keys():
        if (
            prodigy_results_pdb[protein_pair][processsed_pdb] == "biological"
            and processsed_pdb in pre_sel_pdbs
        ):
            biological_pdbs += [processsed_pdb]

    if len(biological_pdbs) == 0:
        biological_pdbs = None

    return biological_pdbs, prodigy_results_pdb


def run_prodigy(pdb_f, chain_i, chain_j):
    """Wrapper for running PRODIGY-CRYSTAL."""

    cmd = f"python {PRODIGY} {pdb_f} --selection {chain_i} {chain_j}"

    p = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )

    out, _ = p.communicate()
    try:
        interface_class = out.decode("utf-8").split(os.linesep)[-2].split()[2]
    except Exception as e:
        logging.debug(f"{pdb_f}, {e}")
        interface_class = "NONE"

    return interface_class


def calc_distances(pdb_f, uniprot_id, folder):
    """Calculate the distances between residues in a molecule, and generate a residue_distance matrix"""

    coord_dic = {}
    distances = {}

    # Obtain the coordinates
    with open(pdb_f) as fh:
        for line in fh:
            if line.startswith("ATOM"):  # Avoid REMARKS, END and other lines
                atom_name = line[12:16].strip()
                # Only uses the CA, other methods can be used here
                if atom_name == "CA":
                    resnum = int(line[22:26])  # residue number
                    # Coordinates:
                    if line[30] == "-":  # Negative coordinate
                        x = float(line[30:38])
                    else:  # Positive coordinate
                        x = float(line[31:38])
                    if line[38] == "-":
                        y = float(line[38:46])
                    else:
                        y = float(line[39:46])
                    if line[46] == "-":
                        z = float(line[46:54])
                    else:
                        z = float(line[47:54])
                    coord_dic[resnum] = (x, y, z)

    # Calculate the distances
    output_matrix = []
    residue_names = []  # Numbers of residues for the later plot of the matrix

    for (
        residue_a
    ) in coord_dic.keys():  # Calculate the distances of residue a with every residue b
        residue_a_matrix = []
        residue_names += [str(residue_a)]
        c_a = coord_dic[residue_a]  # Coordinates of residue a
        for residue_b in coord_dic.keys():
            pair_ab = (
                f"{str(residue_a)}_{str(residue_b)}"  # ResidueA_ResidueB -> e.g. 1_2
            )
            c_b = coord_dic[residue_b]
            distance_a_b = scipy.spatial.distance.euclidean(
                c_a, c_b
            )  # Distance from coordinates a to coordinates b
            distances[pair_ab] = distance_a_b  # Dictionary to be returned
            residue_a_matrix += [distance_a_b]  # To plot the matrix
        output_matrix += [residue_a_matrix]  # To plot the matrix

    # Plot the residue-distance matrix
    mat = pd.DataFrame(output_matrix, columns=residue_names, index=residue_names)
    fig = plt.figure()
    ax = plt.gca()
    im = ax.matshow(mat, interpolation="none", cmap="viridis_r")
    fig.colorbar(im)
    ax.set_title(f"{uniprot_id} residue-residue distances matrix", pad=50)
    ax.set_axis_on()
    fig.tight_layout()
    plt.savefig(f"{folder}/{uniprot_id}_distances_matrix.png")
    plt.close()

    return distances


def partial_d_value(residue_a, residue_b, receptor_distances, sigma):
    residue_pair = f"{str(residue_a)}_{str(residue_b)}"
    a_b = receptor_distances[residue_pair]
    partial_value = math.exp((-a_b) / (4.0 * (sigma**2)))
    return partial_value


def get_d_value_components(
    receptor_distances, ligand1_contacts, ligand2_contacts, sigma
):
    """Calculate the d-value for a pair of ligands for the same receptor"""
    I = 0.0
    J = 0.0
    I_J = 0.0
    for residue_Ii in ligand1_contacts:
        for residue_Ij in ligand1_contacts:
            I += partial_d_value(residue_Ii, residue_Ij, receptor_distances, sigma)
        for residue_Jj in ligand2_contacts:
            I_J += partial_d_value(residue_Ii, residue_Jj, receptor_distances, sigma)
    for residue_Ji in ligand2_contacts:
        for residue_Jj in ligand2_contacts:
            J += partial_d_value(residue_Ji, residue_Jj, receptor_distances, sigma)

    return I, J, I_J


def plot_matrix(
    receptor_d_matrix, ligand_names, receptor_id, matrix_type, matrix_folder
):
    """Plot a D-values matrix"""

    df = pd.DataFrame(receptor_d_matrix, columns=ligand_names, index=ligand_names)
    # print(df) # Print the matrix
    fig = plt.figure()
    ax = plt.gca()
    im = ax.matshow(df, interpolation="none", cmap="viridis_r")
    fig.colorbar(im)
    ax.set_xticks(np.arange(len(ligand_names)))
    ax.set_xticklabels(ligand_names)
    ax.set_yticks(np.arange(len(ligand_names)))
    ax.set_yticklabels(ligand_names)
    # Set ticks on both sides of axes on
    ax.tick_params(axis="x", bottom=True, labelbottom=True)
    # Rotate and align bottom ticklabels
    plt.setp(
        [tick.label1 for tick in ax.xaxis.get_major_ticks()],
        rotation=60,
        ha="right",
        va="center",
        rotation_mode="anchor",
    )
    # Rotate and align top ticklabels
    plt.setp(
        [tick.label2 for tick in ax.xaxis.get_major_ticks()],
        rotation=45,
        ha="left",
        va="center",
        rotation_mode="anchor",
    )
    ax.set_title(f"{receptor_id} {matrix_type} matrix", pad=55)
    fig.tight_layout()
    plt.savefig(f"{str(matrix_folder)}/{receptor_id}_interface_matrix.png")
    plt.close()


def load_negatome_data(negatome_url):
    """Load the negatome data."""
    negatome_f = negatome_url.split("/")[-1]
    if not Path(negatome_f).exists():
        subprocess.run(["wget", negatome_url], check=True)

    negatome_dic = {}
    with open(negatome_f, "r") as fh:
        for line in fh.readlines():
            partner_a, partner_b = line.split()
            if partner_a not in negatome_dic.keys():
                negatome_dic[partner_a] = []
            if partner_b not in negatome_dic[partner_a]:
                negatome_dic[partner_a] += [partner_b]

            if partner_b not in negatome_dic.keys():
                negatome_dic[partner_b] = []
            if partner_a not in negatome_dic[partner_b]:
                negatome_dic[partner_b] += [partner_a]

    return negatome_dic


if __name__ == "__main__":

    # Process the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", help="Introduce benchmark location")
    parser.add_argument(
        "data_folder", help="Introduce folder location to store the downloaded PDBs"
    )
    parser.add_argument(
        "prodigy_results_pdb",
        nargs="?",
        help="Introduce prodigy_results_pdb.txt file if possible",
    )
    parser.add_argument(
        "interactome_db",
        nargs="?",
        help="Introduce interactome_db.txt file if possible",
    )
    parser.add_argument(
        "-outHD",
        action="store_true",
        help="Exclude Homodimers interactions from the dataset",
    )
    parser.add_argument(
        "-outIg",
        action="store_true",
        help="Exclude Immunoglobulin interactions from the dataset",
    )
    parser.add_argument(
        "-outHM",
        action="store_true",
        help="Exclude homomultimers interactions from the dataset",
    )
    parser.add_argument(
        "-AO",
        action="store_true",
        help="All Out, exclude HD, Ig and HM from the dataset",
    )
    parser.add_argument(
        "-t2m",
        action="store_true",
        help="Take only interactions with 2 proteins in a pdb (that's 2 match, that's too much)",
    )
    parser.add_argument(
        "-cpc",
        action="store_true",
        help="Check partners coverage, i.e. both interacting molecules surpass the coverage cutoff",
    )
    parser.add_argument(
        "-cbi",
        action="store_true",
        help="Check biological interaction, i.e. run ProdigyCrystal to check if the interaction is biological",
    )
    args = parser.parse_args()
    bm_directory = Path(args.benchmark)
    data_folder = args.data_folder
    prodigy_results_pdb_file = args.prodigy_results_pdb
    prodigy_results_pdb = {}
    interactome_db_file = args.interactome_db
    interactome_db = {}
    outHD = args.outHD
    outIg = args.outIg
    outHM = args.outHM
    AO = args.AO
    t2m = args.t2m
    cpc = args.cpc
    cbi = args.cbi

    if AO:  # Exclude all possible things
        outHD = outIg = outHM = True
        logging.info(f"Modes outHD, outIg and outHM selected")

    if (
        t2m and not outHM
    ):  # Mode only 2 molecules can only go with the mode out homomultimers
        outHD = True
        logging.info(f"Mode outHM also activated")

    if (
        outHM and not cpc
    ):  # It can be not activated but it is advisable until new pipeline
        logging.warning("<<<<< CPC mode not activated, please consider using it >>>>>")

    # Running Prodigy Crystal to check if the interactions are biological is a slow process, to make it faster,
    # the script can re-use the results of running it before, stored in a file called prodigy_results_pdb.txt.
    # Check if prodigy_results_pdb.txt file provided
    if prodigy_results_pdb_file is None:
        logging.info(f"No prodigy_results_pdb.txt file introduced")
    else:
        if Path(prodigy_results_pdb_file).exists():
            logging.info(
                "<<<<<< prodigy_results_pdb.txt file found, reusing it. >>>>>>"
            )
            read_prodigy_results_pdb = open(prodigy_results_pdb_file, "r").read()
            prodigy_results_pdb = ast.literal_eval(
                read_prodigy_results_pdb
            )  # Make a python dictionary out of the file
        else:  # File not found
            logging.info("<<<<<< prodigy_results_pdb.txt file not found. >>>>>>")

    # The script can run faster if we get the information of the residues inside each PDB structure done in the dataset
    # generation. It not provided the script will get them.
    if interactome_db_file is None:
        logging.info(f"No interactome_db.txt file introduced")
        abini = True  # Ab initio, redo the obtaining of receptor residues.
    else:
        if Path(interactome_db_file).exists():
            logging.info("<<<<<< interactome_db.txt file found, reusing it. >>>>>>")
            read_interactome_db = open(interactome_db_file, "r").read()
            interactome_db = ast.literal_eval(
                read_interactome_db
            )  # Make a python dictionary out of the file
            abini = False
        else:  # File not found
            logging.info("<<<<<< interactome.txt file not found. >>>>>>")
            abini = True

    # Folder to store the generated interface_data
    interface_ana_folder = f"{str(bm_directory)}/interface_analysis"
    if Path(interface_ana_folder).exists():
        shutil.rmtree(interface_ana_folder)
    os.mkdir(interface_ana_folder)

    # Folder to store the downloaded pdbs by ProdigyCrystal
    if not Path(data_folder).exists():
        os.mkdir(data_folder)

    # Dictionaries to store the data
    distance_matrix = {}
    contacts_dictionary = {}
    d_values_matrix = {"standard": {}, "1-cos": {}, "sin2": {}, "sin": {}}
    # Make a folder for every type of matrix. *Jesus: this can be reduced to only sin now, indeed.
    for matrix in d_values_matrix.keys():
        matrix_folder = f"{interface_ana_folder}/{matrix}"
        os.mkdir(matrix_folder)

    """ Take the data from the NEGATOME to discard the found interactions that appear as negative there """
    NEGATOME_DIC = load_negatome_data(NEGATOME_URL)

    """ Iterate every folder in the benchmark """
    for folder in bm_directory.iterdir():
        folder_name = folder.stem
        molecules = folder_name.split("_")
        receptor_id = molecules[0]

        """ Some molecules in the dataset have more than one pdb file that cover different
        regions of the protein. They are written as 'b' or 'c', so to get the real ID
        this lowercase letter has to be removed in case there is. (*) """
        if str(receptor_id[-1]).islower():
            real_receptor_id = receptor_id[:-1]
        else:
            real_receptor_id = receptor_id

        # Avoid other folders and repeated recpetors.
        if (
            folder.is_dir()
            and len(molecules) == 3
            and receptor_id not in distance_matrix.keys()
        ):

            """Find the receptor pdb file inside the folder and calculate the residue distances"""
            for file in folder.iterdir():
                if str(file).endswith("r_b.pdb"):
                    receptor_file = file
                    receptor_pdb = (
                        f"{file.stem[:4]}_{file.stem[4]}"  # pdb_id + chain_id
                    )
            if abini:  # interactome_db.txt not provided or not found
                logging.info(
                    f"Obtaining the structure residues for receptor {receptor_id}_{receptor_pdb}"
                )
                rec_residues = get_pdb_info(
                    receptor_file
                )  # Get a list of residues of the PDB structure
            else:  # already processed
                rec_residues = interactome_db[real_receptor_id][receptor_pdb][
                    "residues"
                ]

            if rec_residues != None:  # If None, there is some error in the PDB file
                """Calculate the inter-residue distances of the receptor structure"""
                logging.info(
                    f"Calculating residue-residue distances for receptor {receptor_id}_{receptor_pdb}"
                )
                distance_matrix[receptor_id] = calc_distances(
                    receptor_file, receptor_id, interface_ana_folder
                )  # Get residue-residue distances

                """ Obtain the interface residues for each ligand for the receptor id """
                logging.info(
                    f"Obtaining interface data for receptor {receptor_id}_{receptor_pdb}"
                )
                (
                    contacts_dictionary[receptor_id],
                    prodigy_results_pdb,
                ) = get_interface_data(
                    receptor_id,
                    real_receptor_id,
                    receptor_pdb,
                    rec_residues,
                    prodigy_results_pdb,
                    data_folder,
                    outHD,
                    outIg,
                    Ig_words,
                    outHM,
                    t2m,
                    cpc,
                    cbi,
                )

                # Since it is a slow process, it is useful to be writing the prodigy_results file at the same time
                # in case of error of abortion and faster usage later.
                with open(prodigy_results_pdb_file, "w") as fh:
                    fh.write(str(prodigy_results_pdb))
                    fh.close()

                if (
                    contacts_dictionary[receptor_id] == None
                ):  # No valid partners for receptor
                    logging.warning(
                        f"No possible ligands found for receptor {receptor_id}, continue"
                    )
                    continue

                """ Obtain the d_vales for the receptor, both as dictionary and matrix """
                for matrix in d_values_matrix.keys():
                    d_values_matrix[matrix][receptor_id] = {}
                receptor_d_matrix = {"standard": [], "1-cos": [], "sin2": [], "sin": []}

                logging.info(f"Processing receptor {receptor_id} interfaces")
                ligand_names = contacts_dictionary[receptor_id].keys()
                logging.info(f"Detected ligands are {ligand_names}")

                if len(ligand_names) > 0:  # There is a similar check in line 941
                    # Calculate D for every pair of ligands (ligand1_ligand2)
                    for ligand_1 in ligand_names:
                        ligand1_d_matrix = {}
                        for (
                            matrix
                        ) in d_values_matrix.keys():  # standard, cos, sin2 and sin
                            ligand1_d_matrix[matrix] = []
                        for ligand_2 in ligand_names:
                            ligand_pair = f"{ligand_1}_{ligand_2}"
                            # Get components of equation ||I||, ||J|| and <I,J>
                            I, J, I_J = get_d_value_components(
                                distance_matrix[receptor_id],
                                contacts_dictionary[receptor_id][ligand_1],
                                contacts_dictionary[receptor_id][ligand_2],
                                sigma=SIGMA,
                            )
                            # Get the different D values for each approach
                            st_d_value = math.sqrt(I + J - 2.0 * I_J)
                            cos_IJ = float(I_J / math.sqrt(I * J))
                            cos_d_value = 1 - cos_IJ
                            sin2_d_value = 1 - cos_IJ**2
                            sin_d_value = math.sqrt(sin2_d_value)

                            # Store values for each approach
                            d_values_matrix["standard"][receptor_id][
                                ligand_pair
                            ] = st_d_value
                            d_values_matrix["1-cos"][receptor_id][
                                ligand_pair
                            ] = cos_d_value
                            d_values_matrix["sin2"][receptor_id][
                                ligand_pair
                            ] = sin2_d_value
                            d_values_matrix["sin"][receptor_id][
                                ligand_pair
                            ] = sin_d_value

                            # Store values for each ligand
                            ligand1_d_matrix["standard"] += [st_d_value]
                            ligand1_d_matrix["1-cos"] += [cos_d_value]
                            ligand1_d_matrix["sin2"] += [sin2_d_value]
                            ligand1_d_matrix["sin"] += [sin_d_value]

                        # Store the value of each D-value approach for each ligand in the dictionary of the receptor for the same approach
                        for d_matrix in ligand1_d_matrix.keys():
                            receptor_d_matrix[d_matrix] += [ligand1_d_matrix[d_matrix]]

                    # Make a matrix plot with the D-values data for each approach and store it
                    for d_matrix in receptor_d_matrix.keys():
                        d_matrix_folder = f"{interface_ana_folder}/{d_matrix}"
                        plot_matrix(
                            receptor_d_matrix[d_matrix],
                            ligand_names,
                            receptor_id,
                            d_matrix,
                            d_matrix_folder,
                        )

                else:
                    logging.warning(f"{receptor_id} has no valid ligands")
            else:
                logging.warning(f"{receptor_id} was discarded due to wrong pdb file")

    # Write the distance_matrix dictionary for reusing it after finishing.
    with open(f"{str(interface_ana_folder)}/distance_matrix.txt", "w") as fh:
        for rec in distance_matrix.keys():
            fh.write(rec + " : ")
            fh.write(str(distance_matrix[rec]) + os.linesep)

    # Write the contact dictionary
    with open(f"{str(interface_ana_folder)}/interface_dictionary.txt", "w") as fh:
        fh.write(
            str(contacts_dictionary)
        )  # *Jesus: probably not the best visual output, but the best way for reusing it later

    # Write the interface matrix
    for int_matrix in d_values_matrix.keys():
        with open(
            f"{str(interface_ana_folder)}/interface_{int_matrix}_matrix.txt", "w"
        ) as fh:
            for rec in d_values_matrix[int_matrix].keys():
                fh.write(rec + " : ")
                fh.write(str(d_values_matrix[int_matrix][rec]) + os.linesep)
