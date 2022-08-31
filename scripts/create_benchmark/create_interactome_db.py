# This script is the first step to create a so-called interactome dataset. This is a dataset which contains positive and negative
# interactions organized in a way that the same receptor protein has a positive partner (interaction found in nature)
# and a negative one (interaction proven that does not happen in nature). In that way, the group receptor-positive-negative can
# be used to evaluate the differences between a real and a non-real interaction after any PPI prediction method. It can also be used
# to train AI models.
#
# Since negative data is the limiting factor and the most difficult to found in literature, it is the starating point of the script.
# It will process the Negatome2.0 database to extract pairs of receptor-negative partners. Then it will try to find good pdb
# structures under the established thresholds:
#   - COV_CUTOFF: Coverage of the structure over which the structure is accepted. E.g: 0.5 = 50% of the structure
#   - RES_CUTOFF: Resolution of the structure under which the structure is accepted. E.g: 3.0 = 3.0 Aº
#   - AA_CUTOFF: Number of aminoacids over which the structure is accepted. E.g: 50 (avoid small peptides)
#
# If the receptor has good pdb structures, then it will be processed in search of positive interacting partners. If there are,
# negative partners and positive partners will undergo the same finding good pdb structure process. Finally, positive partners
# will have a reciprocity test to check that their receptor is also a positive partner for them:
# Receptor -> Positive & Positive -> Receptor
#
# If a positive partner for a certain receptor is also a negative partner for the same receptor, it will be removed from both
# sets to avoid conflicting data.
#
# It can be run in different modes:
#   -nrmd = No Remove Duplicates: only do if really trust the process of the Negatome2.0 database
#   -onia = Off Negative Interface API check: Don't need to check if the negative partners have interactions with other molecules
#   -outHM = Accept only 1:1 PPI, i.e. exclude interactions in which there are multiple molecules of the same protein interacting
#   -t2m = Increased stringency of outHM to only 1:1 interactions but just 2 molecules in a PDB file (very restrictive)
#   -cpc = Check partners coverage: check that both molecules in an interactions have a coverage bigger than COV_THRESHOLD
#   -cbi = Check biological intereaction: check that the interaction between two molecules is biological by the sofware
#          ProdigyCrystal
#
# The file prodigy_results_pdb.txt allows to skip already processed interactions by ProdigyCrystal. While create_benchmark
# is the folder where all the data is going to be stored:
#
# Execution: python create_interactome_db.py [-nrmd] [-onia] [-outHM] [-t2m] [-cpc] [-cbi] location/of/create_benchmark [location/of/prodigy_results_pdb.txt]
# e.g. $ python ./scripts/create_interactome_db.py -outHM -cpc -cbi ./create_benchmark ./create_benchmark/prodigy_results_pdb.txt
#
# MSc. Jesús López Rivera, Bonvin Lab :)

import argparse
import ast
import csv
import json
import logging
import multiprocessing
import os
import subprocess
import time
from pathlib import Path

import matplotlib.pyplot as plt
import requests

FORMAT = " %(asctime)s L%(lineno)d %(levelname)s - %(message)s"
logging.basicConfig(format=FORMAT, level="INFO")

# =====================================================================================
# URLs
PDB_MOLECULES_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules"
INTERFACE_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues"
# Get the list of PDB structures mapping to a UniProt accession sorted by coverage of
#  the protein and, if the same, resolution.
PROTEINS_URL = "https://www.ebi.ac.uk/proteins/api/proteins"
BESTPDB_URL = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures"
ALLPDB_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot"
NEGATOME_URL = (
    "http://mips.helmholtz-muenchen.de/proj/ppi/negatome/combined_stringent.txt"
)

# ======================================================================================
# General databases
QUERY_DB = {}
PDB_DB = []
WRONG_PDBS_DB = []

# ======================================================================================
# Cutoffs for including / descarding pdb structures.
COV_CUTOFF = 0.4  # Coverage >= 40% of the whole strcutre
RES_CUTOFF = 3.0  # Resolution <= 3.0 Aº
AA_CUTOFF = 50  # Number of aminoacids >= 50

# ======================================================================================
# Helper software location

PDBrenum = "/trinity/login/jlopez/software/PDBrenum/PDBrenum.py"
PRODIGY = "/trinity/login/jlopez/software/prodigy-cryst/interface_classifier.py"

# =====================================================================================
# Request functions


def make_request(url, data):
    """Helper function to make the requests."""
    # The requests might be done multiple times, use QUERY_DB dictionary to store the results
    # this way it does not send the same query multiple times
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

    # Wait to avoid get banned out of the API
    time.sleep(0.1)
    return data


def run_PDBrenum(pdb_id, WRONG_PDBS_DB):
    """Wrapper for running PDBrenum."""

    renumbered = False

    if pdb_id not in WRONG_PDBS_DB:  # Avoid redoing a problematic pdb

        # Sometimes there is no need of renumber because the pdb files are already correct, so the name can be different.
        no_renum_destination = f"output_PDB/{pdb_id}.pdb"
        renum_destination = f"output_PDB/{pdb_id}_renum.pdb"

        if not Path(renum_destination).exists():  # Check if already renumbered
            logging.info(f"Renumbering PDB_ID {pdb_id}")

            try:
                subprocess.run(
                    ["python", PDBrenum, "-rfla", pdb_id, "-PDB", "-offz"], check=True
                )
                if not (
                    Path(no_renum_destination).exists()
                    or Path(renum_destination).exists()
                ):  # Sometimes the pdb files are not available
                    logging.warning(
                        f"PDBrenum couldn't download or renum the pdb file, treating as problematic."
                    )
                    WRONG_PDBS_DB += [pdb_id]

                else:
                    logging.info(f"Successful renumbering of {pdb_id}")
                    renumbered = True
                    # Rename in case of no_renum file to unify the output
                    if Path(no_renum_destination).exists():
                        os.rename(no_renum_destination, renum_destination)

            except:
                logging.warning(f"{pdb_id} has a problematic pdb file")
                WRONG_PDBS_DB += [pdb_id]
                return renumbered, WRONG_PDBS_DB
        else:
            logging.info(f"{pdb_id} already renumbered, continue")
            renumbered = True
    else:
        logging.warning(f"{pdb_id} is a problematic pdb, skiping")

    return renumbered, WRONG_PDBS_DB


def check_real_chainID(pdb_id, chain_ID, uniprot_ID):
    "Use the PDB API to find the real chain ID."

    url = f"{PDB_MOLECULES_URL}/{pdb_id}"  # url of the pdb
    url_name = f"{PROTEINS_URL}/{uniprot_ID}"  # url of uniprot where to obtain the names of the molecule
    real_ID = False

    # Get pdb data
    try:
        pdb_data = make_request(url, None)
    except Exception as e:
        logging.info(f"Could not make Molecules request for {pdb_id}, {e}")
        return False

    if pdb_data:
        # Get uniprot data
        try:
            protein_data = make_request(url_name, None)
        except Exception as e:
            logging.info(
                f"Could not make InterfaceResidues request for {uniprot_ID}, {e}"
            )
            return False

        if protein_data:
            # Get different names of the protein
            if len(protein_data) == 0:  # Check if the entry is empty.
                logging.info(f"Proteins API for {uniprot_ID} was empty")
                return False
            else:
                # Add all different names to the list
                all_names = []
                if "recommendedName" in protein_data["protein"].keys():
                    try:  # May fail
                        main_name = protein_data["protein"]["recommendedName"][
                            "fullName"
                        ]["value"]
                        all_names += [main_name]
                    except:
                        pass
                if "alternativeName" in protein_data["protein"].keys():
                    try:  # May fail
                        for alt_name in protein_data["protein"]["alternativeName"]:
                            all_names += [alt_name["fullName"]["value"]]
                    except:
                        pass

                if len(all_names) == 0:  # No valid names
                    logging.info(
                        f"Protein API for {uniprot_ID} did NOT provide a valid name"
                    )
                    return False

                # Check pdb data
                if len(pdb_data[pdb_id]) == 0:  # Check if the entry is empty.
                    logging.info(f"PDB_Moleucles API for {pdb_id} was empty")
                    return False
                else:
                    for entity in pdb_data[pdb_id]:
                        # Check in the chain_ID was correct in the previous API
                        if chain_ID in entity["in_chains"]:
                            if (
                                "polypeptide" in entity["molecule_type"]
                            ):  # Avoid other molecules
                                for molecule_name in entity["molecule_name"]:
                                    if molecule_name in all_names:
                                        logging.info(
                                            f"{pdb_id}, ChainID was correct in the previous API"
                                        )
                                        real_ID = chain_ID
                                        # The chain ID was correct, no need to change it
                                        break
                                else:
                                    continue
                                break

                    if (
                        not real_ID
                    ):  # Couldn't find the protein by the previous procedure
                        for entity in pdb_data[pdb_id]:
                            # Check if the chain_ID was wrong in the previous API
                            if chain_ID in entity["in_struct_asyms"]:
                                if "polypetide" in entity["molecule_type"]:
                                    for molecule_name in entity["molecule_name"]:
                                        if (
                                            molecule_name in all_names
                                        ):  # The chain ID was NOT correct.
                                            logging.info(
                                                f"{pdb_id}, ChainID was NOT correct in the previous API"
                                            )
                                            real_ID = entity["in_chains"][
                                                0
                                            ]  # To avoid repeated structures, just take the first one.
                                            break
                                    else:
                                        continue
                                    break

                if (
                    real_ID == False
                ):  # Couldn't find a matching chain ID in "in_struct_asyms"
                    logging.warning(f"{pdb_id}-{chain_ID} is not a valid chain.")

    return real_ID


def get_pdb_residues(pdb_file, chain):
    """Obtain a list of the resiude numbers of a pdb structure"""
    pdb = open(pdb_file, "r").read()
    pdb_lines = pdb.splitlines()
    residues = []

    # Get residue number
    for atom in pdb_lines:
        atom_features = atom.split()
        if atom_features[0] == "ATOM" and "CA" in atom_features:
            if len(atom_features[4]) == len(chain):
                if atom_features[4] == chain:
                    if atom_features[5][
                        -1
                    ].isalpha():  # E.g: 147 147 147 147A 147A 147B...
                        res = int(atom_features[5][:-1])
                    else:
                        res = int(atom_features[5])
                    if res not in residues:
                        residues += [res]
            # When number of residue > 1000
            else:
                in_chain = atom_features[4][: len(chain)]
                res = atom_features[4][len(chain) :]
                if res[-1].isalpha():  # E.g: 1472 1472 1472 1472A 1472A 1472B...
                    res = int(res[:-1])
                else:
                    res = int(res)
                if (
                    res > 5000
                ):  # There is a problem with the renumbering, better to discard
                    return None

                else:
                    if in_chain == chain:
                        if res not in residues:
                            residues += [res]

    if len(residues) == 0:
        residues = None

    return residues


def get_bestpdb(uniprot_id, WRONG_PDBS_DB):
    "Get the best pdb structures from a given UniprotID."

    logging.info(f"Inspecting pdb structures for {uniprot_id}")
    pdb_f = None
    url = f"{BESTPDB_URL}/{uniprot_id}?pretty=false"

    # The request might fail.
    try:
        bestpdb_data = make_request(url, None)
    except Exception as e:
        logging.debug(f"Could not make BestPDB request for {uniprot_id}, {e}")
        return None, WRONG_PDBS_DB

    if bestpdb_data:

        # Each protein might have different pdb structures with different coverages, so it captures the
        # information of all of those above the cutoffs only if they cover different parts of the protein
        # for later use while getting the interacting partners.
        n_pdb = 0
        max_n_pdb = len(bestpdb_data[uniprot_id])
        conditions = True  # Allows to shorten the search, since structures are ordered by quality.

        pdb_f = {}
        pdb_ids = []

        while (
            conditions and n_pdb < max_n_pdb
        ):  # Switch off if conditions doesn't apply to save some time because structures are ordered by quality:
            # first coverage, and then resolution. If a structure has lower coverage than the cutoff, the rest can
            # be skipped, but if the coverage is greater than the cutoff and so is the resolution, it can exist another
            # structure with smaller coverage than the previous but still greater than the cutoff with lower resolution,
            # and that one can be included.

            # Some structures are not numbered accordingly to their protein residue order, but according to the authors
            # criteria. For later use, it is crucial that proteins are correctly renumbered, for which it is used the
            # software PDBrenum. However, it fails sometimes, so in that cases, the following in quality structure (if it
            # surpasses the cutoffs) will be tried and taken it is succeeds.

            pdb_data = bestpdb_data[uniprot_id][n_pdb]
            coverage = pdb_data["coverage"]
            resolution = pdb_data["resolution"]
            pdb_id = pdb_data["pdb_id"]
            chain_id = pdb_data["chain_id"]

            if pdb_id not in pdb_ids:  # Avoid repeated entries
                pdb_ids += [pdb_id]

                if resolution == None:
                    resolution = 0.0
                if coverage == None:
                    coverage = 0.0

                if coverage >= COV_CUTOFF and resolution > RES_CUTOFF:
                    logging.warning(
                        f">> Found BestPDB: {pdb_id}_{chain_id} cov={coverage:.2f},"
                        f"res={resolution:.1f}A but it is below the cutoffs {COV_CUTOFF:.2f} / {RES_CUTOFF:.2f}A, "
                        "Coverage condition satisfied but not resolution condition, continue"
                    )

                elif coverage < COV_CUTOFF:
                    logging.warning(
                        f">> Found BestPDB: {pdb_id}_{chain_id} cov={coverage:.2f},"
                        f"res={resolution:.1f}A but it is below the cutoffs {COV_CUTOFF:.2f} / {RES_CUTOFF:.2f}A, "
                        "Coverage condition unsatisfied, end"
                    )
                    conditions = False  # No interest in the rest of structures

                elif (
                    coverage >= COV_CUTOFF and resolution <= RES_CUTOFF
                ):  # Valid structure
                    logging.info(
                        f"> Found BestPDB: {pdb_id}_{chain_id}, cov={coverage:.2f}, res={resolution:.1f}A"
                    )

                    if len(pdb_f) == 0:  # First valid structure

                        # Try downloading and renumbering pdb_id
                        renumbered, WRONG_PDBS_DB = run_PDBrenum(pdb_id, WRONG_PDBS_DB)
                        if renumbered:

                            # The chain ID that appears at the best-structures API sometimes is not the real ID
                            # that it is inside the pdb file, so it is needed to find the real one.
                            chain_id = check_real_chainID(pdb_id, chain_id, uniprot_id)
                            if chain_id:

                                # Write down the data
                                full_pdb_id = f"{pdb_id}_{chain_id}"
                                pdb_path = f"data/{full_pdb_id}.pdb"
                                current_pdb_path = Path(
                                    f"output_PDB/{pdb_id}_renum.pdb"
                                )

                                # Residues data will be useful for finding interacting ligands.
                                residues = get_pdb_residues(current_pdb_path, chain_id)
                                if residues == None:
                                    logging.warning(
                                        f"{full_pdb_id} was discarded due to wrong renumbering"
                                    )
                                else:
                                    # Check if total of residues is bigger than the cutoff
                                    if len(residues) >= AA_CUTOFF:
                                        # Write data
                                        logging.info(f"Successful {full_pdb_id}!")
                                        pdb_f[full_pdb_id] = {}
                                        pdb_f[full_pdb_id]["pdb_path"] = pdb_path
                                        pdb_f[full_pdb_id]["residues"] = residues
                                        pdb_f[full_pdb_id]["positives"] = []
                                    else:
                                        logging.warning(
                                            f"{pdb_id} was discarded due to size <= {AA_CUTOFF} residues"
                                        )

                            else:
                                logging.warning(
                                    f"{pdb_id} was discarded due to wrong chain ID"
                                )
                        else:
                            logging.warning(
                                f"{pdb_id} was discarded due to wrong renumbering"
                            )

                    # If not the first one.
                    else:
                        # Try downloading and renumbering pdb_id
                        renumbered, WRONG_PDBS_DB = run_PDBrenum(pdb_id, WRONG_PDBS_DB)
                        if renumbered:

                            # The chain ID that appears at the best-structures API sometimes is not the real ID
                            # that it is in the pdb file, so it needs to find the real one.
                            chain_id = check_real_chainID(pdb_id, chain_id, uniprot_id)
                            if chain_id:
                                full_pdb_id = f"{pdb_id}_{chain_id}"
                                pdb_path = f"data/{full_pdb_id}.pdb"
                                current_pdb_path = Path(
                                    f"output_PDB/{pdb_id}_renum.pdb"
                                )

                                # Residues data will be useful for finding interacting ligands.
                                residues = get_pdb_residues(current_pdb_path, chain_id)

                                if residues == None:
                                    # The renumbering didn't go well.
                                    logging.warning(
                                        f"{full_pdb_id} was discarded due to wrong renumbering"
                                    )

                                else:
                                    # Check if total od residues is bigger than the cutoff
                                    if len(residues) >= AA_CUTOFF:
                                        # Check if the residues in this new structure are different to the rest of the structures.
                                        for pdb_structure in pdb_f.keys():
                                            different = False
                                            for resi in residues:
                                                if (
                                                    resi
                                                    not in pdb_f[pdb_structure][
                                                        "residues"
                                                    ]
                                                ):
                                                    different = True
                                            if not different:
                                                logging.warning(
                                                    f"{full_pdb_id} was discarded for not bringing a different coverage of the structure"
                                                )
                                                break
                                        else:
                                            # Write down the data
                                            logging.info(f"Successful {full_pdb_id}!")
                                            pdb_f[full_pdb_id] = {}
                                            pdb_f[full_pdb_id]["pdb_path"] = pdb_path
                                            pdb_f[full_pdb_id]["residues"] = residues
                                            pdb_f[full_pdb_id]["positives"] = []

                                    else:
                                        logging.warning(
                                            f"{pdb_id} was discarded due to size <= {AA_CUTOFF} residues"
                                        )

                            else:
                                logging.warning(
                                    f"{full_pdb_id} was discarded due to wrong chain ID"
                                )
                        else:
                            logging.warning(
                                f"{pdb_id} was discarded due to wrong renumbering"
                            )
            else:
                logging.info(f"{pdb_id} was already processed")

            n_pdb += 1
            if n_pdb == max_n_pdb:
                logging.info("Last pdb structure found analyzed")

        # Check if any structure has passed the filters.
        if len(pdb_f) == 0:
            pdb_f = None

    return pdb_f, WRONG_PDBS_DB


def check_multimers_receptor(ligand_id, receptor_id, rec_int_pdbs):
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

                if partner_id == receptor_id:  # Receptor-partner interaction
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
    confirmed_pdbs = []
    for pdb in candidate_pdbs:
        pdb_url = f"{PDB_MOLECULES_URL}/{pdb}"

        # The request might fail.
        try:
            pdb_data = make_request(pdb_url, None)
        except Exception as e:
            logging.debug(f"Could not make PDB request for {pdb}, {e}")
            return None

        if pdb_data:
            n_proteins = 0
            for entity in pdb_data[pdb]:
                if "molecule_type" in entity.keys():
                    if "polypeptide" in entity["molecule_type"]:
                        if "in_struct_asyms" in entity.keys():
                            n_proteins += len(entity["in_struct_asyms"])

            if n_proteins == 2:
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


def get_interacting_partners(
    uniprot_id, pdbf, prodigy_results_pdb, outHM, t2m, cpc, cbi
):
    "Use the interface API to find partners for the given Uniprot ID."

    logging.info(f"Inspecting interacting partners for {uniprot_id}")

    # Inputs:
    #   uniprot_id: ID of the molecule
    #   pdbf: dictionary that stores previous information for the molecule
    #   cov_cutoff: coverage threshold
    #   prodigy_results_pdb: dictionary that stores information of previous ProdigyCrystal runs that tell if the interactions are
    #                        biological or not
    #   outHM: flag to exclude homomultimers (i.e. only 1:1 PPI)
    #   t2m: flag to consider only PDBs in which there are only 2 protein molecules
    #   cpc: flag to "check partners coverage", i.e. check that in the PDB where the interaction has been seen, the coverage
    #        of both receptor and ligand is greater than the established threshold.
    #   cbi: flah to "check biological interaction", i.e. run ProdigyCrystal to know if the interaction is biological or not

    # Get info from the API
    url = f"{INTERFACE_URL}/{uniprot_id}"
    # The request might fail
    try:
        interface_data = make_request(url, None)
    except Exception as e:
        logging.debug(f"Could not make InterfaceResidues request for {uniprot_id}, {e}")
        return False

    if interface_data and len(interface_data) != 0:  # Check if the entry is empty
        if len(interface_data[uniprot_id]) == 0:  # Check if the entry is empty.
            logging.debug(f"InterfaceResidues for {uniprot_id} was empty")
            pdbf = False
        else:  # There is information

            rec_len = len(
                interface_data[uniprot_id]["sequence"]
            )  # Lenght of the receptor according to Uniprot
            # Iterate every interacting partner
            for partner in interface_data[uniprot_id]["data"]:
                partner_id = partner["accession"]  # Uniprot_id of the partner
                logging.info(
                    f"Evaluating interaction between {uniprot_id} and {partner_id}"
                )

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
                        logging.info(f"Checking receptor multimers for {uniprot_id}")
                        confirmed_int_pdbs = check_multimers_receptor(
                            partner_id, uniprot_id, included_int_pdbs
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
                                        f"{uniprot_id}_{partner_id} discarded due to no pdb with only 2 proteins"
                                    )
                                    continue
                        else:
                            logging.warning(
                                f"{uniprot_id}_{partner_id} discarded due to receptor only interacts as multimer"
                            )
                            continue
                    else:
                        logging.warning(
                            f"{uniprot_id}_{partner_id} discarded due to ligand only interacts as multimer"
                        )
                        continue
                else:
                    confirmed_int_pdbs = partner["allPDBEntries"]

                if (
                    cpc
                ):  # Check if coverage of the proteins in the interacting PDBs is over the established threshold
                    logging.info(
                        f"Checking partners coverage for {uniprot_id}_{partner_id}"
                    )
                    confirmed_int_pdbs = check_partners_coverage(
                        uniprot_id, rec_len, partner_id, confirmed_int_pdbs
                    )
                    if confirmed_int_pdbs == None:
                        logging.warning(
                            f"{uniprot_id}_{partner_id} discarded due to coverage under cutoff"
                        )
                        continue

                if (
                    cbi
                ):  # Check: if the interactions for the confirmed pdbs are biological by Prodigy Crystal
                    logging.info(
                        f"Checking biological interaction for {uniprot_id}-{partner_id}"
                    )
                    confirmed_int_pdbs, prodigy_results_pdb = validate_interaction(
                        uniprot_id, partner_id, confirmed_int_pdbs, prodigy_results_pdb
                    )
                    if confirmed_int_pdbs == None:
                        logging.warning(
                            f"{partner_id} discarded due to no biological interaction"
                        )
                        continue

                # Get interface residues of the single PPI interaction
                logging.info(
                    f"Checking matching interfaces in the structure of {uniprot_id}"
                )
                interface = get_interface(partner, confirmed_int_pdbs)
                if interface != None:  # There are valid PDBs and valid residues

                    for pdb_structure in pdbf.keys():
                        if [partner_id] not in pdbf[pdb_structure]["positives"]:
                            # Check the interface of each PDB separatedly
                            for int_pdb in interface.keys():

                                # Check if interface residues are present in the receptor structure
                                for residue in interface[int_pdb]:
                                    if residue not in pdbf[pdb_structure]["residues"]:
                                        logging.warning(
                                            f"Interaction between partner {partner_id} and receptor {uniprot_id}-{pdb_structure} \
                                            due to interface not in pdb structure"
                                        )
                                        break

                                else:  # The interface is valid
                                    logging.info(
                                        f"Successful interacting partner {partner_id} with {uniprot_id}-{pdb_structure}"
                                    )
                                    pdbf[pdb_structure]["positives"] += [partner_id]
                                    break

                else:
                    logging.warning(
                        f"{partner_id} discarded due to interface with {uniprot_id}_{pdb_structure} not valid"
                    )

            # Check if the pdb structures have valid interacting partners, and if not removing them
            del_pdbs = []
            for pdb_structure in pdbf.keys():
                if len(pdbf[pdb_structure]["positives"]) == 0:  # Not useful
                    del_pdbs += [pdb_structure]
            if len(del_pdbs) != 0:  # Check if it is needed to delete any structure
                for del_pdb in del_pdbs:
                    del pdbf[del_pdb]
            if len(pdbf) == 0:
                pdbf = None
    else:
        pdbf = False

    return pdbf, prodigy_results_pdb


def validate_interaction(prot_a, prot_b, pre_sel_pdbs, prodigy_results_pdb):
    """Check if the interaction between two proteins is biological."""

    biological_pdbs = []

    url_a = f"{INTERFACE_URL}/{prot_a}"
    try:
        interface_data_a = make_request(url_a, None)
    except Exception as e:
        logging.debug(f"Could not make InterfaceResidues request for {prot_a}, {e}")
        return None

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
            if pre_sel_pdb in prodigy_results_pdb[protein_pair].keys():
                not_redo_pdbs += [pre_sel_pdb]
        # Avoid redoing already done pdbs.
        for not_redo_pdb in not_redo_pdbs:
            if not_redo_pdb in pre_sel_pdbs:
                pre_sel_pdbs.remove(not_redo_pdb)

    # Obtain a dictionary with pdb_id : chain_id for both proteins
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

    # assert pdb_chain_a.keys() == pdb_chain_b.keys()
    # """It gives error because the API contains some mistakes, so it's better to correct them beforehand"""

    for pdb in pdb_chain_a.keys():
        if pdb in pdb_chain_b.keys():
            checked_pdb_chain_a[pdb] = pdb_chain_a[pdb]
            checked_pdb_chain_b[pdb] = pdb_chain_b[pdb]
    assert checked_pdb_chain_a.keys() == checked_pdb_chain_b.keys()

    # Run Prodigy
    for pdb in checked_pdb_chain_a:
        renum_pdb_f = Path(f"output_PDB/{pdb}_renum.pdb")
        if (
            renum_pdb_f.exists()
        ):  # Check if the pdb has already been downloaded and renumbered
            pdb_f = renum_pdb_f
        else:  # Check if the pdb has already been downloaded
            pdb_f = Path(f"output_PDB/{pdb}.pdb")
        if not pdb_f.exists():
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
        # Those pdbs haven't been evaluated as biological before, so they are not
        if pdb not in prodigy_results_pdb[protein_pair].keys():
            prodigy_results_pdb[protein_pair][pdb] = "NOT biological"

    for processsed_pdb in prodigy_results_pdb[protein_pair].keys():
        if prodigy_results_pdb[protein_pair][processsed_pdb] == "biological":
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

    negatome_list = list(negatome_dic.keys())
    return negatome_list, negatome_dic


# ====================================================================================================
# Helper functions


def prepare_pdb(pdb_path):
    """Prepare the pdb file for the pdb id and chain in the path where it will be stored"""

    # pdb_path = Path(f"data/{pdb_id}_{chain}.pdb")
    pdb_filename = str(Path(pdb_path).stem)
    pdb_id = pdb_filename.split("_")[-2]
    chain = pdb_filename.split("_")[-1]
    renum_pdb = f"output_PDB/{pdb_id}_renum.pdb"  # Downloaded and renumbered file

    logging.info(f"Preparing pdb for structure {pdb_id}_{chain}")
    # Check if already processed
    if Path(pdb_path).exists():
        logging.debug(f"{pdb_path} already exists, skipping preparation")
        return
    else:
        # Split models if there are in the pdb structure
        prev_cwd = Path.cwd()
        os.chdir(f"output_PDB")
        renum_file_name = Path(renum_pdb).name
        subprocess.run(["pdb_splitmodel", renum_file_name], check=True)
        os.chdir(prev_cwd)
        split_pdb_name = renum_pdb.rstrip(".pdb") + "_1.pdb"
        if Path(split_pdb_name).exists():
            renum_pdb = split_pdb_name

        # Get files ready to be run in HADDOCK
        cmd = f"grep 'ATOM' {renum_pdb} | pdb_selchain -{chain} | pdb_selaltloc | pdb_tidy"

        p = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )

        out, err = p.communicate()

        # Create the new file
        with open(pdb_path, "w") as fh:
            fh.write(out.decode("utf-8"))

        return


# =====================================================================================


if __name__ == "__main__":
    """Create the Interactome DB."""

    logging.info("============================")

    # Process the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "create_benchmark",
        help="Introduce location to store the data to create the benchmark",
    )
    parser.add_argument(
        "-nrmd",
        action="store_true",
        help="No remove positive/negative partners duplicates",
    )
    parser.add_argument(
        "-onia", action="store_true", help="Off negative interface API check"
    )
    parser.add_argument(
        "-outHM",
        action="store_true",
        help="Exclude homomultimers interactions from the dataset",
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
    parser.add_argument(
        "prodigy_results_pdb",
        nargs="?",
        help="Introduce prodigy_results_pdb.txt file if possible",
    )
    args = parser.parse_args()
    bm_directory = args.create_benchmark
    nrmd = args.nrmd
    onia = args.onia
    outHM = args.outHM
    t2m = args.t2m
    cpc = args.cpc
    cbi = args.cbi

    if t2m and not outHM:
        logging.info(
            f"For a faster running of the script in -t2m mode, please consider also using the -outHM mode."
        )
        exit()

    prodigy_results_pdb_file = args.prodigy_results_pdb

    Path(bm_directory).mkdir(exist_ok=True)

    prodigy_results_pdb = {}

    # Running the script is a slow process with the big bottleneck of running Prodigy Crystal
    # (to check if the interactions are biological). To make it faster, the script can re-use the
    # results of running it before, stored in a file called prodigy_results_pdb.txt.
    if prodigy_results_pdb_file is None:
        logging.info(f"No prodigy_results_pdb.txt file introduced")
    else:
        if Path(prodigy_results_pdb_file).exists():
            logging.info(
                "<<<<<< Prodigy_results_pdb.txt file found, reusing it. >>>>>>"
            )
            read_prodigy_results_pdb = open(prodigy_results_pdb_file, "r").read()
            prodigy_results_pdb = ast.literal_eval(read_prodigy_results_pdb)
        else:
            logging.info("<<<<<< Prodigy_results_pdb.txt file not found. >>>>>>")

    os.chdir(bm_directory)
    Path(f"data").mkdir(exist_ok=True)  # Folder to store processed pdbs

    # Different dictionaries to store processed information.
    interactome_db = (
        {}
    )  # Confirmed groups receptor - positive_ligand(s) - negative_ligand(s)
    all_molecules_pdb = {}  # pdbf of all processed molecules
    interacting_partners = (
        {}
    )  # pdbf after get_interacting_partners function (receptors and positive ligands)
    negative_partners = {}  # Negative partners of a receptor and their pdb structures
    positive_negative_duplicates = (
        []
    )  # Receptor-partner of ligands that appear both as negative and positive

    output_csv = [
        ["Receptor", "Receptor_pdb", "Ligand", "Ligand_pdb", "Interaction_type"]
    ]  # Final valid data to build the benchmark later

    # Exclusion reasons and storage variables.
    exclusion = [["Receptor_ID", "Exclusion reason"]]
    n_no_receptor_structure = 0
    n_no_interface_API = 0
    n_no_negative_partners = 0
    n_no_positive_partners = 0

    # Process negatome database
    negatome_list, negatome_dic = load_negatome_data(NEGATOME_URL)
    total_receptors = len(negatome_list)

    # Iterate the receptors in the negatome database
    for i, receptor_uniprot_id in enumerate(negatome_list, start=1):
        logging.info(
            f"========== {receptor_uniprot_id} ========== ({i}/{total_receptors})"
        )
        interactome_db[receptor_uniprot_id] = {}

        # Check if the pdb structures for the receptor have been processed before
        if receptor_uniprot_id not in all_molecules_pdb.keys():

            # Find pdb structures for the receptor
            receptor_pdbf, WRONG_PDBS_DB = get_bestpdb(
                receptor_uniprot_id, WRONG_PDBS_DB
            )

            # E.g:
            # receptor_pdb_f = {
            #           'a1b2': {
            #                   'pdb_path': 'data/a1b2.pdb',
            #                   'start': 5,
            #                   'end': 167,
            #                   'positives': []
            #                       },
            #           'c3d4': {(...)}, ...
            #                   }

            all_molecules_pdb[receptor_uniprot_id] = receptor_pdbf
        else:
            logging.info(f"{receptor_uniprot_id} structures were processed before")
            receptor_pdbf = all_molecules_pdb[receptor_uniprot_id]

        if receptor_pdbf == None:
            logging.info(
                f"No acceptable structures found for {receptor_uniprot_id}, discarding receptor."
            )
            exclusion += [[receptor_uniprot_id, "No receptor structure"]]
            n_no_receptor_structure += 1

        else:
            logging.info(f"{receptor_uniprot_id} has valid structure(s)")
            # There are valid structures for the receptor.

            # Check if the receptor has already been processed for interacting partners
            if receptor_uniprot_id not in interacting_partners.keys():
                # Get interacting partners
                receptor_pdbf, prodigy_results_pdb = get_interacting_partners(
                    receptor_uniprot_id,
                    receptor_pdbf,
                    prodigy_results_pdb,
                    outHM,
                    t2m,
                    cpc,
                    cbi,
                )
                interacting_partners[receptor_uniprot_id] = receptor_pdbf
            else:
                logging.info(
                    f"{receptor_uniprot_id} interacting partners were processed before"
                )
                receptor_pdbf = interacting_partners[receptor_uniprot_id]

            if receptor_pdbf == False:  # Discard receptor
                logging.warning(
                    f"Receptor {receptor_uniprot_id} discarded due to no interface API"
                )
                exclusion += [[receptor_uniprot_id, "No interface API"]]
                n_no_interface_API += 1

            elif receptor_pdbf == None:  # Discard receptor
                logging.warning(
                    f"Receptor {receptor_uniprot_id} discarded due to no positive partners"
                )
                exclusion += [[receptor_uniprot_id, "No positive partners"]]
                n_no_positive_partners += 1

            else:
                logging.info(
                    "Interacting partners found, proceding to test of duplicity"
                )

                # Test of duplicity
                del_pos_neg_duplicates = []
                for receptor_pdb in receptor_pdbf.keys():
                    for positive in receptor_pdbf[receptor_pdb]["positives"]:
                        if positive in negatome_dic[receptor_uniprot_id]:
                            logging.warning(
                                f"{positive} appears as positive and negative partner, it's a DUPLICATE."
                            )
                            positive_negative_duplicates += [
                                f"{receptor_uniprot_id}_{positive}"
                            ]
                            del_pos_neg_duplicates += [positive]

                # Remove duplicates
                if len(del_pos_neg_duplicates) > 0:
                    if (
                        not nrmd
                    ):  # no remove duplicates (accept them as negative partners)
                        logging.info("Removing duplicates...")
                        for duplicate in del_pos_neg_duplicates:
                            for receptor_pdb in receptor_pdbf.keys():
                                if (
                                    duplicate
                                    in receptor_pdbf[receptor_pdb]["positives"]
                                ):
                                    receptor_pdbf[receptor_pdb]["positives"].remove(
                                        duplicate
                                    )
                            if duplicate in negatome_dic[receptor_uniprot_id]:
                                negatome_dic[receptor_uniprot_id].remove(duplicate)
                    else:
                        logging.info("Removing duplicates from positive partners...")
                        for duplicate in del_pos_neg_duplicates:
                            for receptor_pdb in receptor_pdbf.keys():
                                if (
                                    duplicate
                                    in receptor_pdbf[receptor_pdb]["positives"]
                                ):
                                    receptor_pdbf[receptor_pdb]["positives"].remove(
                                        duplicate
                                    )
                else:
                    logging.info("No positive/negative partners duplicates found")

                # At this point it is better to find the pdb structures of the negative partners first,
                # because it is a less demanding operation than the positive ones in case it is needed
                # to discard the receptor
                negative_matches = []

                # Iterate every negative partner for the receptor
                for negative_uniprot_id in negatome_dic[receptor_uniprot_id]:
                    logging.info(
                        f"--- Inspecting negative candidate {negative_uniprot_id}"
                    )

                    # Check if the negative has already been processed for other receptors.
                    if negative_uniprot_id not in negative_partners.keys():

                        # Check if the pdb structures for the negative have been processed before
                        if negative_uniprot_id not in all_molecules_pdb.keys():

                            # Find pdb structures for the negative
                            negative_pdbf, WRONG_PDBS_DB = get_bestpdb(
                                negative_uniprot_id, WRONG_PDBS_DB
                            )
                            all_molecules_pdb[negative_uniprot_id] = negative_pdbf
                        else:
                            logging.info(
                                f"{negative_uniprot_id} structures were processed before"
                            )
                            negative_pdbf = all_molecules_pdb[negative_uniprot_id]

                        if negative_pdbf == None:  # Discard the negative
                            negative_partners[negative_uniprot_id] = None
                            logging.warning(
                                f"{negative_uniprot_id} discarded due to no good structures"
                            )

                        else:

                            # Check if there is interface API and interacting partners for the negative
                            if not onia:  # off negative interface API check
                                logging.info(
                                    f"Checking interface API and interacting partners for negative {negative_uniprot_id}"
                                )
                                # Check if negative has been processed before for interacting partners
                                if (
                                    negative_uniprot_id
                                    not in interacting_partners.keys()
                                ):
                                    # Get interacting partners
                                    (
                                        negative_pdbf,
                                        prodigy_results_pdb,
                                    ) = get_interacting_partners(
                                        negative_uniprot_id,
                                        negative_pdbf,
                                        prodigy_results_pdb,
                                        outHM,
                                        t2m,
                                        cpc,
                                        cbi,
                                    )
                                    interacting_partners[
                                        negative_uniprot_id
                                    ] = negative_pdbf
                                else:
                                    logging.info(
                                        f"{receptor_uniprot_id} interacting partners were processed before"
                                    )
                                    negative_pdbf = interacting_partners[
                                        negative_uniprot_id
                                    ]

                                if negative_pdbf == False:  # Discard negative
                                    logging.warning(
                                        f"Negative {negative_uniprot_id} discarded due to no interface API"
                                    )
                                    negative_partners[negative_uniprot_id] = None
                                    continue

                                elif negative_pdbf == None:  # Discard negative
                                    logging.warning(
                                        f"Negative {negative_uniprot_id} discarded due to no interacting partners"
                                    )
                                    negative_partners[negative_uniprot_id] = None
                                    continue

                            # Negative is valid, write down data
                            logging.info(
                                f"Negative {negative_uniprot_id} is a valid partner!"
                            )
                            negative_matches += [negative_uniprot_id]
                            negative_partners[negative_uniprot_id] = negative_pdbf

                    else:
                        logging.info(
                            f"{negative_uniprot_id} was fully processed before"
                        )
                        if (
                            negative_partners[negative_uniprot_id] == False
                        ):  # Discard negative
                            logging.warning(
                                f"Negative {negative_uniprot_id} discarded due to no interface API"
                            )
                        elif (
                            negative_partners[negative_uniprot_id] == None
                        ):  # Discard negative
                            logging.warning(
                                f"Negative {negative_uniprot_id} discarded due to no interacting partners"
                            )
                        else:
                            # Negative is valid
                            logging.info(
                                f"Negative {negative_uniprot_id} is a valid partner!"
                            )
                            negative_matches += [negative_uniprot_id]

                # Check if there are valid negative partners.
                if len(negative_matches) > 0:
                    # If so, continue with the positives:
                    logging.info(
                        f"Negative partners found, looking for positive partners"
                    )
                    positive_partners = {}

                    # Iterate every receptor pdb structure
                    for receptor_pdb in receptor_pdbf.keys():
                        logging.info(
                            f"Inspecting positive partners for {receptor_uniprot_id}_{receptor_pdb}"
                        )
                        positive_matches = []
                        # Iterate every positive partner.
                        for positive_uniprot_id in receptor_pdbf[receptor_pdb][
                            "positives"
                        ]:
                            logging.info(
                                f"++ Inspecting positive candidate {positive_uniprot_id}"
                            )

                            # Check if the pdb structures for the positive have been processed before
                            if positive_uniprot_id not in all_molecules_pdb.keys():
                                # Find pdb structures for positive
                                positive_pdbf, WRONG_PDBS_DB = get_bestpdb(
                                    positive_uniprot_id, WRONG_PDBS_DB
                                )
                                all_molecules_pdb[positive_uniprot_id] = positive_pdbf
                            else:
                                logging.info(
                                    f"{positive_uniprot_id} structures were processed before"
                                )
                                positive_pdbf = all_molecules_pdb[positive_uniprot_id]

                            if positive_pdbf == None:
                                logging.warning(
                                    f"Positive {positive_uniprot_id} discarded due to no good structures"
                                )

                            else:
                                logging.info(
                                    f"Acceptable found pdb structure(s) for positive {positive_uniprot_id}, proceding to test of reciprocity."
                                )

                                # Test of reciprocity:
                                # Check if the receptor is also a valid positive partner for the positive partner.

                                # Check if the positive ligand has already been processed for interacting partners.
                                if (
                                    positive_uniprot_id
                                    not in interacting_partners.keys()
                                ):
                                    (
                                        positive_pdbf,
                                        prodigy_results_pdb,
                                    ) = get_interacting_partners(
                                        positive_uniprot_id,
                                        positive_pdbf,
                                        prodigy_results_pdb,
                                        outHM,
                                        t2m,
                                        cpc,
                                        cbi,
                                    )
                                    interacting_partners[
                                        positive_uniprot_id
                                    ] = positive_pdbf
                                else:
                                    logging.info(
                                        f"{receptor_uniprot_id} interacting partners were processed before"
                                    )
                                    positive_pdbf = interacting_partners[
                                        positive_uniprot_id
                                    ]

                                if positive_pdbf == False:
                                    logging.warning(
                                        f"Positive {positive_uniprot_id} discarded due to no interface API"
                                    )

                                elif positive_pdbf == None:
                                    logging.warning(
                                        f"Positive {positive_uniprot_id} discarded due to no positive partners"
                                    )

                                else:
                                    remove_pos_pdb = []
                                    for positive_pdb in positive_pdbf.keys():
                                        logging.info(
                                            f"Evaluating {positive_uniprot_id}_{positive_pdb} reciprocity test"
                                        )
                                        if (
                                            receptor_uniprot_id
                                            in positive_pdbf[positive_pdb]["positives"]
                                        ):
                                            # Successful reciprocity test.
                                            logging.info("Successful reciprocity test")
                                            if (
                                                positive_uniprot_id
                                                not in positive_matches
                                            ):
                                                logging.info(
                                                    f"++ It's a match between {receptor_uniprot_id} and {positive_uniprot_id}!"
                                                )
                                                positive_matches += [
                                                    positive_uniprot_id
                                                ]
                                            else:
                                                logging.info(
                                                    f"{receptor_uniprot_id}_{positive_uniprot_id} was already a match!"
                                                )
                                        else:
                                            logging.info(
                                                f"Failed reciprocity test, discarding positive pdb structure"
                                            )

                                    if positive_uniprot_id in positive_matches:
                                        if len(remove_pos_pdb) > 0:
                                            for pos_pdb in remove_pos_pdb:
                                                del positive_pdbf[pos_pdb]
                                        positive_partners[
                                            positive_uniprot_id
                                        ] = positive_pdbf

                                    else:
                                        logging.info(
                                            f"{positive_uniprot_id} is not a valid positive partner, discarding it."
                                        )

                        if len(positive_matches) > 0:
                            logging.info(
                                f"Positive and negative partners found for {receptor_uniprot_id}-{receptor_pdb}"
                            )

                            # Write database and prepare pdb files.
                            interactome_db[receptor_uniprot_id][receptor_pdb] = {
                                "positives": {},
                                "negatives": {},
                                "residues": receptor_pdbf[receptor_pdb]["residues"],
                                "pdb_path": None,
                            }
                            receptor_pdb_file = receptor_pdbf[receptor_pdb]["pdb_path"]
                            prepare_pdb(receptor_pdb_file)
                            interactome_db[receptor_uniprot_id][receptor_pdb][
                                "pdb_path"
                            ] = receptor_pdb_file
                            for negative in negative_matches:
                                interactome_db[receptor_uniprot_id][receptor_pdb][
                                    "negatives"
                                ][negative] = negative_partners[negative]
                                for negative_pdb in negative_partners[negative].keys():
                                    negative_pdb_path = negative_partners[negative][
                                        negative_pdb
                                    ]["pdb_path"]
                                    prepare_pdb(negative_pdb_path)
                                    output_csv += [
                                        [
                                            receptor_uniprot_id,
                                            receptor_pdb_file,
                                            negative,
                                            negative_pdb_path,
                                            "negative",
                                        ]
                                    ]
                            for positive in positive_matches:
                                interactome_db[receptor_uniprot_id][receptor_pdb][
                                    "positives"
                                ][positive] = positive_partners[positive]
                                for positive_pdb in positive_partners[positive].keys():
                                    positive_pdb_path = positive_partners[positive][
                                        positive_pdb
                                    ]["pdb_path"]
                                    prepare_pdb(positive_pdb_path)
                                    output_csv += [
                                        [
                                            receptor_uniprot_id,
                                            receptor_pdb_file,
                                            positive,
                                            positive_pdb_path,
                                            "positive",
                                        ]
                                    ]

                            """ interactome_db = { 
                                 	Receptor1 {
                                        'rec_pdb1': {
                                            'positives': {
                                                Positive1: {
                                                    pos_pdb1:{
                                                        'positives':[],
                                                        'residues': [1,2,3,4, ...]
                                                        'pdb_path': 'data/a1b2G.pdb'
                                                            },
                                                    pos_pdb2: {}, ...
                                                        },
                                                Positive2: {}, ...
                                                    },
                                                
                                            'negatives': {
                                                Negative1:
                                                    neg_pdb1:{
                                                        'positives':[],
                                                        'residues':[7,8,9,...]
                                                        'pdb_path': 'data/red4A.pdb'
                                                            },
                                                    neg_pdb2: {}, ...
                                                        }
                                                Negative2: {}, ...
                                                    },
                                            'residues':[34,35,36,...],
                                            'pdb_path': 'data/a1b2F.pdb'
                                            },

                                        'rec_pdb2': {...}, ...
                                        },
                                    Receptor2: {}, ...
                                    } """
                        else:
                            logging.info(
                                "No acceptable positive partners, discarding receptor"
                            )
                            exclusion += [[receptor_uniprot_id, "No positive partners"]]
                            n_no_positive_partners += 1
                else:
                    logging.info("No acceptable negative partners, discarding receptor")
                    exclusion += [[receptor_uniprot_id, "No negatives partners"]]
                    n_no_negative_partners += 1

        # Remove not relevant entries
        if len(interactome_db[receptor_uniprot_id]) == 0:
            del interactome_db[receptor_uniprot_id]

        logging.info("============================")

    # Write all the files obtained as output
    with open("interacting_db.csv", "w") as fh:
        writer = csv.writer(fh)
        for element in output_csv:
            writer.writerow(element)
        fh.close()

    with open("interaction_class.tsv", "w") as fh:
        for element in PDB_DB:
            fh.write("\t".join(map(str, element)) + os.linesep)
        fh.close()

    with open("exclusion_reason.csv", "w") as fh:
        writer = csv.writer(fh)
        for receptor in exclusion:
            writer.writerow(receptor)
        fh.close()

    with open("prodigy_results_pdb.txt", "w") as fh:
        fh.write(str(prodigy_results_pdb))
        fh.close()

    with open("interactome_db.txt", "w") as fh:
        fh.write(str(interactome_db))
        fh.close()

    with open("positive_negative_duplicates.txt", "w") as fh:
        fh.write(str(positive_negative_duplicates))
        fh.close()

    # Plot Receptors Exclusion Reason graph
    exclusion_labels = (
        "No receptor structure",
        "No interface API",
        "No negative partners",
        "No positive partners",
    )
    sizes = [
        n_no_receptor_structure,
        n_no_interface_API,
        n_no_negative_partners,
        n_no_positive_partners,
    ]
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=exclusion_labels, autopct="%1.1f%%", startangle=90)
    ax1.axis("equal")
    plt.savefig("exclusion_reason.png")
    plt.close()

    excluded_receptors = (
        n_no_receptor_structure
        + n_no_interface_API
        + n_no_negative_partners
        + n_no_positive_partners
    )

    # Report numerical results
    print()
    print("Dataset creation finished:")
    print(
        f"The number of excluded receptors is {excluded_receptors}: {n_no_receptor_structure} due to no receptor structure, \
    {n_no_interface_API} due to no interface API, {n_no_negative_partners} due to no negative partners, and \
    {n_no_positive_partners} due to no positive partners"
    )
    print()
    print(
        f"The total number of problematic pdbs that PDBRenum couldn't run was {len(WRONG_PDBS_DB)}"
    )
    print(
        f"The number of positive/negative duplicates was {len(positive_negative_duplicates)}"
    )
