import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("benchmark", help="Introduce benchmark location")
args = parser.parse_args()
bm_directory = Path(args.benchmark)

n_pos_interactions = 0
n_neg_interactions = 0
unique_proteins = []
unique_structures = []
unique_receptors = []
unique_receptor_structures = []
unique_ligands = []
unique_ligand_structures = []

unique_positive_ligands = []
unique_positive_ligand_structures = []
unique_negative_ligands = []
unique_negative_ligand_structures = []

for folder in bm_directory.iterdir():
    folder_name = str(folder.stem)
    molecules = folder_name.split("_")
    if len(molecules) == 3 and folder.is_dir():

        receptor = molecules[0]
        if str(receptor[-1]).islower():
            real_receptor = receptor[:-1]
        else:
            real_receptor = receptor

        if real_receptor not in unique_receptors:
            unique_receptors += [real_receptor]
        if receptor not in unique_receptor_structures:
            unique_receptor_structures += [receptor]

        ligand = molecules[1]
        if str(ligand[-1]).islower():
            real_ligand = ligand[:-1]
        else:
            real_ligand = ligand

        if real_ligand not in unique_ligands:
            unique_ligands += [real_ligand]
        if ligand not in unique_ligand_structures:
            unique_ligand_structures += [ligand]

        if molecules[2].split("-")[0] == "positive":
            n_pos_interactions += 1
            if real_ligand not in unique_positive_ligands:
                unique_positive_ligands += [real_ligand]
            if ligand not in unique_positive_ligand_structures:
                unique_positive_ligand_structures += [ligand]
        elif molecules[2].split("-")[0] == "negative":
            n_neg_interactions += 1
            if real_ligand not in unique_negative_ligands:
                unique_negative_ligands += [real_ligand]
            if ligand not in unique_negative_ligand_structures:
                unique_negative_ligand_structures += [ligand]

        for protein in molecules[:2]:
            if str(protein[-1]).islower():
                real_protein = protein[:-1]
            else:
                real_protein = protein

            if real_protein not in unique_proteins:
                unique_proteins += [real_protein]
            if protein not in unique_structures:
                unique_structures += [protein]

print()

if len(unique_receptor_structures) > 0 or len(unique_ligand_structures) > 0:
    print(
        "There are multiple structures for some receptors, please consider using generate_interactome_benchmark.py \
in the -outMS mode"
    )
    print()

print(f"Number of different proteins: {len(unique_proteins)}")
print(f"Number of different protein structures {len(unique_structures)}")
print(f"Number of different receptors: {len(unique_receptors)}")
print(f"Number of different receptor structures: {len(unique_receptor_structures)}")
print(f"Number of different ligands: {len(unique_ligands)}")
print(f"Number of different ligand structures: {len(unique_ligand_structures)}")
print()
print(f"Number of different positive ligands: {len(unique_positive_ligands)}")
print(
    f"Number of different positive ligand structures: {len(unique_positive_ligand_structures)}"
)
print(f"Number of different negative ligands: {len(unique_negative_ligands)}")
print(
    f"Number of different negative ligand structures: {len(unique_negative_ligand_structures)}"
)
print()
print(f"Positive interactions: {n_pos_interactions}")
print(f"Negative interactions: {n_neg_interactions}")
print(
    f"Total of interactions / runs: {int(n_pos_interactions)+int(n_neg_interactions)}"
)
