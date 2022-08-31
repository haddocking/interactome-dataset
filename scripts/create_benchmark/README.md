# create_interactome_dataset.py

""" Requirements: """

- Create a conda environment first for running it correctly:
  - (base) $ git clone https://github.com/Faezov/PDBrenum.git
  - (base) $ cd PDBrenum
  - (base) $ conda create -n PDBrenum python=3.6 numpy=1.17 pandas=0.25.1 biopython=1.76 tqdm=4.36.1 ipython=7.8.0 requests=2.25.1 lxml=4.6.2
  - (base) $ conda activate PDBrenum
- Install pdb-tools and matplotlib
  - (PDBrenum) $ pip install pdb-tools
  - (PDBrenum) $ pip install matplotlib
- Install Prodigy Crystal:
  - (PDBrenum) git clone https://github.com/haddocking/prodigy-cryst.git
  - (PDBrenum) $ cd prodigy-cryst
  - (PDBrenum) $ pip install scipy
  - (PDBrenum) $ pip instal cython
  - (PDBrenum) $ pip install -r requirements
- Remember to change the location of PDBrenum and ProdigyCrystal in the script to match your own location.

""" Execution """

- This script is the first step to create a so-called interactome dataset. This is a dataset which contains positive and negative
  interactions organized in a way that the same receptor protein has a positive partner (interaction found in nature)
  and a negative one (interaction proven that does not happen in nature). In that way, the group receptor-positive-negative can
  be used to evaluate the differences between a real and a non-real interaction after any PPI prediction method. It can also be used
  to train AI models.

- Since negative data is the limiting factor and the most difficult to found in literature, it is the starating point of the script.
  It will process the Negatome2.0 database to extract pairs of receptor-negative partners. Then it will try to find good pdb
  structures under the established thresholds:

  - COV_CUTOFF: Coverage of the structure over which the structure is accepted. E.g: 0.5 = 50% of the structure
  - RES_CUTOFF: Resolution of the structure under which the structure is accepted. E.g: 3.0 = 3.0 Aº
  - AA_CUTOFF: Number of aminoacids over which the structure is accepted. E.g: 50 (avoid small peptides)

- If the receptor has good pdb structures, then it will be processed in search of positive interacting partners. If there are,
  negative partners and positive partners will undergo the same finding good pdb structure process. Finally, positive partners
  will have a reciprocity test to check that their receptor is also a positive partner for them:

* Receptor -> Positive & Positive -> Receptor

- If a positive partner for a certain receptor is also a negative partner for the same receptor, it will be removed from both
  sets to avoid conflicting data.

- It can be run in different modes:
  -nrmd = No Remove Duplicates: only do if really trust the process of the Negatome2.0 database
  -onia = Off Negative Interface API check: Don't need to check if the negative partners have interactions with other molecules
  -outHM = Accept only 1:1 PPI, i.e. exclude interactions in which there are multiple molecules of the same protein interacting
  -t2m = Increased stringency of outHM to only 1:1 interactions but just 2 molecules in a PDB file (very restrictive)
  -cpc = Check partners coverage: check that both molecules in an interactions have a coverage bigger than COV_THRESHOLD
  -cbi = Check biological intereaction: check that the interaction between two molecules is biological by the sofware
  ProdigyCrystal

- The file prodigy_results_pdb.txt allows to skip already processed interactions by ProdigyCrystal. While create_benchmark
  is the folder where all the data is going to be stored:

- Execution: python create_interactome_db.py [-nrmd] [-onia] [-outHM] [-t2m] [-cpc] [-cbi] location/of/create_benchmark [location/of/prodigy_results_pdb.txt]

```
$ python ./scripts/create_interactome_db.py -outHM -cpc -cbi ./create_benchmark ./create_benchmark/prodigy_results_pdb.txt
```

# generate_benchmark.py

- This script creates a folder organized benchmark from the data obtained from create-db.py and
  stored in interacting_db.csv.

- It can be run with different filters:

  - outHD: Exclude homodimers
  - outIg: Exclude Immunoglobulins
  - inIg: Take only Immunoglobulins
  - outCA: Exclude molecules whose pdb structure is only Alpha Carbons
  - outMS: Exclude multiple structures for the same molecule. It takes the one with more number of partners and if tie, the one with more residues.
  - AO: Exclude homodimers, Immunoglobulins, AlphaCarbons and multiple structures.

- The script also stores information of Igs from previous runs, to make it faster to run,
  if wanted to reprocess from the beginning, use the flag -abini (ab initio).

- Execution: python generate_benchmark.py [-outHD] [-outIg] [-inIg] [-outCA] [-outMS] [-AO] [-abini] location/of/interacting_db.csv location/of/create_benchmark_folder location/of/the/benchmark [location/of/interactome_db.txt]

```
$ python generate_benchmark.py -AO - ./create_benchmark/interacting_db.csv ./create_benchmark ./BMAO-230222 ./create_benchmark/interactome_db.txt
```

# print_benchmark_data.py

""" Obtain dataset stats for the generated benchmark """
