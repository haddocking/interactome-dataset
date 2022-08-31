# cluster_interfaces_benchmark.py

""" This is a script to identify possible interfaces in the receptor molecules of a benchmark and cluster
them according to the relative distances from the observed interfaces with different ligands. """

- It gets the information from https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues and can accept different filters:

  - Exclusion of Immunoglobulins: -outIg
  - Exclusion of Homodimers: -outHD
  - Only 1:1 PPI: -outHM (HomoMultimers)
  - Only 2 protein molecules in a PDB: -t2m
  - Check partners coverage: -cpc. Both protein in a PPI, need to have a coverage higher than the established cutoff
  - Check biological interaction: -cbi. Performed by ProdigyCrystal
  - All Out: -AO. Exclude Ig, HD and HM.

- It gets the information of the PDBs in a benchmark to check if the interfaces are present in the strcuture. For that it
  is crucial that those PDBs are correctly renumbered. For doing that, it is strongly recommended to use the Open Software
  PDBrenum.

- For doing the clustering, it first calculates the residue-residue distances in the PDB structure from the euclidean distance
  of their Alpha Carbons. To calculate the D-values for the clustering it does 4 different approaches, of those, the preferred one
  is the sine one.

- The optional files prodigy_results_pdb.txt and interactome_db.txt stores information on biological tests and lists of residues
  for the receptors in a benchmark, respectively.

- Execution: python cluster_interfaces.py location/of/the/benchmark [prodigy_results_pdb] [interactome_db] [-outIg] [-outHD] [-outHM] [-AO] [-cpc] [-cbi]

```
$ python ./scripts/cluster_interfaces.py -AO -cbi ./benchmarks/BM-230222/ ./create_benchmark/prodigy_results_pdb.txt
```

# cluster_distance_matrix.py

""" Provided an interface matrix in input, this script clusters the interfaces and saves the dendrograms as png files """

- Usage:

```
python3 cluster_distance_matrix.py ${PATH_TO_DISTANCE_FILE}
```

# write_restraints.py

""" This scripts creates a clustered benchmark from the clustered interfaces obtained after running cluster_interfaces.py and
cluster_distance_matrix.py. It gets the interfaces from an interface_dictionary.txt file and the different clustered
ligands from the output file of cluster_distance_matrix.py """

- It creates a restraints .tbl file combining the active residues of the interface at the receptors and the passive residues
  from the execution of calc_accessibility.py at the ligands PDBs.

- Execution: python write_restraints.py location/of/the/benchmark location/of/interface_dictionary location/of/clusters_file name_of_new_benchmark

```
$ python ../scripts/write_restraints.py ./OutHM-220622/BM-AO-220622 ./OutHM-220622/BM-AO-220622/interface_analysis/interface_dictionary.txt
./OutHM-220622/BM-AO-220622/cluster_interfaces/clusters_single_thr-0.7071.out NESTOR
```
