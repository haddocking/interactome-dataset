# parse_benchmark.py

- This script identifies different mistakes and errors in a dataset:

  - Sth wrong with the structure: two types of failure, a) Not a protein and b) Not correct structure
  - No positive/negative ligands only positives/negatives: the receptor only has a type of ligands in the benchmark
  - No interface API: there is no interface information for the molecule in the Uniprot API
  - Positive = negative: the same ligand appears as positive and negative for the same receptor

- It is recommended to be used in the mode 'eval' in which it will generate a parse_benchmark.txt file in which
  it can be seen the analysis and evaluation of every molecule in the dataset. It is a way to check that the
  benchmark creation was successful.

- Additionally, it contains the mode 'clean' in which it will create a new benchmark avoiding all the mistakes.
  Note that the preference way must be first fixing the dataset generation, but it can be used as a complement
  in some cases.

- Execution: python parse_bm_linux.py path/to/benchmark mode

```
$ python parse_bm_linux.py ./benchmarks/BM-230222 eval
```

# size_pairs_benchmark.py

- This script calculates the sizes of the the molecules in a given benchmark. It then compares the sizes of positive
  and negatives ligand to check if they are within a similitude threshold set by user. If the comparison is successful
  the pair of ligands will be copied in a new benchmark under the name that the user has decided for it.

- As output files, it generates an histogram for the distribution of sizes of all the molecules and another one
  comparing the sizes of positive ligands vs negative ones. It also creates two files stroing the sizes of the molecules:
  sizes_dic.txt which is a python dictionary that can be used in other scripts and sizes_table.txt that contains the
  same information but with a much better visualization.

- Finally, groups_table.txt is a file generated that contains the information of similar size positive-negative ligands for
  the same receptor. It will be needed later when the analysis to know the pairs of similar size.

- Similitude threshold can accept a percentage or a decimal value: 90 = 90% = 0.9

- This script does not support interface clustered benchmarks as input, so if wanted to filter a clustered benchmark, it is
  advisable to first filter by size and then proceed with the interface clustering.

- Execution: python size_pairs_benchmark.py location/of/the/benchmark name_of_new_benchmark similitude_threshold

```
$ python ./scripts/size_pairs_benchmark.py ./benchmarks/BM-230222/ NEW_BENCHMARK 90
```

# interaction_pair_ana.py

- This script obtains the ranked solutions from Haddock runs stored in /structures/it\*/file.list files,
  from an analysis pair receptor+positive_ligand - receptor+negative_ligand, combines them and sort the solutions
  according to their score and creates a combined file. If size_pairs_benchmark.py has been run
  before, the pairs are stored in the file groups_table.txt. If not, it can determine the pairs by default.

- For the scores, it sets different scenarios: Haddock score (HS), Haddock score by interface size (HSIS)
  and Haddock score by average interface size for different number of solutions (HSIAS_n)

- After that, it takes the previously generated files and calculates the Right Answer Rate (RAR) for different
  thresholds (1, 10, 50, 100, 200). The RAR is defined as the number of positive solutions divided by the number
  of total solutions under a certain threshold. E.g: under threshold 10, there are 6 positive and 4 negative
  solutions, so the RAR = 60%.

- Finally, it calculates the average RAR for each threshold for the whole benchmark as well as its standard
  deviation and plots a bar graph with that information. All the files and graphs will be stored in a new
  folder under the benchmark location with the name analysis/.

- The mode -full creates extra plots for each RAR threshold. These plots have the receptor in the x-axis sorted
  by size from smallest to biggest. On the y-axis there is the RAR for each of the analysis pairs. Each of the
  dot sizes represent the average size of both the ligands (it is useful if they have siilar size, but no point
  if they don't). Finally, the color of the dots represent the average energy (HADDOCK score) of the top threshold
  solutions and the thickness of the border, the variation of energy among solutions.

- Execution: python interaction_pair_ana.py [-full] location/of/the/benchmark [location/of/groups_table.txt] [location/of/sizes_dic.txt]

```
$ python ./scripts/interaction_pairs_ana.py [-full] ./benchmarks/OutHM-220622/BM-AO-220622 [groups_table] [sizes_dic]
```
