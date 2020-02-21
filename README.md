# score-assemblies

A snakemake-wrapper for evaluating *de novo* bacterial genome assemblies, e.g. from Oxford Nanopore (ONT) or Illumina sequencing.

The workflow includes the following programs:
* [pomoxis](https://github.com/nanoporetech/pomoxis) assess_assembly and assess_homopolymers
* dnadiff from the [mummer](https://mummer4.github.io/index.html) package
* [QUAST](http://quast.sourceforge.net/quast)
* [BUSCO](https://busco.ezlab.org/)
* [ideel](https://github.com/mw55309/ideel/)

## Installation
Clone repository, for example:
```
git clone https://github.com/pmenzel/score-assemblies.git /opt/software/score-assemblies
```
Install dependencies into an isolated conda environment
```
conda env create -n score-assemblies --file /opt/software/score-assemblies/environment.yaml
```
and activate environment:
```
source activate score-assemblies
```

## Usage
First, prepare a data folder, which must contain subfolders for the assemblies
and the reference genomes.
```
.
├── assemblies
│   ├── raven.fa
│   ├── raven+medaka.fa
│   └── raven+medaka+pilon.fa
│
└── references
    ├── reference1.fa
    └── reference2.fa
```

Run workflow, e.g. with 20 threads:
```
snakemake -k -s /opt/software/score-assemblies/Snakefile --cores 20
```

### Modules
score-assemblies will run these programs on each assembly:

#### assess_assembly and assess_homopolymers

Each assembly will be scored against each reference genome using the
`assess_assembly` and `assess_homopolymers` scripts from
[pomoxis](https://github.com/nanoporetech/pomoxis).  Additionally to the tables
and plots from pomoxis, summary plots for each reference genome will be plotted
in the files `<reference>_assess_assembly_all_meanQ.pdf` and
`<reference>_assess_homopolymers_all_correct_len.pdf`.

#### BUSCO

Set the lineage via the snakemake call:
```
snakemake -k -s /opt/software/score-assemblies/Snakefile --cores 20 --config busco_lineage=burkholderiales
```
If not set, the default lineage `bacteria` will be used.
Available datasets can be listed with `busco --list-datasets`

The number of complete, fragmented and missing BUSCOs per assembly is plotted in the file `busco/busco_stats.pdf`.


#### dnadiff
Each assembly is compared with each reference and the output files will be
located in `dnadiff/<reference>/<assembly>-dnadiff.report`.  The values for
`AvgIdentity` and `TotalIndels` are extracted from these files and are plotted
for each reference in the files `dnadiff/<reference>_dnadiff_stats.pdf`.

#### QUAST

One QUAST report is generated for each reference genome, containing the results for all assemblies.
The report files are located in `quast/<reference>/report.html`.

#### ideel

Open reading frames are predicted from each assembly via Prodigal and are
search in the Uniprot sprot database with diamond, retaining the best alignment
for each ORF. For each assembly, the distribution of the ratios between length
of the ORF and the matching database sequence are plotted to the file
`ideel/ideel_stats.pdf`.

