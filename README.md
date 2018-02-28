# Primer Blast DX

PCR Primer package written in Python3

## Setup
Make a file called config.json which contains the locations of the blastn and fasta genome files.

This packages uses ncbi-blast-2.7.1 downloaded from [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
It is compiled for linux 64x. If your environment is different, specify the path to blastn in config.json.

Example:
```
{
    "blastn-path": null,
    "maize_v3":"fa/Zea_mays.AGPv3.29.dna.genome.fa"
}
```

## Install
(Optional: virtual environment)
```
virtualenv -p python3 venv
source venv/bin/activate
```

Install package
```
python setup.py install
```

Install package as development (This only creates symlinks)
```
python setup.py develop
```

# How to Run
```
import primer_blast_dx
taskResult = primer_blast_dx.run(task)
```

Take a look at test code for more detail

# Authors
Takao Shibamoto

Kokulapalan (Gokul) Wimalanathan
