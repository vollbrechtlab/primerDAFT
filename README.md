# Primer-Blast-DX

PCR Primer package written in Python3

## Setup
Install blastn from here
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Make a file called genome_config.json which contains the locations to fasta genome files.
Example:
```
{
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
