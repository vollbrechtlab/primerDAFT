# Primer Blast DX

PCR Primer package written in Python3

## Setup
Make a file called config.json which contains the locations of the blastn and fasta genome files.

This packages uses ncbi-blast-2.7.1 downloaded from [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
It is compiled for linux 64x. If your environment is different, specify the path to blastn in config.

Look at `test/primer-dx.conf` for example.

## Install
(Optional: virtual environment)
```
virtualenv -p python3 venv
source venv/bin/activate
```

Install package
```
python3 setup.py install
```

Install package as development (This only creates symlinks)
```
python3 setup.py develop
```

## Useage
```
import primerDAFT
taskResult = primerDAFT.run(task)
```

Take a look at test code for more detail

## Authors
Takao Shibamoto  
Kokulapalan (Gokul) Wimalanathan
