# primerDAFT

Primer Design and Filtering Tool

## Setup
Make a config which contains the locations of the blastn and fasta genome files.

This packages uses ncbi-blast-2.7.1 downloaded from [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
It is compiled for linux 64x. If your environment is different, install the appropriate one and specify the path to blastn in config.

Look at `test/pdaft.conf` for example.

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
taskResult = primerDAFT.run(task, configFile)
```

Take a look at `test/test.py` for more detail

## Authors
Takao Shibamoto  
Kokulapalan (Gokul) Wimalanathan
