# UNDER CONSTRUCTION

# PySNPCheck
PySNPCheck is a command line program that aims to replace the popular web
application SNPCheck. The programs main function is to find variants in
the binding regions of primers and the amplicon region. Input is taken in the
form of two primers and the chromosome they are found on. The program allows
remote access to the human genome assembly 38 and the latest build of dbSNP to achieve ths.
Furthermore users can also specify their own databases for input, allowing
for variant searching across many species. See usage for further details.

## Running the program
There are 3 different ways to run the program:
1. Default : Optimum.
2. Remotely Via NCBI servers.
3. Locally using input data sources.
These are described in more detail bellow.

The program can be run on the command line using the Python prompt
followed by pysnpcheck.py:

```
Python pysnpcheck.py [Input]
```

Usage:

pysnpcheck.py input [-h] [-rb] [-lb BLAST database path] [-sql SQlite database path] [-hg19]

Please read further for dependancies and requirements.

## Input

The program can take two types of input:

1. **A single pair of primers**: This includes typing the input directly
onto the command line. This type of input has four parameters, each must be seperated
by a single space and be in the following order: PRIMER NAME PRIMER1 PRIMER2 CHROMOSOME.
primer1 and 2 are the primer sequence, and chromosome parameter must be a single number or letter.
2. **Batch input**: This is the path to a .txt file containing multiple lines of primers in the
format described in 1. An example of this is the example_primers.txt file included in this repo.

1:
```
Python pysnpcheck.py CEP290_EX19 AGGGAGAAAGTGGGATTAAGATC AGCAAGGCAAATCAACTGGA 12
```
2:
```
Python pysnpcheck.py example_primers.txt
```




## Default run: Optimum

This will initiate a default run which has been identified as the optimum way for searchng for variants. This includes the fastest run, whilst using the least local memory space. Primers are searched against the human genome assembly
hg38, and the latest build of dbSNP. Local BLAST and the NCBIS Entrez is used here.
For this the required dependancies include:
* A local installation of command line BLAST, executables can be found
at the NCBI FTP site by following this link: [Link](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). Note this
dependency can be bypassed if using the -rb remote BLAST argument (see later).
* The Biopython Python package installed by following this link: [Link](http://biopython.org/wiki/Download)
* A local BLAST database formatted for each chromosome. This can be obtained from the following link: [Link](https://www.dropbox.com/sh/5gaaavpp0hxzaou/AADvXGBHRBj0Hxoig0DQCmbva?dl=0)
. Simply place the "blastdb" directory in the same directory as the main
Python script.
* The user should also specify their email on the first line of the email.txt file.
This is necessary when using Entrez to search dbSNP. The email.txt file must also
be placed in the same directory as the Python script.

Once these requirements have been satisfied. Input primers can be SNPChecked.

A default run looks like:

```
Python pysnpcheck.py CEP290_EX19 AGGGAGAAAGTGGGATTAAGATC AGCAAGGCAAATCAACTGGA 12
```

Results are printed to the command line.



### Hg19 liftover

The default run mode returns the results with coordinates relative to the hg38 assembly. The hg19 argument can be used to output the results with the hg19 coordinates instead.

Argument: -hg19

Optional dependancies include:
* The pyliftover Python package which can be obtained by following this link:
[Link](https://pypi.python.org/pypi/pyliftover). Neccassary if using the -hg19 argument.

Example run:
```
Python pysnpcheck.py CEP290_EX19 AGGGAGAAAGTGGGATTAAGATC AGCAAGGCAAATCAACTGGA 12 -hg19
```

## Remote run

This mode works the same way as the default mode, however remote BLAST is used in place of local BLAST. This has several advantages, such as eliminating the need for a local BLAST database to be stored. All databases are coonected to remotely via the NCBI servers. The program returns variants found in the primer regions relative to the hg38 assembly. 

For this the required dependancies include:
* Only Biopython is required.

Argument: -rb

Example of Remote run:
```
Python pysnpcheck.py CEP290_EX19 AGGGAGAAAGTGGGATTAAGATC AGCAAGGCAAATCAACTGGA 12 -rb
```

The hg19 liftover argument can be also be used with a remote run.
Example of Remote run with hg19 liftover:
```
Python pysnpcheck.py CEP290_EX19 AGGGAGAAAGTGGGATTAAGATC AGCAAGGCAAATCAACTGGA 12 -rb -hg19
```

## Local run
## Example usage
## Versioning
## Acknowledgements 
