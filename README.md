[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15111518.svg)](https://doi.org/10.5281/zenodo.15111518)

vcf2diem.py
```
Usage: 
 vcf2diem.py -v <FILE> [-h -n -l <INT> -m -c <INT>] [-f <STR>]... [-s <STR>]...

Options:
 -v, --vcf <FILE>                       VCF file
 -n, --no_annotations                   Suppress writing SNP annotations
 -l, --non_callable_limit <INT>         Maximum number of noncallable genotypes allowed per site (default = no limit)
 -f, --exclude-chromosomes <STR>         Chromosome to exclude (default = None)
 -s, --exclude-samples <STR>             Sample to exclude (default = None)
 -m, --missing-homs                     Include sites missing a ref and/or alt homozygote (default = False)
 -c, --chunks <INT>                     Split diem_files/diem_input/per_chromosome into chunked files (uses /bin/bash)
 -h, --help                             Print this message
"""

Example Command
mamba install -c conda-forge -c bioconda numpy pandas scikit-allel docopt
./vcf2diem.py -v vcf_file.vcf.gz -f contig_1 -s sample_1 -s sample_2
```
