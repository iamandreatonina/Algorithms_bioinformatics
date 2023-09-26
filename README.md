# Project of the course of Algorithms for Bioinformatics

Implementation of Smith-Waterman algorithm in Python as a project for the course Algorithm for Bioinformatics, held by professor Enrico Blanzieri (A.Y 2022-2023).

[Code](https://github.com/iamandreatonina/Algorithms_bioinformatics/blob/main/Code/smith-waterman.py)

### Libraries and modules 
* [NumPy](https://numpy.org/)
* [itertools](https://docs.python.org/3/library/itertools.html)
* [argparse](https://docs.python.org/3/library/argparse.html)


### Cloning the repository
```
git clone https://github.com/iamandreatonina/Algorithms_bioinformatics.git
```
### How to use the algorithm 
```
python smith-waterman.py [-h] [-g GAP_PENALTY] [-m MATCH] [-p MISMATCH] [-o OUTPUT_NAME] seq1 seq2
```
By default without an output name, the code will return also a file text named SW_output

### Options/parameters : 
 * First sequence [required], which corresponds to the rows of the matrix
 * Second sequence [required], which are the columns of the matrix
 * -h, --help [optional] which displays the helper and exit
 * -g, -- gap_penalty [optional] scoring value for the gaps, by default is -1
 * -m, -- match [optional] scoring value for the matches, by default is 2
 * -p, -- mismatch [optional] scoring value for the mismatches, by default is -1
 * -o, --output_name [optional] Name of the output file that will be created (default: SW_output)
 
### Example 
```
python3 smith-waterman.py ATCGGCGATA ATTATACGATA -g -2 -m 3 -p -1 -o alignment_output
```
