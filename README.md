# Project of the course of Algorithms for Bioinformatic

Implementation of Smith-Waterman algorithm in python as project for the course of Algorithm for Bioinformatics, held by professor Enrico Blanzieri,

A Smith Waterman implementation in python, done for the course of algorithms for bioinformatics, of professor Enrico Blanzieri (A.Y 2022-2023).

### Cloning the repository
```
git clone https://github.com/iamandreatonina/Algorithms_bioinformatics.git
```
### How to use the algorithm 
```
python smith-waterman.py first_sequence second_seqeunce [-h HELP] [-g GAP PENALTY] [-m MATCH] [-p MISMATCH] [-o OUTPUT NAME]
```
By default without an output name the code will return also a file txt named SW_output

### Options/parameters : 
 * First seqeunce, which correspond of the rows of the matrix
 * Second sequence, whihc are the columns fo the matrix
 * --h, --help, which disaply the helper 
 * -g , -- gap_penalty, scoring value for the gaps, by defautl is -1
 * -m , -- match, scoring value for the matches, by defautl is 2
 * -p , -- mismatch, scoring value for the misamtches, by defautl is -1
 * -o , --output_name, Name of the output file that will be created (default: SW_output)
 

