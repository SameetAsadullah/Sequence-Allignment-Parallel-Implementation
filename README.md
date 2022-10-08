<h1 align="center">Sequence Allignment Parallel Implementation</h1>

### Description
A `Parallel Implementation` of strings `Sequence Allignment` in `C++` using MPI (Especially used in `Bio-Informatics`). It reads input from `input.txt` and writes the output to file `output.txt`.

### Format of Input
- The first line of the file contains W1 W2 W3 W4 for computation of Alignment Score
- The second line of the file contains a sequence Seq1 (not more than 10000 letters)
- The third line of the file contains a sequence Seq2 (not more than 5000 letters)
- The forth line contains word maximum or minimum which defines the goal of the search

### Format of Output
- The first line of the file contains a Mutant of the Seq2 that produced an answer to the problem
- The second line of the file contains its Offset and the Alignment Score found
