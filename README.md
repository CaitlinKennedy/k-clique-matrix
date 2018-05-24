
# Info:
Code to generate the clique weighted matrix based on github repository maxdan94/kClist
(kClistMatrix.cpp is a heavily modified version of maxdan94/kClist/kClistDens.c)

# kClistMatrix.cpp

This program iterates over all k-cliques and generates the k-clique weighted adjacency matrix.  

## To compile:
"g++ kClistMatrix.cpp -O9 -o kClistMatrix -fopenmp".
Running "make" will compile the code and run some basic tests

## To execute:
"./kClistMatrix p k edgelist.txt".  
- "edgelist.txt" should contain the graph: one edge on each line (i, j, w)  separated by a space.
- k is the size of the k-clique. Add a plus at the end of the number to sum all cliques up to k. Ex "4+" would return the sum of weightings for edge, triangle, and 4 clique"
- p is the number of threads

