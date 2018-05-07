
# Info:
Code to generate the clique weighted matrix based on github repository maxdan94/kClist

# kClistMatrix.cpp

This program iterates over all k-cliques. And generates the k-clique weighted adjacency matrix.  It is highly scallable.

## To compile:
"g++ kClistMatrix.cpp -O9 -o kClistMatrix -fopenmp".

## To execute:
"./kClistDens p k edgelist.txt".  
- "edgelist.txt" should contain the graph: one edge on each line (i, j, w)  separated by a space.
- k is the size of the k-cliques
- p is the number of threads

