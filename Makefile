SHELL := /bin/bash
CC=gcc-7
CFLAGS=-O9
CP=g++-7

all: kClistMatrix kClistMultiple test

kClistMatrix: kClistMatrix.cpp
	$(CP) $(CFLAGS) kClistMatrix.cpp -o kClistMatrix -fopenmp

kClistMultiple: kClistMatrixMultiple.cpp
	$(CP) $(CFLAGS) kClistMatrixMultiple.cpp -o kClistMultiple -fopenmp

test1:
	@./kClistMatrix 1 3 tests/simplegraph.txt 1 > test.txt
	@diff -q <(sort test.txt) tests/simplegraphtest3.txt
	@rm test.txt

test2:
	@./kClistMatrix 1 4 tests/simplegraph.txt 1 > test.txt
	@diff -q <(sort test.txt) tests/simplegraphtest4.txt
	@rm test.txt

test3:
	@./kClistMatrix 1 5 tests/simplegraph.txt 1 > test.txt
	@diff -q <(sort test.txt) tests/simplegraphtest5.txt
	@rm test.txt

test4:
	@./kClistMatrix 1 3 tests/edgelist.txt 1 > test.txt
	@diff -q <(sort test.txt) tests/edgelisttest3.txt
	@rm test.txt

test5:
	@./kClistMatrix 1 4 tests/edgelist.txt 1 > test.txt
	@diff -q <(sort test.txt) tests/edgelisttest4.txt
	@rm test.txt

test6:
	@./kClistMatrix 1 5 tests/edgelist.txt 1 > test.txt
	@diff -q <(sort test.txt) tests/edgelisttest5.txt
	@rm test.txt

test7:
	@./kClistMatrix 2 3 tests/edgelist.txt 1 > test.txt
	@diff -q <(sort test.txt) tests/edgelisttest3.txt
	@rm test.txt

test8:
	@./kClistMatrix 2 2 tests/edgelist.txt 1 > test.txt
	@diff -q <(sort test.txt) <(sort tests/edgelisttest8.txt)
	@rm test.txt

test9:
	@./kClistMatrix 2 2 tests/simplegraph.txt 1 > test.txt
	@diff -q <(sort test.txt) <(sort tests/simplegraphtest9.txt)
	@rm test.txt

test10:
	@./kClistMatrix 2 5+ tests/simplegraph.txt 1 > test.txt
	@diff -q <(sort test.txt) tests/simplegraphtest10.txt
	@rm test.txt

test: test1 test2 test3 test4 test5 test6 test8 test9 test10

clean:
	rm  kClistMatrix kClistMultiple
