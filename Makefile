SHELL := /bin/bash
CC=gcc-7
CFLAGS=-O9
CP=g++-7

all: kClistMatrix test

kClistMatrix: kClistMatrix.cpp
	$(CP) $(CFLAGS) kClistMatrix.cpp -o kClistMatrix -fopenmp

test1:
	@./kClistMatrix 1 3 tests/simplegraph.txt 1 > test1.txt
	@diff -q <(sort test1.txt) tests/simplegraphtest3.txt
	@rm test1.txt

test2:
	@./kClistMatrix 1 4 tests/simplegraph.txt 1 > test2.txt
	@diff -q <(sort test2.txt) tests/simplegraphtest4.txt
	@rm test2.txt

test3:
	@./kClistMatrix 1 5 tests/simplegraph.txt 1 > test3.txt
	@diff -q <(sort test3.txt) tests/simplegraphtest5.txt
	@rm test3.txt

test4:
	@./kClistMatrix 1 3 tests/edgelist.txt 1 > test4.txt
	@diff -q <(sort test4.txt) tests/edgelisttest3.txt
	@rm test4.txt

test5:
	@./kClistMatrix 1 4 tests/edgelist.txt 1 > test5.txt
	@diff -q <(sort test5.txt) tests/edgelisttest4.txt
	@rm test5.txt

test6:
	@./kClistMatrix 1 5 tests/edgelist.txt 1 > test6.txt
	@diff -q <(sort test6.txt) tests/edgelisttest5.txt
	@rm test6.txt

test7:
	@./kClistMatrix 2 3 tests/edgelist.txt 1 > test7.txt
	@diff -q <(sort test7.txt) tests/edgelisttest3.txt
	@rm test7.txt

test: test1 test2 test3 test4 test5 test6 test7

clean:
	rm  kClistMatrix
