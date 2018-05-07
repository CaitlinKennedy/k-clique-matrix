CC=gcc-7
CFLAGS=-O9
CP=g++-7

all: kClistMatrix test

kClistMatrix: kClistMatrix.cpp
	$(CP) $(CFLAGS) kClistMatrix.cpp -o kClistMatrix -fopenmp

test1:
	@./kClistMatrix 1 3 simplegraph.txt 1 > test1.txt
	@diff -q test1.txt simplegraphtest3.txt
	@rm test1.txt

test2:
	@./kClistMatrix 1 4 simplegraph.txt 1 > test2.txt
	@diff -q test2.txt simplegraphtest4.txt
	@rm test2.txt

test3:
	@./kClistMatrix 1 5 simplegraph.txt 1 > test3.txt
	@diff -q test3.txt simplegraphtest5.txt
	@rm test3.txt

test4:
	@./kClistMatrix 1 3 edgelist.txt 1 > test4.txt
	@diff -q test4.txt edgelisttest3.txt
	@rm test4.txt

test5:
	@./kClistMatrix 1 4 edgelist.txt 1 > test5.txt
	@diff -q test5.txt edgelisttest4.txt
	@rm test5.txt

test6:
	@./kClistMatrix 1 5 edgelist.txt 1 > test6.txt
	@diff -q test6.txt edgelisttest5.txt
	@rm test6.txt

test7:
	@./kClistMatrix 3 3 edgelist.txt 1 > test7.txt
	@diff -q test7.txt edgelisttest3.txt
	@rm test7.txt

test: test1 test2 test3 test4 test5 test6 test7

clean:
	rm  kClistMatrix
