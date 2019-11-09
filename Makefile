CC=g++
CFLAGS=-std=c++11 -O2 -Wall -Werror
LINK=-lsdsl -ldivsufsort -ldivsufsort64

quasi:  quasioptimal_partition.cpp algo_utils.o entropy_utils.o  compressor.o \
	partition.o test_lib.o vec_tools.o compressed_vector_rrr.o entropy_coding.o compressed_array.o
	$(CC) quasioptimal_partition.cpp algo_utils.o entropy_utils.o compressor.o \
	partition.o test_lib.o vec_tools.o compressed_vector_rrr.o entropy_coding.o compressed_array.o $(CFLAGS) -o quasi $(LINK)

algo_utils.o: algo_utils.cpp
	$(CC) -c algo_utils.cpp $(CFLAGS)

entropy_utils.o: entropy_utils.cpp
	$(CC) -c entropy_utils.cpp $(CFLAGS)

partition.o: partition.cpp vec_tools.o entropy_utils.o
	$(CC) -c partition.cpp $(CFLAGS)

test_lib.o: test_lib.cpp algo_utils.o entropy_utils.o partition.o vec_tools.o compressed_array.o entropy_utils.o compressor.o
	$(CC) -c test_lib.cpp $(CFLAGS)

vec_tools.o: vec_tools.cpp
	$(CC) -c vec_tools.cpp $(CFLAGS)

compressed_vector_rrr.o: compressed_vector_rrr.cpp
	$(CC) -c compressed_vector_rrr.cpp $(CFLAGS)

entropy_coding.o: entropy_coding.cpp
	$(CC) -c entropy_coding.cpp $(CFLAGS)

compressed_array.o: compressed_array.cpp partition.o entropy_coding.o compressor.o compressed_vector_rrr.o
	$(CC) -c compressed_array.cpp $(CFLAGS)

compressor.o: compressor.cpp partition.o entropy_coding.o
	$(CC) -c compressor.cpp $(CFLAGS)


clean:
	rm *.o quasi
