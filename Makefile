read: main.c
	g++ -Wall -g main.c -O3 -o gtdb rocksdb-7.2.2/librocksdb.a -lpthread -lrt -ldl -lsnappy -lgflags -lz -lbz2 -llz4 -lzstd -lnuma  -ldl -lpthread -lhts
all: gtdb
