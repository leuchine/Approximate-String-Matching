#Makefile for similarity search in complex structure databases
CC=g++ -Wno-deprecated -DNDEBUG -O3
CFLAGS=-c
#Change BIN_HOME to your local binary folder
BIN_HOME=/home/xiaoli/tkde2015/
OFILES=time.o gram.o list.o filtb.o query.o db.o

all:sindex psearch

rebuild:clean all
	
clean:
	rm -f *.o
	rm -f sindex psearch

sindex:${OFILES} sindex.o
	${CC} -o ${BIN_HOME}$@ ${OFILES} sindex.o

psearch:${OFILES} psearch.o
	${CC} -lpthread -o ${BIN_HOME}$@ ${OFILES} psearch.o

# ------------------------------------------
sindex.o:index.cpp
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

psearch.o:pipeSearch.cpp
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
	

# COMMON OBJECT FILES

time.o:Time.cpp Time.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
	
gram.o:Gram.cpp Gram.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

list.o:GramList.cpp GramList.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
	
filtb.o:CountFilter.cpp CountFilter.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

query.o:Query.cpp Query.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

db.o:SeqDB.cpp SeqDB.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
