#
# Makefile for the CMF code
# Hugues Talbot	18 Jun 2008
#

CC=gcc
CFLAGS=-Wall -g -O2

TARGET=cont_max_flow

SOURCE= \
	main.c \
	bimage.c \
	bimage_utils.c \
	contmaxflow.c \
	FMM.c \
	LSTBmsg.c \
	LSTBmsgs.c \
	maxflow.c \
	lbimage_utils.c \
	lcontmaxflow.c \
	lreadLSTBmsgs.c \
	sepgauss.c

OBJECTS=${SOURCE:.c=.o}

.c.o:
	${CC} ${CFLAGS} -c $<

${TARGET}: ${OBJECTS}
	${CC} ${CFLAGS} -o ${TARGET} ${OBJECTS} -lpthread -lm 

all: clean depend 
	${MAKE} ${TARGET}

clean:
	-rm ${TARGET} ${OBJECTS} *~

realclean: clean
	-rm g.pgm result.pgm

test: ${TARGET}
	time ./${TARGET} im48_48.pgm source48_48.pgm sink48_48.pgm result.pgm 48 48 -d

help:
	@echo "Make:"
	@echo "  <default>    : recompiles the executable"
	@echo "  clean        : remove all compiled files"
	@echo "  realclean    : also removes produced images"
	@echo "  depend       : computes dependencies"
	@echo "  test         : runs a test"

depend: ${SOURCE}
	${CC} ${CFLAGS} -M $^ > makedepend

makedepend:
	touch makedepend
	${MAKE} depend

include makedepend