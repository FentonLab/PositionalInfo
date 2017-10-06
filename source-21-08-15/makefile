CFLAGS =-g
CC = gcc

PROG= /usr/local/bin/pos_info_min
OBJS = main.o prob_spk_map.o prob_spk_dir.o


$(PROG) : $(OBJS)
	$(CC) $(CFLAGS)  $(OBJS) -o $(PROG) -lm; rm -f *.o


.c.o:
	$(CC) $(CFLAGS) -c $*.c

