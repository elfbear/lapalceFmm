C = icc -g
FORTRAN = ifort

CFLAGS = $(PROFILE) # -O3 -g
LFLAGS = $(PROFILE)

LAPFLAGS = -lcilkrts -L../lib -ladap_laplace
YUKFLAGS = -nofor_main -lcilkrts -L../lib -ladap_yukawa
EXEC = test_adap_fmm

lap: 
	$(C) $(CFLAGS) -c -I../include $(EXEC).c
	$(C) $(LFLAGS) $(EXEC).o $(LAPFLAGS) -o $(EXEC)

laptest3:
	./$(EXEC) -d 1 -a 3 -b 0
#	./$(EXEC) -d 2 -a 3 -b 0
#	./$(EXEC) -d 3 -a 3 -b 0

laptest6:
	./$(EXEC) -d 1 -a 6 -b 0 
#	./$(EXEC) -d 2 -a 6 -b 0
#	./$(EXEC) -d 3 -a 6 -b 0


yuk:
	$(C) $(CFLAGS) -c -I../include $(EXEC).c
	$(FORTRAN) $(LFLAGS) $(EXEC).o $(YUKFLAGS) -o $(EXEC)

yuktest3:
	./$(EXEC) -d 1 -a 3 -b 0.1
	./$(EXEC) -d 2 -a 3 -b 0.1
	./$(EXEC) -d 3 -a 3 -b 0.1

yuktest6:
	./$(EXEC) -d 1 -a 6 -b 0.1
	./$(EXEC) -d 2 -a 6 -b 0.1
	./$(EXEC) -d 3 -a 6 -b 0.1

clean:
	rm -f *.o
	rm -f *.txt
	rm -f $(EXEC)
	rm -f *~


