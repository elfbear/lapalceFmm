C = icc -g
FORTRAN = ifort
CFLAGS = $(PROFILE) #-O3
FFLAGS = $(PROFILE) # -O3
LFLAGS = $(PROFILE)

COMMON = ./common
OBJDIR = ./build
LAPDIR = ./laplace
YUKDIR = ./yukawa

COMC = $(wildcard $(COMMON)/*.c)
LAPC = $(wildcard $(LAPDIR)/*.c)
YUKC = $(wildcard $(YUKDIR)/*.c)
YUKF = $(wildcard $(YUKDIR)/*.f)

COMOBJ = $(patsubst $(COMMON)/%.c, $(OBJDIR)/%.o, $(COMC))
LAPCOBJ = $(patsubst $(LAPDIR)/%.c, $(OBJDIR)/%.o, $(LAPC))
YUKCOBJ = $(patsubst $(YUKDIR)/%.c, $(OBJDIR)/%.o, $(YUKC))
YUKFOBJ = $(patsubst $(YUKDIR)/%.f, $(OBJDIR)/%.o, $(YUKF))

LIB_LAP = ./build/libadap_laplace.a
LIB_YUK = ./build/libadap_yukawa.a

yuk: $(COMOBJ) $(YUKCOBJ) $(YUKFOBJ)
	$(FORTRAN) $(LFLAGS) $(COMOBJ) $(YUKCOBJ) $(YUKFOBJ)

lib_lap: $(COMOBJ) $(LAPCOBJ)
	ar -cvq $(LIB_LAP) $(COMOBJ) $(LAPCOBJ)

lib_yuk: $(COMOBJ) $(YUKCOBJ) $(YUKFOBJ)
	ar -cvq $(LIB_YUK) $(COMOBJ) $(YUKCOBJ) $(YUKFOBJ)

$(OBJDIR)/%.o: $(COMMON)/%.c
	$(C) $(CFLAGS) -c -I./include -o $@ $<

$(OBJDIR)/%.o: $(LAPDIR)/%.c
	$(C) $(CFLAGS) -c -I./include -o $@ $<

$(OBJDIR)/%.o: $(YUKDIR)/%.c
	$(C) $(CFLAGS) -c -I./include -o $@ $<

$(OBJDIR)/%.o: $(YUKDIR)/%.f
	$(FORTRAN) $(FFLAGS) -c -Iinc -o $@ $<

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(LIB_LAP) $(LIB_YUK)
	rm -f *~ $(LAPDIR)/*~ $(YUKDIR)/*~ $(COMMON)/*~
	rm -f ./include/*~


