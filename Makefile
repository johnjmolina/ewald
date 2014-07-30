ARCH     = macosx
CC	 = clang
CXX      = clang++
#LINKS    = -I/usr/local/include -L/usr/local/lib -lm -stdlib=libc++
#CCOPT    = -O0 -g -stdlib=libc++ 
LINKS    = -I/usr/local/include -L/usr/local/lib -lm -stdlib=libstdc++
CCOPT    = -O0 -g -stdlib=libstdc++ 
CFLAGS   = $(CCOPT)

AUX_OBJS = alloc.o\
       ewald.o\
       ewald_gold.o\
       gen_shell.o

OBJS0 = test_charge.o $(AUX_OBJS)
OBJS1 = test_one.o $(AUX_OBJS)
OBJS2 = test_two.o $(AUX_OBJS)
OBJS3 = test_two_random.o $(AUX_OBJS)
OBJS4 = test_hundred_random.o $(AUX_OBJS)
OBJS5 = test_charge_hundred_random.o $(AUX_OBJS)
OBJS6 = test_fcc.o $(AUX_OBJS)
OBJS7 = test_hundred_random_charge_dipole.o $(AUX_OBJS)
OBJS8 = test_quadrupole.o $(AUX_OBJS)
OBJS9 = test_2quadrupole.o $(AUX_OBJS)

TEST0 = test_zero
TEST1 = test_one
TEST2 = test_two
TEST3 = test_two_random
TEST4 = test_hundred_random
TEST5 = test_charge_hundred_random
TEST6 = test_fcc
TEST7 = test_hundred_q+mu
TEST8 = test_quadrupole
TEST9 = test_2quadrupole

## Implicit rules

.SUFFIXES: .c .cxx .o .out

## Build rules
all: $(TEST0) $(TEST1) $(TEST2) $(TEST3) $(TEST4) $(TEST5) $(TEST6) $(TEST7) $(TEST8) $(TEST9)

$(TEST0): $(OBJS0)
	$(CXX) $(OBJS0) -o $(TEST0).x $(CFLAGS) $(LINKS)

$(TEST1): $(OBJS1)
	$(CXX) $(OBJS1) -o $(TEST1).x $(CFLAGS) $(LINKS)

$(TEST2): $(OBJS2)
	$(CXX) $(OBJS2) -o $(TEST2).x $(CFLAGS) $(LINKS)

$(TEST3): $(OBJS3)
	$(CXX) $(OBJS3) -o $(TEST3).x $(CFLAGS) $(LINKS)

$(TEST4): $(OBJS4)
	$(CXX) $(OBJS4) -o $(TEST4).x $(CFLAGS) $(LINKS)

$(TEST5): $(OBJS5)
	$(CXX) $(OBJS5) -o $(TEST5).x $(CFLAGS) $(LINKS)

$(TEST6): $(OBJS6)
	$(CXX) $(OBJS6) -o $(TEST6).x $(CFLAGS) $(LINKS)

$(TEST7): $(OBJS7)
	$(CXX) $(OBJS7) -o $(TEST7).x $(CFLAGS) $(LINKS)

$(TEST8): $(OBJS8)
	$(CXX) $(OBJS8) -o $(TEST8).x $(CFLAGS) $(LINKS)

$(TEST9): $(OBJS9)
	$(CXX) $(OBJS9) -o $(TEST9).x $(CFLAGS) $(LINKS)

## Compile

.cxx.o: 
	$(CXX) -c $< $(CFLAGS) -o $@

.c.o: 
	$(CC) -c $< $(CFLAGS) -o $@

## Clean
clean:
	rm -f *.o *.x *.dat *~

cleanall:
	rm -f *.o *.x *~
