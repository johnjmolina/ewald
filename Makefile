ARCH     = macosx
CC	 = clang
CXX      = clang++
LINKS    = -I/usr/local/include -L/usr/local/lib -lm -stdlib=libc++
CCOPT    = -O0 -g -stdlib=libc++ 

AUX_OBJS = alloc.o\
       rigid_body.o\
       ewald.o\
       dipole_gold.o\
       gen_shell.o
OBJS1 = test_one.o $(AUX_OBJS)
OBJS2 = test_two.o $(AUX_OBJS)
OBJS3 = test_two_random.o $(AUX_OBJS)
OBJS4 = test_hundred_random.o $(AUX_OBJS)

TEST1 = test_one
TEST2 = test_two
TEST3 = test_two_random
TEST4 = test_hundred_random
## Implicit rules

.SUFFIXES: .c .cxx .o .out

## Build rules
all: $(TEST1) $(TEST2) $(TEST3) $(TEST4)

$(TEST1): $(OBJS1)
	$(CXX) $(OBJS1) -o $(TEST1).x $(CFLAGS) $(LINKS)

$(TEST2): $(OBJS2)
	$(CXX) $(OBJS2) -o $(TEST2).x $(CFLAGS) $(LINKS)

$(TEST3): $(OBJS3)
	$(CXX) $(OBJS3) -o $(TEST3).x $(CFLAGS) $(LINKS)

$(TEST4): $(OBJS4)
	$(CXX) $(OBJS4) -o $(TEST4).x $(CFLAGS) $(LINKS)

## Compile

.cxx.o: 
	$(CXX) -c $< $(CFLAGS) -o $@

.c.o: 
	$(CC) -c $< $(CFLAGS) -o $@

## Clean
clean:
	rm -f *.o *.x *~

cleanall:
	rm -f *.o *.x *~
