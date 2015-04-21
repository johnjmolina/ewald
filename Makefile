ifeq ($(ENV), CLANG)
     CC       = clang
     CXX      = clang++
     COPT     = -O3 -stdlib=libc++
     LINKS    = -I/usr/local/include -L/usr/local/lib -lm -stdlib=libc++
endif
ifeq ($(ENV), CLANG_STD)
     CC       = clang
     CXX      = clang++
     COPT     = -O3 -stdlib=libstdc++
     LINKS    = -I/usr/local/include -L/usr/local/lib -lm -stdlib=libstdc++
endif
ifeq ($(ENV), ICC)
     CC	      = icc
     CXX      = icpc
     CCOPT    = -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 -w0 -lstdc++
     LINKS    = -lm -lstdc++
endif
ifeq ($(ENV), ICC_OMP)
     MKL_DIR  = /home/opt/intel/composer_x2_2013/mkl
     MKL_PATH = $(MKL_DIR)/lib/intel64
     MKL_INCLUDE_PATH = $(MKL_DIR)/include
     CC	      = icc
     CXX      = icpc
     CCOPT    = -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 \
	-ip -openmp -parallel -w0 -L$(MKL_PATH) -I$(MKL_INCLUDE_PATH)
     LINKS    = -lstdc++ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lm
endif

CFLAGS   = $(CCOPT)

AUX_OBJS = alloc.o\
       ewald.o\
       ewald_gold.o\
       gen_shell.o

# charges
OBJSA = test_q_2.o $(AUX_OBJS)
TESTA = test_q_2

OBJSB = test_q_100_rand.o $(AUX_OBJS)
TESTB = test_q_100_rand

OBJSC = test_q_256_fcc.o $(AUX_OBJS)
TESTC = test_q_256_fcc

# dipoles
OBJSD = test_mu_1.o $(AUX_OBJS)
TESTD = test_mu_1

OBJSE = test_mu_2.o $(AUX_OBJS)
TESTE = test_mu_2

OBJSF = test_mu_2_rand.o $(AUX_OBJS)
TESTF = test_mu_2_rand

OBJSG = test_mu_100_rand.o $(AUX_OBJS)
TESTG = test_mu_100_rand

# mixtures
OBJSH = test_q+mu_100_rand.o $(AUX_OBJS)
TESTH = test_q+mu_100_rand

OBJSI = test_q+mu+alpha_2.o $(AUX_OBJS)
TESTI = test_q+mu+alpha_2

## Implicit rules

.SUFFIXES: .c .cxx .o .out

## Build rules
all: $(TESTA) $(TESTB) $(TESTC) $(TESTD) $(TESTE) $(TESTF) $(TESTG) $(TESTH) $(TESTI) $(TESTJ) $(TESTK)

$(TESTA): $(OBJSA)
	$(CXX) $(OBJSA) -o $(TESTA).x $(CFLAGS) $(LINKS)

$(TESTB): $(OBJSB)
	$(CXX) $(OBJSB) -o $(TESTB).x $(CFLAGS) $(LINKS)

$(TESTC): $(OBJSC)
	$(CXX) $(OBJSC) -o $(TESTC).x $(CFLAGS) $(LINKS)

$(TESTD): $(OBJSD)
	$(CXX) $(OBJSD) -o $(TESTD).x $(CFLAGS) $(LINKS)

$(TESTE): $(OBJSE)
	$(CXX) $(OBJSE) -o $(TESTE).x $(CFLAGS) $(LINKS)

$(TESTF): $(OBJSF)
	$(CXX) $(OBJSF) -o $(TESTF).x $(CFLAGS) $(LINKS)

$(TESTG): $(OBJSG)
	$(CXX) $(OBJSG) -o $(TESTG).x $(CFLAGS) $(LINKS)

$(TESTH): $(OBJSH)
	$(CXX) $(OBJSH) -o $(TESTH).x $(CFLAGS) $(LINKS)

$(TESTI): $(OBJSI)
	$(CXX) $(OBJSI) -o $(TESTI).x $(CFLAGS) $(LINKS)

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
