#===============================================================================
# User Options
#===============================================================================

COMPILER = intel
OPTIMIZE = yes
DEBUG    = yes
PROFILE  = no
STATUS   = no
PAPI     = no
MIC      = yes
VERIFY   = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = rsbench

source = \
main.c \
xs_kernel.c \
io.c \
init.c \
material.c \
utils.c

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

CFLAGS = -std=gnu99

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
  CFLAGS += -fopenmp -ffast-math
endif

# Intel compiler
ifeq ($(COMPILER),intel)
  CC = icc
  ifeq ($(MIC),yes)
    #CFLAGS += -openmp -mmic -vec-report=3 -no-opt-prefetch
    CFLAGS += -openmp -mmic -vec-report=3
  else
    CFLAGS += -openmp -vec-report=3
  endif
endif

# Linker Flags
LDFLAGS = -lm

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
  LDFLAGS += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg
  LDFLAGS  += -pg
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3
endif

# Calculation Status % Complete Indicator
ifeq ($(STATUS),yes)
  CFLAGS += -DSTATUS
endif

# PAPI Performance Counters
ifeq ($(PAPI),yes)
  CFLAGS += -DPAPI
  source += papi.c
  LDFLAGS += -lpapi
endif

# Verification of results mode
ifeq ($(VERIFY),yes)
  CFLAGS += -DVERIFICATION
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) rsbench.h
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf rsbench $(obj)

edit:
	vim -p $(source) papi.c rsbench.h

run:
	./rsbench
