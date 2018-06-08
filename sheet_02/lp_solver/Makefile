# # This file is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.

# # It is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.

# # You should have received a copy of the GNU General Public License
# # along with this file.  If not, see <http://www.gnu.org/licenses/>.

# Please set the environment path and uncomment the lines correspondingly.
# You might also need to adopt the  LPINCLUDE & LPLIB  paths further below.
#GUPATH=$(GUROBI_HOME)
#CPLEXPATH=$(CPLEX_HOME)
QSPATH="."

MACH := ${shell uname -srvm}
OS    = ${word 1, $(MACH)}
TYPE  = ${shell uname -m}

ifneq ($(QSPATH),)
LPINCLUDE=$(QSPATH)
LPLIB=$(QSPATH)/qsopt.a.cygwin32
ifeq ($(TYPE),x86_64)
LPLIB=$(QSPATH)/qsopt.a.cygwin64
endif
ifeq ($(OS),Linux)
LPLIB=$(QSPATH)/qsopt.a.linux32
ifeq ($(TYPE),x86_64)
LPLIB=$(QSPATH)/qsopt.a.linux64
endif
endif
LPSOURCE=lpqsopt.o
endif


ifneq ($(GUPATH),)
LPINCLUDE=$(GUPATH)/include
LPLIB=$(GUPATH)/lib/libgurobi.so.4.0.0
LPSOURCE=lpgurobi.o
GRBMWIS=mwis_grb.o
GUROBI_FLAG=-DUSE_GUROBI
endif

ifneq ($(CPLEXPATH),)
PROCESSOR := $(shell uname -p)
LPINCLUDE=$(CPLEXPATH)/include/ilcplex/
LPLIB=$(CPLEXPATH)/lib/x86-64_sles10_4.1/static_pic/libcplex.a
LPSOURCE=lpcplex.o
GRBMWIS=
GUROBI_FLAG=
endif


export CC=gcc
export LD=gcc
CFLAGS+= -g
#CFLAGS+= -O3

export CXX=g++
export LDXX=g++



# For static code analysis with clang we use the clang compiler.
# NOTE: You need to disable optimization otherwise
#export CC=clang

#
# Valgrind does not support fegetround & fesetround. With following compile option
# their use is circumvented. We also recommend to use QSopt as the LP-solver while
# debugging with valgrind, as the commercial solvers impose valgrind errors internally.
#
#CFLAGS+= -DCOMPILE_FOR_VALGRIND



####################################################
# Below this comment changes should be unnecessary.#
####################################################

CFLAGS += -std=c99 -pedantic -Wall -Wextra -Wpedantic -Wextra -Wshadow -W -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wpointer-arith -Wnested-externs -Wundef -Wcast-qual -Wcast-align -Wwrite-strings -I$(LPINCLUDE) -no-pie
export CFLAGS

CXXFLAGS +=  -pedantic -Wall -Wextra -Wpedantic -Wextra -Wshadow -W  -Wmissing-declarations -Wpointer-arith  -Wundef -Wcast-qual -Wcast-align -Wwrite-strings -I$(LPINCLUDE) -no-pie
export CXXFLAGS

MSSFILES=mss.o $(LPSOURCE)
MSSCXXFILES=msscxx.o $(LPSOURCE)

all: mss


mss: $(MSSFILES)
	$(CC) $(CFLAGS) -o mss $(MSSFILES)  $(LPLIB) -lm -lpthread

msscxx: $(MSSCXXFILES)
	$(CXX) $(CXXFLAGS) -o msscxx $(MSSCXXFILES)  $(LPLIB) -lm -lpthread

clean:
	rm -f *.o mss msscxx

lpqsopt.o: lpqsopt.c  lp.h
lpcplex.o: lpcplex.c  lp.h

mss.o:     mss.c lp.h
dummy.o: dummy.c

msscxx.o:   msscxx.C lp.h
dummycxx.o: dummycxx.C
