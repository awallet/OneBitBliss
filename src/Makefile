# Copyright or © or Copr. Leo Ducas and Tancrede Lepoint.
#
# Leo.Ducas@ens.fr and Tancrede.Lepoint@ens.fr
#
# This software is a computer program whose purpose is to provide to the
# research community a proof-of-concept implementation of the BLISS
# digital signature scheme of Ducas, Durmus, Lepoint and Lyubashevsky
# appeared at Crypto 2013.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

CXX	= g++ -no-pie
IFLAGS  =
LDFLAGS = -lgmpxx -lgmp -lntl -lmpfr -lgmp -lcrypto -lboost_program_options
CXXFLAGS= -march=native -mtune=native -O3 -std=gnu++11 -Wall -I/usr/include/eigen3 -g

SRCS    = Entropy.cpp Fft.cpp KeyGen.cpp Sampler.cpp Setup.cpp Sign.cpp Verify.cpp
OBJS    = $(subst .cpp,.o,$(SRCS))

demo: demo.o $(OBJS)
	$(CXX) -o demo demo.o $(OBJS) $(LDFLAGS)

bliss: main.o $(OBJS)
	$(CXX) -o bliss main.o $(OBJS) $(LDFLAGS)

all: demo bliss

clean:
	rm -f *.o bliss histo.dat
