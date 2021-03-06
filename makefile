OBJECTS = hard_sphere.o bop.o g_r.o ini_pos.o vector.o  atom.o
FLAGS   = -std=c++11
CXXFLAGS   = -std=c++11
a.out : $(OBJECTS) 
	g++ -o z.out $(OBJECTS) 
hard_sphere.o : MC_header.hpp ini_pos.hpp g_r.hpp
bop.o : atom.hpp bop.hpp
g_r.o : g_r.hpp atom.hpp
ini_pos.o : ini_pos.hpp atom.hpp
vector.o : atom.hpp vector.hpp
atom.o : atom.hpp

