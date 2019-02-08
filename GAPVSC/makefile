CXX = g++
CXXFLAGS = -g -Wall -O -Wreorder
LDADD =

OBJ = main.o ACO.o BeamSearch.o Benders.o ConstructHeu.o Controller.o Corridor.o Ejection.o ForwardBackward.o GAP.o GeneticAlgorithm.o json.o Kernel.o Lagrangian.o LocalSearch.o LocBranching.o LowerBound.o MetaLocalSearch.o MIPCplex.o MTHG.o Persistence.o Rins.o SimAnnealing.o TabuSearch.o VLSN.o

%.o: ..\%.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

main: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o GAPVS 

clean:
	rm -rf $(OBJ) GAPVS

.PHONY: clean