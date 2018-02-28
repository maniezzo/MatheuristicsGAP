#ifndef BEAM_H
#define BEAM_H
#include "GAP.h"
#include <list>
#include <algorithm>    // std::copy

class BeamSearch
{
   public:
      GeneralizedAssignemnt* GAP;

      BeamSearch(GeneralizedAssignemnt*, int&);
      ~BeamSearch();
      int beamSearch(int** c, int delta, int maxNodes, bool fVerbose); // not nice names with small and big cap, isn't it?

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
      bool  isVerbose;

      // the state of each partial solution
      struct node
      {  int z;      // cost of the partial soluton
         int server; // who the client is assigned to
         int dad;    // which sol was expanded into this
         vector<int> capused;  // array of used capacities
      };

      vector<node> stack;  // stack all nodes expanded during search
      vector<list<int>> fTree;   // forward tree, one list for each level / customer
      vector<list<int>> fList;   // the list of still unexpanded nodes at each level of the forward tree

      int indLastNode;  // aka stack.size(). just it

      int expandNode(int** c, int j, int jlev, int currNode, vector<int> indCost);
      int insertInOrder(list<int> & lst, int ind); // insert a stack index in a list, ordered on a key
      int readSolution(int currNode, vector<int> indCost);  // reads the solutions from the last line of the table
      int findNextNode(int jlev, int newNodes, int openNodes); // finds the next node to expand
};

#endif // BEAM_H
