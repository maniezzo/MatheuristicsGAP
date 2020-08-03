#ifndef FANDB_H
#define FANDB_H
#include "GAP.h"
#include <list>
#include <algorithm>    // std::copy

class FandB
{
   public:
      GeneralizedAssignemnt* GAP;

      FandB(GeneralizedAssignemnt*, int&);
      ~FandB();
      int forwardBackward(int** c, int delta, int maxNodes, bool fVerbose);

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
         int client; // the assigned client
         int server; // who the client is assigned to
         int dad;    // which sol was expanded into this
         vector<int> capused;  // array of used capacities
      };

      vector<node> stack;   // stack all nodes expanded during search
      vector<list<int>> fTree;   // forward tree, one list for each level / customer
      vector<list<int>> bTree;   // backward tree, one list for each level / customer
      vector<list<int>> fList;   // the list of still unexpanded nodes at each level of the forward tree
      vector<list<int>> bList;   // the list of still unexpanded nodes at each level of the backward tree
      vector<int> fTopCost; // max cost of expanded node at each level of the forward tree
      vector<int> bTopCost; // max cost of expanded node at each level of the backward tree

      int indLastNode;  // aka stack.size(). just it
      int numFathomed;  // num fahtomed nodes

      int sweepForward(ofstream&, int** c, int delta, int maxNodes, int openNodes, vector<int> indCost);
      int sweepBackward(ofstream&, int** c, int delta, int maxNodes, int openNodes, vector<int> indCost);
      int expandNode(ofstream&, int** c, int j, int jlev, int currNode, vector<int> indCost, bool isForward);  // generates feasible offspring of a node
      int insertInOrder(list<int> & lst, int ind);              // inserts a stack index in a list, ordered on a key
      int readSolutionF(ofstream&, int currNode, vector<int> indCost);     // reads the solutions from a last node of the forward tree
      int readSolutionB(ofstream&, int currNode, vector<int> indCost);     // reads the solutions from a last node of the backward tree
      int readSolutionFB(ofstream&, int jLevF, int fNode, int bNode, vector<int> indCost); // reads the solution as a mix of forw and a backw partials
      int findNextNodeF(int jlev, int newNodes, int openNodes); // finds the next node to expand forward
      int findNextNodeB(int jlev, int newNodes, int openNodes); // finds the next node to expand backward
      int checkMatch(ofstream&, int jlev, int indLastNode, bool isForward, vector<int> indCost);// checks for matching partial solutions
};

#endif // FANDB_H
