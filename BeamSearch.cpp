#include "BeamSearch.h"

BeamSearch::BeamSearch(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol       = GAP->sol;
   solbest   = GAP->solbest;
   req       = GAP->req;
   isVerbose = false;
}

BeamSearch::~BeamSearch()
{
   //dtor
}

int BeamSearch::beamSearch(int** c, int delta, int maxNodes, bool fVerbose)
{  int i,j,jlev,k,z;
   int currNode,newNodes,openNodes,numExp;

   vector<int> regrets(n),capleft(m);
   vector<int> indCost(n);           // index of order of expansion of the clients
   auto compRegr = [&regrets](int a, int b){ return regrets[a] > regrets[b]; };  // DESC order

   if(GAP->n == NULL)
   {  cout << "Instance undefined. Exiting" << endl;
      return INT_MAX;
   }
   isVerbose = fVerbose;

   for(j=0;j<n;j++)
   {  fTree.push_back(list<int>());
      fList.push_back(list<int>());
   }

   node root;
   root.z       = 0;
   root.dad     = 0;
   root.server  = -1;
   root.capused.resize(m);
   for(i=0;i<m;i++) root.capused[i]=0;
   stack.push_back(root);
   indLastNode = 0;
   zub = INT_MAX;

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   computeRegrets(c,n,m,regrets);

   for(j=0;j<n;j++) indCost[j] = j;       // support for ordering
   std::sort(indCost.begin(), indCost.end(), compRegr);

   z = 0;
   currNode = openNodes = 0;
   openNodes += expandNode(c, indCost[0],-1, currNode, indCost);// stack initialization
   fTree[0].push_back(currNode);                // append to list of expanded nodes

loop:
   jlev = newNodes = 0;
   while(jlev<n)                                // main construction loop
   {  jlev = findNextNode(jlev,newNodes,openNodes);     // backjunmping!
      if(jlev==n) break;
      newNodes = 0;
      
      for(k=0;k<delta && fList[jlev].size()>0;k++)
      {  currNode = fList[jlev].front();        // gets the first element
         if(jlev < n && stack[currNode].z < zub)
         {  j = (jlev == n-1 ? -1 : indCost[jlev+1]);    // client order by regrets
            numExp    = expandNode(c, j, jlev, currNode, indCost);
            openNodes += numExp;
            newNodes  += numExp;
         }
         if(indLastNode > maxNodes) goto end;   // node limit reached
         fTree[jlev].push_back(currNode);       // append to list of expanded nodes
         fList[jlev].pop_front();               // remove from list of nodes to expand
         openNodes--;
      }
      if(isVerbose)
         cout << "Level " << jlev << " expanded " << k << " new nodes " << newNodes << " open nodes " << openNodes << " tot nodes "<< indLastNode << endl;
      if(indLastNode % 1000 < 10)
         cout << "Nodes " << indLastNode << " open " << openNodes << " zub " << zub << endl;
   }
   if(openNodes > 0) goto loop;           // should be unneeded, but just in case

end:
   if(abs(GAP->checkSol(solbest)-zub) > GAP->EPS)
   {  cout << "[beamSearch]: Error, solution cost mismatch" << endl;
      z = INT_MAX;
   }
   else
      cout << "Construction terminated. zub = " << zub << endl;
   return zub;
}

// c: costs; j: current client; jlev; current tree level; corrNode: node currently expanded (will be father)
int BeamSearch::expandNode(int** c, int j, int jlev, int currNode, vector<int> indCost)
{  int i,ii,numNewNodes=0,z=-1;
   vector<int> cost(m),indReq(m);
   auto compCost = [&cost](int a, int b){ return cost[a] < cost[b]; };  // ASC order

   if(j<0)              // reached complete solutions
   {  if(stack[currNode].z < zub)
         z=readSolution(currNode,indCost);
      goto end;
   }

   // cost of expansions
   for(i=0;i<m;i++)
   {  cost[i]= GAP->c[i][j];
      indReq[i] = i;
   }
   std::sort(indReq.begin(), indReq.end(), compCost);

   // expansion of a node
   ii = numNewNodes = 0;
   while(ii<m)
   {  i=indReq[ii];

      if( (stack[currNode].capused[i] + GAP->req[i][j]) <= GAP->cap[i] &&
          (stack[currNode].z + GAP->c[i][j]) < zub)
      {  
         node newNode;
         for(int ii=0;ii<m;ii++) newNode.capused.push_back(stack[currNode].capused[ii]);
         newNode.capused[i] = stack[currNode].capused[i] + GAP->req[i][j];
         newNode.dad = currNode;
         newNode.server = i;
         newNode.z = stack[currNode].z + GAP->c[i][j];
         stack.push_back(newNode);
         indLastNode++;
         insertInOrder(fList[jlev+1], indLastNode);   // inserts the index of the node in the level list
         numNewNodes++;
      }
      ii++;
   }

end:
   return numNewNodes;  // new nodes opened at the level 
}

// level: level where to insert; elem: element to insert (key: node cost)
int BeamSearch::insertInOrder(list<int> & lst, int elem)
{  int res = 0;
   list<int>::iterator it;

   if(lst.size() == 0)
      lst.push_back(elem);
   else
   {  it = lst.begin();
      while(it != lst.end() && stack[*it].z < stack[elem].z)
         ++it;
      lst.insert(it,elem);
   }
   return res;
}

// reads the solutions from the last line of the tree
int BeamSearch::readSolution(int currNode, vector<int> indCost)
{  int res = 0,j,jlev,solNode;

   if(stack[currNode].z < zub)
   {  zub  = stack[currNode].z;
      cout << "New zub: " << zub << endl;
      // reconstruction of the solution
      solNode = currNode;
      jlev = n - 1;
      while(stack[solNode].server > -1)
      {  j = indCost[jlev];
         sol[j] = solbest[j] = stack[solNode].server;
         if(isVerbose)
            cout << "node " << solNode << " j " << j << " i " << sol[j] << " c " << GAP->c[sol[j]][j] << " z " << stack[solNode].z << endl;
         solNode = stack[solNode].dad;
         jlev--;
      }
   }

   return res;
}

// finds the next node to expand
int BeamSearch::findNextNode(int jlev, int newNodes, int openNodes)
{  int jmin;

   if(newNodes>0)    // if there were expansions, go on
      jlev++;
   else              // find the highest level with unexpanded nodes
   {  for(jmin=0;jmin<n;jmin++)
         if(fList[jmin].size()>0)
            break; 
      jlev = jmin;
   }
   return jlev;
}