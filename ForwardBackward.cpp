#include "ForwardBackward.h"

FandB::FandB(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

FandB::~FandB()
{
   //dtor
}

int FandB::forwardBackward(int** c, int delta, int maxNodes, bool fVerbose)
{  int i,j,z,iter,maxIter = 100;
   int currNode,openNodesF,openNodesB;

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
      bTree.push_back(list<int>());
      fList.push_back(list<int>());
      bList.push_back(list<int>());
      fTopCost.push_back(-1);
      bTopCost.push_back(-1);
   }

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   computeRegrets(c,n,m,regrets);

   for(j=0;j<n;j++) indCost[j] = j;       // support for ordering
   std::sort(indCost.begin(), indCost.end(), compRegr);

   node rootF, rootB;
   rootF.z = rootB.z = 0;
   rootF.dad = rootB.dad = 0;
   rootF.server = rootB.server = -1;
   rootF.capused.resize(m);
   rootB.capused.resize(m);
   for(i=0;i<m;i++) 
   {  rootF.capused[i]=0;
      rootB.capused[i]=0;
   }
   stack.push_back(rootF);
   stack.push_back(rootB);
   indLastNode = 1;
   zub = INT_MAX;
   openNodesF = openNodesB = 0;

   z = 0;
   currNode = 0;
   openNodesF += expandNode(c, indCost[0], -1, currNode, indCost, true);// stack initialization, forward
   fTree[0].push_back(currNode);                // append to list of expanded nodes

   z = 0;
   currNode   = 1;
   openNodesB += expandNode(c, indCost[n-1], n, currNode, indCost, false);// stack initialization, backward
   bTree[n-1].push_back(currNode);              // append to list of expanded nodes

   iter = 0;
   while ((openNodesF+openNodesB) > 0 && indLastNode < maxNodes && iter < maxIter)
   {
      cout << "Iter " << iter << " Num nodes: " << indLastNode << " open nodes forw." << openNodesF << " open nodes backw." << openNodesB << endl;
      openNodesF = sweepForward(c, delta,  maxNodes, openNodesF, indCost);
      openNodesB = sweepBackward(c, delta,  maxNodes, openNodesB, indCost);
      iter++;
   }

   if(abs(GAP->checkSol(solbest)-zub) > GAP->EPS)
   {  cout << "[forwardBackward]: Error, solution cost mismatch" << endl;
      z = INT_MAX;
   }
   else
      cout << "Construction terminated. zub = " << zub << endl;
   return zub;
}

// one run of fowrward beam search
int FandB::sweepForward(int** c, int delta, int maxNodes, int openNodes, vector<int> indCost)
{  int j,jlev,k;
   int currNode,newNodes,numExp;

   jlev = newNodes = 0;                        // got to initialize them 
   while(jlev<n)                               // main construction loop, could be while true
   {  jlev = findNextNodeF(jlev,newNodes,openNodes);     // backjunmping!
      if(jlev==n) break;

      newNodes = 0;                            // new nodes at corrent level
      for(k=0;k<delta && fList[jlev].size()>0;k++)       // EXPANSION
      {  currNode = fList[jlev].front();                 // gets the first element
         if(jlev < n && stack[currNode].z < zub)
         {  j = (jlev == n-1 ? -1 : indCost[jlev+1]);    // client order by regrets
            numExp    = expandNode(c, j, jlev, currNode, indCost, true);
            openNodes += numExp;
            newNodes  += numExp;
         }
         if(indLastNode > maxNodes)
         {  cout << "node limit reached" << endl;
            goto end;   
         }
         fTree[jlev].push_back(currNode);       // append to list of expanded nodes
         fList[jlev].pop_front();               // remove from list of nodes to expand
         openNodes--;
         if(stack[currNode].z > fTopCost[jlev]) 
            fTopCost[jlev] = stack[currNode].z; // update max cost of expanded node at the level
         else
            if(isVerbose)
               if(stack[currNode].z < fTopCost[jlev]) cout << "[sweepForward] inner cost insertion" << endl;
      }
      if(isVerbose)
         cout << "[sweepForward] Level " << jlev << " expanded " << k << " new nodes " << newNodes << " open nodes " << openNodes << " tot nodes "<< indLastNode << endl;
   }
end:
   return openNodes;
}

// one run of backward beam search
int FandB::sweepBackward(int** c, int delta, int maxNodes, int openNodes, vector<int> indCost)
{  int j,jlev,k;
   int currNode,newNodes=0,numExp;

   jlev = n-1;
   while(jlev>=0)                                // main construction loop
   {  jlev = findNextNodeB(jlev,newNodes,openNodes); // backjunmping!
      if(jlev < 0) break;

      newNodes = 0;
      for(k=0;k<delta && bList[jlev].size()>0;k++)
      {  currNode = bList[jlev].front();               // gets the first element
         if(jlev >= 0 && stack[currNode].z < zub)
         {  j = (jlev == 0 ? -1 : indCost[jlev-1]);    // client order by regrets
            numExp    = expandNode(c, j, jlev, currNode, indCost, false);
            openNodes += numExp;
            newNodes  += numExp;
         }
         if(indLastNode > maxNodes) goto end;   // node limit reached
         bTree[jlev].push_back(currNode);       // append to list of expanded nodes
         bList[jlev].pop_front();               // remove from list of nodes to expand
         openNodes--;
         if(stack[currNode].z > bTopCost[jlev])
            bTopCost[jlev] = stack[currNode].z; // update max cost of expanded node at the level
         else
            if(isVerbose)
               if(stack[currNode].z < bTopCost[jlev]) cout << "[sweepBackward] inner cost insertion" << endl;
      }
      if(isVerbose)
         cout << "[sweepBackward] Level " << jlev << " expanded " << k << " new nodes " << newNodes << " open nodes " << openNodes << " tot nodes "<< indLastNode << endl;
   }
end:
   return openNodes;
}

// c: costs; j: current client; jlev; current tree level; corrNode: node currently expanded (will be father)
int FandB::expandNode(int** c, int j, int jlev, int currNode, vector<int> indCost, bool isForward)
{  int i,ii,numNewNodes=0,z;
   vector<int> cost(m),indReq(m);
   auto compCost = [&cost](int a, int b){ return cost[a] < cost[b]; };  // ASC order

   if(isForward && j<0)            // reached complete solutions (no dad)
   {  if(stack[currNode].z < zub)
         z=readSolutionF(currNode,indCost);
      goto end;
   }
   if(!isForward && j<0)           // reached complete solutions (no dad)
   {  if(stack[currNode].z < zub)
         z=readSolutionB(currNode,indCost);
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
         if(isForward)
         {  insertInOrder(fList[jlev+1], indLastNode);   // inserts the index of the node in the level list
            checkMatch(jlev+1, indLastNode, isForward, indCost);       // check for feasible completion
            if(isVerbose)   
               cout << "lev." << jlev+1 << " node " << indLastNode << " cost " << newNode.z << endl;
         }
         else
         {  insertInOrder(bList[jlev-1], indLastNode);   // inserts the index of the node in the level list
            checkMatch(jlev-1, indLastNode, isForward, indCost);       // check for feasible completion
            if(isVerbose)   
               cout << "lev." << jlev-1 << " node " << indLastNode << " cost " << newNode.z << endl;
         }
         numNewNodes++;
      }
      ii++;
   }
end:
   return numNewNodes;  // new nodes opened at the level 
}

// level: level where to insert; elem: element to insert (key: node cost)
int FandB::insertInOrder(list<int> & lst, int elem)
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

// check for matching partial solutions, jlev level of indLastNode
int FandB::checkMatch(int jlev, int indLastNode, bool isForward, vector<int> indCost)
{  int i,z=-1,res=0;
   list<int> * lstCompletions;
   std::list<int>::const_iterator iterator;

   if(jlev==n-1 || jlev==0) goto end;

   if(isForward)
      lstCompletions = &bTree[jlev+1];
   else
      lstCompletions = &fTree[jlev-1];

   // Iterate and print values of the completion list
   for (iterator = (*lstCompletions).begin(); iterator != (*lstCompletions).end(); ++iterator) 
   {  if(isVerbose) cout << (isForward ? "F" : "B") << " lev." <<  jlev << " node compl " << *iterator << endl;

      if(stack[*iterator].server < 0) 
         continue;
      for (i=0; i<m; i++)
         if(stack[indLastNode].capused[i] + stack[*iterator].capused[i] > GAP->cap[i])
            goto next;
      z = stack[indLastNode].z + stack[*iterator].z;
      if(isVerbose) cout << "FEASIBLE! cost " << z << endl;

      if(z < zub)
      {  int jLevF = (isForward ? jlev : jlev-1);
         int fNode = (isForward ? indLastNode : *iterator);
         int bNode = (isForward ? *iterator : indLastNode);
         if(isVerbose)
            cout << "f node "<< indLastNode << " z_f " << stack[indLastNode].z << 
               " b node "<< *iterator  << " z_b " << stack[*iterator].z << " zub " << zub << endl;
         z=readSolutionFB(jLevF,fNode,bNode,indCost);
      }
next: ;
   }

end:
   return z;
}

// reads the solutions from a last node of the forward tree
int FandB::readSolutionF(int currNode, vector<int> indCost)
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
         //cout << "node " << solNode << " j " << j << " i " << sol[j] << " c " << GAP->c[sol[j]][j] << " z " << stack[solNode].z << endl;
         solNode = stack[solNode].dad;
         jlev--;
      }
   }

   return res;
}

// reads the solutions from a last node of the backward tree
int FandB::readSolutionB(int currNode, vector<int> indCost)
{  int res = 0,j,jlev,solNode;

   if(stack[currNode].z < zub)
   {  zub  = stack[currNode].z;
      cout << "New zub: " << zub << endl;
      // reconstruction of the solution
      solNode = currNode;
      jlev = 0;
      while(jlev < n)
      {  j = indCost[jlev];
         sol[j] = solbest[j] = stack[solNode].server;
         //cout << "node " << solNode << " j " << j << " i " << sol[j] << " c " << GAP->c[sol[j]][j] << " z " << stack[solNode].z << endl;
         solNode = stack[solNode].dad;
         jlev++;
      }
   }

   return res;
}

// reads a solution as a mix of forw and a backw partials, jLevF last level of forward tree
int FandB::readSolutionFB(int jLevF, int fNode, int bNode, vector<int> indCost)
{  int res = 0,z=0,j,jlev,solNode;

   zub  = stack[fNode].z + stack[bNode].z;
   cout << "New zub: " << zub << endl;

   // reconstruction of the forward part of the solution
   solNode = fNode;
   jlev = jLevF;
   while(jlev > -1 && stack[solNode].server >= 0)
   {  j = indCost[jlev];
      sol[j] = solbest[j] = stack[solNode].server;
      z += GAP->c[sol[j]][j];
      if(isVerbose) cout << "node " << solNode << " lev " << jlev << " j " << j << " i " << sol[j] << " c " << GAP->c[sol[j]][j] << " z " << stack[solNode].z << " ztot " << z << endl;
      solNode = stack[solNode].dad;
      jlev--;
   }

   // reconstruction of the backward part of the solution
   solNode = bNode;
   jlev = jLevF+1;
   while(jlev < n && stack[solNode].server >= 0)
   {  j = indCost[jlev];
      sol[j] = solbest[j] = stack[solNode].server;
      z += GAP->c[sol[j]][j];
      if(isVerbose) cout << "node " << solNode << " lev " << jlev << " j " << j << " i " << sol[j] << " c " << GAP->c[sol[j]][j] << " z " << stack[solNode].z << " ztot " << z << endl;
      solNode = stack[solNode].dad;
      jlev++;
   }
   if(abs(z-zub) > GAP->EPS)
      cout << "[readSolutionFB]: ---------------------------- Error, solution cost mismatch" << endl;

   return z;
}

// finds the next node to expand forward
int FandB::findNextNodeF(int jlev, int newNodes, int openNodes)
{  int jmin;

   if(newNodes>0 || jlev == n-1)    // if there were expansions, go on
      jlev++;
   else              // find the highest level with unexpanded nodes
   {  for(jmin=0;jmin<n;jmin++)
         if(fList[jmin].size()>0)
            break; 
      jlev = jmin;
      //jlev = n;   // just a check
   }
   return jlev;
}

// finds the next node to expand backward
int FandB::findNextNodeB(int jlev, int newNodes, int openNodes)
{  int jmin;

   if(newNodes>0 || jlev == 0)    // if there were expansions, go on
      jlev--;
   else              // find the lowest level with unexpanded nodes
   {  for(jmin=n-1;jmin>0;jmin--)
         if(bList[jmin].size()>0)
            break; 
      if(jlev < n-1) 
         jlev = -1;
      else
         jlev = jmin;
   }
   return jlev;
}
