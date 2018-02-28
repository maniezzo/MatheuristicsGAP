#ifndef CONFIG_H
#define CONFIG_H
#include <string>

using namespace std;

class Config
{
   public:
   class SimAnn
   {  public:
         double maxT;
         double k;
         int maxiter;
   };

   class Tabu
   {  public int Ttenure;
      public int maxiter;
   }

   class LagrAss
   {
      public:
         double alpha;
         double alphastep;
         int maxiter;
         int innerIter;
   }

   class LagrCap
   {
      public:
         double alpha;
         double alphastep;
         int maxiter;
         int innerIter;
   }

   class Corridor
   {  public:
         int delta;
         int maxiter;
   }

   class LocBranching
   {  public:
         int k;
         int maxiter;
   }

   SimAnn* SA;
   Tabu*   TS;
   LagrAss*  lagrAss;
   LagrCap*  lagrCap;
   Corridor* corridor;
   LocBranching*  locBranching;
   string datafile;
   int isverbose;
}

#endif // CONFIG_H
