#ifndef PERSISTENCE_H
#define PERSISTENCE_H
#include <windows.h> // GetModuleFileName
#include <fstream>   // ifstream
#include <sstream>   // stringstream
#include <string>
#include "json.h"
#include "GAP.h"

class Persistence
{
   public:
      Persistence();
      ~Persistence();
      Config* loadConfig();
      int readJSONdata(string);

      GeneralizedAssignemnt* GAP;
   protected:
   private:
};

#endif // PERSISTENCE_H
