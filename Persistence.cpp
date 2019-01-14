#include "Persistence.h"
#include <stdexcept>

Persistence::Persistence()
{
   //ctor
}

Persistence::~Persistence()
{
   //dtor
}

// loading parameters of all supported algorithms
Config* Persistence::loadConfig()
{
   string infile,line;

   infile = "config.json";
   cout << "Opening " << infile << endl;

   ifstream jData (infile.c_str());
   std::stringstream buffer;
   buffer << jData.rdbuf();
   line = buffer.str();
   //std::getline(jData,line);
   jData.close();

   json::Value JSV = json::Deserialize(line);

   json::Object jSA = JSV["SA"];
   GAP->conf->SA->maxT        = (int) jSA["maxT"];
   GAP->conf->SA->k           = jSA["k"];
   GAP->conf->SA->alpha       = jSA["alpha"];
   GAP->conf->SA->iterAnneal  = jSA["iterAnneal"];
   GAP->conf->SA->coefAnneal  = jSA["coefAnneal"];
   GAP->conf->SA->maxiter = jSA["maxiter"];
   GAP->conf->SA->fGoMath = jSA["fGoMath"];

   json::Object jTS = JSV["TS"];
   GAP->conf->TS->Ttenure = jTS["Ttenure"];
   GAP->conf->TS->maxiter = jTS["maxiter"];
   GAP->conf->TS->fGoMath = jTS["fGoMath"];

   json::Object jILS = JSV["ILS"];
   GAP->conf->IterLS->maxiter = jILS["maxiter"];
   GAP->conf->IterLS->alpha = jILS["alpha"];

   json::Object jGRASP = JSV["GRASP"];
   GAP->conf->GRASP->maxiter = jGRASP["maxiter"];
   GAP->conf->GRASP->candNum = jGRASP["candNum"];

   json::Object jVNS = JSV["VNS"];
   GAP->conf->VNS->maxiter = jVNS["maxiter"];
   GAP->conf->VNS->fGoMath = jVNS["fGoMath"];

   json::Object jGA = JSV["GA"];
   GAP->conf->GA->maxiter = jGA["maxiter"];
   GAP->conf->GA->numpop = jGA["numpop"];
   GAP->conf->GA->pc = jGA["pc"];

   json::Object jEC = JSV["EC"];
   GAP->conf->EC->maxiter = jEC["maxiter"];

   json::Object jLagrAss = JSV["lagrAss"];
   GAP->conf->lagrAss->alpha     = jLagrAss["alpha"];
   GAP->conf->lagrAss->alphastep = jLagrAss["alphastep"];
   GAP->conf->lagrAss->minalpha  = jLagrAss["minalpha"];
   GAP->conf->lagrAss->innerIter = jLagrAss["innerIter"];
   GAP->conf->lagrAss->maxiter   = jLagrAss["maxiter"];

   json::Object jlagrCap = JSV["lagrCap"];
   GAP->conf->lagrCap->alpha     = jlagrCap["alpha"];
   GAP->conf->lagrCap->alphastep = jlagrCap["alphastep"];
   GAP->conf->lagrAss->minalpha  = jLagrAss["minalpha"];
   GAP->conf->lagrCap->innerIter = jlagrCap["innerIter"];
   GAP->conf->lagrCap->maxiter   = jlagrCap["maxiter"];

   json::Object jRINS = JSV["RINS"];
   GAP->conf->rinsConf->maxnodes = jRINS["maxnodes"];

   json::Object jbeam = JSV["beam"];
   GAP->conf->beam->delta    = jbeam["delta"];
   GAP->conf->beam->maxnodes = jbeam["maxnodes"];

   json::Object jFB = JSV["FandB"];
   GAP->conf->fbConf->delta    = jFB["delta"];
   GAP->conf->fbConf->maxnodes = jFB["maxnodes"];

   json::Object jCorridor = JSV["corridor"];
   GAP->conf->corridorConf->delta = jCorridor["delta"];
   GAP->conf->corridorConf->maxiter = jCorridor["maxiter"];

   json::Object jLocBranch = JSV["locBranching"];
   GAP->conf->locBranching->k = jLocBranch["k"];
   GAP->conf->locBranching->maxiter = jLocBranch["maxiter"];

   json::Object jBenders = JSV["benders"];
   GAP->conf->bendersConf->maxiter = jBenders["maxiter"];

   json::Object jVLSN = JSV["VLSN"];
   GAP->conf->verylargeConf->maxiter = jVLSN["maxiter"];
   GAP->conf->verylargeConf->k       = jVLSN["k"];

   GAP->conf->datafile  = JSV["datafile"];
   GAP->conf->isVerbose = JSV["isverbose"];
   GAP->conf->aversionf = JSV["aversionf"];
   return GAP->conf;
}

// reads instance data from json formatted files
int Persistence::readJSONdata(string fileList)
{
   string infile,line,path;
   size_t i,j,cont;
   vector <string> arrFiles;

   infile = fileList;
   cout << "Opening " << infile << endl;
   size_t found = fileList.find_last_of("/\\");
   path = fileList.substr(0, found);
   cout << "Data path: " << path << '\n';

   try
   {
      ifstream fList;
      fList.open(infile.c_str(), std::ifstream::in);
      fList.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit);
      string str;
      cont=0;
      while ( std::getline(fList,str) && str.length()>0)
      {  arrFiles.push_back(path+"\\"+str) ;
         cout << cont <<") " << arrFiles[cont] << endl;
         cont++;
      }
      fList.close();
   }
   catch (std::exception const& e)
   {  cout << "Error: " << e.what() << endl;
      return -1;
   }

   // here I got the list of filenames correctly
   try
   {
      cout << "Reading " << arrFiles[0] << endl;
      ifstream jData;
      jData.open (arrFiles[0], std::ifstream::in);
      jData.exceptions ( ifstream::eofbit | ifstream::failbit | ifstream::badbit );
      getline(jData,line);
      //cout << "line:" << line << endl;
      jData.close();
   }
   catch(std::exception const& e)
   {  cout << "Error: " << e.what() << endl;
      return -1;
   }

   json::Value JSV = json::Deserialize(line);
   GAP->name = JSV["name"];
   GAP->n = JSV["numcli"];
   GAP->m = JSV["numserv"];
   GAP->cap = (int*) malloc(GAP->m * sizeof(int));
   for(i=0;i<JSV["cap"].size();i++)
      GAP->cap[i] = JSV["cap"][i];

   GAP->c = (int**) malloc(GAP->m * sizeof(int *));
   for(i=0;i<JSV["cost"].size();i++)
   {  GAP->c[i] = (int*) malloc(GAP->n * sizeof(int));
      for(j=0;j<JSV["cost"][i].size();j++)
         GAP->c[i][j] = JSV["cost"][i][j];
   }

   GAP->req = (int**) malloc(GAP->m * sizeof(int *));
   for(i=0;i<JSV["req"].size();i++)
   {  GAP->req[i] = (int*) malloc(GAP->n * sizeof(int));
      for(j=0;j<JSV["req"][i].size();j++)
         GAP->req[i][j] = JSV["req"][i][j];
   }

   GAP->zub = INT32_MAX;

   cout << "JSON data read" << endl;;
   return 1;
}
