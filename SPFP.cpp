#include "headersMemory.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>


#ifndef GLOBALS
#define GLOBALS

double MIN_X = -124.848974;
double MAX_X = -66.885444;
double MIN_Y = 24.396308;
double MAX_Y = 49.384358;

// double MIN_X = -124.848974;
// double MAX_X = -66.885444;
// double MIN_Y = 10.396308;
// double MAX_Y = 65.384358;

double DATASET_SIZE = 0;
double DELTA_X = 0;
double DELTA_Y = 0;
int MAXSC = 0;                // maximum number of friends plus one
double MAXDIST = 0;          // maximum distance between to points
int MAXT = 0;
double ALPHA = 0;
double GAMMA = 0;
double BETA = 0;
#endif

void check_cooccurrence(SPOs *spos, int u1id, int u2id, int l1id, int l2id, bool validTimeRange,
  set<int>* seenLocations,
  unordered_set< pair<int,int>, PairHasher >* seen_pairs,
  unordered_set< pair<int,int>, PairHasher >* seen_coocc,
  bool isOptimistic){

  bool areFriends = spos->areFriends(u1id, u2id);

  bool notRecorededFriendship = seen_pairs->find(make_pair(u1id, u2id)) == seen_pairs->end();

  bool notRecorededCoOcc = seen_coocc->find(make_pair(u1id, u2id)) == seen_coocc->end();

  bool unseenLocation = (seenLocations->find(u2id) == seenLocations->end());

  bool isDifferent = l1id != l2id && u1id != u2id;

  bool validCoOccurrence = isDifferent && notRecorededCoOcc && (isOptimistic || unseenLocation) && validTimeRange;

  if( validCoOccurrence ){
    if(  areFriends && notRecorededFriendship ){
      seen_pairs->insert(make_pair(u1id, u2id));
      seenLocations->insert( l2id );
    }

    seen_coocc->insert(make_pair(u1id, u2id));
  }
}

void count_cooccurences(SPOs *spos, GPOs* gpos, double radius, double timeRange, bool isOptimistic){
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  double radius_geo_dist = (radius/1000) * 360 / EARTH_CIRCUMFERENCE,
         x=0, y=0;
  unsigned int lid, uid, count=0;

  unsigned int tp = 0, pos = 0, nFriendships=spos->edges;
  double precision, recall;

  boost::posix_time::ptime time;
  boost::posix_time::time_duration td;

  unordered_set< pair<int,int>, PairHasher > seen_pairs, seen_coocc;

  set<int>* seenLocations = new set<int>();
  bool validTimeRange;

  for(auto l = gpos->locations.begin(); l != gpos->locations.end(); l++){
    x   = l->second->getX();
    y   = l->second->getY();
    lid = l->second->getID();
    uid = l->second->getUID();
    time = l->second->getTime();


    // Skipping range search if radius is 0
    if(radius == 0){
      auto points =  gpos->location_to_user.find(lid)->second;
      for(auto u = points->begin(); u != points->end(); u++){
        td = time - (*u)->getTime();
        validTimeRange = (abs(td.total_seconds()) <= timeRange);
        // cout << uid << " " << (*u)->getUID() << " " << lid << " " << (*u)->getID() << " " << validTimeRange << endl;
        check_cooccurrence(spos, uid, (*u)->getUID(), lid, -1,
          validTimeRange,
          seenLocations, &seen_pairs, &seen_coocc,
          isOptimistic);
      }
    } else {
      vector<res_point*>* checkins = gpos->getRange(x, y, radius_geo_dist);
      for(auto c = checkins->begin(); c != checkins->end(); c++){
        td = time - (*c)->time;
        validTimeRange = (abs(td.total_seconds()) <= timeRange);
        check_cooccurrence(spos, uid, (*c)->uid, lid, (*c)->id,
          validTimeRange,
          seenLocations, &seen_pairs, &seen_coocc,
          isOptimistic);
        delete (*c);
      }
      delete checkins;
    }

    count++;

    if(count % 100000==0)
      cout << count << " " << seen_pairs.size() << endl;
  }

  cout << "Number of locations : " << gpos->location_to_user.size() << endl;
  cout << "Number of users : " << gpos->user_to_location.size() << endl;

  tp  = seen_pairs.size();
  pos = seen_coocc.size();

  cout << "Number of friendships that can be inferred from co-occurrences " << tp << endl;
  cout << "Number of co-occurrences " << pos << endl;
  cout << "Number of friendships " << nFriendships << endl;


  precision = tp / (double) pos;
  recall    = tp / (double) nFriendships;

  cout << "Precision : " << precision << endl;
  cout << "Recall : " << recall << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

void testpTools(){

  // int i;
  // double length, entropyTarget;
  double firstEntropy, secondEntropy, thirdEntropy, targetEntropy;
  // double firstMItarget, secondMItarget, thirdMItarget, targetMItarget;
  // int *testFirstVector, *testSecondVector, *testThirdVector, *testMergedVector;
  // struct timeval start,end;

  uint *firstVector = (uint *) calloc(4,sizeof(uint));
  uint *secondVector = (uint *) calloc(4,sizeof(uint));
  uint *thirdVector = (uint *) calloc(4,sizeof(uint));
  uint *targetVector = (uint *) calloc(4,sizeof(uint));

  firstVector[0] = 0;
  firstVector[1] = 0;
  firstVector[2] = 1;
  firstVector[3] = 1;

  secondVector[0] = 0;
  secondVector[1] = 1;
  secondVector[2] = 0;
  secondVector[3] = 1;

  thirdVector[0] = 0;
  thirdVector[1] = 1;
  thirdVector[2] = 1;
  thirdVector[3] = 1;

  targetVector[0] = 0;
  targetVector[1] = 1;
  targetVector[2] = 1;
  targetVector[3] = 0;

  firstEntropy = calcEntropy(firstVector,4);
  secondEntropy = calcEntropy(secondVector,4);
  thirdEntropy = calcEntropy(thirdVector,4);
  targetEntropy = calcEntropy(targetVector,4);
  // targetEntropy = calcRenyiEntropy(0.1,targetVector,4 );

  printf("Entropies - first: %f, second: %f, third: %f, target %f\n",firstEntropy,secondEntropy,thirdEntropy,targetEntropy);

  //verifying section 3.4 EBM

  cout<<" -----------------------Entropy Example ------------------------"<<endl;
  cout<<" ------------------------Section  3.4---------------------------"<<endl;
  uint *CabOcc = (uint *) calloc(5,sizeof(uint));
  uint *CcdOccur = (uint *) calloc(5,sizeof(uint));

  CabOcc[0] = 10;
  CabOcc[1] = 1;
  CabOcc[2] = 0;
  CabOcc[3] = 0;
  CabOcc[4] = 9;

  CcdOccur[0] = 2;
  CcdOccur[1] = 3;
  CcdOccur[2] = 2;
  CcdOccur[3] = 2;
  CcdOccur[4] = 3;

  cout<<" C_ab is : 10 , 1, 0, 0, 9"<<endl;
  cout<<" C_cd is : 2 , 3, 2, 2, 3"<<endl;

  cout<<" ------------------Testing Shannon Entropy------------------------"<<endl;
  double CabEntropy, CcdEntropy, CabDiversity, CcdDiversity;

  CabEntropy = calcEntropyFromCoV(CabOcc,5);
  CcdEntropy = calcEntropyFromCoV(CcdOccur,5);
  CabDiversity = exp(CabEntropy);
  CcdDiversity = exp(CcdEntropy);


  cout<<"CabEntropy = "<<CabEntropy<<" CcdEntropy = "<<CcdEntropy <<" CabEntropy/CcdEntropy = " << CabEntropy/CcdEntropy <<endl;
  cout<<"CabDiversity = "<<CabDiversity<<" CabDiversity = "<<CcdDiversity<<" CabDiversity/CcdDiversity = " << CabDiversity/CcdDiversity <<endl;

  cout<<" ------------------Testing Renyi Entropy------------------------"<<endl;
  double CabRenEntropy, CcdRenEntropy, CabRenDiversity, CcdRenDiversity;

  CabRenEntropy = calcRenyiEntropyFromCoV(0.5,CabOcc,5);
  CcdRenEntropy = calcRenyiEntropyFromCoV(0.5,CcdOccur,5);
  CabRenDiversity = exp(CabRenEntropy);
  CcdRenDiversity = exp(CcdRenEntropy);

  cout<<"CabRenEntropy = "<<CabRenEntropy<<" CcdRenEntropy = "<<CcdRenEntropy<< " CabRenEntropy/CcdRenEntropy = " << CabRenEntropy/CcdRenEntropy<<endl;
  cout<<"CabRenDiversity = "<<CabRenDiversity<<" CcdRenDiversity = "<<CcdRenDiversity<<" CabRenDiversity/CcdRenDiversity = "<< CabRenDiversity/CcdRenDiversity << endl;
}


int main(int argc, char *argv[]){
  cout.precision(15);

  DELTA_X = ((MAX_X - MIN_X)/ (X-1));
  DELTA_Y = ((MAX_Y - MIN_Y)/ (Y-1));

  double r1, r2, tp, tR;
  r1 = atof(argv[3]);
  r2 = atof(argv[4]);
  tR = atof(argv[5]);
  tp = atof(argv[6]);
  bool isOptimistic = (tp == 0);
  ALPHA = 0.480;
  BETA = 0.520;

  // testpTools();
  // exit(-1);

  // Loading social network and checkins
  cout << "------------- Loading SocialGraph ---------------" << endl;
  SPOs* spos = new SPOs();
  spos->load(argv[1]);

  cout << "------------- Loading checkins ---------------" << endl;
  GPOs* gpos = new GPOs(argv[2]);
  cout << "------------- Loading complete ---------------" << endl;


  // cout << "----- Loading Cooccurrence Matrix --- " << endl;
  // gpos->countU2UCoOccurrences();
  // cout << "----- Completed Loading Cooccurrence Matrix --- " << endl;


  // cout << "----- Calculating Location Entropy --- " << endl;
  // gpos->calculateLocationEntropy();
  // cout << "----- Completed Calculating Location Entropy --- " << endl;


  // SimpleQueries* query = new SimpleQueries(gpos, spos);

  // cout << "----- Precomputing matrices --- " << endl;
  // query->buildMatrices(0.5);
  // cout << "----- Completed Precomputing matrices --- " << endl;

  // cout << "----- Calculating Social Strength --- " << endl;
  // query->cacluateSocialStrength();
  // cout << "----- Completed calculating Social Strength --- " << endl;


  // cout << "----- Computing accuracy for threshold --- " << 0.5 <<endl;
  // query->verifySocialStrength(0.5);
  // cout << "--------------------------------------------";

  // cout << "----- Computing accuracy for threshold --- " << 0.6 <<endl;
  // query->verifySocialStrength(0.6);
  // cout << "--------------------------------------------";

  // cout << "----- Computing accuracy for threshold --- " << 0.7 <<endl;
  // query->verifySocialStrength(0.7);
  // cout << "--------------------------------------------";

  // cout << "----- Computing accuracy for threshold --- " << 0.8 <<endl;
  // query->verifySocialStrength(0.8);
  // cout << "--------------------------------------------";

  // cout << "----- Computing accuracy for threshold --- " << 0.9 <<endl;
  // query->verifySocialStrength(0.9);
  // cout << "--------------------------------------------";

  //testing grid snapping
  GPOs* gridSnappiedGPOs = new GPOs();
  gridSnappiedGPOs->createNewGPOsbyGridSnapping(gpos, r1); //second var is the x distance of a cell is km
  cout << "------------- Noise added -------------------" << endl;
  count_cooccurences(spos, gridSnappiedGPOs, 0, tR, isOptimistic);

  // if(r1 == 0){
  //   count_cooccurences(spos, gpos, r2, tR, isOptimistic);
  // } else {
  //   GPOs* purturbedGPOs = new GPOs();
  //   purturbedGPOs->loadPurturbedLocations(gpos, r1);
  //   cout << "------------- Noise added -------------------" << endl;

  //   count_cooccurences(spos, purturbedGPOs, r2, tR, isOptimistic);
  // }
}
