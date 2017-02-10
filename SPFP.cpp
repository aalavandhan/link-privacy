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

// double MIN_X = -180;
// double MAX_X = 180;
// double MIN_Y = -90;
// double MAX_Y = 90;

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

void testDistance(){
  GPOs *g = new GPOs();

  Point p1 = Point(-77.037852, 38.898556, 1);
  Point p2 = Point(-77.043934, 38.897147, 2);
  cout << "Distance between (-77.037852, 38.898556) and (-77.043934, 38.897147) is " << g->distanceBetween(&p1, &p2) << "km" << endl;
}


int main(int argc, char *argv[]){
  cout.precision(15);

  DELTA_X = ((MAX_X - MIN_X)/ (X-1));
  DELTA_Y = ((MAX_Y - MIN_Y)/ (Y-1));

  double r1, r2, tp, time_range_in_seconds;
  r1                    = atof(argv[3]);
  r2                    = atof(argv[4]);
  time_range_in_seconds = atof(argv[5]);
  tp                    = atof(argv[6]);
  bool isOptimistic     = (tp == 0);
  auto locations_of_interest_file = argv[7];

  ALPHA = 0.480;
  BETA = 0.520;

  // testDistance();
  // testpTools();
  // exit(-1);

  // cout << "------------- Loading checkins ---------------" << endl;
  // GPOs* gpos = new GPOs(argv[2]);
  // cout << "------------- Loading complete ---------------" << endl;


  cout << "------------- Loading checkins ---------------" << endl;
  GPOs* g = new GPOs(argv[2]);
  cout << "------------- Loading complete ---------------" << endl;

  // testing grid snapping
  // GPOs* gpos = new GPOs();
  // gpos->createNewGPOsbyGridSnapping(g, r1); //second var is the x distance of a cell is km
  // cout << "Number of locations loaded " << gpos->locations.size() << endl;
  // cout << "------------- Noise added -------------------" << endl;

  GPOs* gg = new GPOs();
  GPOs* gpos = new GPOs();
  gg->loadPurturbedLocations(g, r1);
  cout << "------------- Locations perturbed -------------------" << endl;
  gpos->groupLocationsByRange(gg, r2, isOptimistic);
  cout << "------------- Locations Grouped -------------------" << endl;


  // cout << "----- Loading Cooccurrence Matrix --- " << endl;
  // gpos->countU2UCoOccurrences((uint) time_range_in_seconds);
  // cout << "----- Completed Loading Cooccurrence Matrix --- " << endl;


  // cout << "----- Calculating Location Entropy --- " << endl;
  // gpos->calculateLocationEntropy();
  // cout << "----- Completed Calculating Location Entropy --- " << endl;

  // Loading social network and checkins
  cout << "------------- Loading SocialGraph ---------------" << endl;
  SPOs* spos = new SPOs();
  spos->load(argv[1]);

  // cout << "------------- Loading SocialGraph ---------------" << endl;
  // SPOs* spos = new SPOs(gpos);
  // spos->load(argv[1]);

  // cout << "------------- Computing mean distance between friend pairs ---------------" << endl;
  // cout << "Mean distance between all pairs of friends :" << spos->computeMeanDistanceBetweenAllFriends(gpos) << endl;

  // cout << "------------- Computing node locality of all users ---------------" << endl;
  // spos->computeNodeLocality(gpos);
  // spos->writeNodeLocalityToFile();

  // cout << "------------- Load computed node locality ---------------" << endl;
  // spos->loadNodeLocalityFromFile();

  SimpleQueries* query = new SimpleQueries(gpos, spos);

  cout << "------------- Evaluating range utility ---------------" << endl;
  query->checkUtilityRange(locations_of_interest_file, g, 25);
  query->checkUtilityRange(locations_of_interest_file, g, 50);
  query->checkUtilityRange(locations_of_interest_file, g, 100);
  query->checkUtilityRange(locations_of_interest_file, g, 200);
  query->checkUtilityRange(locations_of_interest_file, g, 400);


  // cout << "------------- Evaluating range proximity ---------------" << endl;
  // query->checkUtilityProximity(locations_of_interest_file, g, 100, 50);
  // query->checkUtilityProximity(locations_of_interest_file, g, 500, 100);
  // query->checkUtilityProximity(locations_of_interest_file, g, 500, 200);

  // cout << "------------- Evaluating query locations ---------------" << endl;
  // query->checkUtilityStats(locations_of_interest_file, 10);
  // query->checkUtilityStats(locations_of_interest_file, 25);
  // query->checkUtilityStats(locations_of_interest_file, 50);
  // query->checkUtilityStats(locations_of_interest_file, 100);
  // query->checkUtilityStats(locations_of_interest_file, 200);
  // query->checkUtilityStats(locations_of_interest_file, 400);
  // query->checkUtilityStats(locations_of_interest_file, 800);
  // query->checkUtilityStats(locations_of_interest_file, 1600);



  // cout << "----- Precomputing matrices --- " << endl;
  // query->buildMatrices(0.1);
  // cout << "----- Completed Precomputing matrices --- " << endl;

  // cout << "----- Calculating Social Strength --- " << endl;
  // query->cacluateSocialStrength();
  // cout << "----- Completed calculating Social Strength --- " << endl;


  // double thresholds[9] = {0.5,1,2,3,4,5,7,10,25};
  // for(double i = 0.8; i < 2.2; i = i + 0.1){
  //   cout << "----- Computing accuracy for threshold --- " << i <<endl;
  //   query->verifySocialStrength(i);
  //   cout << "--------------------------------------------";
  // }

  // cout << "----- Computing accuracy for threshold --- " << 2.1 <<endl;
  // map< int, bool >* users_of_interest = query->getUsersOfInterest(2.1);
  // cout << "--------------------------------------------";

  // map<int, Point*> locations;
  // ofstream output_file;
  // output_file.open("locations-of-interest.csv");

  // for (auto u = users_of_interest->begin(); u != users_of_interest->end(); u++){
  //   auto pts = gpos->getLocations(u->first);
  //   for(auto p=pts->begin(); p!=pts->end(); p++){
  //     double x=(*p)->getX(), y=(*p)->getY();
  //     locations.insert( make_pair((*p)->getID(), *p) );
  //     output_file << x << "," << y <<endl;
  //   }
  // }
  // cout << "------- Locations of interest : " << locations.size() << endl;
  // output_file.close();

  // if(r1 == 0){
  //   count_cooccurences(spos, gpos, r2, tR, isOptimistic);
  // } else {
  //   GPOs* purturbedGPOs = new GPOs();
  //   purturbedGPOs->loadPurturbedLocations(gpos, r1);
  //   cout << "------------- Noise added -------------------" << endl;

  //   count_cooccurences(spos, purturbedGPOs, r2, tR, isOptimistic);
  // }
}
