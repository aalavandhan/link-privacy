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
#endif


void count_cooccrrences(SPOs *spos, GPOs *gpos){
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  // Counting co-occurrences

  // for l in locations
    // for (u1,u2) checkied into l:
    // if are friends of [u1/u2]

  int tp = 0, p = 0, nFriendships=spos->edges;
  double precision, recall;
  float users_per_location=0;

  unordered_set< pair<int,int>, PairHasher > seen_pairs;

  for(auto loc = gpos->location_to_user.begin(); loc != gpos->location_to_user.end(); loc++){
    users_per_location += loc->second->size();
    for(auto it_u1 = loc->second->begin(); it_u1 != loc->second->end(); it_u1++){
      for(auto it_u2 = it_u1; it_u2 != loc->second->end(); it_u2++){
        if(  spos->areFriends(*it_u1, *it_u2) && seen_pairs.find(make_pair(*it_u1, *it_u2)) == seen_pairs.end() ){
          seen_pairs.insert(make_pair(*it_u1,*it_u2));
          tp++;
        }
        p++;
      }
    }
  }

  users_per_location /= gpos->location_to_user.size();
  cout << "Number of locations : " << gpos->location_to_user.size() << endl;
  cout << "Number of users : " << gpos->locationHistory.size() << endl;
  cout << "Mean of users at location : " << users_per_location << endl;
  cout << "Number of (unidirectional) friendships that can be inferred from co-occurrences " << tp << endl;
  cout << "Number of co-occurrences " << p << endl;
  cout << "Number of friendships " << nFriendships << endl;
  precision = tp / (double) p;
  recall    = tp / (double) nFriendships;
  cout << "Precision : " << precision << endl;
  cout << "Recall : " << recall << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

void countPurturbedCoOccurences(SPOs *spos, GPOs* gpos, double radius, bool isOptimistic){
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  double radius_geo_dist = (radius/1000) * 360 / EARTH_CIRCUMFERENCE,
         x=0, y=0;
  unsigned int lid, uid, count=0;

  unsigned int tp = 0, pos = 0, nFriendships=spos->edges;
  double precision, recall;

  unordered_set< pair<int,int>, PairHasher > seen_pairs;

  set<int>* seenLocations = new set<int>();

  for(auto l = gpos->locations.begin(); l != gpos->locations.end(); l++){
    x   = l->second->getX();
    y   = l->second->getY();
    lid = l->second->getID();
    uid = l->second->getUID();

    vector<res_point*>* checkins = gpos->getRange(x, y, radius_geo_dist);
    for(auto c = checkins->begin(); c != checkins->end(); c++){
      bool areFriends = spos->areFriends(uid, (*c)->uid),
         notRecoreded = seen_pairs.find(make_pair(uid, (*c)->uid)) == seen_pairs.end(),
         unseenLocation = (seenLocations->find((*c)->id) == seenLocations->end()),
         validLocation = isOptimistic || unseenLocation;

      if(  areFriends && notRecoreded &&  validLocation ){
        seen_pairs.insert(make_pair(uid, (*c)->uid));
        seenLocations->insert( (*c)->id );
        tp++;
      }

      if(validLocation){
        pos++;
      }

      delete (*c);
    }
    delete checkins;

    count++;

    if(count % 100000==0)
      cout << count << " " << seen_pairs.size() << endl;
  }

  cout << "Number of locations : " << gpos->location_to_user.size() << endl;
  cout << "Number of users : " << gpos->locationHistory.size() << endl;

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

  int i;
  double length, entropyTarget;
  double firstEntropy, secondEntropy, thirdEntropy, targetEntropy;
  // double firstMItarget, secondMItarget, thirdMItarget, targetMItarget;
  int *testFirstVector, *testSecondVector, *testThirdVector, *testMergedVector;
  struct timeval start,end;

  int *firstVector = (int *) calloc(4,sizeof(int));
  int *secondVector = (int *) calloc(4,sizeof(int));
  int *thirdVector = (int *) calloc(4,sizeof(int));
  int *targetVector = (int *) calloc(4,sizeof(int));
  
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
  
  printf("Entropies - first: %f, second: %f, third: %f, target %f\n",firstEntropy,secondEntropy,thirdEntropy,targetEntropy);
    
  // firstMItarget = calcMutualInformation(firstVector,targetVector,4);
  // secondMItarget = calcMutualInformation(secondVector,targetVector,4);
  // thirdMItarget = calcMutualInformation(thirdVector,targetVector,4);
  // targetMItarget = calcMutualInformation(targetVector,targetVector,4);
  
  // printf("MIs - first: %f, second: %f, third: %f, target %f\n",firstMItarget,secondMItarget,thirdMItarget,targetMItarget);
  
  testFirstVector = (int *) calloc(10000,sizeof(int));
  testSecondVector = (int *) calloc(10000,sizeof(int));
  testThirdVector = (int *) calloc(10000,sizeof(int));
  testMergedVector = (int *) calloc(10000,sizeof(int));
  
  for (i = 0; i < 10000; i++)
  {
    testFirstVector[i] = i % 2;
    testSecondVector[i] = i % 4;
    testThirdVector[i] = i % 3;
  }
  /* struct timeval
   * {
   *    time_t         tv_sec      seconds
   *    suseconds_t    tv_usec     microseconds
   * }
   */

  gettimeofday(&start, NULL);
  for (i = 0; i < 1000; i++)
  {
    // miTarget = calcMutualInformation(testFirstVector,testSecondVector,10000);
    entropyTarget = calcEntropy(testFirstVector,10000);
    // cmiTarget = calcConditionalMutualInformation(testFirstVector,testSecondVector,testThirdVector,10000);
    mergeArrays(testFirstVector,testSecondVector,testMergedVector,10000);
  }
  gettimeofday(&end, NULL);
  printf("H(X) = %f\n",entropyTarget);
  
  length = end.tv_sec - start.tv_sec;
  length = length + (end.tv_usec - start.tv_usec) / 1000000.0;
  
  printf("Time taken for a thousand I(X;Y), H(X), I(X;Y|Z), merge(X,Y) is %lf seconds\n",length);

}


int main(int argc, char *argv[]){
  cout.precision(15);

  DELTA_X = ((MAX_X - MIN_X)/ (X-1));
  DELTA_Y = ((MAX_Y - MIN_Y)/ (Y-1));

  double r1, r2, tp;
  r1 = atof(argv[3]);
  r2 = atof(argv[4]);
  tp = atof(argv[5]);
  bool isOptimistic = (tp == 0);

  testpTools();
  exit(-1);

  // Loading social network and checkins
  SPOs* spos = new SPOs();
  spos->load(argv[1]);

  GPOs* gpos = new GPOs(argv[2]);
  cout << "------------- Loading complete ---------------" << endl;


  if(tp == -1){
    count_cooccrrences(spos, gpos);
  } else {
    GPOs* purturbedGPOs = new GPOs();
    purturbedGPOs->loadPurturbedLocations(gpos, r1);
    cout << "------------- Noise added -------------------" << endl;
    countPurturbedCoOccurences(spos, purturbedGPOs, r2, isOptimistic);
  }
}
