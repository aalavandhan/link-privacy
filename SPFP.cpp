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


// Pair hasher
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

//hashes a pair into an unordered set
struct PairHasher
{
    size_t operator()(const pair<int, int> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
};

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
  cout << "Number of users : " << gpos->locationHistory.size() << endl;

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

  // Loading social network and checkins
  SPOs* spos = new SPOs();
  spos->load(argv[1]);

  GPOs* gpos = new GPOs(argv[2]);
  cout << "------------- Loading complete ---------------" << endl;

  gpos->countU2UCoOccurrences();

  // int coocc_count = 0, tot;

  // for(int i=0; i<USER_COUNT; i++){
  //   for(int j=0; j<USER_COUNT; j++){
  //     auto c_list = gpos->coocc_matrix[i][j];
  //     if(c_list != NULL && i != j){
  //       tot=0;
  //       for(auto l=c_list->begin(); l != c_list->end(); l++ ){
  //         tot += l->second;
  //       }
  //       if(tot > 1){
  //         coocc_count++;
  //       }
  //     }
  //   }
  // }

  // cout << "Number of user pairs " << coocc_count << endl;

  // if(r1 == 0){
  //   count_cooccurences(spos, gpos, r2, tR, isOptimistic);
  // } else {
  //   GPOs* purturbedGPOs = new GPOs();
  //   purturbedGPOs->loadPurturbedLocations(gpos, r1);
  //   cout << "------------- Noise added -------------------" << endl;
  //   count_cooccurences(spos, purturbedGPOs, r2, tR, isOptimistic);
  // }
}
