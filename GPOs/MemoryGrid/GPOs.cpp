#include "../../headersMemory.h"

GPOs::GPOs (char* gridFileName){
  kNNExecutions = 0;
  LocationExecutions = 0;
  NextNNExecutions = 0;
  RangeExecutions = 0;
	pureNNexec = 0;
  totalCPUTime = totalTime = 0.0;
  grid = new Grid;
	loadLocations(gridFileName);
  generateFrequencyCache();
  // grid->deleteEmptyCells();
  objects = 0;
  computedNN = returnedNN = finalNextNN = 0;
  nextNNList = new vector<res_point*>();
  flagNextNN = true;
}

GPOs::GPOs(){
  kNNExecutions = 0;
  LocationExecutions = 0;
  NextNNExecutions = 0;
  RangeExecutions = 0;
  pureNNexec = 0;
  totalCPUTime = totalTime = 0.0;
  grid = new Grid;
  objects = 0;
  computedNN = returnedNN = finalNextNN = 0;
  nextNNList = new vector<res_point*>();
  flagNextNN = true;
}

GPOs::GPOs(GPOs *_gpos){
  kNNExecutions = 0;
  LocationExecutions = 0;
  NextNNExecutions = 0;
  RangeExecutions = 0;
  pureNNexec = 0;
  totalCPUTime = totalTime = 0.0;
  grid = new Grid;
  objects = 0;
  computedNN = returnedNN = finalNextNN = 0;

  flagNextNN = true;

  for(auto l = _gpos->locations.begin(); l != _gpos->locations.end(); l++){
    Point *p = *l;
    this->loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
  }
}

GPOs::~GPOs(){
  delete grid;


  for(auto it = user_to_location.begin(); it!= user_to_location.end(); it++){
    delete it->second;
  }

  for(auto it = location_to_user.begin(); it!= location_to_user.end(); it++){
    delete it->second;
  }

  for(auto it = time_to_checkins.begin(); it!= time_to_checkins.end(); it++){
    delete it->second;
  }

  // map< int, map<int,int>* >* location_to_user_to_cooccurrences;
  if(cooccurrences_created){
    for(auto it = location_to_user_to_cooccurrences->begin(); it!= location_to_user_to_cooccurrences->end(); it++){
      delete it->second;
    }
    delete location_to_user_to_cooccurrences;
  }
    //map< int, map<int, pair<int,double> >* >* user_to_order_to_location_displacment;
    if(user_to_order_to_location_displacment != NULL) {
      for(auto it = user_to_order_to_location_displacment->begin(); it!= user_to_order_to_location_displacment->end(); it++){
        delete it->second;
      }
      delete user_to_order_to_location_displacment;
    }

    //map<int, map<int, vector<uint>* >*> locations_users_frequency_map;
    for(auto it = locations_users_frequency_map.begin(); it!= locations_users_frequency_map.end(); it++){
      auto innerMap = it->second;
      for(auto iter = innerMap->begin(); iter!= innerMap->end(); iter++){
        delete iter->second;
      }
      delete innerMap;
    }

    //map<int, map<int, vector< pair<uint, int> >* >*> locations_users_frequency_map_with_order;
    for(auto it = locations_users_frequency_map_with_order.begin(); it!= locations_users_frequency_map_with_order.end(); it++){
      auto innerMap = it->second;
      for(auto iter = innerMap->begin(); iter!= innerMap->end(); iter++){
        delete iter->second;
      }
      delete innerMap;
    }

    // map<int, map<int, vector<pair<int, int> >* >*> cooccurrence_matrix;

    for(auto it = cooccurrence_matrix.begin(); it!= cooccurrence_matrix.end(); it++){
      auto innerMap = it->second;
      for(auto iter = innerMap->begin(); iter!= innerMap->end(); iter++){
        delete iter->second;
      }
      delete innerMap;
    }
}


double GPOs::getTotalCPUTime(){
  return totalCPUTime;
}

double GPOs::getTotalTime(){
  return totalTime;
}

map<int , vector< Point* >*>* GPOs::getLocationToUser(){
  return &location_to_user;
}

vector<Point*>* GPOs::getLocations(){
  return &locations;
}

void GPOs::generateFrequencyCache(){
  cout << "---- GENERATING CACHE ----" << endl;
  boost::posix_time::ptime time_t_epoch(boost::gregorian::date(2000 ,1,1));

  for(auto l_it = location_to_user.begin(); l_it != location_to_user.end(); l_it++){
    map<int, vector<uint>* > * user_frequencies = new map<int, vector<uint>* >();

    int lid = l_it->first;
    auto location_checkins_vector = l_it->second;

    for(auto c_it = location_checkins_vector->begin(); c_it!= location_checkins_vector->end(); c_it++){
      auto point = (*c_it);
      int uid = point->getUID();

      boost::posix_time::time_duration time_difference = point->getTime() - time_t_epoch;
      int time = abs(time_difference.total_seconds());

      auto user_it = user_frequencies->find(uid);
      if(user_it == user_frequencies->end()){
        vector<uint> *time_stamps = new vector<uint>();
        time_stamps->push_back( time );
        user_frequencies->insert(make_pair(uid, time_stamps));
      } else {
        user_it->second->push_back(time);
      }
    }

    for(auto it = user_frequencies->begin(); it!=user_frequencies->end();it++){
      vector<uint> *time_stamps = it->second;
      sort(time_stamps->begin(), time_stamps->end());
    }

    locations_users_frequency_map.insert(make_pair(lid, user_frequencies));
  }
}

void GPOs::generateCooccurrenceCache(){
  cout << "---- GENERATING CACHE ----" << endl;
  boost::posix_time::ptime time_t_epoch(boost::gregorian::date(2000 ,1,1));

  for(auto l_it = location_to_user.begin(); l_it != location_to_user.end(); l_it++){
    map<int, vector< pair<uint, int> >* > * user_frequencies = new map<int, vector< pair<uint, int> >* >();

    int lid = l_it->first;
    auto location_checkins_vector = l_it->second;

    for(auto c_it = location_checkins_vector->begin(); c_it!= location_checkins_vector->end(); c_it++){
      auto point = (*c_it);
      int uid = point->getUID();
      auto user_it = user_frequencies->find(uid);
      boost::posix_time::time_duration time_difference = point->getTime() - time_t_epoch;
      int time = abs(time_difference.total_seconds());

      if(user_it == user_frequencies->end()){
        vector<pair<uint, int>> *time_stamps = new vector<pair<uint, int>>();
        time_stamps->push_back( make_pair( time, point->getOrder() ) );
        user_frequencies->insert(make_pair(uid, time_stamps));
      } else {
        user_it->second->push_back( make_pair( time, point->getOrder() ) );
      }
    }

    for(auto it = user_frequencies->begin(); it!=user_frequencies->end();it++){
      vector<pair<uint, int>> *time_stamps = it->second;
      sort(time_stamps->begin(), time_stamps->end());
    }

    locations_users_frequency_map_with_order.insert(make_pair(lid, user_frequencies));
  }

  cout << "Cache Co-occurrence Generated " << locations_users_frequency_map_with_order.size() << endl;
}

vector< Point* >* GPOs::getLocations(int user_id){
  auto l = user_to_location.find(user_id);
  return l->second;
}

res_point* GPOs::getNextNN(double x, double y, int incrStep){
    clock_t startC, endC;
    struct timeval start, end;
    gettimeofday(&start, NULL);
    startC = clock();


    if(computedNN <= returnedNN && flagNextNN){
        NextNNExecutions++;
        computedNN+=incrStep;
        vector<res_point*>* kNN = grid->getkNN(x, y, computedNN);
        int size = kNN->size();
        for(int i = returnedNN; i < size; i++){
            nextNNList->push_back(util.copy((*kNN)[i]));
        }

        while(!kNN->empty()) {
            delete kNN->back();
            kNN->pop_back();
        }
        delete kNN;

        int newNNsize = nextNNList->size();
        if(computedNN > newNNsize){ // no more!
            //cout<<"here"<<endl;
            flagNextNN = false;
            computedNN = newNNsize;
        }
    }


    if(computedNN > returnedNN){
        endC = clock();
        totalCPUTime += (((double)(endC-startC)*1000.0)/(CLOCKS_PER_SEC));
        gettimeofday(&end, NULL);
        totalTime += util.print_time(start, end);
        return (*nextNNList)[returnedNN++];
    }
    else{
        endC = clock();
        totalCPUTime += (((double)(endC-startC)*1000.0)/(CLOCKS_PER_SEC));
        gettimeofday(&end, NULL);
        totalTime += util.print_time(start, end);

        // cout << "KNN Time : "<< util.print_time(start, end) << " ms" << endl;
        return NULL;
    }
}

void GPOs::getRangeSpatioTemporalBound(Point *p, vector< Point* >* results){
  double radius_geo_dist = (SPATIAL_HARD_BOUND/1000) * 360 / EARTH_CIRCUMFERENCE;
  vector <res_point*> *spatial_candidates = getRange(p->getX(), p->getY(), radius_geo_dist);

  for(auto sc_it=spatial_candidates->begin(); sc_it != spatial_candidates->end(); sc_it++){
    Point *cand = new Point(*sc_it);
    if(p->getTimeDifference(cand) <= TEMPORAL_HARD_BOUND * 3600){
      results->push_back(cand);
    } else
      delete cand;
  }
}

void GPOs::getSkylinePoints(Point *p, unordered_set< Point* > *skylines){
  vector <Point*> points_of_interest;
  getRangeSpatioTemporalBound(p, &points_of_interest);

  for(auto it=points_of_interest.begin(); it != points_of_interest.end(); it++){
    Point *chk = *it;

    if(chk->getID() == p->getID() && chk->getUID() != p->getUID())
      continue;

    if(skylines->size() == 0){
      skylines->insert( chk );
      continue;
    }

    bool chkDominates = false;
    for(auto sk_it = skylines->begin(); sk_it != skylines->end(); ){
      Point *skyline = (*sk_it);
      if( !p->doesSkylineDominate(skyline, chk) ){
        chkDominates = true;
        auto pt_to_delete = sk_it;
        sk_it++;
        skylines->erase(pt_to_delete);
      } else {
        sk_it++;
      }
    }

    if(chkDominates){
      skylines->insert( chk );
    }
  }
}

// Returns knn distance in kilometers
double GPOs::getKNNDistance(Point *p, int k){
  Point *neighbor = getKNN(p, k);
  double distance =  p->computeMinDistInKiloMeters(neighbor->getX(), neighbor->getY());
  delete neighbor;

  if(distance * 1000 <= SPATIAL_HARD_BOUND){
    return distance;
  }

  return std::numeric_limits<double>::infinity();
}

// Returns knn distance in seconds
double GPOs::getTemporalKNNDistance(Point *p, int k){
  map< int, Point* > temporalKNNs;

  getTemporalKNN(p, k, &temporalKNNs);

  if(temporalKNNs.size() > 0){
    auto it = temporalKNNs.rbegin();
    return (double) it->first;
  }

  return std::numeric_limits<double>::infinity();
}

// Returns knn temporal points
void GPOs::getTemporalKNN(Point *p, int k, map< int, Point* > *temporalKNNs){
  map< int, Point* > temporalKNNShortList;

  Point *neighbor = getKNN(p, k);
  int time_block = p->getTimeIndex();

  int left_size = 0;

  for(int i=0; i<=TEMPORAL_HARD_BOUND; i++){
    int cur_time_block = time_block + i;
    auto t_it = time_to_checkins.find(cur_time_block);

    if(t_it != time_to_checkins.end()){
      vector <Point*> *checkins = t_it->second;
      for(auto c_it = checkins->begin(); c_it != checkins->end(); c_it++){
        Point *chk = *c_it;
        if( p->getID() != chk->getID() ){
          temporalKNNShortList.insert( make_pair(p->getTimeDifference(chk), chk) );
        }
      }
    }

    if(temporalKNNShortList.size() >= k)
      break;
  }

  left_size = temporalKNNShortList.size();

  for(int i=1; i<=TEMPORAL_HARD_BOUND; i++){
    int cur_time_block = time_block - i;
    auto t_it = time_to_checkins.find(cur_time_block);

    if(t_it != time_to_checkins.end()){
      vector <Point*> *checkins = t_it->second;
      for(auto c_it = checkins->begin(); c_it != checkins->end(); c_it++){
        Point *chk = *c_it;
        if( p->getID() != chk->getID() ){
          temporalKNNShortList.insert( make_pair(p->getTimeDifference(chk), chk) );
        }
      }
    }

    if(temporalKNNShortList.size() >= k + left_size)
      break;
  }

  for(auto c_it = temporalKNNShortList.begin(); c_it != temporalKNNShortList.end() && temporalKNNs->size() < k; c_it++){
    temporalKNNs->insert(make_pair(c_it->first, c_it->second));
  }

}

// Get's the next nearest checkins discounting the current location
Point* GPOs::getKNN(Point *p, int k){
  int count = 0, incr=0;
  Point* neighbor;

  // cout << "Computing KNN : " << k  << " For " << p->getX() << " " << p->getY() << " " << p->getID() << " " << p->getUID() << endl;

  while(-1){
    res_point* next_neighbor = getNextNN(p->getX(), p->getY(), 25);
    incr++;

    // cout << incr << " " << next_neighbor->x << " " << next_neighbor->y << " " << next_neighbor->id << " " << next_neighbor->uid << endl;

    if( next_neighbor->id != p->getID() && next_neighbor->uid != p->getUID() ){
      // cout << "Found KNN " << endl;
      count++;
      if(count == k){
        neighbor = new Point(next_neighbor);
        clearNextNN();
        break;
      }
    }
  }
  // cout << "Nearest neighbor to " << p->getX() << " " << p->getY() << " is " << neighbor->getX() << " " << neighbor->getY() << endl;
  // cout << "Distance : " << neighbor->computeMinDistInMeters(p->getX(), p->getY()) << endl;
  // vector<Point *> *checkins = location_to_user.find(p->getID())->second;

  return neighbor;
};

vector<res_point*>* GPOs::getkNN(double x, double y, int k){
  // cout << "getKNN start" << endl;
  clock_t startC, endC;
  struct timeval start, end;
  gettimeofday(&start, NULL);
  startC = clock();
  kNNExecutions++;
  // cout << "----" << k << endl;
  vector<res_point*>* res = grid->getkNN(x, y, k);
  // cout << "size = " << res->size() << endl;
  endC = clock();
  totalCPUTime += (((double)(endC-startC)*1000.0)/(CLOCKS_PER_SEC));
  gettimeofday(&end, NULL);
  totalTime += util.print_time(start, end);
  //	cout << "getKNN end" << endl;
  return res;
}

vector<res_point*>* GPOs::getRangeAndDelete(double x, double y, double radius){
  vector<res_point*>* res = grid->getRangeAndDelete(x, y, radius);
  return res;
}

vector<res_point*>* GPOs::getRange(double x, double y, double radius){
  clock_t startC, endC;
  struct timeval start, end;
  gettimeofday(&start, NULL);
  startC = clock();

  RangeExecutions++;

  vector<res_point*>* res = grid->getRange(x, y, radius);

  endC = clock();
  totalCPUTime += (((double)(endC-startC)*1000.0)/(CLOCKS_PER_SEC));
  gettimeofday(&end, NULL);
  totalTime += util.print_time(start, end);

  //    cout << "Num of users (GPOs) in radius = "<<radius<<" are =" << res->size() << endl;
  return res;
}

set<res_point*, res_point_ascending_id>* GPOs::getSetRange(double x, double y, double radius){
  clock_t startC, endC;
  struct timeval start, end;
  gettimeofday(&start, NULL);
  startC = clock();

  RangeExecutions++;

  set<res_point*, res_point_ascending_id>* res = grid->getSetRange(x, y, radius);

  endC = clock();
  totalCPUTime += (((double)(endC-startC)*1000.0)/(CLOCKS_PER_SEC));
  gettimeofday(&end, NULL);
  totalTime += util.print_time(start, end);

  //    cout << "Num of users (GPOs) in radius = "<<radius<<" are =" << res->size() << endl;
  return res;
}


vector<res_point*>* GPOs::getRangeSortedId(double x, double y, double radius){
  clock_t startC, endC;
  struct timeval start, end;
  gettimeofday(&start, NULL);
  startC = clock();

  RangeExecutions++;

  vector<res_point*>* res = grid->getRange(x, y, radius);
  util.sortResPoint_AscendingId(res);

  endC = clock();
  totalCPUTime += (((double)(endC-startC)*1000.0)/(CLOCKS_PER_SEC));
  gettimeofday(&end, NULL);
  totalTime += util.print_time(start, end);
  return res;
}


int GPOs::getkNNExecutions(){
  return kNNExecutions;
}

int GPOs::getLocationExecutions(){
  return LocationExecutions;
}

int GPOs::getNextNNExecutions(){
  return NextNNExecutions;
}

int GPOs::getRangeExecutions(){
  return RangeExecutions;
}

// Only counting user pairs who've co-occurred more than once
unordered_set< pair<int,int>, PairHasher >* GPOs::getCoOccurredUserPairs(){
  return &significantly_cooccured_user_pairs;
}

unordered_map<int, double>* GPOs::getLocationEntropy(){
  return &location_to_H;
}

map<int, map<int, vector<pair<int, int> >* >*>* GPOs::getCooccurrenceMatrix(){
  return &cooccurrence_matrix;
}

map< int, map<int,int>* >* GPOs::getL2U2COOCC(){
  return location_to_user_to_cooccurrences;
}

// Load Point
void GPOs::loadPoint(double x, double y, int lid, int uid, boost::posix_time::ptime time, int order){
  Point* l;

  l = new Point(x, y, lid, uid, time, order);

  locations.push_back(l);
  auto lh_it = user_to_location.find(uid);

  if( lh_it == user_to_location.end() ){
    vector<Point* >* pts = new vector<Point*>();
    pts->push_back(l);
    user_to_location.insert( make_pair(uid, pts) );
  } else {
    vector<Point* >* pts = lh_it->second;
    pts->push_back(l);
  }

  auto lu_it = location_to_user.find(lid);
  if( lu_it == location_to_user.end() ){
    vector<Point* >* pts = new vector<Point*>();
    pts->push_back(l);
    location_to_user.insert( make_pair(lid, pts) );
  } else {
    vector<Point* >* pts = lu_it->second;
    pts->push_back(l);
  }

  int time_block = l->getTimeIndex();
  auto tu_it = time_to_checkins.find(time_block);
  if( tu_it == time_to_checkins.end() ){
    vector<Point* >* pts = new vector<Point*>();
    pts->push_back(l);
    time_to_checkins.insert( make_pair(time_block, pts) );
  } else {
    vector<Point* >* pts = tu_it->second;
    pts->push_back(l);
  }

  grid->addCheckIn(l);

  ids.push_back(lid);
};

// Load locations from file
bool GPOs::loadLocations(const char* fileName){
  cout << "Loading checkins from " << fileName << endl;
  ifstream fin(fileName);

  set<int> *blacklist_set;

  if(DATA_SET == 0 || DATA_SET == 1){
    int location_blacklist[] = {
          672270,    9246,  672328,  672473, 6936696,  671698, 7039899,
          671606, 6886742, 6886818, 6886783,    9247, 6745159, 6886788,
          681993, 6890222,  672035,  702472,  671276, 6887071, 7072522,
            9260,  672401,  672499,    9313,    9326,  671640, 6890228,
         7067634,  671886,    9256,    9227,    9222,  671896,  671812,
          671273,    9221, 7048639, 6886900, 7078149,    9410, 6887608,
          672585,  672070, 7044089, 6886902, 6936591, 6887558,  671665,
          672161,  671223, 6887599, 6887352,  672530, 6386788, 6887242,
         6886743, 6936933,  671570,  671165,   43073,  671728,    9225,
          420315, 6936505,   50782,  671917, 7050076,  671862,    9337,
         6887395,  718353, 6936452,  672406,  672110, 6886745,  688009,
          671967, 6936454,  672336,  712991,  671960,  672523, 6887163,
         5755950, 6886913, 6886773,  672142,  672659,  180328, 4555135,
         6886855, 6886766,   97009,   10305, 6887582, 7030027,   91192,
          693942,  671546,    9618, 6890227, 7198013, 1468206, 7002503,
          700158, 1016127, 6943397,  672400, 7288879,  672104, 6886790,
           10526,  672109,  694184,  689421, 6887100, 4992438,  672193,
         1221889, 4256132,   13075,    9249, 6886757, 6829868, 6886885,
         6887320,    9251, 6864294, 6887435,    9262, 6887532, 6887471,
          671988, 7013789, 6371718,  706040, 6887350,   41994,  672282,
         6887021,    9340,  712502, 6887001,  671781,  672073,    9371,
         6936673, 6937841, 6887473,  671360,  671727,  672017,   28301,
         7319924,  715389, 6936524,  671959,  771475, 4019380,  672200,
         6887112, 6829869,  671231,  694560, 6710875, 7363414, 6890254,
         6886899,   25634,  671485, 6347630, 6886801, 6886772,  671701,
         6887480, 6938064,  692344, 6936458,   17831,   33485,  694223,
         6887065,  582787,  672604, 6936699,  672645, 6886973, 6886760,
         6938075, 7483251, 7100872, 6362928, 6937337,  671358, 6936934,
          671348, 6887333,   14698,  671851, 6887204,   26125, 6890259,
         6886798, 7292568, 6887159,    9510, 7374702, 7380845, 6887063,
         5795051, 6887610,  671421,   43418, 6886595,  696167, 6374375,
         6830983,   17249, 6887127, 7537119, 1594494, 1136304,    9254,
          702615,  672046, 6490365,  671459, 6890223,    9474, 7138792,
         7227399,  672115,  671542,  671197, 6362954, 6829871, 7192395,
         6886938,  682460,  672185,  316637,  671672, 6396141,  671729,
          427712, 6887368,  671322, 7039857};

    blacklist_set = new set<int>(location_blacklist, location_blacklist+249);
  } else {
    int location_blacklist[] = {
      8361,  4079, 13794, 10270, 25065, 20548, 14486, 12174, 11918,
      8298, 33867, 16173, 11279, 14319,  9451,  8442, 25379,  9952,
      8294, 25256, 12120, 11363, 30519, 25608, 17462,  8742, 24384,
      4837, 14034, 10281,  3496,  3682, 27579,  9303, 10024, 19991,
      8338, 17861, 18955,  8137, 19982, 19930, 12561,   488, 11075
    };
    blacklist_set = new set<int>(location_blacklist, location_blacklist+35);
  }

  if (!fin) {
    std::cerr << "Cannot open checkins file " << fileName << std::endl;
    return false;
  }

  int uid, lid, count = 0;
  double x, y;

  std::string date, time, dateTime;

  while (fin){
    fin >> uid >> y >> x >> lid >> date >> time;
    dateTime = date + " " + time;
    boost::posix_time::ptime dtm = boost::posix_time::time_from_string(dateTime);

    // skip newlines
    if (!fin.good()) continue;

    if(blacklist_set->find(lid) == blacklist_set->end()){
      loadPoint(x, y, lid, uid, dtm, count);
      count ++;
    }

    if(count%100000==0)
      cout << count << endl;

  }

  fin.close();


  delete blacklist_set;


  cout << "------- LOADED GPOS from FILE --------- " << endl;
  cout << "Done! Number of checkins: " <<  count << endl;
  cout << "Number of failed checkins: " <<  grid->num_failed << endl;
  cout << "Number of locations: " <<  locations.size() << endl;
  cout << " -------------------------------------- " << endl;
  return true;
}


unordered_map<int, double>* GPOs::calculateLocationEntropy(){
  for(auto it = location_to_user.begin(); it != location_to_user.end(); it++){
    int location_id = it->first;
    vector< Point* >* checkins_at_location = it->second;


    map<int, int> location_counts;
    for(auto u_it = checkins_at_location->begin(); u_it!=checkins_at_location->end(); u_it++){
      int uid = (*u_it)->getUID();

      auto l_count = location_counts.find(uid);
      if(l_count == location_counts.end())
        location_counts.insert(make_pair(uid, 1));
      else
        l_count->second = l_count->second + 1;
    }

    //populating freqVector array with frequencies of users at this location
    uint *freqVector = (uint *) calloc(location_counts.size(), sizeof(uint));
    int i = 0;

    for(auto u_it = location_counts.begin(); u_it!=location_counts.end(); u_it++){
      freqVector[i] = u_it->second;
      i++;
    }

    //converts the frequencies to probabilities before computing the total shannon entropy
    double entropy = calcEntropyFromLocationVector(freqVector , location_counts.size());
    location_to_H.insert(make_pair(location_id, entropy));
  }

  // map <int,int> entropy_histogram;

  // for(auto l_it=location_to_H.begin(); l_it != location_to_H.end(); l_it++){
  //   double entropy = l_it->second;
  //   int bin = (int) floor(entropy * 10);
  //   auto it = entropy_histogram.find(bin);
  //   if(it != entropy_histogram.end()){
  //     it->second = it->second + 1;
  //   }
  //   else{
  //     entropy_histogram.insert(make_pair(bin,1));
  //   }
  // }

  // cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  // cout << "Entropy distribution  : " << endl;
  // for(auto b_it = entropy_histogram.begin(); b_it != entropy_histogram.end(); b_it++){
  //   cout << "Bin : " << b_it->first << "\t Count : " << b_it->second << endl;
  // }
  // cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  return &location_to_H;
}

// Range Queries
map<int, res_point*>* GPOs::getPointsInRange(double x, double y, double radius){
  double radius_geo_dist = (radius/1000) * 360 / EARTH_CIRCUMFERENCE;
  vector<res_point*> *points = getRange(x, y, radius_geo_dist);
  map<int, res_point*> *result = new map<int, res_point*>();
  for(auto p_it = points->begin(); p_it != points->end(); p_it++){
    res_point *p = *p_it;
    result->insert(make_pair(p->oid, p));
  }
  return result;
}

// r1 -> Outer radius, r2 -> inner radius
map<int, res_point*>* GPOs::getPointsInRange(double x, double y, double r1, double r2){
  map<int, res_point*> *r1_points = getPointsInRange(x,y,r1);
  map<int, res_point*> *r2_points = getPointsInRange(x,y,r2);

  map<int, res_point*> *result = new map<int, res_point*>();

  vector<int> r1_orders;
  vector<int> r2_orders;

  for(auto r_it = r1_points->begin(); r_it != r1_points->end(); r_it++)
    r1_orders.push_back( r_it->first );

  for(auto r_it = r2_points->begin(); r_it != r2_points->end(); r_it++)
    r2_orders.push_back( r_it->first );

  sort(r1_orders.begin(), r1_orders.end());
  sort(r2_orders.begin(), r2_orders.end());

  std::vector<int> order_intersection;

  std::set_intersection(r1_orders.begin(), r1_orders.end(),
                        r2_orders.begin(), r2_orders.end(),
                          std::back_inserter(order_intersection));

  for(auto o_it = order_intersection.begin(); o_it != order_intersection.end(); o_it++){
    int oid = (*o_it);

    auto p_it = r1_points->find(oid);
    res_point *p = p_it->second;

    result->insert(make_pair(p->oid, p));
  }

  r1_points->clear();
  delete r1_points;

  r2_points->clear();
  delete r2_points;

  return result;
}

unordered_map< int, vector<int>* >* GPOs::getUsersInRangeByHourBlock(double x, double y, double r1, double r2){
  unordered_map< int, vector<int>* >* r1_list = getUsersInRangeByHourBlock(x, y, r1);
  unordered_map< int, vector<int>* >* r2_list = getUsersInRangeByHourBlock(x, y, r2);

  unordered_map< int, vector<int>* >* user_list = new unordered_map< int, vector<int>* >();

  for(int i=0; i<7;i++){
    vector<int> *u1_list = r1_list->find(i)->second;
    vector<int> *u2_list = r2_list->find(i)->second;
    vector<int> result_list, *result;
    std::set_difference(u1_list->begin(), u1_list->end(), u2_list->begin(), u2_list->end(), std::inserter(result_list, result_list.begin()));
    result = new vector<int>(result_list.begin(), result_list.end());
    user_list->insert(make_pair(i, result));
  }

  return user_list;
}

unordered_map< int, vector<int>* >* GPOs::getUsersInRangeByHourBlock(double x, double y, double r){
  double radius_geo_dist = (r/1000) * 360 / EARTH_CIRCUMFERENCE;
  vector<res_point*>* checkins = getRange(x, y, radius_geo_dist);

  unordered_map< int, vector<int>* >* user_list = new unordered_map< int, vector<int>* >();

  for(int i=0; i<7;i++){
    vector<int> *users = new vector<int>();

    for(auto c = checkins->begin(); c != checkins->end(); c++){
      Point p = Point(*c);
      // Most checked in hour on any day
      if(p.getCheckinHour() == 23 && p.getCheckinDay() == i)
        users->push_back( p.getUID() );
    }

    // Sorting user list
    sort(users->begin(), users->end());
    // Removing duplicates
    users->erase( unique( users->begin(), users->end() ), users->end() );

    user_list->insert(make_pair(i, users));
  }

  for(auto c = checkins->begin(); c != checkins->end(); c++){
    delete (*c);
  }
  delete checkins;

  return user_list;
}

// r1 -> Outer radius, r2 -> inner radius
vector<int>* GPOs::getUsersInRange(double x, double y, double r1, double r2){
  vector<int> *u1_list = getUsersInRange(x, y, r1);
  vector<int> *u2_list = getUsersInRange(x, y, r2);
  vector<int> result_list, *result;
  std::set_difference(u1_list->begin(), u1_list->end(), u2_list->begin(), u2_list->end(), std::inserter(result_list, result_list.begin()));
  result = new vector<int>(result_list.begin(), result_list.end());
  return result;
}

vector<int>* GPOs::getUsersInRange(double x, double y, double radius){
  vector<int> *users = new vector<int>();

  double radius_geo_dist = (radius/1000) * 360 / EARTH_CIRCUMFERENCE;
  vector<res_point*>* checkins = getRange(x, y, radius_geo_dist);

  for(auto c = checkins->begin(); c != checkins->end(); c++){
    users->push_back( (*c)->uid );
    delete (*c);
  }
  delete checkins;

  // Sorting user list
  sort(users->begin(), users->end());
  // Removing duplicates
  users->erase( unique( users->begin(), users->end() ), users->end() );

  return users;
}

vector<int>* GPOs::getUsersInRangeAndTimeBlock(double x, double y, double time_block, int max_checkins, double max_radius){
  vector<int> *users = new vector<int>();

  double max_radius_geo_dist = (max_radius/1000) * 360 / EARTH_CIRCUMFERENCE;

  double radius_bound = estimateNearestDistance( x, y, max_checkins, max_radius_geo_dist);
  vector<res_point*>* checkins = getRange( x, y, radius_bound );

  for(auto c = checkins->begin(); c != checkins->end(); c++){
    res_point *rp = (*c);
    boost::posix_time::ptime p_time = rp->time;
    int p_time_block = (int)( (double)( p_time.time_of_day().hours() * 60 + p_time.time_of_day().minutes() ) / time_block );
    if(p_time_block == time_block)
      users->push_back( rp->uid );

    delete (*c);
  }
  delete checkins;

  // Sorting user list
  sort(users->begin(), users->end());
  // Removing duplicates
  users->erase( unique( users->begin(), users->end() ), users->end() );

  return users;
}

vector<int>* GPOs::getUsersInRange(int source, double radius){
  vector<int> *users = new vector<int>();
  auto source_users_checkins_it = user_to_location.find(source);
  if(source_users_checkins_it != user_to_location.end()){
    vector< Point* >* source_users_checkins = source_users_checkins_it->second;

    for(auto c_it=source_users_checkins->begin(); c_it != source_users_checkins->end(); c_it++){
      Point* checkin = (*c_it);
      vector<int> *n_users = getUsersInRange(checkin->getX(),  checkin->getY(), radius);
      users->insert( users->end(), n_users->begin(), n_users->end() );
      n_users->clear();
      delete n_users;
    }
    // Sorting user list
    sort(users->begin(), users->end());
    // Removing duplicates
    users->erase( unique( users->begin(), users->end() ), users->end() );
  }
  return users;
}


// Spatial Grouping
void GPOs::groupLocationsByKNNDistance(GPOs* gpos, int k, double std_radio){
  double x=0, y=0;
  unsigned int lid, count=0;
  long unsigned int iterations=0;
  unordered_set<int>* seenLocations = new unordered_set<int>();
  boost::posix_time::ptime time;

  GPOs *_duplicate_gpos = new GPOs(gpos);

  // Shuffling the locations
  std::random_shuffle( gpos->locations.begin(), gpos->locations.end() );

  for(auto l = gpos->locations.begin(); l != gpos->locations.end(); l++){
    Point *p = *l;
    x   = p->getX();
    y   = p->getY();
    lid = p->getID();

    // noise distance in meters
    double neighbor_distance = gpos->getKNNDistance(p, k) * 1000;

    vector<res_point*>* checkins = _duplicate_gpos->getRangeAndDelete(x, y, neighbor_distance * std_radio);

    for(auto c = checkins->begin(); c != checkins->end(); c++){
      if( seenLocations->find( (*c)->oid ) == seenLocations->end() ){
        loadPoint(x, y, lid, (*c)->uid, (*c)->time, (*c)->oid);
        seenLocations->insert( (*c)->oid );
        count++;
      }
      iterations++;
      delete (*c);
    }
    delete checkins;

    if(count % 100000==0)
      cout << count << " " << endl;

  };

  delete seenLocations;
  delete _duplicate_gpos;

  generateFrequencyCache();
};

void GPOs::groupLocationsByRange(GPOs* gpos, double radius, bool isOptimistic){
  double radius_geo_dist = (radius/1000) * 360 / EARTH_CIRCUMFERENCE,x=0, y=0;
  unsigned int lid, count=0;
  long unsigned int iterations=0;
  unordered_set<int>* seenLocations = new unordered_set<int>();
  boost::posix_time::ptime time;

  GPOs *_duplicate_gpos = new GPOs(gpos);

  // Shuffling the locations
  std::random_shuffle( gpos->locations.begin(), gpos->locations.end() );

  for(auto l = gpos->locations.begin(); l != gpos->locations.end(); l++){
    Point *p = *l;
    x   = p->getX();
    y   = p->getY();
    lid = p->getID();

    vector<res_point*>* checkins = _duplicate_gpos->getRangeAndDelete(x, y, radius_geo_dist);

    for(auto c = checkins->begin(); c != checkins->end(); c++){
      if(isOptimistic){
        loadPoint(x, y, lid, (*c)->uid, (*c)->time, (*c)->oid);
        count++;
      } else {
        if( seenLocations->find( (*c)->oid ) == seenLocations->end() ){
          loadPoint(x, y, lid, (*c)->uid, (*c)->time, (*c)->oid);
          seenLocations->insert( (*c)->oid );
          count++;
        }
        iterations++;
      }
      delete (*c);
    }
    delete checkins;

    if(count % 100000==0)
      cout << count << " " << endl;

  };
  cout << "Check-ins inserted : " << count << endl;
  cout << "Check-ins failed lookup : " << iterations << endl;

  delete seenLocations;
  delete _duplicate_gpos;

  generateFrequencyCache();
}

// Perturbation

// Time deviation in seconds
void GPOs::loadPurturbedLocationsByTime(GPOs* gpos, uint time_deviation){
  unsigned int lid = LOCATION_NOISE_BOUND, point_count =0;

  purturbed_count = 0;
  total_time_displacement=0;

  for(auto u = gpos->user_to_location.begin(); u != gpos->user_to_location.end(); u++){
    for(auto loc = u->second->begin(); loc != u->second->end(); loc++){
      if(time_deviation != 0){
        Point *p = (*loc);
        boost::posix_time::ptime purtubed_time = util.addTemporalGaussianNoise(p->getTime(), time_deviation);
        loadPoint( p->getX(), p->getY(), lid, u->first, purtubed_time, p->getOrder() );
        double displacement = abs((int)  ((p->getTime() - purtubed_time).total_milliseconds() / 1000.0));
        total_time_displacement+=displacement;
        purturbed_count++;
        temporal_purturbed_count++;
        lid++;
      } else {
        Point *p = (*loc);
        loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
      }
      point_count++;
    }
  }

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"spatially_purtubed_checkins{{"<< spatial_purturbed_count   << "}}" << endl;
  cout<<"temporally_purtubed_checkins{{"<< temporal_purturbed_count << "}}" << endl;

  cout<<"total_temporal_displacement{{"<< total_time_displacement/3600.0 <<"}} hours"<<endl;
  cout<<"average_temporal_displacement{{"<< total_time_displacement  * (1/(float)point_count) <<"}} seconds"<<endl;
  cout<<"average_temporal_displacement_on_purtubed{{"<< total_time_displacement/3600.0  * (1/(float)temporal_purturbed_count) <<"}} hours"<<endl;
}

// Radius in meters Adding Gaussian noise
void GPOs::loadPurturbedLocations(GPOs* gpos, double radius){
  unsigned int point_count = 0, lid=LOCATION_NOISE_BOUND;

  purturbed_count = 0;
  total_spatial_displacement = 0;

  for(auto u = gpos->user_to_location.begin(); u != gpos->user_to_location.end(); u++){
    for(auto loc = u->second->begin(); loc != u->second->end(); loc++){

      if(radius != 0){
        Point *p = (*loc);
        pair<double,double> coordinates_with_noise = util.addGaussianNoise(p->getX(), p->getY(), radius);
        double displacement = util.computeMinimumDistance(p->getX(), p->getY(), coordinates_with_noise.first, coordinates_with_noise.second);
        total_spatial_displacement+=displacement;
        purturbed_count++;
        spatial_purturbed_count++;
        loadPoint( coordinates_with_noise.first, coordinates_with_noise.second, lid, u->first, p->getTime(), p->getOrder() );
        lid++;
      } else {
        Point *p = (*loc);
        loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );

      }
      point_count++;
    }
  }

  total_spatial_displacement = total_spatial_displacement * EARTH_CIRCUMFERENCE / 360;

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"spatially_purtubed_checkins{{"<< spatial_purturbed_count   << "}}" << endl;
  cout<<"temporally_purtubed_checkins{{"<< temporal_purturbed_count << "}}" << endl;
  cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
  cout<<"average_spatial_displacement{{"<< (total_spatial_displacement / point_count) * 1000  <<"}} in meters"<<endl;
  cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / spatial_purturbed_count) * 1000 <<"}} in meters"<<endl;
}

void GPOs::computeSkylineMetrics(bool only_cooccurrences, map< int, map<int,int>* >* _location_to_user_to_cooccurrences){
  ofstream outfile;
  stringstream ss;
  std::string filePath;
  if(!only_cooccurrences){
    ss << "skyline.csv";
  } else {
    ss << "skyline-coocc.csv";
  }
  filePath = ss.str();
  outfile.open( filePath.c_str() );
  int checkin_count = 0;
  for(auto l_it = location_to_user.begin(); l_it != location_to_user.end(); l_it++){
    vector<Point *> *checkins = l_it->second;
    Point *first_point = checkins->at(0);

    if(only_cooccurrences){
      bool location_has_cooccurrences = _location_to_user_to_cooccurrences->find(first_point->getID()) != _location_to_user_to_cooccurrences->end();

      if(!location_has_cooccurrences)
        continue;
    }

    for(auto c_it = checkins->begin(); c_it != checkins->end(); c_it++){
      Point *checkin = (*c_it);

      unordered_set< Point* > skylines;
      getSkylinePoints(checkin, &skylines);

      for(auto sk_it=skylines.begin(); sk_it != skylines.end(); sk_it++){
        Point *skyline = (*sk_it);
        outfile << checkin->getOrder() << " " << skyline->getOrder() << " " << checkin->computeMinDistInKiloMeters(skyline->getX(), skyline->getY()) * 1000 << " " << checkin->getTimeDifference(skyline) << endl;
      }

      checkin_count++;

      if(checkin_count % 10000 == 0)
        cout << checkin_count << endl;

    }


  }
  outfile.close();
}

void GPOs::computeKNNDistances(int k, bool only_cooccurrences, bool compute_spatial, bool compute_temporal, map< int, map<int,int>* >* _location_to_user_to_cooccurrences){
  if(compute_spatial){
    ofstream outfile;
    stringstream ss;
    std::string filePath;
    if(!only_cooccurrences){
      ss << "knn-noise-" << k << ".csv";
    } else {
      ss << "knn-noise-" << k << "-coocc" << ".csv";
    }
    filePath = ss.str();
    outfile.open( filePath.c_str() );
    int location_count = 0;
    for(auto l_it = location_to_user.begin(); l_it != location_to_user.end(); l_it++){
      vector<Point *> *checkins = l_it->second;
      Point *first_point = checkins->at(0);

      if(only_cooccurrences){
        bool location_has_cooccurrences = _location_to_user_to_cooccurrences->find(first_point->getID()) != _location_to_user_to_cooccurrences->end();

        if(!location_has_cooccurrences)
          continue;
      }

      double distance = getKNNDistance(first_point, k);
      outfile << distance * 1000 << endl;

      location_count++;

      if(location_count % 100000 == 0)
        cout << location_count << endl;
    }
    outfile.close();
  }

  if(compute_temporal){
    ofstream outfile;
    stringstream ss;
    std::string filePath;
    if(!only_cooccurrences){
      ss << "knn-noise-temporal-" << k << ".csv";
    } else {
      ss << "knn-noise-temporal-" << k << "-coocc" << ".csv";
    }
    filePath = ss.str();
    outfile.open( filePath.c_str() );
    int location_count = 0;
    for(auto l_it = location_to_user.begin(); l_it != location_to_user.end(); l_it++){
      vector<Point *> *checkins = l_it->second;
      Point *first_point = checkins->at(0);
      double distance = getTemporalKNNDistance(first_point,  k);

      if(only_cooccurrences){
        bool location_has_cooccurrences = _location_to_user_to_cooccurrences->find(first_point->getID()) != _location_to_user_to_cooccurrences->end();

        if(!location_has_cooccurrences)
          continue;
      }

      outfile << distance * 1000 << endl;

      location_count++;

      if(location_count % 100000 == 0)
        cout << location_count << endl;
    }
    outfile.close();
  }
}

void GPOs::loadPurturbedLocationKNNDistance(GPOs* gpos, bool only_cooccurrences, int k, double std_radio, map< int, map<int,int>* >* _location_to_user_to_cooccurrences){
  unsigned int point_count = 0, lid=LOCATION_NOISE_BOUND;
  purturbed_count = 0;
  total_spatial_displacement = 0;

  for(auto l_it = gpos->location_to_user.begin(); l_it != gpos->location_to_user.end(); l_it++){
    vector<Point *> *checkins = l_it->second;
    Point *first_point = checkins->at(0);
    Point *neighbor = gpos->getKNN(first_point, k);

    for(auto loc = checkins->begin(); loc != checkins->end(); loc++){
      Point *p = (*loc);

      bool location_has_cooccurrences = _location_to_user_to_cooccurrences->find(first_point->getID()) != _location_to_user_to_cooccurrences->end();

      if(!only_cooccurrences || (only_cooccurrences && location_has_cooccurrences)){
        // converting noise_radius in meters
        double noise_radius = neighbor->computeMinDistInKiloMeters(p->getX(), p->getY()) * 1000;
        pair<double,double> coordinates_with_noise = util.addGaussianNoise(neighbor->getX(), neighbor->getY(), noise_radius * std_radio);
        double displacement = p->computeMinDistInKiloMeters(coordinates_with_noise.first, coordinates_with_noise.second);
        total_spatial_displacement+=displacement;
        purturbed_count++;
        spatial_purturbed_count++;
        loadPoint( coordinates_with_noise.first, coordinates_with_noise.second, lid, p->getUID(), p->getTime(), p->getOrder() );
        lid++;
      }

      else{
        loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
      }

      point_count++;

      if(point_count % 100000 == 0)
        cout << point_count << endl;
    }
    delete neighbor;
  }

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"spatially_purtubed_checkins{{"<< spatial_purturbed_count   << "}}" << endl;
  cout<<"temporally_purtubed_checkins{{"<< temporal_purturbed_count << "}}" << endl;
  cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
  cout<<"average_spatial_displacement{{"<< (total_spatial_displacement / point_count) * 1000  <<"}} in meters"<<endl;
  cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / spatial_purturbed_count) * 1000 <<"}} in meters"<<endl;
}


void GPOs::loadPurturbedLocationSelectiveKNNDistance(GPOs* gpos, int k, double std_radio, uint time_range_in_seconds, map< int, map<int,int>* >* _location_to_user_to_cooccurrences){
  gpos->generateCooccurrenceCache();

  ofstream outfile;
  outfile.open( "selective-knn-noise.csv" );

  unsigned int point_count = 0, lid=LOCATION_NOISE_BOUND;
  purturbed_count = 0;
  total_spatial_displacement = 0;

  for(auto l_it = gpos->location_to_user.begin(); l_it != gpos->location_to_user.end(); l_it++){
    vector<Point *> *checkins = l_it->second;
    Point *first_point = checkins->at(0);
    Point *neighbor = gpos->getKNN(first_point, k);

    auto loc_coocc_it = gpos->locations_users_frequency_map_with_order.find( first_point->getID() );

    bool location_has_cooccurrences = _location_to_user_to_cooccurrences->find(first_point->getID()) != _location_to_user_to_cooccurrences->end();
    location_has_cooccurrences = location_has_cooccurrences & (loc_coocc_it != locations_users_frequency_map_with_order.end());

    set<int> checkins_to_be_purturbed;

    if(location_has_cooccurrences){
      map<int, vector< pair<uint, int> >* >* loc_coocc = loc_coocc_it->second;
      for(auto u1_it = loc_coocc->begin(); u1_it != loc_coocc->end(); u1_it++){
        int user1 = u1_it->first;
        vector< pair<uint, int> >* u1_timestamps = u1_it->second;
        for(auto u2_it = loc_coocc->begin(); u2_it != loc_coocc->end(); u2_it++){
          int user2 = u2_it->first;
          vector< pair<uint, int> >* u2_timestamps = u2_it->second;

          if(user1 < user2){
            util.getCooccurrencesWithinTimeBlock(u1_timestamps, u2_timestamps, time_range_in_seconds, &checkins_to_be_purturbed);
          }
        }
      }
    }

    // Coorccurrences to be purturbed
    for(auto loc = checkins->begin(); loc != checkins->end(); loc++){
      Point *p = (*loc);

      bool need_to_purturb = checkins_to_be_purturbed.find(p->getOrder()) != checkins_to_be_purturbed.end();

      if(location_has_cooccurrences && need_to_purturb){
        double noise_radius = neighbor->computeMinDistInKiloMeters(p->getX(), p->getY()) * 1000;
        pair<double,double> coordinates_with_noise = util.addGaussianNoise(neighbor->getX(), neighbor->getY(), noise_radius * std_radio);
        double displacement = p->computeMinDistInKiloMeters(coordinates_with_noise.first, coordinates_with_noise.second);
        total_spatial_displacement+=displacement;
        purturbed_count++;
        spatial_purturbed_count++;
        loadPoint( coordinates_with_noise.first, coordinates_with_noise.second, lid, p->getUID(), p->getTime(), p->getOrder() );
        lid++;

        outfile << p->getX() << " " << p->getY() << " " << coordinates_with_noise.first << " " << coordinates_with_noise.second << " " << noise_radius << " " << displacement << endl;

      } else{
        loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
      }

      point_count++;

      if(point_count % 100000 == 0)
        cout << point_count << endl;

    }
  }

  outfile.close();

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"spatially_purtubed_checkins{{"<< spatial_purturbed_count   << "}}" << endl;
  cout<<"temporally_purtubed_checkins{{"<< temporal_purturbed_count << "}}" << endl;
  cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
  cout<<"average_spatial_displacement{{"<< (total_spatial_displacement / point_count) * 1000  <<"}} in meters"<<endl;
  cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / spatial_purturbed_count) * 1000 <<"}} in meters"<<endl;
}

void GPOs::loadPurturbedLocationsBasedOnCombinationFunction(
  GPOs* gpos,
  map< int, map<int, pair<int,double> >* >* user_to_order_to_location_locality,
  map< int, double >*  temporal_locality_map,
  map< int, map<int,int>* >* _location_to_user_to_cooccurrences,
  double radius,
  uint time_deviation,
  bool add_spatial,
  bool add_temporal,
  int noise_function){

  user_to_order_to_location_displacment = new map< int, map<int, pair<int,double> >* >();
  unordered_map<int, double>* HiL_map = gpos->getHiLasMap();
  unordered_map<int, double>* HiJ_map = gpos->getHiJasMap();

  uint lid        = LOCATION_NOISE_BOUND;
  int point_count = 0;


  if(add_spatial){
    if(noise_function == 0)
      cout << "Adding Spatial noise function : 0.5 * ( expHil + expHiJ ) * expHL * checkin_locality " << endl;
    else
      cout << "Adding Spatial noise function : log(1+Cij) * expHL * checkin_locality " << endl;
  }

  if(add_temporal)
    cout << "Adding Temporal noise function : temporal_locality * log(1+Cij) " << endl;


  for(auto u_it = gpos->user_to_location.begin(); u_it != gpos->user_to_location.end(); u_it++){

    int user_id = u_it->first;
    vector< Point* > *checkins = u_it->second;

    double hiL = 0, hiJ = 0;

    auto hil_it = HiL_map->find(user_id);
    if(hil_it != HiL_map->end())
      hiL = hil_it->second;

    auto hij_it = HiJ_map->find(user_id);
    if(hij_it != HiJ_map->end())
      hiL = hij_it->second;

    for(auto loc_it = checkins->begin(); loc_it != checkins->end(); loc_it++){
      Point *p = (*loc_it);
      int order = p->getOrder();

      auto entropy_it = gpos->location_to_H.find(p->getID());
      double entropy = 0;

      if(entropy_it != gpos->location_to_H.end())
        entropy = entropy_it->second;

      double user_cooccurrenes = 0;
      auto iter_outer = _location_to_user_to_cooccurrences->find(p->getID());
      if(iter_outer!= _location_to_user_to_cooccurrences->end()){
        map<int,int>* user_to_cooccurrences_map = iter_outer->second;
        auto iter_inner = user_to_cooccurrences_map->find(user_id);
        if(iter_inner != user_to_cooccurrences_map->end()){
          user_cooccurrenes = iter_inner->second;
        }
      }

      double checkin_locality_value = 0;
      auto iter = user_to_order_to_location_locality->find(user_id);
      if(iter!=user_to_order_to_location_locality->end()){
        map<int,pair<int,double > >* order_to_location_locality_map = iter->second;
        auto iter_inner = order_to_location_locality_map->find(order);
        if(iter_inner!= order_to_location_locality_map->end()){
          pair<int,double> p = iter_inner->second;
          checkin_locality_value = p.second;
        }
      }

      double temporal_locality_value = 0;
      auto tl_it = temporal_locality_map->find(order);
      if(tl_it != temporal_locality_map->end())
        temporal_locality_value = tl_it->second;

      double expHil=0,expHiJ=0;
      if(hiL == 0){
        expHil = 0;
      }else{
        expHil = exp(-hiL / HIL_SCALE);
      }

      if(hiJ == 0){
        expHiJ = 0;
      }else{
        expHiJ = exp(-hiJ / HIJ_SCALE);
      }

      double expHL = exp(-entropy / HL_SCALE);

      // Spatial Noise
      double spatial_noise;

      if(noise_function == 0){
        spatial_noise = 0.5 * ( expHil + expHiJ ) * expHL * checkin_locality_value;
        if(spatial_noise != 0)
          spatial_noise = spatial_noise + 0.45;
        spatial_noise = spatial_noise * radius;
      } else {
        spatial_noise = 0.5 * ( expHil + expHiJ ) * expHL * log( 1 + user_cooccurrenes );
        spatial_noise = spatial_noise * radius;
      }

      // Spatial perturbation
      pair<double,double> coordinates_with_noise = util.addGaussianNoise(p->getX(), p->getY(), spatial_noise);

      // Temporal noise
      double temporal_noise = (temporal_locality_value * log( 1 + user_cooccurrenes )) * (double) time_deviation;

      // Temporal perturbation
      boost::posix_time::ptime purtubed_time = util.addTemporalNoise(p->getTime(), (int) temporal_noise);

      // Loading perturbed data
      double final_x = p->getX(), final_y = p->getY();
      boost::posix_time::ptime final_time = p->getTime();
      int final_lid = p->getID();

      bool to_add_spatial_noise = (add_spatial && spatial_noise > 0);
      bool to_add_temporal_noise = (add_temporal && temporal_noise > 0);

      if(to_add_spatial_noise){
        final_x = coordinates_with_noise.first;
        final_y = coordinates_with_noise.second;
        spatial_purturbed_count++;
      }

      if(to_add_temporal_noise){
        final_time = purtubed_time;
        temporal_purturbed_count++;
      }

      if(to_add_spatial_noise || to_add_temporal_noise){
        final_lid = lid;
        lid++;
        purturbed_count++;
      }

      loadPoint(final_x, final_y, final_lid, p->getUID(), final_time, p->getOrder());
      point_count++;

      // METRICS
      total_spatial_displacement+=util.computeMinimumDistance(p->getX(), p->getY(), coordinates_with_noise.first, coordinates_with_noise.second);
      total_time_displacement+=abs((int)  ((p->getTime() - purtubed_time).total_milliseconds() / 1000.0));
    }
  }
  total_spatial_displacement = total_spatial_displacement * EARTH_CIRCUMFERENCE / 360;

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"spatially_purtubed_checkins{{"<< spatial_purturbed_count   << "}}" << endl;
  cout<<"temporally_purtubed_checkins{{"<< temporal_purturbed_count << "}}" << endl;

  if(add_spatial){
    cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
    cout<<"average_spatial_displacement{{"<< (total_spatial_displacement / point_count) * 1000  <<"}} in meters"<<endl;
    cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / spatial_purturbed_count) * 1000 <<"}} in meters"<<endl;
  }

  if(add_temporal){
    cout<<"total_temporal_displacement{{"<< total_time_displacement/3600.0 <<"}} hours"<<endl;
    cout<<"average_temporal_displacement{{"<< total_time_displacement  * (1/(float)point_count) <<"}} seconds"<<endl;
    cout<<"average_temporal_displacement_on_purtubed{{"<< total_time_displacement/3600.0  * (1/(float)temporal_purturbed_count) <<"}} hours"<<endl;
  }

  generateFrequencyCache();
}


// Compute CoOccurrences
void GPOs::countU2UCoOccurrences(uint time_block){
  int total_cooccurrences=0;

  cout<<"Number of locations: "<<locations_users_frequency_map.size()<<endl;

  for(auto l_it=locations_users_frequency_map.begin(); l_it != locations_users_frequency_map.end(); l_it++){
    int location_id = l_it->first;
    map<int, vector<uint>* > *user_checkin_times = l_it->second;

    if(user_checkin_times->size() == 1){
      continue;
    }


    for(auto u1_it=user_checkin_times->begin(); u1_it != user_checkin_times->end(); u1_it++){
      auto u2_it=u1_it;
      u2_it++;


      for(; u2_it != user_checkin_times->end(); u2_it++){
        int u1Id = u1_it->first;
        auto u1_timestamps = u1_it->second;
        int u2Id = u2_it->first;
        auto u2_timestamps = u2_it->second;

        if(u1Id > u2Id){
          int temp = u2Id;
          u2Id = u1Id;
          u1Id = temp;
        }

        int intersection_count = util.countIntersectionWithinTimeBlock(u1_timestamps,u2_timestamps,time_block, false);

        if(intersection_count > 0){
          auto coV_it = cooccurrence_matrix.find(u1Id);

          if(cooccured_user_pairs.find(make_pair(u1Id, u2Id)) == cooccured_user_pairs.end())
           cooccured_user_pairs.insert(make_pair( u1Id, u2Id ));

          total_cooccurrences += intersection_count;
          if(coV_it != cooccurrence_matrix.end()){ //u1 exist in matrix
            map<int, vector<pair<int, int> >*>* inner_map = coV_it->second;
            //fnd user 2 and update the location list with location_id -> intersection count
            auto flist_it = inner_map->find(u2Id);
            if(flist_it!=inner_map->end()){  //u2 exist in matrix
              vector<pair<int,int> >* location_freqeuncy_list = flist_it->second;
              location_freqeuncy_list->push_back(make_pair(location_id,intersection_count));
            }else{  //u2 does not exist in matrix
              vector<pair<int,int> >* location_freqeuncy_list = new vector<pair<int,int> >();
              location_freqeuncy_list->push_back(make_pair(location_id,intersection_count));
              inner_map->insert(make_pair(u2Id,location_freqeuncy_list));
            }
          }else{  //u1 does not exist in matrix
            map<int, vector<pair<int, int> >* >* inner_map = new map<int, vector<pair<int, int> >* >();
            vector<pair<int,int> >* location_freqeuncy_list = new vector<pair<int,int> >();
            location_freqeuncy_list->push_back(make_pair(location_id,intersection_count));
            inner_map->insert(make_pair(u2Id,location_freqeuncy_list));
            cooccurrence_matrix.insert(make_pair(u1Id,inner_map));
          }
        }

        // if(intersection_count > 1){
        //   cout<<"Location id: "<< location_id <<" uniques count of users checked in: "<< user_checkin_times->size()<<endl;
        //   cout << u1Id<<"  "<<u2Id <<" size: "<<u1_timestamps->size()<<" | "<<u2_timestamps->size()<<endl;
        //   cout<<"intersection_count = "<< intersection_count <<" Cooccurrence matrix size: "<< cooccurrence_matrix.size()<< endl;

        //   cout << "U1 Time stamps : ";
        //   for(auto t=u1_timestamps->begin(); t!= u1_timestamps->end(); t++){
        //     cout << (*t) << " ";
        //   }
        //   cout<<endl;
        //   cout << "U1 Time stamps : ";
        //   for(auto t=u2_timestamps->begin(); t!= u2_timestamps->end(); t++){
        //     cout << (*t) << " ";
        //   }
        //   cout<<endl;

        //   util.countIntersectionWithinTimeBlock(u1_timestamps,u2_timestamps,time_block, true);

        //   auto it = location_to_user.find(location_id);
        //   auto vector_of_points =  it->second;
        //   for(auto v_it = vector_of_points->begin(); v_it != vector_of_points->end(); v_it++){
        //     Point* p = *v_it;
        //     if(p->getUID() == u1Id || p->getUID() ==  u2Id)
        //     cout<<"\t UserID: "<<p->getUID()<<" Time: "<<p->getTime()<<endl;
        //   }

        // }

      }
    }
  }

  //instantiating map for use with combination function based on CiL
  cooccurrences_created = true;
  location_to_user_to_cooccurrences = new map< int, map<int,int>* >();

  for(auto c_it = cooccurrence_matrix.begin(); c_it != cooccurrence_matrix.end(); c_it++){
    int user_1 = c_it->first;
    auto users_location_frequency_map = c_it->second;

    for(auto ulh_it = users_location_frequency_map->begin(); ulh_it != users_location_frequency_map->end(); ulh_it++){
      int user_2 = ulh_it->first, cooccurrence_count=0;
      vector<pair<int, int>>* cooccurrence_counts_vector = ulh_it->second;

      for(auto l_it=cooccurrence_counts_vector->begin(); l_it != cooccurrence_counts_vector->end(); l_it++){
        int location_id = l_it->first;
        int cooccrences_at_l = l_it->second;
        cooccurrence_count += cooccrences_at_l;

        //for cil
        auto iter_outer = location_to_user_to_cooccurrences->find(location_id);
        if(iter_outer != location_to_user_to_cooccurrences->end()){
          //update/create for both users
          map<int,int>* user_to_cooccurrences_map = iter_outer->second;
          auto iter_inner = user_to_cooccurrences_map->find(user_1);
          if(iter_inner != user_to_cooccurrences_map->end()){
            iter_inner->second = iter_inner->second + cooccrences_at_l;
          }else{
            user_to_cooccurrences_map->insert(make_pair(user_1,cooccrences_at_l));
          }

          iter_inner = user_to_cooccurrences_map->find(user_2);
          if(iter_inner != user_to_cooccurrences_map->end()){
            iter_inner->second = iter_inner->second + cooccrences_at_l;
          }else{
            user_to_cooccurrences_map->insert(make_pair(user_2,cooccrences_at_l));
          }
        }else{
          map<int,int>* user_to_cooccurrences_map = new map<int,int>();
          user_to_cooccurrences_map->insert(make_pair(user_1,cooccrences_at_l));
          user_to_cooccurrences_map->insert(make_pair(user_2,cooccrences_at_l));
          location_to_user_to_cooccurrences->insert(make_pair(location_id, user_to_cooccurrences_map));
        }
      }
      //-------------------------

      if(cooccurrence_count>1){
        if(significantly_cooccured_user_pairs.find(make_pair(user_1, user_2)) == significantly_cooccured_user_pairs.end())
          significantly_cooccured_user_pairs.insert(make_pair( user_1, user_2 ));
      } else {
    // remove


        if(insignificantly_cooccured_user_pairs.find(make_pair(user_1, user_2)) == insignificantly_cooccured_user_pairs.end())
          insignificantly_cooccured_user_pairs.insert(make_pair( user_1, user_2 ));
      }
    }
  }


  //To calculate basic cooccurence based metric on noised data
  //  first calculate u1 - > u2 , u3, u4 , u5  cooccurrences.
  // for(auto c_it = cooccurrence_matrix.begin(); c_it != cooccurrence_matrix.end(); c_it++){
  //   int user_1 = c_it->first;
  //   auto users_location_frequency_map = c_it->second;
  //   map<int, vector<pair<int, int> >* >* users_location_frequency_map_insignificant = new map<int, vector<pair<int, int> >* >();

  //   for(auto ulh_it = users_location_frequency_map->begin(); ulh_it != users_location_frequency_map->end(); ulh_it++){
  //     int user_2 = ulh_it->first;
  //     vector<pair<int, int>>* cooccurrence_counts_vector = ulh_it->second;
  //     vector<pair<int, int>>* cooccurrence_counts_vector_insignificant =  new vector<pair<int, int>>();

  //     for(auto l_it=cooccurrence_counts_vector->begin(); l_it != cooccurrence_counts_vector->end(); l_it++){
  //       int location_id = l_it->first;
  //       int cooccrences_at_l = l_it->second;
  //       cooccurrence_counts_vector_insignificant->push_back(make_pair(location_id,cooccrences_at_l));
  //     }
  //     users_location_frequency_map_insignificant->insert(make_pair(user_2,cooccurrence_counts_vector_insignificant));
  //   }
  //   cooccurrence_matrix_insignificant.insert(make_pair(user_1, users_location_frequency_map_insignificant));
  // }

  for(auto p=insignificantly_cooccured_user_pairs.begin(); p!=insignificantly_cooccured_user_pairs.end(); p++){
    int u1=p->first, u2=p->second;

    auto u1_it = cooccurrence_matrix.find(u1);
    if(u1_it != cooccurrence_matrix.end()){
      auto users_location_frequency_map = u1_it->second;
      auto u2_it = users_location_frequency_map->find(u2);
      vector<pair<int, int>>* cooccurrence_counts_vector = u2_it->second;
      // delete cooccurrence_counts_vector;

      if(u2_it != users_location_frequency_map->end()){
        vector<pair<int, int>>* cooccurrence_counts_vector = u2_it->second;
        delete cooccurrence_counts_vector;
        users_location_frequency_map->erase(u2_it);
      }

      if(users_location_frequency_map->size() == 0){
        cooccurrence_matrix.erase(u1_it);
        delete users_location_frequency_map;
      }
    }
  }
  cout << "users_with_atleast_one_cooccurrence{{" << cooccurrence_matrix.size() << "}}" << endl;
  cout << "total_cooccurrence{{" << total_cooccurrences << "}}" << endl;
  cout << "cooccurred_user_paris{{" << cooccured_user_pairs.size() << "}}" << endl;
  cout << "significantly_cooccurred_user_paris{{" << significantly_cooccured_user_pairs.size() << "}}" << endl;
}

int GPOs::getUserCooccurrences(int user_id){
  int number_of_cooccurrences = 0;
  auto it = cooccurrence_matrix.find(user_id);
  if(it != cooccurrence_matrix.end()){
    auto map_of_vectors = it->second;
    for(auto it_map = map_of_vectors->begin(); it_map != map_of_vectors->end(); it_map++){
      auto vector_of_locations = it_map->second;
      for(auto it_vector = vector_of_locations->begin(); it_vector!= vector_of_locations->end();it_vector++){
        number_of_cooccurrences += it_vector->second;
      }
    }
  }
  return number_of_cooccurrences;
}


// HISTOGRAMS
void GPOs::printCooccurrenceMatrix(char *DATASET_PATH){

  ofstream outfile;

  stringstream ss;
  ss << DATASET_PATH << "cooccurrence_matrix.csv";
  const std::string filePath = ss.str();
  outfile.open( filePath.c_str() );

  //for each user user compute his total cooccurrences at all locations and users
  //then iterate again to produce probabilities
  //then compute shanon entropy
  for(auto it = cooccurrence_matrix.begin(); it!= cooccurrence_matrix.end(); it++){

    int user_id_1 = it->first;
    auto map_of_vectors = it->second;

    for(auto it_map = map_of_vectors->begin(); it_map != map_of_vectors->end(); it_map++){
      int user_id_2 = it_map->first;
      auto vector_of_locations = it_map->second;

      for(auto it_vector = vector_of_locations->begin(); it_vector!= vector_of_locations->end();it_vector++){
        // number_of_cooccurrences += it_vector->second;
        int location_id = it_vector->first;

        outfile << user_id_1 << " " << user_id_2<<" "<<location_id <<" "<< it_vector->second << endl;

      }
    }

  }
  outfile.close();
}

map< int, double >* GPOs::computeTemporalLocality(int max_checkins, double max_radius){
  int counter=0;
  map< int, double >* temporal_locality_map = new map< int, double >();
  double max_radius_geo_dist = (max_radius/1000) * 360 / EARTH_CIRCUMFERENCE;

  for(auto l_it=location_to_user.begin(); l_it != location_to_user.end(); l_it++){
    vector< Point* >* checkins_at_l = l_it->second;
    Point *p_sample = checkins_at_l->front();
    double x = p_sample->getX(), y = p_sample->getY();

    double radius_bound = estimateNearestDistance( x, y, max_checkins, max_radius_geo_dist);
    vector<res_point*>* checkins_in_city = getRange( x, y, radius_bound );

    for(auto c_it=checkins_at_l->begin(); c_it != checkins_at_l->end(); c_it++){
      Point *p = (*c_it);
      boost::posix_time::ptime l_time = p->getTime();

      double locality_sum=0, vicinity_count=0, temporal_locality=0;

      for(auto c = checkins_in_city->begin(); c != checkins_in_city->end(); c++){
        boost::posix_time::ptime c_time = (*c)->time;

        boost::posix_time::time_duration time_difference = (c_time - l_time);
        double diff = (double) abs(time_difference.total_seconds());

        if( diff <= 8 * 3600  && (*c)->oid != p->getOrder() ){
          locality_sum += exp( (-diff / (double) (8*3600) ) );
          vicinity_count++;
        }

      }

      if(vicinity_count > 0)
        temporal_locality = locality_sum/vicinity_count;

      temporal_locality_map->insert(make_pair(p->getOrder(), temporal_locality));
    }

    for(auto c = checkins_in_city->begin(); c != checkins_in_city->end(); c++){
      delete (*c);
    }
    delete checkins_in_city;

    counter++;
    if(counter%10000 == 0)
      cout << counter << endl;
  }

  return temporal_locality_map;
}

unordered_map<int, double>* GPOs::getHiLasMap(){
  unordered_map<int, double>* HiL_histogram = new unordered_map<int, double>();

  //for each user user compute his total cooccurrences at all locations and users
  //then iterate again to produce probabilities
  //then compute shanon entropy
  for(auto it = cooccurrence_matrix.begin(); it!= cooccurrence_matrix.end(); it++){
    map <int, int> location_counts;
    int user_id = it->first;
    auto map_of_vectors = it->second;
    for(auto it_map = map_of_vectors->begin(); it_map != map_of_vectors->end(); it_map++){
      auto vector_of_locations = it_map->second;
      for(auto it_vector = vector_of_locations->begin(); it_vector!= vector_of_locations->end();it_vector++){
        // number_of_cooccurrences += it_vector->second;
        int location_id = it_vector->first;
        auto lc_it = location_counts.find(location_id);
        if(lc_it != location_counts.end()){
          lc_it->second = lc_it->second + it_vector->second;
        }else{
          location_counts.insert(make_pair(location_id,it_vector->second));
        }

      }
    }

    uint *freqVector = (uint *) calloc(location_counts.size(), sizeof(uint));
    int i = 0;
    for(auto u_it = location_counts.begin(); u_it!=location_counts.end(); u_it++){
      freqVector[i] = u_it->second;
      i++;
    }
    //converts the frequencies to probabilities before computing the total shannon entropy
    double entropy = calcEntropyFromLocationVector(freqVector , location_counts.size());
    HiL_histogram->insert(make_pair(user_id, entropy));
  }

  return HiL_histogram;
}

unordered_map<int, double>* GPOs::getHiJasMap(){
  unordered_map<int, double>* HiJ_histogram = new unordered_map<int, double>();


  //for each user user compute his total cooccurrences at all locations and users
  //then iterate again to produce probabilities
  //then compute shanon entropy
  for(auto it = cooccurrence_matrix.begin(); it!= cooccurrence_matrix.end(); it++){
    map <int, int> user_counts;
    int user_id = it->first;
    auto map_of_vectors = it->second;
    for(auto it_map = map_of_vectors->begin(); it_map != map_of_vectors->end(); it_map++){
      int user_id_j = it_map->first;
      auto vector_of_locations = it_map->second;
      int cooccurrence_count_for_user_j = 0;
      for(auto it_vector = vector_of_locations->begin(); it_vector!= vector_of_locations->end();it_vector++){
        // number_of_cooccurrences += it_vector->second;
        cooccurrence_count_for_user_j += it_vector->second;
      }
      user_counts.insert(make_pair(user_id_j,cooccurrence_count_for_user_j));
    }

    uint *freqVector = (uint *) calloc(user_counts.size(), sizeof(uint));
    int i = 0; int sum_freq = 0;
    for(auto u_it = user_counts.begin(); u_it!=user_counts.end(); u_it++){
      freqVector[i] = u_it->second;
      sum_freq += u_it->second;
      i++;
    }

    //converts the frequencies to probabilities before computing the total shannon entropy
    double entropy = calcEntropyFromCoV(freqVector , user_counts.size());
    HiJ_histogram->insert(make_pair(user_id, entropy));
  }

  return HiJ_histogram;
}

unordered_map< int, map< int, map<int, double >* >* >* GPOs::getPltMap(double time_block, int max_checkins, double max_radius){
  unordered_map< int, map< int, map<int, double >* >* >* Plt_histogram = new unordered_map< int, map< int, map<int, double >* >* >();

  double max_radius_geo_dist = (max_radius/1000) * 360 / EARTH_CIRCUMFERENCE;
  int processing_count=0;

  boost::posix_time::ptime time_t_epoch(boost::gregorian::date(2000 ,1,1));
  cout << "Computing getPltMap : " << endl;

  for(auto l_it=location_to_user.begin(); l_it != location_to_user.end(); l_it++){
    processing_count++;

    if(processing_count%10000 == 0)
      cout << processing_count << endl;

    int location_id = l_it->first;
    vector< Point* >* checkins_at_l = l_it->second;
    Point *p_sample = checkins_at_l->front();
    double x = p_sample->getX(), y = p_sample->getY();

    double radius_bound = estimateNearestDistance( x, y, max_checkins, max_radius_geo_dist);
    vector<res_point*>* checkins_in_city = getRange( x, y, radius_bound );

    for(auto c_it=checkins_at_l->begin(); c_it != checkins_at_l->end(); c_it++){
      Point *p = (*c_it);
      boost::posix_time::ptime l_time = p->getTime();
      int l_time_block = p->getTimeBlock( time_block );
      int l_day        = l_time.date().day_of_year() + ( l_time.date().year() - 2000 ) * 365;

      double prob_lt = 0;
      int points_in_vicinity = 0;
      double h = 60.0;

      for(auto c = checkins_in_city->begin(); c != checkins_in_city->end(); c++){
        boost::posix_time::ptime c_time = (*c)->time;

        boost::posix_time::time_duration time_difference = (c_time - l_time);
        double diff = (double) time_difference.total_seconds();

        if( abs(diff) <= 3600*12 && (*c)->oid != p->getOrder() ){
          prob_lt = prob_lt + exp( - (diff * diff / h) / (double)( 2* (5*60)^2 ) );
          points_in_vicinity++;
        }
      }

      if(points_in_vicinity > 0)
        prob_lt = 1/((double)points_in_vicinity*h) * prob_lt;

      map< int, map<int, double >* > *day_block_prob_map;
      map<int, double> *block_prob_map;

      auto pl_it = Plt_histogram->find( location_id );
      if(pl_it == Plt_histogram->end()){
        day_block_prob_map = new map< int, map<int, double >* >();
        block_prob_map = new map<int, double>();

        block_prob_map->insert(make_pair(l_time_block, prob_lt));
        day_block_prob_map->insert(make_pair(l_day, block_prob_map));
        Plt_histogram->insert(make_pair(location_id, day_block_prob_map));
      } else{
        day_block_prob_map = pl_it->second;
        auto pldt_it = day_block_prob_map->find(l_day);
        if(pldt_it == day_block_prob_map->end()){
          block_prob_map = new map<int, double>();
        } else {
          block_prob_map = pldt_it->second;
        }
        auto plt_it = block_prob_map->find(l_time_block);
        if(plt_it == block_prob_map->end()){
          block_prob_map->insert(make_pair(l_time_block, prob_lt));
        } else {
          plt_it->second = plt_it->second + prob_lt;
        }
      }
    }

    for(auto c = checkins_in_city->begin(); c != checkins_in_city->end(); c++){
      delete (*c);
    }
    delete checkins_in_city;
  }

  return Plt_histogram;
}

unordered_map<int, double>* GPOs::getHlLasMap(){
  map<int, map<int,int>* > location_to_user_coocc_map;
  unordered_map<int, double>* HlL_histogram = new unordered_map<int, double>();


  ofstream outfile;
  outfile.open("user_to_cooccurrences.csv");

  //for each user user compute his total cooccurrences at all locations and users
  //then iterate again to produce probabilities
  //then compute shanon entropy
  for(auto it = cooccurrence_matrix.begin(); it!= cooccurrence_matrix.end(); it++){

    int user_id_1 = it->first;
    auto map_of_vectors = it->second;

    set<int> number_of_locations;

    for(auto it_map = map_of_vectors->begin(); it_map != map_of_vectors->end(); it_map++){
      int user_id_2 = it_map->first;
      auto vector_of_locations = it_map->second;

      for(auto it_vector = vector_of_locations->begin(); it_vector!= vector_of_locations->end();it_vector++){
        // number_of_cooccurrences += it_vector->second;
        int location_id = it_vector->first;
        number_of_locations.insert(location_id);

        auto lc_it = location_to_user_coocc_map.find(location_id);
        if(lc_it != location_to_user_coocc_map.end()){ //location found in outer map
          auto user_coocc_map = lc_it->second;

          //for user_id_1
          auto inner_user_it = user_coocc_map->find(user_id_1);
          if(inner_user_it!=user_coocc_map->end()){  //user in inner map
            inner_user_it->second = inner_user_it->second + it_vector->second; // update user 1's cocc
          }else{
            user_coocc_map->insert(make_pair(user_id_1,it_vector->second));

          }

          //for user_id_2
          inner_user_it = user_coocc_map->find(user_id_2);
          if(inner_user_it!=user_coocc_map->end()){  //user in inner map
            inner_user_it->second = inner_user_it->second + it_vector->second; // update user 1's cocc
          }else{
            user_coocc_map->insert(make_pair(user_id_2,it_vector->second));

          }

        }else{                                        //location not found in outer map
          map<int,int>* user_coocc_map = new map<int,int>();
          user_coocc_map->insert(make_pair(user_id_1,it_vector->second)); //add value for u1
          user_coocc_map->insert(make_pair(user_id_2,it_vector->second)); //add value for u2
          location_to_user_coocc_map.insert(make_pair(location_id,user_coocc_map)); //add inner map to outermap
        }

      }
    }

    outfile << user_id_1 << " " << number_of_locations.size() << endl;
  }
  outfile.close();

  for(auto it = location_to_user_coocc_map.begin(); it!= location_to_user_coocc_map.end(); it++){
    int location_id = it->first;
    auto user_coocc_map = it->second;
    uint *freqVector = (uint *) calloc(user_coocc_map->size(), sizeof(uint));
    int i = 0;
    for(auto iter = user_coocc_map->begin() ; iter != user_coocc_map->end() ; iter++){
      freqVector[i] = iter->second;
      i++;
    }
    //converts the frequencies to probabilities before computing the total shannon entropy
    double entropy = calcEntropyFromCoV(freqVector , user_coocc_map->size());
    HlL_histogram->insert(make_pair(location_id, entropy));
  }
  return HlL_histogram;
}


// UNUSED
void GPOs::loadLocationsByDayOfWeek(GPOs* gpos, int day_of_week){
  vector<Point*>* o_locations = gpos->getLocations();
  for(auto l_it=o_locations->begin(); l_it != o_locations->end(); l_it++){
    Point *p = (*l_it);
    boost::posix_time::ptime l_time = p->getTime();

    if(l_time.date().day_of_week() == day_of_week)
      loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
  }
  cout << "Loaded checkins with day of the week " << day_of_week << endl;
  cout << "Original data set : " << gpos->location_to_user.size() << endl;
  cout << "Reduced data set : " << location_to_user.size() << endl;
  cout << "-------------------------------------------------------" << endl;
}

//input:  grid cell distance in x direction (say, 100m grid). This value is independent of grid size in headers.h.
//output: all checkins are snapped to the nearest cell corner.
void GPOs::createNewGPOsbyGridSnapping(GPOs* gpos, double grid_distance_on_x_axis_in_km){
  //TODO: create a vistual grid from the input grid dimension = grid_distance_on_x_axis
  //assign location id to each grid corner = n^2 location ids
  double grid_distance_on_x_axis_in_geodist = grid_distance_on_x_axis_in_km * 360 / EARTH_CIRCUMFERENCE;
  double max_x = MAX_X - 0.00001; double max_y = MAX_Y - 0.00001;
  double min_x = MIN_X + 0.00001; double min_y = MIN_Y + 0.00001;


  // Changes the dimensions of the virtual grid to be a square
  // Note: may cause some checkins to fail in the direction in which the grid is extended.
  // this happens because the points may get assisgned to grid corners in the virtual grid
  // that are outside the boundary of the actual grid. (tested on gowalla and found no failed checkins)
  // if(MAX_X - MIN_X > MAX_Y - MIN_Y){  //adjust in y direction by moving the boundaries equally in the north and the south
  //   double seperation = (MAX_X - MIN_X) - (MAX_Y - MIN_Y);
  //   seperation = seperation/2;
  //   min_y =  MIN_Y - seperation;
  //   max_y =  MAX_Y + seperation;

  // }else if(MAX_X - MIN_X < MAX_Y - MIN_Y){  //vice versa
  //   double seperation = (MAX_Y - MIN_Y) - (MAX_X - MIN_X);
  //   seperation = seperation/2;
  //   min_x =  MIN_X - seperation;
  //   max_x =  MAX_X + seperation;
  // }
  //-----------------------------------------------------


  cout<<"Adjusted MIN_X: "<<min_x<<" MAX_X: "<<max_x<<endl;
  cout<<"Adjusted MIN_Y: "<<min_y<<" MAX_Y: "<<max_y<<endl;

  int grid_size = ceil( (max_x - min_x) / grid_distance_on_x_axis_in_geodist ) + 1;
  double delta_on_x = ((max_x - min_x)/ (grid_size));
  double delta_on_y = ((max_y - min_y)/ (grid_size));

  cout<<"Grid Size from input cell size: "<< grid_size<<" X "<<grid_size<<endl;
  cout<<"Adjusted DELTA_X: "<<delta_on_x<<" DELTA_Y : "<<delta_on_y<<endl;
  cout<<"Adjusted (in meters) DELTA_X: "<<delta_on_x*1000*EARTH_CIRCUMFERENCE/360<<" DELTA_Y : "<<delta_on_y*1000*EARTH_CIRCUMFERENCE/360<<endl;

  // int asd = 100;

  //for each checkin do:
  //  find which cell it belongs to
  //  get the cell corner which is closest
  //  load this point to this corner
  for(auto u = gpos->user_to_location.begin(); u != gpos->user_to_location.end(); u++){
    for(auto loc = u->second->begin(); loc != u->second->end(); loc++){

      // asd--;
      // if(asd==0){
      //   exit(-1);
      // }

      Point* p = *loc;
      int q_x = (int)((p->getX() - min_x)/delta_on_x);
      int q_y = (int)((p->getY() - min_y)/delta_on_y);


      // cout<<"Point cell: "<<q_x<<" , " <<q_y<<endl;
      // p->printDetails();

      if(q_x >=0 && q_x < grid_size && q_y >=0 && q_y < grid_size){
        double min_dist = 999;
        double closest_x;
        double closest_y;
        int closest_corner_i,closest_corner_j;
        for(int i = 0;i < 2;i++){
          for(int j = 0;j < 2;j++){
            double distance = util.computeMinimumDistance(p->getX(),p->getY(),min_x+((q_x+i)*delta_on_x), min_y+((q_y+j)*delta_on_y));
            // cout<<"Corner : "<<i*2+j<<" Coordinate: ("<<min_x+((q_x+i)*delta_on_x)<<" , "<<min_y+((q_y+j)*delta_on_y)<<")"<<endl;
            // cout<<"Distance: "<<distance*1000*EARTH_CIRCUMFERENCE/360<<" meters"<<endl;
            if(distance < min_dist){
              min_dist = distance;
              closest_x = min_x+((q_x+i)*delta_on_x);
              closest_y = min_y+((q_y+j)*delta_on_y);
              closest_corner_i = i;
              closest_corner_j = j;
            }
          }
        }
        int location_id = ((closest_corner_i + q_x)*grid_size)+( closest_corner_j + q_y);
        loadPoint(closest_x, closest_y, location_id, u->first, p->getTime(), p->getOrder());
      }
      else{
        cout<<"cell out of bounds...check grid size"<<endl;
      }
    }
  }

  generateFrequencyCache();
}

void GPOs::clearNextNN(){
  //  cout << "clearNextNN nextNNList size = " << nextNNList->size() << endl;
  while(!nextNNList->empty()) {
    delete nextNNList->back();
    nextNNList->pop_back();
    //    objects--;
  }
  //  cout << "clearNextNN nextNNList size = " << nextNNList->size() << " objects = " << objects << endl;
  delete nextNNList;
  //        nextNNList->clear();
  //        delete nextNNList;
  nextNNList = new vector<res_point*>();
  computedNN = returnedNN = finalNextNN = objects = 0;
  flagNextNN = true;
 pureNNexec = 0;
}

double GPOs::estimateNearestDistance(double x, double y, int k, double max_radius){
  return grid->estimateNearestDistance(x,y,k,max_radius);
}

double GPOs::distanceBetween(Point *a, Point *b){
  double lat1r, lon1r, lat2r, lon2r, u, v;
  lat1r = DEG_TO_RAD * a->getY();
  lon1r = DEG_TO_RAD * a->getX();
  lat2r = DEG_TO_RAD * b->getY();
  lon2r = DEG_TO_RAD * b->getX();
  u = sin((lat2r - lat1r)/2);
  v = sin((lon2r - lon1r)/2);
  return 2.0 * EARTH_RADIUS_IN_KILOMETERS * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}


