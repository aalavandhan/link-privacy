#include "../../headersMemory.h"

template <typename It>
double Median(It begin, It end)
{
    using T = typename std::iterator_traits<It>::value_type;
    std::vector<T> data(begin, end);
    std::nth_element(data.begin(), data.begin() + data.size() / 2, data.end());
    return data[data.size() / 2];
}

GPOs::GPOs (char* gridFileName, uint time_range, double spatial_range){
  kNNExecutions = 0;
  LocationExecutions = 0;
  NextNNExecutions = 0;
  RangeExecutions = 0;
	pureNNexec = 0;
  totalCPUTime = totalTime = 0.0;
  grid = new Grid;
	loadLocations(gridFileName);

  coocc_time_range      = time_range;
  coocc_spatial_range   = spatial_range;

  // generateCooccurrenceCache();
  // grid->deleteEmptyCells();
  objects = 0;
  computedNN = returnedNN = finalNextNN = 0;
  nextNNList = new vector<res_point*>();
  flagNextNN = true;
}

GPOs::GPOs(uint time_range, double spatial_range){
  coocc_time_range      = time_range;
  coocc_spatial_range   = spatial_range;

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
  coocc_time_range    = _gpos->coocc_time_range;
  coocc_spatial_range = _gpos->coocc_spatial_range;

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

  for(auto l = _gpos->checkin_list.begin(); l != _gpos->checkin_list.end(); l++){
    Point *p = l->second;
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

    for(auto it = cooccurrence_index.begin(); it!= cooccurrence_index.end(); it++){
      auto innerVector = it->second;
      delete innerVector;
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

vector< pair<int,int> >* GPOs::getCooccurredCheckins(){
  vector< pair<int,int> >* cooccurred_checkins_vector = new vector< pair<int,int> >();
  for(auto c_it = cooccurred_checkins.begin(); c_it != cooccurred_checkins.end(); c_it++){
    int o1 = c_it->first;
    int o2 = c_it->second;
    cooccurred_checkins_vector->push_back(make_pair(o1,o2));
  }
  return cooccurred_checkins_vector;
}

void GPOs::generateCooccurrenceCache(){
  cout << "---- Clearing CACHE ----" << endl;
  for(auto it = locations_users_frequency_map_with_order.begin(); it!= locations_users_frequency_map_with_order.end(); it++){
    auto innerMap = it->second;
    for(auto iter = innerMap->begin(); iter!= innerMap->end(); iter++){
      delete iter->second;
    }
    delete innerMap;
  }

  cout << "---- GENERATING CACHE ----" << endl;

  cout << "del_s for co-occ :"  << coocc_spatial_range  << endl;
  cout << "del_t for co-occ :"  << coocc_time_range     << endl;

  double average_chkns_in_spatial_range = 0, max_chkns_in_spatial_range = 0;

  for(auto l_it = location_to_user.begin(); l_it != location_to_user.end(); l_it++){
    int lid = l_it->first;
    auto checkins_at_l = l_it->second;

    map<int, vector< pair<uint, int> >* > *user_frequencies = new map<int, vector< pair<uint, int> >* >();

    Point *sample = checkins_at_l->front();
    vector<int> *nearby_locations;

    if(coocc_spatial_range != 0){
      nearby_locations = getLocationsInRange(sample->getX(), sample->getY(), coocc_spatial_range);
    } else {
      nearby_locations = new vector<int>();
      nearby_locations->push_back(lid);
    }

    average_chkns_in_spatial_range += nearby_locations->size();
    if(nearby_locations->size() > max_chkns_in_spatial_range){
      max_chkns_in_spatial_range = nearby_locations->size();
    }

    for(auto nl_it = nearby_locations->begin(); nl_it != nearby_locations->end(); nl_it++){
      int n_location = (*nl_it);
      vector< Point* > *checkins_at_n_location = location_to_user.find(n_location)->second;

      for(auto c_it = checkins_at_n_location->begin(); c_it!= checkins_at_n_location->end(); c_it++){
        auto point = (*c_it);
        int uid = point->getUID();
        auto user_it = user_frequencies->find(uid);
        int time = point->getTimeInSeconds();

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
    }
    locations_users_frequency_map_with_order.insert(make_pair(lid, user_frequencies));
    delete nearby_locations;
  }

  average_chkns_in_spatial_range = average_chkns_in_spatial_range / locations_users_frequency_map_with_order.size();

  cout << "Cache Co-occurrence Generated " << locations_users_frequency_map_with_order.size() << endl;
  cout << "Average number of locations in vicinity :  " << average_chkns_in_spatial_range << endl;
  cout << "Max number of locations in vicinity     :  " << max_chkns_in_spatial_range << endl;
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

vector <res_point*>* GPOs::getRangeSpatioTemporalBound(Point *p){
  return getRangeSpatioTemporalBound(p, SPATIAL_HARD_BOUND, TEMPORAL_HARD_BOUND);
}

vector <res_point*>* GPOs::getRangeSpatioTemporalBound(Point *p, double spatial_bound_in_meters, double temporal_bound_in_hours){
  vector <Point*> temporal_candidates;
  vector <res_point*>* candidates = new vector <res_point*>();

  getRangeByTime(p->getTimeIndex() - (temporal_bound_in_hours + 1), p->getTimeIndex() + (temporal_bound_in_hours + 1), &temporal_candidates);

  for(auto cd_it = temporal_candidates.begin(); cd_it != temporal_candidates.end(); cd_it++){
    Point *q = (*cd_it);
    if(p->getUID() != q->getUID() && p->getID() != q->getID()){
      double dist = p->computeMinDistInKiloMeters(q->getX(), q->getY());
      double time_deviation = abs((p->getTime() - q->getTime()).total_seconds());
      if( (dist * 1000) <= spatial_bound_in_meters && ((double)time_deviation/3600.0) <= temporal_bound_in_hours){
        res_point *rp = new res_point();
        rp->id = q->getID();
        rp->uid = q->getUID();
        rp->oid = q->getOrder();
        rp->x = q->getX();
        rp->y = q->getY();
        rp->time = q->getTime();
        rp->dist = dist * 360/EARTH_CIRCUMFERENCE;
        candidates->push_back(rp);
      }
    }
  }

  return candidates;
}


void GPOs::getSkylinePoints(Point *p, vector <res_point*> *candidates, map< int, pair<int, res_point*> > *skylines){
  for(auto it=candidates->begin(); it != candidates->end(); it++){
    res_point *chk = *it;

    if(chk->id == p->getID())
      continue;

    if(chk->id == p->getUID())
      continue;

    if(skylines->size() == 0){
      skylines->insert(make_pair(chk->oid, make_pair(0, chk)));
      continue;
    }

    bool checkinIsDominated = false;
    for(auto sk_it = skylines->begin(); sk_it != skylines->end(); sk_it++){
      int d_count = sk_it->second.first;
      res_point *skyline = sk_it->second.second;
      if(p->doesSkylineDominatePoint(skyline, chk)){
        checkinIsDominated = true;
        sk_it->second.first = d_count + 1;
        break;
      }
    }

    if(!checkinIsDominated){
      int domination_count = 0;
      for(auto sk_it = skylines->begin(); sk_it != skylines->end(); ){
        int skyline_id = sk_it->first;
        int d_count = sk_it->second.first;
        res_point *skyline = sk_it->second.second;
        if(p->doesPointDominateSkyline(skyline, chk)){
          auto pt_to_delete = sk_it;
          sk_it++;
          skylines->erase(skyline_id);
          domination_count = domination_count + d_count;
        } else {
          sk_it++;
        }
      }

      skylines->insert(make_pair(chk->oid, make_pair(domination_count, chk)));
    }
  }
}

double GPOs::getSTKNNDistance(Point *p, int k, vector<res_point*> *candidates, int metric_type){
  priority_queue < pair<double, res_point*>, vector<pair<double, res_point*> > > spatioTemporalKNNs;

  getSpatioTemporalKNN(p, k, &spatioTemporalKNNs, candidates, metric_type);

  if(spatioTemporalKNNs.size() == k){
    return spatioTemporalKNNs.top().first;
  }

  return std::numeric_limits<double>::infinity();
}

// metric_type == 0/3 ( Spatial and Temporal )
// metric_type == 4 ( Spatial and Temporal, dont skip co-locations )
// metric_type == 1 ( Spatial   )
// metric_type == 2 ( Temporal  )
void GPOs::getSpatioTemporalKNN(Point *p, int k,
  priority_queue < pair<double, res_point*>, vector<pair<double, res_point*> > > *spatioTemporalKNNs,
  vector<res_point*> *candidates,
  int metric_type){

  // timeval start_metric,end_metric, start_bound, end_bound;
  // gettimeofday(&start_bound, NULL);

  // gettimeofday(&end_bound, NULL);
  // bound_computation_time+=util.print_time(start_bound, end_bound);
  // gettimeofday(&start_metric, NULL);

  bool coocc_at_p = cooccurrence_index.find(p->getOrder()) != cooccurrence_index.end();

  for(auto it=candidates->begin(); it != candidates->end(); it++){
    res_point *chk = *it;

    if(coocc_at_p && metric_type != 4){
      unordered_set<int>* cooccurred_checkins = cooccurrence_index_indirect.find(p->getOrder())->second;
      // Skip co-occurred check-ins
      if( cooccurred_checkins->find(chk->oid) != cooccurred_checkins->end() ){
        continue;
      }

      // Skip users of co-occurred check-ins
      for(auto coch_it = cooccurred_checkins->begin(); coch_it != cooccurred_checkins->end(); coch_it++){
        int coch_order = (*coch_it);
        Point *ch = checkin_list.find(coch_order)->second;
        if(ch->getUID() == chk->uid)
          continue;
      }
    }

    // Discount the current location/region and current user
    if( p->getID() != chk->id && p->getUID() != chk->uid){
      double distance = 0.0;

      if(metric_type == 0 || metric_type == 3 || metric_type == 4){
        distance = p->getSTDistance(chk, coocc_spatial_range, coocc_time_range);
        if(distance <= 1)
          continue;
      }
      else if(metric_type == 1){
        double spatial_distance  = (chk->dist * EARTH_CIRCUMFERENCE / 360.0) / ( (double)coocc_spatial_range / 1000.0 );
        distance = spatial_distance;
      }
      else if(metric_type == 2){
        double temporal_distance = (double) p->getTimeDifference(chk) / ( (double)coocc_time_range );
        distance = temporal_distance;
      }

      if(spatioTemporalKNNs->size() < k){
        spatioTemporalKNNs->push(make_pair(distance, chk));
      }else{
        if(distance < spatioTemporalKNNs->top().first){
          spatioTemporalKNNs->push(make_pair(distance, chk));
          spatioTemporalKNNs->pop();
        }
      }
    }
  }
  // gettimeofday(&end_metric, NULL);
  // metric_computation_time+=util.print_time(start_metric, end_metric);
}


// Get's the next nearest checkins discounting the current location
Point* GPOs::getKNN(Point *p, int k){
  int count = 0, incr=0;
  Point* neighbor;

  while(-1){
    res_point* next_neighbor = getNextNN(p->getX(), p->getY(), 25);
    incr++;

    if( next_neighbor->id != p->getID() && next_neighbor->uid != p->getUID() ){
      count++;
      if(count == k){
        neighbor = new Point(next_neighbor);
        clearNextNN();
        break;
      }
    }
  }

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

vector<res_point*>* GPOs::getRangeAndDelete(Point *p, double radius, double t_dist){
  vector<res_point*>* res = grid->getRangeAndDelete(p, radius, t_dist);
  return res;
}

vector<res_point*>* GPOs::getRange(Point *original, double radius, double t_dist){
  return grid->getRange(original, radius, t_dist);
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

void GPOs::getRangeByTime(int time_start, int time_end, vector<Point*>* results){
  for(int i= time_start; i<=time_end; i++){
    auto t_it = time_to_checkins.find(i);
    if(t_it != time_to_checkins.end()){
      vector <Point*> *checkins_at_i = t_it->second;
      results->insert(results->end(), checkins_at_i->begin(), checkins_at_i->end());
    }
  }
}

vector<res_point*>* GPOs::getRangeSortedByTime(double x, double y, double radius){
  vector<res_point*>* results = getRange(x, y, radius);
  res_point_checkin_time_comparator_ascending comporator;
  std::sort(results->begin(), results->end(), comporator);
  return results;
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


  checkin_list.insert(make_pair(l->getOrder(), l));

  grid->addCheckIn(l);

  ids.insert(l->getOrder());
};

// Load locations from file
bool GPOs::loadLocations(const char* fileName){
  cout << "Loading checkins from " << fileName << endl;
  ifstream fin(fileName);

  // set<int> *blacklist_set;
  set<int> *blacklist_set = new set<int>();

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

    if(y == 30.2635153076 && x == -97.7401685715)
      continue;

    dateTime = date + " " + time;
    boost::posix_time::ptime dtm = boost::posix_time::time_from_string(dateTime);

    // skip newlines
    if (!fin.good()) continue;

    if(blacklist_set->find(lid) == blacklist_set->end()){
      loadPoint(x, y, lid, uid, dtm, count);
      count ++;
    }

    if(count%1000000==0)
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

vector<int>* GPOs::getCheckinsInRangeByHourBlock(Point *p, double r){
  double radius_geo_dist = (r/1000) * 360 / EARTH_CIRCUMFERENCE;
  vector<res_point*>* checkins = getRange(p->getX(), p->getY(), radius_geo_dist);

  vector<int>* checkin_list = new vector<int>();

  for(auto c = checkins->begin(); c != checkins->end(); c++){
    Point q = Point(*c);
    if( abs((p->getTime() - q.getTime()).total_seconds()) <= 3600 ){
      checkin_list->push_back(q.getOrder());
    }
  }

  sort(checkin_list->begin(), checkin_list->end());
  checkin_list->erase( unique( checkin_list->begin(), checkin_list->end() ), checkin_list->end() );

  for(auto c = checkins->begin(); c != checkins->end(); c++){
    delete (*c);
  }
  delete checkins;

  return checkin_list;
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

vector<int>* GPOs::getLocationsInRange(double x, double y, double radius){
  vector<int> *locations = new vector<int>();
  double radius_geo_dist = (radius/1000) * 360 / EARTH_CIRCUMFERENCE;
  vector<res_point*>* checkins = getRange(x, y, radius_geo_dist);
  for(auto c = checkins->begin(); c != checkins->end(); c++){
    locations->push_back( (*c)->id );
    delete (*c);
  }
  delete checkins;
  // Sorting user list
  sort(locations->begin(), locations->end());
  // Removing duplicates
  locations->erase( unique( locations->begin(), locations->end() ), locations->end() );
  return locations;
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
// void GPOs::groupLocationsByKNNDistance(GPOs* gpos, int k, double std_radio){
//   double x=0, y=0;
//   unsigned int lid, count=0;
//   long unsigned int iterations=0;
//   unordered_set<int>* seenLocations = new unordered_set<int>();
//   boost::posix_time::ptime time;

//   GPOs *_duplicate_gpos = new GPOs(gpos);

//   // Shuffling the locations
//   std::random_shuffle( gpos->locations.begin(), gpos->locations.end() );

//   for(auto l = gpos->locations.begin(); l != gpos->locations.end(); l++){
//     Point *p = *l;
//     x   = p->getX();
//     y   = p->getY();
//     lid = p->getID();

//     // noise distance in meters
//     double neighbor_distance = gpos->getKNNDistance(p, k) * 1000;

//     vector<res_point*>* checkins = _duplicate_gpos->getRangeAndDelete(x, y, neighbor_distance * std_radio);

//     for(auto c = checkins->begin(); c != checkins->end(); c++){
//       if( seenLocations->find( (*c)->oid ) == seenLocations->end() ){
//         loadPoint(x, y, lid, (*c)->uid, (*c)->time, (*c)->oid);
//         seenLocations->insert( (*c)->oid );
//         count++;
//       }
//       iterations++;
//       delete (*c);
//     }
//     delete checkins;

//     if(count % 100000==0)
//       cout << count << " " << endl;

//   };

//   delete seenLocations;
//   delete _duplicate_gpos;

//   generateCooccurrenceCache();
// };

void GPOs::countCoOccurrencesOptimistic(){
  double radius_in_km = coocc_spatial_range / 1000.0;
  double time_deviation_in_hours = coocc_time_range / 3600.0;

  double radius_geo_dist = radius_in_km * 360 / EARTH_CIRCUMFERENCE,x=0, y=0;
  unsigned int count=0, order;

  for(auto l = checkin_list.begin(); l != checkin_list.end(); l++){
    Point *p = l->second;
    order    = l->first;
    x        = p->getX();
    y        = p->getY();

    vector<res_point*>* cooccurrences = getRange(p, radius_geo_dist, time_deviation_in_hours);

    unordered_set<int> *coocc_list = new unordered_set<int>();

    for(auto c = cooccurrences->begin(); c != cooccurrences->end(); c++){
      coocc_list->insert((*c)->oid);

      int o1 = order;
      int o2 = (*c)->oid;
      if(o1 > o2){
        int temp = o2;
        o2 = o1;
        o1 = temp;
      }

      cooccurred_checkins.insert(make_pair(o1, o2));

      delete (*c);
    }
    delete cooccurrences;


    if(coocc_list->size() > 0){
      cooccurrence_index.insert(make_pair(order, coocc_list));
    } else {
      delete coocc_list;
    }

    count++;
    if(count%100000 == 0)
      cout << count << endl;
  }

  double total_indirect_size = 0;
  for(auto c_it = cooccurrence_index.begin(); c_it != cooccurrence_index.end(); c_it++){
    int order = c_it->first;
    unordered_set<int>* direct = c_it->second;
    unordered_set<int>* indirect = new unordered_set<int>();
    for(auto d_it = direct->begin(); d_it != direct->end(); d_it++){
      int other = (*d_it);
      indirect->insert(other);
      unordered_set<int>* others_direct = cooccurrence_index.find(other)->second;
      for(auto od_it = others_direct->begin(); od_it != others_direct->end(); od_it++){
        indirect->insert((*od_it));
      }
    }
    cooccurrence_index_indirect.insert(make_pair(order, indirect));
    total_indirect_size+=indirect->size();
  }
  cout<<"Completed computing cooccurrences in optimistic manner" << endl;
  cout<<"total_cooccurrences{{"<<cooccurred_checkins.size()<<"}}"<<endl;
  cout<<"cooccurrence_index_size{{"<<cooccurrence_index.size()<<"}}"<<endl;
  cout<<"indirect_cooccurrence_index_size{{"<<cooccurrence_index_indirect.size()<<"}}"<<endl;
  cout<<"average_indirect_cooccurrences{{"<<total_indirect_size / cooccurrence_index_indirect.size()<<"}}"<<endl;
}

void GPOs::groupLocationsByST(GPOs* gpos, double radius_in_km, double time_deviation_in_hours){
  double radius_geo_dist = radius_in_km * 360 / EARTH_CIRCUMFERENCE,x=0, y=0;
  unsigned int count=0, order;

  unordered_set<int> seenLocations;
  boost::posix_time::ptime time;
  GPOs *_duplicate_gpos = new GPOs(gpos);

  cout << "Number of original checkins               : " << _duplicate_gpos->checkin_list.size() << endl;

  for(auto l = gpos->checkin_list.begin(); l != gpos->checkin_list.end(); l++){
    Point *p = l->second;
    x        = p->getX();
    y        = p->getY();
    order    = p->getOrder();
    time     = p->getTime();

    if( seenLocations.find( order ) != seenLocations.end() )
      continue;

    loadPoint(x, y, p->getID(), p->getUID(), time, order);
    seenLocations.insert( order );
    count++;

    vector<res_point*>* checkins = _duplicate_gpos->getRangeAndDelete(p, radius_geo_dist, time_deviation_in_hours);

    for(auto c = checkins->begin(); c != checkins->end(); c++){
      if( seenLocations.find( (*c)->oid ) == seenLocations.end() ){
        loadPoint(x, y, p->getID(), (*c)->uid, time, (*c)->oid);
        seenLocations.insert( (*c)->oid );
        count++;
        delete (*c);
      }
    }

    delete checkins;
  };

  cout << "Number of checkins after iteration     : " << count << endl;
  cout << "Check-ins inserted                     : " << seenLocations.size() << endl;
  cout << "Number of total checkins               : " << gpos->checkin_list.size() << endl;

  delete _duplicate_gpos;

  // generateCooccurrenceCache();

}

void GPOs::groupLocationsByDD(GPOs* gpos, set<int> *interested_checkins){
  double radius_geo_dist,x=0, y=0;
  unsigned int order;

  unordered_set<int> seenLocations;
  boost::posix_time::ptime time;

  double co_occurrences = 0;
  for(auto c_it = interested_checkins->begin(); c_it != interested_checkins->end(); c_it++){
    int checkin_order = (*c_it);

    auto p_it = gpos->checkin_list.find(checkin_order);

    Point *p = p_it->second;
    x        = p->getX();
    y        = p->getY();
    order    = p->getOrder();
    time     = p->getTime();

    if(p->getID() >= LOCATION_NOISE_BOUND)
      continue;

    loadPoint(x, y, p->getID(), p->getUID(), time, order);
    seenLocations.insert( order );

    // Get KNNs in SPATIAL_TEMPORAL [ RANGE ]
    priority_queue < pair<double, res_point*>, vector<pair<double, res_point*> > > spatioTemporalKNNs;
    vector<res_point*> *candidates = gpos->getRangeSpatioTemporalBound(p, SPATIAL_SOFT_BOUND, TEMPORAL_SOFT_BOUND);
    gpos->getSpatioTemporalKNN(p, 25, &spatioTemporalKNNs, candidates, 4);

    // Sort KNNs in ascending order
    map< double, res_point* > knnHash;
    while(!spatioTemporalKNNs.empty()){
      double distance = spatioTemporalKNNs.top().first;
      res_point* topK = spatioTemporalKNNs.top().second;
      knnHash.insert(make_pair(distance, topK));
      spatioTemporalKNNs.pop();
    }

    // Create co-occurrences until the first actual location
    vector<res_point*> cooccurrences;
    for(auto knn_it=knnHash.begin(); knn_it!=knnHash.end(); knn_it++){
      res_point *rp = knn_it->second;
      if(rp->id < LOCATION_NOISE_BOUND)
        break;
      if(seenLocations.find(rp->oid) == seenLocations.end()){
        co_occurrences++;
        loadPoint(x, y, p->getID(), rp->uid, time, rp->oid);
        seenLocations.insert(rp->oid);
      }
    }

    // deleting candidates
    for(auto sc_it=candidates->begin(); sc_it != candidates->end(); sc_it++){
      delete *sc_it;
    }
    delete candidates;

    if(checkin_list.size()%100000 == 0)
      cout << checkin_list.size() << " " << co_occurrences << endl;
  }

  for(auto c_it = gpos->checkin_list.begin(); c_it != gpos->checkin_list.end(); c_it++){
    Point *p = c_it->second;
    if( seenLocations.find(p->getOrder()) == seenLocations.end() )
      loadPoint(p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder());
  }

  cout << "Check-ins inserted : " << checkin_list.size()       << endl;
  cout << "Original size      : " << gpos->checkin_list.size() << endl;
}

void GPOs::computeLocationsOfInterest(double spatial_bound_in_meters, double time_bound_in_hours, set<int> *interested_checkins){
  set<int> interested_checkins_unordered;
  vector<pair<double, int>> interested_checkins_ordered;
  for(auto c_it = cooccurred_checkins.begin(); c_it != cooccurred_checkins.end(); c_it++){
    int o1 = c_it->first;
    int o2 = c_it->first;
    interested_checkins_unordered.insert(o1);
    interested_checkins_unordered.insert(o2);
  }

  cout << "Co-occurred check-ins : " << interested_checkins_unordered.size() << endl;

  for(auto i_it = interested_checkins_unordered.begin(); i_it != interested_checkins_unordered.end(); i_it++){
    int order = (*i_it);
    Point *q = checkin_list.find(order)->second;

    auto l_it = location_to_H.find(q->getID());
    double entropy=0;
    if(l_it != location_to_H.end())
      entropy = l_it->second;
    if(q->getID() < LOCATION_NOISE_BOUND){
      interested_checkins_ordered.push_back(make_pair(entropy, q->getOrder()));
    }

    // Get check-ins in vicinity
    double radius_geo_dist = (spatial_bound_in_meters/1000) * 360 / EARTH_CIRCUMFERENCE;
    vector<res_point*>* checkins_in_vicinity = getRange( q, radius_geo_dist, time_bound_in_hours );
    for(auto sc_it=checkins_in_vicinity->begin(); sc_it != checkins_in_vicinity->end(); sc_it++){
      auto l_it = location_to_H.find((*sc_it)->id);
      double entropy=0;
      if(l_it != location_to_H.end())
        entropy = l_it->second;
      if((*sc_it)->id < LOCATION_NOISE_BOUND && entropy != 0)
        interested_checkins_ordered.push_back( make_pair(entropy, (*sc_it)->oid) );
      delete *sc_it;
    }
    delete checkins_in_vicinity;
  }
  std::sort(interested_checkins_ordered.begin(), interested_checkins_ordered.end());
  std::reverse(interested_checkins_ordered.begin(), interested_checkins_ordered.end());

  for(auto k_it=interested_checkins_ordered.begin(); k_it != interested_checkins_ordered.end(); k_it++){
    interested_checkins->insert(k_it->second);
  }
  cout << "Computed check-ins of interest to perform advanced restoration : " << interested_checkins->size() << endl;
}

void GPOs::groupLocationsToTopK(GPOs* gpos, set<int> *interested_checkins, int k, double spatial_bound_in_meters, double temporal_bound_in_hours){
  double x,y,order;
  boost::posix_time::ptime time;
  double co_occurrences=0;
  unordered_set<int> seenLocations;
  for(auto c_it = interested_checkins->begin(); c_it != interested_checkins->end(); c_it++){
    int checkin_order = (*c_it);
    auto p_it = gpos->checkin_list.find(checkin_order);

    Point *p = p_it->second;
    x        = p->getX();
    y        = p->getY();
    order    = p->getOrder();
    time     = p->getTime();

    if(p->getID() >= LOCATION_NOISE_BOUND)
      continue;

    if(seenLocations.find(p->getOrder()) == seenLocations.end()){
      loadPoint(x, y, p->getID(), p->getUID(), time, order);
      seenLocations.insert(p->getOrder());
    }

    priority_queue < pair<double, res_point*>, vector<pair<double, res_point*> > > spatioTemporalKNNs;
    vector<res_point*> *candidates = gpos->getRangeSpatioTemporalBound(p, spatial_bound_in_meters, temporal_bound_in_hours);
    gpos->getSpatioTemporalKNN(p, k, &spatioTemporalKNNs, candidates, 4);

    // KNN in bound
    if(spatioTemporalKNNs.size() == k){
      res_point* topK = spatioTemporalKNNs.top().second;
      if(seenLocations.find(topK->oid) == seenLocations.end() && topK->id >= LOCATION_NOISE_BOUND){
        co_occurrences++;
        loadPoint(x, y, p->getID(), topK->uid, time, topK->oid);
        seenLocations.insert(topK->oid);
      }
    }

    // deleting candidates
    for(auto sc_it=candidates->begin(); sc_it != candidates->end(); sc_it++){
      delete *sc_it;
    }
    delete candidates;

    if(seenLocations.size() % 100000==0)
      cout << seenLocations.size() << " " << co_occurrences << endl;
  }

  for(auto c_it = gpos->checkin_list.begin(); c_it != gpos->checkin_list.end(); c_it++){
    Point *p = c_it->second;
    if( seenLocations.find(p->getOrder()) == seenLocations.end() )
      loadPoint(p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder());
  }

  cout << "Artificially created co-occurrences : " << co_occurrences << endl;
  cout << "Check-ins inserted : " << checkin_list.size()       << endl;
  cout << "Original size      : " << gpos->checkin_list.size() << endl;
}

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

    if(count % 1000000==0)
      cout << count << " " << endl;

  };
  cout << "Check-ins inserted : " << count << endl;
  cout << "Check-ins failed lookup : " << iterations << endl;

  delete seenLocations;
  delete _duplicate_gpos;

  // generateCooccurrenceCache();
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
  cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
  cout<<"average_spatial_displacement{{"<< (total_spatial_displacement / point_count) * 1000  <<"}} in meters"<<endl;
  cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / spatial_purturbed_count) * 1000 <<"}} in meters"<<endl;
}

void GPOs::computeSkylineMetrics(map< int, map<int,int>* >* _location_to_user_to_cooccurrences){
  ofstream outfile;
  stringstream ss;
  std::string filePath;
  ss << "skyline-coocc.csv";

  filePath = ss.str();
  outfile.open( filePath.c_str() );

  int checkin_count = 0;

  set<int> checkins_of_interest;
  pickUniqueCheckinFromCooccurrences(&checkins_of_interest);
  cout << "Checkins of interest : " << checkins_of_interest.size() << endl;

  for(auto c_it = checkins_of_interest.begin(); c_it != checkins_of_interest.end(); c_it++){
    int order = (*c_it);
    Point *p = checkin_list.find(order)->second;
    vector <res_point*> *candidates = getRangeSpatioTemporalBound(p);

    map< int, pair<int, res_point*> > skylines;
    getSkylinePoints(p, candidates, &skylines);

    for(auto sk_it=skylines.begin(); sk_it != skylines.end(); sk_it++){
      int skyline_id = sk_it->first;
      res_point *skyline = sk_it->second.second;
      int domination_count = sk_it->second.first;

      outfile << p->getOrder() << " " << skyline->oid << " " << domination_count << " ";
      outfile << p->computeMinDistInKiloMeters(skyline->x, skyline->y) << " ";
      outfile << (double) p->getTimeDifference(skyline) / 3600.0 << endl;
    }

    for(auto sc_it=candidates->begin(); sc_it != candidates->end(); sc_it++){
      delete *sc_it;
    }
    delete candidates;
    checkin_count++;
    if(checkin_count % 10000 == 0)
      cout << checkin_count << endl;
  }

  outfile.close();
}

void GPOs::pickSingleCheckinFromCooccurrences(set<int> *checkins_of_interest){
  for(auto c_it = cooccurred_checkins.begin(); c_it != cooccurred_checkins.end(); c_it++){
    int o1 = c_it->first;
    int o2 = c_it->second;

    unordered_set<int>* o1_orders = cooccurrence_index.find(o1)->second;
    unordered_set<int>* o2_orders = cooccurrence_index.find(o2)->second;

    int order; // Dominating check-in
    if(o1_orders->size() >= o2_orders->size())
      order = o1;
    else
      order = o2;

    if( checkins_of_interest->find(order) == checkins_of_interest->end() )
      checkins_of_interest->insert(order);
  }
}

void GPOs::pickOtherCheckinFromCooccurrences(set<int> *checkins_of_interest){
  for(auto c_it = cooccurred_checkins.begin(); c_it != cooccurred_checkins.end(); c_it++){
    int o1 = c_it->first;
    int o2 = c_it->second;

    checkins_of_interest->insert(o2);
  }
}

void GPOs::pickUniqueCheckinFromCooccurrences(set<int> *checkins_of_interest){
  for(auto c_it = cooccurred_checkins.begin(); c_it != cooccurred_checkins.end(); c_it++){
    int o1 = c_it->first;
    int o2 = c_it->second;
    checkins_of_interest->insert(o1);
    checkins_of_interest->insert(o2);
  }
}

void GPOs::computeSTKNNDistances(int k, int type){
  ofstream outfile;
  stringstream ss;
  std::string filePath;

  timeval A_start,A_end,B_start,B_end,C_start,C_end,D_start,D_end;

  if(type == 0)
    ss << "knn-noise-combined-" << k <<"-" << coocc_spatial_range << "-" << coocc_time_range << "-coocc" << ".csv";

  if(type == 1)
    ss << "knn-noise-spatial-" << k <<"-" << coocc_spatial_range << "-" << coocc_time_range << "-coocc" << ".csv";

  if(type == 2)
    ss << "knn-noise-temporal-" << k <<"-" << coocc_spatial_range << "-" << coocc_time_range << "-coocc" << ".csv";

  if(type == 3)
    ss << "knn-noise-combined-" << k <<"-" << coocc_spatial_range << "-" << coocc_time_range << "-coocc" << ".csv";

  filePath = ss.str();
  outfile.open( filePath.c_str() );
  int checkin_count = 0;

  cout << "Writing to file : " << ss.str() << endl;

  set<int> *checkins_of_interest;

  if(type != 3){
    checkins_of_interest = new set<int>();
    pickUniqueCheckinFromCooccurrences(checkins_of_interest);
    cout << "Checkins of interest : " << checkins_of_interest->size() << endl;
  } else {
    checkins_of_interest = &ids;
    cout << "Using all checkins : " << checkins_of_interest->size() << endl;
  }

  for(auto c_it = checkins_of_interest->begin(); c_it != checkins_of_interest->end(); c_it++){
    int order = (*c_it);

    Point *p = checkin_list.find(order)->second;
    vector <res_point*> *candidates = getRangeSpatioTemporalBound(p);

    priority_queue < pair<double, res_point*>, vector<pair<double, res_point*> > > spatioTemporalKNNs;
    getSpatioTemporalKNN(p, k, &spatioTemporalKNNs, candidates, type);

    int result_size = spatioTemporalKNNs.size();

    outfile << p->getOrder() << " ";

    for(int i=result_size; i < k; i++ ){
      outfile << -1 << " ";
      outfile << std::numeric_limits<double>::infinity() << " ";
      outfile << std::numeric_limits<double>::infinity() << " ";
      outfile << std::numeric_limits<double>::infinity() << " ";
    }

    while( !spatioTemporalKNNs.empty() ){
      res_point* candidate = spatioTemporalKNNs.top().second;
      outfile << candidate->oid << " ";
      outfile << spatioTemporalKNNs.top().first << " ";
      outfile << p->computeMinDistInKiloMeters(candidate->x, candidate->y) << " ";
      outfile << (double) p->getTimeDifference(candidate) / 3600.0 << " ";
      spatioTemporalKNNs.pop();
    }
    outfile << endl;

    for(auto cd_it = candidates->begin(); cd_it != candidates->end(); cd_it++){
      delete(*cd_it);
    }
    delete candidates;

    checkin_count++;
    if(checkin_count % 10000 == 0)
      cout << checkin_count << endl;

  }

  delete checkins_of_interest;
  outfile.close();
}


void GPOs::hideCoLocations(GPOs* gpos, double prob){
  set<int> checkins_of_interest;
  gpos->pickSingleCheckinFromCooccurrences(&checkins_of_interest);

  cout << "Purtubing only " << prob * checkins_of_interest.size() << " checkins" << endl;

  unsigned int point_count = 0, lid=LOCATION_NOISE_BOUND;

  for(auto c_it = gpos->checkin_list.begin(); c_it != gpos->checkin_list.end(); c_it++){
    int order = c_it->first;
    Point *p = c_it->second;

    if( checkins_of_interest.find(order) != checkins_of_interest.end() && (rand()%101)<=prob*100 ){
      lid++;
      purturbed_count++;
      spatial_purturbed_count++;
      temporal_purturbed_count++;
      continue;
    }

    loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
    point_count++;
  }

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"spatially_purtubed_checkins{{"<< spatial_purturbed_count   << "}}" << endl;
  cout<<"temporally_purtubed_checkins{{"<< temporal_purturbed_count << "}}" << endl;
  cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
  cout<<"average_spatial_displacement{{"<< (total_spatial_displacement / point_count) * 1000  <<"}} in meters"<<endl;
  cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / spatial_purturbed_count) * 1000 <<"}} in meters"<<endl;
  cout<<"total_temporal_displacement{{"<< total_time_displacement <<"}} hours"<<endl;
  cout<<"average_temporal_displacement{{"<< total_time_displacement  * (1/(float)point_count) * 3600 <<"}} seconds"<<endl;
  cout<<"average_temporal_displacement_on_purtubed{{"<< total_time_displacement * (1/(float)temporal_purturbed_count) <<"}} hours"<<endl;
}

void GPOs::dummyCoLocations(GPOs* gpos, int k){
  set<int> checkins_of_interest;
  gpos->pickSingleCheckinFromCooccurrences(&checkins_of_interest);

  cout << "Adding  " << k << "dummy co-locations per co-location" << endl;
  unsigned int point_count = 0, lid=LOCATION_NOISE_BOUND;

  for(auto c_it = gpos->checkin_list.begin(); c_it != gpos->checkin_list.end(); c_it++){
    int order = c_it->first;
    Point *p = c_it->second;

    if( checkins_of_interest.find(order) != checkins_of_interest.end() ){
      for(int j=1; j<=k; j++){
        pair <double,double> pt = util.addGaussianNoise(p->getX(),p->getY(), 50, 0);
        loadPoint( pt.first, pt.second, lid, lid, p->getTime(), lid );
        lid++;
        purturbed_count++;
        spatial_purturbed_count++;
        temporal_purturbed_count++;

        loadPoint( pt.first, pt.second, lid, lid, p->getTime(), lid );
        lid++;
        purturbed_count++;
        spatial_purturbed_count++;
        temporal_purturbed_count++;
      }
    }
    loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
    point_count++;
  }

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"spatially_purtubed_checkins{{"<< spatial_purturbed_count   << "}}" << endl;
  cout<<"temporally_purtubed_checkins{{"<< temporal_purturbed_count << "}}" << endl;
  cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
  cout<<"average_spatial_displacement{{"<< (total_spatial_displacement / point_count) * 1000  <<"}} in meters"<<endl;
  cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / spatial_purturbed_count) * 1000 <<"}} in meters"<<endl;
  cout<<"total_temporal_displacement{{"<< total_time_displacement <<"}} hours"<<endl;
  cout<<"average_temporal_displacement{{"<< total_time_displacement  * (1/(float)point_count) * 3600 <<"}} seconds"<<endl;
  cout<<"average_temporal_displacement_on_purtubed{{"<< total_time_displacement * (1/(float)temporal_purturbed_count) <<"}} hours"<<endl;
}

// Only co-occurrences
void GPOs::loadPurturbedBasedOnSelectiveGaussian(GPOs* gpos, double radius, uint time_deviation){
  set<int> checkins_of_interest;
  gpos->pickSingleCheckinFromCooccurrences(&checkins_of_interest);

  cout << "Purtubing only " << checkins_of_interest.size() << " checkins" << endl;

  unsigned int point_count = 0, lid=LOCATION_NOISE_BOUND;
  double min_spatial_noise_added=std::numeric_limits<double>::infinity(), min_temporal_noise_added=std::numeric_limits<double>::infinity();

  for(auto c_it = gpos->checkin_list.begin(); c_it != gpos->checkin_list.end(); c_it++){
    int order = c_it->first;
    Point *p = c_it->second;

    if( checkins_of_interest.find(order) != checkins_of_interest.end() ){
      pair<double,double> coordinates_with_noise = util.addGaussianNoise( p->getX(), p->getY(), radius,  0 );
      boost::posix_time::ptime purtubed_time = util.addTemporalGaussianNoise( p->getTime(), time_deviation, 0 );

      double sd = p->computeMinDistInKiloMeters(coordinates_with_noise.first, coordinates_with_noise.second);
      double td = (double) abs( (p->getTime() - purtubed_time).total_seconds() ) / 3600.0;
      total_spatial_displacement += sd;
      total_time_displacement += td;

      if(sd < min_spatial_noise_added)
        min_spatial_noise_added  = sd;

      if(td < min_temporal_noise_added)
        min_temporal_noise_added = td;

      loadPoint( coordinates_with_noise.first, coordinates_with_noise.second, lid, p->getUID(), purtubed_time, p->getOrder() );
      lid++;
      purturbed_count++;
      spatial_purturbed_count++;
      temporal_purturbed_count++;
    } else {
      loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
    }

    point_count++;
  }

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"spatially_purtubed_checkins{{"<< spatial_purturbed_count   << "}}" << endl;
  cout<<"temporally_purtubed_checkins{{"<< temporal_purturbed_count << "}}" << endl;
  cout<<"min_spatial_noise_added{{"<<  min_spatial_noise_added <<"}} in km"<<endl;
  cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
  cout<<"average_spatial_displacement{{"<< (total_spatial_displacement / point_count) * 1000  <<"}} in meters"<<endl;
  cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / spatial_purturbed_count) * 1000 <<"}} in meters"<<endl;
  cout<<"min_temporal_noise_added{{"<<  min_temporal_noise_added <<"}} hours"<<endl;
  cout<<"total_temporal_displacement{{"<< total_time_displacement <<"}} hours"<<endl;
  cout<<"average_temporal_displacement{{"<< total_time_displacement  * (1/(float)point_count) * 3600 <<"}} seconds"<<endl;
  cout<<"average_temporal_displacement_on_purtubed{{"<< total_time_displacement * (1/(float)temporal_purturbed_count) <<"}} hours"<<endl;
  cout<<"Locations after perturbation :"<<location_to_user.size()<<endl;
}

void GPOs::anonymizeBasedOnSelectiveSTKNNDistance(GPOs* gpos, int k, bool hide){
  map <int, vector<double>* > st_knn;
  stringstream ss;
  ss << "knn-noise-combined-100-" << coocc_spatial_range << "-" << coocc_time_range << "-coocc" << ".csv";
  ifstream fin(ss.str());

  while(!fin.eof()){
    int order;
    vector<double> *neighbours = new vector<double>();
    fin >> order;

    for(int i = 0; i<100; i++){
      int knn_order;
      double st_distance, s_distance, t_distance;
      string text;
      fin >> knn_order;
      fin >> text;
      st_distance = (text == "inf") ? std::numeric_limits<double>::infinity() : atof(text.c_str());
      fin >> text;
      s_distance = (text == "inf") ? std::numeric_limits<double>::infinity() : atof(text.c_str());
      fin >> text;
      t_distance = (text == "inf") ? std::numeric_limits<double>::infinity() : atof(text.c_str());
      if(knn_order != -1)
        neighbours->push_back(knn_order);
    }

    if(neighbours->size() > 0){
      reverse(neighbours->begin(),neighbours->end());
      st_knn.insert(make_pair(order, neighbours));
    }
    else
      delete neighbours;
  }
  fin.close();
  cout << "Loaded ST_KNN from " << ss.str() << " : " << st_knn.size() << endl;

  unsigned int lid=LOCATION_NOISE_BOUND, cooccurrences_out_of_bound=0, knn_not_added=0;
  double sd, td;

  set <int> seenLocations, checkins_of_interest;
  gpos->pickSingleCheckinFromCooccurrences(&checkins_of_interest);

  for(auto c_it = checkins_of_interest.begin(); c_it != checkins_of_interest.end(); c_it++){
    int order = (*c_it);
    Point *p = gpos->checkin_list.find(order)->second;

    double baseX=p->getX(),
           baseY=p->getY();
    boost::posix_time::ptime baseTime= p->getTime();

    auto knn_it = st_knn.find(order);

    if( seenLocations.find(p->getOrder()) == seenLocations.end() ){
      loadPoint( p->getX(), p->getY(), lid, p->getUID(), p->getTime(), p->getOrder() );
      seenLocations.insert(p->getOrder());
      lid++;
      total_spatial_displacement+=p->computeMinDistInKiloMeters(p->getX(), p->getY());
      total_time_displacement+=(double)abs((p->getTime() - p->getTime()).total_seconds())/3600.0;
      purturbed_count++;
    } else {
      continue;
    }

    if(knn_it == st_knn.end()){ // Hide very sparse
      cooccurrences_out_of_bound++;
      continue;
    }

    vector<double> *neighbours = knn_it->second;
    int kth = k, knn_added = 0;
    for(int i=1; i<=kth && i<=neighbours->size(); i++){
      int neighbor = neighbours->at(i-1);
      Point tp = Point(baseX, baseY, -1);
      Point *q = gpos->checkin_list.find(neighbor)->second;
      // bool no_cooccurrence_at_neighbour = (gpos->cooccurrence_index.find(neighbor) == gpos->cooccurrence_index.end());
      if( seenLocations.find(q->getOrder()) == seenLocations.end() ){
        loadPoint( baseX, baseY, lid, q->getUID(), baseTime, q->getOrder() );
        seenLocations.insert(q->getOrder());
        lid++;
        total_spatial_displacement+=tp.computeMinDistInKiloMeters(q->getX(), q->getY());
        total_time_displacement+=(double) abs( (q->getTime() - baseTime).total_seconds() ) / 3600.0;
        purturbed_count++;
        knn_added++;
      } else {
        kth++; // Check the next NN
      }
    }

    if(knn_added == 0) // Keeping track of co-occurrences which are not anonomized
      knn_not_added++;
  }

  for(auto c_it = gpos->checkin_list.begin(); c_it != gpos->checkin_list.end(); c_it++){
    int order = c_it->first;
    Point *p = c_it->second;

    if(seenLocations.find(order) == seenLocations.end()){
      loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
    }
  }

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"cooccurrences_out_of_bound{{"<< cooccurrences_out_of_bound << "}}" << endl;
  cout<<"knn_not_added{{"<< knn_not_added << "}}" << endl;
  cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
  cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / purturbed_count) * 1000 <<"}} in meters"<<endl;
  cout<<"total_temporal_displacement{{"<< total_time_displacement <<"}} hours"<<endl;
  cout<<"average_temporal_displacement_on_purtubed{{"<< total_time_displacement * (1/(float)purturbed_count) <<"}} hours"<<endl;

  for(auto knn_it=st_knn.begin(); knn_it != st_knn.end(); knn_it++){
    delete knn_it->second;
  }
}

void GPOs::loadPurturbedBasedOnSelectiveSTKNNDistance(GPOs* gpos, int k, bool hide){
  loadPurturbedBasedOnSelectiveSTKNNDistance(gpos, k, true, hide);
}

// Only co-occurrences
void GPOs::loadPurturbedBasedOnSelectiveSTKNNDistance(GPOs* gpos, int k, bool gaussian, bool hide){
  unsigned int point_count = 0, lid=LOCATION_NOISE_BOUND, cooccurrences_out_of_bound=0;

  map <int, vector<double>* > st_knn;
  stringstream ss;
  ss << "knn-noise-combined-100-" << coocc_spatial_range << "-" << coocc_time_range << "-coocc" << ".csv";
  ifstream fin(ss.str());

  while(!fin.eof()){
    int order;
    vector<double> *neighbours = new vector<double>();
    fin >> order;

    for(int i = 0; i<100; i++){
      int knn_order;
      double st_distance, s_distance, t_distance;
      string text;
      fin >> knn_order;
      fin >> text;
      st_distance = (text == "inf") ? std::numeric_limits<double>::infinity() : atof(text.c_str());
      fin >> text;
      s_distance = (text == "inf") ? std::numeric_limits<double>::infinity() : atof(text.c_str());
      fin >> text;
      t_distance = (text == "inf") ? std::numeric_limits<double>::infinity() : atof(text.c_str());
      if(knn_order != -1)
        neighbours->push_back(knn_order);
    }

    if(neighbours->size() > 0){
      reverse(neighbours->begin(),neighbours->end());
      st_knn.insert(make_pair(order, neighbours));
    }
    else
      delete neighbours;
  }
  fin.close();
  cout << "Loaded ST_KNN from " << ss.str() << " : " << st_knn.size() << endl;

  set<int> checkins_of_interest;
  gpos->pickSingleCheckinFromCooccurrences(&checkins_of_interest);

  cout << "Loaded checkins of interest : " << checkins_of_interest.size() << endl;

  for(auto c_it = gpos->checkin_list.begin(); c_it != gpos->checkin_list.end(); c_it++){
    int order = c_it->first;
    Point *p = c_it->second;
    auto knn_it = st_knn.find(order);
    bool checkin_of_interest = (checkins_of_interest.find(order) != checkins_of_interest.end());
    bool knn_out_of_hard_bound = ( knn_it == st_knn.end() );

    if(!checkin_of_interest){
      loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
      point_count++;
      continue;
    }

    Point *base_checkin;
    double noise_radius, time_deviation;

    if(knn_out_of_hard_bound){
      noise_radius = 2 * SPATIAL_HARD_BOUND;
      time_deviation = 2 * TEMPORAL_HARD_BOUND * 3600.0;
      base_checkin = p;
    }
    else {
      vector<double> *neighbours = knn_it->second;
      int k_lim = (neighbours->size() < k) ? neighbours->size() : k;
      int kth = rand() % (k_lim+1);
      int neighbor, neighbours_neighbour;
      if(kth==0){
        neighbor = neighbours->at(0);
        Point *n = gpos->checkin_list.find(neighbor)->second;
        noise_radius = 2 * p->computeMinDistInKiloMeters(n->getX(), n->getY()) * 1000;
        time_deviation = 2 * abs((p->getTime() - n->getTime()).total_seconds());
        base_checkin = p;
      } else {
        neighbor = neighbours->at(kth-1);
        Point *n = gpos->checkin_list.find(neighbor)->second;
        vector <res_point*> *candidates = gpos->getRangeSpatioTemporalBound(n, SPATIAL_SOFT_BOUND, TEMPORAL_SOFT_BOUND);
        if(candidates->size() != 0){
          priority_queue < pair<double, res_point*>, vector<pair<double, res_point*> > > spatioTemporalKNNs;
          gpos->getSpatioTemporalKNN(n, 1, &spatioTemporalKNNs, candidates, 4);
          res_point* q = spatioTemporalKNNs.top().second;
          noise_radius = 2 * n->computeMinDistInKiloMeters(q->x, q->y) * 1000;
          time_deviation = 2 * abs((n->getTime() - q->time).total_seconds());
        } else {
          noise_radius = 2 * SPATIAL_SOFT_BOUND;
          time_deviation = 2 * TEMPORAL_SOFT_BOUND * 3600.0;
        }
        base_checkin = n;
      }
    };

    if(!gaussian){
      noise_radius = 0;
      time_deviation = 0;
    }

    pair<double,double> coordinates_with_noise = util.addGaussianNoise(base_checkin->getX(), base_checkin->getY(), noise_radius);
    boost::posix_time::ptime purtubed_time = util.addTemporalGaussianNoise(base_checkin->getTime(), time_deviation);
    double sd = p->computeMinDistInKiloMeters(coordinates_with_noise.first, coordinates_with_noise.second);
    double td = (double) abs( (p->getTime() - purtubed_time).total_seconds() ) / 3600.0;

    bool knn_out_of_soft_bound = false;
    if(!knn_out_of_hard_bound){
      vector<double> *neighbours = knn_it->second;
      int neighbor = neighbours->at(0);
      Point *q = gpos->checkin_list.find(neighbor)->second;
      knn_out_of_soft_bound = (p->computeMinDistInKiloMeters(q->getX(), q->getY())*1000.0 >= SPATIAL_SOFT_BOUND);
      knn_out_of_soft_bound = knn_out_of_soft_bound || ( abs((p->getTime() - q->getTime()).total_seconds())/3600.0 >= TEMPORAL_SOFT_BOUND );
    }

    if(knn_out_of_soft_bound || knn_out_of_hard_bound){
      total_spatial_displacement_sparse+=sd;
      total_time_displacement_sparse+=td;
      sparse_purturbed_count++;
    } else {
      total_spatial_displacement_dense+=sd;
      total_time_displacement_dense+=td;
      dense_purturbed_count++;
    }

    total_spatial_displacement+=sd;
    total_time_displacement+=td;

    spatial_displacement.push_back(sd);
    temporal_displacement.push_back(td);

    purturbed_count++;
    spatial_purturbed_count++;
    temporal_purturbed_count++;

    loadPoint( coordinates_with_noise.first, coordinates_with_noise.second, lid, p->getUID(), purtubed_time, p->getOrder() );
    lid++;
    point_count++;
  }

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"cooccurrences_out_of_bound{{"<< cooccurrences_out_of_bound << "}}" << endl;
  cout<<"spatially_purtubed_checkins{{"<< spatial_purturbed_count   << "}}" << endl;
  cout<<"temporally_purtubed_checkins{{"<< temporal_purturbed_count << "}}" << endl;
  cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
  cout<<"average_spatial_displacement{{"<< (total_spatial_displacement / point_count) * 1000  <<"}} in meters"<<endl;
  cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / spatial_purturbed_count) * 1000 <<"}} in meters"<<endl;
  cout<<"median_spatial_displacement{{"<< Median(spatial_displacement.begin(), spatial_displacement.end()) * 1000 <<"}} in meters"<<endl;
  cout<<"total_temporal_displacement{{"<< total_time_displacement <<"}} hours"<<endl;
  cout<<"average_temporal_displacement{{"<< total_time_displacement  * (1/(float)point_count) * 3600 <<"}} seconds"<<endl;
  cout<<"average_temporal_displacement_on_purtubed{{"<< total_time_displacement * (1/(float)temporal_purturbed_count) <<"}} hours"<<endl;
  cout<<"median_temporal_displacement{{"<< Median(temporal_displacement.begin(), temporal_displacement.end()) * 3600 <<"}} seconds"<<endl;
  cout<<"dense_purturbed_count{{"<< dense_purturbed_count <<"}}"<<endl;
  cout<<"total_spatial_displacement_dense{{"<< total_spatial_displacement_dense <<"}}"<<endl;
  cout<<"total_time_displacement_dense{{"<< total_time_displacement_dense <<"}}"<<endl;
  cout<<"average_spatial_displacement_dense{{"<< total_spatial_displacement_dense/(double)dense_purturbed_count*1000.0 <<"}}"<<endl;
  cout<<"average_time_displacement_dense{{"<< total_time_displacement_dense/(double)dense_purturbed_count*3600.0 <<"}}"<<endl;
  cout<<"sparse_purturbed_count{{"<< sparse_purturbed_count <<"}}"<<endl;
  cout<<"total_spatial_displacement_sparse{{"<< total_spatial_displacement_sparse <<"}}"<<endl;
  cout<<"total_time_displacement_sparse{{"<< total_time_displacement_sparse <<"}}"<<endl;
  cout<<"average_spatial_displacement_sparse{{"<< total_spatial_displacement_sparse/(double)sparse_purturbed_count*1000.0 <<"}}"<<endl;
  cout<<"average_time_displacement_sparse{{"<< total_time_displacement_sparse/(double)sparse_purturbed_count*3600.0 <<"}}"<<endl;
}

// Only co-occurrences
void GPOs::loadPurturbedBasedOnSelectiveSkyline(GPOs* gpos, int k){
  set<int> checkins_of_interest;
  gpos->pickSingleCheckinFromCooccurrences(&checkins_of_interest);

  unsigned int point_count = 0, lid=LOCATION_NOISE_BOUND;

  map <int, priority_queue< pair<int, int>, vector<pair<int, int>> >* > skyline_map;
  ifstream fin("skyline-coocc.csv");

  while(!fin.eof()){
    int order, skyline_id, d_count;
    double s_distance, t_distance;
    priority_queue< pair<int, int>, vector<pair<int, int>> >* skyline_queue;

    fin >> order >> skyline_id >> d_count >> s_distance >> t_distance;

    auto sk_it = skyline_map.find(order);
    if(sk_it  == skyline_map.end() ){
      skyline_queue = new priority_queue< pair<int, int>, vector<pair<int, int>> >();
      skyline_queue->push(make_pair(d_count, skyline_id));
      skyline_map.insert(make_pair(order, skyline_queue));
    } else {
      skyline_queue = sk_it->second;
      skyline_queue->push(make_pair(d_count, skyline_id));
    }
  }

  fin.close();

  cout << "Loaded skylines " << skyline_map.size() << endl;

  for(auto c_it = gpos->checkin_list.begin(); c_it != gpos->checkin_list.end(); c_it++){
    int order = c_it->first;
    Point *p = c_it->second;

    auto sk_it = skyline_map.find(order);

    if( checkins_of_interest.find(order) != checkins_of_interest.end() && sk_it != skyline_map.end() ){
      priority_queue< pair<int, int>, vector<pair<int, int>>> *skylines = sk_it->second;

      int kth;
      if(k == -1){
        // Pick a Skyline at random
        kth = rand() % skylines->size();
      } else {
        // Pick 1 of top k at random
        int k_lim = (skylines->size() < k) ? skylines->size() : k;
        kth = rand() % k_lim;
      }

      for(int i=0; i<kth; i++)
        skylines->pop();

      int skyline_id = skylines->top().second;
      Point *skyline = gpos->checkin_list.find(skyline_id)->second;

      double noise_radius = p->computeMinDistInKiloMeters(skyline->getX(), skyline->getY()) * 1000;
      double time_deviation = abs((p->getTime() - skyline->getTime()).total_seconds());

      pair<double,double> coordinates_with_noise = util.addGaussianNoise(skyline->getX(), skyline->getY(), noise_radius);
      boost::posix_time::ptime purtubed_time = util.addTemporalGaussianNoise(skyline->getTime(), time_deviation);

      total_spatial_displacement+=p->computeMinDistInKiloMeters(coordinates_with_noise.first, coordinates_with_noise.second);
      total_time_displacement+= (double) abs( (p->getTime() - purtubed_time).total_seconds() ) / 3600.0;

      loadPoint( coordinates_with_noise.first, coordinates_with_noise.second, lid, p->getUID(), purtubed_time, p->getOrder() );

      lid++;
      purturbed_count++;
      spatial_purturbed_count++;
      temporal_purturbed_count++;

    } else {
      loadPoint( p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder() );
    }

    point_count++;
  }

  for(auto sk_it = skyline_map.begin(); sk_it != skyline_map.end(); sk_it++)
    delete sk_it->second;

  cout<<"purtubed_checkins{{"<< purturbed_count << "}}" << endl;
  cout<<"spatially_purtubed_checkins{{"<< spatial_purturbed_count   << "}}" << endl;
  cout<<"temporally_purtubed_checkins{{"<< temporal_purturbed_count << "}}" << endl;

  cout<<"total_spatial_displacement{{"<<  total_spatial_displacement <<"}} in km"<<endl;
  cout<<"average_spatial_displacement{{"<< (total_spatial_displacement / point_count) * 1000  <<"}} in meters"<<endl;
  cout<<"average_spatial_displacement_on_purtubed{{"<< (total_spatial_displacement / spatial_purturbed_count) * 1000 <<"}} in meters"<<endl;

  cout<<"total_temporal_displacement{{"<< total_time_displacement <<"}} hours"<<endl;
  cout<<"average_temporal_displacement{{"<< total_time_displacement  * (1/(float)point_count) * 3600 <<"}} seconds"<<endl;
  cout<<"average_temporal_displacement_on_purtubed{{"<< total_time_displacement * (1/(float)temporal_purturbed_count) <<"}} hours"<<endl;
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

  // generateCooccurrenceCache();
}

// Compute CoOccurrences
void GPOs::countU2UCoOccurrences(){
  int total_cooccurrences=0, count=0;

  int user_checkins_size = 0;

  set <pair<int,int>> cooccurred_checkins; // TODO FIX THIS

  cout<<"Number of locations: "<<locations_users_frequency_map_with_order.size()<<endl;

  for(auto l_it=locations_users_frequency_map_with_order.begin(); l_it != locations_users_frequency_map_with_order.end(); l_it++){
    int location_id = l_it->first;
    map<int, vector< pair<uint, int> >* > *user_checkin_times = l_it->second;

    if(user_checkin_times->size() == 1){
      continue;
    }
    count++;

    user_checkins_size = user_checkins_size + user_checkin_times->size();

    // if(count % 10000 == 0)
    //   cout << count << endl;

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

        int intersection_count = util.getCooccurrencesWithinTimeBlock(u1_timestamps, u2_timestamps, coocc_time_range, &cooccurred_checkins);

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
      }

    }
  }

  cout << "Completed co-occurrences computation " << endl;

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

  cout << "Building co-occurrence index " << endl;

  for(auto c_it = cooccurred_checkins.begin(); c_it != cooccurred_checkins.end(); c_it++){
    int o1 = c_it->first;
    int o2 = c_it->second;
    unordered_set<int> *lst;

    auto o1_it = cooccurrence_index.find(o1);
    if(o1_it == cooccurrence_index.end()){
      lst = new unordered_set<int>();
      cooccurrence_index.insert(make_pair(o1, lst));
    } else {
      lst = o1_it->second;
    }
    lst->insert(o2);

    auto o2_it = cooccurrence_index.find(o2);
    if(o2_it == cooccurrence_index.end()){
      lst = new unordered_set<int>();
      cooccurrence_index.insert(make_pair(o2, lst));
    } else {
      lst = o2_it->second;
    }
    lst->insert(o1);
  }

  cout << "users_with_atleast_one_cooccurrence{{" << cooccurrence_matrix.size() << "}}" << endl;
  cout << "cooccurred_user_paris{{" << cooccured_user_pairs.size() << "}}" << endl;
  cout << "significantly_cooccurred_user_paris{{" << significantly_cooccured_user_pairs.size() << "}}" << endl;
  cout << "total_cooccurrence{{" << cooccurred_checkins.size() << "}}" << endl;
  cout << "cooccurrence_index_size{{" << cooccurrence_index.size() << "}}" << endl;
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

  // generateCooccurrenceCache();
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


