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
      auto user_it = user_frequencies->find(uid);
      boost::posix_time::time_duration time_difference = point->getTime() - time_t_epoch;
      int time = abs(time_difference.total_seconds());

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
        return NULL;
    }
}


vector<res_point*>* GPOs::getkNN(double x, double y, int k){

    //	cout << "getKNN start" << endl;
    clock_t startC, endC;
    struct timeval start, end;
    gettimeofday(&start, NULL);
    startC = clock();

    kNNExecutions++;

//        cout << "----" << k << endl;
    vector<res_point*>* res = grid->getkNN(x, y, k);
//        cout << "size = " << res->size() << endl;


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


bool GPOs::loadLocations(const char* fileName){
  ifstream fin(fileName);


  int location_blacklist[] = { 9246,672270,671606,672473,671698,9247,671886,672328,681993,9313 };
  set<int> blacklist_set (location_blacklist, location_blacklist+10);

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

    if(blacklist_set.find(lid) == blacklist_set.end()){
      loadPoint(x, y, lid, uid, dtm, count);
      count ++;
    }

    if(count%100000==0)
      cout << count << endl;

  }

  fin.close();


  cout << "------- LOADED GPOS from FILE --------- " << endl;
  cout << "Done! Number of checkins: " <<  count << endl;
  cout << "Number of failed checkins: " <<  grid->num_failed << endl;
  cout << "Number of locations: " <<  locations.size() << endl;
  cout << " -------------------------------------- " << endl;
  return true;
}
/*  iterate through each location
**      put each users frequency in a datavector
**      call entropy function on the datavector to obtain the entropy
**      insert entropy into hashmap
*/

unordered_map<int, double>* GPOs::getLocationEntropy(){
  return &location_to_H;
}

map<int, map<int, vector<pair<int, int> >* >*>* GPOs::getCooccurrenceMatrix(){
  return &cooccurrence_matrix;
}

// map<int, map<int, vector<pair<int, int> >* >*>* GPOs::getInsignificantCooccurrenceMatrix(){
//   return &cooccurrence_matrix_insignificant;
// }

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

  grid->addCheckIn(l);

  ids.push_back(lid);
};

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
// TODO FIX this: NEEDS TO BE DIFFERENCE
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

  // cout << "List of users in range : [ ";
  // for (auto i = users->begin(); i != users->end(); ++i)
  //   cout << *i << ',';
  // cout << " ] " << endl;

  return users;
}

void GPOs::groupLocationsByRange(GPOs* gpos, double radius, bool isOptimistic){
  double radius_geo_dist = (radius/1000) * 360 / EARTH_CIRCUMFERENCE,x=0, y=0;
  unsigned int lid, count=0;
  long unsigned int iterations=0;
  unordered_set<int>* seenLocations = new unordered_set<int>();
  boost::posix_time::ptime time;

  GPOs *_duplicate_gpos = new GPOs(gpos);

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
  cout << "Checkins inserted : " << count << endl;
  cout << "Checkins failed lookup : " << iterations << endl;

  delete seenLocations;
  delete _duplicate_gpos;
  generateFrequencyCache();
}

// Radius in meters Adding gaussian noise
// SEED HAS BEEN SET
void GPOs::loadPurturbedLocations(GPOs* gpos, double radius){
  // map <int,int> noise_histogram;
  unsigned int lid = 0;

  for(auto u = gpos->user_to_location.begin(); u != gpos->user_to_location.end(); u++){
    for(auto loc = u->second->begin(); loc != u->second->end(); loc++){

      if(radius != 0){
        Point *p = (*loc);
        pair<double,double> coordinates_with_noise = util.addGaussianNoise(p->getX(), p->getY(), radius);
        double displacement = util.computeMinimumDistance(p->getX(), p->getY(), coordinates_with_noise.first, coordinates_with_noise.second);
        total_displacement+=displacement;
        loadPoint( coordinates_with_noise.first, coordinates_with_noise.second, lid, u->first, p->getTime(), p->getOrder() );
        lid++;
      } else {
        Point *p = (*loc);
        loadPoint( p->getX(), p->getY(), lid, u->first, p->getTime(), p->getOrder() );
        lid++;
      }

      // int bin = (int) floor(noise_distance);
      // auto it = noise_histogram.find(bin);
      // if(it != noise_histogram.end()){
      //   it->second = it->second + 1;
      // }
      // else{
      //   noise_histogram.insert(make_pair(bin,1));
      // }
    }
  }
  cout<<"Total Displacemnt : "<<(((total_displacement*EARTH_CIRCUMFERENCE) /360)/1000) <<" in km"<<endl;
  cout<<"Average Displacemnt : "<<(((total_displacement *EARTH_CIRCUMFERENCE)/360)/lid)*1000 <<" in meters"<<endl;

  // int inRange=0, total=0;
  // double percentage;
  // cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  // cout << "Noise distribution  : " << endl;
  // for(auto b_it = noise_histogram.begin(); b_it != noise_histogram.end(); b_it++){
  //   cout << "Bin : " << b_it->first << "\t Count : " << b_it->second << endl;
  //   if(abs(b_it->first) <= radius ){
  //     inRange+=b_it->second;
  //   }
  //   total+=b_it->second;
  // }
  // percentage = (double) inRange / (double) total;
  // cout << "In Range : " << percentage << endl;
  // cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

void GPOs::loadPurturbedLocationsBasedOnLocationEntropy(GPOs* gpos, double radius, double limit){
  int lid = LOCATION_NOISE_BOUND;
  for(auto l_it = gpos->location_to_user.begin(); l_it != gpos->location_to_user.end(); l_it++){
    int location = l_it->first;
    vector< Point* > *checkins = l_it->second;

    auto h_it = gpos->location_to_H.find(location);
    double entropy = h_it->second;

    Utilities* util = new Utilities();

    if(entropy > limit){
      for(auto loc_it = checkins->begin(); loc_it != checkins->end(); loc_it++){
        Point *p = (*loc_it);
        pair<double,double> coordinates_with_noise = util->addGaussianNoise(p->getX(), p->getY(), radius);
        loadPoint(coordinates_with_noise.first, coordinates_with_noise.second, lid, p->getUID(), p->getTime(), p->getOrder());
        lid++;
      }
    } else {
      for(auto loc_it = checkins->begin(); loc_it != checkins->end(); loc_it++){
        Point *p = (*loc_it);
        loadPoint(p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder());
      }
    }
  }

  generateFrequencyCache();

  cout << "------- Perturbed checkins based on location entropy --------- " << endl;
  cout << "Done! Number of locations: " <<  locations.size() << endl;
  cout << "Number of locations with added noise: " <<  lid - LOCATION_NOISE_BOUND << endl;
  cout << " -------------------------------------- " << endl;
}

void GPOs::loadPurturbedLocationsBasedOnNodeLocality(GPOs* gpos, map<int, double>* node_locality, double radius, double limit){
  // TODO: Pick a random id
  int lid = LOCATION_NOISE_BOUND;
  for(auto u_it = gpos->user_to_location.begin(); u_it != gpos->user_to_location.end(); u_it++){
    int user_id = u_it->first;
    vector< Point* > *user_checkins = u_it->second;

    auto nl_it = node_locality->find(user_id);
    double locality = nl_it->second;

    Utilities* util = new Utilities();

    // Add noise
    if(locality > limit){
      for(auto loc_it = user_checkins->begin(); loc_it != user_checkins->end(); loc_it++){
        Point *p = (*loc_it);
        pair<double,double> coordinates_with_noise = util->addGaussianNoise(p->getX(), p->getY(), radius);
        loadPoint(coordinates_with_noise.first, coordinates_with_noise.second, lid, p->getUID(), p->getTime(),p->getOrder());
        lid++;
      }
    } else {
      for(auto loc_it = user_checkins->begin(); loc_it != user_checkins->end(); loc_it++){
        Point *p = (*loc_it);
        loadPoint(p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(),p->getOrder());
      }
    }
  }

  generateFrequencyCache();

  cout << "------- Perturbed checkins based on node locality --------- " << endl;
  cout << "Done! Number of locations: " <<  locations.size() << endl;
  cout << "Number of locations with added noise: " <<  lid - LOCATION_NOISE_BOUND << endl;
  cout << " -------------------------------------- " << endl;
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

void GPOs::verifyRange(double radius){
  double radius_geo_dist = (radius/1000) * 360 / EARTH_CIRCUMFERENCE,
         meanGroupSize=0;
  unsigned int places=0;

  vector<res_point*>* checkins;

  for(auto l = locations.begin(); l != locations.end(); l++){
    checkins = grid->getRange((*l)->getX(), (*l)->getY(), radius_geo_dist);
    places += checkins->size();
  };

  cout << "Sum of points in range : " << places << endl;
  meanGroupSize = places / locations.size();
  cout << "Mean points in range : " << " " << meanGroupSize << endl;
}


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

      if(cooccurrence_count>1){
        if(significantly_cooccured_user_pairs.find(make_pair(user_1, user_2)) == significantly_cooccured_user_pairs.end())
          significantly_cooccured_user_pairs.insert(make_pair( user_1, user_2 ));
      } else {
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

      if(u2_it != users_location_frequency_map->end())
        users_location_frequency_map->erase(u2_it);

      if(users_location_frequency_map->size() == 0)
        cooccurrence_matrix.erase(u1_it);
    }
  }

  cout << "Number of users with alteast once cooccurrence : " << cooccurrence_matrix.size() << endl;
  cout << "------ Total Cooccrrences ---------  " << total_cooccurrences << endl;
  cout << "------ User pairs who've co-occurred ---------  " << cooccured_user_pairs.size() << endl;
  cout << "------ User pairs who've co-occurred more than once ---------  " << significantly_cooccured_user_pairs.size() << endl;

  // int count=0;

  // for(auto u=cooccured_user_pairs.begin(); u !=cooccured_user_pairs.end(); u++){
  //   count++;
  //   if(count%100000==0)
  //     cout << count << endl;

  //   int u1Id = u->first;
  //   int u2Id = u->second;

  //   unordered_map<int, int> *u1Locations = users_locations_frequency_map.find(u1Id)->second;
  //   unordered_map<int, int> *u2Locations = users_locations_frequency_map.find(u2Id)->second;

  //   vector<pair<int,int> >* location_freqeuncy_list = new vector<pair<int,int> >();
  //   for(auto l = u1Locations->begin(); l != u1Locations->end(); l++){
  //     auto u2Match = u2Locations->find(l->first);
  //     if( u2Match != u2Locations->end()){
  //       int coocc_count = min(u2Match->second, l->second);
  //       location_freqeuncy_list->push_back(make_pair(l->first, coocc_count));
  //     }
  //   }

  //   if(location_freqeuncy_list->size() > 0){
  //     auto coV_it = cooccurrence_matrix.find(u1Id);

  //     if(coV_it != cooccurrence_matrix.end()){
  //       map<int, vector<pair<int, int> >*>* inner_map = coV_it->second;
  //       inner_map->insert(make_pair(u2Id,location_freqeuncy_list));
  //     }else{
  //       map<int, vector<pair<int, int> >* >* inner_map = new map<int, vector<pair<int, int> >* >();
  //       inner_map->insert(make_pair(u2Id,location_freqeuncy_list));
  //       cooccurrence_matrix.insert(make_pair(u1Id,inner_map));
  //     }
  //   }
  // }
}


map< int, map<int,int>* >* GPOs::getL2U2COOCC(){
  return location_to_user_to_cooccurrences;
}

void GPOs::clearNextNN(){

    //	cout << "clearNextNN nextNNList size = " << nextNNList->size() << endl;

    while(!nextNNList->empty()) {
        delete nextNNList->back();
        nextNNList->pop_back();
        //		objects--;
    }

    //	cout << "clearNextNN nextNNList size = " << nextNNList->size() << " objects = " << objects << endl;

    delete nextNNList;

    //        nextNNList->clear();
    //        delete nextNNList;
    nextNNList = new vector<res_point*>();
    computedNN = returnedNN = finalNextNN = objects = 0;
    flagNextNN = true;

	pureNNexec = 0;
}

double GPOs::estimateNearestDistance(double x, double y, int k){
    return grid->estimateNearestDistance(x,y,k);
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


void GPOs::printCooccurrenceMatrix(){

  ofstream outfile;
  outfile.open("cooccurrence_matrix.csv");

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


void GPOs::loadPurturbedLocationsBasedOnCombinationFunction(GPOs* gpos, map< int, map<int, pair<int,double> >* >* user_to_order_to_location_locality , double radius, bool isGaussainNoise, int type){

  user_to_order_to_location_displacment = new map< int, map<int, pair<int,double> >* >();
  unordered_map<int, double>* HiL_map = gpos->getHiLasMap();
  unordered_map<int, double>* HiJ_map = gpos->getHiJasMap();
  // map< int, map<int, pair<int,double> >* >* user_to_order_to_location_locality = spos->getCheckinLocalityMap();

  int lid   = LOCATION_NOISE_BOUND;
  int new_order = 0;
  int purturbed_count = 0;

  ofstream output_file;
  output_file.open("displacement-combination-function.csv");

  cout << "loadPurturbedLocationsBasedOnCombinationFunction : Running type " << type << "\tIsGaussian : " << isGaussainNoise << endl;

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

      double noise;
      // if(type == 0){
      //   noise = 0.5 * ( expHil + expHiJ ) * radius;
      // } else if (type == 1) {
      //   noise = 0.5 * ( (expHil / HIL_SCALE) + (expHiJ / HIJ_SCALE) ) * radius;
      // }else if(type == 2){
      //   noise = checkin_locality_value * (expHil / HIL_SCALE) * radius;
      // } else{
      //   noise = 0.5 * checkin_locality_value * ( (expHil / HIL_SCALE) + (expHiJ / HIJ_SCALE) ) * radius;
      // }
      noise = 0.5 * ( expHil + expHiJ );
      noise = noise * expHL * checkin_locality_value;

      // Offset
      if(noise != 0)
        noise = noise + 0.25;

      noise = noise * radius;

      pair<double,double> coordinates_with_noise;
      if(isGaussainNoise){
        coordinates_with_noise = util.addGaussianNoise(p->getX(), p->getY(), noise);
      } else {
        coordinates_with_noise = util.addNoise(p->getX(), p->getY(), noise);
      }

      double displacement = util.computeMinimumDistance(p->getX(), p->getY(), coordinates_with_noise.first, coordinates_with_noise.second);
      total_displacement+=displacement;

      output_file << user_id << "\t" << p->getID() <<"\t" << order <<"\t"<< displacement <<endl;

      if(noise != 0){
        loadPoint(coordinates_with_noise.first, coordinates_with_noise.second, lid, p->getUID(), p->getTime(), p->getOrder());
        purturbed_count++;
      } else {
        loadPoint(p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), p->getOrder());
      };

      new_order++;
      lid++;
    }
  }

  cout<<"Number of checkins purtubed : "<< purturbed_count << endl;
  cout<<"Total Displacemnt : "<<(((total_displacement*EARTH_CIRCUMFERENCE) /360)/1000) <<" in km"<<endl;
  cout<<"Average Displacemnt : "<<(((total_displacement *EARTH_CIRCUMFERENCE)/360)/new_order)*1000 <<" in meters"<<endl;
  output_file.close();
  generateFrequencyCache();
}

// void GPOs::loadPurturbedLocationsBasedOnCombinationFunctionofCOOCC(GPOs* gpos , map< int, map<int,int>* >* _location_to_user_to_cooccurrences , double radius, bool isGaussainNoise, int type){

//   map<int, double>* HiL_map = gpos->getHiLasMap();
//   map<int, double>* HiJ_map = gpos->getHiJasMap();

//   int lid   = LOCATION_NOISE_BOUND;
//   int new_order = 0;
//   int purturbed_count = 0;

//   cout << "loadPurturbedLocationsBasedOnCombinationFunctionofCOOCC : Running type " << type << "\tIsGaussian : " << isGaussainNoise << endl;

//   for(auto u_it = gpos->user_to_location.begin(); u_it != gpos->user_to_location.end(); u_it++){

//     int user_id = u_it->first;
//     vector< Point* > *checkins = u_it->second;

//     double hiL = 0, hiJ = 0;

//     auto hil_it = HiL_map->find(user_id);
//     if(hil_it != HiL_map->end())
//       hiL = hil_it->second;

//     auto hij_it = HiJ_map->find(user_id);
//     if(hij_it != HiJ_map->end())
//       hiL = hij_it->second;

//     for(auto loc_it = checkins->begin(); loc_it != checkins->end(); loc_it++){
//       Point *p = (*loc_it);
//       double user_cooccurrenes = 0;
//       auto iter_outer = _location_to_user_to_cooccurrences->find(p->getID());
//       if(iter_outer!= _location_to_user_to_cooccurrences->end()){
//         map<int,int>* user_to_cooccurrences_map = iter_outer->second;
//         auto iter_inner = user_to_cooccurrences_map->find(user_id);
//         if(iter_inner != user_to_cooccurrences_map->end()){
//           user_cooccurrenes = iter_inner->second;
//         }
//       }

//       double noise;
//       if(type == 0){
//         noise = log( 1 + user_cooccurrenes ) * radius;
//       } else if(type == 1){
//         noise = log( 1 + user_cooccurrenes ) * ( exp(-hiL) ) * radius;
//       } else if(type == 2){
//         noise = 0.5 * log( 1 + user_cooccurrenes ) * ( exp(-hiL) + exp(-hiJ) ) * radius;
//       } else {
//         noise = 0.5 * log( 1 + user_cooccurrenes ) * ( exp(-hiL/HIL_SCALE) + exp(-hiJ/HIJ_SCALE) ) * radius;
//       };

//       pair<double,double> coordinates_with_noise;
//       if(isGaussainNoise){
//         coordinates_with_noise = util.addGaussianNoise(p->getX(), p->getY(), noise);
//       } else {
//         coordinates_with_noise = util.addNoise(p->getX(), p->getY(), noise);
//       }

//       double displacement = util.computeMinimumDistance(p->getX(), p->getY(), coordinates_with_noise.first, coordinates_with_noise.second);
//       total_displacement+=displacement;

//       if(noise != 0){
//         loadPoint(coordinates_with_noise.first, coordinates_with_noise.second, lid, p->getUID(), p->getTime(), new_order);
//         purturbed_count++;
//       } else {
//         loadPoint(p->getX(), p->getY(), p->getID(), p->getUID(), p->getTime(), new_order);
//       };

//       new_order++;
//       lid++;
//     }
//   }

//   cout<<"Number of checkins purtubed : "<< purturbed_count << endl;
//   cout<<"Total Displacemnt : "<<(((total_displacement*EARTH_CIRCUMFERENCE) /360)/1000) <<" in km"<<endl;
//   cout<<"Average Displacemnt : "<<(((total_displacement *EARTH_CIRCUMFERENCE)/360)/new_order)*1000 <<" in meters"<<endl;
//   generateFrequencyCache();
// }
