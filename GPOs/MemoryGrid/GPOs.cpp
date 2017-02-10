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

GPOs::~GPOs (){
  delete grid;
  delete &locations;
}


double GPOs::getTotalCPUTime(){
  return totalCPUTime;
}

double GPOs::getTotalTime(){
  return totalTime;
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

void GPOs::updateCheckin(Point* p){
  it = locations.find(p->getID());

  if(it!=locations.end()){
    grid->updateCheckIn(p, (*it).second->getX(), (*it).second->getY());
    locations.insert(pair<int, Point*>(p->getID(), p));
  }
}


bool GPOs::loadLocations(const char* fileName){
  ifstream fin(fileName);

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

    loadPoint(x, y, lid, uid, dtm);

    count ++;

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

//TODO: VERIFY FOR CORRECTNESS ACCORDING TO FORMAT OF DATA STRUCTURE
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
  return &location_to_H;
}


void GPOs::loadPoint(double x, double y, int lid, int uid, boost::posix_time::ptime time){
  Point* l;

  l = new Point(x, y, lid, uid, time);

  locations.insert(pair<int, Point*>(lid, l));
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
  users->erase( unique( users->begin(), users->end() ), users->end() );

  return users;
}

void GPOs::groupLocationsByRange(GPOs* gpos, double radius, bool isOptimistic){
double radius_geo_dist = (radius/1000) * 360 / EARTH_CIRCUMFERENCE,x=0, y=0;
  unsigned int lid, uid, count=0;
  set<int>* seenLocations = new set<int>();
  boost::posix_time::ptime time;

  for(auto l = gpos->locations.begin(); l != gpos->locations.end(); l++){
    x   = l->second->getX();
    y   = l->second->getY();
    lid = l->second->getID();
    uid = l->second->getUID();
    time = l->second->getTime();

    vector<res_point*>* checkins = gpos->getRange(x, y, radius_geo_dist);
    for(auto c = checkins->begin(); c != checkins->end(); c++){
      if(isOptimistic){
        loadPoint(x, y, lid, (*c)->uid, (*c)->time);
      } else {
        if(seenLocations->find( (*c)->id ) == seenLocations->end()){
          loadPoint(x, y, lid, (*c)->uid, (*c)->time);
          seenLocations->insert( (*c)->id );
        }
      }
      delete (*c);
    }
    delete checkins;

    count++;
    if(count % 100000==0)
      cout << count << " " << endl;
  }
  generateFrequencyCache();
}

// Radius in meters Adding gaussian noise
// SEED HAS BEEN SET
void GPOs::loadPurturbedLocations(GPOs* gpos, double radius){
  boost::mt19937 rng;
  boost::normal_distribution<> nd(0.0, radius/2);

  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_norormal(rng, nd);

  double lat=0,lon=0, dLat=0, dLon=0, nLat=0, nLon=0, noise_distance=0;

  int R = EARTH_RADIUS_IN_KILOMETERS * 1000;

  map <int,int> noise_histogram;

  unsigned int lid = 0;

  for(auto u = gpos->user_to_location.begin(); u != gpos->user_to_location.end(); u++){
    for(auto loc = u->second->begin(); loc != u->second->end(); loc++){
      noise_distance = var_norormal();

      lon = (*loc)->getX();
      lat = (*loc)->getY();

      dLat = noise_distance / R;
      dLon = noise_distance / (R * cos(PI * lat / 180));

      nLat = lat + dLat * ( 1 / DEG_TO_RAD );
      nLon = lon + dLon * ( 1 / DEG_TO_RAD );

      loadPoint(nLon, nLat, lid, u->first, (*loc)->getTime());
      lid++;

      int bin = (int) floor(noise_distance);
      auto it = noise_histogram.find(bin);
      if(it != noise_histogram.end()){
        it->second = it->second + 1;
      }
      else{
        noise_histogram.insert(make_pair(bin,1));
      }
    }
  }

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

void GPOs::loadPurturbedLocationsBasedOnNodeLocality(GPOs* gpos, map<int, double>* node_locality, double radius){
  // TODO: Pick a random id
  int lid = LOCATION_NOISE_BOUND;
  for(auto u_it = gpos->user_to_location.begin(); u_it != gpos->user_to_location.end(); u_it++){
    int user_id = u_it->first;
    vector< Point* > *user_checkins = u_it->second;

    auto nl_it = node_locality->find(user_id);
    double locality = nl_it->second;

    Utilities* util = new Utilities();

    // Add noise
    if(locality > 0.75){
      for(auto loc_it = user_checkins->begin(); loc_it != user_checkins->end(); loc_it++){
        Point *p = (*loc_it);
        pair<double,double> coordinates_with_noise = util->addGaussianNoise(p->getX(), p->getY(), radius);
        loadPoint(coordinates_with_noise.first, coordinates_with_noise.second, lid, user_id, p->getTime());
        lid++;
      }
    } else {
      for(auto loc_it = user_checkins->begin(); loc_it != user_checkins->end(); loc_it++){
        Point *p = (*loc_it);
        loadPoint(p->getX(), p->getY(), p->getID(), user_id, p->getTime());
      }
    }
  }

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
        loadPoint(closest_x, closest_y, location_id, u->first, p->getTime());
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
    checkins = grid->getRange(l->second->getX(), l->second->getY(), radius_geo_dist);
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

  for(auto c_it = cooccurrence_matrix.begin(); c_it != cooccurrence_matrix.end(); c_it++){
    int user_1 = c_it->first;
    auto users_location_frequency_map = c_it->second;

    for(auto ulh_it = users_location_frequency_map->begin(); ulh_it != users_location_frequency_map->end(); ulh_it++){
      int user_2 = ulh_it->first, cooccurrence_count=0;
      vector<pair<int, int>>* cooccurrence_counts_vector = ulh_it->second;

      for(auto l_it=cooccurrence_counts_vector->begin(); l_it != cooccurrence_counts_vector->end(); l_it++){
        int cooccrences_at_l = l_it->second;
        cooccurrence_count += cooccrences_at_l;
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
