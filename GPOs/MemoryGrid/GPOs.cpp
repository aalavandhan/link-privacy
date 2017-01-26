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
  grid->deleteEmptyCells();
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

void GPOs::getLocation(int id, double* result){
    clock_t startC, endC;
    struct timeval start, end;
    gettimeofday(&start, NULL);
    startC = clock();
    LocationExecutions++;

    it = locations.find(id);

    if(it!=locations.end()){
        result[0]= (*it).second->getX();
        result[1]= (*it).second->getY();
    }
    else{
        result[0] = -1000;
    }

    endC = clock();
    totalCPUTime += (((double)(endC-startC)*1000.0)/(CLOCKS_PER_SEC));
    gettimeofday(&end, NULL);
    totalTime += util.print_time(start, end);

}

res_point* GPOs::getNextNN(double x, double y, int incrStep){
    clock_t startC, endC;
    struct timeval start, end;
    gettimeofday(&start, NULL);
    startC = clock();


    if(computedNN <= returnedNN && flagNextNN){
        NextNNExecutions++;
        computedNN+=incrStep;
//        cout<<"returneNN="<<returnedNN<<endl;
//        cout << "next-" << computedNN << endl;
        vector<res_point*>* kNN = grid->getkNN(x, y, computedNN);
        int size = kNN->size();
//        cout << "size = " << size << endl;

        for(int i = returnedNN; i < size; i++){
            nextNNList->push_back(util.copy((*kNN)[i]));
        }

        while(!kNN->empty()) {
            delete kNN->back();
            kNN->pop_back();
        }
        delete kNN;

        int newNNsize = nextNNList->size();
//        cout << "newNNsize = " << newNNsize << endl;
//        cout << "computted NN = "<<computedNN<<endl;
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

void GPOs::updateCheckin(Point* p){
  it = locations.find(p->getID());

  if(it!=locations.end()){
    grid->updateCheckIn(p, (*it).second->getX(), (*it).second->getY());
    locations.insert(pair<int, Point*>(p->getID(), p));
  }
}

void GPOs::loadPoint(double x, double y, int lid, int uid, boost::posix_time::ptime time){
  Point* l;

  l = new Point(x, y, lid, uid, time);

  locations.insert(pair<int, Point*>(lid, l));
  auto lh_it = locationHistory.find(uid);

  if( lh_it == locationHistory.end() ){
    vector<Point* >* pts = new vector<Point*>();
    pts->push_back(l);
    locationHistory.insert( make_pair(uid, pts) );
  } else {
    vector<Point* >* pts = lh_it->second;
    pts->push_back(l);
  }

  auto lu_it = location_to_user.find(lid);

  if( lu_it == location_to_user.end() ){
    set<int>* users = new set<int>();
    users->insert(uid);
    location_to_user.insert( make_pair(lid, users) );
  } else {
    set<int>* users = lu_it->second;
    users->insert(uid);
  }

  grid->addCheckIn(l);

  ids.push_back(lid);
};


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
  cout << " -------------------------------------- " << endl;
  return true;
}

// Radius in meters
// SEED HAS BEEN SET
void GPOs::loadPurturbedLocations(GPOs* gpos, double radius){
  boost::mt19937 rng;
  boost::normal_distribution<> nd(0.0, radius);

  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_norormal(rng, nd);

  double lat=0,lon=0, dLat=0, dLon=0, nLat=0, nLon=0, noise_distance=0;

  int R = EARTH_RADIUS_IN_KILOMETERS * 1000;

  unsigned int lid = 0;

  for(auto u = gpos->locationHistory.begin(); u != gpos->locationHistory.end(); u++){
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
    }
  }

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
