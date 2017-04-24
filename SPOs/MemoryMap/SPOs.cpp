#include "../../headersMemory.h"

SPOs::SPOs() {
    areFriendsExecutions = getFriendsExecutions = 0;
    totalCPUTime = totalTime = 0.0;
    isGposSet=false;
}

SPOs::SPOs(GPOs * _gpos) {
    areFriendsExecutions = getFriendsExecutions = 0;
    totalCPUTime = totalTime = 0.0;
    gpos = _gpos;
    isGposSet=true;
}

double SPOs::getTotalCPUTime(){
    return totalCPUTime;
}

double SPOs::getTotalTime(){
    return totalTime;
}

int SPOs::getAreFriendsExecutions(){
    return areFriendsExecutions;
}

int SPOs::getGetFriendsExecutions(){
    return getFriendsExecutions;
}

multiset<my_pair, pair_comparator_descending>* SPOs::getDegreeSet(){
    return degreeSet;
}

map< int, double >* SPOs::getNodeLocality(){
  return node_locality_map;
}


// Memorizing Katz computation
double SPOs::getKatzScore(int source, int target){
  if(source > target){
    int temp = target;
    target = source;
    source = temp;
  }

  auto katz_it = katzCache.find(std::make_pair(source,target));
  double score;

  if(katz_it != katzCache.end()){
    score = katz_it->second;
  } else{
    score = 0;
  }

  return score;
}


void SPOs::loadKatzScoreFromMemory(){
  ifstream fin("katz-score.txt");
  if (!fin){
    cout << "Cannot open katz score file katz-score.csv" << endl;
  }
  int u1, u2, count=0;
  double katz_score;

  while(fin){
    fin >> u1 >> u2 >> katz_score;
    katzCache.insert( make_pair(make_pair(u1, u2), katz_score) );
    count++;
  }
  cout << "Processed " << count << " lines " << endl;
  cout << "Loaded Katz score into memory : " << katzCache.size() << endl;
}



void SPOs::precomputeKatzScore(GPOs *gpos, int parts, int part, double dTresh){
  ofstream output_file;

  ostringstream ss;
  ss << "katz-score-" << part << ".txt";

  output_file.open(ss.str());

  int start  = socialgraph_map->size() * (part - 1) / parts;
  int offset = socialgraph_map->size() * part       / parts;

  cout << "Number of users: " << socialgraph_map->size() << endl;
  cout << "Start:  " << start << endl;
  cout << "Offset: " << offset << endl;

  auto u1_it = socialgraph_map->begin();
  std::advance(u1_it, start);

  // int pairs = 0;
  // struct timeval startTime, endTime;
  // gettimeofday(&startTime, NULL);

  for(int i = start; i < offset; i++, u1_it++){
    int source = u1_it->first;
    cout<<"Trying source: "<< source<<endl;

    vector<int>* relavent_users = gpos->getUsersInRange(source, dTresh);
    cout <<"Relevant user list size :"<< relavent_users->size() << endl;


    for(auto u2_it = relavent_users->begin(); u2_it != relavent_users->end(); u2_it++){
      int target = (*u2_it);
      if(source != target){
        if(source < target){
           double score = KatzScore::calculateKatzScore(source, target, socialgraph_map, KATZ_ATTENUATION, KATZ_PATH_LENGTH);
           // pairs++;
           // if(pairs % 1000 == 0){
           //   gettimeofday(&endTime, NULL);
           //   cout << "Number of pairs" << pairs << "\tTime: "<< util.print_time(startTime, endTime)/1000 << endl;
           // }
           if(score > 0){
              output_file << std::fixed << std::setprecision(5) << source << " " << target << " " << score << endl;
              // cout<< source << " " << target << " " << score << endl;
           }
        }else
          break;
      }
    }

    relavent_users->clear();
    delete relavent_users;
  }
  cout << "------- Wrote katz scores to file " << endl;
  output_file.close();
}

void SPOs::writeUserFriendsToFile(char * DATASET_PATH){
  ofstream output_file;

  stringstream ss;
  ss << DATASET_PATH << "user-friendships.csv";
  const std::string filePath = ss.str();
  output_file.open( filePath.c_str() );

  for(auto it=socialgraph_map->begin(); it != socialgraph_map->end(); it++){
    int user = it->first;
    int n_friends = it->second->size();

    if(n_friends > 0)
      output_file << user << " " << n_friends << endl;
  }

  cout << "Users friends written to user-friendships.csv";

  output_file.close();
}

/*
loads the social graph into a map/hashtable containing
(id, value) pair. The "value" data structure stores the
ids of all friends of a given user.
*/
int SPOs::load(const char* file){

	CH_DEG = new map<int,double>();


    ifstream fin(file);
    if (!fin){
        cout << "Cannot open Social graph file " << file << endl;
        return -1;
    }

	//degreeMap =  new map<int, int>():
	degreeSet = new multiset<my_pair, pair_comparator_descending>();
    int id, size;
    unsigned int times = 0;
    int totalFriends = 0;
    cout << "Loading the Social Graph from " << file << endl;
    Value* entry;
    while(fin){ //NUMOFUSERS

        fin >> id >> size;
        if (! fin.good()){
            cout << "fin is not good: id = " << id << endl;
            continue;
        }

		/*
		auto it = CH_DEG->find(size);
		if(it!= CH_DEG->end()){
		it->second = it->second +1;
		}
		else{
			CH_DEG->insert(make_pair(size,1));
		}
		*/

		// maxS.push_back(size);
		// if(size>maxS)
			// maxS = size;

		//degreeMap->insert(make_pair(id,size));

		unordered_set<int>* friend_set;
		if(size!=0)
			friend_set = new unordered_set<int>();
		else
			friend_set = NULL;

		degreeSet->insert(my_pair(id,((double)size/(double)MAXSC)));
        entry = new Value(size, id);

        int* list = (int*) malloc(sizeof(int)*size);
        totalFriends+=sizeof(int)*size;
        for(int i = 0; i<size; i++){
          int friend_id ;
          fin >> friend_id;

          if(!isGposSet || gpos->significantly_cooccured_user_pairs.find(make_pair(id, friend_id)) != gpos->significantly_cooccured_user_pairs.end()){
  			    friend_set->insert(friend_id);
            list[i] = friend_id;
            edges = edges + 1;
          }

        }
        entry->setList(list, size);

		socialgraph_map->insert(make_pair(id, friend_set));
		// cout<<"Inserting user Id; "<<id<<" and size: ";
		// if(size == 0)
			// cout<<"NULL"<<endl;
		// else
			// cout<<friend_set->size()<<endl;


        //ids.push_back(id);

        hashTable.insert(pair<int, Value*>(id, entry));

        times++;
        if(times%1000000 == 0)
            cout << times << endl;
    }

	//print out degreeSet
	// multiset<my_pair,pair_comparator_descending>::iterator itSOC ;
	// for(itSOC = degreeSet->begin(); itSOC != degreeSet->end(); itSOC++){
		// my_pair user = *itSOC;
		// double local_social_score = user.getScore();
		// int local_user_id = user.getId();
		// cout << "Social Degree Set - > User ID: "<<local_user_id<<" | Score: "<<local_social_score<<endl;
	// }


	// sort(maxS.begin(), maxS.end(), std::greater<int>());
	// for(int i=0;i<10;i++){
		// cout<<"MAXSC at i = "<< i <<" is "<<maxS.at(i)<<endl;
	// }
    fin.clear();
    fin.close();

	/*
	// iterate and cumulate the size of users
	auto prev_it = CH_DEG->begin();
	auto it = CH_DEG->begin(); ++it;
	while(it != CH_DEG->end()){
		it->second = (prev_it->second + it->second)*CONV_PERCENTILE;
		++prev_it;
		++it;
	}
	*/

    cout << "Users loaded from Social Graph File : " << times << endl;
    cout << "TotalFriends : " << (totalFriends/(1024)) << "KB" << endl;
    cout << "Edges : " << edges << endl;
    return 0;
}


//void SPOs::addFriendship(int user1, int user2){
//    // add user1 to the list of user2
//    //Value *user2_list = hashTable2[user2];
//    Value *user2_list = hashTable.find(user2)->second;
//    //	cout << "addFriendship one start user = " << user2 << endl;
//    int * tmp2 = (int*) realloc(user2_list->list, (user2_list->size+1)*sizeof(int));

//    if(tmp2 != NULL){
//        user2_list->list = tmp2;
//        user2_list->size = user2_list->size+1;
//        util.addToSortedList(user2_list->list, user2_list->size, user1);
//    }
//    else
//        cout << "Error on memory allocation" << endl;

//    //	cout << "addFriendship one finished" << endl;

//    // add user2 to the list of user1
//    //Value *user1_list = hashTable2[user1];
//    Value *user1_list = hashTable.find(user1)->second;
//    //	cout << "addFriendship for user" << user1 << endl;
//    int * tmp1 = (int*) realloc(user1_list->list, (user1_list->size+1)*sizeof(int));
//    //	cout << "memory located" << endl;
//    if(tmp1 != NULL){
//        user1_list->list = tmp1;
//        user1_list->size = user1_list->size+1;
//        util.addToSortedList(user1_list->list, user1_list->size, user2);
//    }
//    else
//        cout << "Error on memory allocation" << endl;
//    //	cout << "addFriendship two finished" << endl;
//}


//int SPOs::loadfriends(const char* fileName){

//    FILE *file;
//    int id, size;
//    int times = 0;


//    file = fopen(fileName,"r");
//    cout << "Loading the Social Graph from " << file << endl;
//    if(file==NULL)
//        cout << "Cannot open Social graph file " << file << endl;

//    cout << "Loading the Social Graph from " << fileName << endl;
//    int totalFriends = 0;
//    while(!feof(file)){ //NUMOFUSERS
//        //fin >> id >> size;
//        fscanf(file, "%d %d", &id, &size);

//        Value* entry = new Value(size, id);

//        int* list = (int*) malloc(sizeof(int)*size);
//        totalFriends+=sizeof(int)*size;
//        for(int i = 0; i<size; i++){
//            fscanf(file, " %d", &(list[i]));
//            //fin >> list[i];
//        }
//        //cout << "test : " << list[1] << endl;

//        entry->setList(list, size);

//        hashTable2[id] = entry;
//        //hashTable.insert(pair<int, Value*>(id, entry));
//        //delete list;

//        times++;
//        if(times%100000 == 0)
//            cout << times << endl;
//    }
//    fclose(file);
//    cout << "Done!" << endl;
//    cout << "totalFriends = " << (totalFriends/(1024)) << "KB" << endl;
//    return 0;
//}


/*
returns the size of the list of friends of a given user (id)
as the integer "numOfFriends" along with the list in integer
array "friends"
*/
void SPOs::getFriends(int id, int*& friends, unsigned int &numOfFriends){
    clock_t startC, endC;
    struct timeval start, end;
    gettimeofday(&start, NULL);
    startC = clock();
    getFriendsExecutions++;

    Value* v = NULL;
    map<int, Value*>::iterator it = hashTable.find(id);
    if(it != hashTable.end()){
        //element found;
        v = it->second;
    }
    else {
        numOfFriends = 0;
        friends = (int*)malloc(sizeof(int)*numOfFriends);
    }

    if(v!=NULL){
        numOfFriends = v->getListSize();
        int* tmp = v->getList();
        friends = (int*)malloc(sizeof(int)*numOfFriends);

        for(unsigned int i=0; i< numOfFriends; i++){
            friends[i] = tmp[i];
        }
        free(tmp);
    }

    endC = clock();
    totalCPUTime += (((double)(endC-startC)*1000.0)/(CLOCKS_PER_SEC));
    gettimeofday(&end, NULL);
    totalTime += util.print_time(start, end);
}


unordered_set<int>* SPOs::getFriends(int id){
	auto it = socialgraph_map->find(id);
	if(it!=socialgraph_map->end()){
		return it->second;
  }
	else{
    unordered_set<int>* fset = new unordered_set<int>();
    socialgraph_map->insert(make_pair(id, fset));
		return fset;
  }
}

/*
iterates through the list of friends of user1 to find
user2. Returns boolean variable denoting the same.
*/
// bool SPOs::areFriends(int user1, int user2){
//     clock_t startC, endC;
//     struct timeval start, end;
//     gettimeofday(&start, NULL);
//     startC = clock();

//     areFriendsExecutions++;

//     Value* v = NULL;
//     map<int, Value*>::iterator it = hashTable.find(user1);
//     if(it != hashTable.end()){
//         //element found;
//         v = it->second;
//     }

//     if(v != NULL){
//         int size = v->getListSize();
//         int* f = v->getList();

//         for(int j = 0; j < size; j++){
//             if (f[j] == user2){
//                 endC = clock();
//                 totalCPUTime += (((double)(endC-startC)*1000.0)/(CLOCKS_PER_SEC));
//                 gettimeofday(&end, NULL);
//                 totalTime += util.print_time(start, end);
//                 return true;
//             }
//         }
//         free(f);
//     }

//     endC = clock();
//     totalCPUTime += (((double)(endC-startC)*1000.0)/(CLOCKS_PER_SEC));
//     gettimeofday(&end, NULL);
//     totalTime += util.print_time(start, end);
//     return false;
// }

bool SPOs::areFriends(int user1, int user2){
    auto it = socialgraph_map->find(user1);
    if(it!=socialgraph_map->end()){
        auto friends_set = it->second;

        // cout<<"Friend Set size: "<<friends_set->size()<<" friends: ";
        // for(auto it_fset = friends_set->begin();it_fset!=friends_set->end();it_fset++){
        //     cout<<*it_fset<<" ";
        // }
        // cout<<endl;

        auto it_f = friends_set->find(user2);
        if(it_f != friends_set->end()){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        return false;
    }

}

//double SPOs::computeDensity(Group G){

//}
void SPOs::printTriangles( int id, vector<int> friendList){
	int num = 0;

	cout << "PRINTING TRIANGLES of user "<<id<<endl;
	for (unsigned int i = 0; i< friendList.size();i++){
		for (unsigned int j = i; j< friendList.size();j++){
			if(areFriends(friendList.at(i),friendList.at(j))){
				cout << "{"<<friendList.at(i)<<" , "<<friendList.at(j)<<"}" << " ";
				num ++;
			}
		}
	}
	cout<<"Total triangles = "<<num<<endl;
}

int SPOs::getUserDegree(int id){
	int degree = 0;
	Value* v = NULL;
    map<int, Value*>::iterator it = hashTable.find(id);
    if(it != hashTable.end()){
        //element found;
        v = it->second;
		if(v!=NULL){
			degree = v->getListSize();
		}
    }

	return degree;
}

int SPOs::getNumberOfFriends(){
    return edges;
}

double SPOs::computeDistanceBetweenFriends(vector< Point* >* source_checkins, vector< Point* >* friend_checkins){

  double closestDistance = std::numeric_limits<double>::infinity();

  for(auto s_it=source_checkins->begin(); s_it != source_checkins->end(); s_it++){
    Point * source_checkin = (*s_it);
    auto source_checkin_time = source_checkin->getTime();

    for(auto f_it=friend_checkins->begin(); f_it != friend_checkins->end(); f_it++){
      Point * friend_checkin = (*f_it);
      auto friend_checkin_time = friend_checkin->getTime();

      boost::posix_time::time_duration time_difference = source_checkin_time - friend_checkin_time;
      int time_diff_seconds = abs(time_difference.total_seconds());

      if(time_diff_seconds <= TIME_RANGE_IN_SECONDS){

        double distSq = util.computeMinimumDistance(source_checkin->getX(), source_checkin->getY(), friend_checkin->getX(), friend_checkin->getY());
        if(distSq < closestDistance){
          closestDistance = distSq;
        }
      }

    }
  }

  return closestDistance;
}

/*
vector<double>* SPOs::computeDistancesBetweenUserFriends(GPOs* gpos, int source, unordered_set<int>* friends){
  double inf=std::numeric_limits<double>::infinity();
  vector<double>* distances = new vector<double>();

  vector< Point* >* source_checkins;

  auto source_checkins_it = gpos->user_to_location.find(source);

  // Ensuring check-ins are present
  if(source_checkins_it != gpos->user_to_location.end())
    source_checkins = source_checkins_it->second;
  else
    source_checkins = new vector< Point* >();

  for(auto f_it = friends->begin(); f_it != friends->end(); f_it++){
    int fid = (*f_it);

    auto friend_checkins_it = gpos->user_to_location.find(fid);
    vector< Point* >* friend_checkins;

    // Ensuring check-ins are present
    if(friend_checkins_it != gpos->user_to_location.end())
      friend_checkins = friend_checkins_it->second;
    else
      friend_checkins = new vector< Point* >();

    double distance = computeDistanceBetweenFriends(source_checkins, friend_checkins);

    distances->push_back(distance);
  }

  return distances;
}
*/

vector<double>* SPOs::computeDistancesBetweenUserFriends(GPOs* gpos, int source, unordered_set<int>* friends){
  vector<double>* distances = new vector<double>();

  vector< Point* >* source_checkins;

  auto source_checkins_it = gpos->user_to_location.find(source);

  // Ensuring check-ins are present
  if(source_checkins_it != gpos->user_to_location.end())
    source_checkins = source_checkins_it->second;
  else
    source_checkins = new vector< Point* >();

  for(auto f_it = friends->begin(); f_it != friends->end(); f_it++){
    int fid = (*f_it);

    auto friend_checkins_it = gpos->user_to_location.find(fid);
    vector< Point* >* friend_checkins;

    // Ensuring check-ins are present
    if(friend_checkins_it != gpos->user_to_location.end()){
      friend_checkins = friend_checkins_it->second;
      double distance = computeDistanceBetweenFriends(source_checkins, friend_checkins);
      distances->push_back(distance);
    }
  }

  return distances;
}

double SPOs::computeMeanDistanceBetweenAllFriends(GPOs* gpos){
  unsigned int count=0, sum_distance=0, f_counter=0;


  for(auto u_it = socialgraph_map->begin(); u_it != socialgraph_map->end(); u_it++){
    int source = u_it->first;
    unordered_set<int>* friends = u_it->second;

    vector<double>* distances = computeDistancesBetweenUserFriends(gpos, source, friends);

    for(auto d_it=distances->begin(); d_it != distances->end(); d_it++){
      double dist = (*d_it);
      sum_distance += (unsigned int) dist;
      count++;
    }

    if(f_counter%1000 == 0)
      cout << "Processed Edges : " << count << "\tSum distance : " << sum_distance << endl;

    // Freeing memory
    distances->clear();
    delete distances;

    f_counter++;
  }

  double mean_distance = (double) sum_distance / (double) count;
  return mean_distance;
}


// Distance matters: Geo-social metrics for online social networks
double SPOs::computeNodeLocality(GPOs* gpos, int source){
  auto f_it = socialgraph_map->find(source);

  unordered_set<int>* friends = f_it->second;

  double locality_sum=0, node_locality=0;

  vector<double>* distances = computeDistancesBetweenUserFriends(gpos, source, friends);

  for(auto d_it=distances->begin(); d_it != distances->end(); d_it++){
    double dist = (*d_it);

    locality_sum += exp( -dist / NODE_LOCALITY_BETA );
  }

  if(distances->size() != 0)
    node_locality = 1/(double)distances->size() * locality_sum;
  else
    node_locality = 0;

  return node_locality;
}

map< int, double >* SPOs::computeNodeLocality(GPOs* gpos){
  node_locality_map  = new map< int, double >();

  int count=0;

  cout << "Computing node locality" << endl;

  for(auto u_it = socialgraph_map->begin(); u_it != socialgraph_map->end(); u_it++){
    int source = u_it->first;
    double locality = computeNodeLocality(gpos, source);
    node_locality_map->insert(make_pair(source, locality));

    // cout << source << " " << locality << endl;
    if(count%1000 == 0)
      cout << "User Count : " << count << endl;

    count++;
  }

  cout << "Completed node locality" << endl;

  return node_locality_map;
}

// max radius = infinity => Unbounded spatial difference
// Friends checked in within an 8 hour time difference
pair<int,double> SPOs::computeMinimumdistanceToFriend(GPOs* gpos, Point* point_source, vector< Point* >* friend_checkins, double max_radius){
  double closestDistance = std::numeric_limits<double>::infinity();
  int location_id = -1;
  auto point_source_time = point_source->getTime();

  for(auto f_it=friend_checkins->begin(); f_it != friend_checkins->end(); f_it++){
    Point * friend_checkin = (*f_it);

    auto friend_checkin_time = friend_checkin->getTime();
    boost::posix_time::time_duration time_difference = point_source_time - friend_checkin_time;
    int time_diff_seconds = abs(time_difference.total_seconds());

    // Ranges HARD-CODED !!
    if( time_diff_seconds <= ( 8 * 3600 ) ){
      double distSq = util.computeMinimumDistance(point_source->getX(), point_source->getY(), friend_checkin->getX(), friend_checkin->getY());
      double distInKm = 1000 * distSq * EARTH_CIRCUMFERENCE / 360;

      if( distSq < closestDistance && distInKm <= max_radius){
        closestDistance = distSq;
        location_id = friend_checkin->getID();
      }

    }
  }
  return make_pair(location_id,closestDistance);
}

vector<pair<int,double>>* SPOs::computeDistancesToCheckinFriends(GPOs* gpos, Point* point_source, unordered_set<int>* friends, double max_radius){
  vector<pair<int,double>>* location_distances = new vector<pair<int,double> >();

  for(auto f_it = friends->begin(); f_it != friends->end(); f_it++){
    int fid = (*f_it);

    auto friend_checkins_it = gpos->user_to_location.find(fid);
    vector< Point* >* friend_checkins;

    // Ensuring check-ins are present
    if(friend_checkins_it != gpos->user_to_location.end()){
      friend_checkins = friend_checkins_it->second;
      pair<int,double> location_distance_pair = computeMinimumdistanceToFriend(gpos, point_source, friend_checkins, max_radius);

      // Only of users are within range
      // if( location_distance_pair.second != std::numeric_limits<double>::infinity() )
      location_distances->push_back(location_distance_pair);
    }
  }

  return location_distances;
}

// Distance matters: Geo-social metrics for online social networks
double SPOs::computeCheckinLocality(GPOs* gpos, Point* point_source, unordered_set<int>* friends,double max_radius){
  double locality_sum=0, checkin_locality=0;
  unordered_map<int, double>* location_to_H =  gpos->getLocationEntropy();

  vector<pair<int,double>>* location_distances = computeDistancesToCheckinFriends(gpos, point_source, friends, max_radius);

  for(auto d_it=location_distances->begin(); d_it != location_distances->end(); d_it++){
    int location_id = (*d_it).first;
    double dist = (*d_it).second;
    double location_entropy, partial_sum=0;

    auto it_ltH = location_to_H->find(location_id);
    if(it_ltH !=  location_to_H->end()){
      location_entropy = it_ltH->second;
    } else {
      location_entropy = 0;
    }

    partial_sum = exp( (-dist / (double)NODE_LOCALITY_BETA) + (-location_entropy/(double)7) );

    locality_sum += partial_sum;
  }

  if(location_distances->size() != 0)
    checkin_locality = 1/(double)location_distances->size() * locality_sum;
  else
    checkin_locality = 0;

  location_distances->clear();
  delete(location_distances);

  return checkin_locality;
}


map< int, map<int, pair<int,double>>* >* SPOs::computeCheckinLocalityMap(GPOs* gpos, double max_radius){
  checkin_locality_map  = new map<int, map<int, pair<int,double> >*>();

  int count=0;

  cout << "Computing checkin locality" << endl;

  for(auto it = gpos->user_to_location.begin(); it != gpos->user_to_location.end(); it++){
    int source_user_id = it->first;
    auto location_vector = it->second;
    map<int, pair<int,double> >* order_locality_map = new map<int,pair<int,double> >();

    auto f_it = socialgraph_map->find(source_user_id);

    if(f_it==socialgraph_map->end())
      continue;

    unordered_set<int>* friends = f_it->second;

    // Debug:
    // cout << source_user_id << " " << friends->size() << " " << location_vector->size() << endl;

    for(auto loc = location_vector->begin(); loc != location_vector->end(); loc++){
      Point* p = *loc;
      double locality = computeCheckinLocality(gpos, p, friends, max_radius);
      // cout << p->getUID() << " " << p->getOrder() << " " << locality << endl;
      if(locality != 0)
        order_locality_map->insert(make_pair(p->getOrder(),make_pair(p->getID(),locality)));
    }

    if(order_locality_map->size() != 0)
      checkin_locality_map->insert(make_pair(source_user_id, order_locality_map));

    if(count%1000 == 0)
      cout << "User Count : " << count << endl;

    count++;
  }

  cout << "Completed checkin locality of size: "<<checkin_locality_map->size() << endl;

  return checkin_locality_map;
}

// Unbounded temporal difference
// Friends checked in within an 1Km spatial difference ( max_radius, max_checkins )
map< int, double >* SPOs::computeTemporalLocality(int max_checkins, double max_radius, GPOs *gpos){
  int counter=0;
  temporal_locality_map = new map< int, double >();
  double max_radius_geo_dist = (max_radius/1000) * 360 / EARTH_CIRCUMFERENCE;
  for(auto l_it=gpos->location_to_user.begin(); l_it != gpos->location_to_user.end(); l_it++){
    vector< Point* >* checkins_at_l = l_it->second;
    Point *p_sample = checkins_at_l->front();
    double x = p_sample->getX(), y = p_sample->getY();

    double radius_bound = gpos->estimateNearestDistance( x, y, max_checkins, max_radius_geo_dist);
    vector<res_point*>* checkins_in_city = gpos->getRange( x, y, radius_bound );

    for(auto c_it=checkins_at_l->begin(); c_it != checkins_at_l->end(); c_it++){
      Point *p = (*c_it);
      boost::posix_time::ptime l_time = p->getTime();

      map<int, double> user_dist;

      double locality_sum=0, temporal_locality=0, number_of_friends = getFriends( p->getUID() )->size();

      for(auto c = checkins_in_city->begin(); c != checkins_in_city->end(); c++){
        boost::posix_time::ptime c_time = (*c)->time;
        boost::posix_time::time_duration time_difference = (c_time - l_time);
        double diff = (double) abs(time_difference.total_seconds());

        if( areFriends(p->getUID(), (*c)->uid) ){
          auto ud_it = user_dist.find( (*c)->uid );
          if( ud_it == user_dist.end() ){
            user_dist.insert(make_pair((*c)->uid, diff));
          } else if(ud_it->second > diff){
            ud_it->second = diff;
          }
        }
      }

      for(auto ud_it = user_dist.begin(); ud_it != user_dist.end(); ud_it++){
        locality_sum += exp( (-ud_it->second / (double) ( 8 * 3600 )) );
      }

      if(number_of_friends > 0)
        temporal_locality = locality_sum * 1/(double)number_of_friends;

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

void SPOs::writeCheckinLocalityToFile(char* DATASET_PATH){
  ofstream output_file;

  stringstream ss;
  ss << DATASET_PATH << "checkin-locality.csv";
  const std::string filePath = ss.str();
  output_file.open( filePath.c_str() );

  int counter = 0;
  for (auto u_it = checkin_locality_map->begin(); u_it != checkin_locality_map->end(); u_it++){
    double user_id = u_it->first;
    map<int, pair<int,double> >* order_locality_map = u_it->second;
    for(auto iter =order_locality_map->begin(); iter!=order_locality_map->end(); iter++){
      int order = iter->first;
      auto location_nl_pair =  iter->second;
      output_file << user_id << "\t" << location_nl_pair.first <<"\t" << order <<"\t"<< location_nl_pair.second <<endl;
      counter++;
    }
  }

  cout << "------- Wrote checkin locality to files " << counter << endl;
  output_file.close();
}

void SPOs::writeTemporalLocalityToFile(char* DATASET_PATH){
  ofstream output_file;
  stringstream ss;
  ss << DATASET_PATH << "temporal-locality.csv";
  const std::string filePath = ss.str();
  output_file.open( filePath.c_str() );

  int counter = 0;
  for (auto o_it = temporal_locality_map->begin(); o_it != temporal_locality_map->end(); o_it++){
    int order_id = o_it->first;
    double locality = o_it->second;
    if(locality > 0){
      output_file << order_id << "\t" << locality << endl;
      counter++;
    }
  }
  cout << "------- Wrote checkin locality to files " << counter << endl;
  output_file.close();
}

map< int, map<int, pair<int,double> >* >* SPOs::getCheckinLocalityMap(){
  return checkin_locality_map;
}

void SPOs::writeNodeLocalityToFile(){
  ofstream output_file;
  output_file.open("node-locality.csv");
  for (auto u_it = node_locality_map->begin(); u_it != node_locality_map->end(); u_it++){
    double uid = u_it->first;
    double locality = u_it->second;
    output_file << uid << "\t" << locality <<endl;
  }
  cout << "------- Wrote node locality to files " << node_locality_map->size() << endl;
  output_file.close();
}


map< int, double >* SPOs::loadNodeLocalityFromFile(){
  node_locality_map  = new map< int, double >();
  ifstream fin("node-locality.csv");
  if (!fin){
    cout << "Cannot open node locality file node-locality.csv" << endl;
  }
  double user_id, locality;
  while(fin){
    fin >> user_id >> locality;
    node_locality_map->insert(make_pair(user_id, locality));
  }
  return node_locality_map;
}

map<int, double >* SPOs::loadTemporalLocalityFromFile(char* DATASET_PATH){
  temporal_locality_map = new map< int, double >();
  stringstream ss;
  ss << DATASET_PATH << "temporal-locality.csv";
  const std::string filePath = ss.str();
  ifstream fin(filePath.c_str());
  if (!fin){
    cout << "Cannot open node locality file temporal-locality.csv" << endl;
  }
  double order_id, locality;
  while(fin){
    fin >> order_id >> locality;
    temporal_locality_map->insert(make_pair(order_id, locality));
  }
  return temporal_locality_map;
}

map< int, map<int, pair<int,double> >* >* SPOs::loadCheckinLocalityFromFile(char* DATASET_PATH){
  int count=0;
  map< int, map<int, pair<int,double> >* >* user_to_order_to_location_locality  = new map<int, map<int, pair<int,double> >*>();

  stringstream ss;
  ss << DATASET_PATH << "checkin-locality.csv";
  const std::string filePath = ss.str();
  ifstream fin(filePath.c_str());
  if (!fin){
    cout << "Cannot open checkin locality file checkin-locality.csv" << endl;
  }
  double user_id,location_id,order, locality_value;

  while(fin){
    fin >> user_id >> location_id >> order >> locality_value;
    count++;
    auto it = user_to_order_to_location_locality->find(user_id);
    if(it!=user_to_order_to_location_locality->end()){  //user id found
      map<int, pair<int,double> >* order_locality_map = it->second;
      order_locality_map->insert(make_pair(order,make_pair(location_id,locality_value)));
    }else{  //create new map and add user checkin order record
      map<int, pair<int,double> >* order_locality_map = new map<int, pair<int,double> >();
      order_locality_map->insert(make_pair(order,make_pair(location_id,locality_value)));
      user_to_order_to_location_locality->insert(make_pair(user_id,order_locality_map));
    }

  }
  checkin_locality_map = user_to_order_to_location_locality;

  cout<<"Loaded user_to_order_to_location_locality map of unique users: "<<user_to_order_to_location_locality->size()<<" and total values: "<<count<<endl;
  return user_to_order_to_location_locality;
}
