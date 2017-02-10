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

double SPOs::getKatzScore(int source, int target){
    return KatzScore::calculateKatzScore(source, target, socialgraph_map, KATZ_ATTENUATION, KATZ_PATH_LENGTH);
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

    cout << "Users loaded from Social Graph " << times << endl;
    cout << "totalFriends = " << (totalFriends/(1024)) << "KB" << endl;
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
	if(it!=socialgraph_map->end())
		return it->second;
	else
		return NULL;
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
    for(auto f_it=friend_checkins->begin(); f_it != friend_checkins->end(); f_it++){
      Point * friend_checkin = (*f_it);
      double distSq = gpos->distanceBetween(source_checkin, friend_checkin);
      if(distSq < closestDistance){
        closestDistance = distSq;
      }
    }
  }

  return closestDistance;
}

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
    if(friend_checkins_it != gpos->user_to_location.end())
      friend_checkins = friend_checkins_it->second;
    else
      friend_checkins = new vector< Point* >();

    double distance = computeDistanceBetweenFriends(source_checkins, friend_checkins);

    distances->push_back(distance);
  }

  return distances;
}


double SPOs::computeMeanDistanceBetweenAllFriends(GPOs* gpos){
  double inf=std::numeric_limits<double>::infinity();
  unsigned int count=0, sum_distance=0, f_counter=0;


  for(auto u_it = socialgraph_map->begin(); u_it != socialgraph_map->end(); u_it++){
    int source = u_it->first;
    unordered_set<int>* friends = u_it->second;

    vector<double>* distances = computeDistancesBetweenUserFriends(gpos, source, friends);

    for(auto d_it=distances->begin(); d_it != distances->end(); d_it++){
      double dist = (*d_it);
      if(dist != inf){
        sum_distance += (unsigned int) dist;
        count++;
      }
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
// double SPOs::computeNodeLocality(GPOs* gpos, int source){
//   int *friends;
//   unsigned int numOfFriends;

//   auto source_checkins_it = gpos->user_to_location->find(source);
//   vector< Point* >* source_checkins = source_checkins_it->second;

//   getFriends(source, &friends, &numOfFriends);

//   double locality_sum=0, node_locality;

//   for(auto f = friends->begin(); f != friends->end(); f++){
//     int fid = (*f);
//     auto friend_checkins_it = gpos->user_to_location->find(fid);
//     vector< Point* >* friend_checkins = friend_checkins_it->second;

//     double closestDistanceSq = std::numeric_limits<double>::infinity();

//     for(auto s_it=source_checkins->begin(); s_it != source_checkins->end(); s_it++){
//       Point * source_checkin = (*s_it);
//       for(auto f_it=friend_checkins->begin(); f_it != friend_checkins->end(); f_it++){
//         Point * friend_checkin = (*f_it);
//         double distSq = gpos->distanceBetween(source_checkin, friend_checkin);
//         if(distSq < closestDistanceSq){
//           closestDistanceSq = distSq;
//         }
//       }
//     }
//     delete (*f);

//     locality_sum += exp( -sqrt(closestDistanceSq) / NODE_LOCALITY_BETA );
//   }
//   delete friends;

//   node_locality = (1/numOfFriends) * locality_sum;

//   return node_locality;
// }
