#include "../headers.h"

SimpleQueries::SimpleQueries(IGPOs *spatialIndex, ISPOs *socialGraph){
    spos = socialGraph;
    gpos = spatialIndex;
}


SimpleQueries::~SimpleQueries(){}

// map<int, vector<my_pair>*> social_strength_matrix;
// map<int, map<int, double>*> diversity_matrix;
// map<int, map<int, double>*> weighted_frequency_matrix;
// map<int, map<int, vector<pair<int, int>>*>*> cooccurrence_matrix;
int SimpleQueries::countCooccurredFriends(){
  int friends = 0;
  auto cooccured_user_pairs = gpos->getCoOccurredUserPairs();
  for(auto p=cooccured_user_pairs->begin(); p!=cooccured_user_pairs->end(); p++){
    int u1=p->first, u2=p->second;
    if(spos->areFriends(u1, u2))
      friends++;
  }
  return friends;
}

void SimpleQueries::checkUtilityStats(const char* fileName, double radius, double noise_distance){
  ifstream fin(fileName);
  double x,y, avg_users;
  int total_users=0, locations_with_users=0;
  int count=0;

  if (!fin) {
    std::cerr << "Cannot open locations of interest file file " << fileName << std::endl;
  }

  vector<int> *u_set;

  while (fin){
    fin >> y >> x;
    u_set = gpos->getUsersInRange(x, y, radius, noise_distance);

    if(u_set->size() > 0)
      locations_with_users++;

    total_users += u_set->size();
    count++;
  }

  avg_users = (double) total_users / (double) count;


  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Number of locations " << count << " | Range " << radius << "m "  << " | Noise Radius " << noise_distance << endl;
  cout << "Locations with users :" << locations_with_users << endl;
  cout << "Locations without users :" << count - locations_with_users << endl;
  cout << "Mean users around a location :" << avg_users << endl;
  cout << "Total users across all locations :" << total_users << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

// Given a set of locations of interest and a range; this utility compares the usersInRange from each location
// between base_gpos and this->gpos
void SimpleQueries::checkUtilityRange(const char* fileName, IGPOs *base_gpos, double radius, double noise_distance){
  ifstream fin(fileName);
  double x,y, precision, recall, avg_precision=0, avg_recall=0;
  int count=-1;

  if (!fin) {
    std::cerr << "Cannot open locations of interest file file " << fileName << std::endl;
  }

  vector<int> *u1_set, *u2_set;

  while (fin){
    fin >> y >> x;
    u1_set = base_gpos->getUsersInRange(x, y, radius, radius-noise_distance);
    u2_set = gpos->getUsersInRange(x, y, radius, radius-noise_distance);

    // cout << x << "\t" << y <<"\t" << u1_set->size() << "\t" <<  u2_set->size() << endl;

    std::vector<int> v_intersection;

    std::set_intersection(u1_set->begin(), u1_set->end(),
                          u2_set->begin(), u2_set->end(),
                          std::back_inserter(v_intersection));

    if(u2_set->size() != 0){
      precision = (double) v_intersection.size() / (double) u2_set->size();

      if(v_intersection.size() != 0)
        recall    = (double) v_intersection.size() / (double) u1_set->size();
      else
        recall    = 0;

      avg_precision += precision;
      avg_recall    += recall;

      count++;
    }
  }

  avg_precision /= count;
  avg_recall    /= count;

  double f1 = 2 * avg_precision * avg_recall / ( avg_precision + avg_recall );

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Utility [ RANGE QUERY ]" << endl;
  cout << "Number of locations " << count << " | Range " << radius << "m " << endl;
  cout << "Precision :" << avg_precision << endl;
  cout << "Recall    :" << avg_recall    << endl;
  cout << "F1 Score  :" << f1            << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

// Given a set of locations of interest and a range  this utility counts the pairs of users with
// katz score greater than the defined treshold
void SimpleQueries::checkUtilityProximity(const char* fileName, IGPOs *base_gpos, double radius, double tresh, double noise_distance){
  cout << "Computing user proximity list for base gpos " << endl;
  // Before adding noise
  vector< unordered_set< pair<int,int>, PairHasher >* >* base_proximity_list = SimpleQueries(base_gpos, spos).computeProximityUserList(fileName, radius, tresh, noise_distance);
  cout << "Computed user proximity list for base gpos : " << base_proximity_list->size() << endl;

  // With noise
  cout << "Computing user proximity list for other gpos : " << endl;
  vector< unordered_set< pair<int,int>, PairHasher >* >* cmp_proximity_list  = computeProximityUserList(fileName, radius, tresh, noise_distance);
  cout << "Computed user proximity list for other gpos : " << cmp_proximity_list->size() << endl;

  double precision, recall, avg_precision=0, avg_recall=0;
  int count=base_proximity_list->size();

  for(int i=0; i<count; i++){
    unordered_set< pair<int,int>, PairHasher >* base_proximate_users = base_proximity_list->at(i);
    unordered_set< pair<int,int>, PairHasher >* cmp_proximate_users  = cmp_proximity_list->at(i);

    cout << "Location : " << i << endl;

    int positive = cmp_proximate_users->size(), tp=0, gt=base_proximate_users->size();

    for(auto pair_it=base_proximate_users->begin(); pair_it != base_proximate_users->end(); pair_it++){
      int u1=pair_it->first;
      int u2=pair_it->second;

      // cout << "Checking pair " << u1 << " " << u2 << endl;

      if(cmp_proximate_users->find(make_pair(u1, u2)) != cmp_proximate_users->end())
        tp++;
    }

    if(positive != 0)
      precision = (double) tp / (double) positive;
    else
      precision = 0;

    recall    = (double) tp / (double) gt;

    cout << "\tPrecision : " << precision <<"\tRecall : "<< recall << endl;

    avg_precision += precision;
    avg_recall    += recall;
  }

  avg_precision /= count;
  avg_recall    /= count;

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Utility [ Proximity QUERY ]" << endl;
  cout << "Number of locations " << count << " | Range " << radius << "m " << endl;
  cout << "Precision :" << avg_precision << endl;
  cout << "Recall    :" << avg_recall    << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}


// Threshold defines
vector< unordered_set< pair<int,int>, PairHasher >* >* SimpleQueries::computeProximityUserList(const char* fileName, double radius, double tresh, double noise_distance){
  ifstream fin(fileName);

  double x,y;

  if (!fin) {
    std::cerr << "Cannot open locations of interest file file " << fileName << std::endl;
  }

  vector<int> *u1_set, *u2_set, *u_set;
  int lno=-1;

  vector< unordered_set< pair<int,int>, PairHasher >* >* proximate_users_list = new vector< unordered_set< pair<int,int>, PairHasher >* >();

  while (fin){
    fin >> y >> x;
    lno++;

    cout << "Processing location number: " << lno << endl;

    // u1_set = gpos->getUsersInRange(x, y, radius);
    // u2_set = gpos->getUsersInRange(x, y, radius);
    u_set = gpos->getUsersInRange(x, y, radius, noise_distance);
    u1_set = u_set;
    u2_set = u_set;

    multiset<ranked_pair, ranked_pair_comparator_descending>* proximate_users_set = new multiset<ranked_pair, ranked_pair_comparator_descending>();
    unordered_set< pair<int,int>, PairHasher >* proximate_users = new unordered_set< pair<int,int>, PairHasher >();
    unordered_set< pair<int,int>, PairHasher >* ranked_proximate_users = new unordered_set< pair<int,int>, PairHasher >();

    cout << "number of users around this location : " << u_set->size() << endl;
    if(u1_set->size() > 1){
      for(auto u1_it=u1_set->begin(); u1_it != u1_set->end(); u1_it++){
        for(auto u2_it=u2_set->begin(); u2_it != u2_set->end(); u2_it++){
          int u1id = (*u1_it);
          int u2id = (*u2_it);
          if(u1id != u2id){
            double KatzScore = spos->getKatzScore(u1id, u2id);
            if(proximate_users->find(make_pair(u1id, u2id)) == proximate_users->end() && KatzScore > 0){
              proximate_users->insert(make_pair(u1id, u2id));
              proximate_users_set->insert(ranked_pair(u1id, u2id, KatzScore));
            }
          }
        }
      }
      // cout << "computed ranked pair : " << proximate_users->size() << endl;
      double min_score;
      auto rk_it=proximate_users_set->begin();
      for(int count=0; count < tresh && rk_it != proximate_users_set->end(); count++, rk_it++){
        int u1id = rk_it->getId1();
        int u2id = rk_it->getId2();
        ranked_proximate_users->insert(make_pair(u1id, u2id));
        min_score = rk_it->getScore();
      }
      // cout << "\tMinimum pair score is " << min_score << endl;
      // if(ranked_proximate_users->size() > 0 && min_score > 0.01){
      //   cout << y << " " << x << " "  << min_score << " " << ranked_proximate_users->size() << endl;
      // }
      cout << "Persisting top " << ranked_proximate_users->size() << " friendships" << endl;
    }

    // Delete
    // proximate_users_set
    proximate_users_set->clear();
    delete proximate_users_set;
    // proximate_users
    proximate_users->clear();
    delete proximate_users;
    // cout << "keeping top " << ranked_proximate_users->size() <<" user_paris " << endl;


    proximate_users_list->push_back(ranked_proximate_users);
  }

  return proximate_users_list;
}

map< int, bool >* SimpleQueries::getUsersOfInterest(double tresh){
  map< int, bool >* users_of_interest = new map< int, bool >();

  for (auto s_it = social_strength_matrix.begin(); s_it != social_strength_matrix.end(); s_it++){
    int user_1 = s_it->first;
    auto user_ss_list = s_it->second;

    for(auto ss_it = user_ss_list->begin(); ss_it!= user_ss_list->end();ss_it++){
      if(ss_it->second >= tresh){
        if(users_of_interest->find(user_1) == users_of_interest->end())
          users_of_interest->insert(make_pair(user_1, true));
      }
    }
  }

  return users_of_interest;
}

bool SimpleQueries::areEBMFriends(int source, int target, double tresh){
  if(source > target){
    int temp = target;
    target = source;
    source = temp;
  }

  auto iter = social_strength_matrix.find(source);

  if(iter == social_strength_matrix.end())
    return false;

  auto user_social_strengths_map = iter->second;

  auto user_social_strength_it = user_social_strengths_map->find(target);

  if(user_social_strength_it == user_social_strengths_map->end())
    return false;

  double social_strength = user_social_strength_it->second;

  if(social_strength >= tresh)
    return true;
  else
    return false;
}

bool SimpleQueries::areTrueEBMFriends(int source, int target, double tresh){
  return ( areEBMFriends(source, target, tresh) && spos->areFriends(source, target) );
  // return areEBMFriends(source, target, tresh);
}

void SimpleQueries::verifySocialStrength(double tresh){

  int postitive=0, true_positive=0, gt = countCooccurredFriends(), total_score=0, nFriends=spos->getNumberOfFriends();

  for (auto s_it = social_strength_matrix.begin(); s_it != social_strength_matrix.end(); s_it++){
    int user_1 = s_it->first;
    auto user_ss_list = s_it->second;

    for(auto ss_it = user_ss_list->begin(); ss_it!= user_ss_list->end();ss_it++){
      int user_2 = ss_it->first;

      if(ss_it->second >= tresh){
        if(spos->areFriends(user_1, user_2)){
          true_positive++;
        }
        postitive++;
        total_score += ss_it->second;
      }
    }
  }

  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Number of correct friendships inferred " << true_positive << endl;
  cout << "Number of  friendships inferred " << postitive << endl;
  cout << "Number of friendships with more than one cooccrrences " << gt << endl;
  cout << "Number of friendships " << nFriends << endl;
  double precision = true_positive / (double) postitive;
  double recall    = true_positive / (double) RECALL_BOUND;
  double mean_score    = total_score / (double) postitive;
  double f1 = 2 * precision * recall / ( precision + recall );
  cout << "Precision : " << precision << endl;
  cout << "Recall : " << recall << endl;
  cout << "F1 Score  :" << f1   << endl;
  cout << "Mean Score : " << mean_score << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

void SimpleQueries::buildMatrices(double q){
  unordered_map<int, double>* location_to_H =  gpos->getLocationEntropy();
  map<int, map<int, vector<pair<int, int> >* >*>* cooccurrence_matrix = gpos->getCooccurrenceMatrix();


  for(auto c_it = cooccurrence_matrix->begin(); c_it != cooccurrence_matrix->end(); c_it++){
    int user_1 = c_it->first;
    auto users_location_frequency_map = c_it->second;

    for(auto ulh_it = users_location_frequency_map->begin(); ulh_it != users_location_frequency_map->end(); ulh_it++){
      int user_2 = ulh_it->first;
      vector<pair<int, int>>* cooccurrence_counts_vector = ulh_it->second;

      uint *cooccVector = (uint *) calloc(cooccurrence_counts_vector->size(), sizeof(uint));
      int i=0;
      double weighted_frequency=0;

      for(auto u_it = cooccurrence_counts_vector->begin(); u_it!=cooccurrence_counts_vector->end(); u_it++){
        int location_id = u_it->first;

        double location_entropy;

        auto it_ltH = location_to_H->find(location_id);
        if(it_ltH !=  location_to_H->end()){
          location_entropy = it_ltH->second;
        } else {
          location_entropy = 0;
        }

        double partial_weighted_frequency = u_it->second * exp(-1 * location_entropy);
        weighted_frequency += partial_weighted_frequency;

        cooccVector[i] = u_it->second;
        i++;
      }

      double entropy = calcRenyiEntropyFromCoV(q, cooccVector , cooccurrence_counts_vector->size());
      double diversity = exp(entropy);


      if(user_1 > user_2){
        int temp = user_2;
        user_2 = user_1;
        user_1 = temp;
      }


      // Insert into diversity matrix
      auto dm_it = diversity_matrix.find(user_1);
      if(dm_it == diversity_matrix.end()){
        map<int, double>* dmlist = new map<int, double>();
        dmlist->insert(make_pair(user_2, diversity));
        diversity_matrix.insert(make_pair(user_1, dmlist));
      }else{
        map<int, double>* dmlist = dm_it->second;
        dmlist->insert(make_pair(user_2, diversity));
      }

      // Insert into weighted frequency matrix
      auto wt_it = weighted_frequency_matrix.find(user_1);
      if(wt_it == weighted_frequency_matrix.end()){
        map<int, double>* wtlist = new map<int, double>();
        wtlist->insert(make_pair(user_2, weighted_frequency));
        weighted_frequency_matrix.insert(make_pair(user_1, wtlist));
      }else{
        map<int, double>* wtlist = wt_it->second;
        wtlist->insert(make_pair(user_2, weighted_frequency));
      }
    }
  }
}


map<int, map<int, double>*> SimpleQueries::cacluateSocialStrength(){
    for (auto d_it = diversity_matrix.begin(); d_it != diversity_matrix.end(); d_it++){
        int user_1 = d_it->first;
        auto user_diversity_list = d_it->second;

        auto fm_fit = weighted_frequency_matrix.find(user_1);
        if(fm_fit == weighted_frequency_matrix.end())
            cout<<"ERROR---user not found--------"<<endl;
        auto user_wfreq_list = fm_fit->second;

        for(auto list_it = user_diversity_list->begin(); list_it!= user_diversity_list->end();list_it++){
            int user_2 = list_it->first;
            double diversity_value = list_it->second;
            auto f_it = user_wfreq_list->find(user_2);
            double weighted_frequency_value = f_it->second;
            double social_strength_value = (ALPHA * diversity_value) + (BETA * weighted_frequency_value);


            //swap values of user_1 and user_2 to ensure user1 is less than user2
            if(user_1 > user_2){
              int temp = user_2;
              user_2 = user_1;
              user_1 = temp;
            }

            //find in matrix
            //if not found create new list
            // insert list as pair
            //if found
            //insert into existing list
            auto ss_it = social_strength_matrix.find(user_1);
            if(ss_it == social_strength_matrix.end()){
                map<int, double> *ssmap = new map<int, double>();
                ssmap->insert(make_pair(user_2,social_strength_value));
                social_strength_matrix.insert(make_pair(user_1, ssmap));
            }else{
                map<int, double> *ssmap = ss_it->second;
                ssmap->insert(make_pair(user_2,social_strength_value));
            }
        }
    }
    return social_strength_matrix;
}

void SimpleQueries::cacluateCooccurrenceDistributionBasedOnLocationEntropy(){
  cout << "Computing correlation between location_entropy and cooccurrences : " << endl;

  unordered_map<int, double>* location_to_H = gpos->getLocationEntropy();
  map<int , vector< Point* >*>* location_to_user = gpos->getLocationToUser();
  double diff=0.5;

  for(double i=0; i <= 7.5; i = i + diff){
    int cooccurrence_count = 0;
    for(auto it=location_to_H->begin(); it != location_to_H->end(); it++){
      int location = it->first;
      double entropy = it->second;
      if(entropy >= i && entropy < i+diff){
        auto l_it = location_to_user->find(location);
        vector <Point*>* checkins = l_it->second;
        for(auto checkin=checkins->begin(); checkin != checkins->end(); checkin++){
          Point *p = *checkin;
          cooccurrence_count += gpos->getUserCooccurrences(p->getUID());
        }
      }
    }
    cout << i << " " << cooccurrence_count << endl;
  }
}

void SimpleQueries::cacluateCooccurrenceDistributionBasedOnNodeLocality(){
  cout << "Computing correlation between node_locality and cooccurrences : " << endl;

  map< int, double >* node_locality = spos->getNodeLocality();
  double diff = 0.1;

  for(double i=0; i <= 1; i = i + diff){
    int cooccurrence_count = 0;
    for(auto it=node_locality->begin(); it != node_locality->end(); it++){
      int user = it->first;
      double locality = it->second;
      if(locality >= i && locality < i + diff){
        cooccurrence_count += gpos->getUserCooccurrences(user);
      }
    }
    cout << i << " " << cooccurrence_count << endl;
  }
}

void SimpleQueries::cacluateCooccurrenceDistribution(vector <int> *users){
  map <int,int> co_occ_histogram;

  for(auto u_it=users->begin(); u_it != users->end(); u_it++){
    int c_user = (*u_it);
    int count = gpos->getUserCooccurrences(c_user);

    auto h_it = co_occ_histogram.find(count);
    int bin = (int) floor(count/10);
    if(h_it != co_occ_histogram.end()){
      h_it->second = h_it->second + 1;
    }
    else{
      co_occ_histogram.insert(make_pair(bin,1));
    }
  }

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Cooccurrence distribution  : " << endl;
  for(auto b_it = co_occ_histogram.begin(); b_it != co_occ_histogram.end(); b_it++){
    cout << "Bin : " << b_it->first << "\t Count : " << b_it->second << endl;
  }
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}


void SimpleQueries::writeHistogramstoFile(double tresh){

  ofstream outfile;
  outfile.open("HiL.csv");
  map<int, double>* HiL_map = gpos->getHiLasMap();
  for(auto it = HiL_map->begin(); it !=HiL_map->end(); it++){
    int user_id = it->first;
    double H = it->second;
    outfile<< std::fixed << setprecision(10) << user_id << " "<<H<<endl;
  }
  outfile.close();

  cout<<"Printing HiL complete. size:"<<HiL_map->size()<<endl;

  outfile.open("HiJ.csv");
  map<int, double>* HiJ_map =gpos->getHiJasMap();
  for(auto it = HiJ_map->begin(); it !=HiJ_map->end(); it++){
    int user_id = it->first;
    double H = it->second;
    outfile<< std::fixed << setprecision(10) << user_id << " "<<H<<"\n";
  }
  outfile.close();

  cout<<"Printing HiJ complete.size:"<<HiJ_map->size()<<endl;

  outfile.open("HlL.csv");
  map<int, double>* HlL_map =gpos->getHlLasMap();
  for(auto it = HlL_map->begin(); it !=HlL_map->end(); it++){
    int location_id = it->first;
    double H = it->second;
    outfile<< std::fixed << setprecision(10) << location_id << " "<<H<<"\n";
  }
  outfile.close();

  cout<<"Printing HlL complete. size:"<<HlL_map->size()<<endl;


  //-------------------------------------
  //-- for a given user how many friends ebm truly inferred
  //-- this is not the true detected freinds
  map<int,int> user_to_trueEBMinferred;
  for (auto s_it = social_strength_matrix.begin(); s_it != social_strength_matrix.end(); s_it++){
    int user_1 = s_it->first;
    auto user_ss_list = s_it->second;

    for(auto ss_it = user_ss_list->begin(); ss_it!= user_ss_list->end();ss_it++){
      int user_2 = ss_it->first;

      if(ss_it->second >= tresh){
        if(spos->areFriends(user_1, user_2)){

          auto it = user_to_trueEBMinferred.find(user_1);
          if(it!=user_to_trueEBMinferred.end()){
            it->second = it->second+1;
          }else{
            user_to_trueEBMinferred.insert(make_pair(user_1,1));
          }

          it = user_to_trueEBMinferred.find(user_2);
          if(it!=user_to_trueEBMinferred.end()){
            it->second = it->second+1;
          }else{
            user_to_trueEBMinferred.insert(make_pair(user_2,1));
          }

        }
      }
    }
  }

  outfile.open("user_to_trueEBMinferred.csv");
  for(auto it = user_to_trueEBMinferred.begin(); it !=user_to_trueEBMinferred.end(); it++){
    int user_id = it->first;
    int frequency = it->second;
    outfile << std::fixed << std::setprecision(10) << user_id << " "<<frequency<<"\n";
  }
  outfile.close();

  cout<<"Printing user_to_trueEBMinferred complete. size: "<<user_to_trueEBMinferred.size()<<endl;


  //number of freinds ebm inferred that may not be true
  {
      map<int,int> user_to_EBMinferred;
      for (auto s_it = social_strength_matrix.begin(); s_it != social_strength_matrix.end(); s_it++){
        int user_1 = s_it->first;
        auto user_ss_list = s_it->second;

        for(auto ss_it = user_ss_list->begin(); ss_it!= user_ss_list->end();ss_it++){
          int user_2 = ss_it->first;

          if(ss_it->second >= tresh){

              auto it = user_to_EBMinferred.find(user_1);
              if(it!=user_to_EBMinferred.end()){
                it->second = it->second+1;
              }else{
                user_to_EBMinferred.insert(make_pair(user_1,1));
              }

              it = user_to_EBMinferred.find(user_2);
              if(it!=user_to_EBMinferred.end()){
                it->second = it->second+1;
              }else{
                user_to_EBMinferred.insert(make_pair(user_2,1));
              }
          }
        }
      }

      outfile.open("user_to_EBMinferred.csv");
      for(auto it = user_to_EBMinferred.begin(); it !=user_to_EBMinferred.end(); it++){
        int user_id = it->first;
        int frequency = it->second;
        outfile << std::fixed << std::setprecision(10) << user_id << " "<<frequency<<"\n";
      }
      outfile.close();

      cout<<"Printing user_to_EBMinferred complete. size: "<<user_to_EBMinferred.size()<<endl;
  }


  {
    //-----------------------------------
    //-- locations to the number of true ebm freinds inferred at the location
    map<int,int> location_to_EMBinferred;
    auto cooccurrence_matrix = gpos->getCooccurrenceMatrix();
    for(auto it = cooccurrence_matrix->begin(); it!= cooccurrence_matrix->end(); it++){

      int user_id_1 = it->first;
      auto map_of_vectors = it->second;
      for(auto it_map = map_of_vectors->begin(); it_map != map_of_vectors->end(); it_map++){
        int user_id_2 = it_map->first;
        auto vector_of_locations = it_map->second;
        bool arefriend = areEBMFriends(user_id_1, user_id_2, tresh);
        if(arefriend){
          for(auto it_vector = vector_of_locations->begin(); it_vector!= vector_of_locations->end();it_vector++){
            // number_of_cooccurrences += it_vector->second;

            int location_id = it_vector->first;
            auto lc_it = location_to_EMBinferred.find(location_id);
            if(lc_it != location_to_EMBinferred.end()){
              lc_it->second = lc_it->second + 1;
            }else{
              location_to_EMBinferred.insert(make_pair(location_id,1));
            }
          }
        }
      }
    }

    outfile.open("location_to_EMBinferred.csv");
    for(auto it = location_to_EMBinferred.begin(); it !=location_to_EMBinferred.end(); it++){
      int user_id = it->first;
      int frequency = it->second;
      outfile<< std::fixed << setprecision(10) << user_id << " "<<frequency<<"\n";
    }
    outfile.close();

    cout<<"Printing location_to_EMBinferred complete. size: "<<location_to_EMBinferred.size()<<endl;

  }
  {
    //-----------------------------------
    //-- locations to the number of true ebm freinds inferred at the location
    map<int,int> location_to_trueEBMinferred;
    auto cooccurrence_matrix = gpos->getCooccurrenceMatrix();
    for(auto it = cooccurrence_matrix->begin(); it!= cooccurrence_matrix->end(); it++){

      int user_id_1 = it->first;
      auto map_of_vectors = it->second;
      for(auto it_map = map_of_vectors->begin(); it_map != map_of_vectors->end(); it_map++){
        int user_id_2 = it_map->first;
        auto vector_of_locations = it_map->second;
        bool arefriend = areTrueEBMFriends(user_id_1, user_id_2, tresh);
        if(arefriend){
          for(auto it_vector = vector_of_locations->begin(); it_vector!= vector_of_locations->end();it_vector++){
            // number_of_cooccurrences += it_vector->second;

            int location_id = it_vector->first;
            auto lc_it = location_to_trueEBMinferred.find(location_id);
            if(lc_it != location_to_trueEBMinferred.end()){
              lc_it->second = lc_it->second + 1;
            }else{
              location_to_trueEBMinferred.insert(make_pair(location_id,1));
            }
          }
        }
      }
    }

    outfile.open("location_to_trueEBMinferred.csv");
    for(auto it = location_to_trueEBMinferred.begin(); it !=location_to_trueEBMinferred.end(); it++){
      int user_id = it->first;
      int frequency = it->second;
      outfile<< std::fixed << setprecision(10) << user_id << " "<<frequency<<"\n";
    }
    outfile.close();
    cout<<"Printing location_to_trueEBMinferred complete. size: "<<location_to_trueEBMinferred.size()<<endl;
  }


  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}