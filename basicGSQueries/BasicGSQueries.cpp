#include "../headersMemory.h"

SimpleQueries::SimpleQueries(GPOs *spatialIndex, ISPOs *socialGraph){
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


void SimpleQueries::getInterestingQueryPoints(double radius, const char* query_file, int DATA_SET){
  ofstream outfile;
  outfile.open( query_file );
  double x,y;
  int day, count=0;

  set<int> checkins_of_interest, selected_checkins;
  gpos->pickSingleCheckinFromCooccurrences(&checkins_of_interest);
  for(auto c_it = checkins_of_interest.begin(); c_it != checkins_of_interest.end(); c_it++){
    int order = (*c_it);
    Point *p = gpos->checkin_list.find(order)->second;
    vector<int>* user_list = gpos->getCheckinsInRangeByHourBlock(p, radius);
    outfile << p->getX() << " " << p->getY() << " " << p->getTime() << " " << p->getOrder() << " " << user_list->size() << endl;
    delete user_list;
  }

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

void SimpleQueries::checkUtilityStats(const char* fileName, double radius, double noise_distance){
  ifstream fin(fileName);

  double x,y, avg_users;
  int total_users=0, locations_with_users=0;
  int count=0;
  int day;

  if (!fin) {
    std::cerr << "Cannot open locations of interest file file " << fileName << std::endl;
  }

  // while (fin){
  //   fin >> y >> x >> day;
  //   unordered_map< int, vector<int>* >* user_list = gpos->getUsersInRangeByHourBlock(x,y,radius,radius-noise_distance);
  //   vector<int> *u_set = user_list->find(day)->second;
  //   locations_with_users++;
  //   total_users += u_set->size();
  //   count++;
  //   delete user_list;
  // }

  // avg_users = (double) total_users / (double) locations_with_users;

  // cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  // cout << "Location time blocks with users : " << locations_with_users << " | Range " << radius << "m "  << " | Noise Radius " << noise_distance << endl;
  // cout << "Total location time_blocks : " << count << endl;
  // cout << "Mean users around a location at an 1 hour interval : " << avg_users << endl;
  // cout << "Total users across all location time_blocks :" << total_users << endl;
  // cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

void SimpleQueries::checkUtilityBasic(GPOs *base_gpos){
  vector< pair<int,int> >* base_cooccurrences = base_gpos->getCooccurredCheckins();
  vector< pair<int,int> >* purturbed_cooccurrences = gpos->getCooccurredCheckins();

  cout << "Base Co-occurrences      :" << base_cooccurrences->size() << endl;
  cout << "Perturbed Co-occurrences :" << purturbed_cooccurrences->size() << endl;

  unordered_set< pair<int,int>, PairHasher > base_cooccurrences_hash, purturbed_cooccurrences_hash;

  for(auto c_it = base_cooccurrences->begin(); c_it != base_cooccurrences->end(); c_it++){
    int o1 = c_it->first;
    int o2 = c_it->second;
    base_cooccurrences_hash.insert(make_pair(o1, o2));
  }

  for(auto c_it = purturbed_cooccurrences->begin(); c_it != purturbed_cooccurrences->end(); c_it++){
    int o1 = c_it->first;
    int o2 = c_it->second;
    purturbed_cooccurrences_hash.insert(make_pair(o1, o2));
  }

  cout << "Base Co-occurrences      :" << base_cooccurrences_hash.size() << endl;
  cout << "Perturbed Co-occurrences :" << purturbed_cooccurrences_hash.size() << endl;

  {
    int true_positive = 0, gt = base_cooccurrences_hash.size(), positive = 0;
    for(auto c_it = purturbed_cooccurrences_hash.begin(); c_it != purturbed_cooccurrences_hash.end(); c_it++){
      int o1 = c_it->first;
      int o2 = c_it->second;
      if(o1 > o2){
        int temp = o2;
        o2 = o1;
        o1 = temp;
      }
      if(base_cooccurrences_hash.find(make_pair(o1, o2)) != base_cooccurrences_hash.end()){
        true_positive++;
      }
      positive++;
    }
    double precision = (double) true_positive / (double) positive;
    double recall    = (double) true_positive / (double) gt;
    double f1        = 2 * precision * recall / ( precision + recall );
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Utility [ BASIC METRIC ]" << endl;
    cout << "Number of co-occurrences before : " << gt << endl;
    cout << "Number of after                 : " << purturbed_cooccurrences->size() << endl;
    cout << "utility_basic_precision{{" << precision  << "}}" << endl;
    cout << "utility_basic_recall{{" << recall  << "}}" << endl;
    cout << "utility_basic_f1{{" << f1  << "}}" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  }

  // Basic utility : Sparse vs Dense
  map <int, vector<double>* > st_knn;
  stringstream ss;
  ss << "knn-noise-combined-100-25-1200-coocc.csv";
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
        neighbours->push_back(st_distance);
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

  double buckets[] = { 0, 0.53208299999999997, 0.58567060000000004, 0.64515509999999998, 0.70223800000000003, 0.75724099999999994, 0.80847699999999989, 0.8628517, 0.91566639999999999, 0.97168499999999991, 1.0265299999999999, 1.0890230000000001, 1.157, 1.2290700000000001, 1.305666, 1.3864299999999998, 1.47228, 1.5558350000000001, 1.6438799999999998, 1.7360349999999998, 1.8352200000000001, 1.9357689999999994, 2.0341900000000002, 2.13984, 2.2491719999999997, 2.35955, 2.4794160000000001, 2.5962529999999999, 2.7206040000000007, 2.8508929999999997, 2.9753799999999999, 3.1026599999999998, 3.2446519999999999, 3.379346, 3.5142840000000004, 3.6633449999999996, 3.8194480000000008, 3.9711199999999995, 4.1262600000000003, 4.2805239999999998, 4.4452800000000003, 4.6173569999999993, 4.8019539999999994, 4.9884849999999998, 5.19076, 5.3870250000000004, 5.5763959999999999, 5.7756300000000005, 5.9883299999999977, 6.2125680000000001, 6.4404900000000005, 6.6806920000000005, 6.9297880000000003, 7.1815429999999987, 7.4403980000000027, 7.6929750000000006, 7.9469760000000012, 8.2380029999999991, 8.542764, 8.8573359999999948, 9.1903600000000001, 9.5360929999999975, 9.9122579999999996, 10.314669999999998, 10.720680000000003, 11.171049999999999, 11.636740000000005, 12.134170000000005, 12.6488, 13.187069999999999, 13.7812, 14.385189999999998, 15.080300000000003, 15.838469999999997, 16.639820000000004, 17.5017, 18.434159999999999, 19.47711, 20.563219999999987, 21.719450000000005, 23.030799999999999, 24.4726, 26.019639999999967, 27.705359999999999, 29.481639999999985, 31.457100000000001, 33.506659999999975, 35.59911000000001, 37.76133999999999, 40.319940000000003, 43.094999999999999, 46.378520000000002, 50.193020000000004, 55.179540000000017, 61.402519999999875, 68.65625, 77.564479999999975, 89.588520000000003, 111.13499999999999, 153.13430000000071, 695 };
  int buckets_size = (sizeof(buckets)/sizeof(*buckets));

  vector<double> bucket_vector;
  map<int, unordered_set<pair<int,int>, PairHasher>*> bucket_hash;
  for(int i=0; i<buckets_size; i++){
    unordered_set<pair<int,int>, PairHasher>* b_hash = new unordered_set<pair<int,int>, PairHasher>();
    bucket_vector.push_back(buckets[i]);
    bucket_hash.insert(make_pair( i, b_hash ));
  }
  sort(bucket_vector.begin(), bucket_vector.end());
  cout << "STEP 1: Generated bucket bounds to calculate accuracy." << endl;

  vector<int> true_positive_vector, gt_vector, positive_vector;
  for(auto c_it = base_cooccurrences_hash.begin(); c_it != base_cooccurrences_hash.end(); c_it++){
    int o1 = c_it->first;
    int o2 = c_it->second;
    double knn_dist;
    if(st_knn.find(o1) != st_knn.end())
      knn_dist = st_knn.find(o1)->second->at(0);
    else if(st_knn.find(o2) != st_knn.end())
      knn_dist = st_knn.find(o2)->second->at(0);
    else
      knn_dist = buckets[buckets_size-1] - 1;

    auto b_val_it = lower_bound(bucket_vector.begin(), bucket_vector.end(), knn_dist);
    double b_val = (*b_val_it);
    int b_val_pos = (b_val_it - bucket_vector.begin());
    auto bset_it = bucket_hash.find( b_val_pos );
    if(bset_it == bucket_hash.end()){
      cout << "BOUND ERROR : " << b_val << " " << b_val_pos << endl;
      continue;
    }
    unordered_set<pair<int,int>, PairHasher>* b_hash = bset_it->second;
    b_hash->insert(make_pair(o1, o2));
  }
  cout << "STEP 2: Split co-occurrences per bucket." << endl;

  for(auto b_it = bucket_hash.begin(); b_it != bucket_hash.end(); b_it++){
    int bucket = b_it->first;
    unordered_set<pair<int,int>, PairHasher>* co_occurred_checkins = b_it->second;
    unordered_set<pair<int,int>, PairHasher> p_co_occurred_checkins;
    set<int> checkins;

    // cout << "\tSTEP 3: Processing bucket " << bucket << endl;
    for(auto co_it = co_occurred_checkins->begin(); co_it != co_occurred_checkins->end(); co_it++){
      int o1 = co_it->first;
      int o2 = co_it->second;
      // True Positives
      if(gpos->cooccurred_checkins.find(make_pair(o1, o2)) != gpos->cooccurred_checkins.end()){
        if(p_co_occurred_checkins.find(make_pair(o1, o2)) == p_co_occurred_checkins.end()){
          p_co_occurred_checkins.insert(make_pair(o1, o2));
        }
      }

      // False Positives
      auto co_index_it = gpos->cooccurrence_index.find(o1);
      if(co_index_it != gpos->cooccurrence_index.end()){
        unordered_set<int>* co_index_set = co_index_it->second;
        for(auto co_ch_it = co_index_set->begin(); co_ch_it != co_index_set->end(); co_ch_it++){
          int order = o1;
          int other = (*co_ch_it);
          if(other < order){
            int temp = order;
            order = other;
            other = temp;
          }
          if(gpos->cooccurred_checkins.find(make_pair(order, other)) == gpos->cooccurred_checkins.end() ){
            if(p_co_occurred_checkins.find(make_pair(order, other)) == p_co_occurred_checkins.end()){
              p_co_occurred_checkins.insert(make_pair(order, other));
            }
          }
        }
      }

      co_index_it = gpos->cooccurrence_index.find(o2);
      if(co_index_it != gpos->cooccurrence_index.end()){
        unordered_set<int>* co_index_set = co_index_it->second;
        for(auto co_ch_it = co_index_set->begin(); co_ch_it != co_index_set->end(); co_ch_it++){
          int order = o2;
          int other = (*co_ch_it);
          if(other < order){
            int temp = order;
            order = other;
            other = temp;
          }
          if(gpos->cooccurred_checkins.find(make_pair(order, other)) == gpos->cooccurred_checkins.end() ){
            if(p_co_occurred_checkins.find(make_pair(order, other)) == p_co_occurred_checkins.end()){
              p_co_occurred_checkins.insert(make_pair(order, other));
            }
          }
        }
      }
    }

    // cout << "\tSTEP 3: Computed co_occ list for bucket" << endl;
    int true_positive = 0, gt = co_occurred_checkins->size(), positive = p_co_occurred_checkins.size();
    for(auto c_it = co_occurred_checkins->begin(); c_it != co_occurred_checkins->end(); c_it++){
      int o1 = c_it->first;
      int o2 = c_it->second;

      auto f_it = p_co_occurred_checkins.find(make_pair(o1,o2));
      if(f_it != p_co_occurred_checkins.end())
        true_positive++;
    }
    delete co_occurred_checkins;

    true_positive_vector.push_back(true_positive);
    gt_vector.push_back(gt);
    positive_vector.push_back(positive);

    cout << "\tSTEP 3: Computed accuracy" << endl;
  }

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "true_positive_count{{";
  for(auto it = true_positive_vector.begin(); it != true_positive_vector.end(); it++){
    cout << *it << ",";
  }
  cout << "}}" << endl;
  cout << "ground_truth_count{{";
  for(auto it = gt_vector.begin(); it != gt_vector.end(); it++){
    cout << *it << ",";
  }
  cout << "}}" << endl;
  cout << "positive_count{{";
  for(auto it = positive_vector.begin(); it != positive_vector.end(); it++){
    cout << *it << ",";
  }
  cout << "}}" << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  for(auto knn_it=st_knn.begin(); knn_it != st_knn.end(); knn_it++){
    delete knn_it->second;
  }
}

// Given a set of locations of interest and a range; this utility compares the usersInRange from each location
// between base_gpos and this->gpos
void SimpleQueries::checkUtilityRange(const char* fileName, GPOs *base_gpos, double radius){
  ifstream fin(fileName);
  double precision, recall, avg_precision=0, avg_recall=0, total_results=0;
  int count=0, order;

  if (!fin) {
    std::cerr << "Cannot open locations of interest file file " << fileName << std::endl;
  }

  while (fin){
    fin >> order;
    auto p_it = base_gpos->checkin_list.find(order);
    if(p_it == base_gpos->checkin_list.end()){
      cout << "Count not find checking in query file" << endl;
      continue;
    }

    Point *p = p_it->second;
    vector<int> *u1_set, *u2_set;

    u1_set = base_gpos->getCheckinsInRangeByHourBlock(p, radius);
    u2_set = gpos->getCheckinsInRangeByHourBlock(p, radius);

    total_results += u1_set->size();

    cout << "Result size : "<< u1_set->size() << endl;

    std::vector<int> v_intersection;

    std::set_intersection(u1_set->begin(), u1_set->end(),
                          u2_set->begin(), u2_set->end(),
                          std::back_inserter(v_intersection));

    if(u1_set->size() > 0){
      precision = (double) v_intersection.size() / (double) u2_set->size();

      if(v_intersection.size() != 0)
        recall    = (double) v_intersection.size() / (double) u1_set->size();
      else
        recall    = 0;

      avg_precision += precision;
      avg_recall    += recall;

      count++;
    }

    delete u1_set;
    delete u2_set;
  }

  avg_precision /= count;
  avg_recall    /= count;

  double f1 = 2 * avg_precision * avg_recall / ( avg_precision + avg_recall );

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Utility [ RANGE QUERY ]" << endl;
  cout << "Number of locations time blocks" << count << " | Range " << radius << "m" << endl;
  cout << "Average number of results per query" << total_results/count << endl;
  cout << "utility_range_precision{{" << avg_precision <<"}}" << endl;
  cout << "utility_range_recall{{" << avg_recall    <<"}}"<< endl;
  cout << "utility_range_f1{{" << f1            <<"}}"<< endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

// Given a set of locations of interest and a range  this utility counts the pairs of users with
// katz score greater than the defined treshold
void SimpleQueries::checkUtilityProximity(const char* fileName, GPOs *base_gpos, double radius, double tresh, double noise_distance){
  ifstream fin(fileName);

  if (!fin) {
    std::cerr << "Cannot open locations of interest file file " << fileName << std::endl;
  }

  double precision, recall, avg_precision=0, avg_recall=0;
  int count = 0;
  int valid_query_points = 0;

  while (fin){
    double x,y;
    fin >> y >> x;

    count++;

    unordered_set< pair<int,int>, PairHasher >*  base_proximate_users = SimpleQueries(base_gpos, spos).computeProximityUserList(x, y, radius, tresh, noise_distance);

    if(base_proximate_users->size() == tresh){
      valid_query_points++;

      unordered_set< pair<int,int>, PairHasher >*  cmp_proximate_users  = computeProximityUserList(x, y, radius, tresh, noise_distance);

      int positive = cmp_proximate_users->size(), tp=0, gt=base_proximate_users->size();

      for(auto pair_it=base_proximate_users->begin(); pair_it != base_proximate_users->end(); pair_it++){
        int u1=pair_it->first;
        int u2=pair_it->second;

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

      cmp_proximate_users->clear();
      delete cmp_proximate_users;
    }

    base_proximate_users->clear();
    delete base_proximate_users;
  }

  avg_precision /= valid_query_points;
  avg_recall    /= valid_query_points;

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Utility [ Proximity QUERY ]" << endl;
  cout << "Number of locations " << valid_query_points << " | Range " << radius << "m " << " | Noise " << noise_distance << "m " << endl;
  cout << "Number of invalid locations " << (count - valid_query_points) << endl;
  cout << "Precision :" << avg_precision << endl;
  cout << "Recall    :" << avg_recall    << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}


// Threshold defines
unordered_set< pair<int,int>, PairHasher >* SimpleQueries::computeProximityUserList(double x, double y, double radius, double tresh, double noise_distance){
  vector<int> *u1_set, *u2_set, *u_set;

  u_set = gpos->getUsersInRange(x, y, radius, radius-noise_distance);
  u1_set = u_set;
  u2_set = u_set;

  multiset<ranked_pair, ranked_pair_comparator_descending>* proximate_users_set = new multiset<ranked_pair, ranked_pair_comparator_descending>();
  unordered_set< pair<int,int>, PairHasher >* proximate_users = new unordered_set< pair<int,int>, PairHasher >();
  unordered_set< pair<int,int>, PairHasher >* ranked_proximate_users = new unordered_set< pair<int,int>, PairHasher >();

  // cout << "number of users around this location : " << u_set->size() << endl;

  if(u1_set->size() > 1){
    for(auto u1_it=u1_set->begin(); u1_it != u1_set->end(); u1_it++){
      for(auto u2_it=u2_set->begin(); u2_it != u2_set->end(); u2_it++){

        int u1id = (*u1_it);
        int u2id = (*u2_it);

        if(u1_it != u2_it){
          double KatzScore = spos->getKatzScore(u1id, u2id);

          // cout << u1id << " " << u2id << " " << KatzScore << endl;

          if( proximate_users->find( make_pair(u1id, u2id) ) == proximate_users->end() && KatzScore > 0.005 ){
            proximate_users->insert(make_pair(u1id, u2id));
            proximate_users_set->insert(ranked_pair(u1id, u2id, KatzScore));
          }

        }

      }
    }

    // cout << "computed ranked pair : " << proximate_users->size() << endl;
    double min_score=0;
    auto rk_it=proximate_users_set->begin();
    for(int count=0; count < tresh && rk_it != proximate_users_set->end(); count++, rk_it++){
      int u1id = rk_it->getId1();
      int u2id = rk_it->getId2();
      ranked_proximate_users->insert(make_pair(u1id, u2id));
      min_score = rk_it->getScore();
    }

    proximate_users->clear();
    delete proximate_users;

    proximate_users_set->clear();
    delete proximate_users_set;

    // cout << "\tMinimum pair score is " << min_score << endl;

    // if(ranked_proximate_users->size() > 0 && min_score > 0.01){
    //   cout << y << " " << x << " "  << min_score << " " << ranked_proximate_users->size() << endl;
    // }
    // cout << "Persisting top " << ranked_proximate_users->size() << " friendships" << endl;
  }

  // cout << "keeping top " << ranked_proximate_users->size() <<" user_paris " << endl;

  return ranked_proximate_users;
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

double SimpleQueries::getSocialStrength(int source, int target){
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

  return social_strength;
}

bool SimpleQueries::areEBMFriends(int source, int target, double tresh){
  return ( getSocialStrength(source, target) >= tresh );
}

bool SimpleQueries::areTrueEBMFriends(int source, int target, double tresh){
  return ( areEBMFriends(source, target, tresh) && spos->areFriends(source, target) );
}

void SimpleQueries::computeAccuracyOfSocialStrength(double precision_limit){

  double prev_precision = 0, precision, recall, mean_score, f1, i;
  int postitive, true_positive, gt, nFriends;
  double tresh, total_score;

  gt = countCooccurredFriends();

  for(i = 0; i < 100; i = i + 0.05){
    tresh = i;
    postitive=0;
    true_positive=0;
    total_score=0;
    nFriends=spos->getNumberOfFriends();

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

    precision = true_positive / (double) postitive;
    recall    = true_positive / (double) gt;
    mean_score    = total_score / (double) postitive;
    f1 = 2 * precision * recall / ( precision + recall );

    if(prev_precision <= precision_limit && precision >= precision_limit)
      break;

    prev_precision = precision;
  }

  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "final_number_of_correct_friendships_inferred{{" << true_positive << "}}" << endl;
  cout << "final_number_of_friendships_inferred{{" << postitive << "}}" << endl;
  cout << "final_number_of_friendships_with_2_or_more_cooccurrences{{" << gt << "}}" << endl;
  cout << "final_number_of_friendships{{" << nFriends << "}}" << endl;
  cout << "final_threshold{{" << i << "}}" << endl;
  cout << "final_precision{{" << precision << "}}" << endl;
  cout << "final_recall{{" << recall << "}}" << endl;
  cout << "final_f1{{" << f1   << "}}" << endl;
  cout << "final_mean_score{{" << mean_score << "}}" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

void SimpleQueries::countEBMInferredFriendships(double tresh){
  int postitive=0;
  double total_score=0;

  for (auto s_it = social_strength_matrix.begin(); s_it != social_strength_matrix.end(); s_it++){
    auto user_ss_list = s_it->second;
    for(auto ss_it = user_ss_list->begin(); ss_it!= user_ss_list->end();ss_it++){
      if(ss_it->second >= tresh){
        postitive++;
        total_score += ss_it->second;
      }
    }
  }

  double mean_score    = total_score / (double) postitive;

  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "number_of_friendships_inferred{{" << postitive << "}}" << endl;
  cout << "threshold{{" << tresh << "}}" << endl;
  cout << "mean_score{{" << mean_score << "}}" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

void SimpleQueries::generateSocialGraph(char *DATASET_PATH, double tresh){
  ofstream outfile;

  stringstream ss;
  ss << DATASET_PATH << "socialGraph.txt";
  const std::string filePath = ss.str();
  outfile.open( filePath.c_str() );

  for (auto s_it = social_strength_matrix.begin(); s_it != social_strength_matrix.end(); s_it++){
    int user_1 = s_it->first;
    auto user_ss_list = s_it->second;

    vector<int> friends;

    for(auto ss_it = user_ss_list->begin(); ss_it!= user_ss_list->end();ss_it++){
      int user_2 = ss_it->first;
      if(ss_it->second >= tresh){
        friends.push_back(user_2);
      }
    }

    if(friends.size() == 0)
      continue;

    outfile << user_1 << " ";
    outfile << friends.size();

    for(auto u_it = friends.begin(); u_it != friends.end(); ++u_it) {
      outfile << " " << *u_it;
    }

    outfile << endl;
  }

  outfile.close();
}

void SimpleQueries::verifySocialStrength(double tresh){

  int postitive=0, true_positive=0, gt = countCooccurredFriends(), nFriends=spos->getNumberOfFriends();
  double total_score=0;

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

  double precision = true_positive / (double) postitive;
  double recall    = true_positive / (double) gt;
  double mean_score    = total_score / (double) postitive;
  double f1 = 2 * precision * recall / ( precision + recall );

  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "number_of_correct_friendships_inferred{{" << true_positive << "}}" << endl;
  cout << "number_of_friendships_inferred{{" << postitive << "}}" << endl;
  cout << "number_of_friendships_with_2_or_more_cooccurrences{{" << gt << "}}" << endl;
  cout << "number_of_friendships{{" << nFriends << "}}" << endl;
  cout << "threshold{{" << tresh << "}}" << endl;
  cout << "precision{{" << precision << "}}" << endl;
  cout << "recall{{" << recall << "}}" << endl;
  cout << "f1{{" << f1   << "}}" << endl;
  cout << "mean_score{{" << mean_score << "}}" << endl;
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
            cout<<"ERROR---user not found---mismatch in diversity and weighted frequency matrices-----"<<endl;
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

    //Memory Cleanup for map<int, map<int, double>*>
    for(auto it = diversity_matrix.begin(); it!= diversity_matrix.end();it++){
        auto user_diversity_list = it->second;
        user_diversity_list->clear();
        delete(user_diversity_list);
    }
    diversity_matrix.clear();

    for(auto it = weighted_frequency_matrix.begin(); it!= weighted_frequency_matrix.end();it++){
        auto user_diversity_list = it->second;
        user_diversity_list->clear();
        delete(user_diversity_list);
    }
    weighted_frequency_matrix.clear();

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

void SimpleQueries::writeHistogramstoFile(char *DATASET_PATH, double tresh, double time_block, map< int, double >* temoral_locality_map){
  printPartialDiversityAndWeightedFrequencyValues(DATASET_PATH, RENY_Q);
  unordered_map<int, double>* location_to_H =  gpos->getLocationEntropy();

  ofstream outfile;
  cout << "Writing location entropy to file : " << endl;

  {
    ofstream outfile;
    stringstream ss;
    std::string filePath;
    ss << DATASET_PATH << "locationEntropy.txt";
    filePath = ss.str();
    outfile.open( filePath.c_str() );
    for(auto it = location_to_H->begin();it!=location_to_H->end();it++){
      outfile<<it->first<<" "<<it->second<<endl;
    }
    outfile.close();
  }

  {
    ofstream outfile;
    stringstream ss;
    std::string filePath;
    ss << DATASET_PATH << "HiL.csv";
    filePath = ss.str();
    outfile.open( filePath.c_str() );
    unordered_map<int, double>* HiL_map = gpos->getHiLasMap();
    for(auto it = HiL_map->begin(); it !=HiL_map->end(); it++){
      int user_id = it->first;
      double H = it->second;
      outfile<< std::fixed << setprecision(10) << user_id << " "<<H<<endl;
    }
    outfile.close();
    cout<<"Printing HiL complete. size:"<<HiL_map->size()<<endl;
  }

  {
    ofstream outfile;
    stringstream ss;
    std::string filePath;
    ss << DATASET_PATH << "HiJ.csv";
    filePath = ss.str();
    outfile.open( filePath.c_str() );
    unordered_map<int, double>* HiJ_map =gpos->getHiJasMap();
    for(auto it = HiJ_map->begin(); it !=HiJ_map->end(); it++){
      int user_id = it->first;
      double H = it->second;
      outfile<< std::fixed << setprecision(10) << user_id << " "<<H<<"\n";
    }
    outfile.close();
    cout<<"Printing HiJ complete.size:"<<HiJ_map->size()<<endl;
  }

  {
    ofstream outfile;
    stringstream ss;
    std::string filePath;
    ss << DATASET_PATH << "HlL.csv";
    filePath = ss.str();
    outfile.open( filePath.c_str() );
    unordered_map<int, double>* HlL_map =gpos->getHlLasMap();
    for(auto it = HlL_map->begin(); it !=HlL_map->end(); it++){
      int location_id = it->first;
      double H = it->second;
      outfile<< std::fixed << setprecision(10) << location_id << " "<<H<<"\n";
    }
    outfile.close();
    cout<<"Printing HlL complete. size:"<<HlL_map->size()<<endl;
  }


  //-------------------------------------
  //-- for a given user how many friends ebm truly inferred
  //-- this is not the true detected freinds
  {
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

    ofstream outfile;
    stringstream ss;
    ss << DATASET_PATH << "user_to_trueEBMinferred.csv";
    const std::string filePath = ss.str();
    outfile.open( filePath.c_str() );
    for(auto it = user_to_trueEBMinferred.begin(); it !=user_to_trueEBMinferred.end(); it++){
      int user_id = it->first;
      int frequency = it->second;
      outfile << std::fixed << std::setprecision(10) << user_id << " "<<frequency<<"\n";
    }
    outfile.close();

    cout<<"Printing user_to_trueEBMinferred complete. size: "<<user_to_trueEBMinferred.size()<<endl;
  }


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


    ofstream outfile;
    stringstream ss;
    ss << DATASET_PATH << "user_to_EBMinferred.csv";
    const std::string filePath = ss.str();
    outfile.open( filePath.c_str() );
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

    ofstream outfile;
    stringstream ss;
    ss << DATASET_PATH << "location_to_EMBinferred.csv";
    const std::string filePath = ss.str();
    outfile.open( filePath.c_str() );
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


    ofstream outfile;
    stringstream ss;
    ss << DATASET_PATH << "location_to_trueEBMinferred.csv";
    const std::string filePath = ss.str();
    outfile.open( filePath.c_str() );
    for(auto it = location_to_trueEBMinferred.begin(); it !=location_to_trueEBMinferred.end(); it++){
      int user_id = it->first;
      int frequency = it->second;
      outfile<< std::fixed << setprecision(10) << user_id << " "<<frequency<<"\n";
    }
    outfile.close();
    cout<<"Printing location_to_trueEBMinferred complete. size: "<<location_to_trueEBMinferred.size()<<endl;
  }

  //OUT user - location - ebm - truebm - precision
  {
    //-----------------------------------
    //-- user and location to the number of ebm freinds inferred at the location
    map<int, map<int,int> *> location_user_to_EBMinferred;
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
            int location_id = it_vector->first;

            auto lc_it = location_user_to_EBMinferred.find(location_id);
            if(lc_it != location_user_to_EBMinferred.end()){
              map<int,int>* user_toEBMinferred = lc_it->second;

              //check for both locations and update/create map as necessary

              auto iter_u1 = user_toEBMinferred->find(user_id_1);
              if(iter_u1 != user_toEBMinferred->end()){
                iter_u1->second = iter_u1->second+1;
              }else{
                user_toEBMinferred->insert(make_pair(user_id_1 , 1));
              }

              auto iter_u2 = user_toEBMinferred->find(user_id_2);
              if(iter_u2 != user_toEBMinferred->end()){
                iter_u2->second = iter_u2->second+1;
              }else{
                user_toEBMinferred->insert(make_pair(user_id_2 , 1));
              }

            }else{
              map<int,int>* user_toEBMinferred = new map<int,int>();
              user_toEBMinferred->insert(make_pair(user_id_1,1));
              user_toEBMinferred->insert(make_pair(user_id_2,1));
              location_user_to_EBMinferred.insert(make_pair(location_id,user_toEBMinferred));
            }
          }
        }
      }
    }

    ofstream outfile;
    stringstream ss;
    ss << DATASET_PATH << "location_user_to_EBMinferred.csv";
    const std::string filePath = ss.str();
    outfile.open( filePath.c_str() );
    for(auto it = location_user_to_EBMinferred.begin(); it !=location_user_to_EBMinferred.end(); it++){
      int location_id = it->first;
      map<int,int>* user_toEBMinferred = it->second;
      for(auto iter = user_toEBMinferred->begin(); iter!= user_toEBMinferred->end(); iter++){
        int user_id = iter->first;
        int frequency = iter->second;
        outfile<< std::fixed << setprecision(10) << user_id << " " <<location_id<<" "<<frequency<<"\n";
      }
    }
    outfile.close();

    cout<<"Printing location_user_to_EBMinferred complete. size: "<<location_user_to_EBMinferred.size()<<endl;

  }

  //-----------------------------------
  //-- user and location to the number of ebm freinds inferred at the location
  {

    map<int, map<int,int> *> location_user_to_trueEBMinferred;
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
            int location_id = it_vector->first;

            auto lc_it = location_user_to_trueEBMinferred.find(location_id);
            if(lc_it != location_user_to_trueEBMinferred.end()){
              map<int,int>* user_toEBMinferred = lc_it->second;

              //check for both locations and update/create map as necessary

              auto iter_u1 = user_toEBMinferred->find(user_id_1);
              if(iter_u1 != user_toEBMinferred->end()){
                iter_u1->second = iter_u1->second+1;
              }else{
                user_toEBMinferred->insert(make_pair(user_id_1 , 1));
              }

              auto iter_u2 = user_toEBMinferred->find(user_id_2);
              if(iter_u2 != user_toEBMinferred->end()){
                iter_u2->second = iter_u2->second+1;
              }else{
                user_toEBMinferred->insert(make_pair(user_id_2 , 1));
              }

            }else{
              map<int,int>* user_toEBMinferred = new map<int,int>();
              user_toEBMinferred->insert(make_pair(user_id_1,1));
              user_toEBMinferred->insert(make_pair(user_id_2,1));
              location_user_to_trueEBMinferred.insert(make_pair(location_id,user_toEBMinferred));
            }
          }
        }
      }
    }

    ofstream outfile;
    stringstream ss;
    ss << DATASET_PATH << "location_user_to_trueEBMinferred.csv";
    const std::string filePath = ss.str();
    outfile.open( filePath.c_str() );
    for(auto it = location_user_to_trueEBMinferred.begin(); it !=location_user_to_trueEBMinferred.end(); it++){
      int location_id = it->first;
      map<int,int>* user_totrueEBMinferred = it->second;

      for(auto iter = user_totrueEBMinferred->begin(); iter!= user_totrueEBMinferred->end(); iter++){
        int user_id = iter->first;
        int frequency = iter->second;
        outfile<< std::fixed << setprecision(10) << user_id << " " <<location_id<<" "<<frequency<<"\n";
      }


    }
    outfile.close();

    cout<<"Printing location_user_to_trueEBMinferred complete. size: "<<location_user_to_trueEBMinferred.size()<<endl;

  }

  // Location Time Locality
  {
    map<int , vector< Point* >*>* location_to_user = gpos->getLocationToUser();

    ofstream outfile;
    stringstream ss;
    ss << DATASET_PATH << "locationTime_Locality.csv";
    const std::string filePath = ss.str();
    outfile.open( filePath.c_str() );

    for(auto l_it=location_to_user->begin(); l_it != location_to_user->end(); l_it++){
      int location_id = l_it->first;
      vector< Point* >* checkins_at_l = l_it->second;

      for(auto c_it=checkins_at_l->begin(); c_it != checkins_at_l->end(); c_it++){
        Point *p = (*c_it);
        int user_id = p->getUID(), order = p->getOrder();
        boost::posix_time::ptime l_time = p->getTime();
        int l_time_block = p->getTimeBlock( time_block );
        int l_day        = l_time.date().day_of_year() + ( l_time.date().year() - 2000 ) * 365;


        double locality=0;

        auto loc_it = temoral_locality_map->find(order);
        if(loc_it != temoral_locality_map->end())
          locality = loc_it->second;

        outfile<< std::fixed << setprecision(10) << location_id << " " <<  l_day << " " << l_time_block <<" " << user_id << " " <<  locality << endl;
      }
    }

    outfile.close();
  }


  //OUT location - time - ebm - truebm/ebm - precision
  {

    map<int , vector< Point* >*>* location_to_user = gpos->getLocationToUser();

    ofstream outfile;
    stringstream ss;
    ss << DATASET_PATH << "locationTime_EBM.csv";
    const std::string filePath = ss.str();
    outfile.open( filePath.c_str() );

    for(auto l_it=location_to_user->begin(); l_it != location_to_user->end(); l_it++){
      int location_id = l_it->first;
      vector< Point* >* checkins_at_l = l_it->second;

      unordered_map<int, map< int, vector<int>*>* >* day_time_block_ulist_map = new unordered_map<int, map< int, vector<int>*>* >();
      map< int, vector<int>*>* time_block_ulist_map;
      vector<int> *u_list;

      for(auto c_it=checkins_at_l->begin(); c_it != checkins_at_l->end(); c_it++){
        Point *p = (*c_it);
        int user_id = p->getUID();
        boost::posix_time::ptime l_time = p->getTime();
        int l_time_block = p->getTimeBlock( time_block );
        int l_day        = l_time.date().day_of_year() + ( l_time.date().year() - 2000 ) * 365;

        auto dtb_it = day_time_block_ulist_map->find(l_day);
        if(dtb_it == day_time_block_ulist_map->end()){
          time_block_ulist_map = new map< int, vector<int>*>();
          u_list = new vector<int>();
          time_block_ulist_map->insert(make_pair(l_time_block, u_list));
          day_time_block_ulist_map->insert(make_pair(l_day, time_block_ulist_map));
        } else {
          time_block_ulist_map = dtb_it->second;
          auto tb_it = time_block_ulist_map->find(l_time_block);
          if(tb_it == time_block_ulist_map->end()){
            u_list = new vector<int>();
            time_block_ulist_map->insert(make_pair(l_time_block, u_list));
          } else {
            u_list = tb_it->second;
          }
        }
        u_list->push_back(user_id);
      }

      for(auto d_it=day_time_block_ulist_map->begin(); d_it!=day_time_block_ulist_map->end(); d_it++){
        time_block_ulist_map = d_it->second;
        for(auto t_it=time_block_ulist_map->begin(); t_it != time_block_ulist_map->end(); t_it++){
          vector<int> *u_list = t_it->second;
          sort(u_list->begin(), u_list->end());
          u_list->erase( unique( u_list->begin(), u_list->end() ), u_list->end() );
        }
      }

      for(auto d_it=day_time_block_ulist_map->begin(); d_it!=day_time_block_ulist_map->end(); d_it++){
        int l_day = d_it->first;
        time_block_ulist_map = d_it->second;

        for(auto t_it=time_block_ulist_map->begin(); t_it != time_block_ulist_map->end(); t_it++){
          int l_time_block = t_it->first;
          vector<int> *u_list = t_it->second;

          int trueEbmFriends=0,ebmFriends=0, actualFriends=0;

          for(auto u1_it = u_list->begin(); u1_it != u_list->end(); u1_it++){
            for(auto u2_it = u1_it; u2_it != u_list->end(); u2_it++){
              if(u1_it != u2_it){
                int user_1 = (*u1_it);
                int user_2 = (*u2_it);
                if(areTrueEBMFriends(user_1, user_2, tresh))
                  trueEbmFriends++;
                if(areEBMFriends(user_1, user_2, tresh))
                  ebmFriends++;
                if(spos-> areFriends(user_1, user_2))
                  actualFriends++;
              }
            }
          }

          outfile<< std::fixed << setprecision(10) << location_id << " " <<  l_day << " " << l_time_block <<" "<< trueEbmFriends << " " << ebmFriends << " " << actualFriends << endl;
        }
      }

      delete day_time_block_ulist_map;
    }
    outfile.close();
    cout<<"Printing locationTime_EBM complete. " << endl;
  }


  //OUT location - time - u1,u2 - number of co-occ
  {

    map<int , vector< Point* >*>* location_to_user = gpos->getLocationToUser();

    ofstream outfile;
    stringstream ss;
    ss << DATASET_PATH << "locationTime_CoOcc.csv";
    const std::string filePath = ss.str();
    outfile.open( filePath.c_str() );

    for(auto l_it=location_to_user->begin(); l_it != location_to_user->end(); l_it++){
      int location_id = l_it->first;
      vector< Point* >* checkins_at_l = l_it->second;

      unordered_map<int, map< int, vector<int>*>* >* day_time_block_ulist_map = new unordered_map<int, map< int, vector<int>*>* >();
      map< int, vector<int>*>* time_block_ulist_map;
      vector<int> *u_list;

      for(auto c_it=checkins_at_l->begin(); c_it != checkins_at_l->end(); c_it++){
        Point *p = (*c_it);
        int user_id = p->getUID();
        boost::posix_time::ptime l_time = p->getTime();
        int l_time_block = p->getTimeBlock( time_block );
        int l_day        = l_time.date().day_of_year() + ( l_time.date().year() - 2000 ) * 365;

        auto dtb_it = day_time_block_ulist_map->find(l_day);
        if(dtb_it == day_time_block_ulist_map->end()){
          time_block_ulist_map = new map< int, vector<int>*>();
          u_list = new vector<int>();
          time_block_ulist_map->insert(make_pair(l_time_block, u_list));
          day_time_block_ulist_map->insert(make_pair(l_day, time_block_ulist_map));
        } else {
          time_block_ulist_map = dtb_it->second;
          auto tb_it = time_block_ulist_map->find(l_time_block);
          if(tb_it == time_block_ulist_map->end()){
            u_list = new vector<int>();
            time_block_ulist_map->insert(make_pair(l_time_block, u_list));
          } else {
            u_list = tb_it->second;
          }
        }
        u_list->push_back(user_id);
      }

      for(auto d_it=day_time_block_ulist_map->begin(); d_it!=day_time_block_ulist_map->end(); d_it++){
        time_block_ulist_map = d_it->second;
        for(auto t_it=time_block_ulist_map->begin(); t_it != time_block_ulist_map->end(); t_it++){
          vector<int> *u_list = t_it->second;
          sort(u_list->begin(), u_list->end());
          u_list->erase( unique( u_list->begin(), u_list->end() ), u_list->end() );
        }
      }

      for(auto d_it=day_time_block_ulist_map->begin(); d_it!=day_time_block_ulist_map->end(); d_it++){
        int l_day = d_it->first;
        time_block_ulist_map = d_it->second;

        for(auto t_it=time_block_ulist_map->begin(); t_it != time_block_ulist_map->end(); t_it++){
          int l_time_block = t_it->first;
          vector<int> *u_list = t_it->second;

          for(auto u1_it = u_list->begin(); u1_it != u_list->end(); u1_it++){
            for(auto u2_it = u1_it; u2_it != u_list->end(); u2_it++){
              if(u1_it != u2_it){
                int user_1 = (*u1_it);
                int user_2 = (*u2_it);
                outfile<< std::fixed << setprecision(10) << location_id << " " <<  l_day << " " << l_time_block <<" "<< user_1 << " " << user_2 << endl;
              }
            }
          }

        }
      }

      delete day_time_block_ulist_map;
    }
    outfile.close();
    cout<<"Printing locationTime_CoOcc complete. " << endl;
  }

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

void SimpleQueries::printPartialDiversityAndWeightedFrequencyValues(char * DATASET_PATH,double alpha){
  unordered_map<int, double>* location_to_H =  gpos->getLocationEntropy();
  map<int, map<int, vector<pair<int, int> >* >*>* cooccurrence_matrix = gpos->getCooccurrenceMatrix();

  ofstream outfile;

  stringstream ss;
    ss << DATASET_PATH << "partialDWFmatrix";
    const std::string filePath = ss.str();
    outfile.open( filePath.c_str() );

  for(auto c_it = cooccurrence_matrix->begin(); c_it != cooccurrence_matrix->end(); c_it++){
    int user_1 = c_it->first;
    auto users_location_frequency_map = c_it->second;

    for(auto ulh_it = users_location_frequency_map->begin(); ulh_it != users_location_frequency_map->end(); ulh_it++){
      int user_2 = ulh_it->first;
      vector<pair<int, int>>* cooccurrence_counts_vector = ulh_it->second;

      uint *cooccVector = (uint *) calloc(cooccurrence_counts_vector->size(), sizeof(uint));

      vector<double> temporary_entropy_value_vector;

      int i=0;

      for(auto u_it = cooccurrence_counts_vector->begin(); u_it!=cooccurrence_counts_vector->end(); u_it++){
        cooccVector[i] = u_it->second;
        i++;
      }

      double *stateProbs;
      double stateLength = cooccurrence_counts_vector->size();

      stateProbs = (double *) checkedCalloc(stateLength,sizeof(double));
      int sumOfCo = sumState(cooccVector, stateLength);
      for (i = 0; i < stateLength; i++) {
         stateProbs[i] = cooccVector[i]/(double)sumOfCo;
      }

      double power_factor = 1/(double)(1.0-alpha);
      double diversity_from_new = 0;

      for (i = 0; i < stateLength; i++) {
        double tempValue = stateProbs[i];
        if (tempValue > 0) {
          diversity_from_new += pow(tempValue,alpha);
        }
      }
      diversity_from_new = pow(diversity_from_new,power_factor);
      for (i = 0; i < stateLength; i++) {
        double tempValue = stateProbs[i];
        temporary_entropy_value_vector.push_back(diversity_from_new * tempValue);
      }

      i=0;
      for(auto u_it = cooccurrence_counts_vector->begin(); u_it!=cooccurrence_counts_vector->end(); u_it++){
        int location_id = u_it->first;

        double location_entropy;

        auto it_ltH = location_to_H->find(location_id);
        if(it_ltH !=  location_to_H->end()){
          location_entropy = it_ltH->second;
        } else {
          location_entropy = 0;
        }

        double partial_weighted_frequency;
        if(DATA_SET == 2){
          partial_weighted_frequency = u_it->second * exp(-1 * location_entropy / 2.5);
        } else {
          partial_weighted_frequency = u_it->second * exp(-1 * location_entropy);
        }

        outfile<<user_1<<" "<<user_2 <<" "<< location_id<<" "<<temporary_entropy_value_vector.at(i)<< " " <<partial_weighted_frequency<<endl;
        i++;
      }

    }
  }
}
