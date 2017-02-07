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


// Given a set of locations of interest and a range; this utility compares the usersInRange from each location
// between base_gpos and this->gpos
void SimpleQueries::checkUtilityRange(const char* fileName, IGPOs *base_gpos, double radius){
  ifstream fin(fileName);
  double x,y, precision, recall, avg_precision=0, avg_recall=0;
  int count=-1;

  if (!fin) {
    std::cerr << "Cannot open locations of interest file file " << fileName << std::endl;
  }

  vector<int> *u1_set, *u2_set;

  while (fin){
    fin >> y >> x;
    u1_set = base_gpos->getUsersInRange(x, y, radius);
    u2_set = gpos->getUsersInRange(x, y, radius);

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


  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Utility [ RANGE QUERY ]" << endl;
  cout << "Number of locations " << count << " | Range " << radius << "m " << endl;
  cout << "Precision :" << avg_precision << endl;
  cout << "Recall    :" << avg_recall    << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

// Given a set of locations of interest and a range  this utility counts the pairs of users with
// katz score greater than the defined treshold
void SimpleQueries::checkUtilityProximity(const char* fileName, IGPOs *base_gpos, double radius, double tresh){
  cout << "Computing user proximity list for base gpos : " << endl;
  // Before adding noise
  vector< unordered_set< pair<int,int>, PairHasher >* >* base_proximity_list = SimpleQueries(base_gpos, spos).computeProximityUserList(fileName, radius, tresh);
  // With noise

  cout << "Computing user proximity list for other gpos : " << endl;
  vector< unordered_set< pair<int,int>, PairHasher >* >* cmp_proximity_list  = computeProximityUserList(fileName, radius, tresh);

  double precision, recall, avg_precision=0, avg_recall=0;
  int count=base_proximity_list->size();

  for(int i=0; i<count; i++){
    unordered_set< pair<int,int>, PairHasher >* base_proximate_users = base_proximity_list->at(i);
    unordered_set< pair<int,int>, PairHasher >* cmp_proximate_users  = cmp_proximity_list->at(i);

    int positive = cmp_proximate_users->size(), tp=0, gt=base_proximate_users->size();

    for(auto pair=base_proximate_users->begin(); pair != base_proximate_users->end(); base_proximate_users++){
      int u1=pair->first;
      int u2=pair->second;

      if(cmp_proximate_users->find(make_pair(u1, u2)) != cmp_proximate_users->end())
        tp++;
    }

    precision = (double) tp / (double) positive;
    recall    = (double) tp / (double) gt;

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
vector< unordered_set< pair<int,int>, PairHasher >* >* SimpleQueries::computeProximityUserList(const char* fileName, double radius, double tresh){
  ifstream fin(fileName);

  double x,y;

  if (!fin) {
    std::cerr << "Cannot open locations of interest file file " << fileName << std::endl;
  }

  vector<int> *u_set;
  int lno=-1;

  vector< unordered_set< pair<int,int>, PairHasher >* >* proximate_users_list = new vector< unordered_set< pair<int,int>, PairHasher >* >();

  while (fin){
    fin >> y >> x;
    lno++;

    cout << "Processing location number: " << lno << endl;

    u_set = gpos->getUsersInRange(x, y, radius);

    multiset<ranked_pair, ranked_pair_comparator_descending>* proximate_users_set = new multiset<ranked_pair, ranked_pair_comparator_descending>();
    unordered_set< pair<int,int>, PairHasher >* proximate_users = new unordered_set< pair<int,int>, PairHasher >();

    cout << "number of users around this location : " << u_set->size() << endl;
    cout << "computing proximate_users_set  for location" << endl;

    for(auto u1_it=u_set->begin(); u1_it != u_set->end(); u1_it++){
      auto u2_it=u1_it;
      u2_it++;
      for(; u2_it != u_set->end(); u2_it++){
        int u1id = (*u1_it);
        int u2id = (*u2_it);
        if(u1id > u2id){
          int temp = u2id;
          u2id = u1id;
          u1id = temp;
        }

        double KatzScore = spos->getKatzScore(u1id, u2id);
        // cout << u1id << " " << u2id << " " << KatzScore << endl;
        if(KatzScore>= tresh && proximate_users->find(make_pair(u1id, u2id)) == proximate_users->end()){
          proximate_users->insert(make_pair(u1id, u2id));
          proximate_users_set->insert(ranked_pair(u1id, u2id, KatzScore));
        }
      }
    }

    unordered_set< pair<int,int>, PairHasher >* ranked_proximate_users = new unordered_set< pair<int,int>, PairHasher >();

    cout << "keeping top " << tresh <<" user_paris " << endl;
    int count=0;
    for(auto rk_it=proximate_users_set->begin(); rk_it != proximate_users_set->end(); rk_it++){
      int u1id = rk_it->getId1();
      int u2id = rk_it->getId2();
      if(u1id > u2id){
        int temp = u2id;
        u2id = u1id;
        u1id = temp;
      }
      ranked_proximate_users->insert(make_pair(u1id, u2id));
      count++;
      if(count > tresh)
        break;
    }
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
      if(ss_it->getScore() >= tresh){
        if(users_of_interest->find(user_1) == users_of_interest->end())
          users_of_interest->insert(make_pair(user_1, true));
      }
    }
  }

  return users_of_interest;
}

void SimpleQueries::verifySocialStrength(double tresh){

  int postitive=0, true_positive=0, gt = countCooccurredFriends(), total_score=0, nFriends=spos->getNumberOfFriends();

  for (auto s_it = social_strength_matrix.begin(); s_it != social_strength_matrix.end(); s_it++){
    int user_1 = s_it->first;
    auto user_ss_list = s_it->second;

    for(auto ss_it = user_ss_list->begin(); ss_it!= user_ss_list->end();ss_it++){
      int user_2 = ss_it->getId();

      if(ss_it->getScore() >= tresh){
        if(spos->areFriends(user_1, user_2)){
          true_positive++;
        }
        postitive++;
        total_score += ss_it->getScore();
      }
    }
  }

  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Number of correct friendships inferred " << true_positive << endl;
  cout << "Number of  friendships inferred " << postitive << endl;
  // cout << "Number of friendships with more than one cooccrrences " << gt << endl;
  cout << "Number of friendships " << nFriends << endl;
  double precision = true_positive / (double) postitive;
  double recall    = true_positive / (double) nFriends;
  double mean_score    = total_score / (double) postitive;
  cout << "Precision : " << precision << endl;
  cout << "Recall : " << recall << endl;
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


map<int, vector<my_pair>*> SimpleQueries::cacluateSocialStrength(){
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
                vector<my_pair>* sslist = new vector<my_pair>();
                sslist->push_back(my_pair(user_2,social_strength_value));
                social_strength_matrix.insert(make_pair(user_1, sslist));
            }else{
                vector<my_pair>* sslist = ss_it->second;
                sslist->push_back(my_pair(user_2,social_strength_value));
            }
        }
    }
    return social_strength_matrix;
}
