#include "headersMemory.h"
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>


#ifndef GLOBALS
#define GLOBALS

int DATA_SET = 0;
double MIN_X = -124.848974;
double MAX_X = -66.885444;
double MIN_Y = 24.396308;
double MAX_Y = 49.384358;
double DATASET_SIZE = 3661488;
char *DATASET_PATH;

double DELTA_X = 0;
double DELTA_Y = 0;
int MAXSC = 0;                // maximum number of friends plus one
double MAXDIST = 50;        // maximum distance between to points
int MAXT = 0;

// EBM CONSTANTS FOR GOWALLA FULL
// double ALPHA = 0.55;
// double GAMMA = 0;
// double BETA  = 0.45;

// EBM CONSTANTS FOR GOWALLA
double ALPHA = 0.480;
double GAMMA = 0;
double BETA  = 0.520;

double RENY_Q=0.1; // Order of diversity
// int TIME_RANGE_IN_SECONDS = 1200;

// Global parameters
char *checkins_file, *graph_file, *query_file, *query_test_file;
int iteration_type = 1;
bool run_utilties = false;
bool is_noise_method = false;
int data_set = 0, k;
double noise_radius, group_radius, grid_size_in_km, locality_treshold, entropy_treshold,
         total_parts, part_number, distance_treshold, social_strength_tresh, combination_type,
         function_type, p1, p2, p3, p4, p5, p6,
         time_deviation, group_time_radius,
         coocc_time_range = TIME_RANGE_IN_SECONDS,
         coocc_spatial_range = SPATIAL_RANGE_IN_METERS,
         day_of_week, time_block, noise_type, noise_function;
#endif

void printParameters(){
  cout << "------------- Input Parameters -------------" << endl;
  cout << "noise_radius         : " << noise_radius << endl;
  cout << "group_radius         : " << group_radius << endl;
  cout << "grid_size_in_km      : " << grid_size_in_km << endl;
  cout << "locality_treshold    : " << locality_treshold << endl;
  cout << "entropy_treshold     : " << entropy_treshold << endl;
  cout << "total_parts          : " << total_parts << endl;
  cout << "part_number          : " << part_number << endl;
  cout << "distance_treshold    : " << distance_treshold << endl;
  cout << "social_strength_tresh: " << social_strength_tresh << endl;
  cout << "combination_type     : " << combination_type << endl;
  cout << "function_type        : " << function_type << endl;
  cout << "time_deviation       : " << time_deviation << endl;
  cout << "group_time_radius    : " << group_time_radius << endl;
  cout << "coocc_time_range     : " << coocc_time_range << endl;
  cout << "coocc_spatial_range  : " << coocc_spatial_range << endl;
  cout << "day_of_week          : " << day_of_week << endl;
  cout << "time_block           : " << time_block << endl;
  cout << "noise_type           : " << noise_type << endl;
  cout << "K                    : " << k << endl;
  cout << "------------------------------------------" << endl;
}

GPOs* loadCheckins(char* checkins_file){
  cout << "------------- Loading checkins ---------------" << endl;
  GPOs* gpos = new GPOs(checkins_file, (uint) coocc_time_range, coocc_spatial_range);
  return gpos;
}

GPOs* loadCheckins(char* checkins_file, bool preload_LE, bool preload_OCC){
  GPOs* gpos = loadCheckins(checkins_file);

  if(preload_LE)
    gpos->calculateLocationEntropy();

  if(preload_OCC){
    gpos->generateCooccurrenceCache();
    gpos->countU2UCoOccurrences();
  }

  return gpos;
}

SPOs* loadSocialGraph(char* graph_file, GPOs *gpos){
  cout << "------------- Loading SocialGraph ---------------" << endl;
  SPOs* spos = new SPOs(gpos);
  spos->load(graph_file);
  return spos;
}

SPOs* loadSocialGraph(char* graph_file){
  cout << "------------- Loading SocialGraph ---------------" << endl;
  SPOs* spos = new SPOs();
  spos->load(graph_file);
  return spos;
}

void runRangeUtility(GPOs *purturbedGPOs, GPOs *baseGPOs, SPOs *spos){
  SimpleQueries* query = new SimpleQueries(purturbedGPOs, spos);

  cout << "------------- Evaluating range utility ---------------" << endl;

  double base_radius = (double) max((int)noise_radius, (int)800);
  query->checkUtilityRange(query_file, baseGPOs, base_radius, noise_radius);
}

void runProximityUtility(GPOs *purturbedGPOs, GPOs *baseGPOs, SPOs *spos){
  spos->loadKatzScoreFromMemory();

  SimpleQueries* query = new SimpleQueries(purturbedGPOs, spos);

  cout << "------------- Evaluating range proximity ---------------" << endl;

  double base_radius = (double) max((int)noise_radius, (int)800);
  query->checkUtilityProximity(query_file, baseGPOs, base_radius, 100, noise_radius);
}

void runBasicUtility(GPOs *purturbedGPOs, GPOs *baseGPOs, SPOs *spos){
  SimpleQueries* query = new SimpleQueries(purturbedGPOs, spos);

  cout << "------------- Evaluating basic utility ---------------" << endl;

  double base_radius = (double) max((int)noise_radius, (int)800);
  query->checkUtilityBasic(baseGPOs);
}

void runUtilities(GPOs *purturbedGPOs, GPOs *baseGPOs, SPOs *spos){
  runRangeUtility(purturbedGPOs, baseGPOs, spos);
  // runProximityUtility(purturbedGPOs, baseGPOs, spos);
}

void runEBMWithoutGroundTruth(){
  bool preload_LE  = true;
  bool preload_OCC = true;

  GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = new SPOs(); // Dummy spos

  SimpleQueries* query = new SimpleQueries(gpos, spos);

  cout << "----- Precomputing matrices --- " << endl;
  query->buildMatrices(RENY_Q);

  cout << "----- Calculating Social Strength --- " << endl;
  query->cacluateSocialStrength();

  // Compute stats: Number of Friendships inferred vs Threshold
  for(double i = 0; i < 10; i = i + 0.1){
    query->countEBMInferredFriendships(i);
  }
}

void runEBM(GPOs *gpos, SPOs *spos){
  SimpleQueries* query = new SimpleQueries(gpos, spos);

  cout << "----- Precomputing matrices --- " << endl;
  query->buildMatrices(RENY_Q);

  cout << "----- Calculating Social Strength --- " << endl;
  query->cacluateSocialStrength();

  // for(double i = 0; i < 10; i = i + 0.1){
  //  query->verifySocialStrength(i);
  // }
  query->computeAccuracyOfSocialStrength(0.80);
}

// void runBasicOnNoised(GPOs *baseGPOs, GPOs *cmpGPOs, SPOs *spos, bool areFriends){
//   // cmpGPOs->countU2UCoOccurrences();
//   // cmpGPOs->calculateLocationEntropy();
//   // runEBM(cmpGPOs, spos);

//   //grab cooccurence Matrix from baseGPOs
//   //grab cooccurence matrix from cmpGPOs
//   //compare them both.
//   map<int, pair<double,double>> user_to_precision_recall;

//   auto cooccurence_matrix_base = baseGPOs->getCooccurrenceMatrix();
//   auto cooccurence_matrix_cmp = cmpGPOs->getCooccurrenceMatrix();

//   unordered_set<int> seen_users_in_base;


//   //printing--------
//   ofstream outfile;
//   ostringstream ss;
//   ss<<"RES_Basic" << p1<<"_"<<p2<<"_"<<p3<<"_"<<p4;
//   outfile.open(ss.str());

//   double precision=0,recall=0,precision_in_friends=0,recall_in_friends=0;
//   int sum_of_weights_precision=0;
//   int sum_of_weights_recall=0;

//   for(auto c_it = cooccurence_matrix_base->begin(); c_it != cooccurence_matrix_base->end(); c_it++){
//     int user_1 = c_it->first;
//     seen_users_in_base.insert(user_1);
//     auto users_location_frequency_map_base = c_it->second;
//     map<int, vector<pair<int, int> >* >* users_location_frequency_map_cmp = NULL;


//     auto iter_outer = cooccurence_matrix_cmp->find(user_1);
//     if(iter_outer != cooccurence_matrix_cmp->end()){
//       users_location_frequency_map_cmp = iter_outer->second;
//       int cooccurrences_across_users_base = 0;
//       int cooccurrences_across_users_cmp = 0;
//       int cooccurrences_across_users_min = 0;
//       for(auto ulh_it = users_location_frequency_map_base->begin(); ulh_it != users_location_frequency_map_base->end(); ulh_it++){
//         int user_2 = ulh_it->first;

//         if(areFriends && !spos->areFriends(user_1, user_2)){
//           continue;
//         }

//         vector<pair<int, int>>* cooccurrence_counts_vector_base = ulh_it->second;
//         vector<pair<int, int>>* cooccurrence_counts_vector_cmp = NULL;

//         int cooccurrence_count_cmp = 0;
//         int cooccurrence_count_base = 0;
//         auto iter_inner = users_location_frequency_map_cmp->find(user_2);
//         if(iter_inner != users_location_frequency_map_cmp->end()){
//           cooccurrence_counts_vector_cmp = iter_inner->second;
//           for(auto l_it=cooccurrence_counts_vector_cmp->begin(); l_it != cooccurrence_counts_vector_cmp->end(); l_it++){
//             int cooccrences_at_l = l_it->second;
//             cooccurrence_count_cmp += cooccrences_at_l;
//           }
//         }else{
//           // do nothing because cooccurences are zero
//           // cout<<"ERROR ? baseGPOs and cmpGPOs do not have similar cooccurrence_counts_vector"<<endl;
//         }

//         for(auto l_it=cooccurrence_counts_vector_base->begin(); l_it != cooccurrence_counts_vector_base->end(); l_it++){
//           int cooccrences_at_l = l_it->second;
//           cooccurrence_count_base += cooccrences_at_l;
//         }

//         cooccurrences_across_users_base += cooccurrence_count_base;
//         cooccurrences_across_users_cmp += cooccurrence_count_cmp;
//         cooccurrences_across_users_min += min(cooccurrence_count_cmp,cooccurrence_count_base);
//       }

//       if(cooccurrences_across_users_cmp != 0){
//         double temp_precision = cooccurrences_across_users_min/(double)cooccurrences_across_users_cmp;
//         double temp_recall = cooccurrences_across_users_min/(double)cooccurrences_across_users_base;
//         if(temp_recall > 1) {temp_recall=1;}
//         precision += temp_precision * users_location_frequency_map_base->size();
//         recall += temp_recall * users_location_frequency_map_base->size();
//         sum_of_weights_precision += users_location_frequency_map_base->size();
//         sum_of_weights_recall += users_location_frequency_map_base->size();
//       }
//       else{
//         sum_of_weights_precision += users_location_frequency_map_base->size();
//         sum_of_weights_recall += users_location_frequency_map_base->size();
//       }
//       // outfile<<it->first<< " "<<precision<<" "<<recall<<" "<<cooccurrence_counts_vector_base->size()<<endl;
//       // user_to_precision_recall.insert(make_pair(user_1,make_pair(precision,recall)));
//     }else{
//       //case when user is not there in cmp cooccurrence matrix
//       sum_of_weights_precision += users_location_frequency_map_base->size();
//       sum_of_weights_recall += users_location_frequency_map_base->size();
//       // user_to_precision_recall.insert(make_pair(user_1,make_pair(0,0)));
//       // cout<<"ERROR ? baseGPOs and cmpGPOs do not have similar cooccurence matrices"<<endl;
//     }
//   }

//   for(auto c_it = cooccurence_matrix_cmp->begin(); c_it != cooccurence_matrix_cmp->end(); c_it++){
//     int user = c_it->first;
//     auto users_location_frequency_map_base = c_it->second;
//     if(seen_users_in_base.find(user) == seen_users_in_base.end()){
//       sum_of_weights_precision += users_location_frequency_map_base->size();
//       // user_to_precision_recall.insert(make_pair(user,make_pair(0,0)));
//     }
//   }
//   double res_precision = precision/(double)sum_of_weights_precision;
//   double res_recall    = recall/(double)sum_of_weights_recall;
//   double f1 = 2 * res_precision * res_recall / (res_precision + res_recall);

//   if(areFriends){
//     cout << "are_freinds_basic_metric_precision{{" << res_precision << "}}" << endl;
//     cout << "are_freinds_basic_metric_recall{{" << res_recall << "}}" << endl;
//     cout << "are_freinds_basic_metric_f1{{" << f1 << "}}" << endl;
//   } else {
//     cout << "basic_metric_precision{{" << res_precision << "}}" << endl;
//     cout << "basic_metric_recall{{" << res_recall << "}}" << endl;
//     cout << "basic_metric_f1{{" << f1 << "}}" << endl;
//   }


//   // for(auto it = user_to_precision_recall.begin(); it != user_to_precision_recall.end(); it++){
//   //   auto precision_recall_pair = it->second;
//   //   double precision = precision_recall_pair.first;
//   //   double recall = precision_recall_pair.second;
//   //   outfile<<it->first<< " "<<precision<<" "<<recall<<endl;
//   // }
//   // outfile.close();
//   //------------------
// }


void runEBMOnNoised(GPOs *baseGPOs, GPOs *purturbedGPOs, GPOs *cmpGPOs, SPOs *spos){
  cmpGPOs->countU2UCoOccurrences();
  cmpGPOs->calculateLocationEntropy();
  runEBM(cmpGPOs, spos);

  if(run_utilties){
    runBasicUtility(cmpGPOs, baseGPOs, spos);
    runUtilities(purturbedGPOs, baseGPOs, spos);
  }
}

void plainEBM(){
  bool preload_LE  = true;
  bool preload_OCC = true;

  GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, gpos);

  runEBM(gpos, spos);
}

void selectiveGaussianNoiseIdealGrouping(int isOptimistic, int f){
  bool preload_LE  = false;
  bool preload_OCC = false;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  GPOs* fixedGPOs = baseGPOs;
  fixedGPOs->countCoOccurrencesOptimistic();

  double spatial_radi[] =  { 1.35 };
  double temporal_radi[] = { 1.15 };

  double noise_radius   = 100 * f;
  double time_deviation = 1200 * f;

  cout << "Using spatial noise : (m)"  << noise_radius << endl;
  cout << "Using time    noise : (mi)" << time_deviation/60 << endl;

  GPOs* purturbedGPOs = new GPOs(coocc_time_range, coocc_spatial_range);
  purturbedGPOs->loadPurturbedBasedOnSelectiveGaussian(fixedGPOs, noise_radius, time_deviation);
  if(run_utilties){
    runUtilities(purturbedGPOs, fixedGPOs, spos);
  }

  double mean_radius_spatial  = (double) purturbedGPOs->total_spatial_displacement / (double) purturbedGPOs->purturbed_count;
  double mean_radius_temporal = (double) purturbedGPOs->total_time_displacement    / (double) purturbedGPOs->purturbed_count;

  cout << "Mean Radius Spatial  :" << mean_radius_spatial  << endl;
  cout << "Mean Radius Temporal :" << mean_radius_temporal << endl;

  if(isOptimistic == 1)
    cout << "OPTIMISTIC GROUPING STRATEGY" << endl;
  else
    cout << "PESIMISTIC GROUPING STRATEGY" << endl;

  for(int i=0; i<3; i++){
    for(int j=0; j<2; j++){
      double sg = spatial_radi[i]  * mean_radius_spatial;
      double tg = temporal_radi[j] * mean_radius_temporal;

      cout << "Using Spatial  Grouping (m):  " << sg * 1000 << endl;
      cout << "Using Temporal Grouping (mi): " << tg * 60   << endl;

      GPOs* cmpGPOs;
      if(!isOptimistic){
        coocc_spatial_range   = 0;
        coocc_time_range      = 1;
        cmpGPOs  = new GPOs(coocc_time_range, coocc_spatial_range);
        cmpGPOs->groupLocationsByST(purturbedGPOs, sg, tg);
        cmpGPOs->countCoOccurrencesOptimistic();
      } else {
        cmpGPOs  = new GPOs(purturbedGPOs);
        cmpGPOs->coocc_spatial_range   = sg * 1000;
        cmpGPOs->coocc_time_range      = tg * 3600;
        cmpGPOs->countCoOccurrencesOptimistic();
      }

      if(run_utilties){
        runBasicUtility(cmpGPOs, fixedGPOs, spos);
      }

      delete cmpGPOs;
    }
  }

}

void selectiveGaussianNoise(int isOptimistic){
  bool preload_LE  = false;
  bool preload_OCC = false;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  GPOs* fixedGPOs = baseGPOs;
  fixedGPOs->countCoOccurrencesOptimistic();

  for(int i=1; i<=7; i++){
    double noise_radius   = 100 * i;
    double time_deviation = 1200 * i;

    cout << "Using spatial noise : (m)"  << noise_radius << endl;
    cout << "Using time    noise : (mi)" << time_deviation/60 << endl;

    GPOs* purturbedGPOs = new GPOs(coocc_time_range, coocc_spatial_range);
    purturbedGPOs->loadPurturbedBasedOnSelectiveGaussian(fixedGPOs, noise_radius, time_deviation);
    if(run_utilties){
      runUtilities(purturbedGPOs, fixedGPOs, spos);
    }
    double mean_radius_spatial  = (double) purturbedGPOs->total_spatial_displacement / (double) purturbedGPOs->purturbed_count;
    double mean_radius_temporal = (double) purturbedGPOs->total_time_displacement    / (double) purturbedGPOs->purturbed_count;

    cout << "Mean Radius Spatial  :" << mean_radius_spatial  << endl;
    cout << "Mean Radius Temporal :" << mean_radius_temporal << endl;

    if(isOptimistic == 1)
      cout << "OPTIMISTIC GROUPING STRATEGY" << endl;
    else
      cout << "PESIMISTIC GROUPING STRATEGY" << endl;

    double sg = mean_radius_spatial * group_radius;
    double tg = mean_radius_temporal * group_time_radius;

    cout << "Using Spatial  Grouping (m):  " << sg * 1000 << endl;
    cout << "Using Temporal Grouping (mi): " << tg * 60   << endl;

    GPOs* cmpGPOs;
    if(!isOptimistic){
      coocc_spatial_range   = 0;
      coocc_time_range      = 1;
      cmpGPOs  = new GPOs(coocc_time_range, coocc_spatial_range);
      cmpGPOs->groupLocationsByST(purturbedGPOs, sg, tg);
      cmpGPOs->countCoOccurrencesOptimistic();
    } else {
      cmpGPOs  = new GPOs(purturbedGPOs);
      cmpGPOs->coocc_spatial_range   = sg * 1000;
      cmpGPOs->coocc_time_range      = tg * 3600;
      cmpGPOs->countCoOccurrencesOptimistic();
    }

    if(run_utilties){
      runBasicUtility(cmpGPOs, fixedGPOs, spos);
    }

    delete purturbedGPOs;
    delete cmpGPOs;
  }
}

void selectiveGaussianNoiseDDAdversary(int k, int isOptimistic){
  bool preload_LE  = false;
  bool preload_OCC = false;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);
  baseGPOs->countCoOccurrencesOptimistic();

  for(int i=1; i<=7; i++){
    double noise_radius   = 100 * i;
    double time_deviation = 1200 * i;

    GPOs* purturbedGPOs = new GPOs(coocc_time_range, coocc_spatial_range);
    purturbedGPOs->loadPurturbedBasedOnSelectiveGaussian(baseGPOs, noise_radius, time_deviation);

    GPOs* cmpGPOs;
    if(!isOptimistic){
      cmpGPOs       = new GPOs(coocc_time_range,coocc_spatial_range);
      cmpGPOs->groupLocationsByDD(purturbedGPOs, k);
      cmpGPOs->countCoOccurrencesOptimistic();
    } else {
      cmpGPOs  = new GPOs(purturbedGPOs);
      cmpGPOs->countCoOccurrencesOptimisticDD(k);
    }

    if(run_utilties){
      runBasicUtility(cmpGPOs, baseGPOs, spos);
      runUtilities(purturbedGPOs, baseGPOs, spos);
    }

    delete purturbedGPOs;
    delete cmpGPOs;
  }
}


void selectiveSTKNNNoise(int k){
  bool preload_LE  = false;
  bool preload_OCC = false;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);
  baseGPOs->countCoOccurrencesOptimistic();

  GPOs* purturbedGPOs = new GPOs(coocc_time_range, coocc_spatial_range);
  purturbedGPOs->loadPurturbedBasedOnSelectiveSTKNNDistance(baseGPOs, k);
  if(run_utilties){
    runUtilities(purturbedGPOs, baseGPOs, spos);
  }
  purturbedGPOs->countCoOccurrencesOptimistic();

  double mean_radius_spatial  = (double) purturbedGPOs->total_spatial_displacement / (double) purturbedGPOs->purturbed_count;
  double mean_radius_temporal = (double) purturbedGPOs->total_time_displacement    / (double) purturbedGPOs->purturbed_count;

  cout << "Mean Radius Spatial  :" << mean_radius_spatial  << endl;
  cout << "Mean Radius Temporal :" << mean_radius_temporal << endl;

  double sg = mean_radius_spatial  * 1.35;
  double tg = mean_radius_temporal * 1.15;

  cout << "Using Spatial  Grouping (m):  " << sg * 1000 << endl;
  cout << "Using Temporal Grouping (mi): " << tg * 60   << endl;

  // {
  //   GPOs* cmpGPOs       = new GPOs(coocc_time_range, coocc_spatial_range);
  //   cmpGPOs->groupLocationsByST(purturbedGPOs, sg, tg);
  //   cmpGPOs->countCoOccurrencesOptimistic();
  //   if(run_utilties){
  //     runBasicUtility(cmpGPOs, baseGPOs, spos);
  //   }
  //   delete cmpGPOs;
  // }

  {
    GPOs* cmpGPOs       = new GPOs(coocc_time_range, coocc_spatial_range);
    cmpGPOs->groupLocationsByDD(purturbedGPOs, std::max(k, 1));
    cmpGPOs->countCoOccurrencesOptimistic();
    if(run_utilties){
      runBasicUtility(cmpGPOs, baseGPOs, spos);
    }
    delete cmpGPOs;
  }
}

void selectiveSkylineNoise(int k){
  bool preload_LE  = false;
  bool preload_OCC = true;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  double spatial_grouping[]  = { 0.05, 0.10, 0.25, 0.5, 0.75 };
  double temporal_grouping[] = { 0.05, 0.10, 0.25, 0.5, 0.75 };

  GPOs* purturbedGPOs = new GPOs(coocc_time_range,coocc_spatial_range);
  purturbedGPOs->loadPurturbedBasedOnSelectiveSkyline(baseGPOs, k);
  if(run_utilties){
    runUtilities(purturbedGPOs, baseGPOs, spos);
  }

  double group_radius_spatial  = (double) purturbedGPOs->total_spatial_displacement / (double) purturbedGPOs->purturbed_count;
  double group_radius_temporal = (double) purturbedGPOs->total_time_displacement    / (double) purturbedGPOs->purturbed_count;

  cout << "Mean Radius Spatial  :" << group_radius_spatial  << endl;
  cout << "Mean Radius Temporal :" << group_radius_temporal << endl;

  for(int i=0; i<5;i++){
    for(int j=0; j<5;j++){
      double sg = group_radius_spatial * 1000.0 * spatial_grouping[ i ] + 3.3;
      double tg = group_radius_temporal * 3600.0 * temporal_grouping[ j ] + 180;

      cout << "Using Spatial  Grouping : " << sg << endl;
      cout << "Using Temporal Grouping : " << tg << endl;

      GPOs* cmpGPOs       = new GPOs(coocc_time_range,coocc_spatial_range);
      cmpGPOs->groupLocationsByST(purturbedGPOs, sg, tg);
      cmpGPOs->countU2UCoOccurrences();

      if(run_utilties){
        runBasicUtility(cmpGPOs, baseGPOs, spos);
      }

      delete cmpGPOs;
    }
  }

}

void gaussianNoiseVsEBM(double noise_radius, double group_radius, double time_deviation, bool add_spatial, bool add_temporal){
  bool preload_LE  = false;
  bool preload_OCC = true;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  GPOs* spatiallyPurturbedGPOs = new GPOs(coocc_time_range,coocc_spatial_range);
  GPOs* spatiallyAndTemporallyPurturbedGPOs = new GPOs(coocc_time_range,coocc_spatial_range);
  GPOs* cmpGPOs  = new GPOs(coocc_time_range,coocc_spatial_range);

  if( add_spatial ){
    spatiallyPurturbedGPOs->loadPurturbedLocations(baseGPOs, noise_radius);
    cout << "------------- Locations spatially perturbed -------------------" << endl;
  }

  if( add_temporal ){
    spatiallyAndTemporallyPurturbedGPOs->loadPurturbedLocationsByTime(spatiallyPurturbedGPOs, time_deviation);
    cout << "------------- Locations temporally perturbed -------------------" << endl;
  }

  cmpGPOs->groupLocationsByRange(spatiallyAndTemporallyPurturbedGPOs, group_radius, false);
  cout << "------------- Locations Grouped -------------------" << endl;

  runEBMOnNoised(baseGPOs, spatiallyAndTemporallyPurturbedGPOs, cmpGPOs, spos);
}

void CombinationNoiseVsEBM(double noise_radius, double time_deviation, bool add_spatial, bool add_temporal, int noise_function){
  bool preload_LE  = false;
  bool preload_OCC = true;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  cout << "------------- Load computed checkin locality ---------------" << endl;
  SPOs *tmp_spos = new SPOs();

  GPOs* purturbedGPOs = new GPOs(coocc_time_range,coocc_spatial_range);
  GPOs* cmpGPOs       = new GPOs(coocc_time_range,coocc_spatial_range);

  purturbedGPOs->loadPurturbedLocationsBasedOnCombinationFunction(
    baseGPOs,
    tmp_spos->loadCheckinLocalityFromFile(DATASET_PATH),
    tmp_spos->loadTemporalLocalityFromFile(DATASET_PATH),
    baseGPOs->getL2U2COOCC(),
    noise_radius,
    (uint) time_deviation,
    add_spatial,
    add_temporal,
    noise_function);
  cout << "------------- Locations perturbed -------------------" << endl;

  cmpGPOs->groupLocationsByRange(purturbedGPOs, group_radius, false);
  cout << "------------- Locations Grouped -------------------" << endl;

  runEBMOnNoised(baseGPOs, purturbedGPOs, cmpGPOs, spos);
}


void checkQueryFileStats(){
  cout << "------------- Evaluating query locations ---------------" << endl;

  GPOs* baseGPOs = loadCheckins(checkins_file);
  SPOs* spos = loadSocialGraph(graph_file);

  SimpleQueries* query = new SimpleQueries(baseGPOs, spos);

  query->getInterestingQueryPoints(query_test_file, 800, 100, query_file, DATA_SET);

  query->checkUtilityStats(query_file, 800, 100);
  query->checkUtilityStats(query_file, 800, 200);
  query->checkUtilityStats(query_file, 800, 300);
  query->checkUtilityStats(query_file, 800, 400);
  query->checkUtilityStats(query_file, 800, 500);
  query->checkUtilityStats(query_file, 800, 600);
  query->checkUtilityStats(query_file, 800, 700);
}

void testpTools(){
  // int i;
  // double length, entropyTarget;
  double firstEntropy, secondEntropy, thirdEntropy, targetEntropy;
  // double firstMItarget, secondMItarget, thirdMItarget, targetMItarget;
  // int testFirstVector, testSecondVector, testThirdVector, testMergedVector;
  // struct timeval start,end;

  uint *firstVector = (uint *) calloc(4,sizeof(uint));
  uint *secondVector = (uint *) calloc(4,sizeof(uint));
  uint *thirdVector = (uint *) calloc(4,sizeof(uint));
  uint *targetVector = (uint *) calloc(4,sizeof(uint));

  firstVector[0] = 0;
  firstVector[1] = 0;
  firstVector[2] = 1;
  firstVector[3] = 1;

  secondVector[0] = 0;
  secondVector[1] = 1;
  secondVector[2] = 0;
  secondVector[3] = 1;

  thirdVector[0] = 0;
  thirdVector[1] = 1;
  thirdVector[2] = 1;
  thirdVector[3] = 1;

  targetVector[0] = 0;
  targetVector[1] = 1;
  targetVector[2] = 1;
  targetVector[3] = 0;

  firstEntropy = calcEntropy(firstVector,4);
  secondEntropy = calcEntropy(secondVector,4);
  thirdEntropy = calcEntropy(thirdVector,4);
  targetEntropy = calcEntropy(targetVector,4);
  // targetEntropy = calcRenyiEntropy(0.1,targetVector,4 );

  printf("Entropies - first: %f, second: %f, third: %f, target %f\n",firstEntropy,secondEntropy,thirdEntropy,targetEntropy);

  //verifying section 3.4 EBM

  cout<<" -----------------------Entropy Example ------------------------"<<endl;
  cout<<" ------------------------Section  3.4---------------------------"<<endl;
  uint *CabOcc = (uint *) calloc(5,sizeof(uint));
  uint *CcdOccur = (uint *) calloc(5,sizeof(uint));

  CabOcc[0] = 10;
  CabOcc[1] = 1;
  CabOcc[2] = 0;
  CabOcc[3] = 0;
  CabOcc[4] = 9;

  CcdOccur[0] = 2;
  CcdOccur[1] = 3;
  CcdOccur[2] = 2;
  CcdOccur[3] = 2;
  CcdOccur[4] = 3;

  cout<<" C_ab is : 10 , 1, 0, 0, 9"<<endl;
  cout<<" C_cd is : 2 , 3, 2, 2, 3"<<endl;

  cout<<" ------------------Testing Shannon Entropy------------------------"<<endl;
  double CabEntropy, CcdEntropy, CabDiversity, CcdDiversity;

  CabEntropy = calcEntropyFromCoV(CabOcc,5);
  CcdEntropy = calcEntropyFromCoV(CcdOccur,5);
  CabDiversity = exp(CabEntropy);
  CcdDiversity = exp(CcdEntropy);


  cout<<"CabEntropy = "<<CabEntropy<<" CcdEntropy = "<<CcdEntropy <<" CabEntropy/CcdEntropy = " << CabEntropy/CcdEntropy <<endl;
  cout<<"CabDiversity = "<<CabDiversity<<" CabDiversity = "<<CcdDiversity<<" CabDiversity/CcdDiversity = " << CabDiversity/CcdDiversity <<endl;

  cout<<" ------------------Testing Renyi Entropy------------------------"<<endl;
  double CabRenEntropy, CcdRenEntropy, CabRenDiversity, CcdRenDiversity;

  CabRenEntropy = calcRenyiEntropyFromCoV(0.1,CabOcc,5);
  CcdRenEntropy = calcRenyiEntropyFromCoV(0.1,CcdOccur,5);
  CabRenDiversity = exp(CabRenEntropy);
  CcdRenDiversity = exp(CcdRenEntropy);

  cout<<"CabRenEntropy = "<<CabRenEntropy<<" CcdRenEntropy = "<<CcdRenEntropy<< " CabRenEntropy/CcdRenEntropy = " << CabRenEntropy/CcdRenEntropy<<endl;
  cout<<"CabRenDiversity = "<<CabRenDiversity<<" CcdRenDiversity = "<<CcdRenDiversity<<" CabRenDiversity/CcdRenDiversity = "<< CabRenDiversity/CcdRenDiversity << endl;
}

int main(int argc, char *argv[]){
  cout.precision(15);

  // testpTools();

  if (argc != 10){
    cout << "Usage: data_set [ebm|grouped-ebm|gaussian-v-ebm|snapping-v-ebm|nl-v-ebm|le-v-ebm|occ-hist|compute-katz|ebm-without-gt|compute-social-graph] run_utilties p1 p2 p3 p4 p5 p6 " << endl;
    return -1;
  }

  data_set = atoi(argv[1]);

  if(data_set == 0){
    cout << "Using dataset GOWALLA_SMALL" << endl;

    checkins_file    = "data_Gowalla/checkins.txt";
    graph_file       = "data_Gowalla/socialGraph.txt";
    query_file       = "data_Gowalla/queries.txt";
    query_test_file  = "data_Gowalla/queries-initial.txt";

    MIN_X = -124.848974;
    MAX_X = -66.885444;
    MIN_Y = 24.396308;
    MAX_Y = 49.384358;

    DATA_SET = 0;
    DATASET_SIZE = 3669249;
    DATASET_PATH = "data_Gowalla/";

    cout << "Using path: " << DATASET_PATH << endl;

  } else if(data_set == 1){
    cout << "Using dataset GOWALLA_LARGE" << endl;

    checkins_file    = "data_GowallaFull/checkins.txt";
    graph_file       = "data_GowallaFull/socialGraph.txt";
    query_file       = "data_GowallaFull/queries.txt";
    query_test_file  = "data_GowallaFull/queries-initial.txt";


    MIN_X = -124.848974;
    MAX_X = -66.885444;
    MIN_Y = 24.396308;
    MAX_Y = 49.384358;

    DATA_SET = 1;
    DATASET_SIZE = 18290199;
    DATASET_PATH     = "data_GowallaFull/";

    cout << "Using path: " << DATASET_PATH << endl;

  } else if(data_set == 2){
    cout << "Using dataset SHANGHAI"      << endl;

    checkins_file    = "data_Shanghai/checkins.txt";
    graph_file       = "data_Shanghai/socialGraph.txt";
    query_file       = "data_Shanghai/queries.txt";
    query_test_file  = "data_Shanghai/queries-initial.txt";


    MIN_X = 120.711564;
    MAX_X = 122.222184;
    MIN_Y = 30.343011;
    MAX_Y = 32.417996;

    DATA_SET = 2;
    DATASET_SIZE = 7858442;
    DATASET_PATH     = "data_Shanghai/";

    cout << "Using path: " << DATASET_PATH << endl;

  } else if(data_set == 3){
    cout << "Using dataset SIMULATED"      << endl;

    checkins_file    = "data_Simulated/checkins.txt";
    graph_file       = "data_Simulated/socialGraph.txt";
    query_file       = "data_Simulated/queries.txt";
    query_test_file  = "data_Simulated/queries-initial.txt";


    MIN_X = -98.15;
    MAX_X = -97;
    MIN_Y = 30.02;
    MAX_Y = 30.70;

    DATA_SET = 2;
    DATASET_SIZE = 7858442;
    DATASET_PATH     = "data_Simulated/";

    cout << "Using path: " << DATASET_PATH << endl;

  } else {
    cout << "Invalid data_set"            << endl;
    exit(-1);
  }

  DELTA_X = ((MAX_X - MIN_X)/ (X-1));
  DELTA_Y = ((MAX_Y - MIN_Y)/ (Y-1));

  if (strcmp(argv[2], "ebm") == 0)
    iteration_type = 1;
  else if (strcmp(argv[2], "grouped-ebm") == 0)
    iteration_type = 2;
  else if (strcmp(argv[2], "gaussian-v-ebm") == 0)
    iteration_type = 3;
  else if (strcmp(argv[2], "comb-v-ebm") == 0)
    iteration_type = 4;

  else if (strcmp(argv[2], "ebm-without-gt") == 0)
    iteration_type = 5;

  else if (strcmp(argv[2], "selective-gaussian") == 0)
    iteration_type = 6;

  else if (strcmp(argv[2], "selective-skyline") == 0)
    iteration_type = 7;

  else if (strcmp(argv[2], "selective-stknn") == 0)
    iteration_type = 8;

  else if (strcmp(argv[2], "selective-gaussian-dd") == 0)
    iteration_type = 9;

  else if (strcmp(argv[2], "selective-gaussian-ideal-grouping") == 0)
    iteration_type = 10;

  else if (strcmp(argv[2], "occ-hist") == 0)
    iteration_type = 90;
  else if (strcmp(argv[2], "compute-katz") == 0)
    iteration_type = 91;
  else if (strcmp(argv[2], "compute-locality") == 0)
    iteration_type = 92;
  else if (strcmp(argv[2], "compute-histograms") == 0)
    iteration_type = 93;
  else if (strcmp(argv[2], "query-metrics") == 0)
    iteration_type = 94;

  else if (strcmp(argv[2], "compute-social-graph") == 0)
    iteration_type = 95;

  else if (strcmp(argv[2], "compute-knn-metrics") == 0)
    iteration_type = 96;

  else if (strcmp(argv[2], "compute-skyline-metrics") == 0)
    iteration_type = 97;

  else if (strcmp(argv[2], "compute-stknn-metrics") == 0)
    iteration_type = 98;

  else if (strcmp(argv[2], "compute-all-stknn-metrics") == 0)
    iteration_type = 99;

  else if (strcmp(argv[2], "compute-coocc-metrics") == 0)
    iteration_type = 100;

  else
    iteration_type = -1;

  run_utilties       = atoi(argv[3]); // [ 0 / 1]

  p1 = atof(argv[4]); // ( Parameter 1 )
  p2 = atof(argv[5]); // ( Parameter 2 )
  p3 = atof(argv[6]); // ( Parameter 3 )
  p4 = atof(argv[7]); // ( Parameter 4 )
  p5 = atof(argv[8]); // ( Parameter 5 )
  p6 = atof(argv[9]); // ( Parameter 6 )


  switch(iteration_type){
    case 1:{
      cout << "ITRATION: Running EBM" << endl;
      noise_radius            = 0;
      time_deviation          = 0;
      group_radius            = 0;
      coocc_time_range   = p4;
      printParameters();
      plainEBM();
      break;
    }

    case 2:{
      cout << "ITRATION: Running EBM with grouping" << endl;
      noise_radius            = 0;
      time_deviation          = 0;
      group_radius            = p3;
      coocc_time_range        = p4;
      printParameters();
      gaussianNoiseVsEBM(noise_radius, group_radius, time_deviation, true, true);
      break;
    }

    case 3:{
      cout << "ITRATION: Running EBM with gaussian noise" << endl;
      noise_radius            = p1;
      time_deviation          = p2;
      group_radius            = p3;
      coocc_time_range   = p4;
      noise_type              = p5;
      bool add_temporal=false, add_spatial=false;
      printParameters();
      if( noise_type == 0 ){
        add_spatial = true;
        add_temporal = true;
        cout << "Adding both spatial and temporal noise" << endl;
      } else if( noise_type == 1){
        add_spatial = true;
        add_temporal = false;
        cout << "Adding spatial noise" << endl;
      } else if( noise_type == 2){
        add_spatial = false;
        add_temporal = true;
        cout << "Adding temporal noise" << endl;
      } else {
        cout << "Invalid option" << endl;
      }
      gaussianNoiseVsEBM(noise_radius, group_radius, time_deviation, add_spatial, add_temporal);
      break;
    }

    case 4:{
      cout << "ITRATION: Running EBM with combination noise" << endl;
      noise_radius            = p1;
      time_deviation          = p2;
      group_radius            = p3;
      coocc_time_range        = p4;
      noise_type              = p5;
      noise_function          = p6;
      bool add_temporal=false, add_spatial=false;

      if(noise_function == 0)
        cout << "Using Noise function HIj and HiL" << endl;
      else
        cout << "Using Noise function Cij" << endl;

      printParameters();
      if( noise_type == 0 ){
        add_spatial = true;
        add_temporal = true;
        cout << "Adding both spatial and temporal noise" << endl;
      } else if( noise_type == 1){
        add_spatial = true;
        add_temporal = false;
        cout << "Adding spatial noise" << endl;
      } else if( noise_type == 2){
        add_spatial = false;
        add_temporal = true;
        cout << "Adding temporal noise" << endl;
      } else {
        cout << "Invalid option" << endl;
      }
      CombinationNoiseVsEBM(noise_radius, time_deviation, add_spatial, add_temporal, noise_function);
      break;
    }

    case 5:{
      cout << "ITRATION: Running EBM without ground truth" << endl;
      noise_radius            = 0;
      time_deviation          = 0;
      group_radius            = 0;
      coocc_time_range        = p4;
      printParameters();
      runEBMWithoutGroundTruth();
      break;
    }

    case 6:{
      cout << "ITRATION: Selective Gaussian Noise" << endl;

      coocc_spatial_range     = p1;
      coocc_time_range        = p2;
      group_radius            = p3;
      group_time_radius       = p4;

      int isOptimistic        = p5;

      printParameters();
      selectiveGaussianNoise(isOptimistic);

      break;
    }

    case 7:{
      cout << "ITRATION: Selective Noise based on Skyline" << endl;
      k                   = p1;

      if(k == -1){
        cout << "Picking skylines at random " << endl;
      } else {
        cout << "TopK" << endl;
      }

      group_radius            = p3;
      time_block              = p4;

      printParameters();
      selectiveSkylineNoise(k);

      break;
    }

    case 8:{
      cout << "ITRATION: Selective Noise based on st-knn" << endl;

      k                   = p1;

      coocc_spatial_range = p2;
      coocc_time_range    = p3;

      printParameters();
      selectiveSTKNNNoise(k);

      break;
    }

    case 9:{
      cout << "ITRATION: Selective Noise vs DD Adversary" << endl;


      coocc_spatial_range     = p1;
      coocc_time_range        = p2;
      k                       = p3;
      int isOptimistic        = p4;

      printParameters();
      selectiveGaussianNoiseDDAdversary(k, isOptimistic);

      break;
    }

    case 10:{
      cout << "ITRATION: Selective Gaussian Noise Ideal grouping" << endl;

      coocc_spatial_range = p1;
      coocc_time_range    = p2;
      int f               = p3;
      int isOptimistic    = p4;

      printParameters();
      selectiveGaussianNoiseIdealGrouping(isOptimistic, f);

      break;
    }

    // case 4:
    //   cout << "ITRATION: Running EBM with grid snapping noise" << endl;
    //   grid_size_in_km         = p1;
    //   noise_radius            = grid_size_in_km * 1000;
    //   time_deviation          = p2;
    //   group_radius            = p3;
    //   coocc_time_range   = p4;

    //   printParameters();
    //   gridSnappingVsEBM(grid_size_in_km);
    //   break;

    case 90:{
      cout << "METRICS: Calculating Co occurrence distribution" << endl;
      bool preload_LE  = true;
      bool preload_OCC = true;

      printParameters();
      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = loadSocialGraph(graph_file);
      spos->loadNodeLocalityFromFile();

      SimpleQueries* query = new SimpleQueries(gpos, spos);

      query->cacluateCooccurrenceDistributionBasedOnNodeLocality();
      query->cacluateCooccurrenceDistributionBasedOnLocationEntropy();
      break;
    }

    case 91:{
      cout << "METRICS: Pre-Compute KATZ score" << endl;
      GPOs* gpos = loadCheckins(checkins_file);
      SPOs* spos = loadSocialGraph(graph_file);
      total_parts = p1;
      part_number = p2;
      distance_treshold = p3;

      printParameters();
      spos->precomputeKatzScore(gpos, total_parts, part_number, distance_treshold);
      break;
    }

    case 92:{
      cout << "METRICS: node_locality and temporal_locality" << endl;
      bool preload_LE  = false;
      bool preload_OCC = true;

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = loadSocialGraph(graph_file, gpos);

      int max_checkins          = 250;  // checkins in vicinity
      double max_radius         = 1000; // max_radius 2km

      printParameters();
      cout << "Data driven search Max checkins : " << max_checkins << endl;
      cout << "Data driven search Max radius   : " << max_radius << endl;

      spos->computeCheckinLocalityMap(gpos, std::numeric_limits<double>::infinity());
      spos->writeCheckinLocalityToFile(DATASET_PATH);

      spos->computeTemporalLocality(max_checkins, max_radius, gpos);
      spos->writeTemporalLocalityToFile(DATASET_PATH);

      // cout << "------------- Computing mean distance between friend pairs ---------------" << endl;
      // cout << "Mean distance between all pairs of friends :" << spos->computeMeanDistanceBetweenAllFriends(gpos) << endl;
      // spos->computeNodeLocality(gpos);
      break;
    }

    case 93:{
      cout << "METRICS: Compute histograms" << endl;

      social_strength_tresh = p1;
      time_block            = p2;

      printParameters();

      bool preload_LE  = true;
      bool preload_OCC = true;

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = loadSocialGraph(graph_file, gpos);

      SimpleQueries* query = new SimpleQueries(gpos, spos);
      query->buildMatrices(RENY_Q);
      query->cacluateSocialStrength();
      cout << "Using Threshold" << social_strength_tresh << endl;

      gpos->printCooccurrenceMatrix(DATASET_PATH);

      map< int, double >* temoral_locality_map = spos->loadTemporalLocalityFromFile(DATASET_PATH);
      query->writeHistogramstoFile(DATASET_PATH, social_strength_tresh, time_block, temoral_locality_map);

      spos->writeUserFriendsToFile(DATASET_PATH);
      break;
    }

    case 94:{
      printParameters();
      checkQueryFileStats();
      break;
    }

    case 95:{
      cout << "Additional: Generate Social Graph" << endl;

      social_strength_tresh = p1;

      bool preload_LE  = true;
      bool preload_OCC = true;

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = new SPOs(); // Dummy spos

      SimpleQueries* query = new SimpleQueries(gpos, spos);

      cout << "----- Precomputing matrices --- " << endl;
      query->buildMatrices(RENY_Q);

      cout << "----- Calculating Social Strength --- " << endl;
      query->cacluateSocialStrength();

      cout << "--- Computing social graph -- " << endl;
      query->generateSocialGraph(DATASET_PATH, social_strength_tresh);

      break;
    }

    case 96:{
      cout << "METRICS: Compute KNN Metrics" << endl;
      cout << "KNN for co_occurred locations " << endl;

      bool compute_spatial  = (p1 == 1);
      bool compute_temporal = (p2 == 1);

      if(compute_spatial)
        cout << "Computing spatial KNN " << endl;

      if(compute_temporal)
        cout << "Computing temporal KNN " << endl;

      printParameters();

      bool preload_LE  = false;
      bool preload_OCC = false;

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = loadSocialGraph(graph_file, gpos);
      gpos->countCoOccurrencesOptimistic();

      if(compute_spatial)
        gpos->computeSTKNNDistances(10, 1);

      if(compute_temporal)
        gpos->computeSTKNNDistances(10, 2);

      break;
    }

    case 97:{
      cout << "METRICS: Compute Skyline Metrics" << endl;
      cout << "Skyline for co_occurred locations " << endl;

      printParameters();

      bool preload_LE  = false;
      bool preload_OCC = true;

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = loadSocialGraph(graph_file, gpos);

      gpos->computeSkylineMetrics(gpos->getL2U2COOCC());

      break;
    }

    case 98:{
      cout << "METRICS: Compute ST KNN Metrics" << endl;
      cout << "ST KNN for co_occurred locations " << endl;

      coocc_spatial_range = p1;
      coocc_time_range    = p2;

      printParameters();

      bool preload_LE  = false;
      bool preload_OCC = false;

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = loadSocialGraph(graph_file, gpos);

      gpos->countCoOccurrencesOptimistic();
      gpos->computeSTKNNDistances(10, 0);
      break;
    }

    case 99:{
      cout << "METRICS: Compute ST KNN Metrics" << endl;
      cout << "ST KNN for all locations " << endl;

      coocc_spatial_range = p1;
      coocc_time_range    = p2;

      printParameters();

      bool preload_LE  = false;
      bool preload_OCC = false;

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = loadSocialGraph(graph_file, gpos);

      gpos->countCoOccurrencesOptimistic();
      gpos->computeSTKNNDistances(10, 3);
      break;
    }

    case 100:{
      cout << "METRICS: Co-Occurrence Metrics" << endl;

      bool preload_LE  = false;
      bool preload_OCC = false;

      double coocc_spatial_radius[] = { 0, 25, 50, 100, 200 };
      double coocc_temporal_radius[] = { 1, 20, 40, 60, 120 };

      printParameters();

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      GPOs* fixed = new GPOs(coocc_time_range, coocc_spatial_range);
      fixed->groupLocationsByRange(gpos, 10, false);
      delete gpos;

      for (int i = 0; i < 5; ++i){
        for (int j = 0; j < 5; ++j){
          GPOs* test_gpos = new GPOs(fixed);
          test_gpos->coocc_spatial_range = coocc_spatial_radius[ i ];
          test_gpos->coocc_time_range    = coocc_temporal_radius[ j ] * 60;
          test_gpos->generateCooccurrenceCache();
          test_gpos->countU2UCoOccurrences();
          delete test_gpos;
        }
      }

      break;
    }

    default:
      cout << "Invalid iteration type; Enter a valid option" << endl;
  }

  // Plain EBM


  // cout << "------------- Loading checkins ---------------" << endl;
  // GPOs* gpos = new GPOs(argv[2]);
  // cout << "------------- Loading complete ---------------" << endl;

  // cout << "------------- Loading checkins ---------------" << endl;
  // GPOs* g = new GPOs(argv[2]);
  // g->countU2UCoOccurrences();
  // GPOs* gpos = new GPOs(argv[2]);
  // gpos->groupLocationsByRange(g, r2, isOptimistic);
  // cout << "------------- Loading complete ---------------" << endl;

  // cout << "------------- Loading checkins ---------------" << endl;
  // GPOs* g = new GPOs(argv[2]);
  // g->countU2UCoOccurrences();
  // cout << "------------- Loading complete ---------------" << endl;


  // testing grid snapping
  // GPOs* gpos = new GPOs(coocc_time_range,coocc_spatial_range);
  // gpos->createNewGPOsbyGridSnapping(g, r1); //second var is the x distance of a cell is km
  // cout << "Number of locations loaded " << gpos->locations.size() << endl;
  // cout << "------------- Noise added -------------------" << endl;


  // test Gaussian noise with grouping
  // GPOs* gg = new GPOs(coocc_time_range,coocc_spatial_range);
  // GPOs* gpos = new GPOs(coocc_time_range,coocc_spatial_range);
  // gg->loadPurturbedLocations(g, r1);
  // cout << "------------- Locations perturbed -------------------" << endl;
  // gpos->groupLocationsByRange(gg, r2, isOptimistic);
  // cout << "------------- Locations Grouped -------------------" << endl;

  // test Gaussian for high node locality  noise without grouping
  // cout << "------------- Load computed node locality ---------------" << endl;
  // SPOs* tmp_spos = new SPOs();
  // map< int, double >* node_locality = tmp_spos->loadNodeLocalityFromFile();
  // GPOs* gpos = new GPOs(coocc_time_range,coocc_spatial_range);
  // gpos->loadPurturbedLocationsBasedOnNodeLocality(g, node_locality, r1, 0.75);
  // gpos->generateFrequencyCache();
  // cout << "------------- Locations perturbed Based on  Node locality -------------------" << endl;

  // test Gaussian for high location entropy noise without grouping
  // GPOs* gpos = new GPOs(coocc_time_range,coocc_spatial_range);
  // g->calculateLocationEntropy();
  // gpos->loadPurturbedLocationsBasedOnLocationEntropy(g, r1, 1);
  // gpos->generateFrequencyCache();
  // cout << "------------- Locations perturbed Based on  Location Entropy -------------------" << endl;


  // test gaussian noise without grouping
  // GPOs* gpos = new GPOs(coocc_time_range,coocc_spatial_range);
  // gpos->loadPurturbedLocations(g, r1);
  // cout << "------------- Locations perturbed -------------------" << endl;

  // cout << "----- Loading Cooccurrence Matrix --- " << endl;
  // gpos->countU2UCoOccurrences();
  // cout << "----- Completed Loading Cooccurrence Matrix --- " << endl;


  // cout << "----- Calculating Location Entropy --- " << endl;
  // gpos->calculateLocationEntropy();
  // cout << "----- Completed Calculating Location Entropy --- " << endl;

  // // Loading social network and checkins
  // cout << "------------- Loading SocialGraph ---------------" << endl;
  // SPOs* spos = new SPOs();
  // spos->load(argv[2]);

  // cout << "------------- Loading SocialGraph ---------------" << endl;
  // SPOs* spos = new SPOs(g);
  // spos->load(argv[2]);


  // cout << "------------- Loading SocialGraph ---------------" << endl;
  // SPOs* spos = new SPOs(gpos);
  // spos->load(argv[2]);

  // spos->loadKatzScoreFromMemory();

  // spos->precomputeKatzScore(gpos, r1, r2, coocc_time_range);

  // cout << "------------- Computing mean distance between friend pairs ---------------" << endl;
  // cout << "Mean distance between all pairs of friends :" << spos->computeMeanDistanceBetweenAllFriends(gpos) << endl;

  // cout << "------------- Computing node locality of all users ---------------" << endl;
  // spos->computeNodeLocality(gpos);
  // spos->writeNodeLocalityToFile();
  // spos->loadNodeLocalityFromFile();

  // SimpleQueries* query = new SimpleQueries(gpos, spos);
  // query->cacluateCooccurrenceDistributionBasedOnNodeLocality();
  // query->cacluateCooccurrenceDistributionBasedOnLocationEntropy();

  // cout << "------------- Evaluating range utility ---------------" << endl;
  // query->checkUtilityRange(locations_of_interest_file, g, 10);
  // query->checkUtilityRange(locations_of_interest_file, g, 25);
  // query->checkUtilityRange(locations_of_interest_file, g, 50);
  // query->checkUtilityRange(locations_of_interest_file, g, 100);
  // query->checkUtilityRange(locations_of_interest_file, g, 200);
  // query->checkUtilityRange(locations_of_interest_file, g, 400);


  // cout << "------------- Evaluating range proximity ---------------" << endl;
  // query->checkUtilityProximity(locations_of_interest_file, g, 1000, 100);
  // query->checkUtilityProximity(locations_of_interest_file, g, 500, 100);
  // query->checkUtilityProximity(locations_of_interest_file, g, 500, 200);


  // cout << "------------- Evaluating query locations ---------------" << endl;
  // query->checkUtilityStats(locations_of_interest_file, 10);
  // query->checkUtilityStats(locations_of_interest_file, 25);
  // query->checkUtilityStats(locations_of_interest_file, 50);
  // query->checkUtilityStats(locations_of_interest_file, 100);
  // query->checkUtilityStats(locations_of_interest_file, 200);
  // query->checkUtilityStats(locations_of_interest_file, 400);
  // query->checkUtilityStats(locations_of_interest_file, 800);
  // query->checkUtilityStats(locations_of_interest_file, 1600);


  // {
  //   cout << "----- Precomputing matrices --- " << endl;
  //   query->buildMatrices(0.1);
  //   cout << "----- Completed Precomputing matrices --- " << endl;

  //   cout << "----- Calculating Social Strength --- " << endl;
  //   query->cacluateSocialStrength();
  //   cout << "----- Completed calculating Social Strength --- " << endl;

  //   // double thresholds[9] = {0.5,1,2,3,4,5,7,10,25};

  //   for(double i = 1; i < 3; i = i + 0.1){
  //     cout << "----- Computing accuracy for threshold --- " << i <<endl;
  //     query->verifySocialStrength(i);
  //     cout << "--------------------------------------------";
  //   }

  //   // cout << "----- Computing accuracy for threshold --- 1.6 " <<endl;
  //   // query->verifySocialStrength(1.6);
  //   // cout << "--------------------------------------------";
  // }


  // map<int, Point*> locations;
  // ofstream output_file;
  // output_file.open("locations-of-interest.csv");

  // for (auto u = users_of_interest->begin(); u != users_of_interest->end(); u++){
  //   auto pts = gpos->getLocations(u->first);
  //   for(auto p=pts->begin(); p!=pts->end(); p++){
  //     double x=(*p)->getX(), y=(*p)->getY();
  //     locations.insert( make_pair((*p)->getID(), *p) );
  //     output_file << x << "," << y <<endl;
  //   }
  // }
  // cout << "------- Locations of interest : " << locations.size() << endl;
  // output_file.close();

  // if(r1 == 0){
  //   count_cooccurences(spos, gpos, r2, tR, isOptimistic);
  // } else {
  //   GPOs* purturbedGPOs = new GPOs(coocc_time_range,coocc_spatial_range);
  //   purturbedGPOs->loadPurturbedLocations(gpos, r1);
  //   cout << "------------- Noise added -------------------" << endl;

  //   count_cooccurences(spos, purturbedGPOs, r2, tR, isOptimistic);
  // }
}
