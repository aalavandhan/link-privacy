#include "headersMemory.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>


#ifndef GLOBALS
#define GLOBALS

double MIN_X = -124.848974;
double MAX_X = -66.885444;
double MIN_Y = 24.396308;
double MAX_Y = 49.384358;

double DATASET_SIZE = 0;
double DELTA_X = 0;
double DELTA_Y = 0;
int MAXSC = 0;                // maximum number of friends plus one
double MAXDIST = 0;          // maximum distance between to points
int MAXT = 0;


// EBM CONSTANTS
double ALPHA = 0.480;
double GAMMA = 0;
double BETA = 0.520;
double Q=0.1; // Order of diversity
// int TIME_RANGE_IN_SECONDS = 1200;

// Global parameters
char *checkins_file, *graph_file, *query_file;
int iteration_type = 1;
bool run_utilties = false;
bool is_noise_method = false;

double noise_radius, group_radius, grid_size_in_km, locality_treshold, entropy_treshold,
         total_parts, part_number, distance_treshold, social_strength_tresh, combination_type,
         function_type, p1, p2, p3, p4;
#endif

GPOs* loadCheckins(char* checkins_file){
  cout << "------------- Loading checkins ---------------" << endl;
  GPOs* gpos = new GPOs(checkins_file);
  return gpos;
}

GPOs* loadCheckins(char* checkins_file, bool preload_LE, bool preload_OCC){
  GPOs* gpos = loadCheckins(checkins_file);

  if(preload_LE)
    gpos->calculateLocationEntropy();

  if(preload_OCC)
    gpos->countU2UCoOccurrences((uint) TIME_RANGE_IN_SECONDS);

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

void runRangeUtility(GPOs *purturbedGPOs, GPOs *cmpGPOs, SPOs *spos){
  SimpleQueries* query = new SimpleQueries(cmpGPOs, spos);

  cout << "------------- Evaluating range utility ---------------" << endl;

  if(noise_radius < 50)
    query->checkUtilityRange(query_file, purturbedGPOs, 50, noise_radius);

  if(noise_radius < 100)
    query->checkUtilityRange(query_file, purturbedGPOs, 100, noise_radius);

  if(noise_radius < 200)
    query->checkUtilityRange(query_file, purturbedGPOs, 200, noise_radius);

  if(noise_radius < 400)
    query->checkUtilityRange(query_file, purturbedGPOs, 400, noise_radius);

  if(noise_radius < 800)
    query->checkUtilityRange(query_file, purturbedGPOs, 800, noise_radius);

  if(noise_radius < 1600)
    query->checkUtilityRange(query_file, purturbedGPOs, 1600, noise_radius);
}

void runProximityUtility(GPOs *purturbedGPOs, GPOs *cmpGPOs, SPOs *spos){
  spos->loadKatzScoreFromMemory();

  SimpleQueries* query = new SimpleQueries(cmpGPOs, spos);

  cout << "------------- Evaluating range proximity ---------------" << endl;
  query->checkUtilityProximity(query_file, purturbedGPOs, 1000, 500, noise_radius);
}

void runUtilities(GPOs *purturbedGPOs, GPOs *cmpGPOs, SPOs *spos){
  runRangeUtility(purturbedGPOs, cmpGPOs, spos);
  runProximityUtility(purturbedGPOs, cmpGPOs, spos);
}

void runEBM(GPOs *gpos, SPOs *spos){
  SimpleQueries* query = new SimpleQueries(gpos, spos);

  cout << "----- Precomputing matrices --- " << endl;
  query->buildMatrices(Q);

  cout << "----- Calculating Social Strength --- " << endl;
  query->cacluateSocialStrength();

  for(double i = 0.3; i < 2.5; i = i + 0.1){
    cout << "----- Computing accuracy for threshold --- " << i <<endl;
    query->verifySocialStrength(i);
    cout << "--------------------------------------------";
  }
}

void runBasicOnNoised(GPOs *baseGPOs, GPOs *purturbedGPOs, GPOs *cmpGPOs, SPOs *spos){
  // cmpGPOs->countU2UCoOccurrences((uint) TIME_RANGE_IN_SECONDS);
  // cmpGPOs->calculateLocationEntropy();
  // runEBM(cmpGPOs, spos);

  //grab cooccurence Matrix from baseGPOs
  //grab cooccurence matrix from cmpGPOs
  //compare them both.
  map<int, pair<double,double>> user_to_precision_recall;

  auto cooccurence_matrix_base = baseGPOs->getInsignificantCooccurrenceMatrix();
  auto cooccurence_matrix_cmp = cmpGPOs->getInsignificantCooccurrenceMatrix();

  unordered_set<int> seen_users_in_base;

  for(auto c_it = cooccurence_matrix_base->begin(); c_it != cooccurence_matrix_base->end(); c_it++){
    int user_1 = c_it->first;
    seen_users_in_base.insert(user_1);
    auto users_location_frequency_map_base = c_it->second;
    map<int, vector<pair<int, int> >* >* users_location_frequency_map_cmp = NULL;


    auto iter_outer = cooccurence_matrix_cmp->find(user_1);
    if(iter_outer != cooccurence_matrix_cmp->end()){
      users_location_frequency_map_cmp = iter_outer->second;
      int cooccurrences_across_users_base = 0;
      int cooccurrences_across_users_cmp = 0;
      int cooccurrences_across_users_min = 0;
      for(auto ulh_it = users_location_frequency_map_base->begin(); ulh_it != users_location_frequency_map_base->end(); ulh_it++){
        int user_2 = ulh_it->first;
        vector<pair<int, int>>* cooccurrence_counts_vector_base = ulh_it->second;
        vector<pair<int, int>>* cooccurrence_counts_vector_cmp = NULL;


        int cooccurrence_count_cmp = 0;
        int cooccurrence_count_base = 0;
        auto iter_inner = users_location_frequency_map_cmp->find(user_2);
        if(iter_inner != users_location_frequency_map_cmp->end()){
          cooccurrence_counts_vector_cmp = iter_inner->second;

          for(auto l_it=cooccurrence_counts_vector_cmp->begin(); l_it != cooccurrence_counts_vector_cmp->end(); l_it++){
          int cooccrences_at_l = l_it->second;
          cooccurrence_count_cmp += cooccrences_at_l;
          }
        }else{
          // do nothing because cooccurences are zero
          // cout<<"ERROR ? baseGPOs and cmpGPOs do not have similar cooccurrence_counts_vector"<<endl;
        }

        for(auto l_it=cooccurrence_counts_vector_base->begin(); l_it != cooccurrence_counts_vector_base->end(); l_it++){
          int cooccrences_at_l = l_it->second;
          cooccurrence_count_base += cooccrences_at_l;
        }
        cooccurrences_across_users_base += cooccurrence_count_base;
        cooccurrences_across_users_cmp += cooccurrence_count_cmp;

        cooccurrences_across_users_min += min(cooccurrence_count_cmp,cooccurrence_count_base);
      }
      double precision = cooccurrences_across_users_min/(double)cooccurrences_across_users_cmp;
      double recall = cooccurrences_across_users_min/(double)cooccurrences_across_users_base;
      if(recall > 1) {recall=1;}
      user_to_precision_recall.insert(make_pair(user_1,make_pair(precision,recall)));
    }else{
      //case when user is not there in cmp cooccurrence matrix
      user_to_precision_recall.insert(make_pair(user_1,make_pair(0,0)));
      // cout<<"ERROR ? baseGPOs and cmpGPOs do not have similar cooccurence matrices"<<endl;
    }
  }

  for(auto c_it = cooccurence_matrix_base->begin(); c_it != cooccurence_matrix_base->end(); c_it++){
    int user = c_it->first;
    if(seen_users_in_base.find(user) == seen_users_in_base.end()){
      user_to_precision_recall.insert(make_pair(user,make_pair(0,-1)));
    }
  }

  //printing--------
  ofstream outfile;
  ostringstream ss;
  ss<<"Basic" << p1<<"_"<<p2<<" "<<p3<<"_"<<p4;
  outfile.open(ss.str());
  for(auto it = user_to_precision_recall.begin(); it != user_to_precision_recall.end(); it++){
    auto precision_recall_pair = it->second;
    double precision = precision_recall_pair.first;
    double recall = precision_recall_pair.second;
    outfile<<it->first<< " "<<precision<<" "<<recall<<endl;
  }
  outfile.close();
  //------------------
}

void runEBMOnNoised(GPOs *baseGPOs, GPOs *purturbedGPOs, GPOs *cmpGPOs, SPOs *spos){
  cmpGPOs->countU2UCoOccurrences((uint) TIME_RANGE_IN_SECONDS);
  cmpGPOs->calculateLocationEntropy();
  runEBM(cmpGPOs, spos);
  runBasicOnNoised(baseGPOs, purturbedGPOs, cmpGPOs, spos);
  runUtilities(purturbedGPOs, cmpGPOs, spos);
}

void plainEBM(){
  bool preload_LE  = true;
  bool preload_OCC = true;

  GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, gpos);

  runEBM(gpos, spos);
}

void gridSnappingVsEBM(double grid_size_in_km){
  bool preload_LE  = false;
  bool preload_OCC = true;

  GPOs* baseGPOs   = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  GPOs* cmpGPOs    = new GPOs();
  cmpGPOs->createNewGPOsbyGridSnapping(baseGPOs, grid_size_in_km);
  cout << "Number of locations loaded " << cmpGPOs->locations.size() << endl;
  cout << "------------- Noise added -------------------" << endl;

  runEBMOnNoised(baseGPOs, baseGPOs, cmpGPOs, spos);
}

void gaussianNoiseVsEBM(double noise_radius, double group_radius){
  bool preload_LE  = false;
  bool preload_OCC = true;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  GPOs* purturbedGPOs = new GPOs();
  GPOs* cmpGPOs  = new GPOs();

  purturbedGPOs->loadPurturbedLocations(baseGPOs, noise_radius);
  cout << "------------- Locations perturbed -------------------" << endl;

  cmpGPOs->groupLocationsByRange(purturbedGPOs, group_radius, false);
  cout << "------------- Locations Grouped -------------------" << endl;

  runEBMOnNoised(baseGPOs, purturbedGPOs, cmpGPOs, spos);
}

void CombinationNoiseVsEBM(double noise_radius){
  bool preload_LE  = false;
  bool preload_OCC = true;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  cout << "------------- Load computed checkin locality ---------------" << endl;
  SPOs *tmp_spos = new SPOs();
  map< int, map<int, pair<int,double> >* >* user_to_order_to_location_locality = tmp_spos->loadCheckinLocalityFromFile();

  GPOs* purturbedGPOs = new GPOs();
  GPOs* cmpGPOs  = new GPOs();

  if(combination_type == 0){
    purturbedGPOs->loadPurturbedLocationsBasedOnCombinationFunction(baseGPOs, user_to_order_to_location_locality, noise_radius, false, function_type);
  } else if(combination_type == 1) {
    purturbedGPOs->loadPurturbedLocationsBasedOnCombinationFunction(baseGPOs, user_to_order_to_location_locality, noise_radius, true, function_type);
  } else if(combination_type == 2) {
    purturbedGPOs->loadPurturbedLocationsBasedOnCombinationFunctionofCOOCC(baseGPOs, baseGPOs->getL2U2COOCC(), noise_radius, false, function_type);
  } else if(combination_type == 3) {
    purturbedGPOs->loadPurturbedLocationsBasedOnCombinationFunctionofCOOCC(baseGPOs, baseGPOs->getL2U2COOCC(), noise_radius, true, function_type);
  }

  cout << "------------- Locations perturbed -------------------" << endl;
  purturbedGPOs->countU2UCoOccurrences((uint) TIME_RANGE_IN_SECONDS);

  cmpGPOs->groupLocationsByRange(purturbedGPOs, group_radius, false);
  cout << "------------- Locations Grouped -------------------" << endl;

  runEBMOnNoised(baseGPOs, purturbedGPOs, cmpGPOs, spos);
}


void nodeLocalityNoiseVsEBM(double noise_radius, double locality_treshold){
  bool preload_LE  = false;
  bool preload_OCC = true;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  GPOs* cmpGPOs  = new GPOs();

  cout << "------------- Load computed node locality ---------------" << endl;
  SPOs *tmp_spos = new SPOs();
  map< int, double >* node_locality = tmp_spos->loadNodeLocalityFromFile();

  cmpGPOs->loadPurturbedLocationsBasedOnNodeLocality(baseGPOs, node_locality, noise_radius, locality_treshold);
  cout << "------------- Locations perturbed Based on  Node locality -------------------" << endl;

  runEBMOnNoised(baseGPOs, baseGPOs, cmpGPOs, spos);
}

void locationEntropyNoiseVsEBM(double noise_radius, double entropy_treshold){
  bool preload_LE  = true;
  bool preload_OCC = true;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  GPOs* cmpGPOs  = new GPOs();

  cmpGPOs->loadPurturbedLocationsBasedOnLocationEntropy(baseGPOs, noise_radius, entropy_treshold);
  cout << "------------- Locations perturbed Based on  Location Entropy -------------------" << endl;

  runEBMOnNoised(baseGPOs, baseGPOs, cmpGPOs, spos);
}

void checkQueryFileStats(){
  cout << "------------- Evaluating query locations ---------------" << endl;

  GPOs* baseGPOs = loadCheckins(checkins_file);
  SPOs* spos = loadSocialGraph(graph_file);

  SimpleQueries* query = new SimpleQueries(baseGPOs, spos);

  query->checkUtilityStats(query_file, 10, 5);
  query->checkUtilityStats(query_file, 25, 12.5);
  query->checkUtilityStats(query_file, 50, 25);
  query->checkUtilityStats(query_file, 100, 50);
  query->checkUtilityStats(query_file, 200, 100);
  query->checkUtilityStats(query_file, 400, 200);
  query->checkUtilityStats(query_file, 800, 400);
  query->checkUtilityStats(query_file, 1600, 800);
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
// void test(){
//   Utilities util;
//   pair<double,double> point;
//   for(int i = 0; i < 2000; i++){
//      point = util.addGaussianNoise(0,0,200);
//      cout << util.computeMinimumDistance(point.first, point.second);
//   }
// }

int main(int argc, char *argv[]){
  cout.precision(15);

  testpTools();

  DELTA_X = ((MAX_X - MIN_X)/ (X-1));
  DELTA_Y = ((MAX_Y - MIN_Y)/ (Y-1));

  if (argc != 10){
    cout << "Usage: " << argv[0] << " graph_file checkins_file query_file [ebm|grouped-ebm|gaussian-v-ebm|snapping-v-ebm|nl-v-ebm|le-v-ebm|occ-hist|compute-katz] run_utilties p1 p2 p3 p4" << endl;
    return -1;
  }

  graph_file    = argv[1];
  checkins_file = argv[2];
  query_file    = argv[3];


  if (strcmp(argv[4], "ebm") == 0)
    iteration_type = 1;
  else if (strcmp(argv[4], "grouped-ebm") == 0)
    iteration_type = 2;
  else if (strcmp(argv[4], "gaussian-v-ebm") == 0)
    iteration_type = 3;
  else if (strcmp(argv[4], "snapping-v-ebm") == 0)
    iteration_type = 4;
  else if (strcmp(argv[4], "nl-v-ebm") == 0)
    iteration_type = 5;
  else if (strcmp(argv[4], "le-v-ebm") == 0)
    iteration_type = 6;
  else if (strcmp(argv[4], "occ-hist") == 0)
    iteration_type = 7;
  else if (strcmp(argv[4], "compute-katz") == 0)
    iteration_type = 8;
  else if (strcmp(argv[4], "compute-node-locality") == 0)
    iteration_type = 9;
  else if (strcmp(argv[4], "compute-histograms") == 0)
    iteration_type = 10;
  else if (strcmp(argv[4], "query-metrics") == 0)
    iteration_type = 11;
  else if (strcmp(argv[4], "comb-v-ebm") == 0)
    iteration_type = 12;
  else
    iteration_type = -1;

  run_utilties       = atoi(argv[5]); // [ 0 / 1]

  p1 = atof(argv[6]); // ( Parameter 1 )
  p2 = atof(argv[7]); // ( Parameter 2 )
  p3 = atof(argv[8]); // ( Parameter 3 )
  p4 = atof(argv[9]); // ( Parameter 4 )


  switch(iteration_type){
    case 1:
      cout << "ITRATION: Running EBM" << endl;
      plainEBM();
      break;

    case 2:
      cout << "ITRATION: Running EBM with grouping" << endl;
      noise_radius = 0;
      group_radius = p2;
      gaussianNoiseVsEBM(noise_radius, group_radius);
      break;

    case 3:
      cout << "ITRATION: Running EBM with gaussian noise" << endl;
      noise_radius = p1;
      group_radius = p2;

      gaussianNoiseVsEBM(noise_radius, group_radius);
      break;

    case 4:
      cout << "ITRATION: Running EBM with grid snapping noise" << endl;
      grid_size_in_km = p1;
      noise_radius = grid_size_in_km * 1000;

      gridSnappingVsEBM(grid_size_in_km);
      break;

    case 5:
      cout << "ITRATION: Running EBM with node locality noise" << endl;
      noise_radius = p1;
      locality_treshold = p2;

      nodeLocalityNoiseVsEBM(noise_radius, locality_treshold);
      break;

    case 6:
      cout << "ITRATION: Running EBM with location entropy noise noise" << endl;
      noise_radius = p1;
      entropy_treshold = p2;

      locationEntropyNoiseVsEBM(noise_radius, entropy_treshold);
      break;

    case 7:{
      cout << "METRICS: Calculating Cooccurrence distribution" << endl;
      bool preload_LE  = true;
      bool preload_OCC = true;

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = loadSocialGraph(graph_file);
      spos->loadNodeLocalityFromFile();

      SimpleQueries* query = new SimpleQueries(gpos, spos);

      query->cacluateCooccurrenceDistributionBasedOnNodeLocality();
      query->cacluateCooccurrenceDistributionBasedOnLocationEntropy();
      break;
    }

    case 8:{
      cout << "METRICS: Pre-Compute KATZ score" << endl;
      GPOs* gpos = loadCheckins(checkins_file);
      SPOs* spos = loadSocialGraph(graph_file);
      total_parts = p1;
      part_number = p2;
      distance_treshold = p3;
      spos->precomputeKatzScore(gpos, total_parts, part_number, distance_treshold);
      break;
    }

    case 9:{
      cout << "METRICS: Pre-Compute node_locality" << endl;
      bool preload_LE  = false;
      bool preload_OCC = true;

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = loadSocialGraph(graph_file, gpos);

      spos->computeCheckinLocalityMap(gpos);
      spos->writeCheckinLocalityToFile();
      // cout << "------------- Computing mean distance between friend pairs ---------------" << endl;
      // cout << "Mean distance between all pairs of friends :" << spos->computeMeanDistanceBetweenAllFriends(gpos) << endl;

      // spos->computeNodeLocality(gpos);

      break;
    }

    case 10:{
      cout << "METRICS: Compute histograms" << endl;
      bool preload_LE  = true;
      bool preload_OCC = true;

      GPOs* gpos = loadCheckins(checkins_file, preload_LE, preload_OCC);
      SPOs* spos = loadSocialGraph(graph_file, gpos);

      social_strength_tresh = p1;

      SimpleQueries* query = new SimpleQueries(gpos, spos);
      query->buildMatrices(Q);
      query->cacluateSocialStrength();

      cout << "Using Threshold" << social_strength_tresh << endl;

      query->writeHistogramstoFile(social_strength_tresh);
      spos->writeUserFriendsToFile();
      break;
    }

    case 11:{
      checkQueryFileStats();
    }

    case 12:
      cout << "ITRATION: Running EBM with combination noise" << endl;
      noise_radius     = p1;
      group_radius     = p2;
      combination_type = p3;
      function_type    = p4;
      CombinationNoiseVsEBM(noise_radius);
      break;

    default:
      cout << "Invalid iteration type; Enter a valid option" << endl;
  }

  // Plain EBM


  // cout << "------------- Loading checkins ---------------" << endl;
  // GPOs* gpos = new GPOs(argv[2]);
  // cout << "------------- Loading complete ---------------" << endl;

  // cout << "------------- Loading checkins ---------------" << endl;
  // GPOs* g = new GPOs(argv[2]);
  // g->countU2UCoOccurrences((uint) time_range_in_seconds);
  // GPOs* gpos = new GPOs(argv[2]);
  // gpos->groupLocationsByRange(g, r2, isOptimistic);
  // cout << "------------- Loading complete ---------------" << endl;

  // cout << "------------- Loading checkins ---------------" << endl;
  // GPOs* g = new GPOs(argv[2]);
  // g->countU2UCoOccurrences((uint) time_range_in_seconds);
  // cout << "------------- Loading complete ---------------" << endl;


  // testing grid snapping
  // GPOs* gpos = new GPOs();
  // gpos->createNewGPOsbyGridSnapping(g, r1); //second var is the x distance of a cell is km
  // cout << "Number of locations loaded " << gpos->locations.size() << endl;
  // cout << "------------- Noise added -------------------" << endl;


  // test Gaussian noise with grouping
  // GPOs* gg = new GPOs();
  // GPOs* gpos = new GPOs();
  // gg->loadPurturbedLocations(g, r1);
  // cout << "------------- Locations perturbed -------------------" << endl;
  // gpos->groupLocationsByRange(gg, r2, isOptimistic);
  // cout << "------------- Locations Grouped -------------------" << endl;

  // test Gaussian for high node locality  noise without grouping
  // cout << "------------- Load computed node locality ---------------" << endl;
  // SPOs* tmp_spos = new SPOs();
  // map< int, double >* node_locality = tmp_spos->loadNodeLocalityFromFile();
  // GPOs* gpos = new GPOs();
  // gpos->loadPurturbedLocationsBasedOnNodeLocality(g, node_locality, r1, 0.75);
  // gpos->generateFrequencyCache();
  // cout << "------------- Locations perturbed Based on  Node locality -------------------" << endl;

  // test Gaussian for high location entropy noise without grouping
  // GPOs* gpos = new GPOs();
  // g->calculateLocationEntropy();
  // gpos->loadPurturbedLocationsBasedOnLocationEntropy(g, r1, 1);
  // gpos->generateFrequencyCache();
  // cout << "------------- Locations perturbed Based on  Location Entropy -------------------" << endl;


  // test gaussian noise without grouping
  // GPOs* gpos = new GPOs();
  // gpos->loadPurturbedLocations(g, r1);
  // cout << "------------- Locations perturbed -------------------" << endl;

  // cout << "----- Loading Cooccurrence Matrix --- " << endl;
  // gpos->countU2UCoOccurrences((uint) time_range_in_seconds);
  // cout << "----- Completed Loading Cooccurrence Matrix --- " << endl;


  // cout << "----- Calculating Location Entropy --- " << endl;
  // gpos->calculateLocationEntropy();
  // cout << "----- Completed Calculating Location Entropy --- " << endl;

  // // Loading social network and checkins
  // cout << "------------- Loading SocialGraph ---------------" << endl;
  // SPOs* spos = new SPOs();
  // spos->load(argv[1]);

  // cout << "------------- Loading SocialGraph ---------------" << endl;
  // SPOs* spos = new SPOs(g);
  // spos->load(argv[1]);


  // cout << "------------- Loading SocialGraph ---------------" << endl;
  // SPOs* spos = new SPOs(gpos);
  // spos->load(argv[1]);

  // spos->loadKatzScoreFromMemory();

  // spos->precomputeKatzScore(gpos, r1, r2, time_range_in_seconds);

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
  //   GPOs* purturbedGPOs = new GPOs();
  //   purturbedGPOs->loadPurturbedLocations(gpos, r1);
  //   cout << "------------- Noise added -------------------" << endl;

  //   count_cooccurences(spos, purturbedGPOs, r2, tR, isOptimistic);
  // }
}
