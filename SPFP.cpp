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

void runRangeUtility(GPOs *baseGPOs, GPOs *cmpGPOs, SPOs *spos){
  SimpleQueries* query = new SimpleQueries(cmpGPOs, spos);

  cout << "------------- Evaluating range utility ---------------" << endl;
  query->checkUtilityRange(query_file, baseGPOs, 50);
  query->checkUtilityRange(query_file, baseGPOs, 100);
  query->checkUtilityRange(query_file, baseGPOs, 200);
  query->checkUtilityRange(query_file, baseGPOs, 400);
  query->checkUtilityRange(query_file, baseGPOs, 800);
}

void runProximityUtility(GPOs *baseGPOs, GPOs *cmpGPOs, SPOs *spos){
  spos->loadKatzScoreFromMemory();

  SimpleQueries* query = new SimpleQueries(cmpGPOs, spos);

  cout << "------------- Evaluating range proximity ---------------" << endl;
  query->checkUtilityProximity(query_file, baseGPOs, 1000, 100);
}

void runUtilities(GPOs *baseGPOs, GPOs *cmpGPOs, SPOs *spos){
  if(is_noise_method && run_utilties){
    runRangeUtility(baseGPOs, cmpGPOs, spos);
    runProximityUtility(baseGPOs, cmpGPOs, spos);
  }
}

void runEBM(GPOs *gpos, SPOs *spos){
  SimpleQueries* query = new SimpleQueries(gpos, spos);

  cout << "----- Precomputing matrices --- " << endl;
  query->buildMatrices(Q);

  cout << "----- Calculating Social Strength --- " << endl;
  query->cacluateSocialStrength();

  for(double i = 1; i < 3; i = i + 0.1){
    cout << "----- Computing accuracy for threshold --- " << i <<endl;
    query->verifySocialStrength(i);
    cout << "--------------------------------------------";
  }
}

void runEBMOnNoised(GPOs *baseGPOs, GPOs *cmpGPOs, SPOs *spos){
  cmpGPOs->countU2UCoOccurrences((uint) TIME_RANGE_IN_SECONDS);
  cmpGPOs->calculateLocationEntropy();
  runEBM(cmpGPOs, spos);

  runUtilities(baseGPOs, cmpGPOs, spos);
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

  runEBMOnNoised(baseGPOs, cmpGPOs, spos);
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

  runEBMOnNoised(baseGPOs, cmpGPOs, spos);
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

  runEBMOnNoised(baseGPOs, cmpGPOs, spos);
}

void locationEntropyNoiseVsEBM(double noise_radius, double entropy_treshold){
  bool preload_LE  = true;
  bool preload_OCC = true;

  GPOs* baseGPOs = loadCheckins(checkins_file, preload_LE, preload_OCC);
  SPOs* spos = loadSocialGraph(graph_file, baseGPOs);

  GPOs* cmpGPOs  = new GPOs();

  cmpGPOs->loadPurturbedLocationsBasedOnLocationEntropy(baseGPOs, noise_radius, entropy_treshold);
  cout << "------------- Locations perturbed Based on  Location Entropy -------------------" << endl;

  runEBMOnNoised(baseGPOs, cmpGPOs, spos);
}

void checkQueryFileStats(){
  cout << "------------- Evaluating query locations ---------------" << endl;

  GPOs* baseGPOs = loadCheckins(checkins_file);
  SPOs* spos = loadSocialGraph(graph_file);

  SimpleQueries* query = new SimpleQueries(baseGPOs, spos);

  query->checkUtilityStats(query_file, 10);
  query->checkUtilityStats(query_file, 25);
  query->checkUtilityStats(query_file, 50);
  query->checkUtilityStats(query_file, 100);
  query->checkUtilityStats(query_file, 200);
  query->checkUtilityStats(query_file, 400);
  query->checkUtilityStats(query_file, 800);
  query->checkUtilityStats(query_file, 1600);
}

int main(int argc, char *argv[]){
  cout.precision(15);

  DELTA_X = ((MAX_X - MIN_X)/ (X-1));
  DELTA_Y = ((MAX_Y - MIN_Y)/ (Y-1));

  if (argc != 9){
    cout << "Usage: " << argv[0] << " graph_file checkins_file query_file [ebm|grouped-ebm|gaussian-v-ebm|snapping-v-ebm|nl-v-ebm|le-v-ebm|occ-hist|compute-katz] run_utilties p1 p2 p3" << endl;
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

  else
    iteration_type = -1;

  run_utilties       = atoi(argv[5]); // [ 0 / 1]
  is_noise_method    = false;

  double p1 = atof(argv[6]); // ( Parameter 1 )
  double p2 = atof(argv[7]); // ( Parameter 2 )
  double p3 = atof(argv[8]); // ( Parameter 3 )

  double noise_radius, group_radius, grid_size_in_km, locality_treshold, entropy_treshold,
         total_parts, part_number, distance_treshold, social_strength_tresh;


  switch(iteration_type){
    case 1:
      cout << "ITRATION: Running EBM" << endl;
      plainEBM();
      break;

    case 2:
      cout << "ITRATION: Running EBM with grouping" << endl;
      group_radius = p2;
      gaussianNoiseVsEBM(0, group_radius);
      is_noise_method = true;
      break;

    case 3:
      cout << "ITRATION: Running EBM with gaussian noise" << endl;
      noise_radius = p1;
      group_radius = p2;

      gaussianNoiseVsEBM(noise_radius, group_radius);
      is_noise_method = true;
      break;

    case 4:
      cout << "ITRATION: Running EBM with grid snapping noise" << endl;
      grid_size_in_km = p1;

      gridSnappingVsEBM(grid_size_in_km);
      is_noise_method = true;
      break;

    case 5:
      cout << "ITRATION: Running EBM with node locality noise" << endl;
      noise_radius = p1;
      locality_treshold = p2;

      nodeLocalityNoiseVsEBM(noise_radius, locality_treshold);
      is_noise_method = true;
      break;

    case 6:
      cout << "ITRATION: Running EBM with location entropy noise noise" << endl;
      noise_radius = p1;
      entropy_treshold = p2;

      locationEntropyNoiseVsEBM(noise_radius, entropy_treshold);
      is_noise_method = true;
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
      GPOs* gpos = loadCheckins(checkins_file);
      SPOs* spos = loadSocialGraph(graph_file);
      spos->computeNodeLocality(gpos);
      spos->writeNodeLocalityToFile();
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
      query->writeHistogramstoFile(social_strength_tresh);
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
