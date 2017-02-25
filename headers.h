#ifndef HEADER_H
#define HEADER_H


#include <iostream>
#include <sys/time.h>
#include <vector>
#include <queue>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <cmath>

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <list>
#include <set>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>

#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time_io.hpp"
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/random.hpp>
#include <boost/math/distributions/normal.hpp>

#include <stdint.h>
#include <stdarg.h>

#include <errno.h>
#include <bitset>
// #include "boost/algorithm/string.hpp"


// --------------------- Server

//#include <sys/types.h>
//#include <sys/socket.h>
//#include <winsock2.h>
//#include <netinet/in.h>

// --------------------- Server
#define EARTH_CIRCUMFERENCE 40075.017

using namespace std;

//VLDB dataset
//#define MAXGROUPSIZE 12652
//#define MAXSC 420                  // maximum number of friends plus one
//#define MAXDIST 52.1636             // maximum distance between to points
//#define MAXDIST 0.021636

//gowalla dataset

// #define MAXGROUPSIZE 1141
// #define MAXSC 71                // maximum number of friends plus one
// #define MAXDIST 39.352          // maximum distance between to points

//austin dataset from Nick
#define MAXGROUPSIZE 5868
extern int MAXSC;                // maximum number of friends plus one
extern double MAXDIST;          // maximum distance between to points
extern int MAXT;
extern double ALPHA;
extern double BETA;
extern double GAMMA;

//scallability

// d maxS maxT
// 10 1239 570
// 30 2011 4311
// 50 2931 14354
// 70 3215 21460
// 90 3695 34098


// Scalability
// #define MAXGROUPSIZE 1500000

//d=10
// #define MAXSC 1240
// #define MAXT 570
//d=30
// #define MAXSC 2012
// #define MAXT 4311
//d=50
// #define MAXSC 2932
// #define MAXT 14354
//d=70
// #define MAXSC 3216
// #define MAXT 21460
//d=90
// #define MAXSC 3696
// #define MAXT 34098

// #define MAXDIST 0.89831527707  //correspons to 100km


#define EULERS_NUMBER 2.718281828459045
#define BOUNDARY_ERROR 0.00000000000001


#define KATZ_ATTENUATION 0.1
#define KATZ_PATH_LENGTH 2

#define NODE_LOCALITY_BETA 384


// Number of friendships in gowalla dataset with more than one
// occccurrences
#define RECALL_BOUND 10505

#define LOCATION_NOISE_BOUND 9999999


#define TIME_RANGE_IN_SECONDS 1200 // defines the time difference between 2 checkins to be considered a co-occurrence

#include "utilities/res_point.h"
#include "utilities/Utilities.h"
#include "utilities/my_pair.h"
#include "utilities/ranked_pair.h"
#include "utilities/PairHasher.h"

#include "pTools/MIToolbox.h"
#include "pTools/ArrayOperations.h"
#include "pTools/CalculateProbability.h"
#include "pTools/Entropy.h"
#include "pTools/RenyiEntropy.h"
#include "pTools/WeightedEntropy.h"

#include "GPOs/MemoryGrid/grid/Point.h"

#include "GPOs/IGPOs.h"
#include "SPOs/ISPOs.h"

#include "utilities/Group.h"
#include "basicGSQueries/BasicGSQueries.h"

#endif