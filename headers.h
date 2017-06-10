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
#include <boost/math/distributions.hpp>
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
extern double RENY_Q;
extern int DATA_SET;

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


#define KATZ_ATTENUATION 0.05
#define KATZ_PATH_LENGTH 2

#define NODE_LOCALITY_BETA 0.1

#define HIJ_SCALE 7
#define HIL_SCALE 7
#define HL_SCALE 0.68


// Number of friendships in gowalla dataset with more than one
// occccurrences
#define RECALL_BOUND 50932

#define LOCATION_NOISE_BOUND 59999999


#define TIME_RANGE_IN_SECONDS 1200 // defines the time difference between 2 checkins to be considered a co-occurrence
#define SPATIAL_RANGE_IN_METERS 25 // defines the spatial distance between 2 checkins to be considered a co-occurrence

// dataset scalability
extern double DATASET_SIZE;
#define X 15000 // table[X][Y]
#define Y 15000
extern double MIN_X;
extern double MIN_Y;
extern double MAX_X;
extern double MAX_Y;

//#define DATASET_SIZE 10000

/// @brief The usual PI/180 constant
#define DEG_TO_RAD 0.017453292519943295769236907684886

/// @brief Earth's quatratic mean radius for WGS-84
#define EARTH_RADIUS_IN_KILOMETERS 6371
#define EARTH_CIRCUMFERENCE 40075.017

extern double DELTA_X ;
extern double DELTA_Y ;

#define SPATIAL_HARD_BOUND 5000
#define TEMPORAL_HARD_BOUND 12 // hours

#define SPATIAL_SOFT_BOUND 1000
#define TEMPORAL_SOFT_BOUND 13.33 // hours

/*
 // DENSE

// Grid Set Up
#define X 300 // table[X][Y]
#define Y 300
#define MIN_X 40
#define MIN_Y -73
#define MAX_X 41
#define MAX_Y -72
#define DELTA_X ((MAX_X - MIN_X)/ (X-1))
#define DELTA_Y ((MAX_Y - MIN_Y)/ (Y-1))

*/

/*
 //SPARSE

// Grid Set Up
#define X 300 // table[X][Y]
#define Y 300
#define MIN_X 30
#define MIN_Y -84
#define MAX_X 50
#define MAX_Y -62
#define DELTA_X ((MAX_X - MIN_X)/ (X-1))
#define DELTA_Y ((MAX_Y - MIN_Y)/ (Y-1))


// Dataset
#define DATASET_SIZE 10000
*/

// Cell type, it is used for the heaps
#define CELL 0
#define RECTANGLE 1
#define POINT 2

// CPM directions
#define UP 0
#define RIGHT 1
#define DOWN 2
#define LEFT 3

// PI
#define PI 3.1428


//#include "SPOs/MemoryMap/pch.h"
//#include "GPOs/MemoryGrid/headers.h"

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
#include "GPOs/MemoryGrid/grid/Cell.h"
#include "GPOs/MemoryGrid/grid/Grid.h"
#include "GPOs/MemoryGrid/GPOs.h"

#include "SPOs/ISPOs.h"

#include "utilities/Group.h"
#include "basicGSQueries/BasicGSQueries.h"
#include "COCqueries/COCqueries.h"
#include "PurturbQ/PurturbQ.h"

#endif
