
#ifndef HEADERSMEMORY_H
#define HEADERSMEMORY_H

#include "headers.h"

// REAL

// Grid Set Up
//scallability
// #define DATASET_SIZE 2000000
// #define X 500 // table[X][Y]
// #define Y 500

//250 km
// #define MAX_X -72.79393586019292
// #define MAX_Y 41.82401035104755
// #define MIN_X -75.30216129645595
// #define MIN_Y 39.579425624113

//64km
// #define MAX_X -73.67130332202277
// #define MAX_Y 41.066482183460764
// #define MIN_X -74.33533012375237
// #define MIN_Y 40.40050397157456

// 100km - official
// #define MAX_X -73.11029266745638
// #define MAX_Y 41.27306291471666
// #define MIN_X -75.41071675619725
// #define MIN_Y 39.36143608785296

// AUSTIN
// #define DATASET_SIZE 5868
// #define X 100 // table[X][Y]
// #define Y 100
// #define MAX_X -97.5564050674
// #define MAX_Y 30.4098181886
// #define MIN_X -97.8997278214
// #define MIN_Y 30.1313942671

//GOWALLA
// #define DATASET_SIZE 1141
// #define X 1000 // table[X][Y]
// #define Y 1000

// #define MAX_X -97.63822406
// #define MAX_Y 30.4087502667
// #define MIN_X -97.8879422167
// #define MIN_Y 30.1541876667


//VLDB Experiments star-graph
//#define DATASET_SIZE 12652
//#define X 3000 // table[X][Y]
//#define Y 1000

//#define MIN_Y -37.8
//#define MIN_X -158
//#define MAX_Y 64
//#define MAX_X 175.0


//#define MIN_Y -37.8
//#define MIN_X -158
//#define MAX_Y 64
//#define MAX_X 175.0

//Experiments NNTopK
//#define DATASET_SIZE 5202
//#define X 1000 // table[X][Y]
//#define Y 1000

//#define MIN_Y 40.5598
//#define MIN_X -74.1994
//#define MAX_Y 40.8996
//#define MAX_X -73.691

// Dataset
//#define DATASET_SIZE 12652

#include "SPOs/MemoryMap/Value.h"
#include "SPOs/MemoryMap/SPOs.h"

//#include "SPOs/MemoryMapWeighted/Pair.h"
#include "SPOs/MemoryMapWeighted/User.h"
#include "SPOs/MemoryMapWeighted/Graph.h"

#include "pTools/KatzScore.h"

#endif
