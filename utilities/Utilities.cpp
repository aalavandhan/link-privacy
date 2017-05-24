#include "../headersMemory.h"
#include <cmath>


template<class T>
struct distribution
{ // general case, assuming T is of integral type
  typedef boost::uniform_int<> type;
};

template<>
struct distribution<float>
{ // float case
  typedef boost::uniform_real<> type;
};

template<>
struct distribution<double>
{ // double case
  typedef boost::uniform_real<> type;
};
template <typename N> N getRandom(N min, N max, N sd)
{
  typedef typename distribution<N>::type distro_type;
  boost::mt19937 seed( (int) sd );
  distro_type dist(min,max);
  boost::variate_generator<boost::mt19937&, distro_type > random(seed, dist);
  return random();
};

Utilities::Utilities(){
  srand((unsigned)time(NULL));
}

Utilities::~Utilities(){}

double Utilities::distanceBetween(double lat1, double lon1, double lat2, double lon2){
  double lat1r, lon1r, lat2r, lon2r, u, v;
  lat1r = DEG_TO_RAD * lat1;
  lon1r = DEG_TO_RAD * lon1;
  lat2r = DEG_TO_RAD * lat2;
  lon2r = DEG_TO_RAD * lon2;
  u = sin((lat2r - lat1r)/2);
  v = sin((lon2r - lon1r)/2);
  return 2.0 * EARTH_RADIUS_IN_KILOMETERS * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

double Utilities::angleFromCoordinate(double lat1, double long1, double lat2, double long2) {
  double dLon = (long2 - long1);
  double y = sin(dLon) * cos(lat2);
  double x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon);
  double brng = atan2(y, x);
  brng = brng * (1/DEG_TO_RAD);
  brng = ((int)brng + 360) % 360;
  brng = 360 - brng; // count degrees counter-clockwise - remove to make clockwise
  return brng;
}

pair<double,double> Utilities::addNoise(double x, double y, double radius){
  static unsigned int seed = rand() % 10000;
  ++seed;

  int R = EARTH_RADIUS_IN_KILOMETERS * 1000;
  int direction = getRandom<int>(0,360, seed);

  double nLat = y  + (radius * cos( direction ) / R) * (1/DEG_TO_RAD);
  double nLon = x  + (radius * sin( direction ) / R) * (1/DEG_TO_RAD) / cos(y * DEG_TO_RAD);

  // cout << lat << " " << lon << " " << nLat << " " << nLon << endl;
  // double newDirection = angleFromCoordinate(nLat, nLon, lat, lon);
  // cout << distanceBetween(nLat, nLon, y, x) * 1000 << " " << abs(noise_distance)  << " " << direction << " " << newDirection <<endl;

  return make_pair(nLon, nLat);
}

pair<double,double> Utilities::addGaussianNoise(double x, double y, double radius){
  return addGaussianNoise(x, y, radius, 0);
}

// Noise radius in meters
pair<double,double> Utilities::addGaussianNoise(double x, double y, double radius, double offset){
  boost::mt19937 rng;
  static unsigned int seed = rand() % 10000;
  rng.seed(++seed);

  boost::normal_distribution<> nd(0.0, radius/2);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_norormal(rng, nd);

  double lat=0,lon=0, nLat=0, nLon=0, noise_distance=0;
  int R = EARTH_RADIUS_IN_KILOMETERS * 1000;

  noise_distance = abs(var_norormal()) + offset;
  int direction = getRandom<int>(0,360, seed);

  lon = x;
  lat = y;

  // nLat = lat  + noise_distance * cos( direction );
  // nLon = lon  + noise_distance * sin( direction );

  nLat = lat  + (noise_distance * cos( direction ) / R) * (1/DEG_TO_RAD);
  nLon = lon  + (noise_distance * sin( direction ) / R) * (1/DEG_TO_RAD) / cos(lat * DEG_TO_RAD);

  // cout << lat << " " << lon << " " << nLat << " " << nLon << endl;
  // double newDirection = angleFromCoordinate(nLat, nLon, lat, lon);
  // cout << distanceBetween(nLat, nLon, y, x) * 1000 << " " << abs(noise_distance)  << " " << direction << " " << newDirection <<endl;

  return make_pair(nLon, nLat);
}

boost::posix_time::ptime Utilities::addTemporalGaussianNoise(boost::posix_time::ptime time, uint deviation_in_seconds){
  return addTemporalGaussianNoise(time, deviation_in_seconds, 0);
}

boost::posix_time::ptime Utilities::addTemporalGaussianNoise(boost::posix_time::ptime time, uint deviation_in_seconds, int offset){
  boost::mt19937 rng;
  static unsigned int seed = rand() % 10000;
  rng.seed(++seed);

  boost::normal_distribution<> nd(0.0, deviation_in_seconds/2);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_norormal(rng, nd);

  int noise_in_seconds = (int) var_norormal() + offset;

  return ( time + boost::posix_time::seconds(noise_in_seconds) );
}

boost::posix_time::ptime Utilities::addTemporalNoise(boost::posix_time::ptime time, int deviation_in_seconds){
  static unsigned int seed = rand() % 10000;
  ++seed;

  int direction = getRandom<int>(0, 1, seed);

  if(direction == 0){
    return ( time + boost::posix_time::seconds(deviation_in_seconds) );
  } else {
    return ( time - boost::posix_time::seconds(deviation_in_seconds) );
  }
}


//change
//for Diversity between topk results
double Utilities::getDistanceBetween(int friendsO[], int friendsK[], int friendsSizeO, int friendsSizeK, string f){

	//if(f == "A"){
		// Jaccard Similarity
		return 1 - computeJaccard( friendsO, friendsK, friendsSizeO, friendsSizeK );
	//}
	/*
	else if(f == "B"){
		//Normalized by min
		if(friendsSizeK != 0 && friendsSizeO != 0)
			return countIntersection( friendsO, friendsK, friendsSizeO, friendsSizeK ) / (double) (friendsSizeO > friendsSizeK ? friendsSizeK : friendsSizeO) ;
		else
			return 0;
	}
	else if(f == "C"){
		// Exponentially smoothed - need to experiment with accumulating exponential smoothing
		return (exp(countIntersection( friendsO, friendsK, friendsSizeO, friendsSizeK ) / (double) (friendsSizeO > friendsSizeK ? friendsSizeK : friendsSizeO))) - 1 ;
	}
	*/
}

double Utilities::getSocialDistanceBetween( unordered_set<int>* _f1, unordered_set<int>* _f2 ){
	int intersections = 0, union_size = 0;
	if(_f1 != NULL && _f2 != NULL){
		if(_f1->size() == 0 || _f2->size() == 0)
			return 0;
		if(_f1->size() > _f2->size() ) {
			for(auto it = _f2->begin();it!=_f2->end();it++){
				int id = *it;
				if(_f1->find(id) != _f1->end()){
					intersections++;
				}
			}
		}
		else{
			for(auto it = _f1->begin();it!=_f1->end();it++){
				int id = *it;
				if(_f2->find(id) != _f2->end()){
					intersections++;
				}
			}
		}
		union_size = _f2->size() + _f1->size() - intersections;
		return (union_size - intersections) / (double) union_size;
	}
	return 0;
}

double Utilities::getTextDistanceBetween( unordered_set<string>* user_profile_1, unordered_set<string>* user_profile_2 ){
	int intersections = 0, union_size = 0;
	if(user_profile_1 != NULL && user_profile_2 != NULL){
		if(user_profile_1->size() == 0 || user_profile_2->size() == 0)
			return 0;
		if(user_profile_1->size() > user_profile_2->size() ) {
			for(auto it = user_profile_2->begin();it!=user_profile_2->end();it++){
				string word = *it;
				if(user_profile_1->find(word) != user_profile_1->end()){
					intersections++;
				}
			}
		}
		else{
			for(auto it = user_profile_1->begin(); it != user_profile_1->end() ;it++){
				string word = *it;
				if(user_profile_2->find(word) != user_profile_2->end()){
					intersections++;
				}
			}
		}
		union_size = user_profile_2->size() + user_profile_1->size() - intersections;
		return (union_size - intersections) / (double) union_size;
	}
	return 0;
}



double Utilities::computeSetIntersection( unordered_set<int>* _f1, unordered_set<int>* _f2 ){
	int intersections = 0;
	if(_f1 != NULL && _f2 != NULL){
		if(_f1->size() == 0 || _f2->size() == 0)
			return 0;
		if(_f1->size() > _f2->size() ) {
			for(auto it = _f2->begin();it!=_f2->end();it++){
				int id = *it;
				if(_f1->find(id) != _f1->end()){
					intersections++;
				}
			}
		}
		else{
			for(auto it = _f1->begin();it!=_f1->end();it++){
				int id = *it;
				if(_f2->find(id) != _f2->end()){
					intersections++;
				}
			}
		}

		return intersections;
	}
	return 0;
}


double Utilities::getTextDistanceBetween( unordered_map<string, double>* user_profile_1, unordered_set<string>* user_profile_2 ){
	int intersections = 0, union_size = 0;
	if(user_profile_1 != NULL && user_profile_2 != NULL){
		if(user_profile_1->size() == 0 || user_profile_2->size() == 0)
			return 0;
		if(user_profile_1->size() > user_profile_2->size() ) {
			for(auto it = user_profile_2->begin();it!=user_profile_2->end();it++){
				string word = *it;
				if(user_profile_1->find(word) != user_profile_1->end()){
					intersections++;
				}
			}
		}
		else{
			for(auto it = user_profile_1->begin(); it != user_profile_1->end() ;it++){
				string word = it->first;
				if(user_profile_2->find(word) != user_profile_2->end()){
					intersections++;
				}
			}
		}
		union_size = user_profile_2->size() + user_profile_1->size() - intersections;
		return (union_size - intersections) / (double) union_size;
	}

	return 0;
}



double Utilities::computeJaccard(int arr1[], int arr2[], int m, int n){
	if(m==0 || n==0){
		return 1;
	}

	int intersections = countIntersection(arr1, arr2, m, n);
	int union_size = m + n - intersections;

	// cout<<"Array 1: ";
	// for(int i = 0; i < m; i++)cout<<" "<<arr1[i];
	// cout<<" | Array 2: ";
	// for(int i = 0; i < n; i++)cout<<" "<<arr2[i];
	// cout<<" intersections: "<<intersections<<" union_size: "<<union_size<<endl;

	return  intersections/(double) union_size;
}

/* Function prints union of arr1[] and arr2[]
   m is the number of elements in arr1[]
   n is the number of elements in arr2[] */
int Utilities::countUnion(int arr1[], int arr2[], int m, int n){
  int count = 0;
  int i = 0, j = 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j]){
	++count;
     ++i;
	 }
    else if (arr2[j] < arr1[i]){
      ++count;
	  ++j;
	  }
    else
    {
      ++count;
	  ++j;
	  ++i;
    }
  }

  /* Print remaining elements of the larger array */
  while(i < m){
  ++count;
     ++i;
	}
  while(j < n){
   ++count;
	  ++j;
  }
  return count;
}



int Utilities::countIntersection(int arr1[], int arr2[], int m, int n)
{
	int count=0;

  int i = 0, j = 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j])
      i++;
    else if (arr2[j] < arr1[i])
      j++;
    else /* if arr1[i] == arr2[j] */
    {
	++count;
      ++j;
      i++;
    }
  }

  return count;
}

int Utilities::getCooccurrencesWithinTimeBlock(vector<pair<uint, int>>* arr1, vector<pair<uint, int>>* arr2, uint time_block, set< pair<int,int> >* orders){
  int count=0;
  int m = arr1->size();
  int n = arr2->size();
  int i = 0, j = 0;

  while (i < m && j < n)
  {
    int diff = arr1->at(i).first - arr2->at(j).first;

    if(abs(diff) <= (int) time_block){
      int o1 = arr1->at(i).second;
      int o2 = arr2->at(j).second;
      if(o1 > o2){
        int temp = o2;
        o2 = o1;
        o1 = temp;
      }

      int size_before = orders->size();
      orders->insert(make_pair(o1, o2));
      // New co-occurrence
      if(orders->size() - size_before == 1)
        count++;

      ++j;
      i++;
    }
    else if(arr1->at(i).first > arr2->at(j).first){
        j++;
    } else {
        i++;
    }
  }

  return count;
}

void Utilities::pickRandomCooccurrencesWithinTimeBlock(vector<pair<uint, int>>* arr1, vector<pair<uint, int>>* arr2, uint time_block, set<int>* orders)
{
  int m = arr1->size();
  int n = arr2->size();
  int i = 0, j = 0;

  while (i < m && j < n)
  {
    int diff = arr1->at(i).first - arr2->at(j).first;


    if(abs(diff) <= (int) time_block) {

      if( rand()%2 == 0)
        orders->insert(arr1->at(i).second);
      else
        orders->insert(arr2->at(j).second);

      ++j;
      i++;
    }
    else if(arr1->at(i).first > arr2->at(j).first){
        j++;
    } else {
        i++;
    }
  }

}

int Utilities::countIntersectionWithinTimeBlock(vector<pair<uint, int>>* arr1, vector<pair<uint, int>>* arr2, uint time_block){
  int count=0;
  int m = arr1->size();
  int n = arr2->size();
  int i = 0, j = 0;
  while (i < m && j < n){

    int diff = arr1->at(i).first - arr2->at(j).first;

    if(abs(diff) <= (int) time_block) {
      ++count;
      ++j;
      i++;
    }
    else if(arr1->at(i) > arr2->at(j)){
        j++;
    } else {
        i++;
    }
  }


  return count;
}

int Utilities::countIntersectionWithinTimeBlock(vector<uint>* arr1, vector<uint>* arr2, uint time_block, bool debug)
{
  int count=0;
  int m = arr1->size();
  int n = arr2->size();
  int i = 0, j = 0;
  while (i < m && j < n)
  {
    int diff = arr1->at(i) - arr2->at(j);

    if(debug)
        cout<<"diff: "<<  diff <<" ";

    if(abs(diff) <= (int) time_block) {
      ++count;
      ++j;
      i++;
    }
    else if(arr1->at(i) > arr2->at(j)){
        j++;
    } else {
        i++;
    }
  }

  if(debug)
    cout<<endl;

  return count;
}



//time in microseconds
double Utilities::print_time(struct timeval &start, struct timeval &end){
    double usec;

    usec = (end.tv_sec*1000000 + end.tv_usec) - (start.tv_sec*1000000 + start.tv_usec);
    return usec;
}


int Utilities::getRandomInteger(int min, int max){
    //srand((unsigned)time(NULL));
    return (rand()%(max-min+1))+min;

}

double Utilities::getX(int p, int mod){

    int x = rand()%mod;
    int x2 = rand()%1000;

    double f2 = (double)x2/1000000.0f;
    double f  = (double)x/1000.0f;

    return p+f2+f;
}



void Utilities::addToSortedList(int *list,int n,int x){

    int i;
    for(i = n-2; i > 1 && list[i] > x; i--)
        list[i+1] = list[i];

    list[i+1] = x;
}

void Utilities::sortResPoint_AscendingId(vector<res_point*>* vc){

    struct res_point_ascending_id tmp;
    sort(vc->begin(), vc->end(), tmp);

}

void Utilities::sortResPoint_AscendingDist(vector<res_point*>* vc){

    struct res_point_ascending_dist tmp;
    sort(vc->begin(), vc->end(), tmp);
}

// erase from the vector the res_points whose dist is greater than thres
void Utilities::updateUsersInValidRange(vector<res_point*>* res, double thres){
    unsigned int i = 0;
    while(i < res->size()){
        if((*res)[i]->dist > thres)
            res->erase(res->begin()+i);
        else
            i++;
    }
}

/*
Group* Utilities::computeMyGroup(vector<res_point*>* usersInRange, int* friends, int friendsSize, res_point* p, int m){

    Group* result = new Group(p);
    priority_queue<res_point*, vector<res_point*>, res_point_ascending_dist> userFriends;

    int users = usersInRange->size();

    if(friendsSize > users){
        for(int i = 0; i< users; i++){
            if(binarySearch(friends, 0, friendsSize, (*usersInRange)[i]->id))
                userFriends.push(copy((*usersInRange)[i]));
        }
    }
    else{
        for(int i = 0; i< friendsSize; i++){
            res_point* tmp = binarySearch_res_point(usersInRange, 0, users, friends[i]);
            if(tmp != NULL)
                userFriends.push(copy(tmp));
        }
    }

    int j = 0;
    double adist = p->dist;
    while(!userFriends.empty()){
        res_point*tmp = userFriends.top();
        if(j < m){
            adist+=tmp->dist;
            j++;
        }
        result->addFriend(tmp);
        userFriends.pop()
    }
    result->adist = adist;

    return result;

}
*/

bool Utilities::binarySearch(int sortedArray[], int first, int last, int key) {
    // function:
    //   Searches sortedArray[first]..sortedArray[last] for key.
    // returns: index of the matching element if it finds key,
    //         otherwise  -(index where it could be inserted)-1.
    // parameters:
    //   sortedArray in  array of sorted (ascending) values.
    //   first, last in  lower and upper subscript bounds
    //   key         in  value to search for.
    // returns:
    //   index of key, or -insertion_position -1 if key is not
    //                 in the array. This value can easily be
    //                 transformed into the position to insert it.

    while (first <= last) {
        int mid = (first + last) / 2;  // compute mid point.
        if (key > sortedArray[mid])
            first = mid + 1;  // repeat search in top half.
        else if (key < sortedArray[mid])
            last = mid - 1; // repeat search in bottom half.
        else
            return true;     // found it. return position /////
    }
    return false;    // failed to find key
}

bool Utilities::binarySearch(vector<int>* sortedArray, int first, int last, int key) {
    // function:
    //   Searches sortedArray[first]..sortedArray[last] for key.
    // returns: index of the matching element if it finds key,
    //         otherwise  -(index where it could be inserted)-1.
    // parameters:
    //   sortedArray in  array of sorted (ascending) values.
    //   first, last in  lower and upper subscript bounds
    //   key         in  value to search for.
    // returns:
    //   index of key, or -insertion_position -1 if key is not
    //                 in the array. This value can easily be
    //                 transformed into the position to insert it.

    while (first <= last) {
        int mid = (first + last) / 2;  // compute mid point.
        if (key > (*sortedArray)[mid])
            first = mid + 1;  // repeat search in top half.
        else if (key < (*sortedArray)[mid])
            last = mid - 1; // repeat search in bottom half.
        else
            return true;     // found it. return position /////
    }
    return false;    // failed to find key
}


res_point* Utilities::binarySearch_res_point(vector<res_point*>* sortedArray, int first, int last, int key) {
    // function:
    //   Searches sortedArray[first]..sortedArray[last] for key.
    // returns: index of the matching element if it finds key,
    //         otherwise  -(index where it could be inserted)-1.
    // parameters:
    //   sortedArray in  array of sorted (ascending) values.
    //   first, last in  lower and upper subscript bounds
    //   key         in  value to search for.
    // returns:
    //   index of key, or -insertion_position -1 if key is not
    //                 in the array. This value can easily be
    //                 transformed into the position to insert it.
    while (first <= last) {
        int mid = (first + last) / 2;  // compute mid point.
        if (key > (*sortedArray)[mid]->id)
            first = mid + 1;  // repeat search in top half.
        else if (key < (*sortedArray)[mid]->id)
            last = mid - 1; // repeat search in bottom half.
        else{
            return (*sortedArray)[mid];     // found it. return position /////
        }
    }

    return NULL;    // failed to find key
}



double Utilities::computeMinimumDistance(double x1, double y1, double x2, double y2){
    return sqrt(((x1-x2)*(x1-x2))+ ((y1 - y2)*(y1 - y2)));
}

//double Utilities::computeMinimumDistance(double x, double y, double x2, double y2){

//        // pi/180 = 0.0174532925199433 (precise to double precision)

//        double dLong=(y-y2)*0.0174532925199433;
//        double dLat=(x-x2)*0.0174532925199433;

//        double aHarv= (sin(dLat/2.0)*sin(dLat/2.0))+(cos(x*0.01745329251994333)*cos(x2*0.01745329251994333)*sin(dLong/2)*sin(dLong/2));
//        double cHarv=2*atan2(sqrt(aHarv),sqrt(1.0-aHarv));

//        return 6378.137*cHarv;
//}


res_point* Utilities::copy(res_point* toBeCopied){

    res_point* rp = new res_point();
    rp->id = toBeCopied->id;
    rp->x = toBeCopied->x;
    rp->y = toBeCopied->y;
    rp->oid = toBeCopied->oid;
    rp->uid = toBeCopied->uid;
    rp->dist = toBeCopied->dist;

    return rp;
}

res_point* Utilities::createResultPoint(int id, double x, double y, double distance){

    res_point* rp = new res_point();
    rp->id = id;
    rp->x = x;
    rp->y = y;
    rp->dist = distance;

    return rp;
}





string Utilities::int_to_string(int i) {
	stringstream out;
	out << i;
	return out.str();
}
