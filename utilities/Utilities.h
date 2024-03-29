#include <unordered_set>
#include <unordered_map>

class Utilities{

public:
    Utilities();
    ~Utilities();

    double print_time(struct timeval &start, struct timeval &end);
    double getX(int p, int mod);

    bool binarySearch(int sortedArray[], int first, int last, int key);
    bool binarySearch(vector<int>* sortedArray, int first, int last, int key);
    res_point* binarySearch_res_point(vector<res_point*>* sortedArray, int first, int last, int key);
    void sortResPoint_AscendingId(vector<res_point*>* vc);
    void sortResPoint_AscendingDist(vector<res_point*>* vc);

    double computeMinimumDistance(double x1, double y1, double x2, double y2);

    void updateUsersInValidRange(vector<res_point*>* res, double thres);

    //		Group* computeMyGroup(vector<res_point*>* usersInRange, int* friends, int friendsSize, res_point* p, int m);

    res_point* copy(res_point* toBeCopied);
    res_point* createResultPoint(int id, double x, double y, double distance);
    int getRandomInteger(int min, int max);
	string int_to_string(int);

    void addToSortedList(int *list, int listSize, int newElement);
	double getTextDistanceBetween(unordered_set<string>* user_profile_1, unordered_set<string>* user_profile_2);
	double getTextDistanceBetween( unordered_map<string, double>* user_profile_1, unordered_set<string>* user_profile_2 );
	double getSocialDistanceBetween(unordered_set<int>* _f1, unordered_set<int>* _f2);
	double getDistanceBetween(int arr1[], int arr2[], int m, int n, string f);
	int countIntersection(int arr1[], int arr2[], int m, int n);
	int countUnion(int arr1[], int arr2[], int m, int n);
	double computeJaccard(int arr1[], int arr2[], int m, int n);
	double computeSetIntersection( unordered_set<int>* _f1, unordered_set<int>* _f2 );

    int countIntersectionWithinTimeBlock(vector<uint>* arr1, vector<uint>* arr2, uint time_block, bool debug);
    int countIntersectionWithinTimeBlock(vector<pair<uint, int>>* arr1, vector<pair<uint, int>>* arr2, uint time_block);
    void pickRandomCooccurrencesWithinTimeBlock(vector<pair<uint, int>>* arr1, vector<pair<uint, int>>* arr2, uint time_block, set<int>* orders);

    int getCooccurrencesWithinTimeBlock(vector<pair<uint, int>>* arr1, vector<pair<uint, int>>* arr2, uint time_block, set< pair<int,int> >* orders);

    pair<double,double> addGaussianNoise(double x, double y, double radius);
    pair<double,double> addGaussianNoise(double x, double y, double radius, double offset);
    pair<double,double> addNoise(double x, double y, double radius);

    double calulateTauB(vector<int> *A, vector<int> *B);
    double distanceBetween(double lat1, double lon1, double lat2, double lon2);
    double angleFromCoordinate(double lat1, double long1, double lat2, double long2);
    boost::posix_time::ptime addTemporalGaussianNoise(boost::posix_time::ptime time, uint deviation_in_seconds);
    boost::posix_time::ptime addTemporalGaussianNoise(boost::posix_time::ptime time, uint deviation_in_seconds, int offset);
    boost::posix_time::ptime addTemporalNoise(boost::posix_time::ptime time, int deviation_in_seconds);

    void geometricMedian(vector<Point*> *points);
    int computeNodesToAddForKAnon(int v, int e, int k);
};
