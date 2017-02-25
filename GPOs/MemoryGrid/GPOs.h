class GPOs : public IGPOs
{

private:

    // Point* locations2[DATASET_SIZE];
    map<int, Point*>::iterator it;

    //int* pages;
    int kNNExecutions;
    int LocationExecutions;
    int NextNNExecutions;
    int RangeExecutions;
    double totalCPUTime, totalTime;

    // nextNN without the incremental approach of NN
    vector<res_point*>* nextNNList;
    int computedNN, returnedNN, finalNextNN;
    bool flagNextNN;
    int objects;
    Utilities util;
    int pureNNexec;


public:
    Grid *grid;
    GPOs(char* gridFileName);
    GPOs();
    ~GPOs();
    vector<int> ids;
    map<int, Point*> locations;
    //user id to checkins
    map<int , vector< Point* >*> user_to_location;
    map<int , vector< Point* >*> location_to_user;

    // map<int, unordered_map<int, vector<uint>* >*> users_locations_frequency_map;
    map<int, map<int, vector<uint>* >*> locations_users_frequency_map;

    unordered_set< pair<int,int>, PairHasher > cooccured_user_pairs;
    unordered_set< pair<int,int>, PairHasher > significantly_cooccured_user_pairs;
    unordered_set< pair<int,int>, PairHasher > insignificantly_cooccured_user_pairs;

    unordered_map<int, double> location_to_H;
    map<int, map<int, vector<pair<int, int> >* >*> cooccurrence_matrix;

    void createNewGPOsbyGridSnapping(GPOs* gpos, double grid_distance_on_x_axis_in_km);
    virtual vector< Point* >* getLocations(int user_id);
    virtual res_point* getNextNN(double x, double y, int incrStep);
    virtual vector<res_point*>* getkNN(double x, double y, int k);
    virtual vector<res_point*>* getRange(double x, double y, double radius);
    virtual set<res_point*, res_point_ascending_id>* getSetRange(double x, double y, double radius);
    virtual vector<res_point*>* getRangeSortedId(double x, double y, double radius);
    virtual double estimateNearestDistance(double x, double y, int k);
    virtual void clearNextNN();
    virtual unordered_map<int, double>* getLocationEntropy();
    virtual map<int, map<int, vector<pair<int, int> >* >*>* getCooccurrenceMatrix();
    virtual unordered_set< pair<int,int>, PairHasher >* getCoOccurredUserPairs();
    virtual vector<int>* getUsersInRange(double x, double y, double r1, double r2);
    virtual vector<int>* getUsersInRange(double x, double y, double radius);
    virtual vector<int>* getUsersInRange(int source, double radius);
    virtual int getUserCooccurrences(int user_id);
    virtual map<int , vector< Point* >*>* getLocationToUser();

    unordered_map<int, double>* calculateLocationEntropy();

    virtual map<int, double>* getHiLasMap();
    virtual map<int, double>* getHiJasMap();
    virtual map<int, double>* getHlLasMap();

    void generateFrequencyCache();


    // nextNN without the incremental approach of NN
    //virtual res_point* getNextNN(double x, double y, int incrStep);

    int getPageAccesses();
    void resetPageAccesses();

    int getkNNExecutions();
    int getLocationExecutions();
    int getNextNNExecutions();
    int getRangeExecutions();

    bool loadLocations(const char* fileName);
	bool loadLocations(const char* fileName, int numOfFiles);
    double getTotalCPUTime();
    double getTotalTime();

    void updateCheckin(Point* p);

    void loadPoint(double x, double y, int lid, int uid, boost::posix_time::ptime time, int order);
    void groupLocationsByRange(GPOs* gpos, double radius, bool isOptimistic);
    void loadPurturbedLocations(GPOs* gpos, double radius);
    void verifyRange(double radius);
    void countU2UCoOccurrences(uint time_block);
    double distanceBetween(Point *a, Point *b);

    void loadPurturbedLocationsBasedOnNodeLocality(GPOs* gpos, map<int, double>* node_locality, double radius, double limit);
    void loadPurturbedLocationsBasedOnLocationEntropy(GPOs* gpos, double radius, double limit);
};
