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
    GPOs(GPOs *gpos);
    ~GPOs();
    vector<int> ids;
    vector<Point*> locations;
    double total_spatial_displacement = 0;
    double total_time_displacement    = 0;
    int purturbed_count = 0, spatial_purturbed_count=0, temporal_purturbed_count=0;

    //user id to checkins
    map<int , vector< Point* >*> user_to_location;
    map<int , vector< Point* >*> location_to_user;
    map<int , vector< Point* >*> time_to_checkins;

    bool cooccurrences_created=false;
    map< int, map<int,int>* >* location_to_user_to_cooccurrences;
    map< int, map<int, pair<int,double> >* >* user_to_order_to_location_displacment = NULL;

    // map<int, unordered_map<int, vector<uint>* >*> users_locations_frequency_map;
    map<int, map<int, vector<uint>* >*> locations_users_frequency_map;
    map<int, map<int, vector< pair<uint, int> >* >*> locations_users_frequency_map_with_order;

    unordered_set< pair<int,int>, PairHasher > cooccured_user_pairs;
    unordered_set< pair<int,int>, PairHasher > insignificantly_cooccured_user_pairs;
    unordered_set< pair<int,int>, PairHasher > significantly_cooccured_user_pairs;

    unordered_map<int, double> location_to_H;
    map<int, map<int, vector<pair<int, int> >* >*> cooccurrence_matrix;
    // map<int, map<int, vector<pair<int, int> >* >*> cooccurrence_matrix_insignificant;

    void createNewGPOsbyGridSnapping(GPOs* gpos, double grid_distance_on_x_axis_in_km);
    virtual vector< Point* >* getLocations(int user_id);
    virtual res_point* getNextNN(double x, double y, int incrStep);
    virtual vector<res_point*>* getkNN(double x, double y, int k);
    virtual vector<res_point*>* getRange(double x, double y, double radius);
    virtual vector<res_point*>* getRangeAndDelete(double x, double y, double radius);
    virtual set<res_point*, res_point_ascending_id>* getSetRange(double x, double y, double radius);
    virtual vector<res_point*>* getRangeSortedId(double x, double y, double radius);
    virtual double estimateNearestDistance(double x, double y, int k, double max_radius);
    virtual void clearNextNN();
    virtual unordered_map<int, double>* getLocationEntropy();
    virtual map<int, map<int, vector<pair<int, int> >* >*>* getCooccurrenceMatrix();
    // virtual map<int, map<int, vector<pair<int, int> >* >*>* getInsignificantCooccurrenceMatrix();
    virtual unordered_set< pair<int,int>, PairHasher >* getCoOccurredUserPairs();
    virtual vector<int>* getUsersInRange(double x, double y, double r1, double r2);

    virtual vector<int>* getUsersInRange(double x, double y, double radius);
    virtual vector<int>* getUsersInRange(int source, double radius);

    virtual unordered_map< int, vector<int>* >* getUsersInRangeByHourBlock(double x, double y, double r1, double r2);
    virtual unordered_map< int, vector<int>* >* getUsersInRangeByHourBlock(double x, double y, double radius);

    virtual map<int, res_point*>* getPointsInRange(double x, double y, double radius);
    virtual map<int, res_point*>* getPointsInRange(double x, double y, double r1, double r2);
    virtual int getUserCooccurrences(int user_id);
    virtual map<int , vector< Point* >*>* getLocationToUser();

    unordered_map<int, double>* calculateLocationEntropy();

    void printCooccurrenceMatrix(char *DATA_SET);

    virtual vector<Point*>* getLocations();
    virtual unordered_map<int, double>* getHiLasMap();
    virtual unordered_map<int, double>* getHiJasMap();
    virtual unordered_map<int, double>* getHlLasMap();

    virtual unordered_map< int, map< int, map<int, double >* >* >* getPltMap(double time_block, int max_checkins, double max_radius);
    virtual vector<int>* getUsersInRangeAndTimeBlock(double x, double y, double time_block, int max_checkins, double max_radius);

    virtual map< int, map<int,int>* >* getL2U2COOCC();

    void generateFrequencyCache();
    void generateCooccurrenceCache();


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
    void loadLocationsByDayOfWeek(GPOs* gpos, int day_of_week);
    double getTotalCPUTime();
    double getTotalTime();

    void loadPoint(double x, double y, int lid, int uid, boost::posix_time::ptime time, int order);
    void groupLocationsByRange(GPOs* gpos, double radius, bool isOptimistic);
    void groupLocationsByKNNDistance(GPOs* gpos, int k, double std_radio);
    void loadPurturbedLocations(GPOs* gpos, double radius);
    void loadPurturbedLocationsByTime(GPOs* gpos, uint time_deviation);
    void verifyRange(double radius);
    void countU2UCoOccurrences(uint time_block);
    double distanceBetween(Point *a, Point *b);

    void loadPurturbedLocationsBasedOnNodeLocality(GPOs* gpos, map<int, double>* node_locality, double radius, double limit);
    void loadPurturbedLocationsBasedOnLocationEntropy(GPOs* gpos, double radius, double limit);
    void loadPurturbedLocationsBasedOnCombinationFunction(GPOs* gpos,
        map< int, map<int, pair<int,double> >* >* user_to_order_to_location_locality,
        map< int, double >* order_to_temporal_locality,
        map< int, map<int,int>* >* _location_to_user_to_cooccurrences,
        double radius,
        uint time_deviation,
        bool add_spatial,
        bool add_temporal,
        int noise_function);
    // void loadPurturbedLocationsBasedOnCombinationFunctionofCOOCC(GPOs* gpos, map< int, map<int, int>* >* _location_to_user_to_cooccurrences , double radius, bool isGaussainNoise, int function_type);

    virtual map< int, double >* computeTemporalLocality(int max_checkins, double max_radius);

    Point* getKNN(Point *p, int k);
    double getKNNDistance(Point *p, int k);
    void loadPurturbedLocationKNNDistance(GPOs* gpos, bool only_cooccurrences, int k, double std_radio, map< int, map<int,int>* >* _location_to_user_to_cooccurrences);
    void loadPurturbedLocationSelectiveKNNDistance(GPOs* gpos, int k, double std_radio, uint time_range_in_seconds, map< int, map<int,int>* >* _location_to_user_to_cooccurrences);

    void computeKNNDistances(int k, bool only_cooccurrences, bool compute_spatial, bool compute_temporal, map< int, map<int,int>* >* _location_to_user_to_cooccurrences);

    double getTemporalKNNDistance(Point *p, int k);
    void getTemporalKNN(Point *p, int k, map< int, Point* > *temporalKNNs);

    void getRangeSpatioTemporalBound(Point *p, vector< Point* >* results);

    void getSkylinePoints(Point *p,
            multiset<Point, point_checkin_time_comparator_ascending>::iterator lb_it,
            multiset<Point, point_checkin_time_comparator_ascending>::iterator ub_it,
            unordered_set< Point* > *skylines);

    void computeSkylineMetrics(bool only_cooccurrences, map< int, map<int,int>* >* _location_to_user_to_cooccurrences);
    void getSpatioTemporalKNN(Point *p, int k, map< double, Point* > *spatioTemporalKNNs);
    double getSTKNNDistance(Point *p, int k);
    void computeSTKNNDistances(int k, bool only_cooccurrences, map< int, map<int,int>* >* _location_to_user_to_cooccurrences);
};
