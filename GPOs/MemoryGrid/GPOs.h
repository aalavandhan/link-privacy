class GPOs
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
    GPOs(char* gridFileName, uint time_range_in_seconds, double spatial_range);
    GPOs(uint time_range, double spatial_range);
    GPOs(GPOs *gpos);
    ~GPOs();
    set<int> ids;
    vector<Point*> locations;
    double total_spatial_displacement = 0;
    double total_time_displacement    = 0;
    int purturbed_count = 0, spatial_purturbed_count=0, temporal_purturbed_count=0;

    uint coocc_time_range;
    double coocc_spatial_range;

    //user id to checkins
    map<int , vector< Point* >*> user_to_location;
    map<int , vector< Point* >*> location_to_user;
    map<int , vector< Point* >*> time_to_checkins;
    map< int, Point* > checkin_list;

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
    vector< Point* >* getLocations(int user_id);
    res_point* getNextNN(double x, double y, int incrStep);
    vector<res_point*>* getkNN(double x, double y, int k);
    vector<res_point*>* getRange(double x, double y, double radius);
    vector<res_point*>* getRange(Point *original, double radius, double t_dist);
    vector<res_point*>* getRange(Point *original, double radius, double t_dist, double st_distance);
    vector<res_point*>* getRangeAndDelete(double x, double y, double radius);
    vector<res_point*>* getRangeAndDelete(Point *p, double radius, double t_dist);
    vector<res_point*>* getRangeSortedByTime(double x, double y, double radius);
    set<res_point*, res_point_ascending_id>* getSetRange(double x, double y, double radius);
    vector<res_point*>* getRangeSortedId(double x, double y, double radius);
    void getRangeByTime(int time_start, int time_end, vector<Point*>* results);
    double estimateNearestDistance(double x, double y, int k, double max_radius);
    void clearNextNN();
    unordered_map<int, double>* getLocationEntropy();
    map<int, map<int, vector<pair<int, int> >* >*>* getCooccurrenceMatrix();
    // map<int, map<int, vector<pair<int, int> >* >*>* getInsignificantCooccurrenceMatrix();
    unordered_set< pair<int,int>, PairHasher >* getCoOccurredUserPairs();
    vector<int>* getUsersInRange(double x, double y, double r1, double r2);

    vector<int>* getUsersInRange(double x, double y, double radius);
    vector<int>* getLocationsInRange(double x, double y, double radius);

    vector<int>* getUsersInRange(int source, double radius);

    unordered_map< int, vector<int>* >* getUsersInRangeByHourBlock(double x, double y, double r1, double r2);
    unordered_map< int, vector<int>* >* getUsersInRangeByHourBlock(double x, double y, double radius);

    map<int, res_point*>* getPointsInRange(double x, double y, double radius);
    map<int, res_point*>* getPointsInRange(double x, double y, double r1, double r2);
    int getUserCooccurrences(int user_id);
    map<int , vector< Point* >*>* getLocationToUser();

    unordered_map<int, double>* calculateLocationEntropy();

    void printCooccurrenceMatrix(char *DATA_SET);

    vector<Point*>* getLocations();
    unordered_map<int, double>* getHiLasMap();
    unordered_map<int, double>* getHiJasMap();
    unordered_map<int, double>* getHlLasMap();

    unordered_map< int, map< int, map<int, double >* >* >* getPltMap(double time_block, int max_checkins, double max_radius);
    vector<int>* getUsersInRangeAndTimeBlock(double x, double y, double time_block, int max_checkins, double max_radius);

    map< int, map<int,int>* >* getL2U2COOCC();

    void generateCooccurrenceCache();

    set< pair<int,int> >* getCooccurredCheckins();
    set< pair<int,int> > cooccurred_checkins;
    unordered_map< int, unordered_set<int>* > cooccurrence_index;
    set<int> unpurtrubed_cooccurrences;

    // nextNN without the incremental approach of NN
    //virtual res_point* getNextNN(double x, double y, int incrStep);

    double bound_computation_time =0 , metric_computation_time=0;
    double range_query_computation = 0, cooccurence_selection_computation = 0;
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
    void groupLocationsByDD(GPOs* gpos, int k);
    void groupLocationsByDD(GPOs* gpos, int k, double factor);

    void groupLocationsByKNNDistance(GPOs* gpos, int k, double std_radio);
    void groupLocationsByST(GPOs* gpos, double radius, double time_deviation);
    void loadPurturbedLocations(GPOs* gpos, double radius);
    void loadPurturbedLocationsByTime(GPOs* gpos, uint time_deviation);
    void verifyRange(double radius);
    void countU2UCoOccurrences();
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

    void loadPurturbedBasedOnSelectiveGaussian(GPOs* gpos, double radius, uint time_deviation);
    void loadPurturbedBasedOnSelectiveSTKNNDistance(GPOs* gpos, int k);
    void loadPurturbedBasedOnSelectiveSkyline(GPOs* gpos, int k);

    virtual map< int, double >* computeTemporalLocality(int max_checkins, double max_radius);

    Point* getKNN(Point *p, int k);
    double getKNNDistance(Point *p, int k);
    void loadPurturbedLocationKNNDistance(GPOs* gpos, bool only_cooccurrences, int k, double std_radio, map< int, map<int,int>* >* _location_to_user_to_cooccurrences);
    void loadPurturbedLocationSelectiveKNNDistance(GPOs* gpos, int k, double std_radio, map< int, map<int,int>* >* _location_to_user_to_cooccurrences);

    void pickSingleCheckinFromCooccurrences(set<int> *checkins_of_interest);
    void pickOtherCheckinFromCooccurrences(set<int> *checkins_of_interest);
    void pickUniqueCheckinFromCooccurrences(set<int> *checkins_of_interest);

    vector <res_point*>* getRangeSpatioTemporalBound(Point *p);

    void computeSkylineMetrics(map< int, map<int,int>* >* _location_to_user_to_cooccurrences);
    void getSkylinePoints(Point *p, vector <res_point*> *spatial_candidates, map< int, pair<int, res_point*> > *skylines);

    void computeSTKNNDistancesOptimized(int k, int type);
    void computeSTKNNDistances(int k, int type);
    double getSTKNNDistance(Point *p, int k, vector<res_point*> *spatial_candidates, int type);
    void getSpatioTemporalKNN(Point *p, int k,
        priority_queue < pair<double, res_point*>, vector<pair<double, res_point*> > > *spatioTemporalKNNs,
        vector<res_point*> *spatial_candidates,
        int type);

    void countCoOccurrencesOptimistic();
    void countCoOccurrencesOptimisticDD(int k);

    pair<double, double> maxDistanceOutsideCooccurrence(Point *p);
};
