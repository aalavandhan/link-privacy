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
    map<int , vector< Point* >*> locationHistory;
    map<int , vector< Point* >*> location_to_user;

    unordered_map<int, double> location_to_H;
    map<int, map<int, vector<int, int>*>*> cooccurrence_matrix;

    virtual void getLocation(int id, double* result);
    virtual res_point* getNextNN(double x, double y, int incrStep);
    virtual vector<res_point*>* getkNN(double x, double y, int k);
    virtual vector<res_point*>* getRange(double x, double y, double radius);
    virtual set<res_point*, res_point_ascending_id>* getSetRange(double x, double y, double radius);
    virtual vector<res_point*>* getRangeSortedId(double x, double y, double radius);
    virtual double estimateNearestDistance(double x, double y, int k);
    virtual void clearNextNN();
    virtual unordered_map<int, double>* calculateLocationEntropy(map<int , set<int>*> location_History);

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

    void loadPoint(double x, double y, int lid, int uid, boost::posix_time::ptime time);
    void groupByRange(GPOs* gpos, double radius);
    void loadPurturbedLocations(GPOs* gpos, double radius);
    void verifyRange(double radius);
    void countU2UCoOccurrences();
};
