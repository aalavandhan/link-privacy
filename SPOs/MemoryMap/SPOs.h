#include <unordered_set>
#include <unordered_map>

class SPOs : public ISPOs
{
private:
    double totalCPUTime;
    double totalTime;
    map<int, Value*> hashTable;

	map< int, unordered_set<int>* >* socialgraph_map = new map< int, unordered_set<int>* >();

//    Value* hashTable2[DATASET_SIZE];
    int areFriendsExecutions, getFriendsExecutions;
    Utilities util;
    GPOs *gpos;
    bool isGposSet;
    map< int, double >* node_locality_map;
    map< int, double >* temporal_locality_map;
    map< int, map<int, pair<int,double> >* >* checkin_locality_map;

    unordered_map<std::pair<int,int>, double, PairHasher> katzCache;

public:
    SPOs();
    SPOs(GPOs *gpos);
    ~SPOs() {}
	//map<int, int>* degreeMap;
	multiset<my_pair, pair_comparator_descending>* degreeSet;
    vector<int> ids;

	//cumulative  degree historgram
	map<int, double>* CH_DEG;

    int load(const char* file);
//    int loadfriends(const char* file);
//    void addFriendship(int user1, int user2);

    // to be implemented by ISPOs
    virtual void getFriends(int id, int*& friends,unsigned int& numOfFriends);
	virtual unordered_set<int>* getFriends(int id);
    virtual bool areFriends(int user1, int user2);
	virtual int getUserDegree(int id);
    virtual int getNumberOfFriends();
	virtual multiset<my_pair, pair_comparator_descending>* getDegreeSet();
    virtual map< int, double >* getNodeLocality();

	void printTriangles( int id, vector<int> friendList);

    int getAreFriendsExecutions();
    int getGetFriendsExecutions();

    double getTotalCPUTime();
    double getTotalTime();

    virtual double getKatzScore(int source, int target);

    int edges=0;

    //checkin locality functions
    map< int, map<int, pair<int,double>>* >* computeCheckinLocalityMap(GPOs* gpos, double max_radius);
    pair<int,double> computeMinimumdistanceToFriend(GPOs* gpos, Point* point_source, vector< Point* >* friend_checkins, double max_radius);
    double computeCheckinLocality(GPOs* gpos, Point* point_source, unordered_set<int>* friends, double max_radius);
    vector<pair<int,double>>* computeDistancesToCheckinFriends(GPOs* gpos, Point* point_source, unordered_set<int>* friends, double max_radius);
    //----------------------------


    double computeDistanceBetweenFriends(vector< Point* >* source_checkins, vector< Point* >* friend_checkins);
    double computeMeanDistanceBetweenAllFriends(GPOs* gpos);
    vector<double>* computeDistancesBetweenUserFriends(GPOs* gpos, int source, unordered_set<int>* friends);
    double computeNodeLocality(GPOs* gpos, int source);
    map< int, double >* computeNodeLocality(GPOs* gpos);
    map< int, double >* computeTemporalLocality(int max_checkins, double max_radius, GPOs *gpos);
    map< int, map<int, pair<int,double> >* >* getCheckinLocalityMap();
    map< int, double >* loadNodeLocalityFromFile();
    map< int, map<int, pair<int,double> >* >* loadCheckinLocalityFromFile(char* DATASET_PATH);
    map< int, double >* loadTemporalLocalityFromFile(char* DATASET_PATH);

    void writeNodeLocalityToFile();
    void writeCheckinLocalityToFile(char* DATASET_PATH);
    void writeTemporalLocalityToFile(char* DATASET_PATH);

    void precomputeKatzScore(GPOs* gpos, int parts, int part, double dTresh);
    void loadKatzScoreFromMemory();

    void writeUserFriendsToFile(char* DATASET_PATH);
};
