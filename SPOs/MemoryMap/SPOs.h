#include<unordered_set>
#include<unordered_map>

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
    map< int, map<int, double>* >* checkin_locality_map;

    unordered_map<std::pair<int,int>, double, PairHasher> katzCache = unordered_map<std::pair<int,int>, double, PairHasher>();

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

    double computeMinimumdistanceToFriend(GPOs* gpos, Point* point_source, vector< Point* >* friend_checkins);
    double computeCheckinLocality(GPOs* gpos, Point* point_source, unordered_set<int>* friends);
    vector<double>* computeDistancesToCheckinFriends(GPOs* gpos, Point* point_source, unordered_set<int>* friends);



    double computeDistanceBetweenFriends(vector< Point* >* source_checkins, vector< Point* >* friend_checkins);
    double computeMeanDistanceBetweenAllFriends(GPOs* gpos);
    vector<double>* computeDistancesBetweenUserFriends(GPOs* gpos, int source, unordered_set<int>* friends);
    double computeNodeLocality(GPOs* gpos, int source);
    map< int, double >* computeNodeLocality(GPOs* gpos);
    map< int, map<int, double>* >* computeCheckinLocalityMap(GPOs* gpos);

    map< int, double >* loadNodeLocalityFromFile();
    void writeNodeLocalityToFile();
    void writeCheckinLocalityToFile();
    void precomputeKatzScore(GPOs* gpos, int parts, int part, double dTresh);
    void loadKatzScoreFromMemory();
};
