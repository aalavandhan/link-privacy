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
    map< int, double >* node_locality;

    map<std::pair<int,int>, double> katzCache = map<std::pair<int,int>, double>();

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

	void printTriangles( int id, vector<int> friendList);

    int getAreFriendsExecutions();
    int getGetFriendsExecutions();

    double getTotalCPUTime();
    double getTotalTime();

    virtual double getKatzScore(int source, int target);

    int edges=0;

    double computeDistanceBetweenFriends(vector< Point* >* source_checkins, vector< Point* >* friend_checkins);
    double computeMeanDistanceBetweenAllFriends(GPOs* gpos);
    vector<double>* computeDistancesBetweenUserFriends(GPOs* gpos, int source, unordered_set<int>* friends);
    double computeNodeLocality(GPOs* gpos, int source);
    map< int, double >* computeNodeLocality(GPOs* gpos);
    map< int, double >* loadNodeLocalityFromFile();
    void writeNodeLocalityToFile();
    double precomputeKatzScore(int parts, int part);
};
