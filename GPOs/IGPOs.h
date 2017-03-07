class IGPOs
{
public:
    virtual vector<res_point*>* getkNN(double x, double y, int k) = 0;
    virtual vector<res_point*>* getRange(double x, double y, double radius) = 0;
    virtual set<res_point*, res_point_ascending_id>* getSetRange(double x, double y, double radius) = 0;
    virtual res_point* getNextNN(double x, double y, int incrStep) = 0;
    virtual void clearNextNN() = 0;
    virtual vector<res_point*>* getRangeSortedId(double x, double y, double radius) = 0;
    virtual double estimateNearestDistance(double x, double y, int k) = 0;
    virtual unordered_map<int, double>* getLocationEntropy() =0;
    virtual map<int, map<int, vector<pair<int, int> >* >*>* getCooccurrenceMatrix()=0;
    virtual map<int, map<int, vector<pair<int, int> >* >*>* getInsignificantCooccurrenceMatrix()=0;
    virtual unordered_set< pair<int,int>, PairHasher >* getCoOccurredUserPairs()=0;
    virtual vector<int>* getUsersInRange(double x, double y, double r1, double r2)=0;
    virtual vector<int>* getUsersInRange(double x, double y, double radius)=0;
    virtual vector<int>* getUsersInRange(int source, double radius)=0;
    virtual map<int, res_point*>* getPointsInRange(double x, double y, double radius)=0;
    virtual map<int, res_point*>* getPointsInRange(double x, double y, double r1, double r2)=0;
    virtual int getUserCooccurrences(int user_id)=0;
    virtual map<int , vector< Point*>* >* getLocationToUser()=0;
    virtual map<int, double>* getHiLasMap() = 0;
    virtual map<int, double>* getHiJasMap() = 0;
    virtual map<int, double>* getHlLasMap() = 0;
    virtual map< int, map<int,int>* >* getL2U2COOCC() = 0;
    virtual void printCooccurrenceMatrix() = 0;

};

