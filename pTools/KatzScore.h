#include <unordered_set>
#include <unordered_map>
#include <map>
using namespace std;

class KatzScore{

public:
    KatzScore();
    ~KatzScore();


    static double calculateKatzScore(int source ,int target, map<int, unordered_set<int>* >* social_graph_map, double katz_epsilon, int katz_path_limit);

};

