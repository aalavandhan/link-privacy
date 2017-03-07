class SimpleQueries{

private:
    IGPOs *gpos;
    ISPOs *spos;
    Utilities util;

    map<int, map<int, double>*> social_strength_matrix;
    map<int, map<int, double>*> diversity_matrix;
    map<int, map<int, double>*> weighted_frequency_matrix;

public:
    SimpleQueries(IGPOs *spatialIndex, ISPOs *socialGraph);
    ~SimpleQueries();

    map<int, map<int, double>*> cacluateSocialStrength();
    void verifySocialStrength(double tresh);
    void buildMatrices(double q);
    int countCooccurredFriends();
    map< int, bool >* getUsersOfInterest(double tresh);
    void checkUtilityRange(const char* fileName, IGPOs *base_gpos, double radius, double noise_distance);
    void checkUtilityProximity(const char* fileName, IGPOs *base_gpos, double radius, double tresh, double noise_distance);
    unordered_set< pair<int,int>, PairHasher >* computeProximityUserList(double x, double y, double radius, double tresh, double noise_distance);
    void checkUtilityStats(const char* fileName, double radius, double noise_radius);
    void cacluateCooccurrenceDistribution(vector <int> *users);
    void cacluateCooccurrenceDistributionBasedOnNodeLocality();
    void cacluateCooccurrenceDistributionBasedOnLocationEntropy();

    bool areEBMFriends(int source, int target, double tresh);
    bool areTrueEBMFriends(int source, int target, double tresh);
    void writeHistogramstoFile(double tresh);

    void printPartialDiversityAndWeightedFrequencyValues(double alpha);
};

