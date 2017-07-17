class SimpleQueries{

private:
    GPOs *gpos;
    ISPOs *spos;
    Utilities util;

    map<int, map<int, double>*> social_strength_matrix;
    map<int, map<int, double>*> diversity_matrix;
    map<int, map<int, double>*> weighted_frequency_matrix;

public:
    SimpleQueries(GPOs *spatialIndex, ISPOs *socialGraph);
    ~SimpleQueries();

    map<int, map<int, double>*> cacluateSocialStrength();
    void countEBMInferredFriendships(double tresh);
    void verifySocialStrength(double tresh);
    void buildMatrices(double q);
    int countCooccurredFriends();
    map< int, bool >* getUsersOfInterest(double tresh);
    void checkUtilityRange(const char* fileName, GPOs *base_gpos, double radius);
    void checkUtilityProximity(const char* fileName, GPOs *base_gpos, double radius, double tresh, double noise_distance);
    unordered_set< pair<int,int>, PairHasher >* computeProximityUserList(double x, double y, double radius, double tresh, double noise_distance);
    void checkUtilityStats(const char* fileName, double radius, double noise_radius);
    void getInterestingQueryPoints(double radius, const char* query_file, int DATA_SET);

    void cacluateCooccurrenceDistribution(vector <int> *users);
    void cacluateCooccurrenceDistributionBasedOnNodeLocality();
    void cacluateCooccurrenceDistributionBasedOnLocationEntropy();

    double getSocialStrength(int source, int target);
    bool areEBMFriends(int source, int target, double tresh);
    bool areTrueEBMFriends(int source, int target, double tresh);
    void writeHistogramstoFile(char *DATASET_PATH, double tresh, double time_block, map< int, double >* temoral_locality_map);

    void printPartialDiversityAndWeightedFrequencyValues(char *DATASET_PATH, double alpha);

    void computeAccuracyOfSocialStrength(double precision_limit);
    void generateSocialGraph(char *DATASET_PATH, double tresh);

    void checkUtilityBasic(GPOs *base_gpos);
    void checkUtilityKNN(const char* fileName, GPOs *base_gpos);
    void checkUtilityRKNN(int k);
};

