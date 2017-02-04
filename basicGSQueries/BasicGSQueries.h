class SimpleQueries{

private:
    IGPOs *gpos;
    ISPOs *spos;
    Utilities util;

    map<int, vector<my_pair>*> social_strength_matrix;
    map<int, map<int, double>*> diversity_matrix;
    map<int, map<int, double>*> weighted_frequency_matrix;

public:
    SimpleQueries(IGPOs *spatialIndex, ISPOs *socialGraph);
    ~SimpleQueries();

    map<int, vector<my_pair>*> cacluateSocialStrength();
    void verifySocialStrength(double tresh);
    void buildMatrices(double q);
    int countCooccurredFriends();
    map< int, bool >* getUsersOfInterest(double tresh);
    void checkUtilityRange(const char* fileName, IGPOs *base_gpos, double radius);
};

