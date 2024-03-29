class Grid {

private:
    Cell* table[X+1][Y+1];
    int x, y;
    map<int , Point*>* locations;

public:
    Grid ();
    ~Grid ();

    Cell* getCell(double x, double y);
    Cell* getCell(int i, int j);
    bool addCheckIn(Point* user);
    void updateCheckIn(Point* p, double old_x, double old_y);
    Point* getPoint(int id);

    vector<res_point*>* getkNN(double x, double y, int k);
    vector<res_point*>* getRange(double x, double y, double radius);
    vector<res_point*>* getRange(Point *original, double radius, double t_dist);
    vector<res_point*>* getRangeAndDelete(double x, double y, double radius);
    vector<res_point*>* getRangeAndDelete(Point *original, double radius, double t_dist);
    set<res_point*, res_point_ascending_id>* getSetRange(double x, double y, double radius);
    multiset<res_point*, res_point_checkin_time_comparator_ascending>* getSetRangeByTime(double x, double y, double radius);

    void getRectangle(int direction, int level, double x, double y, Cell& c);

    bool loadFromFile(const char* fileName);
	bool loadFromFile(const char* fileName, int numOfFiles);

    list<Cell*>* getIntersectedCellsWithRectangle(double x1, double y1, double x2, double y2);
    double estimateNearestDistance(double x, double y, int k, double max_radius);
    void deleteEmptyCells();

    Cell* makeCell(double x, double y);

    int num_failed=0;

    double getCellCornerCoordinates(int c_i, int c_j, int corner_id);
};

