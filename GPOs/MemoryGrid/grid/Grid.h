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
    set<res_point*, res_point_ascending_id>* getSetRange(double x, double y, double radius);

    void getRectangle(int direction, int level, double x, double y, Cell& c);

    bool loadFromFile(const char* fileName);
	bool loadFromFile(const char* fileName, int numOfFiles);

    list<Cell*>* getIntersectedCellsWithRectangle(double x1, double y1, double x2, double y2);
    double estimateNearestDistance(double x, double y, int k);
    void deleteEmptyCells();

    Cell* makeCell(double x, double y);

    int num_failed=0;
};

