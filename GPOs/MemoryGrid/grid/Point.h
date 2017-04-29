class Point {

public:
    Point ();
    Point (res_point *p);
    Point (double x, double y, int id);
    Point (double x, double y, int id, int uid, boost::posix_time::ptime time, int order);
    virtual ~Point ();

    virtual double getX();
    virtual double getY();
    virtual bool updateXY(double xx, double yy);
    virtual int getID();
    virtual int getUID();
    virtual boost::posix_time::ptime getTime();
    virtual int getOrder();
    virtual double computeMinDist(double x, double y);
    virtual double computeMinDistInMeters(double x, double y);
    virtual void printDetails();
    virtual double getMinDist();
    virtual int getTimeBlock(int time_block_duration);
    virtual int getWeekTimeBlock();
    virtual int getCheckinHour();
    virtual int getCheckinDay();

    //Overload the < operator.
    bool operator< (const Point& p) const;
    //Overload the > operator.
    bool operator> (const Point& p) const;

    struct ascending : public std::binary_function<Point*, Point*, bool>
    {
        bool operator()(const Point* __x, const Point* __y) const { return __x->p_minDist > __y->p_minDist; }
    };

private:
    int p_id, p_uid, p_order;
    double p_minDist;
    double p_x, p_y;
    boost::posix_time::ptime p_time;
};


