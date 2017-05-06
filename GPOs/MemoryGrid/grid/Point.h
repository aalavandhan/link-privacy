class Point {

public:
    Point ();
    Point (res_point *p);
    Point (int time);
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
    virtual double computeMinDistInKiloMeters(double x, double y);
    virtual void printDetails();
    virtual double getMinDist();
    virtual int getTimeBlock(int time_block_duration);
    virtual int getTimeInSeconds() const;
    virtual int getWeekTimeBlock();
    virtual int getCheckinHour();
    virtual int getCheckinDay();
    virtual int getTimeIndex();
    virtual int getTimeDifference(Point *q);
    virtual int getTimeDifference(res_point *q);
    virtual bool doesSkylineDominatePoint(Point *skyline, Point *other);
    virtual bool doesPointDominateSkyline(Point *skyline, Point *other);

    //Overload the < operator.
    bool operator< (const Point& p) const;
    //Overload the > operator.
    bool operator> (const Point& p) const;

    struct ascending : public std::binary_function<Point*, Point*, bool>
    {
        bool operator()(const Point* __x, const Point* __y) const { return __x->p_minDist > __y->p_minDist; }
    };

    static boost::posix_time::ptime START_DATE_TIME;

private:
    int p_id, p_uid, p_order;
    double p_minDist;
    double p_x, p_y;
    boost::posix_time::ptime p_time;
};

struct point_checkin_time_comparator_ascending : public binary_function<Point, Point, bool> {

    bool operator()(const Point  &x, const Point  &y) const {
        return y.getTimeInSeconds() > x.getTimeInSeconds();
    }

};
