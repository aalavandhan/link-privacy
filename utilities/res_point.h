struct res_point{
    int id;
    int uid;
    int oid;
    double x;
    double y;
    double dist;
    int cached_time;
    boost::posix_time::ptime time;

    int getTimeInSeconds() const{
        boost::posix_time::ptime time_t_epoch(boost::gregorian::date(2000 ,1,1));
        boost::posix_time::time_duration time_difference = time - time_t_epoch;
        return (int)( time_difference.total_seconds() );
    }
};

struct res_point_equal : public binary_function<res_point*, res_point*, bool>
{
    bool operator()(const res_point* __x, const res_point* __y) const { return __x->id == __y->id; }
};

struct res_point_ascending_dist : public binary_function<res_point*, res_point*, bool>
{
    bool operator()(const res_point* __x, const res_point* __y) const { return __x->dist < __y->dist; }
};

struct res_point_descending_dist : public binary_function<res_point*, res_point*, bool>
{
    bool operator()(const res_point* __x, const res_point* __y) const { return __x->dist > __y->dist; }
};

struct res_point_ascending_id : public binary_function<res_point*, res_point*, bool>
{
    bool operator()(const res_point* __x, const res_point* __y) const { return __x->id < __y->id; }
};

struct res_point_descending_id : public binary_function<res_point*, res_point*, bool>
{
    bool operator()(const res_point* __x, const res_point* __y) const { return __x->id > __y->id; }
};

struct res_point_checkin_time_comparator_ascending : public binary_function<res_point*, res_point*, bool>
{
    bool operator()(const res_point* __x, const res_point* __y) const { return __x->time < __y->time; }
};
