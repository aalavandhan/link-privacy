#include "../../../headersMemory.h"

Point::Point (){}

Point::Point (double x, double y, int id){
    p_x = x;
    p_y = y;
    p_id = id;
}

Point::Point (double x, double y, int id, int uid, boost::posix_time::ptime time, int order){
    p_x = x;
    p_y = y;
    p_id = id;
    p_uid = uid;
    p_time = time;
    p_order = order;
}

Point::Point (res_point *p){
    p_x = p->x;
    p_y = p->y;
    p_id = p->id;
    p_uid = p->uid;
    p_order = p->oid;
    p_time = p->time;
}

Point::Point (int time){
    p_x = -1;
    p_y = -1;
    p_id = -1;
    p_uid = -1;
    p_order = -1;
    boost::posix_time::ptime time_t_epoch(boost::gregorian::date(2000 ,1,1));
    p_time = time_t_epoch + boost::posix_time::seconds(time);
}

Point::~Point (){
    // delete p_x;
    // delete p_y;
    // delete p_id;
}

boost::posix_time::ptime time_t_epoch(boost::gregorian::date(2000 ,1,1));
boost::posix_time::ptime Point::START_DATE_TIME = time_t_epoch;


double Point::getX(){return p_x;}

double Point::getY(){return p_y;}

int Point::getID(){return p_id;}

int Point::getUID(){return p_uid;}

int Point::getOrder(){return p_order;}

boost::posix_time::ptime Point::getTime(){return p_time;}


// Euclidean distance
double Point::computeMinDist(double x, double y){
    p_minDist = sqrt((x - p_x)*(x - p_x) + (y - p_y)*(y - p_y));
    return p_minDist;
}

// Distance in Meters
double Point::computeMinDistInKiloMeters(double x, double y){
    double lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = DEG_TO_RAD * y;
    lon1r = DEG_TO_RAD * x;
    lat2r = DEG_TO_RAD * p_y;
    lon2r = DEG_TO_RAD * p_x;
    u = sin((lat2r - lat1r)/2);
    v = sin((lon2r - lon1r)/2);
    return 2.0 * EARTH_RADIUS_IN_KILOMETERS * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

/*
// Haversine distance
double Point::computeMinDist(double x, double y){

    // pi/180 = 0.0174532925199433 (precise to double precision)

    double dLong=(y-p_y)*0.0174532925199433;
    double dLat=(x-p_x)*0.0174532925199433;

    double aHarv = (sin(dLat/2.0)*sin(dLat/2.0))+(cos(x*0.01745329251994333)*cos(p_x*0.01745329251994333)*sin(dLong/2)*sin(dLong/2));
    double cHarv = 2*atan2(sqrt(aHarv),sqrt(1.0-aHarv));

    p_minDist = 6378.137*cHarv;

    return p_minDist;
}
*/

void Point::printDetails(){
    cout << "ID = " << p_id << "\t (" << p_x << ", " << p_y << ")" << " dist = " << p_minDist << endl;
}

double Point::getMinDist(){
    return p_minDist;
}

bool Point::updateXY(double xx, double yy){
    p_x = xx;
    p_y = yy;
    return true;
}

//Overload the < operator.
bool Point::operator< (const Point& p) const
{
    return p_minDist > p.p_minDist;
}

//Overload the > operator.
bool Point::operator> (const Point& p) const
{
    return p_minDist < p.p_minDist;
}

int Point::getTimeBlock(int time_block_duration){
    int p_time_block = (int)( (double)( p_time.time_of_day().hours() * 60 + p_time.time_of_day().minutes() ) / time_block_duration );
    return p_time_block;
}

bool Point::doesSkylineDominatePoint(res_point *skyline, res_point *other){
    int sk_del_t    = getTimeDifference(skyline);
    double sk_del_s = computeMinDistInKiloMeters(skyline->x, skyline->y);

    int ot_del_t    = getTimeDifference(other);
    double ot_del_s = computeMinDistInKiloMeters(other->x, other->y);

    // cout << "Skyline : " << skyline->getOrder() << " "<< sk_del_s << " " << sk_del_t << endl;
    // cout << "Point   : " << other->getOrder() << " " << ot_del_s << " " << ot_del_t << endl;
    // cout << "Domination : " << (sk_del_s <= ot_del_s && sk_del_t <= ot_del_t) << endl;

    return (sk_del_s < ot_del_s + BOUNDARY_ERROR && sk_del_t <= ot_del_t);
}

bool Point::doesPointDominateSkyline(res_point *skyline, res_point *other){
    int sk_del_t    = getTimeDifference(skyline);
    double sk_del_s = computeMinDistInKiloMeters(skyline->x, skyline->y);

    int ot_del_t    = getTimeDifference(other);
    double ot_del_s = computeMinDistInKiloMeters(other->x, other->y);

    return (ot_del_s < sk_del_s + BOUNDARY_ERROR && ot_del_t < sk_del_t);
}

int Point::getTimeInSeconds() const{
    boost::posix_time::time_duration time_difference = p_time - START_DATE_TIME;
    return (int)( time_difference.total_seconds() );
}

// 1 Hour time blocks
int Point::getTimeIndex(){
    return (int)( (double) getTimeInSeconds() / (double) (3 * 60 * 60) );
}

int Point::getTimeDifference(Point *q){
    return abs((p_time - q->getTime()).total_seconds());
}

int Point::getTimeDifference(res_point *q){
    return abs((p_time - q->time).total_seconds());
}

int Point::getWeekTimeBlock(){
    int p_time_block = p_time.date().day_of_week() * 24 + p_time.time_of_day().hours();
    return p_time_block;
}

int Point::getCheckinHour(){
    return p_time.time_of_day().hours();
}

int Point::getCheckinDay(){
    return p_time.date().day_of_week();
}
