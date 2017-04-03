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

Point::~Point (){
    // delete p_x;
    // delete p_y;
    // delete p_id;
}


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
