#include "../../../headersMemory.h"

Grid::Grid (){
    locations       = new map<int , Point*>();

    cout << "MIN_X = " << MIN_X << " MAX_X = " << MAX_X << " DELTA_X = " << DELTA_X << endl;
    cout << "MIN_Y = " << MIN_Y << " MAX_Y = " << MAX_Y << " DELTA_Y = " << DELTA_Y << endl;

    // memset(table, 0, sizeof table);

    // for (int i=0; i<X; i++){
    //     for (int j=0; j<Y; j++){
    //         // Cell* c = new Cell;
    //         // c->setDimensions(MIN_X+(i*DELTA_X), MIN_Y+(j*DELTA_Y), MIN_X+((i+1)*DELTA_X), MIN_Y+((j+1)*DELTA_Y));
    //         // //            cout << MIN_X+(i*DELTA_X) << " , " <<  MIN_Y+(j*DELTA_Y) << " , " <<  MIN_X+((i+1)*DELTA_X) << " , " <<  MIN_Y+((j+1)*DELTA_Y) << endl;
    //         // c->setType(CELL);
    //         // c->setIndex(i, j);
    //         // table[i][j] = c;
    //         table[i][j] = NULL;
    //     }
    // }

    cout << "GRID SETUP " << endl;

}


Grid::~Grid (){
    for (int i=0; i<X; i++){
        for (int j=0; j<Y; j++){
            delete table[i][j];
        }
    }
    delete locations;
}

void Grid::deleteEmptyCells(){
  int emptyCount=0;
  for(int i=0; i<X;i++){
    for(int j=0; j<Y;j++){
      if(table[i][j] != NULL && table[i][j]->getCheckinCount() == 0){
        emptyCount++;
        delete table[i][j];
      }
    }
  }
  cout << "GRID EMPTY : " << emptyCount << endl;
  cout << "Percentage EMPTY : " << (double) emptyCount / (X*Y) << endl;
}

Cell* Grid::getCell(double x, double y){
    int q_x = (int)((x - MIN_X)/DELTA_X);
    int q_y = (int)((y - MIN_Y)/DELTA_Y);

    // cout << "cell[" << q_x << ", " << q_y << "]" << endl;
    // cout << x << " " << MIN_X << " " << y << " " << MIN_Y << endl;

    if(q_x >=0 && q_x < X && q_y >=0 && q_y < Y)
        return table[q_x][q_y];
    else
        return NULL;
}

Cell* Grid::getCell(int i, int j){
    if(i >= 0 && j >= 0 && i < X && j < Y)
        return table[i][j];
    else
        return NULL;

}

Cell* Grid::makeCell(double x, double y){
    int q_x = (int)((x - MIN_X)/DELTA_X);
    int q_y = (int)((y - MIN_Y)/DELTA_Y);

    if(q_x >=0 && q_x < X && q_y >=0 && q_y < Y){
        if (table[q_x][q_y] != NULL){
            return table[q_x][q_y];
        }

        Cell* c = new Cell;
        c->setDimensions(MIN_X+(q_x*DELTA_X), MIN_Y+(q_y*DELTA_Y), MIN_X+((q_x+1)*DELTA_X), MIN_Y+((q_y+1)*DELTA_Y));
        //            cout << MIN_X+(q_x*DELTA_X) << " , " <<  MIN_Y+(q_y*DELTA_Y) << " , " <<  MIN_X+((q_x+1)*DELTA_X) << " , " <<  MIN_Y+((q_y+1)*DELTA_Y) << endl;
        c->setType(CELL);
        c->setIndex(q_x, q_y);
        table[q_x][q_y] = c;

        return c;
    }
    else{
        return NULL;
    }
}


bool Grid::addCheckIn(Point* user){
    Cell * c = makeCell(user->getX(),user->getY());

    if(c != NULL){
        c->newCheckIn(user);
        int id = user->getID();
        //cout << "Add checkin id = " << id << endl;
        locations->insert(pair<int, Point*>(id, user));
        //cout << "locations size = " << locations->size() << endl;
        return true;
    }
    else{
        // cout << "checkin failed! (" << user->getX() <<", " << user->getY() << endl;
        num_failed++;
        return false;
    }
}


void Grid::updateCheckIn(Point* p, double old_x, double old_y){

    //go to the previous cell and delete him
    //Point *oldP = getPoint(p->getID());
    //oldP->printDetails();
    Cell * c = getCell(old_x, old_y);
    if(c != NULL){
        c->deleteCheckIn(p);
    }
    //cout << "go to the previous cell and delete him finished " << endl;
    //the update of the hashtable has already be done in the gpos! this is only responsible for the grid
    //add him to the list of the new cell
    c = getCell(p->getX(), p->getY());
    if(c != NULL){
        c->newCheckIn(p);
    }
    //cout << "add him to the list of the new cell finished " << endl;
}


Point* Grid::getPoint(int id){

    map<int, Point*>::iterator it;
    it = locations->find(id);
    if(it != locations->end())
        return (*it).second;

    return NULL;
}

//returns the rectangle for a direction and level
// query cell (i, j)
void Grid::getRectangle(int direction, int level, double x, double y, Cell& c){

    int i = (int)((x - MIN_X)/DELTA_X);
    int j = (int)((y - MIN_Y)/DELTA_Y);

    double x1=0, y1=0, x2=0, y2=0;

    if (direction == UP){

        x1 = (i - level)*DELTA_X;
        y1 = (j + 1 + level)*DELTA_Y;
        x2 = (i + 2 + level)*DELTA_X;
        y2 = (j + 2 + level)*DELTA_Y;

    }
    else if (direction == DOWN){

        x1 = (i - level-1)*DELTA_X;
        y1 = (j - 1 - level)*DELTA_Y;
        x2 = (i + 1 + level)*DELTA_X;
        y2 = (j - level)*DELTA_Y;

    }
    else if (direction == RIGHT){

        x1 = (i + level + 1)*DELTA_X;
        y1 = (j - level - 1)*DELTA_Y;
        x2 = (i + 2 + level)*DELTA_X;
        y2 = (j + level + 1)*DELTA_Y;

    }
    else if (direction == LEFT){

        x1 = (i - 1 - level)*DELTA_X;
        y1 = (j - level)*DELTA_Y;
        x2 = (i - level)*DELTA_X;
        y2 = (j + level + 2)*DELTA_Y;

    }
    else
        cout << "Error, there is no such a direction" << endl;


    c.setDimensions(x1+MIN_X, y1+MIN_Y, x2+MIN_X, y2+MIN_Y);
    c.setDirection(direction);
    c.computeMinDist(x, y);
    c.setType(RECTANGLE);
}



list<Cell*>* Grid::getIntersectedCellsWithRectangle(double x1, double y1, double x2, double y2){

    list <Cell*> *result = new list<Cell*>;

    if(x1 < MIN_X)
        x1 = MIN_X;
    else if(x1 > MAX_X)
        x1 = MAX_X;

    if(x2 < MIN_X)
        x2 = MIN_X;
    else if(x2 > MAX_X)
        x2 = MAX_X;

    if(y1 < MIN_Y)
        y1 = MIN_Y;
    else if(y1 > MAX_Y)
        y1 = MAX_Y;

    if(y2 < MIN_Y)
        y2 = MIN_Y;
    else if(y2 > MAX_Y)
        y2 = MAX_Y;

    //cout << "x1 = " << x1 << "\t y1 = " << y1 << "\t x2 = " << x2 << "\t y2 = " << y2 << endl;

    int x_start = (int)((x1 - MIN_X)/DELTA_X);
    int x_end = (int)((x2 - MIN_X)/DELTA_X);
    int y_start = (int)((y1 - MIN_Y)/DELTA_Y);
    int y_end = (int)((y2 - MIN_Y)/DELTA_Y);

    for(int x = x_start; x <= x_end; x++){
        for(int y = y_start; y <= y_end; y++){
            Cell *c = getCell(x, y);
            //cout << "intersected cell: (" << x << ", " << y << ") => " << "x1 = " << c->getX1() << "\t y1 = " << c->getY1() << "\t x2 = " << c->getX2() << "\t y2 = " << c->getY2() << endl;
            if(c != NULL)
                result->push_front(c);
        }
    }

    return result;
}


vector<res_point*>* Grid::getkNN(double x, double y, int k){
    //	k = 10;
    //	cout << "getKNN start k = " << k << endl;
    if(k > DATASET_SIZE)
        k = DATASET_SIZE;

    //	cout << "getKNN start k = " << k << endl;
    //double x = v.getX();
    //double y = v.getY();
    vector<res_point*>* result = new vector<res_point*>();
    res_point* rp;
    double best_dist = 9999999;
    //Point* inn = NULL;
    priority_queue<Cell*, vector<Cell*>, Cell::ascending>* heap = new priority_queue<Cell*, vector<Cell*>, Cell::ascending>;
    //list<Point*>* best_NN = v.getPointersToResults();
    int levels[] = {0, 0, 0, 0};

    /*
    if(v.getResult()->size() >= DATASET_SIZE){
        v.finished == true;
        return;
    }
*/
    //cout<<"query location cell: (" << x << ", " << y <<")"<< endl;
    Cell* cq = getCell(x, y);
    cq->setMinDist(0.0);
    heap->push(cq);

    Cell* c;
    for(int i=0; i<4; i++){
        c = new Cell;
        getRectangle(i, 0, x, y, *c);
        if(c != NULL){
            heap->push(c);
            levels[i] = levels[i] + 1;
        }
    }

    Cell* tmp = heap->top();
    int count = 0;

    while ( !heap->empty() &&  tmp->getMinDist() < best_dist && count < k) {

        heap->pop();

        if (tmp->getType() == CELL){                               //it is a cell
            list<Point*>* L = tmp->getCheckIns();
            list<Point*>::iterator it;
            Point* p = NULL;
            for(it=L->begin(); it != L->end(); ++it){
                p = *it;
                Cell * c = new Cell;
                c->setType(POINT);
                c->setID(p->getID());
                c->setDimensions(p->getX(), p->getY(), p->getX(), p->getY());
                c->computeMinDist(x, y);
                heap->push(c);
            }
        }
        else if(tmp->getType() == RECTANGLE){                      //it is a rectangle for a specific direction
            int direction = tmp->getDirection();
            double x1 = tmp->getX1();
            double y1 = tmp->getY1();
            double x2 = tmp->getX2();
            double y2 = tmp->getY2();

            if (direction == UP || direction == DOWN){
                for(double i = (x1+DELTA_X/2.0); i < x2; i+=DELTA_X){
                    Cell* cr = getCell(i, (y1+DELTA_Y/2.0));
                    if(cr != NULL){
                        cr->computeMinDist(x, y);
                        heap->push(cr);
                    }
                }
            }else if (direction == RIGHT || direction == LEFT){
                for(double j = (y1+DELTA_Y/2); j < y2; j+=DELTA_Y){
                    Cell* cr = getCell((x1+DELTA_X/2), j);
                    if(cr != NULL){
                        cr->computeMinDist(x, y);
                        heap->push(cr);
                    }
                }
            }else{
                cout << "Grid::getkNN::Error, there is no such direction: " << direction << endl;
            }

            Cell* c = new Cell;
            getRectangle(direction, levels[direction], x, y, *c);
            heap->push(c);

            levels[direction] = levels[direction] + 1;
            delete tmp;
        }
        else{                                                                 // it is a point
            //v.setResult(tmp->getID(), tmp->getMinDist(), tmp->getX1(), tmp->getY1());
            //v.setBestDist(best_dist);

            rp = new res_point();
            rp->id = tmp->getID();
            rp->x = tmp->getX1();
            rp->y = tmp->getY1();
            rp->dist = tmp->getMinDist();

            //			cout << "new result = " << rp->id << endl;

            result->push_back(rp);
            count++;
            delete tmp;
        }

        tmp = heap->top();

    }
    // delete the rest heap
    while(!heap->empty()) {
        Cell* c = heap->top();
        if(c->getType() != CELL)
            delete c;
        heap->pop();
    }
    delete heap;

    //	cout << "getKNN end size = " << result->size() << endl;

    return result;
}


vector<res_point*>* Grid::getRangeAndDelete(double x, double y, double radius){
    Point* p = NULL;
    Cell* c = NULL;
    vector<res_point*>* result = new vector<res_point*>();

    int x_start = (int)(((x-radius) - MIN_X)/DELTA_X);
    int x_end = (int)(((x+radius) - MIN_X)/DELTA_X);
    int y_start = (int)(((y-radius) - MIN_Y)/DELTA_Y);
    int y_end = (int)(((y+radius) - MIN_Y)/DELTA_Y);

    if(x_end >= X)
        x_end = X-1;

    if(y_end >= Y)
        y_end = Y-1;

    if(x_start < 0)
        x_start = 0;

    if(y_start < 0)
        y_start = 0;

    int visits = 0;
    int totalcells=0;

    for(int i = x_start; i <= x_end; i++){
        for(int j = y_start; j <= y_end; j++){
            totalcells++;
            c = getCell(i, j);

            if(c != NULL && c->getCheckIns() != NULL && c->intersectsWithCircle(x, y, radius)){
                visits++;
                list<Point*>* L = c->getCheckIns();
                list<Point*>::iterator it;

                it=L->begin();
                while(it != L->end()){
                    p = *it;
                    //count++;
                    if(p->computeMinDist(x, y) <= radius){
                        res_point* rp = new res_point();
                        rp->id = p->getID();
                        rp->uid = p->getUID();
                        rp->x = p->getX();
                        rp->y = p->getY();
                        rp->dist = p->getMinDist();
                        rp->time = p->getTime();
                        rp->oid = p->getOrder();
                        result->push_back(rp);
                        auto del_iter = it;
                        ++it;
                        L->erase(del_iter);
                        delete p;
                    }else{
                        ++it;
                    }
                }
            }

        }
    }
    return result;
}


// 1. get the cells which are in the circle's MBR.
// 2. for each one of these, check if it intersects the circle,
// 3. if yes, get the checkins that are in the circle.
vector<res_point*>* Grid::getRange(double x, double y, double radius){

    //    Cell* qPoint = NULL;
    //    qPoint=getCell(x,y);
    //    cout<<"query point is in cell: ("<<qPoint->getIndexI()<<" , "<<qPoint->getIndexJ()<<")"<<endl;

    Point* p = NULL;
    Cell* c = NULL;

    vector<res_point*>* result = new vector<res_point*>();

    int x_start = (int)(((x-radius) - MIN_X)/DELTA_X);
    int x_end = (int)(((x+radius) - MIN_X)/DELTA_X);
    int y_start = (int)(((y-radius) - MIN_Y)/DELTA_Y);
    int y_end = (int)(((y+radius) - MIN_Y)/DELTA_Y);

    if(x_end >= X)
        x_end = X-1;

    if(y_end >= Y)
        y_end = Y-1;

    if(x_start < 0)
        x_start = 0;

    if(y_start < 0)
        y_start = 0;

    // int count = 0;
    int visits = 0;
    int totalcells=0;

    //    cout<<"xstart ="<<x_start<<"  ,  "<<"xend ="<<x_end<<endl;
    //    cout<<"ystart ="<<y_start<<"  ,  "<<"yend ="<<y_end<<endl;

    for(int i = x_start; i <= x_end; i++){
        for(int j = y_start; j <= y_end; j++){
            totalcells++;
            c = getCell(i, j);

            if(c != NULL && c->getCheckIns() != NULL && c->intersectsWithCircle(x, y, radius)){
                // count+=c->getCheckIns()->size();
                //cout<<"cell: "<<"("<<i<<","<<j<<")"<<endl;
                visits++;
                list<Point*>* L = c->getCheckIns();
                list<Point*>::iterator it;
                for(it=L->begin(); it != L->end(); ++it){
                    p = *it;
                    //count++;
                    if(p->computeMinDist(x, y) <= radius){
                        res_point* rp = new res_point();
                        rp->id = p->getID();
                        rp->uid = p->getUID();
                        rp->x = p->getX();
                        rp->y = p->getY();
                        rp->dist = p->getMinDist();
                        rp->time = p->getTime();
                        rp->oid = p->getOrder();
                        result->push_back(rp);
                    }
                }
            }

        }
    }

    //    v.setStatistics(visits, count);
    //    cout << "total users = " << count << endl;
    //    cout << "total cells = " << totalcells << endl;
    //    cout << "result Size = " << result->size() << endl;
    //    cout << "intersected cells = " << visits << endl;


    return result;
}



set<res_point*, res_point_ascending_id>* Grid::getSetRange(double x, double y, double radius){

    //    Cell* qPoint = NULL;
    //    qPoint=getCell(x,y);
    //    cout<<"query point is in cell: ("<<qPoint->getIndexI()<<" , "<<qPoint->getIndexJ()<<")"<<endl;

    Point* p = NULL;
    Cell* c = NULL;

    set<res_point*, res_point_ascending_id>* result = new set<res_point*, res_point_ascending_id>();

    int x_start = (int)(((x-radius) - MIN_X)/DELTA_X);
    int x_end = (int)(((x+radius) - MIN_X)/DELTA_X);
    int y_start = (int)(((y-radius) - MIN_Y)/DELTA_Y);
    int y_end = (int)(((y+radius) - MIN_Y)/DELTA_Y);

    if(x_end >= X)
        x_end = X-1;

    if(y_end >= Y)
        y_end = Y-1;

    if(x_start < 0)
        x_start = 0;

    if(y_start < 0)
        y_start = 0;

    int count = 0;
    int visits = 0;
    int totalcells=0;

    //    cout<<"xstart ="<<x_start<<"  ,  "<<"xend ="<<x_end<<endl;
    //    cout<<"ystart ="<<y_start<<"  ,  "<<"yend ="<<y_end<<endl;

    for(int i = x_start; i <= x_end; i++){
        for(int j = y_start; j <= y_end; j++){
            totalcells++;
            c = getCell(i, j);
            count+=c->getCheckIns()->size();
            if(c != NULL && c->intersectsWithCircle(x, y, radius)){
                //cout<<"cell: "<<"("<<i<<","<<j<<")"<<endl;
                visits++;
                list<Point*>* L = c->getCheckIns();
                list<Point*>::iterator it;
                for(it=L->begin(); it != L->end(); ++it){
                    p = *it;
                    //count++;
                    if(p->computeMinDist(x, y) <= radius){
                        res_point* rp = new res_point();
                        rp->id = p->getID();
                        rp->x = p->getX();
                        rp->y = p->getY();
                        rp->dist = p->getMinDist();
                        result->insert(rp);
                    }
                }
            }

        }
    }

    //    v.setStatistics(visits, count);
    //    cout << "total users = " << count << endl;
    //    cout << "total cells = " << totalcells << endl;
    //    cout << "result Size = " << result->size() << endl;
    //    cout << "intersected cells = " << visits << endl;

    return result;
}

bool Grid::loadFromFile(const char* fileName){
    //read file and add to the grid

    ifstream fin(fileName);
    if (! fin)
    {
        std::cerr << "Cannot open data file " << fileName << "." << std::endl;
        return false;
    }

    cout << "Loading Grid ..." << endl;


    double uid, x, y, lid;
    int count = 0;
    Point* l;

    while (fin)
    {
        fin >> uid >> x >> y >> lid;

        // skip newlines, etc.
        if (! fin.good()) continue;


        if(x> MIN_X && x < MAX_X && y > MIN_Y && y < MAX_Y){
            l = new Point(x, y, lid);

            addCheckIn(l);
            count ++;

            if(count%100000==0)
                cout << count << endl;
        }
    }

    cout << "Done" << endl;

    return true;
}

//returns single coordinate x1 (1) , x2(2), y1(3), y2(4)
double Grid::getCellCornerCoordinates(int c_i, int c_j, int corner_id){
    if(corner_id == 1){
        return MIN_X+(c_i*DELTA_X);
    }else if(corner_id == 2){
        return MIN_X+((c_i+1)*DELTA_X);
    }else if(corner_id == 3){
        return MIN_Y+((c_j+0)*DELTA_Y);
    }else if(corner_id == 4){
        return MIN_Y+((c_j+1)*DELTA_Y);
    }
    return 0;
}


/*
TKDE 2004 Tao, papadias
An efficient cost model for optimization of nearest neighbor search in low and medium dimensional
*/
double Grid::estimateNearestDistance(double x, double y, int k, double max_radius){
    //double lx_plus, lx_minus, ly_plus, ly_minus;
    //priority_queue<double, vector<double>, greater<double>> HP;
    int En_old = 0;
    double l_old = 0;
    double Dk, Lr, L;
    int u=0;
    double min;
    double HP[4]; // pairs (position,direction): (0, lx_plus), (1, lx_minus), (2, ly_plus), (3, ly_minus)
    int direction=0;
    int HP_lastCells[4][2];
    Point* query = new Point(x, y, -1);
    bool validity[4] = {true, true, true, true};
    Cell *queryCell = getCell(x, y);

    for(int i = 0; i < 4; i++){
        HP_lastCells[i][0] = queryCell->getIndexI();
        HP_lastCells[i][1] = queryCell->getIndexJ();
    }

    double Cv = 1.7724538509;            // Cv = sprt(pi) for 2 dimensional scale;

    if(DATASET_SIZE<k)  {
        cout << "BOUNDARY_ERROR DATASET_SIZE is less than k" << endl;
        exit(1);

        HP_lastCells[0][0] = X-1;   HP_lastCells[0][1] = Y-1;
        HP_lastCells[1][0] = X-1;   HP_lastCells[1][1] = 0;
        HP_lastCells[2][0] = 0;     HP_lastCells[2][1] = 0;
        HP_lastCells[3][0] = 0;     HP_lastCells[3][1] = Y-1;

        HP[0] = query->computeMinDist(getCell(HP_lastCells[0][0], HP_lastCells[0][1])->getX2(), getCell(HP_lastCells[0][0], HP_lastCells[0][1])->getY2());
        HP[1] = query->computeMinDist(getCell(HP_lastCells[1][0], HP_lastCells[1][1])->getX2(), getCell(HP_lastCells[1][0], HP_lastCells[1][1])->getY1());
        HP[2] = query->computeMinDist(getCell(HP_lastCells[2][0], HP_lastCells[2][1])->getX1(), getCell(HP_lastCells[2][0], HP_lastCells[2][1])->getY1());
        HP[3] = query->computeMinDist(getCell(HP_lastCells[3][0], HP_lastCells[3][1])->getX1(), getCell(HP_lastCells[3][0], HP_lastCells[3][1])->getY2());

        double max = -1;
        for(int i = 0; i < 4; i++){
            if(max < HP[i]){
                max = HP[i];
            }
        }
        Dk=2*max/Cv;
        /*
        cout <<"--------------------------------------" << endl;
        cout<<"Data Dependent Estimation Completed with"<<endl;
        cout<<"Estimated radius (Dk) of circle = " << Dk <<endl;
        cout<< "Actual Users in the estimated radius = " << getRange(x, y, Dk)->size() << endl;
        cout<< "--------------------------------------" << endl;
        */
        return Dk;
    }

    do{

        HP[0] = query->computeMinDist(getCellCornerCoordinates(HP_lastCells[0][0], HP_lastCells[0][1], 2), y);
        HP[1] = query->computeMinDist(getCellCornerCoordinates(HP_lastCells[1][0], HP_lastCells[1][1], 1), y);
        HP[2] = query->computeMinDist(x, getCellCornerCoordinates(HP_lastCells[2][0], HP_lastCells[2][1], 4));
        HP[3] = query->computeMinDist(x, getCellCornerCoordinates(HP_lastCells[3][0], HP_lastCells[3][1], 3));

        min = MAXDIST;
        for(int i = 0; i < 4; i++){
            if(min > HP[i] && validity[i]){
                direction = i;
                min = HP[i];
            }
        }

        // cout << "Min("<< direction << ") = " << min << endl;

        L = 2*min;

        list<Cell*>* intersectedCells = getIntersectedCellsWithRectangle(x-L/2+BOUNDARY_ERROR, y-L/2+BOUNDARY_ERROR, x+L/2-BOUNDARY_ERROR, y+L/2-BOUNDARY_ERROR);

        int usersInIntersectedCells = 0;
        u=0;
        double intersectedVolPercentage = 0;

        for(list<Cell*>::iterator it = intersectedCells->begin(); it!= intersectedCells->end(); ++it){

            Cell *tmp = *it;

            intersectedVolPercentage = tmp->intersectedVolumeWithRectangle(x-L/2, y-L/2, x+L/2, y+L/2)/(DELTA_X*DELTA_Y);
            if(intersectedVolPercentage > 1)
                intersectedVolPercentage = 1;

            u+= ceil(tmp->getCheckIns()->size()*intersectedVolPercentage);
            usersInIntersectedCells+=tmp->getCheckIns()->size();
        }

        delete intersectedCells;

        //u=(int)usersCount;

        if(k <= u || min >= max_radius){
            // compute Lr
            Lr = pow((l_old*l_old*(k-u)-L*L*(k-En_old))/(En_old-u), 0.5);
            // cout << "Total # of intersected cells = " << intersectedCells->size()<<endl;
            // cout<< "--------------------------------------" << endl;
            // cout<< "Total number of users in the intersected cells = " << usersInIntersectedCells << endl;

            // compute Dk
            Dk = Lr/Cv;

            // cout <<"--------------------------------------" << endl;
            // cout<<"Data Dependent Estimation Completed with"<<endl;
            // cout<<"Estimated radius (Dk) of circle = " << Dk <<endl;
            // cout<< "--------------------------------------" << endl;
            // cout<<"Estimated number of users enclosed(En) in the rectangle"<<endl;
            // cout<<"with extent Lr = "<<Lr<<" are within " << En_old <<" and "<<u<<"."<< endl;
            // cout<< "--------------------------------------" << endl;
            // cout<< "Actual Users in the estimated radius = " << getRange(x, y, Dk)->size() << endl;
            // cout<< "--------------------------------------" << endl;

            return Dk;
        }
        else{
            l_old = L;
            En_old = u;
            bool flag = false;

            if(direction == 0){
                if(HP_lastCells[0][0] > 0 && HP_lastCells[0][0] < X-1){
                    HP_lastCells[0][0] = HP_lastCells[0][0] +1;
                    flag =true;
                }
            }
            else if(direction == 1){
                if(HP_lastCells[1][0] > 0 && HP_lastCells[1][0] < X-1){
                    HP_lastCells[1][0] = HP_lastCells[1][0] -1;
                    flag =true;
                }
            }
            else if(direction == 2){
                if(HP_lastCells[2][1] > 0 && HP_lastCells[2][1] < Y-1){
                    HP_lastCells[2][1] = HP_lastCells[2][1] +1;
                    flag =true;
                }
            }
            else {
                if(HP_lastCells[3][1] > 0 && HP_lastCells[3][1] < Y-1){
                    HP_lastCells[3][1] = HP_lastCells[3][1] -1;
                    flag =true;
                }
            }
            if(!flag)
                validity[direction] = false;
        }
    }while(true);
}



