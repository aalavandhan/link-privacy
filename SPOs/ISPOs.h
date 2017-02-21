class ISPOs
{
public:
    virtual void getFriends(int id, int*& friends,unsigned int& size) = 0;
    virtual bool areFriends(int user1, int user2) = 0;
	  virtual int getUserDegree(int id ) =0;
    virtual int getNumberOfFriends( ) =0;
    virtual double getKatzScore(int source, int target) = 0;
    virtual map< int, double >* getNodeLocality()=0;
};
