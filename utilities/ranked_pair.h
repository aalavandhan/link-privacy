class ranked_pair{

private:
	int id1, id2;
	double score;

public:
	ranked_pair(int id1, int id2, double score): id1(id1), id2(id2), score(score){}

	int getId1() const {return id1;}
	int getId2() const {return id2;}
	double getScore() const {return score;}

};

struct ranked_pair_comparator_descending{

	bool operator()(const ranked_pair  &x, const ranked_pair  &y){
		return y.getScore() < x.getScore();
	}

	bool operator()(const ranked_pair  *x, const ranked_pair  *y){
		return y->getScore() < x->getScore();
	}

};
