#include "../headers.h"

SimpleQueries::SimpleQueries(IGPOs *spatialIndex, ISPOs *socialGraph){
    spos = socialGraph;
    gpos = spatialIndex;
}


SimpleQueries::~SimpleQueries(){}


map<int, vector<my_pair>*> SimpleQueries::cacluateSocialStrength(){
    for (auto d_it = diversity_matrix.begin(); d_it != diversity_matrix.end(); d_it++){
        int user_1 = d_it->first;
        auto user_diversity_list = d_it->second;

        auto fm_fit = weighted_frequency_matrix.find(user_1);
        if(fm_fit == weighted_frequency_matrix.end())
            cout<<"ERROR---user not found--------"<<endl;
        auto user_wfreq_list = fm_fit->second;

        for(auto list_it = user_diversity_list->begin(); list_it!= user_diversity_list->end();list_it++){
            int user_2 = list_it->first;
            double diversity_value = list_it->second;
            auto f_it = user_wfreq_list->find(user_2);
            double weighted_frequency_value = f_it->second;
            double social_strength_value = (ALPHA * diversity_value) + (BETA * weighted_frequency_value) + GAMMA;

            //swap values of user_1 and user_2 to ensure user1 is less than user2
            if(user_1 > user_2){
                int temp = user_1;
                user_2 = user_1;
                user_1 = temp;
            }

            //find in matrix
            //if not found create new list
            // insert list as pair
            //if found
            //insert into existing list
            auto ss_it = social_strength_matrix.find(user_1);
            if(ss_it == social_strength_matrix.end()){
                vector<my_pair>* sslist = new vector<my_pair>();
                sslist->push_back(my_pair(user_2,social_strength_value));
                social_strength_matrix.insert(make_pair(user_1, sslist));
            }else{
                vector<my_pair>* sslist = ss_it->second;
                sslist->push_back(my_pair(user_2,social_strength_value));
            }
        }
    }
    return social_strength_matrix;
}
