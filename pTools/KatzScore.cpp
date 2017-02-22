#include "../pTools/KatzScore.h"

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <vector>
#include <queue>
#include <iostream>
#include <unordered_set>
#include <map>
#include <math.h>

using namespace std;

map< int, unordered_set<int>* > GRAPH;

KatzScore::KatzScore(){
}
KatzScore::~KatzScore(){
}

inline void print_path(vector<int>path)
{
    cout<<"[ ";
    for(int i=0;i<path.size();++i)
    {
        cout<<path[i]<<" ";
    }
    cout<<"]"<<endl;
}
bool isadjacency_node_not_present_in_current_path(int node,vector<int>path)
{
    for(int i=0;i<path.size();++i)
    {
        if(path[i]==node)
        return false;
    }
    return true;
}

int findpaths(int source ,int target, map <int, int >* pathLen_to_freq, map<int, unordered_set<int>* >* social_graph_map, int katz_path_limit)
{
    vector<int>path;
    path.push_back(source);
    queue<vector<int> >q;
    q.push(path);

    while(!q.empty())
    {
        path=q.front();
        q.pop();

        int last_nodeof_path=path[path.size()-1];
        if(last_nodeof_path==target)
        {
            auto it = pathLen_to_freq->find(path.size());
            if(it!=pathLen_to_freq->end()){
                it->second = it->second + 1;
            }else{
                pathLen_to_freq->insert(make_pair(path.size(),1));
            }

            // if(source == 0 && target == 4){
                // cout<<"The Required path is:: ";
                // print_path(path);
            // }
        }
        else
        {
            // if(source == 0 && target == 4){
                print_path(path);
            // }
        }

        if(path.size() > katz_path_limit)
            continue;

        auto it = social_graph_map->find(last_nodeof_path);
        if(it!= social_graph_map->end()){
            unordered_set<int>* friend_set = it->second;
            // if(source == 0 && target == 4){
                // cout<<"Adding number of paths to heap: "<<friend_set->size()<<endl;
            // }
            for(auto set_it = friend_set->begin(); set_it != friend_set->end(); set_it++)
            {
                int friend_id = *set_it;
                if(isadjacency_node_not_present_in_current_path(friend_id,path))
                {
                    vector<int>new_path(path.begin(),path.end());
                    new_path.push_back(friend_id);
                    q.push(new_path);

                    // if(source == 0 && target == 4){
                        // cout << "\tAdding path to heap: ";
                        // print_path(new_path);
                    // }
                }
            }

        }
    }
    return 1;
}

double KatzScore::calculateKatzScore(int source ,int target, map< int, unordered_set<int>* >* social_graph_map, double katz_epsilon, int katz_path_limit){
    map <int, int > pathLen_to_freq;
    findpaths(source, target, & pathLen_to_freq, social_graph_map, katz_path_limit);
    double katz_score = 0;
    for(auto it = pathLen_to_freq.begin(); it!= pathLen_to_freq.end(); ++it){
        int path_length = it->first;
        int path_freq = it->second;
        katz_score += pow(katz_epsilon, path_length - 1) * path_freq;
    }
    return katz_score;
}

// int main()
// {
//     // freopen("out.txt","w",stdout);
//     int T,N,M,u,v,source,target;
//     scanf("%d",&T);
//     while(T--)
//     {
//         printf("Enter Total Nodes & Total Edges\n");
//         scanf("%d%d",&N,&M);
//         for(int i=1;i<=M;++i)
//         {
//             scanf("%d%d",&u,&v);
//             auto it = GRAPH.find(u);
//             if(it == GRAPH.end()){
//                 unordered_set<int>* friend_set = new unordered_set<int>();
//                 friend_set->insert(v);
//                 GRAPH.insert(make_pair(u,friend_set));
//             }else{
//                 unordered_set<int>* friend_set = it->second;
//                 friend_set->insert(v);
//                 GRAPH.insert(make_pair(u,friend_set));
//             }
//         }
//         printf("(Source, target)\n");
//         scanf("%d%d",&source,&target);
//         cout<<KatzScore::calculateKatzScore(source,target, &GRAPH, 0.01, 3)<<endl;;
//     }
//     //system("pause");
//     return 0;
// }
