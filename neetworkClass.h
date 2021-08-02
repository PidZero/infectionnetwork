#ifndef NETCLASS1_H
#define NETCLASS1_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <numeric>

class cl_network{
    public:
        cl_network(){};
        ~cl_network(){};

        std::vector < std::vector < double > > adjacency;
        std::vector < std::vector < int > > adList;
        std::vector < bool > S;
        std::vector < bool > R;
        std::vector < bool > I;
        int numberR;
        int numberI;
        int numberS;

    public:
        // als naechses Netzwerk mit poissonverteilten kanten um den gew'hlten knoten
        void initNet(int arg_size, int arg_meanDegree){
            N = arg_size;
            std::poisson_distribution<int> poissonDist(arg_meanDegree);
            adjacency.clear();
            adjacency.resize(arg_size);
            for(int ii = 0; ii < arg_size; ii++){
                adjacency.at(ii).resize(arg_size);
            }
            for(int ii = 0; ii < arg_size-1; ii++){
                std::uniform_int_distribution<int> uniDist(ii+1, arg_size-1);
                int N = poissonDist(gen);
                for(int jj = 0; jj < N; jj++){
                    int k = uniDist(gen);
                    adjacency.at(ii).at(k) = 1;
                    adjacency.at(k).at(ii) = 1;
                }
            }

            makeList();
        }

    public: 
        void initInfection(int arg_I, double arg_infectionProb){
            int k; 
            p = arg_infectionProb;
            std::uniform_int_distribution<int> uniDist(0, N-1);
            S.resize(N, 1);
            I.resize(N, 0);
            R.resize(N, 0);
            for(int ii = 0; ii < arg_I; ii++){
                k = uniDist(gen);
                if(I.at(k) == 0){
                    I.at(k) = 1;
                    S.at(k) = 0;
                }else{
                    ii--;
                }
            }
        }

    public:
        void timestep(){
            std::uniform_real_distribution<> chance(0., 1.);
            std::vector <bool> oldI;
            oldI.resize(N, 0);
            for(int ii = 0; ii < N; ii++){
                oldI.at(ii) = I.at(ii);
                if(oldI.at(ii)){
                    R.at(ii) = 1;
                    I.at(ii) = 0;
                }
            }
            for(int ii = 0; ii < N; ii++){
                if(oldI.at(ii)){
                    for(int jj = 0; jj < adList.at(ii).size(); jj++){

                        if((chance(gen) < p) and (S.at(adList.at(ii).at(jj)))){
                            I.at(adList.at(ii).at(jj)) = 1;
                            S.at(adList.at(ii).at(jj)) = 0;
                        }
                    }
                }
            }
            numberR = std::accumulate(R.begin(), R.end(), 0);
            numberI = std::accumulate(I.begin(), I.end(), 0);
            numberS = std::accumulate(S.begin(), S.end(), 0);
        }




    public:
        void printState(){
            int iS, iI, iR;
            iS = std::accumulate(S.begin(), S.end(), 0);
            iI = std::accumulate(I.begin(), I.end(), 0);
            iR = std::accumulate(R.begin(), R.end(), 0);
            std::cout<<iS<<"\t"<<iI<<"\t"<<iR<<"\t"<<N<<"\t"<<iS+iI+iR<<std::endl;
        }

    public:
        void printAdjacency(std::string arg_filename){
            std::ofstream myfile;
            myfile.open(arg_filename);
            myfile << ", ";
            for(int ii = 0; ii<adjacency.size(); ii++){
                if(ii < adjacency.size()){
                    myfile << ii << ", ";
                }
            }
            myfile << std::endl;

            for(int ii = 0; ii<adjacency.size(); ii++){
                myfile << ii << ", ";
                for(int jj= 0; jj < adjacency.size(); jj++){
                    myfile << adjacency.at(ii).at(jj);
                    if(jj < adjacency.size()-1){
                        myfile << ", ";
                    }
                }
                myfile << std::endl;
            }
            myfile.close();
        }


    private:
        std::mt19937_64 gen{std::random_device{}()};
        int N{0};
        double p{0};


    private:
        void makeList(){
            adList.resize(N);
            for(int ii = 0; ii < N; ii++){
                for(int jj = 0; jj < N; jj++){
                    if(adjacency.at(ii).at(jj)){
                        adList.at(ii).push_back(jj);
                    }
                }
            }
        }

};

#endif
