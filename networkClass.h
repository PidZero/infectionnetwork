// //////////////////////////////////////////////// //
// (c) J. Neidhart, 2021                            //
// Basic Network Class for simple SIR dynamics      //
// //////////////////////////////////////////////// //
#ifndef NETCLASS1_H
#define NETCLASS1_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <numeric>
#include <cmath>

class cl_network{
    public:
        cl_network(){};
        ~cl_network(){};
        
        // Adjaceny Matrix
        std::vector < std::vector < double > > adjacency;

        // Adjacency List
        std::vector < std::vector < int > > adList;

        // Vector to store the labels of the susceptible
        std::vector < bool > S;

        // Vector to store the labels of the resistant
        std::vector < bool > R;

        // Vector to store the labels of the infected
        std::vector < bool > I;

        // Numbers of the resistant, infected and susceptible
        int numberR;
        int numberI;
        int numberS;

    public:
        // create a network with Poisson distributed number of random edges
        void initNetRandom(int arg_size, int arg_meanDegree){
            N = arg_size;
            std::poisson_distribution<int> poissonDist(arg_meanDegree);
            adjacency.clear();
            adjacency.resize(arg_size);
            for(int ii = 0; ii < arg_size; ii++){
                adjacency.at(ii).resize(arg_size);
            }
            for(int ii = 0; ii < arg_size-1; ii++){
                std::uniform_int_distribution<int> uniDist(ii+1, arg_size-1);
                int neighbours = poissonDist(gen);
                for(int jj = 0; jj < neighbours; jj++){
                    int k = uniDist(gen);
                    adjacency.at(ii).at(k) = 1;
                    adjacency.at(k).at(ii) = 1;
                }
            }

            makeList();
        }

    public:
        // create a network with normal distributed edges 
        void initNetNormal(int arg_size, int arg_meanDegree, double arg_width){
            N = arg_size;
            std::poisson_distribution<int> poissonDist(arg_meanDegree);
            std::normal_distribution<double> normalDist(0., arg_width);
            adjacency.clear();
            adjacency.resize(arg_size);
            for(int ii = 0; ii < arg_size; ii++){
                adjacency.at(ii).resize(arg_size);
            }
            for(int ii = 0; ii < arg_size-1; ii++){
                std::uniform_int_distribution<int> uniDist(ii+1, arg_size-1);
                int neighbours = poissonDist(gen);
                for(int jj = 0; jj < neighbours; jj++){
                    int k = ceil(abs(normalDist(gen)));
                    if(ii+k < arg_size){
                        adjacency.at(ii).at(ii+k) = 1;
                        adjacency.at(ii+k).at(ii) = 1;
                    }
                }
            }

            makeList();
        }

    public:
        // create an undirected Price model
        void initNetUndirectedPrice(int arg_size, int c, int a){
            N = arg_size;
            double r = double(c)/(double(c) + double(a));
            std::uniform_real_distribution<> real_uniDist(0., 1.);
            std::vector < int > targetList;
            adjacency.clear();
            adjacency.resize(N);
            for(int ii = 0; ii < arg_size; ii++){
                adjacency.at(ii).resize(N);
            }

            targetList.push_back(0);
            for(int ii = 0; ii < N; ii++){
                for(int jj = 0; jj < c; jj++){
                    if(real_uniDist(gen) < r){
                        std::uniform_int_distribution<int> uniDist(0, targetList.size()-1);
                        int k = targetList.at(uniDist(gen));
                        adjacency.at(ii).at(k) = 1;
                        adjacency.at(k).at(ii) = 1;
                        targetList.push_back(ii);
                        targetList.push_back(k);
                    }else{
                        std::uniform_int_distribution<int> uniDist(0, N-1);
                        int k = uniDist(gen);
                        adjacency.at(ii).at(k) = 1;
                        adjacency.at(k).at(ii) = 1;
                        targetList.push_back(ii);
                        targetList.push_back(k);
                    }
                }
            }

            makeList();
        }


    public:
        // create an undirected Price model with poisson distributed parameter c
        void initNetUndirectedPoissonPrice(int arg_size, int poissonParameter, int a){
            std::poisson_distribution<int> poissonDist(poissonParameter);
            N = arg_size;
            double r;
            int c;            
            std::uniform_real_distribution<> real_uniDist(0., 1.);
            std::vector < int > targetList;
            adjacency.clear();
            adjacency.resize(N);
            for(int ii = 0; ii < arg_size; ii++){
                adjacency.at(ii).resize(N);
            }

            targetList.push_back(0);
            for(int ii = 0; ii < N; ii++){
                c = poissonDist(gen);
                r =  double(c)/(double(c) + double(a));
                for(int jj = 0; jj < c; jj++){
                    if(real_uniDist(gen) < r){
                        std::uniform_int_distribution<int> uniDist(0, targetList.size()-1);
                        int k = targetList.at(uniDist(gen));
                        adjacency.at(ii).at(k) = 1;
                        adjacency.at(k).at(ii) = 1;
                        targetList.push_back(ii);
                        targetList.push_back(k);
                    }else{
                        std::uniform_int_distribution<int> uniDist(0, N-1);
                        int k = uniDist(gen);
                        adjacency.at(ii).at(k) = 1;
                        adjacency.at(k).at(ii) = 1;
                        targetList.push_back(ii);
                        targetList.push_back(k);
                    }
                }
            }

            makeList();
        }


    public: 
        // Initialize the starting setup for SIR dynamic
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
        // generate time step
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
        // print the current state of the infectious dynamic
        void printState(){
            int iS, iI, iR;
            iS = std::accumulate(S.begin(), S.end(), 0);
            iI = std::accumulate(I.begin(), I.end(), 0);
            iR = std::accumulate(R.begin(), R.end(), 0);
            std::cout<<iS<<"\t"<<iI<<"\t"<<iR<<"\t"<<N<<"\t"<<iS+iI+iR<<std::endl;
        }

    public:
        // print the adjacency matrix to a file
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
        // the ramdom number generator is of the Mersenne Twister type
        std::mt19937_64 gen{std::random_device{}()};
        // The number of nodes is N
        int N{0};
        // The infection probability is p
        double p{0};


    private:
        // create an adjacency list from the adjacency matrix
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
