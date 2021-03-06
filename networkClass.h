// /////////////////////////////////////////////////////////////// //
// Network Class. Basic Network Class for simple dynamics          //
// (c) 2021, J. Neidhart                                           //
// /////////////////////////////////////////////////////////////// //
#ifndef NETCLASS1_H
#define NETCLASS1_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <numeric>
#include <cmath>

#include "sirClass.h"


class cl_network{
    public:
        cl_network(){};
        ~cl_network(){};

        // Adjaceny Matrix
        std::vector < std::vector < double > > adjacency;

        // Adjacency List
        std::vector < std::vector < int > > adList;

        // Degree Vector
        std::vector < int > degrees;

        // Mean Degree
        double meanDegree;

    public:
        // create a network with Poisson distributed number of random edges
        void initNetRandom(int arg_size, int arg_meanDegree){
            N = arg_size;
            std::poisson_distribution<int> poissonDist(arg_meanDegree);
            adjacency.clear();
            adjacency.resize(arg_size);
            for(auto & a: adjacency){
                a.resize(arg_size);
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
            for(auto & a: adjacency){
                a.resize(arg_size);
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
            for(auto & a: adjacency){
                a.resize(arg_size);
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
            for(auto & a: adjacency){
                a.resize(arg_size);
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
        // create a Watts-Strogatz Graph Model
        void initNetWattsStrogatz(int arg_size, int arg_edges, double arg_rewiringProbability){
            N = arg_size;
            std::uniform_real_distribution<> real_uniDist(0., 1.);
            std::uniform_int_distribution<int> uniDist(0, N-1);
            int k = arg_edges;
            double beta = arg_rewiringProbability;
            adjacency.clear();
            adjacency.resize(N);
            for(auto & a: adjacency){
                a.resize(arg_size);
            }


            // connect eacch node to the k next nodes (for strict Watts Strogatz it should be connected to k/2 before
            // and k/2 after)
            for(int ii = 0; ii < N; ii++){
                for(int jj = 1; jj <= k; jj++){
                    int h = (ii+jj)%N;
                    adjacency.at(ii).at(h) = 1;
                    adjacency.at(h).at(ii) = 1;
                }
            }
            // Go through all edges and rewire each edge with probability beta
            for(int ii = 0; ii < N; ii++){
                for(int jj = ii; jj < N; jj++){
                    if(adjacency.at(ii).at(jj)){
                        // rewire with probability beta
                        if(real_uniDist(gen) < beta){
                            //reset edge
                            adjacency.at(ii).at(jj) = 0;
                            adjacency.at(jj).at(ii) = 0;
                            int h = ii;
                            while(h == ii){
                                h = uniDist(gen);
                            }
                            adjacency.at(ii).at(h) = 1;
                            adjacency.at(h).at(ii) = 1;
                        }
                    }
                }
            }
            makeList();
        }


    public:
        // calculate the degree vector and the mean degree
        void calculateDegrees(){
            meanDegree = 0;
            degrees.clear();
            degrees.resize(N);
            for(int ii = 0; ii < N; ii++){
                for(int jj = 0; jj < N; jj++){
                    degrees.at(ii) += adjacency.at(ii).at(jj);
                }
                meanDegree += degrees.at(ii);
            }
            meanDegree /= (double)N;
        }


    private:
        // calculate distance of all nodes to given node. This is a direct implementation of Newman, Ch.10.3.3
        std::vector < int > distances(int initialNode){
            std::vector < int > queue;
            std::vector < int > distList;
            int readPointer = 0;
            int writePointer = 1;
            queue.clear();
            queue.resize(N);
            distList.clear();
            queue.at(0) = initialNode;
            distList.resize(N, -1);
            distList.at(initialNode) = 0;
            while(true){
                if(readPointer == writePointer){
                    break;
                }
                int label = queue.at(readPointer);
                readPointer++;
                int d = distList.at(label);
                for(const auto & a: adList.at(label)){
                    if(distList.at(a) == -1){
                        distList.at(a) = d+1;
                        queue.at(writePointer) = a;
                        writePointer++;
                    }
                }
            }
            return(distList);
        }


    public:
        // Implement Breadth First Search for finding shortest paths and Diameter, which is the longest shortest path
        void calculatePathStatistics(){
            //    for(auto &a: distances(0)){
            double sum{0.};
            int max{0};
            for(int ii = 0; ii < N; ii++){ 
                for(auto &a: distances(ii)){
                    sum += double(a);
                    if(a > max){
                        max = a;        
                    }
                }
            }
            meanDistance = (sum/double(N))/double(N);
            maxDistance = max;
            std::cout<<"Mean Distance: "<<meanDistance<<";  Diameter: "<<maxDistance<<std::endl;
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
        double p{0.};
        // mean Distance between nodes
        double meanDistance{0.};
        // Diameter i.e. maximum distance between nodes
        int maxDistance{0};

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
