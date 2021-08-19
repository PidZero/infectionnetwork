// /////////////////////////////////////////////////////////////// //
// SIR Class. Basic Class for simple SIR dynamics on a network     //
// (c) 2021, J. Neidhart                                           //
// /////////////////////////////////////////////////////////////// //

# ifndef SIRCLASS_H
# define SIRCLASS_H

#include "networkClass.h"
#include <vector>

class cl_sir{
    public:

        // Constructor
        cl_sir(std::vector < std::vector < int > > arg_adList) : adList(arg_adList){
            N = adList.size();
        };

        // destructor
        ~cl_sir(){};

    private:

        // Adjacency List as Vector of vectors
        std::vector < std:: vector <int> > adList;

        // the ramdom number generator is of the Mersenne Twister type
        std::mt19937_64 gen{std::random_device{}()};
        
        // The number of nodes is N
        int N{0};
        
        // The infection probability is p
        double p{0.};


        // Vector to store the labels of the susceptible
        std::vector < bool > S;

        // Vector to store the labels of the resistant
        std::vector < bool > R;

        // Vector to store the labels of the infected
        std::vector < bool > I;

        // Numbers of the resistant, infected and susceptible
        int numberR{0};
        int numberI{0};
        int numberS{0};


    public: 
        // Initialize the starting setup for SIR dynamic
        void initInfection(int arg_I, double arg_infectionProb){
            int k; 
            p = arg_infectionProb;https://manjaro.org/
                std::uniform_int_distribution<int> uniDist(0, N-1);
            S.clear();
            I.clear();
            R.clear();
            S.resize(N, 1);
            I.resize(N, 0);
            R.resize(N, 0);
            numberI = arg_I;
            numberS = N-arg_I;
            numberR = 0;
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
        // do statistics on a previously initialized Network
        void statistics(int arg_repetitions){
            std::vector < double > timeline;
            int time{0};
            for(int ii = 0; ii < arg_repetitions; ii++){
                initInfection(3, 0.3);
                time = 0;
                while(numberI > 0){
                    if(timeline.size()<= time){
                        timeline.push_back(numberI);
                    }else{
                        timeline.at(time) += numberI;
                    }
                    timestep();
                    time++;
                }
            }
            for(auto & a:timeline){
                std::cout<<double(a)/double(arg_repetitions)<<std::endl;
            }
        }



};

#endif
