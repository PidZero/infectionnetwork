#include "networkClass.h"

int main(){
    cl_network net;
   net.initNetWattsStrogatz(1000, 3, .25);
    net.initInfection(8, 0.3);
    while(net.numberI > 0){
        net.printState();
        net.timestep();
    }
    net.printState();
    net.printAdjacency("net1.csv");
    net.calculateDegrees();
    net.calculatePathStatistics();
    std::cout<<"Mean Degree: "<<net.meanDegree<<std::endl;

    return(0);
}
