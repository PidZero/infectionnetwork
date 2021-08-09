#include "networkClass.h"

int main(){
    cl_network net;
   net.initNetWattsStrogatz(100, 4, .25);
    net.initInfection(3, 0.3);
    while(net.numberI > 0){
        net.printState();
        net.timestep();
    }
    net.printState();
    net.printAdjacency("net1.csv");
    net.calculateDegrees();
    std::cout<<"Mean Degree: "<<net.meanDegree<<std::endl;
    return(0);
}
