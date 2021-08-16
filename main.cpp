// /////////////////////////////////////////////////////////////// //
// Network Class example main.                                     //
// (c) 2021, J. Neidhart                                           //
// /////////////////////////////////////////////////////////////// //

#include "networkClass.h"



int main(){
    cl_network net;
   net.initNetWattsStrogatz(100, 4, .25);
    /*while(net.numberI > 0){
        net.printState();
        net.timestep();
    }
    net.printState();
    net.printAdjacency("net1.csv");
    net.calculateDegrees();
    net.calculatePathStatistics();
    std::cout<<"Mean Degree: "<<net.meanDegree<<std::endl;

*/
   net.statistics(1000);
    return(0);
}
