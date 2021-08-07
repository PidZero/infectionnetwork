#include "networkClass.h"

int main(){
    cl_network net;
    net.initNetWattsStrogatz(100, 2, .1);
    net.initInfection(3, 0.3);
    while(net.numberI > 0){
        net.printState();
        net.timestep();
    }
    net.printState();
    net.printAdjacency("net1.csv");
    return(0);
}
