#include "networkClass.h"

int main(){
    cl_network net;
    net.initNetNormal(1000, 5, 15);
    net.initInfection(3, 0.1);
    while(net.numberI > 0){
        net.printState();
        net.timestep();
    }
    net.printState();
    net.printAdjacency("net1.csv");
    return(0);
}
