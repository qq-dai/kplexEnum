#include "plexenum.h"

int main(int argc, char *argv[])
{
    Algorithm *mc = new KplexEnum();
    mc->read_graph(argv[1]);
    mc->setParameters(argc, argv);
    //mc->testprintGraph();
    mc->run();
    delete mc;
    return 1;
}