#include "parallel_plex.h"

int main(int argc, char *argv[])
{   
    Algorithm *mc = new Parallel_Enum();

    mc->read_graph(argv[1]);
    mc->setParameters(argc, argv);
    mc->run();
    
    delete mc;
    return 1;
}