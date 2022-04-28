#pragma once
#include "algorithms.h"
#include <cstring>
#include <omp.h>

typedef struct _Graph
{
    int32 nodes;
    vector<vector<int32>> adj;
    vector<int32> excluded;
    vector<int32> nonbrcnts;
}Graph;

typedef struct _Vertxs
{
    vector<int32> counters;
    vector<int32> noadj;
    vector<int32> maps;
    vector<int32> subD;
    int32 totalcnt = 0;
    long kplexnums = 0;
    long iterations = 0;
}Vertxs;

#define ParTASK 5

class Parallel_Enum : public Algorithm
{
private:
    int32 k = 2;
    int32 maxPlex = 0;
    int32 minsize = 3;
    int32 maxcore = 0;
    vector<vector<int32>> fwdadj;
    //Table tbl;
    vector<CuckooHash> cuhash;
    vector<int32> VtxPos;

    int32 threads = 1;
    int32 tasks = 0;
    int32 overtasks = 0;
    vector<long> lcnt_plx, lcnt_its;
    bool division = false;

public:
    Parallel_Enum(/* args */);
    ~Parallel_Enum();

    void setParameters(int32 argc, char *argv[]);
    void initIndexAndFwd(vector<int> &nodeset);
    void initRoot(Graph &subG, Vertxs &sVtx, int32 v, vector<int32> &R, Node &P, vector<boolean> &visited);
    void setSubadj(Graph &subG, Vertxs &sVtx, vector<int32> &R, int32 rsize, Node &P);
    boolean isMaximal(Graph &subG, Vertxs &sVtx, vector<int32> &R, Node &P, int32 rsize);
    void addToPlex(Vertxs &sVtx, Node &P, int32 &id, vector<int32> &R, int32 rsize);
    void resetCounts(Vertxs &sVtx, Node &P, int32 &id);
    boolean updateSet(Graph &subG, Vertxs &sVtx, vector<int32> &R, int32 rsize, vector<int32> &_counters, Node &P, int32 id, Node &_P, vector<int32> &delnodes);
    void recursionCalls(Graph &subG, Vertxs &sVtx, vector<int32> &R, int32 rsize, Node &P, int32 pi);
    void recursion(Graph &subG, Vertxs &sVtx, vector<int32> &R, int32 rsize, Node &P);
    void recursion_seq(Graph &subG, Vertxs &sVtx, vector<int32> &R, int32 rsize, Node &P);
    void delNodes(vector<int32> &nbrs, vector<int32> &subD) {
        for (auto &u : nbrs) subD[u]--;
    }
    void addNodes(vector<int32> &nbrs, vector<int32> &subD) {
        for (auto &u : nbrs) subD[u]++;
    }
    void exchanges(Node &P, int32 s, int32 t) {
        if (s == t) return;
        upair a = P.nodes[s];
        P.nodes[s] = P.nodes[t];
        P.nodes[t] = a;
    }

    void penum();
    void run();
    inline bool is_nbr(int32 u, int32 v) {return cuhash[u].find(v);}
};
