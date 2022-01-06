#include "algorithms.h"
#include <cstring>

class KplexEnum : public Algorithm
{
private:
    /* data */
    int32 k         = 2;
    int32 minsize   = 3;
    int32 maxcore   = 0;
    int32 maxPlex   = 0;
    long64  kplexnums = 0;
    long64 iterations = 0;

    vector<vector<int32>> fwdadj;

    //Table tbl;
    vector<int32> counters;
    vector<int32> _counters;
    vector<int32> excluded;
    vector<int32> noadj;
    vector<boolean> exclusion;
    vector<CuckooHash> cuhash;

    vector<int32> pushnonbrs;

    vector<int32> sbin;
    vector<int32> ta, tb;

    vector<int32> subdeg;
    vector<vector<int32>> subadj;

    vector<int32> nonbrcnts;
    int32 totalcnt;

public:
    KplexEnum(/* args */);
    ~KplexEnum();
    void setParameters(int32 argc, char8 *argv[]);
    void initIndexAndFwd(int32 *nodeset);

    void run();
    void enums();
    void initRoot(int32 v, vector<int32> &R, Node &P, vector<boolean> &visited);
    void recursion(vector<int32> &R, int32 rsize, Node &P);
    void addToPlex(Node &P, int32 &id, vector<int32> &R, int32 rsize);
    void setCounts(vector<int32> &cnts, int32 &size);
    void resetCounts(Node &P, int32 &id);
    boolean updateSet(vector<int32> &R, int32 rsize, vector<int32> &_counters, Node &P, int32 id, Node &_P, vector<int32> &Rnodes);
    boolean isMaximal(vector<int32> &R, int32 rsize);
    boolean isMaximal(vector<int32> &R, Node &P, int32 rsize);
    inline boolean is_nbr(int32 u, int32 v);

    void sortNodes(Node &P, int32 s, int32 t);
    void setSubadj(vector<int32> R, int32 rsize, Node &P);

    void removeNodes(const int32 &u) {
        for (auto u : subadj[u])
            subdeg[u]--;
    }
    void addNodes(const int32 &u) {
        for (auto u : subadj[u])
            subdeg[u]++;
    }

    void exchanges(Node &P, int32 s, int32 t) {
        if (s == t) return;
        //assert(s < P.size && s >= P.exsize);
        //if (t >= P.size || t < P.exsize) {
        //    printf("P.size=%d, s=%d, t=%d\n",P.size,s,t);
        //}
        //assert(t < P.size && t >= P.exsize);
        upair a = P.nodes[s];
        P.nodes[s] = P.nodes[t];
        P.nodes[t] = a;
    }
};
