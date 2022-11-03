#include <assert.h>
#include "algorithms.h"

Algorithm::Algorithm(/* args */)
{
    n = m = md = 0;
}

Algorithm::~Algorithm()
{
    //printf("algorithm=%d, mincliquesize=%d\n", algorithm, mincliquesize);
}

void Algorithm::testprintGraph()
{
    for (int32 i = 0; i < n; ++i) {
        printf("nbr[%d]: deg=%d\n", i, deg[i]);
        int32 d = deg[i];
        for (int32 j = 0; j < d; ++j) {
            printf("\t%d\n", adj[i][j]);
        }
    }
}

void Algorithm::read_graph(const char *str)
{
    printf("file: %s\n", str);
    bool is_bin = false;
    clock_t stm = clock();
    if (strstr(str,".bin")) is_bin = true;
    if (is_bin) {
        FILE *in = fopen(str, "rb");
        if (in == NULL) {
            printf("No such file: %s\n", str);
            exit(1);
        }
        if(fread(&n, sizeof(int32), 1, in)!=1) printf("err: read n!\n");
        if(fread(&m, sizeof(int32), 1, in)!=1) printf("err: read m!\n");
        deg.resize(n); adj.resize(n);
        for (int32 i = 0, s = 0; i < n; ++i)
            if (fread(&deg[i], sizeof(int32), 1, in)!=1)
                printf("err: read deg[%d]!\n",i);
        int u, d;
        for (int32 i = 0, s = 0; i < n; ++i) {
            d = deg[i];
            adj[i].reserve(d);
            for (int32 j = 0; j < d; j++) {
                if (fread(&u, sizeof(int32), 1, in)!=1)
                    printf("err: read adj[%d][%d]!\n",i, j);
                adj[i].emplace_back(u);
            }
            md = d > md ? d : md;
        }
        fclose(in);
        printf("n = %d, m = %d, maxdeg = %d\n", n, m, md);
    }
    else {
        FILE *in = fopen(str, "r");
        if (in == NULL) {
            printf("No such file: %s\n", str);
            exit(1);
        }
        char8 line[128];
        fgets(line, 128, in);
        if (sscanf(line, "%d %d", &n, &m) != 2) exit(1);
        assert(n > 0); assert(m > 0);
        vector<pair<int32,int32>> tempE;
        tempE.reserve(m);

        deg.resize(n);

        int32 u, v, cnt = 0;
        int32 maxid = 0;
        for (long i = 0; i < m; ++i) {
            char *r = fgets(line, 128, in);
            sscanf(line, "%d %d", &u, &v);
            if (u >= v) continue;
            assert(u < n && u >= 0);
            assert(v < n && v >= 0);
            tempE.emplace_back(u,v);
            deg[u]++; deg[v]++;
            if (feof(in)) break;
        }
        fclose(in);
        m = tempE.size();
        adj.resize(n);
        for (long i = 0; i < m; ++i) {
            u = tempE[i].first;
            v = tempE[i].second;
            adj[u].emplace_back(v);
            adj[v].emplace_back(u);
        }
        for (int32 i = 0; i < n; ++i)
            md = max(md, deg[i]);

        printf("n = %d, m = %d, maxdeg = %d\n", n, m * 2, md);
    }
    printf("Reading time: %lf s\n", double(clock()-stm)/CLOCKS_PER_SEC);
}

int Algorithm::core_decompsition(vector<int> &nodeset, int32 nodesize)
{
    bool flag = nodeset.empty() ? true : false;
    int32 maxcore = 0;
    int32 len     = nodesize;
    vector<int32> bin, pos, curdeg, sequence;
    pos.resize(n);
    bin.resize(md+1, 0);
    sequence.resize(n);
    curdeg = deg;
    if (core.empty()) core.resize(n,0);
    for (int32 i = 0; i < len; ++i) flag ? bin[curdeg[i]]++ : bin[curdeg[nodeset[i]]]++;
    //for (int32 i = 0; i <= md; ++i) printf("bin[%d]=%d\n", i, bin[i]);
    for (int32 i = 1; i <= md; ++i) bin[i] += bin[i-1];
    for (int32 i = md; i > 0; --i) bin[i] = bin[i-1]; bin[0] = 0;

    //for (int32 i = 0; i <= md; ++i) printf("offs[%d]=%d\n", i, bin[i]);

    for (int32 i = 0; i < len; ++i) {
        int32 v = flag ? i : nodeset[i];
        int32 posv = bin[deg[v]]++;
        sequence[posv] = v;
        pos[v] = posv;
    }

    int32 k = 0;
    for (int32 i = 0; i < len; ++i) {
        int32 v = sequence[i];
        int32 d = deg[v];
        k = max(k, curdeg[v]);
        maxcore = max(k,maxcore);
        flag ? i : nodeset[i] = v;
        core[v] = k;
        for (int32 j = 0; j < d; ++j) {
            int32 w = adj[v][j];
            int32 dw = curdeg[w]--;
            if (dw > k) {
                int32 posw = pos[w];
                int32 pdws = bin[dw-1]++;
                if (posw != pdws) {
                    sequence[posw] = sequence[pdws];
                    pos[sequence[posw]] = posw;
                    sequence[pdws] = w;
                    pos[w] = pdws;
                }
            }
        }
    }
    return maxcore;
}

int Algorithm::coloring(vector<int> &nodeset, int nodesize)
{
    int color_nums = 0, len = nodesize;
    vector<int> bin(md+1), sequence(len);

    //core_decompsition(nodeset, nodesize);
    //for (int i = 0; i < len; ++i) sequence[i] = nodeset[i];

    for (int i = 0; i < len; ++i) bin[deg[nodeset[i]]]++;
    for (int i = 1; i <= md; ++i) bin[i] += bin[i-1];
    for (int i = md; i > 0; --i) bin[i] = bin[i-1]; bin[0] = 0;
    for (int i = 0; i < len; ++i) {
        int v = nodeset[i];
        int posi = bin[deg[v]]++;
        sequence[posi] = v;
    }
    if (colors.empty()) colors.resize(n, 0);
    for (int i = len-1; i >= 0; --i) {
        int maxc = -1, curc = -1;
        int u = sequence[i];
        int d = deg[u];
        for (int j = 0; j < d; ++j) {
            int v = adj[u][j];
            int cv = colors[v];
            if (cv >= 0) {bin[cv]++; maxc = max(maxc, cv);}
        }
        for (int j = 0; j <= maxc; ++j) {
            if (bin[j] == 0 && curc == -1) curc = j;
            bin[j] = 0;
        }
        if (curc == -1) {
            color_nums = max(color_nums, ++maxc + 1);
            colors[u] = maxc;
        }
        else colors[u] = curc;
    }
    printf("Color nums: %d\n", color_nums);
    return color_nums;
}
