#include <assert.h>
#include "algorithms.h"

Algorithm::Algorithm(/* args */)
{
    n = m = md = 0;
    deg = NULL;
	adj = NULL;
	datas = NULL;
}

Algorithm::~Algorithm()
{
    //printf("algorithm=%d, mincliquesize=%d\n", algorithm, mincliquesize);
    if (deg != NULL) delete[] deg; deg = NULL;
	if (adj != NULL) delete[] adj; adj = NULL;
	if (datas != NULL) delete[] datas; datas = NULL;
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

void Algorithm::read_graph(const char8 *str)
{
    printf("file: %s\n", str);
    boolean is_bin = false;
	clock_t stm = clock();
    if (strstr(str,".bin")) is_bin = true;
    if (is_bin) {
        FILE *in = fopen(str, "rb");
        if (in == NULL) {
            printf("No such file: %s\n", str);
            exit(1);
        }

		size_t FRead = 0;
		FRead = fread(&n, sizeof(int32), 1, in);
		FRead = fread(&m, sizeof(int32), 1, in);
		deg = new int32[n]();
		adj = new int32*[n];
		datas = new int32[m]();
		FRead = fread(deg, sizeof(int32), n, in);
        FRead = fread(datas, sizeof(int32), m, in);
		fclose(in);

		for (int32 i = 0, s = 0; i < n; ++i) { // Construct offs of all vertices
			adj[i] = datas + s;
			s += deg[i];
			md = deg[i] > md ? deg[i] : md;
		}
		printf("n=%d, m=%d, md=%d\n", n, m, md);
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
        //printf("n=%d, m=%d\n", n, m);
        assert(n > 0); assert(m > 0);
        vector<pair<int32,int32>> tempE; tempE.reserve(m);

        if (deg != NULL) exit(1);
        deg = new int32[n]();

        int32 u, v, cnt = 0;
        int32 maxid = 0;
        for (int32 i = 0; i < m; ++i) {
            char8 *r = fgets(line, 128, in);
            if (feof(in)) break;
            sscanf(line, "%d %d", &u, &v);
            //printf("u=%d, v=%d\n", u, v);
            if (u >= v) continue;
            assert(u < n && u >= 0);
            assert(v < n && v >= 0);
            tempE.emplace_back(u,v);
            deg[u]++; deg[v]++;
        }
        fclose(in);
        m = tempE.size();
        for (int32 i = 0; i < n; ++i) md = max(md, deg[i]);
        for (int32 i = 1; i < n; ++i) deg[i] += deg[i-1];
        for (int32 i = n-1; i > 0; --i) deg[i] = deg[i-1]; deg[0] = 0;

        adj = new int32*[n];
        datas = new int32[m*2]();
        for (int32 i = 0; i < m; ++i) {
            u = tempE[i].first;
            v = tempE[i].second;
            
            int32 pu = deg[u]++, pv = deg[v]++; 
            datas[pu] = v; datas[pv] = u;
        }
        for (int32 i = n-1; i > 0 ; --i) deg[i] -= deg[i-1];
        *adj = datas;
        for (int32 i = 1; i < n; ++i) adj[i] = adj[i-1] + deg[i-1]; 

        printf("n=%d, m=%d, md=%d\n", n, m * 2, md);
    }
	printf("Reading time: %lf s\n", double(clock()-stm)/CLOCKS_PER_SEC);

}


int32 Algorithm::core_decompsition(int32 *nodeset, int32 nodesize)
{
    boolean flag = nodeset == NULL ? true : false;
    if (core.empty()) core.resize(n, 0);
    int32 maxcore   = 0;
    int32 len       = nodesize;
    int32 *sequence = new int32[nodesize];
    int32 *bin      = new int32[md+1]();
    int32 *pos      = new int32[n];
    int32 *curdeg   = new int32[n];
    memcpy(curdeg, deg, sizeof(int32) * n);

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

    delete[] sequence;
    delete[] bin;
    delete[] pos;
    delete[] curdeg;
    return maxcore;
}
