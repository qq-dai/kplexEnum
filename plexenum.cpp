#include "plexenum.h"

KplexEnum::KplexEnum(/* args */)
{
}

KplexEnum::~KplexEnum()
{
}

void KplexEnum::setParameters(int32 argc, char *argv[]) {
    for (int32 i = 1; i < argc; ++i) {
        char *p = argv[i];
        if (strstr(p, "-q=")) {
            minsize = atoi(p+3);
        }
        if (strstr(p, "-k=")) {
            k = atoi(p+3);
        }
    }
    printf("minsize=%d, k=%d\n", minsize, k);
}

void KplexEnum::initIndexAndFwd(int32 *nodeset)
{
    cuhash.resize(n);
    fwdadj.resize(n);
    vector<boolean> visited(n, false);
    for (int32 i = 0; i < n; ++i) {
        int32 v = nodeset[i];
        int32 d = deg[v];
        int32 cv = core[v];
        cuhash[v].reserve(d);
        fwdadj[v].reserve(cv);
        for (int32 j = 0; j < d; ++j) {
            int32 u = adj[v][j];
            if (!visited[u]) 
                fwdadj[v].emplace_back(u);
            cuhash[v].insert(u);
        }
        visited[v] = true;
    }
}

void KplexEnum::initRoot(int32 v, vector<int32> &R, Node &P, vector<boolean> &visited)
{
    P.clear(); visited[v] = true;
    excluded.clear();
    vector<vector<int32>> temp(2);
    int32 p = 0, q = 1;
    boolean it = true;
    temp[0].reserve(fwdadj[v].size());
    temp[1].reserve(fwdadj[v].size());
    for (auto u : fwdadj[v]) {
        visited[u] = true;
        temp[0].emplace_back(u);
    }
    while (it) {
        vector<int32> &atemp = temp[p];
        vector<int32> &btemp = temp[q];
        btemp.clear();
        it = false;
        for (size_t i = 0; i < atemp.size(); ++i) {
            int32 u = atemp[i];
            int32 d = 0;
            int32 maxc = 0, clnums = 0;
            for (size_t j = 0; j < atemp.size(); ++j) {
                int32 w = atemp[j];
                if (u != w && is_nbr(u,w)) d++;
            }
            if (d + 2*k >= minsize)
                btemp.emplace_back(u);
        }
        p = q; q = 1-p;
        if (btemp.size() != atemp.size()) it = true;
    }
    if (temp[p].size()+k+1 < minsize) {
        for (auto u : fwdadj[v]) visited[u] = false;
        return;
    }
    for (auto u : temp[p]) P.cinsert(u,0);
    int32 len = P.size;
    for (int32 i = 0; i < deg[v]; ++i) {
        int32 u = adj[v][i];
        if (core[u]+k >= minsize && !visited[u]){
            int32 d = 0;
            visited[u] = true;
            for (int32 j = 0; j < len; ++j) {
                int32 w = P.nodes[j].first;
                if (is_nbr(u,w)) d++;
            }
            if (d + 2*k >= minsize)
                excluded.emplace_back(u);
        }
    }
    int32 exlen = excluded.size();
    for (int32 i = 0; i < len; ++i) {
        int32 u = P.nodes[i].first;
        int32 du = deg[u];
        for (int32 j = 0; j < du; ++j) {
            int32 w = adj[u][j];
            if (core[w] + k < minsize) continue;
            if (!visited[w] && !exclusion[w]) {
                visited[w] = true;
                P.cinsert(w, 1);
                noadj[w*k] = 0;
            }
            if (!visited[w] && exclusion[w]) {
                excluded.emplace_back(w);
                visited[w] = true;
            }
        }
    }
    for (int32 i = exlen; i < excluded.size(); ++i) {
        int32 u = excluded[i];
        int32 d = 0;
        visited[u] = false;
        for (int32 j = 0; j < len; ++j) {
            int32 w = P.nodes[j].first;
            if (is_nbr(u,w)) d++;
        }
        if (d + 2*k >= minsize+2)
            excluded[exlen++] = u;
    }
    excluded.resize(exlen);

    int32 len1 = len;
    for (int32 i = len1; i < P.size; ++i) {
        int32 u = P.nodes[i].first;
        int32 d = 0;
        visited[u] = false;
        for (int32 j = 0; j < len1; ++j) {
            int32 w = P.nodes[j].first;
            if (is_nbr(u,w)) d++;
        }
        if (d + 2*k >= minsize + 2){
            //P.counts[len]  = d;
            P.nodes[len++].first = u;
        }
    }
    P.size = len;
    for (int32 i = 0; i < len1; ++i) {
        int32 u = P.nodes[i].first;
        visited[u] = false;
    }

    visited[v] = false;
    R[0] = v; counters[0] = 1;
    //for (auto u : fwdadj[v]) visited[u] = false;
    for (int32 i = 0; i < deg[v]; ++i) visited[adj[v][i]] = false;
}

void KplexEnum::run() {
    clock_t tm = clock();
    enums();
    printf("All time: %f sec\n", double(clock()-tm)/CLOCKS_PER_SEC);
}

void KplexEnum::enums() {
    int32 *nodeset = new int32[n];
    for (int32 i = 0; i < n; ++i) nodeset[i] = i;
    maxcore = core_decompsition(nodeset, n);
    initIndexAndFwd(nodeset);
    Node P; P.reserve(n);
    vector<int32>  R(maxcore+k);
    vector<boolean> visited(n,false);

    counters.resize(maxcore+k,0);
    _counters.reserve(maxcore+k);
    noadj.resize(n*k);
    exclusion.resize(n, false);
    excluded.reserve(n);

    nonbrcnts.resize(k+1);
    pushnonbrs.resize(2*k);

    subdeg.resize(n,0);
    subadj.resize(n);
    
    clock_t st = clock();
    int32 computed=0;
    int32 maxp = 0;
    for (int32 i = 0; i < n; ++i) {
        int32 v = nodeset[i];
        if (fwdadj[v].size() + k < minsize) continue;
        initRoot(v, R, P, visited);
        exclusion[v] = true;
        if (P.size+1 < minsize) continue;
        maxp = max(maxp, P.size);
        setSubadj(R, 1, P);
        if (subdeg[v]+k < minsize) continue;
        //recursion(R, 1, P);
        totalcnt = k-1;
        recursion(R, 1, P);
        if (computed++ == n) break;
    }
    printf("Number of plexes    : %ld \n", kplexnums);
    printf("Number of iterations: %ld \n", iterations);
    printf("Runtime             : %f sec\n", double(clock() - st) / CLOCKS_PER_SEC);
    //printf("MaxPlex             : %d \n", maxPlex);
    //printf("MaxCandidate        : %d \n", maxp);
    delete[] nodeset;
}

void KplexEnum::recursion(vector<int32> &R, int32 rsize, Node &P)
{
    if (rsize + P.size < minsize + P.exsize) return;
    iterations++;
    if (P.size == 0 && rsize >= 2 * k - 1) {
        if (isMaximal(R, P, rsize))
        {
            maxPlex = max(maxPlex,rsize);
            kplexnums++;
        }
        return;
    }
    if (P.size - P.exsize == 0) return;

    int32 mini = -1;
    int32 mind = n;
    for (int32 i = 0; i < rsize; ++i) {
        int32 u = R[i];
        if (mind == n || subdeg[u] < mind) {
            mind = subdeg[u];
            mini = i;
        }
    }
    if (mind+k < minsize) return;
    if (mind+k+P.exsize < rsize+P.size) {
        int32 v = R[mini];
        int32 pi = -1;
        int32 cnt = 0;
        mind = n;
        for (int32 i = P.exsize; i < P.size; ++i) {
            int32 u = P.nodes[i].first;
            int32 c = P.nodes[i].second;
            int32 d = subdeg[u];
            if (!is_nbr(u,v)) {
                if (mind == d && cnt < c) {
                    pi = i; cnt = c;
                }
                else if (mind > d) {
                    mind = d;
                    pi = i; cnt = c;
                }
            }
        }
        Node _P; _P.resize(P.size);
        vector<int32> test; 
        test.reserve(P.size - P.exsize);
        addToPlex(P, pi, R, rsize);
        _counters.clear();
        for (int32 j = 0; j < rsize; ++j) 
            _counters.emplace_back(counters[j]);
        if (updateSet(R, rsize+1, _counters, P, P.exsize, _P, test)){
            for (auto u : test) removeNodes(u);
            recursion(R, rsize+1, _P);
            for (auto u : test) addNodes(u);
        }
        resetCounts(P, pi);

        int32 u = P.nodes[pi].first; 
        for (int32 j = 0; j < P.size; ++j)
            _P.nodes[j] = P.nodes[j];
        _P.exsize = P.exsize+1;
        _P.size = P.size;
        exchanges(_P, P.exsize, pi);
        removeNodes(u);
        recursion(R, rsize, _P);
        addNodes(u);
        return;
    }

    int32 minicnt = 0;
    for (int32 i = P.exsize; i < P.size; ++i) {
        int32 u = P.nodes[i].first;
        int32 c = P.nodes[i].second;
        if (mind == n || subdeg[u] < mind) {
            mind = subdeg[u];
            mini = i + rsize;
            minicnt = c;
        }
        else if (subdeg[u] == mind && minicnt < c) {
            mini = i + rsize;
            minicnt = c;
        }
    }

    //printf("deg[%d]=%d, mind=%d\n", P.nodes[mini].first, subdeg[P.nodes[mini].first], mind);
    if (mind+k+P.exsize >= rsize+P.size) {
        if (rsize+P.size-P.exsize >= 2*k-1) {
            if (isMaximal(R, P, rsize))
            {
                maxPlex = max(maxPlex,rsize+P.size-P.exsize);
                kplexnums++;
            }
        }
        return;
    }

    Node _P; _P.resize(P.size);
    int32 i = mini-rsize;
    if (mind+k >= minsize) {
        vector<int32> test; 
        test.reserve(P.size - P.exsize);
        addToPlex(P, i, R, rsize);
        _counters.clear();
        for (int32 j = 0; j < rsize; ++j) 
            _counters.emplace_back(counters[j]);
        if (updateSet(R, rsize+1, _counters, P, P.exsize, _P, test)){
            for (auto u : test) removeNodes(u);
            recursion(R, rsize+1, _P);
            for (auto u : test) addNodes(u);
        }
        resetCounts(P, i);
    }

    int32 u = P.nodes[i].first;
    for (int32 j = 0; j < P.size; ++j)
        _P.nodes[j] = P.nodes[j];
    _P.exsize = P.exsize+1;
    _P.size = P.size;
    exchanges(_P, P.exsize, i);
    removeNodes(u);
    recursion(R, rsize, _P);
    addNodes(u);
}


boolean KplexEnum::updateSet(vector<int32> &R, int32 rsize, vector<int32> &_counters, Node &P, int32 id, Node &_P, vector<int32> &Rnodes)
{
    int32 rl  = rsize-1;
    int32 crl = counters[rl];
    int32 u   = R[rl];
    int32 cnt = 0;
    int32 len1 = 0;
    for (int32 i = 0; i < k; ++i) nonbrcnts[i] = 0;
    for (int32 i = 0; i < P.exsize; ++i) {
        int32 v = P.nodes[i].first;
        int32 c = P.nodes[i].second;
        int32 cn = c;
        boolean ok = true;

        if (!is_nbr(v,u)) {
            if (crl+1 > k || c+1 >= k) continue;
            noadj[v*k+cn] = rl;
            cn += 1;
        }

        for (int32 j = v*k; j < v*k + c; ++j) {
            int32 wid = noadj[j];
            if (counters[wid]+1 > k) {ok = false; break;}
        }
        if (ok)  {
            _P.nodes[len1].first = v;
            _P.nodes[len1].second = cn;
            len1++;
            //_P.einsert(v, cn);
        }
    }
    int32 len2 = len1;
    if (len1 > id) {
        printf("id=%d, P.exsize=%d\n",id, len1);
    }
    int32 cnt1 = 0;
    for (int32 i = id; i < P.size; ++i) {
        int32 v = P.nodes[i].first;
        int32 c = P.nodes[i].second;
        int32 cn = c;
        boolean ok = true;
        if (v == u) continue;
        if (!is_nbr(v,u)) {
            if (crl+1 > k || c+1 >= k)  {
                Rnodes.emplace_back(v);
                continue;
            }
            noadj[v*k+cn] = rl;
            cn += 1;
        }
        for (int32 j = v*k; j < v*k + c; ++j) {
            if (counters[noadj[j]]+1 > k) {
            	ok = false; break;
            }
        }
        if (ok) {
            _P.nodes[len2].first = v;
            _P.nodes[len2++].second = cn;
            _P.cinsert(v, cn);
            // prunning
            if (cn == c) {
                cnt1++; 
                nonbrcnts[c]++;
                if (c == 0) {cnt++;}
                else {
                    int32 maxcounts = -1;
                    int32 maxwid = 0;
                    for (int32 j = v*k; j < v*k + c; ++j) {
                        int32 wid = noadj[j];
                        int32 ct = _counters[wid];
                        if (ct > maxcounts) {
                            maxwid = wid; 
                            maxcounts = ct;
                        }
                    }
                    if (maxcounts < k) {
                        _counters[maxwid]++;
                        cnt++;
                    }
                }
            }
        } else Rnodes.emplace_back(v); 
    }
    _P.exsize = len1;
    _P.size = len2;

    if (rsize+k+cnt < minsize+crl) return false;
    int32 nbrs = nonbrcnts[0];
    int32 s = totalcnt + (crl-k);
    for (int32 i = 1; i < k; ++i) {
        int32 c = nonbrcnts[i];
        if (s < i * c) {
            nbrs += s/i;
            break;
        }
        nbrs += c;
        s -= i*c;
    }
    if (nbrs+rsize+k < minsize+crl) return false;

    return true;
}

boolean KplexEnum::isMaximal(vector<int32> &R, Node &P, int32 rsize) {
    for (int32 i = P.exsize; i < P.size; ++i)
        R[rsize++] = P.nodes[i].first;
    //printf("rsize=%d\n", rsize);
    for (int32 i = 0; i < excluded.size(); ++i) {
        int32 u = excluded[i];
        int32 cnt = 0;
        boolean ok = true;
        for (int32 j = 0; j < rsize; ++j) {
            int32 w = R[j];
            if (!is_nbr(u, w))
                if (++cnt >= k || subdeg[w]+k <= rsize) {
                    ok = false; break;
                }
        }
        if (ok) return false;
    }
    for (int32 i = 0; i < P.exsize; ++i) {
        int32 u = P.nodes[i].first;
        int32 cnt = 0;
        boolean ok = true;
        for (int32 j = 0; j < rsize; ++j) {
            int32 w = R[j];
            if (!is_nbr(u, w))
                if (++cnt >= k || subdeg[w]+k <= rsize) {
                    ok = false; break;
                }
        }
        if (ok) return false;
    }
    return true;
}

boolean KplexEnum::isMaximal(vector<int32> &R, int32 rsize)
{
    for (int32 i = 0; i < excluded.size(); ++i) {
        int32 u = excluded[i];
        int32 cnt = 0;
        boolean ok = true;
        for (int32 j = 0; j < rsize; ++j) {
            if (!is_nbr(u, R[j]))
                if (++cnt >= k || counters[j]+1 > k) {
                    ok = false; break;
                }
        }
        if (ok) return false;
    }
    return true;
}

void KplexEnum::addToPlex(Node &P, int32 &id, vector<int32> &R, int32 rsize)
{
    int32 v = P.nodes[id].first;
    int32 c = P.nodes[id].second;
    R[rsize] = v;
    for (int32 i = v*k; i < v*k + c; ++i) {
        int32 uid = noadj[i];
        counters[uid] += 1;
    }
    counters[rsize] = c+1;
    totalcnt += (k-2*c-1);
}

void KplexEnum::resetCounts(Node &P, int32 &id)
{
    int32 v = P.nodes[id].first;
    int32 c = P.nodes[id].second;
    for (int32 i = v*k; i < v*k + c; ++i) {
        int32 uid = noadj[i];
        counters[uid] -= 1;
    }
    totalcnt -= (k-2*c-1);
}

inline boolean KplexEnum::is_nbr(int32 u, int32 v) {
    return cuhash[u].find(v);
}

void KplexEnum::sortNodes(Node &P, int32 s, int32 t) {
    ta.clear(); tb.clear();
    ta.resize(t-s); tb.resize(t-s);
    for (int32 i = s; i < t; ++i) {
        ta[i-s] = P.nodes[i].first;
        tb[i-s] = P.nodes[i].second;
    }

    sbin.clear();
    sbin.resize(k+1, 0);
    for (auto c : tb) sbin[c]++;

    for (int32 j ,i = k, l = P.exsize; i >= 0 ; --i) {
        j = sbin[i];
        sbin[i] = l;
        l += j;
    }
    assert(sbin[k] <= t);
    for (int32 i = 0; i < ta.size(); ++i) {
        int32 u = ta[i];
        int32 c = tb[i];
        int32 pos = sbin[c]++;
        P.nodes[pos].first = u;
        P.nodes[pos].second = c;
    }
}

void KplexEnum::setSubadj(vector<int32> R, int32 rsize, Node &P)
{
    int32 v = R[0];
    subadj[v].clear();
    for (int32 i = P.exsize; i < P.size; ++i) {
        int32 u = P.nodes[i].first;
        subadj[u].clear();
    }
    for (int32 i = P.exsize; i < P.size; ++i) {
        int32 u = P.nodes[i].first;
        if (P.nodes[i].second == 0) {
            subadj[u].emplace_back(v);
            subadj[v].emplace_back(u);
        }
        for (int32 j = i+1; j < P.size; ++j) {
            int32 w = P.nodes[j].first;
            if (is_nbr(u,w)) {
                subadj[u].emplace_back(w);
                subadj[w].emplace_back(u);
            }
        }
        subdeg[u] = subadj[u].size();
    }
    subdeg[v] = subadj[v].size();

    // test
    for (int32 i = P.exsize; i < P.size; ++i) {
        int32 u = P.nodes[i].first;
        int32 d = 0;
        for (int32 j = P.exsize; j < P.size; ++j) {
            int32 w = P.nodes[j].first;
            if (is_nbr(u,w)) ++d;
        }
        if (is_nbr(u,v)) ++d;
        if (d != subdeg[u]) {
            printf("v=%d, u=%d\n",v,u);
        }
        assert(d == subdeg[u]);
    }
}