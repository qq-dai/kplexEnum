#include "parallel_plex.h"
thread_local long cntT = 0;
thread_local long plexcntT = 0;

Parallel_Enum::Parallel_Enum(/* args */)
{
}

Parallel_Enum::~Parallel_Enum()
{
}

void Parallel_Enum::setParameters(int32 argc, char *argv[]) {
    for (int32 i = 1; i < argc; ++i) {
        char *p = argv[i];
        if (strstr(p, "-q=")) {
            minsize = atoi(p+3);
        }
        else if (strstr(p, "-k=")) {
            k = atoi(p+3);
        }
        else if (strstr(p, "-t=")) {
            threads = atoi(p+3);
        }
    }
    printf("k=%d, minsize=%d, threads=%d\n", k, minsize, threads);
    assert(minsize >= 2*k-1);
}

void Parallel_Enum::initIndexAndFwd(vector<int> &nodeset)
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

void Parallel_Enum::initRoot(Graph &subG, Vertxs &sVtx, int32 v, vector<int32> &R, Node &P, vector<boolean> &visited)
{
    P.clear(); visited[v] = true;
    vector<vector<int32>> temp(2);
    int32 p = 0, q = 1;
    int32 vpid = VtxPos[v];
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
    subG.excluded.clear();
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
               subG.excluded.emplace_back(u);
        }
    }
    int32 exlen = subG.excluded.size();
    for (int32 i = 0; i < len; ++i) {
        int32 u = P.nodes[i].first;
        int32 du = deg[u];
        for (int32 j = 0; j < du; ++j) {
            int32 w = adj[u][j];
            if (core[w] + k < minsize) continue;
            if (!visited[w] && VtxPos[w]>vpid) {
                visited[w] = true;
                P.cinsert(w, 1);
            }
            if (!visited[w] && VtxPos[w]<vpid) {
                subG.excluded.emplace_back(w);
                visited[w] = true;
            }
        }
    }
    for (int32 i = exlen; i < subG.excluded.size(); ++i) {
        int32 u = subG.excluded[i];
        int32 d = 0;
        visited[u] = false;
        for (int32 j = 0; j < len; ++j) {
            int32 w = P.nodes[j].first;
            if (is_nbr(u,w)) d++;
        }
        if (d + 2*k >= minsize+2)
            subG.excluded[exlen++] = u;
    }
    subG.excluded.resize(exlen);

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
            P.nodes[len++].first = u;
        }
    }
    P.size = len;
    for (int32 i = 0; i < len1; ++i) {
        int32 u = P.nodes[i].first;
        visited[u] = false;
    }

    visited[v] = false;
    R[0] = v; sVtx.counters[0] = 1;
    for (int32 i = 0; i < deg[v]; ++i) visited[adj[v][i]] = false;
}

void Parallel_Enum::setSubadj(Graph &subG, Vertxs &sVtx,vector<int32> &R, int32 rsize, Node &P)
{
    int32 v = R[0];
    int32 vid = P.size;
    R[0] = vid;
    subG.adj[vid].clear();
    sVtx.maps[vid] = v;
    sVtx.subD[vid] = 0;
    subG.nodes = P.size+1;
    for (int32 i = 0; i < P.size; ++i) {
        int32 u = P.nodes[i].first;
        P.nodes[i].first = i;
        subG.adj[i].clear();
        sVtx.maps[i] = u;
        sVtx.subD[i] = 0;
        sVtx.noadj[k*i] = 0;
    }
    for (int32 i = P.exsize; i < P.size; ++i) {
        int32 u = sVtx.maps[i]; 
        if (P.nodes[i].second == 0) {
            subG.adj[i].emplace_back(vid);
            subG.adj[vid].emplace_back(i);
            sVtx.subD[i]++; sVtx.subD[vid]++;
        }
        for (int32 j = i+1; j < P.size; ++j) {
            int32 w = sVtx.maps[j];
            if (is_nbr(u,w)) {
                subG.adj[i].emplace_back(j);
                subG.adj[j].emplace_back(i);
                sVtx.subD[i]++; sVtx.subD[j]++;
            }
        }
    }
}

boolean Parallel_Enum::isMaximal(Graph &subG, Vertxs &sVtx, vector<int32> &R, Node &P, int32 rsize) {
    for (int32 i = P.exsize; i < P.size; ++i)
        R[rsize++] = P.nodes[i].first;
    //printf("rsize=%d\n", rsize);
    for (int32 i = 0; i < subG.excluded.size(); ++i) {
        int32 u = subG.excluded[i];
        int32 cnt = 0;
        boolean ok = true;
        for (int32 j = 0; j < rsize; ++j) {
            int32 wid = R[j];
            int32 w = sVtx.maps[wid];
            if (!is_nbr(u, w))
                if (++cnt >= k || sVtx.subD[wid]+k <= rsize) {
                    ok = false; break;
                }
        }
        if (ok) return false;
    }
    for (int32 i = 0; i < P.exsize; ++i) {
        int32 uid = P.nodes[i].first;
        int32 cnt = 0;
        boolean ok = true;
        for (int32 j = 0; j < rsize; ++j) {
            int32 wid = R[j];
            if (!is_nbr(sVtx.maps[uid], sVtx.maps[wid]))
                if (++cnt >= k || sVtx.subD[wid]+k <= rsize) {
                    ok = false; break;
                }
        }
        if (ok) return false;
    }
    return true;
}

void Parallel_Enum::addToPlex(Vertxs &sVtx, Node &P, int32 &id, vector<int32> &R, int32 rsize)
{
    assert(id < P.size);
    assert(id >= P.exsize);
    int32 vid = P.nodes[id].first;
    int32 c = P.nodes[id].second;
    R[rsize] = vid;
    for (int32 i = vid*k; i < vid*k + c; ++i) {
        int32 ucid = sVtx.noadj[i];
        sVtx.counters[ucid] += 1;
    }
    sVtx.counters[rsize] = c+1;
    sVtx.totalcnt += (k-2*c-1);
}

void Parallel_Enum::resetCounts(Vertxs &sVtx, Node &P, int32 &id)
{
    int32 vid = P.nodes[id].first;
    int32 c = P.nodes[id].second;
    for (int32 i = vid*k; i < vid*k + c; ++i) {
        int32 ucid = sVtx.noadj[i];
        sVtx.counters[ucid] -= 1;
    }
    sVtx.totalcnt -= (k-2*c-1);
}

boolean Parallel_Enum::updateSet(Graph &subG, Vertxs &sVtx, vector<int32> &R, int32 rsize, vector<int32> &_counters, Node &P, int32 id, Node &_P, vector<int32> &delnodes)
{
    int32 rl  = rsize-1;
    int32 crl = sVtx.counters[rl];
    int32 uid = R[rl];
    int32 cnt = 0;
    int32 len1 = 0;
    int32 u = sVtx.maps[uid];
    for (int32 i = 0; i < k; ++i) subG.nonbrcnts[i] = 0;
    for (int32 i = 0; i < P.exsize; ++i) {
        int32 vid = P.nodes[i].first;
        int32 c = P.nodes[i].second;
        int32 v = sVtx.maps[vid];
        int32 cn = c;
        boolean ok = true;
        
        if (!is_nbr(u,v)) {
            if (crl+1 > k || c+1 >= k) continue;
            sVtx.noadj[vid*k+cn] = rl;
            cn += 1;
        }

        for (int32 j = vid*k; j < vid*k + c; ++j) {
            int32 wcid = sVtx.noadj[j];
            if (sVtx.counters[wcid]+1 > k) {ok = false; break;}
        }
        if (ok)  {
            _P.nodes[len1].first = vid;
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
        int32 vid = P.nodes[i].first;
        int32 c = P.nodes[i].second;
        int32 v = sVtx.maps[vid];
        int32 cn = c;
        boolean ok = true;
        if (vid == uid) continue;
        if (!is_nbr(v,u)) {
            if (crl+1 > k || c+1 >= k)  {
                delnodes.emplace_back(vid);
                continue;
            }
            sVtx.noadj[vid*k+cn] = rl;
            cn += 1;
        }
        for (int32 j = vid*k; j < vid*k + c; ++j) {
            if (sVtx.counters[sVtx.noadj[j]]+1 > k) {
                ok = false; break;
            }
        }
        if (ok) {
            _P.nodes[len2].first = vid;
            _P.nodes[len2++].second = cn;

            // prunning
            if (cn == c) {
                cnt1++; 
                subG.nonbrcnts[c]++;
                if (c == 0) {cnt++;}
                else {
                    int32 maxcounts = -1;
                    int32 maxwid = 0;
                    for (int32 j = vid*k; j < vid*k + c; ++j) {
                        int32 wcid = sVtx.noadj[j];
                        int32 ct = _counters[wcid];
                        if (ct > maxcounts) {
                            maxwid = wcid; 
                            maxcounts = ct;
                        }
                    }
                    if (maxcounts < k) {
                        _counters[maxwid]++;
                        cnt++;
                    }
                }
            }
        } else delnodes.emplace_back(vid); 
    }
    _P.exsize = len1;
    _P.size = len2;

    if (rsize+k+cnt < minsize+crl) return false;

    return true;
}

void Parallel_Enum::recursion_seq(Graph &subG, Vertxs &sVtx,vector<int32> &R, int32 rsize, Node &P)
{
    if (rsize + P.size < minsize + P.exsize) return;
    cntT++;
    if (P.size == 0 && rsize >= 2 * k - 1) {
        if (isMaximal(subG,sVtx,R,P,rsize)) {
            plexcntT++;
        }
        return;
    }
    if (P.size - P.exsize == 0) return;

    int32 mini = -1;
    int32 mind = n;
    for (int32 i = 0; i < rsize; ++i) {
        int32 uid = R[i];
        if (mind == n || sVtx.subD[uid] < mind) {
            mind = sVtx.subD[uid];
            mini = i;
        }
    }
    if (mind+k < minsize) return;
    if (mind+k+P.exsize < rsize+P.size) {
        int32 vid = R[mini];
        int32 pi = -1;
        int32 cnt = 0;
        int32 v = sVtx.maps[vid];
        mind = n;
        for (int32 i = P.exsize; i < P.size; ++i) {
            int32 uid = P.nodes[i].first;
            int32 c = P.nodes[i].second;
            int32 ud = sVtx.subD[uid];
            if (!is_nbr(sVtx.maps[uid], v)) {
                if (mind == ud && cnt <= c) {
                    pi = i; cnt = c;
                }
                else if (mind > ud) {
                    mind = ud;
                    pi = i; cnt = c;
                }
            }
        }
        Node _P; _P.resize(P.size);
        vector<int32> test, _counters; 
        test.reserve(P.size - P.exsize);
        _counters.reserve(rsize);
        if (pi < P.exsize) {
            printf("pi=%d\n", pi);
            printf("v=%d, rsize=%d, P:",vid, rsize);
            for (int i = 0; i < P.exsize; ++i) printf(" [%d, %d]",P.nodes[i].first,P.nodes[i].second);
            printf("\t");
            for (int i = P.exsize; i < P.size; ++i) printf(" [%d, %d]",P.nodes[i].first,P.nodes[i].second);
            printf("\n");
        }
        assert(pi >= P.exsize);
        addToPlex(sVtx, P, pi, R, rsize);
        _counters.clear();
        for (int32 j = 0; j < rsize; ++j) 
            _counters.emplace_back(sVtx.counters[j]);
        if (updateSet(subG,sVtx,R,rsize+1, _counters, P, P.exsize, _P, test)){
            for (auto uid : test) delNodes(subG.adj[uid], sVtx.subD);
            recursion(subG,sVtx,R,rsize+1,_P);
            for (auto uid : test) addNodes(subG.adj[uid], sVtx.subD);
        }
        resetCounts(sVtx, P, pi);

        int32 uid = P.nodes[pi].first; 
        for (int32 j = 0; j < P.size; ++j)
            _P.nodes[j] = P.nodes[j];
        _P.exsize = P.exsize+1;
        _P.size = P.size;
        exchanges(_P, P.exsize, pi);
        delNodes(subG.adj[uid],sVtx.subD);
        recursion(subG,sVtx,R, rsize, _P);
        addNodes(subG.adj[uid],sVtx.subD);
        return;
    }

    int32 minicnt = 0;
    for (int32 i = P.exsize; i < P.size; ++i) {
        int32 uid = P.nodes[i].first;
        int32 c = P.nodes[i].second;
        if (mind == n || sVtx.subD[uid] < mind) {
            mind = sVtx.subD[uid];
            mini = i + rsize;
            minicnt = c;
        }
        else if (sVtx.subD[uid] == mind && minicnt < c) {
            mini = i + rsize;
            minicnt = c;
        }
    }

    //printf("deg[%d]=%d, mind=%d\n", P.nodes[mini-rsize].first, sVtx.subD[P.nodes[mini-rsize].first], mind);
    if (mind+k+P.exsize >= rsize+P.size) {
        if (rsize+P.size-P.exsize >= 2*k-1) {
            if (isMaximal(subG,sVtx, R, P, rsize))
            {
                plexcntT++;
            }
        }
        return;
    }

    Node _P; _P.resize(P.size);
    int32 i = mini-rsize;
    if (mind+k >= minsize) {
        vector<int32> test, _counters; 
        test.reserve(P.size - P.exsize);
        _counters.reserve(rsize+1);
        assert(i >= P.exsize);
        addToPlex(sVtx,P, i, R, rsize);
        _counters.clear();
        for (int32 j = 0; j < rsize; ++j) 
            _counters.emplace_back(sVtx.counters[j]);
        if (updateSet(subG,sVtx,R, rsize+1, _counters, P, P.exsize, _P, test)){
           for (auto uid : test) delNodes(subG.adj[uid], sVtx.subD);
           recursion(subG,sVtx,R, rsize+1, _P);
           for (auto uid : test) addNodes(subG.adj[uid], sVtx.subD);
        }
        resetCounts(sVtx,P, i);
    }

    int32 uid = P.nodes[i].first;
    for (int32 j = 0; j < P.size; ++j)
        _P.nodes[j] = P.nodes[j];
    _P.exsize = P.exsize+1;
    _P.size = P.size;
    exchanges(_P, P.exsize, i);
    delNodes(subG.adj[uid], sVtx.subD);
    recursion(subG,sVtx,R, rsize, _P);
    addNodes(subG.adj[uid], sVtx.subD);
}

void Parallel_Enum::recursionCalls(Graph &subG, Vertxs &sVtx,vector<int32> &R, int32 rsize, Node &P, int pi)
{
    assert(pi >= P.exsize);
    if (sVtx.subD[P.nodes[pi].first]+k >= minsize) {
        Node _P; _P.resize(P.size);
        vector<int32> test, _counters; 
        test.reserve(P.size - P.exsize);
        _counters.reserve(rsize);
        addToPlex(sVtx, P, pi, R, rsize);
        _counters.clear();
        for (int32 j = 0; j < rsize; ++j) 
            _counters.emplace_back(sVtx.counters[j]);
        if (updateSet(subG,sVtx,R,rsize+1, _counters, P, P.exsize, _P, test)){
            if (threads>=2 && division && (tasks-overtasks)<threads && (rsize+_P.size)>minsize+_P.exsize+ParTASK) {
                Vertxs *_sVtx = new Vertxs[1];
                int subsize = subG.nodes;
                int minc = min(maxcore+k,rsize+P.size-P.exsize);
                (*_sVtx).maps.assign(sVtx.maps.begin(),sVtx.maps.begin()+subsize);
                (*_sVtx).subD.assign(sVtx.subD.begin(),sVtx.subD.begin()+subsize);
                (*_sVtx).noadj.assign(sVtx.noadj.begin(),sVtx.noadj.begin()+subsize*k);
                (*_sVtx).counters.assign(sVtx.counters.begin(),sVtx.counters.begin()+minc);
                (*_sVtx).totalcnt = sVtx.totalcnt;
                #pragma omp atomic
                tasks++;
                for (auto uid : test) delNodes(subG.adj[uid], (*_sVtx).subD);
                #pragma omp task shared(subG,_sVtx)
                {
                    (*_sVtx).kplexnums = 0; (*_sVtx).iterations = 0;
                    recursion(subG,(*_sVtx),R,rsize+1,_P);
                    #pragma omp atomic update
                    overtasks++;

                    int tid = omp_get_thread_num();
                    lcnt_plx[tid] += (*_sVtx).kplexnums;
                    lcnt_its[tid] += (*_sVtx).iterations;

                    delete[] _sVtx;
                }
            }
            else
            {
                for (auto uid : test) delNodes(subG.adj[uid], sVtx.subD);
                recursion(subG,sVtx,R,rsize+1,_P);
                for (auto uid : test) addNodes(subG.adj[uid], sVtx.subD);

            }
        }
        resetCounts(sVtx, P, pi);

    }
    Node _P; _P.resize(P.size);
    int32 uid = P.nodes[pi].first; 
    for (int32 j = 0; j < P.size; ++j)
        _P.nodes[j] = P.nodes[j];
    _P.exsize = P.exsize+1;
    _P.size = P.size;
    exchanges(_P, P.exsize, pi);
    delNodes(subG.adj[uid],sVtx.subD);
    recursion(subG,sVtx,R, rsize, _P);
    addNodes(subG.adj[uid],sVtx.subD);
}

void Parallel_Enum::recursion(Graph &subG, Vertxs &sVtx,vector<int32> &R, int32 rsize, Node &P)
{
    if (rsize + P.size < minsize + P.exsize) return;
    sVtx.iterations++;
    if (P.size == 0 && rsize >= minsize) {
        if (isMaximal(subG,sVtx,R,P,rsize)) {
            sVtx.kplexnums++;
        }
        return;
    }
    if (P.size - P.exsize == 0) return;

    int32 mini = -1;
    int32 mind = n;
    for (int32 i = 0; i < rsize; ++i) {
        int32 uid = R[i];
        if (mind == n || sVtx.subD[uid] < mind) {
            mind = sVtx.subD[uid];
            mini = i;
        }
    }
    if (mind+k < minsize) return;
    if (mind+k+P.exsize < rsize+P.size) {
        int32 vid = R[mini];
        int32 pi = -1;
        int32 cnt = 0;
        int32 v = sVtx.maps[vid];
        mind = n;
        for (int32 i = P.exsize; i < P.size; ++i) {
            int32 uid = P.nodes[i].first;
            int32 c = P.nodes[i].second;
            int32 ud = sVtx.subD[uid];
            if (!is_nbr(sVtx.maps[uid], v)) {
                if (mind == ud && cnt <= c) {
                    pi = i; cnt = c;
                }
                else if (mind > ud) {
                    mind = ud;
                    pi = i; cnt = c;
                }
            }
        }
        recursionCalls(subG,sVtx,R,rsize,P,pi);
        return;
    }

    int32 minicnt = 0;
    for (int32 i = P.exsize; i < P.size; ++i) {
        int32 uid = P.nodes[i].first;
        int32 c = P.nodes[i].second;
        if (mind == n || sVtx.subD[uid] < mind) {
            mind = sVtx.subD[uid];
            mini = i + rsize;
            minicnt = c;
        }
        else if (sVtx.subD[uid] == mind && minicnt < c) {
            mini = i + rsize;
            minicnt = c;
        }
    }

    if (mind+k+P.exsize >= rsize+P.size) {
        if (rsize+P.size-P.exsize >= minsize) {
            if (isMaximal(subG,sVtx, R, P, rsize)) {
                sVtx.kplexnums++;
            }
        }
        return;
    }
    recursionCalls(subG,sVtx,R,rsize,P,mini-rsize);
}


void Parallel_Enum::penum()
{
    double dtm = omp_get_wtime();
    omp_set_num_threads(threads);
    vector<int32> nodeset(n);
    #pragma omp parallel for
    for (int32 i = 0; i < n; ++i) nodeset[i] = i;
    maxcore = core_decompsition(nodeset, n);
    initIndexAndFwd(nodeset);
    //printf("Maxcore Number:%d\nDecomp time:%f sec\n", maxcore, (omp_get_wtime()-dtm));
    printf("Decomp time:%f sec\n", (omp_get_wtime()-dtm));

    VtxPos.resize(n);
    double tm = omp_get_wtime();
    int computed = 0;
    long kplexnums = 0;
    long iterations = 0;
    lcnt_plx.resize(threads,0);
    lcnt_its.resize(threads,0);
    int sn =  min(n, maxcore*md+maxcore);
    if (2*k < minsize) sn = min(n, maxcore*md/(minsize-2*k)+maxcore);
    #pragma omp parallel reduction(+:computed,iterations,kplexnums)
    {
        //printf("Test: thread_num:%d\n", omp_get_thread_num());
        Graph subG;
        Vertxs sVtx;
        Node P;
        vector<int> R;
        vector<bool> vis;
        R.resize(maxcore+k);
        vis.resize(n, 0);
        P.reserve(md);
        subG.adj.resize(sn);
        subG.excluded.reserve(md);
        subG.nonbrcnts.resize(maxcore+k);
        sVtx.counters.resize(maxcore+k);
        sVtx.noadj.resize(k*sn);
        sVtx.maps.resize(sn);
        sVtx.subD.resize(sn);
        double wtime = omp_get_wtime();
        
        #pragma omp for
        for (int i = 0; i < n; ++i) 
            VtxPos[nodeset[i]] = i;
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            int v = nodeset[i];
            if (core[v] + k < minsize) continue;
            initRoot(subG,sVtx,v,R,P,vis);
            if (P.size+1 < minsize) continue;
            if (core[v]>=maxcore && !division) division = true;
            setSubadj(subG,sVtx,R,1,P);

            sVtx.totalcnt = k-1;
            sVtx.iterations = 0;
            sVtx.kplexnums = 0;
            #pragma omp taskgroup
            recursion(subG,sVtx,R,1,P);

            computed++;
            iterations += sVtx.iterations;
            kplexnums += sVtx.kplexnums;
        }
        wtime = omp_get_wtime() - wtime;
    }
    for (int i = 0; i < threads; ++i) {
        kplexnums+=lcnt_plx[i];
        iterations+=lcnt_its[i];
    }
    //printf("computed : %d \n", computed);
    printf("Number of plexes    : %ld \n", kplexnums);
    printf("Number of iterations: %ld \n", iterations);
    //printf("Tasks=%d, Overtasks=%d \n", tasks, overtasks);
    printf("Enum time: %f sec\n", (omp_get_wtime() - tm));
}

void Parallel_Enum::run()
{
    //double tm = omp_get_wtime();
    penum();
    //printf("All time:%f sec\n", (omp_get_wtime()-tm));
}
