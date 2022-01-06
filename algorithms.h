#pragma once
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <queue>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cassert>
#include <cmath>
#include <utility>
#include <immintrin.h>

using namespace std;

typedef int int32;
typedef char char8;
typedef long long64;
typedef bool boolean;
typedef unsigned int uint32;
typedef pair<int32, int32> upair;

class Algorithm
{
protected:
	/* data */
	int32 n, m, md;
    int32 *deg, *datas, **adj;
	int32 mincliquesize;
    int32 algorithm;
    vector<int32> core, topcore, colors;
public:
    Algorithm(/* args */);
    virtual ~Algorithm();
    virtual void run() {}
    virtual void setParameters(int32 argc, char **argv) {}

    void read_graph(const char *str);
    void testprintGraph();
    int32 core_decompsition(int32 *nodeset, int32 nodesize);
};

#define unfilled -1
class CuckooHash
{
private:
	/* data */
	int32 capacity;
	int32 mask;
	int32 size;
	int32 buff_size = sizeof(int32);
	int32 *hashtable;

	void rehash(int32 **_table) {
		int32 oldcapacity = capacity;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		capacity = (mask + 1) * buff_size;
		int32 *newhash = new int32[capacity];
		memset((newhash), unfilled, sizeof(int32) * capacity);
		for (int32 i = 0; i < oldcapacity; ++i){
			if ((*_table)[i] != unfilled) insert((*_table)[i], &newhash);
		}
		swap((*_table), newhash);
		delete[] newhash;
	}
	void insert(const int32 &_u, int32 **_table) {
		
		int32 hs = hash1(_u);
		for (int32 i = 0; i < buff_size; ++i) {
			if ((*_table)[hs * buff_size + i] == unfilled){
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}
		hs = hash2(_u);
		for (int32 i = 0; i < buff_size; ++i) {
			if ((*_table)[hs * buff_size + i] == unfilled){
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}

		boolean use_hash1 = true;
		int32 u = _u;
		for (int32 i = 0; i < mask; ++i) {
			int32 replaced;
			if (use_hash1) hs = hash1(u);
			else hs = hash2(u);
			int32 j = 0;
			for (; j < buff_size; ++j) {
				if ((*_table)[hs * buff_size + j] == unfilled) break;
			}
			if (buff_size == j) {
				replaced = move((*_table)[hs * buff_size]);
				j = 1;
				for (; j < buff_size; j++) {
					(*_table)[hs * buff_size + j - 1] =
						move((*_table)[hs * buff_size + j]);
				}
				(*_table)[hs * buff_size + j - 1] = u;
			}
			else {
				replaced = move((*_table)[hs * buff_size + j]);
				(*_table)[hs * buff_size + j] = u;
			}
			use_hash1 = hs == hash2(replaced);
			u = move(replaced);
			if (u == unfilled) return;
		}
		rehash(_table);
		insert(u, _table);
	}

	int32 hash1(const int32 &x) { return x & mask;}
	int32 hash2(const int32 &x) { return ~x & mask;}

public:
	CuckooHash(/* args */) {
		capacity = 0;
		hashtable = NULL;
		mask = 0;
		size = 0;
	}
	~CuckooHash() {
		if (hashtable) delete[] hashtable;
	}

	void reserve(int32 _size) {
		if (capacity >= _size) return;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		while (_size >= mask * buff_size) mask = (mask << 1) | 1;
		capacity = (mask + 1) * buff_size;
		if (hashtable) delete[] hashtable;
		hashtable = new int32[capacity];
		memset(hashtable, unfilled, sizeof(int32) * capacity);
	}

	void insert(const int32 &_u) {
		if (find(_u)) return;
		insert(_u, &hashtable);
		size++;
	}

	boolean find(const int32 &_u) {
		int32 hs1 = hash1(_u);
		int32 hs2 = hash2(_u);

		assert(buff_size == 4 && sizeof (int32) == 4);
		__m128i cmp = _mm_set1_epi32(_u);
		__m128i b1 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs1]);
		__m128i b2 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs2]);
		__m128i flag = _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));

		return _mm_movemask_epi8(flag) != 0;
	}
	int32 getcapacity() {return capacity;}
	int32 getmask() {return mask;}
	int32 *gethashtable() {return hashtable;}
};

typedef struct _Node
{
    int32 size, exsize; 
    vector<upair> nodes;

	_Node() {
		size 	= 0;
		exsize  = 0;
	}
	_Node(int32 n) {
		nodes.resize(n);
		size 	= n;
		exsize  = 0;
	}
	void resize(int32 n) {
		nodes.resize(n);
		size 	= n;
		exsize  = 0;
	}
	void reserve(int32 n) {
		nodes.reserve(n);
		size 	= 0;
		exsize  = 0;
	}
	void cinsert(const int32 &x, const int32 &c) {
		nodes.emplace_back(x,c);
		size++;
	}
	void einsert(const int32 &x, const int32 &c) {
		nodes.emplace_back(x, c);
		size++;
		exsize++;
	}
	void cinsert(const int32 &id, const int32 &x, const int32 &c) {
		nodes[id].first  = x;
		nodes[id].second = c;
		size++;
	}
	void einsert(const int32 &id, const int32 &x, const int32 &c) {
		nodes[id].first  = x;
		nodes[id].second = c;
		size++;
		exsize++;
	}
	void clear() {
		size	= 0;
		exsize  = 0;
		nodes.clear();
	}
} Node;