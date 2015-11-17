/*
 * SeqDB.h
 *
 *  Created on: Sep 14, 2010
 *      Author: xiaoliwang
 */

#ifndef _SEQDB_H_
#define _SEQDB_H_

#include "Time.h"
#include "Gram.h"
#include "GramList.h"
#include "CountFilter.h"
#include "Query.h"

#include <vector>
#include <set>
#include <algorithm>
#include <stdint.h>
#include <string>
#include <fstream>
#include <queue>

using namespace std;

#ifdef _WIN32
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

class Candi {
public:
  unsigned id;
  unsigned count;

  Candi(unsigned id, unsigned count) : id(id), count(count) {}
  ~Candi(){}
};

struct queue_entry{
	unsigned m_sid;
	unsigned m_dist;

	queue_entry(unsigned sid, unsigned dist){ m_sid = sid; m_dist = dist; };

	bool operator< (const queue_entry &other) const
	{
		return (m_dist<other.m_dist);
	};
};

struct min_queue_entry{
	unsigned m_sid;
	unsigned m_dist;

	min_queue_entry(unsigned sid, unsigned dist){ m_sid = sid; m_dist = dist; };

	bool operator< (const min_queue_entry &other) const
	{
		return (m_dist>other.m_dist);
	};
};

template <class InvList>
class CSeqDB;

template <class InvList = vector<unsigned> >
class CSeqDB
{
protected:
#ifdef _WIN32
  typedef typename unordered_map<unsigned, CGramList<InvList>* > GramList;
  typedef typename unordered_map<unsigned, string> GramIndex;
  typedef typename unordered_map<unsigned, unsigned> FQueue;
#else 
  typedef typename tr1::unordered_map<unsigned, CGramList<InvList>* > GramList;
  typedef typename tr1::unordered_map<unsigned, string> GramIndex;
  typedef typename tr1::unordered_map<unsigned, unsigned> FQueue;
#endif

	vector<string>* data;
	vector<queue_entry> datasizes;
	int* datasize_groupindex;

	GramList gramListUpper;
	vector<unsigned> gramCodeUpper;
	vector<string> gramStrUpper;
	GramList gramListLower;

	bool* processedData;
	unsigned** dataCount;
	unsigned fk;
	int** edcost;

	vector<int> post_cand_low;
	vector<int> post_cand_up;
	
	// For statistics
public:
	priority_queue<queue_entry> m_queue;
	CQuery* theQuery;
	CCountFilter *filter;
	CGram* gramGenUpper;
	CGram* gramGenLower;
	unsigned gramGenUpper_len;
	unsigned gramGenLower_len;
	bool stop;
	unsigned max_listlen;
	unsigned min_listlen;
	unsigned processed;
	unsigned data_maxlen;
	unsigned gram_maxed;

public:
	CSeqDB(const unsigned& maxlen, const unsigned max_gramed, vector<string>* sequences = NULL, CGram* gramgenupper = NULL, CGram* gramgenlower = NULL, CCountFilter *fltable = NULL) : 
	  data(sequences), gramGenUpper(gramgenupper), gramGenLower(gramgenlower), filter(fltable)
	{
		this->datasize_groupindex = new int[maxlen+1];
		for (unsigned i = 0; i < sequences->size(); i++) {
			datasizes.push_back(queue_entry(i, sequences->at(i).length()));
			this->datasize_groupindex[sequences->at(i).length()] = i;
		}
		this->processedData = new bool[sequences->size()];
		this->gram_maxed = max_gramed;
		this->dataCount = NULL;
		this->data_maxlen = maxlen+1;
		this->edcost = new int *[this->data_maxlen];
		for (unsigned i = 0; i < this->data_maxlen; i++) 
			this->edcost[i] = new int[this->data_maxlen];
	}
	~CSeqDB(void)
	{
		for (unsigned i = 0; i < this->data_maxlen; i++) 
		{
			if (this->edcost[i] != NULL)
			{
				delete [] this->edcost[i];
				this->edcost[i] = NULL;
			}
		}
		delete [] this->edcost;
		this->edcost = NULL;

		if (this->dataCount != NULL) {
			for (unsigned i = 0; i < this->gram_maxed; i++) {
				if (this->dataCount[i] != NULL)
				{
					delete [] this->dataCount[i];
					this->dataCount[i] = NULL;
				}
			}
			delete [] this->dataCount;
			this->dataCount = NULL;
		}

		if (this->processedData != NULL) {
			delete [] this->processedData;
			this->processedData = NULL;
		}

		if (this->datasize_groupindex != NULL) {
			delete [] this->datasize_groupindex;
			this->datasize_groupindex = NULL;
		}

		this->data = NULL;
		this->gramGenUpper = NULL;
		this->gramGenLower = NULL;
		this->filter = NULL;
	}

	void initParas(const unsigned& kf, CGram* gramgenupper, CGram* gramgenlower, CCountFilter *fltable)
	{
		this->fk = kf;
		this->gramGenUpper = gramgenupper;
		this->gramGenLower = gramgenlower;
		this->filter = fltable;
	}

	void reset()
	{
		for (unsigned i = 0; i < this->data->size(); i++) {
			this->processedData[i] = false;
		}
		if(!this->m_queue.empty()) {
			priority_queue<queue_entry> tmp;
			swap(this->m_queue, tmp);
		}
		if(this->post_cand_low.size() > 0) {
			vector<int> tmp;
			swap(this->post_cand_low, tmp);
		}
		if(this->post_cand_up.size() > 0) {
			vector<int> tmp;
			swap(this->post_cand_up, tmp);
		}
		this->processed = 0;
		this->stop = false;
	}

	void buildindex();
	void insertSequenceIntoIndex(const string& current, const unsigned sid, GramIndex& gramIndex);
	void insertNgramIntoIndex(const string& current, const unsigned gid);
	void print(FILE* pf);
	void save(const char* indexName);
	void saveIndex(ofstream& pf, GramList& theGramList);
	void saveIndexGramCode(ofstream& pf);
	void saveIndexGramStr(ofstream& pf);
	void load(const char* indexName);
	void loadIndexGramCode(ifstream& pf);
	void loadIndexGramStr(ifstream& pf);
	void loadIndex(ifstream& pf, GramList& theGramList);

	bool getEditDistance(const string &tString, const string &qString, int constraint);
	int getRealEditDistance_ns(const string &tString, const string &qString, int constraint);
	int getRealEditDistance_nsd(const string &tString, const string &qString);

	void getGramLists(const vector<unsigned>& gramCodes, GramList& theGramList, vector<InvList*> &invLists);
	static bool cmpInvList(const InvList* a, const InvList* b) 
	{
		return a->size() < b->size();
	}
	void expProbe(const typename InvList::iterator start, 
		const typename InvList::iterator end, 
		typename InvList::iterator& lbound, 
		typename InvList::iterator& ubound, 
		unsigned value) const;
	void merge(vector<InvList*>& invLists, const unsigned threshold, vector<Candi>& candiscur);
	void gramQuery(const CQuery& query, const unsigned ged, vector< vector<unsigned> >& similarGrams);
	void getIDBounds(const int querysize, const int threshold, int& idlow, int& idhigh);

	// The fuction for pipe knn search
	void init_threshold(CQuery& query);
	void knn_postprocess();
	void accumulateFrequency(long ged);

	// The functions for round robin processing
	

	// Interfaces for GPU processing testing
	void queryGramLists(const CQuery& query);
};

template <class InvList>
void CSeqDB<InvList>::buildindex() 
{
	GramIndex upperGrams;
	unsigned gramCode;
	for(unsigned sid = 0; sid < this->data->size(); sid++) {
		string current = this->data->at(sid);
		insertSequenceIntoIndex(current, sid, upperGrams);
	}

	typename GramIndex::iterator iter;
	unsigned i = 0;
	for(iter = upperGrams.begin(); iter != upperGrams.end(); iter++, i++)
	{
		this->gramCodeUpper.push_back(iter->first);
		this->gramStrUpper.push_back(iter->second);
		insertNgramIntoIndex(iter->second, i);
	}
}

template <class InvList>
void CSeqDB<InvList>::insertSequenceIntoIndex(const string& current, const unsigned sid, GramIndex& uppergrams) 
{
	// add sid to gramListUpper
	vector<string> grams;
	vector<unsigned> gramCodes;
	this->gramGenUpper->decompose(current, grams, gramCodes);
	for(unsigned i = 0; i < gramCodes.size(); i++) {
		unsigned gramCode = gramCodes[i];

		if (this->gramListUpper.find(gramCode) == this->gramListUpper.end()) {
			// a new gram
			CGramList<InvList>* newGramList = new CGramList<InvList>();
			this->gramListUpper[gramCode] = newGramList;
			newGramList->getArray()->push_back(sid);
		}
		else {
			// an existing gram
			CGramList<InvList>* theGramList = this->gramListUpper[gramCode];
			// avoid adding duplicate sequences
			if(theGramList->getArray()->back() != sid)
				theGramList->getArray()->push_back(sid);
		}

		if (uppergrams.find(gramCode) == uppergrams.end()) {
			uppergrams[gramCode] = grams[i];
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::insertNgramIntoIndex(const string& current, const unsigned gid)
{
	// add sid to gramListUpper
	vector<unsigned> gramCodes;
	this->gramGenLower->decompose(current, gramCodes);
	for(unsigned i = 0; i < gramCodes.size(); i++) {
		unsigned gramCode = gramCodes[i];

		if (this->gramListLower.find(gramCode) == this->gramListLower.end()) {
			// a new gram
			CGramList<InvList>* newGramList = new CGramList<InvList>();
			this->gramListLower[gramCode] = newGramList;
			newGramList->getArray()->push_back(gid);
		}
		else {
			// an existing gram
			CGramList<InvList>* theGramList = this->gramListLower[gramCode];
			// avoid adding duplicate sequences
			if(theGramList->getArray()->back() != gid)
				theGramList->getArray()->push_back(gid);
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::print(FILE* pf) 
{
	fprintf(pf, "The upper-level index:\n");
	typename GramList::iterator iter;
	for(iter = gramListUpper.begin(); iter != gramListUpper.end(); iter++) {
		fprintf(pf, "%u:\t", iter->first);
		InvList* theList = iter->second->getArray();
		for(unsigned i = 0; i < theList->size(); i++) {
			fprintf(pf, "%u ", theList->at(i));
		}
		fprintf(pf, "\n");
	}
	fprintf(pf, "The upper-level gram codes:\n");
	for(unsigned i = 0; i<this->gramCodeUpper.size(); i++) {
		fprintf(pf, "%u\n", this->gramCodeUpper[i]);
	}
	fprintf(pf, "The upper-level grams:\n");
	for(unsigned i = 0; i<this->gramStrUpper.size(); i++) {
		fprintf(pf, "%s\n", this->gramStrUpper[i].c_str());
	}
	fprintf(pf, "The lower-level index:\n");
	for(iter = gramListLower.begin(); iter != gramListLower.end(); iter++) {
		fprintf(pf, "%u:\t", iter->first);
		InvList* theList = iter->second->getArray();
		for(unsigned i = 0; i < theList->size(); i++) {
			fprintf(pf, "%u ", theList->at(i));
		}
		fprintf(pf, "\n");
	}
}

template <class InvList>
void CSeqDB<InvList>::save(const char* indexName)
{
	ofstream pf;
	pf.open(indexName, ios::out);
	if(pf.is_open()) {
		// Write the gram length for the upper and lower level
		unsigned gl_size = this->gramGenUpper->getGramLength();
		pf.write((const char*)&gl_size, sizeof(unsigned));
		gl_size = this->gramGenLower->getGramLength();
		pf.write((const char*)&gl_size, sizeof(unsigned));
		saveIndex(pf, this->gramListUpper);
		saveIndexGramCode(pf);
		saveIndexGramStr(pf);
		saveIndex(pf, this->gramListLower);
		pf.close();
	}
	else {
		fprintf(stderr, "Error: Failed to open the index file - %s\n", indexName);
	}
}

template <class InvList>
void CSeqDB<InvList>::saveIndex(ofstream& pf, GramList& theGramList)
{
	unsigned val;
	// write the gram list
	val = theGramList.size();
	pf.write((const char*)&val, sizeof(unsigned));
	for(typename GramList::iterator iter = theGramList.begin(); iter != theGramList.end(); iter++) {
		val = iter->first;
		pf.write((const char*)&val, sizeof(unsigned));
		InvList* theList = iter->second->getArray();
		if (theList) {
			val = theList->size();
			pf.write((const char*)&val, sizeof(unsigned));
			for(unsigned i = 0; i < theList->size(); i++) {
				val = theList->at(i);
				pf.write((const char*)&val, sizeof(unsigned));
			}
		}
		else {
			val = 0;
			pf.write((const char*)&val, sizeof(unsigned));
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::saveIndexGramCode(ofstream& pf)
{
	unsigned val;
	// write the gram list
	val = this->gramCodeUpper.size();
	pf.write((const char*)&val, sizeof(unsigned));
	for(unsigned i = 0; i < this->gramCodeUpper.size(); i++) {
		val = this->gramCodeUpper[i];
		pf.write((const char*)&val, sizeof(unsigned));
	}
}

template <class InvList>
void CSeqDB<InvList>::saveIndexGramStr(ofstream& pf)
{
	unsigned val;
	// write the gram list
	val = this->gramStrUpper.size();
	pf.write((const char*)&val, sizeof(unsigned));
	string gram;
	size_t gramLen = this->gramGenUpper->getGramLength();
	for(unsigned i = 0; i < this->gramStrUpper.size(); i++) {
		gram = this->gramStrUpper[i];
		char *a=new char[gram.size()+1];
		strcpy(a, gram.c_str());
		pf.write(a, gramLen);
	}
}

template <class InvList>
void CSeqDB<InvList>::load(const char* indexName)
{
	ifstream pf;
	pf.open(indexName, ios::in);
	if(pf.is_open()) {
		// Read the gram length for upper and lower level
		unsigned gl_size;
		pf.read((char*)&gl_size, sizeof(unsigned));
		this->gramGenUpper_len = gl_size;
		pf.read((char*)&gl_size, sizeof(unsigned));
		this->gramGenLower_len = gl_size;
		this->max_listlen = 0;
		this->min_listlen = this->data->size();
		loadIndex(pf, this->gramListUpper);
		const unsigned maxlen = this->max_listlen;
		const unsigned minlen = this->min_listlen;
		loadIndexGramCode(pf);
		loadIndexGramStr(pf);
		loadIndex(pf, this->gramListLower);
		this->max_listlen = maxlen;
		this->min_listlen = minlen;
		pf.close();

		// Initialize required parameters
		this->gram_maxed = this->gram_maxed < ((this->gramGenUpper_len-this->gramGenLower_len)/this->gramGenLower_len+1) ? this->gram_maxed : ((this->gramGenUpper_len-this->gramGenLower_len)/this->gramGenLower_len+1);
		this->dataCount = new unsigned *[this->gram_maxed];
		for (unsigned i = 0; i < this->gram_maxed; i++) {
			this->dataCount[i] = new unsigned[this->data->size()];
			for (unsigned j = 0; j < this->data->size(); j++)
				this->dataCount[i][j] = 0;
		}
	}
	else {
		fprintf(stderr, "Error: Failed to open the index file - %s\n", indexName);
	}
}

template <class InvList>
void CSeqDB<InvList>::loadIndex(ifstream& pf, GramList& theGramList)
{
	// read the gram list
	unsigned gramListSize;
	unsigned gramCode;
	unsigned listSize;
	pf.read((char*)&gramListSize, sizeof(unsigned));
	for(unsigned i = 0; i < gramListSize; i++) {
		pf.read((char*)&gramCode, sizeof(unsigned));
		pf.read((char*)&listSize, sizeof(unsigned));
		if (this->max_listlen < listSize)
		this->max_listlen = listSize;
		if (this->min_listlen > listSize)
		this->min_listlen = listSize;
		if (listSize > 0) {
			CGramList<InvList>* newGramList = new CGramList<InvList>();
			theGramList[gramCode] = newGramList;
			InvList* theList = newGramList->getArray();
			for(unsigned j = 0; j < listSize; j++) {
				unsigned sid;
				pf.read((char*)&sid, sizeof(unsigned));
				theList->push_back(sid);
			}
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::loadIndexGramCode(ifstream& pf)
{
	// read the gram index
	unsigned gramListSize;
	unsigned gramCode;
	pf.read((char*)&gramListSize, sizeof(unsigned));
	for(unsigned i = 0; i < gramListSize; i++) {
		pf.read((char*)&gramCode, sizeof(unsigned));
		this->gramCodeUpper.push_back(gramCode);
	}
}

template <class InvList>
void CSeqDB<InvList>::loadIndexGramStr(ifstream& pf)
{
	// read the gram index
	unsigned gramListSize;
	pf.read((char*)&gramListSize, sizeof(unsigned));
	const size_t& gramLen = this->gramGenUpper_len;
	for(unsigned i = 0; i < gramListSize; i++) {
		char *gramArr = new char[gramLen+1];
		pf.read(gramArr, gramLen);
		string gramStr = gramArr;
		this->gramStrUpper.push_back(gramStr);
	}
}

// Origin verify function from Flamingo
template <class InvList>
bool CSeqDB<InvList>::getEditDistance(const string &tString, const string &qString, int constraint)
{
	int qStrLen = (int)qString.length();
	int tStrLen = (int)tString.length();
	int i, j, min_cost;
	for(i = 0; i < qStrLen+1; i++)
		edcost[0][i] = i;
	int startPos, endPos;
	for(i = 1; i < tStrLen+1; i++) {
		edcost[i][0] = i;
		min_cost = qStrLen;
		startPos = (i - (constraint + 1)) > 1 ? (i - (constraint + 1)) : 1;
		endPos = (i + (constraint + 1)) < qStrLen ? (i + (constraint + 1)) : qStrLen;
		if(startPos > 1) {
			edcost[i][startPos-1] = edcost[i-1][startPos-1] < edcost[i-1][startPos-2] ? (edcost[i-1][startPos-1] + 1) : (edcost[i-1][startPos-2] + 1);
		}
		for(j = startPos; j < endPos+1; j++) {
			edcost[i][j] = std::min(edcost[i-1][j-1]+((tString[i-1]==qString[j-1])?0:1), std::min(edcost[i-1][j]+1, edcost[i][j-1]+1));
			min_cost = edcost[i][j] < min_cost ? edcost[i][j] : min_cost;
		}
		if(j < qStrLen+1) {
			edcost[i][j] = edcost[i][j-1]<edcost[i-1][j-1]?(edcost[i][j-1]+1):(edcost[i-1][j-1]+1);
		}

		if(min_cost > constraint)
			return false;
	}	
	return edcost[tStrLen][qStrLen] <= constraint;
}


// Origin verify function from Flamingo
template <class InvList>
int CSeqDB<InvList>::getRealEditDistance_ns(const string &tString, const string &qString, int constraint)
{
	const int qStrLen=(int)qString.length(), tStrLen=(int)tString.length();
	int i, j, min_cost;
	for(i = 0; i < qStrLen+1; i++)
		edcost[0][i] = i;
	int startPos, endPos;
	for(i = 1; i < tStrLen+1; i++) {
		edcost[i][0] = i;
		min_cost = qStrLen;
		startPos = (i - (constraint + 1)) > 1 ? (i - (constraint + 1)) : 1;
		endPos = (i + (constraint + 1)) < qStrLen ? (i + (constraint + 1)) : qStrLen;
		if(startPos > 1) {
			edcost[i][startPos-1] = edcost[i-1][startPos-1] < edcost[i-1][startPos-2] ? (edcost[i-1][startPos-1] + 1) : (edcost[i-1][startPos-2] + 1);
		}
		for(j = startPos; j < endPos+1; j++) {
			edcost[i][j] = std::min(edcost[i-1][j-1]+((tString[i-1]==qString[j-1])?0:1), std::min(edcost[i-1][j]+1, edcost[i][j-1]+1));
			min_cost = edcost[i][j] < min_cost ? edcost[i][j] : min_cost;
		}
		if(j < qStrLen+1) {
			edcost[i][j] = edcost[i][j-1]<edcost[i-1][j-1]?(edcost[i][j-1]+1):(edcost[i-1][j-1]+1);
		}
		if(min_cost > constraint+1) {
			return min_cost;
		}

	}
	return edcost[tStrLen][qStrLen];
}

template <class InvList>
int CSeqDB<InvList>::getRealEditDistance_nsd(const string &tString, const string &qString)
{
	const int qStrLen=(int)qString.length(), tStrLen=(int)tString.length();
	int i, j;
	for(i=0; i < qStrLen+1; i++) edcost[0][i] = i;
	for(i=1; i < tStrLen+1; i++) {
		edcost[i][0] = i;
		for(j = 1; j < qStrLen+1; j++)
			edcost[i][j] = std::min(edcost[i-1][j-1]+((tString[i-1]==qString[j-1])?0:1), std::min(edcost[i-1][j]+1, edcost[i][j-1]+1));
	}
	return edcost[tStrLen][qStrLen];
}

template <class InvList>
void CSeqDB<InvList>::getGramLists(const vector<unsigned>& gramCodes, GramList& theGramList, vector<InvList*> &invLists)
{
	unsigned gramCode;
	for(unsigned i = 0; i < gramCodes.size(); i++) {
		gramCode = gramCodes[i];
		if(theGramList.find(gramCode) != theGramList.end()) {
		  InvList* tmp = theGramList[gramCode]->getArray();
		  invLists.push_back(tmp);
		}
	} 
}

template <class InvList>
void CSeqDB<InvList>::expProbe(const typename InvList::iterator start, 
	 const typename InvList::iterator end, 
	 typename InvList::iterator& lbound, 
	 typename InvList::iterator& ubound, 
	 unsigned value) const {
  
	unsigned c = 0;    
	lbound = start;
	for(;;) {      
		ubound = start + (1 << c);
		if(ubound >= end) {
			ubound = end;
			return;
		}
		else if(*ubound >= value) return;
		c++;
		lbound = ubound;
	}
}

template <class InvList>
void CSeqDB<InvList>::merge(vector<InvList*>& invLists, const unsigned threshold, vector<Candi>& candiscur)
{
	sort(invLists.begin(), invLists.end(), CSeqDB<InvList>::cmpInvList);
  
	unsigned numShortLists = invLists.size() - threshold + 1;
  
	// process the short lists using the algorithm from
	// Naoaki Okazaki, Jun'ichi Tsujii
	// "Simple and Efficient Algorithm for Approximate Dictionary Matching", COLING 2010

	vector<Candi> candis;
	for(unsigned i = 0; i < numShortLists; i++) {
		vector<Candi> tmp;
		tmp.reserve(candis.size() + invLists[i]->size());
		vector<Candi>::const_iterator candiIter = candis.begin();

		typename InvList::iterator invListIter = invLists[i]->begin();
		typename InvList::iterator invListEnd = invLists[i]->end();
    
		while(candiIter != candis.end() || invListIter != invListEnd) {
			if(candiIter == candis.end() || (invListIter != invListEnd && candiIter->id > *invListIter)) {
				tmp.push_back(Candi(*invListIter, 1));
				invListIter++;
			} 
			else if (invListIter == invListEnd || (candiIter != candis.end() && candiIter->id < *invListIter)) {
				tmp.push_back(Candi(candiIter->id, candiIter->count));
				candiIter++;
			} 
			else {
				tmp.push_back(Candi(candiIter->id, candiIter->count + 1));
				candiIter++;
				invListIter++;
			}
		}
		std::swap(candis, tmp);
	}
  
	// process long lists
	unsigned stop = invLists.size();
	for(unsigned j = numShortLists; j < stop; j++) {
		if(candis.size() < stop - j) break; // termination heuristic
    
		vector<Candi> tmp;
		tmp.reserve(candis.size());
		typename InvList::iterator listIter = invLists[j]->begin();
        
		for(unsigned i = 0; i < candis.size(); i++) {
			typename InvList::iterator start, end;
			expProbe(listIter, invLists[j]->end(), start, end, candis[i].id);      
			listIter = lower_bound(start, end, candis[i].id);
			if(listIter != invLists[j]->end() && *listIter == candis[i].id) {
				candis[i].count++;
				if(candis[i].count >= threshold) candiscur.push_back(candis[i]);
				else tmp.push_back(candis[i]);
			}
			else {
				if(candis[i].count + stop - j > threshold)
					tmp.push_back(candis[i]);
			}
		}
    
		std::swap(candis, tmp);       
	}

	for (unsigned i = 0; i < candis.size(); i++) {
		candiscur.push_back(candis[i]);
	}
}

template <class InvList>
void CSeqDB<InvList>::gramQuery(const CQuery& query, const unsigned ged, vector< vector<unsigned> >& similarGrams)
{
	const vector<string>& gramStrs = query.ngrams;
	vector<unsigned> gramCodes;
	vector<Candi> candis;
	vector<unsigned> results;
	int threshold = (int)this->gramGenLower->getGramLength();
	threshold = (int)(this->gramGenUpper->getGramLength()) - threshold + 1 - (int)ged*threshold;
	for (unsigned i = 0; i < gramStrs.size(); i++) {
		gramCodes.clear();
		this->gramGenLower->decompose(gramStrs[i], gramCodes);
		vector<InvList*> lists;
		this->getGramLists(gramCodes, this->gramListLower, lists);
		results.clear();
		if ((unsigned)threshold > lists.size()) {
			similarGrams.push_back(results);
			continue;
		}

		candis.clear();
		this->merge(lists, (unsigned)threshold, candis);
		for(unsigned j = 0; j < candis.size(); j++) {
			if (this->getEditDistance(this->gramStrUpper[candis[j].id], gramStrs[i], (int)ged))
				results.push_back(this->gramCodeUpper[candis[j].id]);
		}
		similarGrams.push_back(results);
	}
}

template <class InvList>
void CSeqDB<InvList>::init_threshold(CQuery& query) {
	int idlow, idhigh, ed;
	this->getIDBounds((int)query.length, 0, idlow, idhigh);
	idhigh = idhigh > (idlow + (int)query.constraint) ? idhigh : ((idlow + (int)query.constraint) < this->data->size() ? (idlow + (int)query.constraint) : this->data->size());
	idlow = idlow < (idhigh-(int)query.constraint) ? idlow : ((idhigh-(int)query.constraint)>0 ? (idhigh-(int)query.constraint) : 0);
	for (; idlow<idhigh; idlow++) {
		if (this->m_queue.size() < query.constraint) {
			ed = this->getRealEditDistance_nsd(this->data->at(idlow), query.sequence);
			this->processedData[idlow] = true;
			this->processed++;
			this->m_queue.push(queue_entry(idlow, (unsigned)ed));
			continue;
		}
		query.threshold = this->m_queue.top().m_dist;
		if (query.threshold == 0) {
			return;
		}	
		ed = this->getRealEditDistance_ns(this->data->at(idlow), query.sequence, (int)query.threshold);
		this->processedData[idlow] = true;
		this->processed++;
		if (ed < (int)query.threshold) {
			this->m_queue.pop();
			this->m_queue.push(queue_entry(idlow, (unsigned)ed));
			query.threshold = this->m_queue.top().m_dist;
			if (query.threshold == 0) {
				return;
			}
		}
	}
	query.threshold = this->m_queue.top().m_dist;
}

template <class InvList>
void CSeqDB<InvList>::knn_postprocess() {
	bool checked;
	int ed, i, j, appgram_count, idlow, idhigh, idnew;
	this->getIDBounds((int)this->theQuery->length, this->theQuery->threshold, idlow, idhigh);
	const int& querylen_index = this->datasize_groupindex[this->theQuery->length];
	if (idlow<idhigh) {
/*		if (this->fk > 1) {
			// @TODO should be testing the correctness
			// If use the f-queue
			vector<min_queue_entry> the_fqueue;
			bool checked_bound;
			// For those sizes no larger than query
			for (; idlow<=querylen_index; idlow++) {
				if (!this->processedData[idlow]) {
					if (this->dataCount[0][idlow]==0) {
						this->post_cand_low.push_back(idlow);
						continue;
					}
					checked = true;
					appgram_count = 0;
					for (i = 0; i < this->gram_maxed; i++) {
						appgram_count += this->dataCount[i][idlow];
						if (appgram_count < this->filter->tabUpQuery[this->theQuery->threshold][i]) {
							checked = false;
							break;
						}
					}
					if (checked)
					the_fqueue.push_back(min_queue_entry(idlow, this->dataCount[0][idlow]));
					if (the_fqueue.size() == this->fk || idlow == querylen_index) {
						sort(the_fqueue.begin(), the_fqueue.end());
						checked_bound = false;
						for (vector<min_queue_entry>::iterator iter=the_fqueue.begin(); iter!=the_fqueue.end(); iter++) {
							const unsigned& dataid = iter->m_sid;
							appgram_count = iter->m_dist;
							checked = true;
							if (checked_bound) {
								if (appgram_count < this->filter->tabUpQuery[this->theQuery->threshold][0]) {
									break;
								}
								for (i = 1; i < this->gram_maxed; i++) {
									appgram_count += this->dataCount[i][dataid];
									if (appgram_count < this->filter->tabUpQuery[this->theQuery->threshold][i]) {
										checked = false;
										break;
									}
								}
							}
							if (checked) {
								ed = this->getRealEditDistance_ns(this->data->at(dataid), this->theQuery->sequence, (int)this->theQuery->threshold);
								this->processed++;
								checked_bound = false;
								if (ed < (int)this->theQuery->threshold) {
									this->m_queue.pop();
									this->m_queue.push(queue_entry(dataid, (unsigned)ed));
									if (this->theQuery->threshold > this->m_queue.top().m_dist) {
										this->theQuery->threshold = this->m_queue.top().m_dist;
										if (this->theQuery->threshold == 0) {
											return;
										}
										if (this->filter->tabUpQuery[this->theQuery->threshold][this->gram_maxed-1] > 0)
										checked_bound = true;
										this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idnew, idhigh);
										if (idnew > querylen_index)
											break;
										idlow = idlow < (idnew-1) ? (idnew-1) : idlow;
									}
								}
							}
						}
						vector<min_queue_entry> tmp;
						swap(the_fqueue, tmp);
					}
				}
			}
			// For those sizes larger than query
			for (; idlow<idhigh; idlow++) {
				if (!this->processedData[idlow]) {
					if (this->dataCount[0][idlow]==0) {
						this->post_cand_up.push_back(idlow);
						continue;
					}
					checked = true;
					appgram_count = 0;
					const int& tmp_datasize = (int)this->data->at(idlow).size();
					for (i = 0; i < this->gram_maxed; i++) {
						appgram_count += this->dataCount[i][idlow];
						if ((appgram_count+this->filter->tabUp[this->theQuery->threshold][i]) < tmp_datasize) {
							checked = false;
							break;
						}
					}
					if (checked)
					the_fqueue.push_back(min_queue_entry(idlow, this->dataCount[0][idlow]));
					if (the_fqueue.size() == this->fk || idlow == idhigh-1) {
						sort(the_fqueue.begin(), the_fqueue.end());
						checked_bound = false;
						for (vector<min_queue_entry>::iterator iter=the_fqueue.begin(); iter!=the_fqueue.end(); iter++) {
							const unsigned& dataid = iter->m_sid;
							checked = true;
							if (checked_bound) {
								const int& _datasize = (int)this->data->at(dataid).size();
								if ((this->dataCount[0][dataid]+this->filter->tabUp[this->theQuery->threshold][0]) < _datasize) {
									break;
								}
								appgram_count = this->dataCount[0][dataid];
								for (i = 1; i < this->gram_maxed; i++) {
									appgram_count += this->dataCount[i][dataid];
									if ((appgram_count+this->filter->tabUp[this->theQuery->threshold][i]) < _datasize) {
										checked = false;
										break;
									}
								}
							}
							if (checked) {
								ed = this->getRealEditDistance_ns(this->data->at(dataid), this->theQuery->sequence, (int)this->theQuery->threshold);
								this->processed++;
								checked_bound = false;
								if (ed < (int)this->theQuery->threshold) {
									this->m_queue.pop();
									this->m_queue.push(queue_entry(dataid, (unsigned)ed));
									if (this->theQuery->threshold > this->m_queue.top().m_dist) {
										this->theQuery->threshold = this->m_queue.top().m_dist;
										if (this->theQuery->threshold == 0) {
											return;
										}
										checked_bound = true;
										this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idnew, idhigh);
										if (idnew >= idhigh)
											break;
										idlow = idlow < (idnew-1) ? (idnew-1) : idlow;
									}
								}
							}
						}
						vector<min_queue_entry> tmp;
						swap(the_fqueue, tmp);
					}
				}
			}
		}*/
//		else {
			// For those sizes no larger than query
			for (; idlow<=querylen_index; idlow++) {
				if (!this->processedData[idlow]) {
					if (this->dataCount[0][idlow]==0) {
						// Even do not have a similar gram
						if (this->dataCount[1][idlow]==0)
						continue;
						this->post_cand_low.push_back(idlow);
						continue;
					}
					// Check the bounds on approximate ngrams
					checked = true;
					appgram_count = 0;
					for (i = 0; i < this->gram_maxed; i++) {
						appgram_count += this->dataCount[i][idlow];
						if (appgram_count < this->filter->tabUpQuery[this->theQuery->threshold][i]) {
							checked = false;
							break;
						}
					}
					if (checked) {
						ed = this->getRealEditDistance_ns(this->data->at(idlow), this->theQuery->sequence, (int)this->theQuery->threshold);
						this->processed++;
						if (ed < (int)this->theQuery->threshold) {
							this->m_queue.pop();
							this->m_queue.push(queue_entry(idlow, (unsigned)ed));
							if (this->theQuery->threshold > this->m_queue.top().m_dist) {
								this->theQuery->threshold = this->m_queue.top().m_dist;
								if (this->theQuery->threshold == 0) {
									return;
								}
								this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idnew, idhigh);
								if (idnew > querylen_index)
									break;
								idlow = idlow < (idnew-1) ? (idnew-1) : idlow;
							}
						}
					}
				}
			}
			// For those sizes larger than query
			for (; idlow<idhigh; idlow++) {
				if (!this->processedData[idlow]) {
					if (this->dataCount[0][idlow]==0) {
						// Even do not have a similar gram
						if (this->dataCount[1][idlow]==0)
						continue;
						this->post_cand_up.push_back(idlow);
						continue;
					}
					// Check the bounds on approximate ngrams
					checked = true;
					appgram_count = 0;
					const int& tmp_datasize = (int)this->data->at(idlow).size();
					for (i = 0; i < this->gram_maxed; i++) {
						appgram_count += this->dataCount[i][idlow];
						if ((appgram_count + this->filter->tabUp[this->theQuery->threshold][i]) < tmp_datasize) {
							checked = false;
							break;
						}
					}
					if (checked) {
						ed = this->getRealEditDistance_ns(this->data->at(idlow), this->theQuery->sequence, (int)this->theQuery->threshold);
						this->processed++;
						if (ed < (int)this->theQuery->threshold) {
							this->m_queue.pop();
							this->m_queue.push(queue_entry(idlow, (unsigned)ed));
							if (this->theQuery->threshold > this->m_queue.top().m_dist) {
								this->theQuery->threshold = this->m_queue.top().m_dist;
								if (this->theQuery->threshold == 0) {
									return;
								}
								this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idnew, idhigh);
								if (idnew >= idhigh)
									break;
								idlow = idlow < (idnew-1) ? (idnew-1) : idlow;
							}
						}
					}
				}
			}
//		}
		// Check approximate ngrams
		this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idlow, idhigh);
		for (j = 0; j < this->post_cand_low.size(); j++) {
			const unsigned& dataid = this->post_cand_low[j];
			if (dataid >= idlow && dataid < idhigh) {
				checked = true;
				appgram_count = 0;
				for (i = 1; i < this->gram_maxed; i++) {
					appgram_count += this->dataCount[i][dataid];
					if (appgram_count < this->filter->tabUpQuery[this->theQuery->threshold][i]) {
						checked = false;
						break;
					}
				}
				if (checked) {
					ed = this->getRealEditDistance_ns(this->data->at(dataid), this->theQuery->sequence, (int)this->theQuery->threshold);
					this->processed++;
					if (ed < (int)this->theQuery->threshold) {
						this->m_queue.pop();
						this->m_queue.push(queue_entry(dataid, (unsigned)ed));
						if (this->theQuery->threshold > this->m_queue.top().m_dist) {
							this->theQuery->threshold = this->m_queue.top().m_dist;
							if (this->theQuery->threshold == 0) {
								return;
							}
							this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idlow, idhigh);
						}
					}
				}
			}
		}
		for (j = 0; j < this->post_cand_up.size(); j++) {
			const unsigned& dataid = this->post_cand_up[i];
			if (dataid >= idlow && dataid < idhigh) {
				checked = true;
				appgram_count = 0;
				const int& tmp_datasize = (int)this->data->at(dataid).size();
				for (i = 1; i < this->gram_maxed; i++) {
					appgram_count += this->dataCount[i][dataid];
					if ((appgram_count + this->filter->tabUp[this->theQuery->threshold][i]) < tmp_datasize) {
						checked = false;
						break;
					}
				}
				if (checked) {
					ed = this->getRealEditDistance_ns(this->data->at(dataid), this->theQuery->sequence, (int)this->theQuery->threshold);
					this->processed++;
					if (ed < (int)this->theQuery->threshold) {
						this->m_queue.pop();
						this->m_queue.push(queue_entry(dataid, (unsigned)ed));
						if (this->theQuery->threshold > this->m_queue.top().m_dist) {
							this->theQuery->threshold = this->m_queue.top().m_dist;
							if (this->theQuery->threshold == 0) {
								return;
							}
							this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idlow, idhigh);
						}
					}
				}
			}
		}
	}
}

// Accumulate the frequency of approximate ngrams
// the largest value of ged is decided by upper gram length according to the CA strategy
template <class InvList>
void CSeqDB<InvList>::accumulateFrequency(long ged)
{
	vector<InvList*> lists;
	if (ged == 0) {
		this->getGramLists(this->theQuery->gramCodes, this->gramListUpper, lists);
	}
	else {
		vector< vector<unsigned> > similarGrams;
		this->gramQuery(*(this->theQuery), ged, similarGrams);
		for (unsigned i = 0; i < similarGrams.size(); i++) {
			this->getGramLists(similarGrams[i], this->gramListUpper, lists);
		}
	}
	int idlow,  idhigh;
	this->getIDBounds((int)this->theQuery->length, (int)this->theQuery->threshold, idlow, idhigh);
	typename InvList::iterator iterlow, iterhigh;
	for(unsigned i = 0; i < lists.size(); i++) {
		iterlow = lists[i]->begin();
		iterhigh = lists[i]->end();
		if (idlow < idhigh) {
			iterlow = lower_bound(lists[i]->begin(), lists[i]->end(), idlow);
			iterhigh = upper_bound(lists[i]->begin(), lists[i]->end(), idhigh);
			for(; iterlow < iterhigh; ++iterlow) {
				if (!this->processedData[*iterlow])
					this->dataCount[ged][*iterlow]++;
			}
		}
	}
}

template <class InvList>
void CSeqDB<InvList>::getIDBounds(const int querysize, const int threshold, int& idlow, int& idhigh)
{
	int lowlen = querysize- threshold > 0 ? querysize - threshold : 0;
	int highlen = querysize + threshold;
	if ((unsigned)lowlen > this->datasizes.front().m_dist) {
		vector<queue_entry>::iterator _iterLow = lower_bound(this->datasizes.begin(), this->datasizes.end(), queue_entry(0, (unsigned)lowlen));
		idlow = (int)(_iterLow->m_sid);
	}
	else {
		idlow = this->datasizes.front().m_sid;
	}
	if ((unsigned)highlen < this->datasizes.back().m_dist) {
		vector<queue_entry>::iterator _iterHigh = upper_bound(this->datasizes.begin(), this->datasizes.end(), queue_entry(0, (unsigned)highlen));
		idhigh = (int)(_iterHigh->m_sid);
	}else {
		idhigh = this->datasizes.back().m_sid;
	}
}

template <class InvList>
void CSeqDB<InvList>::queryGramLists(const CQuery& query)
{
	vector<InvList*> invLists;
	this->getGramLists(query.gramCodes, this->gramListUpper, invLists);
	fprintf(stdout, "%u\n", invLists.size());
	typename InvList::iterator iterlow;
	for(unsigned i = 0; i < invLists.size(); i++) {
		for(iterlow = invLists[i]->begin(); iterlow != invLists[i]->end(); iterlow++) {
			fprintf(stdout, "%u\t", *iterlow);
		}
		fprintf(stdout, "\n", *iterlow);
	}
}

#endif /* _SEQDB_H_ */
