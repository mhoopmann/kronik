/*
Copyright 2007-2015, Michael R. Hoopmann, University of Washington

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef _CKRONIK_H
#define _CKRONIK_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

using namespace std;

#define FWHMCONST 2.3548200450309493820231386529194
#define SQRTTWO 1.4142135623730950488016887242097

typedef struct sPep{
  char charge;
	char mods;
  float intensity;
  double monoMass;
	double basePeak;
  double xCorr;
} sPep;

typedef struct sScan {

	char file;
	int scanNum;
	float rTime;
	vector<sPep> *vPep;

	//Constructors & Destructor
	sScan(){vPep = new vector<sPep>;}
	sScan(const sScan& s){
		vPep = new vector<sPep>;
		for(unsigned int i=0;i<s.vPep->size();i++) vPep->push_back(s.vPep->at(i));
		scanNum = s.scanNum;
		file = s.file;
		rTime=s.rTime;
	}
	~sScan(){delete vPep;}

  //Copy operator
	sScan& operator=(const sScan& s){
		if(this!=&s){
			delete vPep;
			vPep = new vector<sPep>;
			for(unsigned int i=0;i<s.vPep->size();i++) vPep->push_back(s.vPep->at(i));
			scanNum = s.scanNum;
			file = s.file;
			rTime=s.rTime;
		}
		return *this;
	}

  //Clear
	void clear(){
		delete vPep;
		vPep = new vector<sPep>;
	}

  static int compareIntRev1(const void *p1, const void *p2){
    const sPep d1 = *(sPep *)p1;
    const sPep d2 = *(sPep *)p2;
    if(d1.intensity>d2.intensity) return -1;
    else if(d1.intensity<d2.intensity) return 1;
    else return 0;
  }

	static int compareMonoMass1(const void *p1, const void *p2){
    const sPep d1 = *(sPep *)p1;
    const sPep d2 = *(sPep *)p2;
    if(d1.monoMass<d2.monoMass) return -1;
    else if(d1.monoMass>d2.monoMass) return 1;
    else return 0;
  }

  void sortIntRev(){
    if(vPep->size()==0) return;
    qsort(&vPep->at(0),vPep->size(),sizeof(sPep),compareIntRev1);
  }

	void sortMonoMass(){
    if(vPep->size()==0) return;
    qsort(&vPep->at(0),vPep->size(),sizeof(sPep),compareMonoMass1);
  }

} sScan;

typedef struct sProfileData {
  bool interpolated;

  int scanNum;

  float intensity;
  float rTime;

  double monoMass;
  double xCorr;

  //Constructor
  sProfileData(){
    interpolated=false;
    scanNum=0;
    intensity=0.0f;
    rTime=0.0f;
    monoMass=0.0;
    xCorr=0.0;
  }
} sProfileData;

typedef struct sPepProfile {
  char charge;
	char mods;

	int lowScan;
  int highScan;
  int bestScan;
  int MS2Events;

  unsigned int datapoints;

  float rTime;
  float firstRTime;
  float lastRTime;
  float intensity;
  float sumIntensity;

  double monoMass;
  double basePeak;
  double xCorr;
	double gaussian[4];
	double gaussR2;

	//string gene;
	//string sequence;
	//char* gene;
	//char* sequence;
  char gene[32];
  char sequence[64];

  sProfileData* profile;

  //Constructors
  sPepProfile() {
    profile = new sProfileData[1];
    sProfileData d;
    profile[0]=d;
    datapoints=1;
    lowScan=0;
    highScan=0;
    bestScan=0;
    MS2Events=0;
    charge=0;
    rTime=0.0f;
    firstRTime=0.0f;
    lastRTime=0.0f;
    intensity=0.0f;
    sumIntensity=0.0f;
    monoMass=0.0;
    basePeak=0.0;
    xCorr=0.0;
		gaussian[0]=gaussian[1]=gaussian[2]=gaussian[3]=0.0;
		gaussR2=0.0;
    mods=0;
    //sequence=new char[1];
		//gene=new char[1];
		sequence[0]='\0';
		gene[0]='\0';
  }
  sPepProfile(const unsigned int i) {
    profile = new sProfileData[i];
    sProfileData d;
    unsigned int j;
    for(j=0;j<i;j++) profile[j]=d;
    datapoints=i;
    lowScan=0;
    highScan=0;
    bestScan=0;
    charge=0;
    MS2Events=0;
    rTime=0.0f;
    firstRTime=0.0f;
    lastRTime=0.0f;
    intensity=0.0f;
    sumIntensity=0.0f;
    monoMass=0.0;
    basePeak=0.0;
    xCorr=0.0;
		gaussian[0]=gaussian[1]=gaussian[2]=gaussian[3]=0.0;
		gaussR2=0.0;
    mods=0;
		//sequence=new char[1];
		//gene=new char[1];
		sequence[0]='\0';
		gene[0]='\0';
  }
  sPepProfile(const sPepProfile& p){
    unsigned int i;
    profile = new sProfileData[p.datapoints];
    for(i=0;i<p.datapoints;i++) profile[i]=p.profile[i];
    lowScan = p.lowScan;
    highScan = p.highScan;
    bestScan = p.bestScan;
    charge = p.charge;
    MS2Events = p.MS2Events;
    datapoints = p.datapoints;
    rTime = p.rTime;
    firstRTime = p.firstRTime;
    lastRTime = p.lastRTime;
    intensity = p.intensity;
    sumIntensity = p.sumIntensity;
    monoMass = p.monoMass;
    basePeak = p.basePeak;
    xCorr = p.xCorr;
		gaussian[0]=p.gaussian[0];
		gaussian[1]=p.gaussian[1];
		gaussian[2]=p.gaussian[2];
		gaussian[3]=p.gaussian[3];
		gaussR2=p.gaussR2;
    mods=p.mods;
		//delete [] sequence;
		//sequence = new char[strlen(p.sequence)+1];
		strcpy(sequence,p.sequence);
		//delete [] gene;
		//gene = new char[strlen(p.gene)+1];
		strcpy(gene,p.gene);
  }

  //Destructor
  ~sPepProfile(){
		//delete [] gene;
		//delete [] sequence;
    delete [] profile;
  }

  //Copy operator
  sPepProfile& operator=(const sPepProfile& p){
    if(this!=&p){
      unsigned int i;
      delete [] profile;
      profile = new sProfileData[p.datapoints];
      for(i=0;i<p.datapoints;i++) profile[i]=p.profile[i];
      lowScan = p.lowScan;
      highScan = p.highScan;
      bestScan = p.bestScan;
      MS2Events = p.MS2Events;
      charge = p.charge;
      datapoints = p.datapoints;
      rTime = p.rTime;
      firstRTime = p.firstRTime;
      lastRTime = p.lastRTime;
      intensity = p.intensity;
      sumIntensity = p.sumIntensity;
      monoMass = p.monoMass;
      basePeak = p.basePeak;
      xCorr = p.xCorr;
			gaussian[0]=p.gaussian[0];
			gaussian[1]=p.gaussian[1];
			gaussian[2]=p.gaussian[2];
			gaussian[3]=p.gaussian[3];
			gaussR2=p.gaussR2;
      mods=p.mods;
			//delete [] sequence;
			//sequence = new char[strlen(p.sequence)+1];
			strcpy(sequence,p.sequence);
			//delete [] gene;
			//gene = new char[strlen(p.gene)+1];
			strcpy(gene,p.gene);
    }
    return *this;
  }

  void setPoints(const unsigned int& i){
    delete [] profile;
    profile = new sProfileData[i];
    datapoints=i;
    sProfileData d;
    for(unsigned int j=0;j<i;j++) profile[j]=d;
  }

	static int compareMonoMass2(const void *p1, const void *p2){
    const sProfileData d1 = *(sProfileData *)p1;
    const sProfileData d2 = *(sProfileData *)p2;
    if(d1.monoMass<d2.monoMass) return -1;
    else if(d1.monoMass>d2.monoMass) return 1;
    else return 0;
  }

  static int compareScanNum2(const void *p1, const void *p2){
    const sProfileData d1 = *(sProfileData *)p1;
    const sProfileData d2 = *(sProfileData *)p2;
    if(d1.scanNum<d2.scanNum) return -1;
    else if(d1.scanNum>d2.scanNum) return 1;
    else return 0;
  }

  static int compareIntRev2(const void *p1, const void *p2){
    const sProfileData d1 = *(sProfileData *)p1;
    const sProfileData d2 = *(sProfileData *)p2;
    if(d1.intensity>d2.intensity) return -1;
    else if(d1.intensity<d2.intensity) return 1;
    else return 0;
  }

  void sortIntRev(){
    if(datapoints==0) return;
    qsort(profile,datapoints,sizeof(sProfileData),compareIntRev2);
  }

  void sortMonoMass(){
    if(datapoints==0) return;
    qsort(profile,datapoints,sizeof(sProfileData),compareMonoMass2);
  }

  void sortScanNum(){
    if(datapoints==0) return;
    qsort(profile,datapoints,sizeof(sProfileData),compareScanNum2);
  }

} sPepProfile;

typedef struct iTwo{
  int scan;
  int pep;
} iTwo;

typedef struct intenIndex{
	float intensity;
	int scan;
	int pep;
} intenIndex;

class CKronik {
public:

  //Constructors and destructors
  CKronik();
  CKronik(const CKronik& c);
  ~CKronik();

  //Operator overrides
  CKronik& operator=(const CKronik& c);
  sPepProfile& operator[ ](const unsigned int& i);

  //Statistics and math functions  
  bool		pearson(int i1, int i2, bool byScan, bool interpolate, double& rval, double& pval, double& slope, double& intercept, float& sum1, float& sum2);
	double	polynomialBestFit(sProfileData* profile, int sz, float max, double* coeff, int degree=2);
  
  //Parameter accessors and modifiers
	void setGaussFit(bool b);
  void setPPMTol(double d);
  void setMatchTol(int i);
  void setGapTol(int i);

  //Automation
  int getPercent();
  bool loadHK(const char* in);
  bool processHK(const char* in, const char* out="\0");

  //Tools
  bool getRT(int scanNum, float& rt);

  //Filters
  void filterRT(float rt1, float rt2);
  void removeContaminants(float rt);
	void removeMass(double m1, double m2);

  //Data manipulation functions
  sPepProfile& at(unsigned int i);
  void add(sPepProfile& p);
  void clear();
  void clearHK();
  void erase(unsigned int i);
  unsigned int size();

  //Sorting functions
  void sortBasePeak();
  void sortBestScan();
  void sortMonoMass();
  void sortFirstRTime();
  void sortIntensityRev();
  void sortRTime();
	void sortSumIntensityRev();

	//Unsorted
	vector<string> vFile;
	vector<string> vMods;

protected:
private:
	int binarySearch(vector<sScan>& v, int index, double mass, int charge);
  bool findMax(vector<sScan>& v, int& s, int& p);
  double interpolate(int x1, int x2, double y1, double y2, int x);
  
  //Statistics functions
  double betai(double a, double b, double x);
  double betacf(double a, double b, double x);
  double gammaln(double xx);

  //Data Members: data analysis
  vector<sPepProfile> vPeps;
  vector<sScan> hkData;

  //Data Members: parameters
	bool		bGaussFit;	//default false
  double	dPPMTol;		//default 10.0
  int			iGapTol;    //default 1
  int			iMatchTol;  //Default 3
  int			iPercent;

  //Sorting Functions
  void sortPeptide();
  static int compareBP3(const void *p1, const void *p2);
  static int compareBS3(const void *p1, const void *p2);
  static int compareMM3(const void *p1, const void *p2);
  static int compareRT3(const void *p1, const void *p2);
  static int compareFRT3(const void *p1, const void *p2);
	static int compareIntRev3(const void *p1, const void *p2);
  static int compareIRev3(const void *p1, const void *p2);
	static int compareSumIRev3(const void *p1, const void *p2);

};

#endif
