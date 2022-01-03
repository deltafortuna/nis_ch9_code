#ifndef PARAMS_H
#define PARAMS_H

// GLOBAL PARAMETERS
// must be constant at compile time
extern const int bitlength = 4000;

// the following are defined by values in the parameters file
extern int popsize;
extern double mutrate;
extern int seqlength;
extern int sampsize;
extern int sampfreq;
extern vector<bool> useMS;
extern vector<string> mscommand;
extern vector<int> demography;
extern vector<double> dem_parameter;
extern vector<int> dem_start_gen;
extern vector<int> dem_end_gen;
extern vector<int> carrying_cap;
extern map<int, vector<int> > pop_schedule;
extern double recrate;
extern double hotrecrate;
extern bool useRec;
extern bool useHotRec;
extern int hotrecStart;
extern int hotrecStop;
extern bool getWindowStats;
extern int windowSize;
extern int windowStep;
extern int printhapfreq;
extern int diploid_sample;
extern bool trackAlleleBirths;
extern bool modelMigration;
extern map<int, vector<int> > splitgenesis;
extern map<int, vector<int> > mergegenesis;
extern map<int, vector<double> > sellocus, possel, nfdsel, negsel;
extern bool polygenicTrait;
extern double baseTraitValue;
extern vector<double> alleleEffects;
extern vector<double> dominanceEffects;
extern double enviroSD;
extern vector<int> epistaticTypes;
extern int epistaticLoci;
extern bool epistatic;
extern Matrix<int> epiD;
extern Matrix<int> epiR;
extern Matrix<int> epiAA;
extern Matrix<int> epiAD;
extern Matrix<int> epiDA;
extern Matrix<int> epiDD;
extern vector<bool> useHapStartingData;
extern vector<string> hapFile;
extern map<int, vector<double> > polysel;
#endif
