#include "params.h"  // access to declarations of global parameter values

map<int, map<string, string> > read_parameters_file(const string &parameters_fn)
{
	map<int, map<string, string> > params_by_block;
	map<string,string> params;
	ifstream paramfile(parameters_fn.c_str());
	string line;
	int block = -1;
	regex query("DEME");
	while(getline(paramfile, line)) {
		istringstream iss(line); // does it need to be line.c_str()?????????
		string key, nextone, value;
		iss >> key;
		if (regex_search(key, query)) {  // check if entering next deme block
			params_by_block[block] = params;
			params.clear();
			++block;
		} else {
			while (iss >> nextone)
				value += nextone + " ";
			params[key] = value;
		}
	}
	params_by_block[block] = params; // assign values for final deme
	return  params_by_block;
}

map<int, map<string, string> > p = read_parameters_file("parameters"); // change name if parameters file not named parameters

vector<int> get_multi_int_param(const string &key, map<string, string> &parameters)
{
	vector<int> vec;
	istringstream iss(parameters[key].c_str());
	string param;
	while(getline(iss, param, ' '))
		vec.push_back(atoi(param.c_str()));
	return vec;
}

vector<double> get_multi_double_param(const string &key, map<string, string> &parameters)
{
	vector<double> vec;
	istringstream iss(parameters[key].c_str());
	string param;
	while(getline(iss, param, ' '))
		vec.push_back(atof(param.c_str()));
	return vec;
}

vector<int> create_pop_schedule()
{
	vector<int> ps;
	int i=0;
	int cursize = popsize;
	for (int step = 0; step < demography.size(); ++step) {
		for (; i<dem_start_gen[step]; ++i)
		{
			ps.push_back(cursize);
		}
		for (; i <= dem_end_gen[step]; ++i) {
			switch(demography[step]) {
				case 0: ps.push_back(cursize);  // no size change
					    break;
				case 1: if (i == dem_start_gen[step])
							cursize += dem_parameter[step];
						ps.push_back(cursize); // instantaneous
						break;
				case 2: cursize += dem_parameter[step];
						ps.push_back(cursize);  // linear
						break;
				case 3: cursize *= exp(dem_parameter[step]);  // exponential
						ps.push_back(cursize);
						break;
				case 4: cursize = (carrying_cap[step] * cursize) /
									(cursize + (carrying_cap[step] - cursize)*exp(-1*dem_parameter[step]));  // logistic
						ps.push_back(cursize);
			}
		}
	}
	return ps;
}

// variable names
int popsize, sampsize, seqlength, sampfreq, hotrecStart, hotrecStop, windowSize, windowStep, pop_num, runlength, printhapfreq, diploid_sample, epistaticLoci;
double mutrate, recrate, hotrecrate;
bool useRec, useHotRec, getWindowStats, modelMigration, trackAlleleBirths;
vector<int> demography;
vector<double> dem_parameter;
vector<int> dem_start_gen;
vector<int> dem_end_gen;
vector<int> carrying_cap;
vector<double> m;
vector<string> mscommand, hapFile;
vector<bool> useMS, useHapStartingData;
vector<int> birthgen, extinctgen;
map<int, vector<int> > pop_schedule, splitgenesis, mergegenesis;
map<int, vector<double> > sellocus, possel, nfdsel, negsel, polysel;
double baseTraitValue, enviroSD;
bool  polygenicTrait, epistatic;
vector<double> alleleEffects, dominanceEffects;
vector<int> epistaticTypes;

int process_parameters() { // replaces old version of function
	for (auto iter = p.begin(); iter!=p.end(); ++iter ) {
		map<string, string> parameters = iter->second;
		if (iter->first == -1) {  // block of global parameters
			mutrate = atof(parameters["mutrate"].c_str());
			recrate = atof(parameters["recrate"].c_str());
			hotrecrate = atof(parameters["hotrecrate"].c_str());
			useRec = atoi(parameters["useRec"].c_str());
			useHotRec = atoi(parameters["useHotRec"].c_str());
			hotrecStart = atoi(parameters["hotrecStart"].c_str());
			hotrecStop = atoi(parameters["hotrecStop"].c_str());
			sampsize = atoi(parameters["sampsize"].c_str());
			seqlength = atof(parameters["seqlength"].c_str()); // covernsion using atof() enables use of e notation in parameters file
			sampfreq = atoi(parameters["sampfreq"].c_str());
			getWindowStats = atoi(parameters["getWindowStats"].c_str());
			windowSize = atoi(parameters["windowSize"].c_str());
			windowStep = atoi(parameters["windowStep"].c_str());
			pop_num = atoi(parameters["pop_num"].c_str());
			runlength = atoi(parameters["runlength"].c_str());
			m = get_multi_double_param("migration_rates", parameters);
			printhapfreq = atoi(parameters["printhapfreq"].c_str());
			diploid_sample = atoi(parameters["diploid_sample"].c_str());
			trackAlleleBirths = atoi(parameters["trackAlleleBirths"].c_str());
			modelMigration = atoi(parameters["modelMigration"].c_str());
			polygenicTrait = atoi(parameters["polygenicTrait"].c_str());
			baseTraitValue = atof(parameters["baseTraitValue"].c_str());
			alleleEffects = get_multi_double_param("alleleEffects", parameters);
			if (polygenicTrait && seqlength > 1 && alleleEffects.size() == 1) // then set all effect sizes to this value
				for (int i=1; i<seqlength; ++i)
					alleleEffects.push_back(alleleEffects[0]);
			dominanceEffects = get_multi_double_param("dominanceEffects", parameters);
			if (seqlength > 1 && dominanceEffects.size() == 1)
				for (int i=1; i<seqlength; ++i)
					dominanceEffects.push_back(dominanceEffects[0]);
			enviroSD = atof(parameters["enviroSD"].c_str());
			epistatic =  atoi(parameters["epistatic"].c_str());
			epistaticTypes = get_multi_int_param("epistaticTypes", parameters);
			epistaticLoci = epistaticTypes.size();
			if (epistatic)   {
				seqlength = seqlength + epistaticLoci;
				for (int i = 0; i<epistaticLoci; ++i)
					alleleEffects.push_back(0.);
			} else {
				epistaticLoci = 0;
			}
		} else {   // block of deme parameters
			popsize = atoi(parameters["popsize"].c_str());
			demography =  get_multi_int_param("demography", parameters);
			dem_parameter = get_multi_double_param("dem_parameter", parameters);
			dem_start_gen =  get_multi_int_param("dem_start_gen", parameters);
			dem_end_gen = get_multi_int_param("dem_end_gen", parameters);
			carrying_cap = get_multi_int_param("carrying_cap", parameters);
			birthgen.push_back( atoi(parameters["birthgen"].c_str()) );
			extinctgen.push_back( atoi(parameters["extinctgen"].c_str()) );
			useMS.push_back( atoi(parameters["useMS"].c_str()) );
			mscommand.push_back( parameters["mscommand"] );
			pop_schedule[iter->first] = create_pop_schedule();
			splitgenesis[iter->first] = get_multi_int_param("splitgenesis", parameters);
			mergegenesis[iter->first] = get_multi_int_param("mergegenesis", parameters);
			sellocus[iter->first] = get_multi_double_param("sellocus", parameters);
			possel[iter->first] = get_multi_double_param("possel", parameters);
			nfdsel[iter->first] = get_multi_double_param("nfdsel", parameters);
			negsel[iter->first] = get_multi_double_param("negsel", parameters);
			useHapStartingData.push_back( atoi(parameters["useHapStartingData"].c_str()) );
			hapFile.push_back( parameters["hapFile"] );
			polysel[iter->first] = get_multi_double_param("polysel", parameters);
		}
	}
	return 1;
}

int good_parameters = process_parameters();
double* a = &m[0]; // matrix takes array argument
Matrix<double> mig(pop_num, pop_num, a);
int ed[] = {0, 0, 0, 1, 0, 0, 2, 0, 0};
Matrix<int> epiD(3, 3, ed); // dominance epistasis
int er[] = {0, 0, 0, 1, 1, 0, 2, 2, 0};
Matrix<int> epiR(3, 3, er);   // recessive epistasis
int aa[] = {0, 1, 2, 1, 1, 1, 2, 1, 0};
Matrix<int> epiAA(3, 3, aa); // AxA epistasis
int ad[] = {2, 0, 2, 1, 1, 1, 0, 2, 0};
Matrix<int> epiAD(3, 3, ad); // AxD epistasis
int da[] = {0, 2, 0, 1, 1, 1, 2, 0, 2};
Matrix<int> epiDA(3, 3, da); // DxA epistasis
int dd[] = {0, 2, 0, 2, 0, 2, 0, 2, 0};
Matrix<int> epiDD(3, 3, dd); // DxD epistasis
