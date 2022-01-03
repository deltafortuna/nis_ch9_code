#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "allele.h"

class Individual {

private:
	vector<vector<int>> sequences;
	double trait_value;
	void calculate_trait_value() {
		normal_distribution<double> Ve(0, enviroSD); // environmental variance
		trait_value = baseTraitValue;
		vector<int> genotypes(seqlength); // initialize with seqlength entries set to zero
		for (int i=0; i<2; ++i)
			for (auto iter = sequences[i].begin(); iter != sequences[i].end(); ++iter)
				genotypes[*iter]++;
		for (int i=epistaticLoci; i<seqlength-epistaticLoci; ++i) { // calculate contributions of non-interacting loci
			if (genotypes[i] == 1)
				trait_value += (alleleEffects[i] + dominanceEffects[i]);
			if (genotypes[i] == 2)
				trait_value += (2 * alleleEffects[i]);
		}
		if (epistatic) {
			int mainLocusIndex = 0;
			// calcluate contributions of epistatic pairs
			for (int i=seqlength-epistaticLoci; i < seqlength; ++i) {  // index of epistatic locus
				double traitMod;
				switch(epistaticTypes[mainLocusIndex]) {
					case 0:	// dominance epistasis as in figure 9.5B
							traitMod = epiD[genotypes[mainLocusIndex]][genotypes[i]] * alleleEffects[mainLocusIndex];
							trait_value +=  traitMod;
								break;
					case 1: // recessive epistasis as in figure 9.5A
							traitMod = epiR[genotypes[mainLocusIndex]][genotypes[i]] * alleleEffects[mainLocusIndex];
							trait_value += traitMod;
								break;
					case 2: // AxA epistasis
							traitMod= epiAA[genotypes[mainLocusIndex]][genotypes[i]] * alleleEffects[mainLocusIndex];
							trait_value += traitMod;
								break;
					case 3: // AxD epistasis
							traitMod = epiAD[genotypes[mainLocusIndex]][genotypes[i]] * alleleEffects[mainLocusIndex];
							trait_value += traitMod;
							break;
					case 4: // DxA epistasis
							traitMod = epiDA[genotypes[mainLocusIndex]][genotypes[i]] * alleleEffects[mainLocusIndex];
							trait_value += traitMod;
							break;
					case 5: // DxD epistasis
							traitMod = epiDD[genotypes[mainLocusIndex]][genotypes[i]] * alleleEffects[mainLocusIndex];
							trait_value += traitMod;
				}
				mainLocusIndex++;
			}
		}
		trait_value += Ve(n); // add environmental effect
		if (trait_value < baseTraitValue)
				trait_value = baseTraitValue;
	}

	void remove_allele_by_position (int seqnum, int position) {
		auto pos = find(sequences[seqnum].begin(), sequences[seqnum].end(), position);
		if (pos != sequences[seqnum].end())  // ensures position in vector
			sequences[seqnum].erase(pos);
	}

	void resolve_crossover (vector<int>& breaks) {
		int numbreaks = breaks.size();
		vector<int>::iterator lower, upper;
		map<int, vector<int> > segments;
		vector<int> newvec;

		for (int seq = 0; seq < 2; ++seq) {
			lower = sequences[seq].begin();
			for (int i=0; i<numbreaks; ++i) {
				upper  = upper_bound(sequences[seq].begin(), sequences[seq].end(), breaks[i]);
				newvec.assign(lower, upper);
				lower = upper;
				segments[i + seq*(numbreaks+1)] = newvec;
			}
			newvec.assign(lower, sequences[seq].end());
			segments[numbreaks + seq*(numbreaks+1)] = newvec;
		}

		for(int i=0; i<= numbreaks; ++i) {
			if (i%2 == 0) {
				if (i == 0) {sequences[0] = segments[0]; sequences[1] = segments[numbreaks+1];}
				else {
					sequences[0].insert(sequences[0].end(), segments[i].begin(), segments[i].end());
					sequences[1].insert(sequences[1].end(), segments[i+numbreaks+1].begin(), segments[i+numbreaks+1].end());
				}
			} else{
				sequences[0].insert(sequences[0].end(), segments[i+numbreaks+1].begin(), segments[i+numbreaks+1].end());
				sequences[1].insert(sequences[1].end(), segments[i].begin(), segments[i].end());
			}
		}
	}

public:
	inline vector<int> get_sequence(int whichseq) { return sequences[whichseq]; }
	inline vector<vector<int> > get_sequences() { return sequences; }
	inline double get_trait_value() {return trait_value;}

	const vector<int> & get_seq(int whichseq) { return sequences[whichseq]; }

	double get_polygenic_fitness(int popn) {
		double fitness;
		double theta = polysel[popn][1];
		double s = polysel[popn][2];
		if (polysel[popn][0] == 1) { // Gaussian fitness function
			fitness = exp(-1 * s * (pow (trait_value - theta, 2)));
		} else if (polysel[popn][0] == 2) { // Evolutionary constraints
			if (trait_value < polysel[popn][1])
				fitness = 1 - ((polysel[popn][1] - trait_value) * polysel[popn][3]); // [3] sel. strength
			else if (trait_value > polysel[popn][2])
				fitness = 1 - ((trait_value - polysel[popn][2]) * polysel[popn][3]);
			else
				fitness = 1;
		}
		if (fitness < 0) fitness = 0;
		if (fitness > 1) fitness = 1;
		return fitness;
	}

	void remove_fixed_allele(int to_remove) {
		for (int i = 0; i<2; ++i) {
			vector<int>::iterator p = find(sequences[i].begin(), sequences[i].end(), to_remove);
			if (p != sequences[i].end()) // i.e., element not found
		    		sequences[i].erase(p);
		 }
	}

	int get_genotype(int pos) { // returns number of ancestral alleles (i.e., 0s)
		int geno = 2;
		if ( find( sequences[0].begin(), sequences[0].end(), pos) != sequences[0].end() ) --geno;
		if ( find( sequences[1].begin(), sequences[1].end(), pos) != sequences[1].end() ) --geno;
		return geno;
	}

	int get_negsel_genotype(vector<int>& sites) { // returns number of seqs with a deleterious allele
		int geno = 0;
		for (int i=0; i<2; ++i) {
			vector<int> intersectionality;
			set_intersection(sites.begin(), sites.end(), sequences[i].begin(), sequences[i].end(), back_inserter(intersectionality));
			if (intersectionality.size() >0) ++geno;
		}
		return geno;
	}

/*	vector<int> get_alleles() {
		vector<int> a = sequences[0];
		a.insert(a.end(), sequences[1].begin(), sequences[1].end());
		return(a);
	}
*/

	// generation 0 constructor
	Individual (vector<vector<int>> seqs): sequences(seqs) {
		if (polygenicTrait)
			calculate_trait_value();
	}

	// intra-simulation constructor
	Individual (Individual *p1, Individual *p2, vector<vector<int> > mutation_results, vector<int> breakpoints) {
		if (polygenicTrait) {
			sequences = mutation_results;
			calculate_trait_value();
		} else {
			sequences.push_back((*p1).get_sequence(mutation_results[2][0]));
			sequences.push_back((*p2).get_sequence(mutation_results[3][0]));
			for (int i=0; i<2; ++i) {
				for (int j=1; j<mutation_results[i].size(); ++j)  {
					if (mutation_results[i][j] > 0) {
						sequences[i].push_back(mutation_results[i][j]);
						sort(sequences[i].begin(), sequences[i].end());
					} else
						remove_allele_by_position(i, -1 * mutation_results[i][j]);
				}
			}
			if (breakpoints.size() > 0) {
				resolve_crossover(breakpoints);
			}
		}
	}

	~Individual() {} // destructor

	static mt19937 n;
};

#endif
