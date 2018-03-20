/*=============================================================================================================
Main.cpp
===============================================================================================================

Population genetics growth models

Ailene McPearson and Freek de Haas

Program version
xx/xx/xxxx	:

=============================================================================================================*/

#include "random.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <bitset>
#include <math.h>
#include <assert.h>
#include <iomanip>
#include <ctime>

enum type { type1, type2 };
enum fitness { logistic_r, logistic_k };
enum direction {stay, left, right, up, down };

const int NROW = 3;
const int NCOL = 3;
const int PACCURACY = 100;
const int BITSETACCURACY = 1000;

typedef boost::dynamic_bitset<> dynamicbitset;

class Individual {
public:
	Individual(const int &NLOCI, const int &type) : haplotype(NLOCI, type) { ; }

	Individual(const Individual &Parent1, const Individual &Parent2, dynamicbitset r_localbit, dynamicbitset m_localbit) {

		dynamicbitset temp1 = Parent1.haplotype & r_localbit;
		dynamicbitset temp2 = Parent2.haplotype & r_localbit.flip();
		haplotype = temp1 | temp2;

		temp1 = haplotype & m_localbit;
		temp2 = haplotype | m_localbit;

		haplotype = temp2 & temp1.flip();
	}

	inline int haplotypeint() {
		return haplotype.to_ulong();
	}

	inline bool returnlocus(int locus) {
		return haplotype[locus];
	}

private:
	boost::dynamic_bitset<> haplotype;
};

typedef std::vector<Individual*> subpopulation;

subpopulation::iterator it;
std::vector<rnd::discrete_distribution>::iterator itdist;

struct DataSet {
	DataSet(const int &NLOCI, const int &NGEN) : NLOCI(NLOCI), NGEN(NGEN) {
		// [NGEN][LOCUS][REP]
		ploc.resize(NGEN, std::vector<std::vector<double>>(NLOCI, std::vector<double>()));
		MetaOut.resize(NGEN, std::vector<std::vector<int>>(NLOCI, std::vector<int>(PACCURACY, 0)));
		SubOut.resize(NROW*NCOL, std::vector<std::vector<std::vector<int>>>(NGEN, std::vector<std::vector<int>>(NLOCI, std::vector<int>(PACCURACY, 0))));
		plocsub.resize(NROW*NCOL, std::vector<std::vector<std::vector<double>>>(NGEN, std::vector<std::vector<double>>(NLOCI, std::vector<double>())));
		MetaOutTotal.resize(NGEN, std::vector<double>(NLOCI, 0.0));
	}

	void AddMetaData(subpopulation Metapopulation[], const int &GEN) {

		int sumpopsize = 0;
		std::vector<int> nloc(NLOCI);	// [LOCUS]

		for (int i = 0; i < NROW * NCOL; ++i) {
			sumpopsize += Metapopulation[i].size();
			for (it = Metapopulation[i].begin(); it < Metapopulation[i].end(); ++it) {
				for (int j = 0; j < NLOCI; ++j) {
					nloc[j] += (int)(*it)->returnlocus(j);
				}
			}
		}

		if (sumpopsize != 0) {
			// Calculate p
			for (int loc = 0; loc < NLOCI; ++loc) {
				double temp = (double)nloc[loc] / (double)sumpopsize;
				ploc[GEN][loc].push_back(temp);
			}
		}
	}

	void AddSubData(subpopulation Metapopulation[], const int &GEN) {
		for (int sub = 0; sub < NROW*NCOL; ++sub) {
			int sumpopsize = Metapopulation[sub].size();
			std::vector<int> nloc(NLOCI);

			for (it = Metapopulation[sub].begin(); it < Metapopulation[sub].end(); ++it) {
				for (int j = 0; j < NLOCI; ++j) {
					nloc[j] += (int)(*it)->returnlocus(j);
				}
			}
		
			if (sumpopsize != 0) {
				// Calculate p
				for (int loc = 0; loc < NLOCI; ++loc) {
					double temp = (double)nloc[loc] / (double)sumpopsize;
					plocsub[sub][GEN][loc].push_back(temp);
				}
			}
		}
	}

	void Analysis(std::vector<std::ofstream> &MetaStream, std::vector<std::vector<std::ofstream>> &SubStream) {

		// Metadata

		 // [NGEN][LOCUS][PACCURACY]

		for (int loc = 0; loc < NLOCI; ++loc) {
			for (int gen = 0; gen < NGEN; ++gen) {
				for (unsigned int i = 0; i < ploc[gen][loc].size(); ++i) {
					double intpart, fracpart = modf(ploc[gen][loc][i], &intpart);
					intpart = static_cast<int>(fracpart*PACCURACY);
					++MetaOut[gen][loc][intpart];
					++MetaOutTotal[gen][loc];
				}
			}
		}

		for (int loc = 0; loc < NLOCI; ++loc) {
			for (int gen = 0; gen < NGEN; ++gen) {
				for (int i = 0; i < PACCURACY; ++i) {
					MetaStream[loc] << (double)MetaOut[gen][loc][i] / MetaOutTotal[gen][loc] << MetaStream[loc].fill();
				}
				MetaStream[loc] << std::endl;
			}
		}

		// Subdata
		for (int sub = 0; sub < NROW*NCOL; ++sub) {
			for (int loc = 0; loc < NLOCI; ++loc) {
				for (int gen = 0; gen < NGEN; ++gen) {
					for (unsigned int i = 0; i < plocsub[sub][gen][loc].size(); ++i) {
						double intpart, fracpart = modf(plocsub[sub][gen][loc][i], &intpart);
						intpart = static_cast<int>(fracpart*PACCURACY);
						++SubOut[sub][gen][loc][intpart];
					}
				}
			}
		}

		for (int sub = 0; sub < NROW*NCOL; ++sub) {
			for (int loc = 0; loc < NLOCI; ++loc) {
				for (int gen = 0; gen < NGEN; ++gen) {
					for (int i = 0; i < PACCURACY; ++i) {
						SubStream[sub][loc] << SubOut[sub][gen][loc][i] << SubStream[sub][loc].fill();
					}
					SubStream[sub][loc] << std::endl;
				}
			}
		}

	}

private:
	const int NLOCI;
	const int NGEN;

	std::vector<std::vector<std::vector<double>>> ploc;	// [NGEN][LOCUS][REP]
	std::vector<std::vector<std::vector<int>>> MetaOut; // [NGEN][LOCUS][PACCURACY]
	std::vector<std::vector<double>> MetaOutTotal;

	std::vector<std::vector<std::vector<std::vector<double>>>> plocsub;	// [Sub][NGEN][LOCUS][REP]
	std::vector<std::vector<std::vector<std::vector<int>>>>  SubOut;	//[Sub][NGEN][LOCUS][PACCURACY]
};

void Setrandombitset(dynamicbitset &recombit, dynamicbitset &mutationbit, const double &MUTATIONRATE, const double &RECOMBINATIONRATE) {
	recombit.resize(BITSETACCURACY);
	mutationbit.resize(BITSETACCURACY);

	for (int i = 0; i < BITSETACCURACY; ++i) {
		mutationbit[i] = rnd::uniform() < MUTATIONRATE ? true : false;
	}

	bool focal = rnd::bernoulli();
	for (int i = 0; i < BITSETACCURACY; ++i) {
		if (rnd::uniform() < RECOMBINATIONRATE) { focal = focal == true ? false : true; }
		recombit[i] = focal;
	}

}

inline void PopulationbyType (subpopulation &pop, int types[], const int &NHAPLOTYPES) {
	for (int i = 0; i < NHAPLOTYPES; ++i) {
		types[i] = 0;
	}
	
	for (it = pop.begin(); it != pop.end(); ++it) {
		++types[(*it)->haplotypeint()];
	}
}

void MigrationMatrix(std::vector<rnd::discrete_distribution*> &distarray) {
	
	const int NTOT = NROW * NCOL;	
	for (int i = 0; i < NTOT; ++i) {
		distarray.push_back(new rnd::discrete_distribution(NTOT));
	}
	
	/* distarray */
	for (int i = 0; i < NTOT; ++i) {
		for (int j = 0; j < NTOT; ++j) {
			(*distarray[i])[j] = 1.0;  // this does not make any sense now}
		}
	}
}

void getbitsets(dynamicbitset &r_localbit, dynamicbitset &m_localbit, const dynamicbitset &r_globalbit, const dynamicbitset &m_globalbit, const int &NLOCI) {
	
	r_localbit.resize(NLOCI);
	m_localbit.resize(NLOCI);

	int rand1, rand2;
	do { rand1 = rnd::integer(BITSETACCURACY); } while (rand1 > BITSETACCURACY - NLOCI);
	do { rand2 = rnd::integer(BITSETACCURACY); } while (rand1 > BITSETACCURACY - NLOCI);

	for (int i = 0; i < NLOCI; ++i) {
		r_localbit[i] = r_globalbit[rand1 + i];
		m_localbit[i] = m_globalbit[rand2 + i];
	}
};

void Mating(subpopulation Metapopulation[], const int &NLOCI, const int &NHAPLOTYPES, const dynamicbitset &r_globalbit, const dynamicbitset &m_globalbit, int n[], const double z[], const double &r, const double &k) {
	subpopulation MetapopulationAfter[NROW*NCOL];

	for (int pop = 0; pop < NCOL*NROW; ++pop) {
			subpopulation focalBefore = Metapopulation[pop];
			const int popsize = focalBefore.size();
			PopulationbyType(focalBefore, n, NHAPLOTYPES);
			for (it = focalBefore.begin(); it != focalBefore.end(); ++it) {
				int nind = rnd::poisson((1.0 + (r * z[(*it)->haplotypeint()] * (1.0 - (double)popsize / k))));
				for (int j = 0; j < nind; ++j) {
					dynamicbitset r_localbit, m_localbit;
					getbitsets(r_localbit, m_localbit, r_globalbit, m_globalbit,NLOCI);
					MetapopulationAfter[pop].push_back(new Individual(**it, *focalBefore[rnd::integer(popsize)], r_localbit, m_localbit));
				}
			}
			/*clear memory*/
			for (it = focalBefore.begin(); it != focalBefore.end(); ++it) { delete *it; }
	}
	
	for (int pop = 0; pop < NROW*NCOL; ++pop) {
		Metapopulation[pop] = MetapopulationAfter[pop];
	}
} 

void Migration(subpopulation Metapopulation[], const std::vector<rnd::discrete_distribution*> &distribution ) {
	
	subpopulation MetapopulationAfter[NROW*NCOL];

	for (int i = 0; i < NROW*NCOL; ++i) {
			subpopulation focalsubopulation = Metapopulation[i];
			for (it = focalsubopulation.begin(); it < focalsubopulation.end(); ++it) {
				int next = (*distribution[i]).sample();
				MetapopulationAfter[next].push_back(*it);
			}
		}

	for(int i = 0; i < NROW*NCOL; ++i) {
		Metapopulation[i] = MetapopulationAfter[i];
	}
}

void RunSimulation(subpopulation InitialMetapopulation[], const std::vector<rnd::discrete_distribution*> &MigrationMatrix, const dynamicbitset &r_globalbit, const dynamicbitset &m_globalbit, const int &NGEN, const int &NHAPLOTYPES, const int &NLOCI, int n[], const double z[], const double &r, const double &k, DataSet &Save) {
	assert(NHAPLOTYPES > 0 && NHAPLOTYPES < 10);

	subpopulation FocalMetapopulation[NROW*NCOL];
	for (int i = 0; i < NROW*NCOL; ++i) {
		FocalMetapopulation[i] = InitialMetapopulation[i];
	}

	for (int i = 0; i < NGEN; ++i) {
		Save.AddMetaData(FocalMetapopulation,i);
		Save.AddSubData(FocalMetapopulation, i);

		Migration(FocalMetapopulation, MigrationMatrix);
		Mating(FocalMetapopulation, NLOCI, NHAPLOTYPES, r_globalbit, m_globalbit, n, z, r, k);
	}

	for (int i = 0; i < NROW*NCOL; ++i) {
		for (it = FocalMetapopulation[i].begin(); it != FocalMetapopulation[i].end(); ++it) { delete *it; }
	}
}

std::string ReturnTimeStamp(const std::string &CurrentDirectory) {
	std::string out;

	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
	std::string str = oss.str();
	
	out = CurrentDirectory;
	out.append("/");
	out.append(str);

	return out;
}

void OutputParameters(std::ofstream &ofstream, std::vector<rnd::discrete_distribution*> dist, const double &r, const double &k, const int &NGEN, const int &NLOCI, const int &NHAPLOTYPES, const double z[]) {
	ofstream.fill(',');
	ofstream << "r"					<< ofstream.fill() << r << std::endl;
	ofstream << "k"					<< ofstream.fill() << k << std::endl;
	for (int i = 0; i < NHAPLOTYPES; ++i) {
		ofstream << "z" << i		<< ofstream.fill() << z[i] << std::endl;
	}
	ofstream << "NGEN"				<< ofstream.fill() << NGEN << std::endl;
	ofstream << "NLOCI"				<< ofstream.fill() << NLOCI << std::endl;
	ofstream << "NHAPLOTYPES"		<< ofstream.fill() << NHAPLOTYPES << std::endl;
	ofstream << "NROW"				<< ofstream.fill() << NROW << std::endl;
	ofstream << "NCOL"				<< ofstream.fill() << NCOL << std::endl;

	ofstream << std::endl << std::endl;
	ofstream << "MigrationMatrix" << std::endl << std::endl;

	for (unsigned int i = 0; i < dist.size(); ++i) {
		for (int j = 0; j < dist[i]->size(); ++j) {
			ofstream << (*dist[i])[j] << ofstream.fill();
		}
		ofstream << std::endl;
	}


}

void CreateOutputStreams(std::ofstream &ParametersOfstream, std::vector<std::ofstream> &MetaOfstream, std::vector<std::vector<std::ofstream>> &SubOfStream) {
	std::string CurrentWorkingDirectory = "/home/freek/PopGenDensityDependence/";
	//std::string Current = boost::filesystem::current_path();
	std::string MainOutputFolder = ReturnTimeStamp(CurrentWorkingDirectory);
	std::string MetaOutputFolder = MainOutputFolder;
	std::string SubOutputFolder = MainOutputFolder;
	MetaOutputFolder.append("/MetaPopulation");
	SubOutputFolder.append("/SubPopulations");
	boost::filesystem::create_directories(MainOutputFolder.c_str());
	boost::filesystem::create_directories(MetaOutputFolder.c_str());
	boost::filesystem::create_directories(SubOutputFolder.c_str());

	ParametersOfstream.open(MainOutputFolder.append("/Parameters.csv"));

	// MetaFiles:
	
	for (unsigned int i = 0; i < MetaOfstream.size(); ++i) {
		std::string focallocus = MetaOutputFolder;
		focallocus.append("/Locus");
		focallocus = focallocus + std::to_string(i);
		focallocus = focallocus.append(".csv");
		MetaOfstream[i].open(focallocus);
		assert(MetaOfstream[i].is_open());
		MetaOfstream[i].fill(',');
	}

	// SubFiles:
	for (unsigned int sub = 0; sub < SubOfStream.size(); ++sub) {
		std::string focal = SubOutputFolder;
		focal.append("/Sub");
		focal = focal + std::to_string(sub);
		boost::filesystem::create_directories(focal.c_str());
		for (unsigned loc = 0; loc < SubOfStream[sub].size(); ++loc) {
			std::string temp = focal;
			temp.append("/Locus");
			temp = temp + std::to_string(loc);
			temp = temp.append(".csv");
			SubOfStream[sub][loc].open(temp);
			assert(SubOfStream[sub][loc].is_open());
			SubOfStream[sub][loc].fill(',');
		}
	}
}

int main() {

	/* PARAMETERS */
	const double r = 0.5, k = 100.0, RECOMBINATIONRATE = 0.5, MUTATIONRATE = 0.00001;
	const int NGEN = 30, NLOCI = 2, NHAPLOTYPES = 4, NREP = 200; 
	const double z[NHAPLOTYPES] = { 0.1, 0.1, 0.1, 0.1 };
	std::vector<rnd::discrete_distribution*> Migrationdist;
	MigrationMatrix(Migrationdist);

	/* OFSTREAM */
	std::ofstream ParametersOfstream;
	std::vector<std::ofstream> MetaOfstream(NLOCI);
	std::vector<std::vector<std::ofstream>> SubOfstream(NROW*NCOL);

	for (int j = 0; j < NROW*NCOL; ++j) {
		for (int i = 0; i < NLOCI; ++i) {
			SubOfstream[j].push_back(std::ofstream());
		}
	}

	CreateOutputStreams(ParametersOfstream, MetaOfstream, SubOfstream);
	assert(ParametersOfstream.is_open());

	/* INITIALIZE */
	rnd::set_seed();
	dynamicbitset r_globalbit, m_globalbit;
	Setrandombitset(r_globalbit, m_globalbit, MUTATIONRATE,RECOMBINATIONRATE);
	OutputParameters(ParametersOfstream, Migrationdist, r, k, NGEN, NLOCI, NHAPLOTYPES, z);
	DataSet Save(NLOCI,NGEN);

	for (int i = 0; i < NREP; ++i) {
		subpopulation initialpop;
		for (int i = 0; i < 15; ++i) {
			initialpop.push_back(new Individual(NLOCI, rnd::integer(NHAPLOTYPES)));
		}
		subpopulation initialmetapopulation[NROW*NCOL];
		initialmetapopulation[0] = initialpop;

		std::cout << "Replicate: " << i << std::endl;

		int n[NHAPLOTYPES];
		RunSimulation(initialmetapopulation, Migrationdist, r_globalbit, m_globalbit, NGEN, NHAPLOTYPES, NLOCI, n, z, r, k, Save);
	}

	/* ANALYSIS */
	Save.Analysis(MetaOfstream, SubOfstream);
	
	return 0;
}