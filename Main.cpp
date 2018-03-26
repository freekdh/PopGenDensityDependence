/*=============================================================================================================
Main.cpp
===============================================================================================================

Population genetics growth models

Ailene McPearson and Freek de Haas

Program version
xx/xx/xxxx	:

=============================================================================================================*/

/// Serial founding events for migration.
/// Deme x time x 'mean'-frequency

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

const int BITSETACCURACY = 100000;
typedef boost::dynamic_bitset<> dynamicbitset;
dynamicbitset r_globalbit, m_globalbit;
std::vector<rnd::discrete_distribution>::iterator itdist;

class Individual {
	public:
	Individual(const int &NLOCI, const int &type) : haplotype(NLOCI, type) {}

	Individual(const Individual &Parent1, const Individual &Parent2) {
		
		const int NLOCI = Parent1.haplotype.size();
		assert(NLOCI == Parent2.haplotype.size());
		
		dynamicbitset m_localbit, r_localbit;
		getbitsets(r_localbit,m_localbit,NLOCI);
		assert(r_localbit.size() == Parent1.haplotype.size());

		dynamicbitset temp1 = Parent1.haplotype & r_localbit;
		dynamicbitset temp2 = Parent2.haplotype & r_localbit.flip();
		haplotype = temp1 | temp2;

		temp1 = haplotype & m_localbit;
		temp2 = haplotype | m_localbit;

		haplotype = temp2 & temp1.flip();
	}

	inline int haplotypeint() {
		int type = haplotype.to_ulong();
		return type;
	}

	inline bool returnlocus(int locus) {
		return haplotype[locus];
	}

	private:
	boost::dynamic_bitset<> haplotype;

	void getbitsets(dynamicbitset &r_localbit, dynamicbitset &m_localbit, const int &NLOCI) {
		
		r_localbit.resize(NLOCI);
		m_localbit.resize(NLOCI);

		int rand1, rand2;
		do { rand1 = rnd::integer(BITSETACCURACY); } while (rand1 > BITSETACCURACY - NLOCI);
		do { rand2 = rnd::integer(BITSETACCURACY); } while (rand1 > BITSETACCURACY - NLOCI);

		for (int i = 0; i < NLOCI; ++i) {
			r_localbit[i] = r_globalbit[rand1 + i];
			m_localbit[i] = m_globalbit[rand2 + i];
		}

		assert(r_localbit.size() == NLOCI);
		assert(m_localbit.size() == NLOCI);
	}
};

typedef std::vector<Individual*> Population;
Population::iterator it;

struct Parameters {
	//default parameters
	double RECOMBINATIONRATE = 0.5;
	double MUTATIONRATE = 0.000001;
	int NGEN = 100;
	int NLOCI = 1;
	int NREP = 1;
	int NMETA = 1;
	std::vector<double> r;
	std::vector<double> k;
	std::vector<int> NINIT;
	std::vector<rnd::discrete_distribution*> Migrationdist;
};

struct DataSet {
	DataSet(Parameters* parspointer) : pars(parspointer) {
		// [NGEN][LOCUS][REP]
		ploc.resize(NGEN, std::vector<std::vector<double>>(NLOCI, std::vector<double>()));
		MetaOut.resize(NGEN, std::vector<std::vector<int>>(NLOCI, std::vector<int>(PACCURACY, 0)));
		SubOut.resize(NROW*NCOL, std::vector<std::vector<std::vector<int>>>(NGEN, std::vector<std::vector<int>>(NLOCI, std::vector<int>(PACCURACY, 0))));
		plocsub.resize(NROW*NCOL, std::vector<std::vector<std::vector<double>>>(NGEN, std::vector<std::vector<double>>(NLOCI, std::vector<double>())));
		MetaOutTotal.resize(NGEN, std::vector<double>(NLOCI, 0.0));
	}

	void AddMetaData(std::vector<Population> Metapopulation, const int &GEN) {

		int sumpopsize = 0;
		std::vector<int> nloc(NLOCI);	// [LOCUS]

		for (int i = 0; i < pars.NMETA; ++i) {
			sumpopsize += Metapopulation[i].size();
			for (it = Metapopulation[i].begin(); it < Metapopulation[i].end(); ++it) {
				for (int j = 0; j < pars.NLOCI; ++j) {
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

	void AddSubData(std::vector<Population> Metapopulation, const int &GEN) {
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
	const Parameters* pars;

	std::vector<std::vector<std::vector<double>>> ploc;	// [NGEN][LOCUS][REP]
	std::vector<std::vector<std::vector<int>>> MetaOut; // [NGEN][LOCUS][PACCURACY]
	std::vector<std::vector<double>> MetaOutTotal;

	std::vector<std::vector<std::vector<std::vector<double>>>> plocsub;	// [Sub][NGEN][LOCUS][REP]
	std::vector<std::vector<std::vector<std::vector<int>>>>  SubOut;	//[Sub][NGEN][LOCUS][PACCURACY]

};

void InitializeRandomBitsets(const double &MUTATIONRATE, const double &RECOMBINATIONRATE) {
	r_globalbit.resize(BITSETACCURACY);
	m_globalbit.resize(BITSETACCURACY);

	for (int i = 0; i < BITSETACCURACY; ++i) {
		m_globalbit[i] = rnd::uniform() < MUTATIONRATE ? true : false;
	}

	bool focal = rnd::bernoulli();
	for (int i = 0; i < BITSETACCURACY; ++i) {
		if (rnd::uniform() < RECOMBINATIONRATE) { if(focal == true) {focal = false;} else {focal = true;};}
		r_globalbit[i] = focal;
	}

}

void Mating(std::vector<Population> &Metapopulation, const Parameters &pars) {
	assert(Metapopulation.size() == pars.NMETA);
	std::vector<Population> MetapopulationAfter(pars.NMETA);

	for (int i = 0; i < pars.NMETA; ++i) {
		// Mating: determin nr of offspring -> randomly assign partner per offspring
		Population FocalPopulation = Metapopulation[i];
		const int FocalPopulationSize = FocalPopulation.size();
		for (it = FocalPopulation.begin(); it != FocalPopulation.end(); ++it) {
			int NOffspring = rnd::poisson((1.0 + (pars.r[(*it)->haplotypeint()] * (1.0 - (double)FocalPopulationSize / pars.k[(*it)->haplotypeint()]))));
			for (int j = 0; j < NOffspring; ++j) {
				MetapopulationAfter[i].push_back(new Individual(**it, *FocalPopulation[rnd::integer(FocalPopulationSize)]));
			}
		}
		/*clear memory*/
		for (it = FocalPopulation.begin(); it != FocalPopulation.end(); ++it) { delete *it; }
	}
	
	Metapopulation = MetapopulationAfter;
} 

void Migration(std::vector<Population> &Metapopulation, const Parameters &pars) {
	assert(Metapopulation.size() == pars.NMETA);
	std::vector<Population> MetapopulationAfter(pars.NMETA);

	for (int i = 0; i < Metapopulation.size(); ++i) {
			for (it = Metapopulation[i].begin(); it < Metapopulation[i].end(); ++it) {
				int next = (*pars.Migrationdist[i]).sample();
				assert(0 <= next <= pars.NMETA);
				MetapopulationAfter[next].push_back(*it);
			}
		}

	Metapopulation = MetapopulationAfter;
	assert(Metapopulation.size() == pars.NMETA);
}

void RunSimulation(const Parameters &pars, DataSet &data) {

	Population InitialPopulation;
	for (int type = 0; type < pars.NINIT.size(); ++type) {
		for(int n = 0; n < pars.NINIT[type]; ++n){
			InitialPopulation.push_back(new Individual(pars.NLOCI, type));
		}
	}

	std::vector<Population> FocalMetapopulation(pars.NMETA);
	FocalMetapopulation[0] = InitialPopulation;

	for (int i = 0; i < pars.NGEN; ++i) {
		data.AddMetaData(FocalMetapopulation, i);
		data.AddSubData(FocalMetapopulation, i);

		Migration(FocalMetapopulation, pars.Migrationdist);
		Mating(FocalMetapopulation, pars);
	}

	for (int i = 0; i < FocalMetapopulation.size(); ++i) {
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

void OutputParameters(std::ofstream &ofstream, const Parameters &pars) {
	ofstream.fill(',');
	ofstream << "r"					<< ofstream.fill() << pars.r << std::endl;
	for (int i = 0; i < pars.r.size(); ++i) {
		ofstream << "r" << i		<< ofstream.fill() << pars.r[i] << std::endl;
	}
	for (int i = 0; i < pars.k.size(); ++i) {
		ofstream << "k" << i		<< ofstream.fill() << pars.k[i] << std::endl;
	}
	for (int i = 0; i < pars.NINIT.size(); ++i) {
		ofstream << "NINIT" << i	<< ofstream.fill() << pars.NINIT[i] << std::endl;
	}
	ofstream << "NGEN"				<< ofstream.fill() << pars.NGEN << std::endl;
	ofstream << "NLOCI"				<< ofstream.fill() << pars.NLOCI << std::endl;
	ofstream << "NMETA"				<< ofstream.fill() << pars.NMETA << std::endl;

	ofstream << std::endl << std::endl;
	ofstream << "MigrationMatrix" << std::endl << std::endl;

	for (unsigned int i = 0; i < pars.Migrationdist.size(); ++i) {
		for (int j = 0; j < pars.Migrationdist[i]->size(); ++j) {
			ofstream << (*pars.Migrationdist[i])[j] << ofstream.fill();
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

	/* PARAMETERS HAPLOID SIMULATION */
	Parameters pars;
	
	pars.RECOMBINATIONRATE = 0.5, pars.MUTATIONRATE = 0.00001;
	pars.NGEN = 30, pars.NLOCI = 1, pars.NREP = 50, pars.NMETA = 1;
	assert(pars.NMETA > 0); assert(pars.NREP > 0); assert(pars.NLOCI > 0); assert(pars.NGEN > 0);
	
	// growth: n1[t+1] = n1[t] + n1[t] z[1] (1 - (n1[t]+n2[t])/k) 
	pars.r.push_back(0.1); pars.r.push_back(0.1);
	pars.k.push_back(100.0); pars.k.push_back(100.0);
	pars.NINIT.push_back(10); pars.NINIT.push_back(10);			
	assert(pars.r.size() == pow(2,pars.NLOCI));
	assert(pars.NINIT.size() == pow(2,pars.NLOCI));
	assert(pars.k.size() == pow(2,pars.NLOCI));

	std::vector<rnd::discrete_distribution*> Migrationdist;
	for (int i = 0; i < pars.NMETA; ++i) {
		pars.Migrationdist.push_back(new rnd::discrete_distribution(pars.NMETA));
		for (int j = 0; j < pars.NMETA; ++j) {
			(*pars.Migrationdist[i])[j] = 1.0;
		}
	}

	/* OFSTREAM */
	std::ofstream ParametersOfstream;
	std::vector<std::ofstream> MetaOfstream(pars.NLOCI);
	std::vector<std::vector<std::ofstream>> SubOfstream(pars.NMETA);
	for (int j = 0; j < pars.NMETA; ++j) {
		for (int i = 0; i < pars.NLOCI; ++i) {
			SubOfstream[j].push_back(std::ofstream());
		}
	}
	CreateOutputStreams(ParametersOfstream, MetaOfstream, SubOfstream);

	/* INITIALIZE */
	rnd::set_seed();
	InitializeRandomBitsets(pars.MUTATIONRATE, pars.RECOMBINATIONRATE);
	OutputParameters(ParametersOfstream, pars);
	DataSet data(pars.NLOCI,pars.NGEN);

	for (int i = 0; i < pars.NREP; ++i) {
		std::cout << "Replicate: " << i << std::endl;
		RunSimulation(pars, data);	// input, output
	}

	/* ANALYSIS */
	data.Analysis(MetaOfstream, SubOfstream);
	
	return 0;
}