#include <iostream>
#include <vector>
#include <time.h>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include "mtrandom.h"

using namespace std;

struct Individual
{
	int number;
	double fitness;
	double dGenes;
	vector<int> bGenes;
	vector<double> solution;
};

void readFile(int& populationSize, int& fitnessFunction, int& c, int& d, int& e, int& f);
void initialisePop(vector<Individual>& population);
double fitnessCheck(int options, vector<double> population2);
pair<Individual, Individual> crossover();
Individual tournamentSelection();
double f1(const vector<double>& xs);
double f2(const vector<double>& xs);
vector<int> to_binary(double x, const pair<double,double>& prange, unsigned int num_bits, bool is_gray_coded);
double from_binary(const vector<int>& bits, const pair<double,double>& prange, bool is_gray_coded);

const int _f1_ = 1;
const int _f2_ = 2;
int populationSize;

vector<Individual> population;
int c, d, e, f;
int fitnessFunction;

mtrandom rnd;

pair<double,double> range;
pair<Individual, Individual> children;

int main()
{
	cout << "Evolutionary computing!" << endl;

	rnd.seed_random(time(NULL));

	readFile(populationSize, fitnessFunction, c, d, e, f);
	initialisePop(population);
	children = crossover();

	for(int i=0; i<populationSize; i++)
	{
		cout << "Individual number: " << population[i].number << " -> ";
		cout << "fitness = " << population[i].fitness << " ";
		for(unsigned int j=0; j<population[i].solution.size(); j++)
		{
			cout << "solution : " << population[i].solution[j] << endl;
		}
	}

	//TODO: Select individuals for mating
//	cout << "fitness: " << fitnessCheck(fitnessFunction, population) << endl;

	//TODO: Mutate offspring

	//TODO: Select offsprings to survive

//	double f1pop = f1(population);

//	cout << "f1: " << f1pop << endl;

    return 0;
}

void readFile(int& populationSize, int& fitnessFunction, int& c, int& d, int& e, int& f)
{
	std::fstream myfile("data", std::ios_base::in);

    if(myfile >> populationSize >> fitnessFunction >> c >> d >> e >> f)
    {
		switch(fitnessFunction)
		{
			case 1:
				range = make_pair(-5.12 , 5.11);
				break;
			default:
				break;
		}
    }
    else cout << "Couldn't read from file" << endl;

	cout << "Population size: " << populationSize << endl;
    cout << populationSize << "\t" << fitnessFunction << "\t" << c << "\t" << d << "\t" << e << "\t" << f << endl;

    getchar();
}

void initialisePop(vector<Individual>& population)
{
	for(int i=0; i<populationSize; i++)	// initialize population
	{
		double rndtmp = rnd.random();
		Individual individual;
		individual.solution.push_back(rndtmp);
		individual.fitness = f1(individual.solution);
		individual.number = i;
		individual.bGenes = to_binary(rndtmp, range, 10, false);
		population.push_back(individual);
	}
}

double fitnessCheck(int fitness, vector<double> population2)
{
	switch(fitness)
	{
		case _f1_:
			return f1(population2);
			break;
		default:
			return f1(population2);
	}
}

pair<Individual, Individual> crossover()
{
	pair<Individual, Individual> children;
	Individual parent1 = tournamentSelection();
	Individual parent2 = tournamentSelection();
	while(parent1.number == parent2.number)
		parent2 = tournamentSelection();

	cout << "pair first " << range.first << endl;
	cout << "pair second " << range.second << endl;

	cout << "Parent 1 : " << parent1.number << endl;
	cout << "Parent 2 : " << parent2.number << endl;

	cout << "parent 1 binary: " << endl;
	for(unsigned int i=0; i<parent1.bGenes.size(); i++)
	{
		cout << parent1.bGenes[i] << " ";
	}

	cout << "\nparent 2 binary: " << endl;
	for(unsigned int i=0; i<parent2.bGenes.size(); i++)
	{
		cout << parent2.bGenes[i] << " ";
	}

	cout << endl;

	Individual child1, child2;
	for(unsigned int i=0; i<parent1.bGenes.size()/2; i++)
	{
		child1.bGenes.push_back(parent1.bGenes[i]);
		child2.bGenes.push_back(parent2.bGenes[i]);
	}
	for(unsigned int i=parent2.bGenes.size()/2; i<10; i++)
	{
		child1.bGenes.push_back(parent2.bGenes[i]);
		child2.bGenes.push_back(parent1.bGenes[i]);
	}

	child1.solution.push_back(from_binary(child1.bGenes, range, false));
	child2.solution.push_back(from_binary(child2.bGenes, range, false));

	cout << endl;
	cout << "child 1 solution: " << endl;
	for(unsigned int i=0; i<child1.solution.size(); i++)
		cout << child1.solution[i];

	cout << endl;

	cout << "child 2 solution: " << endl;
	for(unsigned int i=0; i<child2.solution.size(); i++)
		cout << child2.solution[i];

	cout << endl;

	cout << "Child 1 binary: " << endl;
	for(unsigned int i=0; i<child1.bGenes.size(); i++)
	{
		cout << child1.bGenes[i] << " ";
	}
	cout << endl;

	cout << "Child 2 binary: " << endl;
	for(unsigned int i=0; i<child2.bGenes.size(); i++)
	{
		cout << child2.bGenes[i] << " ";
	}
	cout << endl;

	children = make_pair(child1, child2);

	return children;
}

Individual tournamentSelection()
{
	int rand1 = rnd.random(populationSize);
	int rand2 = rnd.random(populationSize);

	while(rand1 == rand2)
		rand2 = rnd.random(populationSize);

	cout << "rand1 : " << rand1 << endl;
	cout << "rand2 : " << rand2 << endl;

	Individual candidate1 = population[rand1];
	Individual candidate2 = population[rand2];

	cout << "candidate 1 fitness : " << candidate1.fitness << endl;
	cout << "candidate 2 fitness : " << candidate2.fitness << endl;

	if(candidate1.fitness < candidate2.fitness) return candidate1;

	else return candidate2;
}

// Dejong’s F1 (sphere) function N 2
// Minimize f (x̄)
// N = any positive integer
// −5.12 ≤ xi ≤ 5.11
// Optimum is f (x) = 0 at xi = 0
double f1(const vector<double>& xs)
{
	double sum = 0.0;

	for(unsigned int i=0; i<xs.size(); ++i)
	{
		sum += xs[i] * xs[i];
	}

	return sum;
}

// Dejong’s F2 function
// Minimize f ( ̄) = 100(x2 − x2 )2 + (1.0 − x1 )2 .x1
// N = 2
// −2.048 ≤ xi ≤ 2.047
// Optimum is f (x) = 0 at (1, 1)
double f2(const vector<double>& xs)
{
	double term1 = (xs[0] * xs[0]) - xs[1];
	double term2 = 1.0 - xs[0];

	return 100 * term1 * term1 + term2 * term2;
}

vector<int> to_binary(double x, const pair<double,double>& prange, unsigned int num_bits, bool is_gray_coded)
{
	// map the double onto the appropriate integer range
	double range = prange.second - prange.first;
	double max_bit_val = pow(2.0, static_cast<double>(num_bits)) - 1;
	int int_val = static_cast<int>((x - prange.first) * max_bit_val / range + 0.5);

	// convert the integer to binary
	std::vector<int> result(num_bits);
	for(unsigned int b=0; b<num_bits; ++b)
	{
		result[b] = int_val % 2;
		int_val/=2;
	}

	if(is_gray_coded)
	{
		for(unsigned int b=0; b<num_bits-1; ++b)
		{
			result[b] = !(result[b] == result[b+1]);
		}
	}

	return result;
}

double from_binary(const std::vector<int>& bits, const std::pair<double,double>& prange, bool is_gray_coded)
{
	unsigned int num_bits = bits.size();
	unsigned int intval = 0;

	if(is_gray_coded)
	{
	// convert from gray to binary
		std::vector<int> binary(num_bits);
		binary[num_bits-1] = bits[num_bits-1];
		intval = binary[num_bits-1];

		for(int i=num_bits-2; i>=0; i--)
		{
			binary[i] = !(binary[i+1] == bits[i]);
			intval += intval + binary[i];
		}
	}
	else
	{
		// convert from binary encoding to integer
		for(int i=num_bits-1; i>=0; i--)
			intval += intval + bits[i];
}
	// convert from integer to double in the appropriate range
	double range = prange.second - prange.first;
	double m = range / (pow(2.0, double(num_bits)) - 1.0);

	return m * double(intval) + prange.first;
}
