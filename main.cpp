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
	double fitness;
	double dGenes;
	vector<int> bGenes;
	vector<double> solution;
};

void readFile(int& populationSize, int& fitnessFunction, int& c, int& d, int& e, int& f);
void initialisePop(vector<Individual>& population);
double fitnessCheck(int options, vector<double> population2);
double f1(const vector<double>& xs);
double f2(const vector<double>& xs);
vector<int> to_binary(double x, const pair<double,double>& prange, unsigned int num_bits, bool is_gray_coded);

const int _f1_ = 1;
const int _f2_ = 2;
int populationSize;

vector<Individual> population;
int c, d, e, f;
int fitnessFunction;

mtrandom test;

int main()
{
	cout << "Evolutionary computing!" << endl;


	readFile(populationSize, fitnessFunction, c, d, e, f);
	initialisePop(population);

	for(int i=0; i<populationSize; i++)
	{
		cout << "Individual number: " << i << " : ";
		cout << " fitness = " << population[i].fitness << " ";
		for(int j=0; j<1; j++)
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

    if(myfile >> populationSize >> fitnessFunction >> c >> d >> e >> f){}
    else cout << "Couldn't read from file" << endl;

	cout << "Population size: " << populationSize << endl;
    cout << populationSize << "\t" << fitnessFunction << "\t" << c << "\t" << d << "\t" << e << "\t" << f << endl;

    getchar();
}

void initialisePop(vector<Individual>& population)
{
	for(int i=0; i<populationSize; i++)	// initialize population
	{
		double rndtmp = test.random();
		Individual individual;
		individual.solution.push_back(rndtmp);
		individual.fitness = f1(individual.solution);
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
