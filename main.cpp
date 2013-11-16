#include <iostream>
#include <vector>
#include <time.h>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "mtrandom.h"

using namespace std;

struct Individual
{
	int number;
	double fitness;
	vector<int> chromosome;
	vector<double> solution;
};

void display();
void readFile(int& populationSize, int& fitnessFunction);
void initialize(vector<Individual>& population);
pair<Individual,Individual> crossover();
Individual tournamentSelection();
void mutate(pair<Individual,Individual>& children);
void replacement(pair<Individual,Individual> children);
Individual weakestLink();
double evaluate();
double f1(const vector<double>& xs);
double f2(const vector<double>& xs);
vector<int> to_binary(double x, const pair<double,double>& prange, unsigned int num_bits, bool is_gray_coded);
double from_binary(const vector<int>& bits, const pair<double,double>& prange, bool is_gray_coded);

int populationSize;
vector<Individual> population;
int fitnessFunction;
int chromosomeSize = 10;

mtrandom rnd;

pair<double,double> range;
pair<Individual,Individual> children;

bool done = false;
int bestI;
int weakestI;

int main()
{
	cout << "Evolutionary computing!" << endl;

	rnd.seed_random(time(NULL));
	readFile(populationSize, fitnessFunction);
	initialize(population);
	//display();

	while(!done)
	{
		children = crossover();
		mutate(children);
		//TODO: Select offsprings to survive
		replacement(children);
//		cout << setprecision(10) << "Best solution: " << evaluate() << "\r";
		cout << std::fixed;
		cout << setprecision(10) << "Best solution: " << evaluate() << "\r";
		//cout << setprecision(9) << f << '\n';
		if(evaluate() == 0)
			break;
	}
	cout << endl;

    return 0;
}

void display()
{
	for(int i=0; i<populationSize; i++)
	{
		cout << "Individual #" << population[i].number << " -> ";
		cout << "fitness = " << population[i].fitness << " ";
		for(unsigned int j=0; j<population[i].solution.size(); j++)
		{
			cout << "solution : " << population[i].solution[j] << endl;
		}
	}
}

void readFile(int& populationSize, int& fitnessFunction)
{
	fstream myfile("data", ios_base::in);

    if(myfile >> populationSize >> fitnessFunction)
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
//	cout << populationSize << "\t" << fitnessFunction << "\t" << c << "\t" << d << "\t" << e << "\t" << f << endl;

    getchar();
}

void initialize(vector<Individual>& population)
{
	double best;
	double weakest;

	for(int i=0; i<populationSize; i++)	// initialize population
	{
		double tmprnd = rnd.random(range.first, range.second);

		Individual individual;
		individual.solution.push_back(tmprnd);
		individual.fitness = f1(individual.solution);

		if(i == 0)
		{
			weakestI = i;
//			weakest = individual.fitness;
			weakest = best = individual.fitness;
			bestI = i;
		}

		else if(individual.fitness > weakest)
		{
			weakest = individual.fitness;
			weakestI = i;
		}

		else if(individual.fitness < best)
		{
			best = individual.fitness;
			bestI = i;
		}

		individual.number = i;
		individual.chromosome = to_binary(tmprnd, range, chromosomeSize, false);
		population.push_back(individual);
/*
		cout << "Best fitness : \t\t" << best << endl;
		cout << "Best index : \t\t" << bestI  << endl;
		cout << "Weakest fitness : \t" << weakest << endl;
		cout << "Weakest index : \t" << weakestI << endl;
*/
	}
}

pair<Individual,Individual> crossover()
{
	pair<Individual, Individual> children;
	Individual parent1 = tournamentSelection();
	Individual parent2 = tournamentSelection();
	while(parent1.number == parent2.number)
		parent2 = tournamentSelection();

	Individual child1, child2;
	for(int i=0; i<chromosomeSize/2; i++)
	{
		child1.chromosome.push_back(parent1.chromosome[i]);
		child2.chromosome.push_back(parent2.chromosome[i]);
	}
	for(int i=chromosomeSize/2; i<10; i++)
	{
		child1.chromosome.push_back(parent2.chromosome[i]);
		child2.chromosome.push_back(parent1.chromosome[i]);
	}
/*
	cout << "Child 1 binary: " << endl;
	for(int i=0; i<chromosomeSize; i++)
	{
		cout << child1.chromosome[i] << " ";
	}
	cout << endl;

	cout << "Child 2 binary: " << endl;
	for(int i=0; i<chromosomeSize; i++)
	{
		cout << child2.chromosome[i] << " ";
	}
	cout << endl;
*/
	children = make_pair(child1, child2);

	return children;
}

Individual tournamentSelection()
{
	int rand1 = rnd.random(populationSize);
	int rand2 = rnd.random(populationSize);

	while(rand1 == rand2)
		rand2 = rnd.random(populationSize);

	Individual candidate1 = population[rand1];
	Individual candidate2 = population[rand2];

	if(candidate1.fitness < candidate2.fitness) return candidate1;

	else return candidate2;
}

void mutate(pair<Individual,Individual>& children)
{
	int flipBit = rnd.random(chromosomeSize);			// select random "gene" to mutate
//	cout << "bit to flip: " << flipBit << endl;

	if(children.first.chromosome[flipBit] == 0) children.first.chromosome[flipBit] = 1;		// mutate first child
	else children.first.chromosome[flipBit] = 0;

	flipBit = rnd.random(chromosomeSize);				// select new random "gene" to mutate for child 2

	if(children.second.chromosome[flipBit] == 0) children.second.chromosome[flipBit] = 1;		// mutate second child
	else children.second.chromosome[flipBit] = 0;

	children.first.solution.push_back(from_binary(children.first.chromosome, range, false));
	children.second.solution.push_back(from_binary(children.second.chromosome, range, false));

	children.first.fitness = f1(children.first.solution);
	children.second.fitness = f1(children.second.solution);
/*
	cout << "child 1 solution: " << endl;
	for(unsigned int i=0; i<children.first.solution.size(); i++)
		cout << children.first.solution[i];
	cout << endl;

	cout << "child 2 solution: " << endl;
	for(unsigned int i=0; i<children.second.solution.size(); i++)
		cout << children.second.solution[i];
	cout << endl;
*/
}

Individual weakestLink()
{
	double best;
	double weakest;

	for(int i=0; i<populationSize; i++)
	{
		if(i == 0)
		{
			weakestI = i;
			weakest = best = population[i].fitness;
			bestI = i;
		}

		else if(population[i].fitness > weakest)
		{
			weakest = population[i].fitness;
			weakestI = i;
		}

		else if(population[i].fitness < best)
		{
			best = population[i].fitness;
			bestI = i;
		}
	}

	return population[weakestI];
}

void replacement(pair<Individual,Individual> children)
{
	short child;
	Individual weakest = weakestLink();
/*
	cout << "Before replacement: " << endl;
	display();

	cout << "Weakest link : " << weakest.number << " " << weakest.fitness << endl;
	cout << "Child 1 fitness : " << children.first.fitness << endl;
	cout << "Child 2 fitness : " << children.second.fitness << endl;
*/
	if(children.first.fitness < children.second.fitness)	// Which child has lower fitness
	{
		child = 0;
		if(weakest.fitness > children.first.fitness)
		{
//			cout << "Replace weakest with child 1" << endl;
			children.first.number = population[weakestI].number;
			population[weakestI] = children.first;
		}
	}
	else
	{
		child = 1;
		if(weakest.fitness > children.second.fitness)
		{
//			cout << "Replace weakest with child 2" << endl;
			children.second.number = population[weakestI].number;
			population[weakestI] = children.second;
		}
	}

	weakest = weakestLink();

	if(child == 0)
	{
		if(children.second.fitness < weakest.fitness)
		{
			children.second.number = population[weakestI].number;
			population[weakestI] = children.second;
		}
	}
	else
	{
		if(children.first.fitness < weakest.fitness)
		{
			children.first.number = population[weakestI].number;
			population[weakestI] = children.first;
		}
	}
/*
	cout << "After replacement: " << endl;
	display();
	cout << "Weakest link : " << weakest.number << " " << weakest.fitness << endl;
*/
}

double evaluate()
{
	double best;

	for(int i=0; i<populationSize; i++)
	{
		if(i == 0)
			best = population[i].fitness;
		else if(population[i].fitness < best)
			best = population[i].fitness;
	}

	return best;
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
