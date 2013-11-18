#include <iostream>
#include <vector>
#include <time.h>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include "mtrandom.h"

using namespace std;

int N;
struct Individual
{
	int age;
	int number;
	double fitness;
	vector<int> chromosomes[2];
	vector<double> solution;
};

void display();
void readFile(int& populationSize, int& fitnessFunction);
void readFile(string file);
void initialize(vector<Individual>& population);
pair<Individual, Individual> parentSelection();
Individual randomSelection();
Individual tournamentSelection();
Individual tournamentSelectionWeaker();
pair<Individual,Individual> crossover(pair<Individual,Individual> parents);
pair<Individual,Individual> fixedPointCrossover(pair<Individual,Individual> parents);
pair<Individual,Individual> randomPointCrossover(pair<Individual,Individual> parents);
pair<Individual,Individual> twoPointCrossover(pair<Individual,Individual> parents);
void mutate(pair<Individual,Individual>& children);
void mutate2(pair<Individual,Individual>& children);
void replacement(pair<Individual,Individual> children);
void replaceWithTournament(pair<Individual,Individual> children);
void replaceWeakest(pair<Individual,Individual> children);
void replaceRandom(pair<Individual,Individual> children);
void replaceOldest(pair<Individual,Individual> children);
Individual oldest();
Individual weakestLink();
Individual fittest();
double evaluate();
double fitnessCheck(vector<double> solution);
void updateAge(int cycles);

double f1(const vector<double>& xs);
double f2(const vector<double>& xs);
double f3(const vector<double>& xs);
vector<int> to_binary(double x, const pair<double,double>& prange, unsigned int num_bits, bool is_gray_coded);
double from_binary(const vector<int>& bits, const pair<double,double>& prange, bool is_gray_coded);

int populationSize;
vector<Individual> population;
int fitnessFunction;
int chromosomeSize = 10;

mtrandom rnd;

pair<double,double> range;
pair<Individual,Individual> children;
pair<Individual,Individual> parents;

bool useAge = false;
bool done = false;
int bestI;
int weakestI;
string file;

bool gray;
char replacementStrategy;
int iterations = 0;
double mutationRate;
char selectionStrategy;
int crossoverStrategy;
int maxIterations;

double currentBest;
double lastBest;

int main(int argc, char *argv[])
{
	rnd.seed_random(time(NULL));
	cout << "Evolutionary computing!" << endl;
	if(argc != 2)
	{
		cout << "usage" << endl;
		exit(1);
	}
	else
	{
		file = argv[1];
		readFile(file);
	}

	initialize(population);

//	display();

	do
	{
		parents = parentSelection();
		children = crossover(parents);
		mutate2(children);
		if(useAge) updateAge(iterations);
		replacement(children);
		cout << std::fixed;
		cout << setprecision(10) << "Best solution: " << evaluate() << " x: " << population[bestI].solution[0] << " y:" << population[bestI].solution[1] << " iteration#: " << iterations << "\r";
//		cout << "\n Iteration: " << iterations << "\r";
//		cout << setprecision(10) << "Best solution: " << evaluate() << endl;
		currentBest = evaluate();
		if(iterations == 0)
			lastBest = currentBest;
		if(iterations % 10000 == 0 && iterations > 0)
		{
			cout << "10000 iterations" << endl;
			cout << "currentBest: " << currentBest << endl;
			cout << "lastBest: " << lastBest << endl;
			if(currentBest >= lastBest)
				done = true;
			else lastBest = currentBest;
		}
/*
		if(evaluate() == 0.000000000)
			done = true;
*/
		iterations++;
	}while(!done && iterations < maxIterations);
	cout << endl;
	cout << "currentBest: " << currentBest << endl;
	cout << "lastBest: " << lastBest << endl;

//	display();

    return 0;
}

void display()
{
	for(int i=0; i<populationSize; i++)
	{
		cout << "#" << population[i].number << " -> " << endl;
		cout << "\tfitness = " << population[i].fitness << endl;
		for(unsigned int j=0; j<population[i].solution.size(); j++)
		{
			cout << "\tsolution : #" << j << " " << population[i].solution[j] << endl;
		}
	}
}

void readFile(string file)
{
	fstream myfile(file.c_str(), ios_base::in);

    if(myfile >> populationSize >> fitnessFunction >> gray >> replacementStrategy >> mutationRate >> selectionStrategy >> crossoverStrategy >> maxIterations)
    {
		switch(fitnessFunction)						// problem to solve
		{
			case 1:									// problem 1
				range = make_pair(-5.12 , 5.11);
				N = 1;
				break;
			case 2:									// problem 2
				range = make_pair(-2.048, 2.047);
				N = 2;
				break;
			case 3:
				range = make_pair(-65.536, 65.535);
				N = 2;
				break;
			case 4:
				range = make_pair(-512.0, 511.0);
				N = 10;
			default:
				break;
		}
		switch(replacementStrategy)
		{
			case 'O':
				useAge = true;
				break;
			default:
				break;
		}
    }
    else
    {
		cout << "Couldn't read from file" << endl;
		exit(1);
	}

	cout << "Population size: " << populationSize << endl;
	cout << "problem: " << fitnessFunction << endl;
	cout << "Gray coded: " << gray << endl;
	cout << "Replacement strategy: " << replacementStrategy <<  endl;
	cout << "Mutation rate: " << mutationRate << endl;
	cout << "Selection strategy: " << selectionStrategy << endl;
	cout << "Crossover strategy: " << crossoverStrategy << endl;
	cout << "Max iterations: " << maxIterations << endl;

//    getchar();
}

void initialize(vector<Individual>& population)
{
	double best;
	double weakest;

	for(int i=0; i<populationSize; i++)	// initialize population
	{
		Individual individual;
		individual.number = i;

		for(int j=0; j<N; j++)
		{
			vector<int> chromosome;
			double tmprnd = rnd.random(range.first, range.second);
			individual.chromosomes[j] = chromosome;
			individual.chromosomes[j] = to_binary(tmprnd, range, chromosomeSize, gray);
			individual.solution.push_back(tmprnd);
		}

		individual.fitness = fitnessCheck(individual.solution);

		if(i == 0)								// if it's the first individual
		{
			weakestI = i;
			weakest = individual.fitness;
//			weakest = best = individual.fitness;
			bestI = i;
		}

		/*else if(individual.fitness > weakest)
		{
			weakest = individual.fitness;
			weakestI = i;
		}*/

		else if(individual.fitness < best)
		{
			best = individual.fitness;
			bestI = i;
		}

		population.push_back(individual);
	}
}

pair<Individual, Individual> parentSelection()
{
	Individual parent1, parent2;
	pair<Individual,Individual> parents;

	switch(selectionStrategy)							// decide on what selection function to use
	{
		case 'R':											// selectionStrategy 1 means parents are selected randomly

			parent1 = randomSelection();				// select parent 1 randomly
			parent2 = randomSelection();				// select parent 2 randomly

			while(parent1.number == parent2.number)		// while we have the same parents
				parent2 = randomSelection();			// select a new random parent 2

			parents = make_pair(parent1, parent2);		// make the pair
			return parents;								// return parents pair

		case 'T':											// selectionStrategy 2 means parents are selected with tournament selection

			parent1 = tournamentSelection();			// select parent 1 with tournament selection
			parent2 = tournamentSelection();			// select parent 2 with tournament selection

			while(parent1.number == parent2.number)		// while we have the same parents
				parent2 = tournamentSelection();		// select a new parent 2 with tournament selection

			parents = make_pair(parent1, parent2);		// make the pair
			return parents;								// return parents pair

		default:										// default is tournament selection

			parent1 = tournamentSelection();			// select parent 1 with tournament selection
			parent2 = tournamentSelection();			// select parent 2 with tournament selection

			while(parent1.number == parent2.number)		// while we have the same parents
				parent2 = tournamentSelection();		// select a new parent 2 with tournament selection

			parents = make_pair(parent1, parent2);		// make the pair
			return parents;								// return parents pair
	}
}

Individual randomSelection()
{
	Individual individual;
	individual = population[rnd.random(populationSize)];

	return individual;
}

Individual tournamentSelection()
{
	int rand1 = rnd.random(populationSize);
	int rand2 = rnd.random(populationSize);

	while(rand1 == rand2) rand2 = rnd.random(populationSize);

	Individual candidate1 = population[rand1];
	Individual candidate2 = population[rand2];

	if(candidate1.fitness < candidate2.fitness) return candidate1;

	else return candidate2;
}

Individual tournamentSelectionWeaker()
{
	int rand1 = rnd.random(populationSize);
	int rand2 = rnd.random(populationSize);

	while(rand1 == rand2) rand2 = rnd.random(populationSize);

	Individual candidate1 = population[rand1];
	Individual candidate2 = population[rand2];

	if(candidate1.fitness > candidate2.fitness) return candidate1;

	else return candidate2;
}

pair<Individual,Individual> crossover(pair<Individual,Individual> parents)
{
	switch(crossoverStrategy)
	{
		case 1:
			return fixedPointCrossover(parents);
		case 2:
			return randomPointCrossover(parents);
        case 3:
            return twoPointCrossover(parents);
		default:
			return randomPointCrossover(parents);
	}
}

pair<Individual,Individual> fixedPointCrossover(pair<Individual,Individual> parents)
{
	pair<Individual, Individual> children;
	Individual child1, child2;

	for(int i=0; i<N; i++)
		for(int j=0; j<chromosomeSize/2; j++)
		{
			child1.chromosomes[i].push_back(parents.first.chromosomes[i][j]);
			child2.chromosomes[i].push_back(parents.second.chromosomes[i][j]);
		}

	for(int i=0; i<N; i++)
		for(int j=chromosomeSize/2; j<chromosomeSize; j++)
		{
			child1.chromosomes[i].push_back(parents.second.chromosomes[i][j]);
			child2.chromosomes[i].push_back(parents.first.chromosomes[i][j]);
		}

	children = make_pair(child1, child2);

	return children;
}

pair<Individual,Individual> randomPointCrossover(pair<Individual,Individual> parents)
{
	pair<Individual, Individual> children;
	Individual child1, child2;

	for(int i=0; i<N; i++)
	{
		int point = rnd.random(0, chromosomeSize);
		for(int j=0; j<point; j++)
		{
			child1.chromosomes[i].push_back(parents.first.chromosomes[i][j]);
			child2.chromosomes[i].push_back(parents.second.chromosomes[i][j]);
		}
		for(int j=point; j<chromosomeSize; j++)
		{
			child1.chromosomes[i].push_back(parents.first.chromosomes[i][j]);
			child2.chromosomes[i].push_back(parents.second.chromosomes[i][j]);
		}
	}

	children = make_pair(child1, child2);

	return children;
}

pair<Individual,Individual> twoPointCrossover(pair<Individual,Individual> parents)
{
    pair<Individual, Individual> children;
	Individual child1, child2;

	for(int i=0; i<N; i++)
	{
		int point1 = rnd.random(0, chromosomeSize);
		int point2;
		do
        {
            point2 = rnd.random(0, chromosomeSize);
        }while(point2 == point1);

        if(point1 > point2)
        {
            int temp = point1;
            point1 = point2;
            point2 = temp;
        }
        //cout << "\nP1: " << point1 << " P2: " << point2 << endl;
        if(point1 < point2)
        {
            for(int j=0; j<point1; j++)
            {
                child1.chromosomes[i].push_back(parents.first.chromosomes[i][j]);
                child2.chromosomes[i].push_back(parents.second.chromosomes[i][j]);
            }
            for(int j=point1; j<point2+1; j++)
            {
                child1.chromosomes[i].push_back(parents.second.chromosomes[i][j]);
                child2.chromosomes[i].push_back(parents.first.chromosomes[i][j]);
            }
            for(int j=point2; j<chromosomeSize; j++)
            {
                child1.chromosomes[i].push_back(parents.first.chromosomes[i][j]);
                child2.chromosomes[i].push_back(parents.second.chromosomes[i][j]);
            }
        }
	}
/*
        cout << "\nP1: ";
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < chromosomeSize; j++)
                cout << parents.first.chromosomes[i][j];

            cout << "    ";
        }
        cout << "\nP2: ";
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < chromosomeSize; j++)
                cout << parents.second.chromosomes[i][j];

             cout << "    ";
        }
        cout << "\nC1: ";
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < chromosomeSize; j++)
                cout << child1.chromosomes[i][j];

             cout << "    ";
        }
        cout << "\nC2: ";
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < chromosomeSize; j++)
                cout << child2.chromosomes[i][j];

             cout << "    ";
        }
        cout << endl;
*/
	children = make_pair(child1, child2);

	return children;

}

void mutate2(pair<Individual,Individual>& children)
{
	for(int i=0; i<N; i++)
		for(int j=0; j<chromosomeSize; j++)
		{
			double mut = rnd.random();
			if(mut < mutationRate)
			{
				if(children.first.chromosomes[i][j] == 0)
					children.first.chromosomes[i][j] = 1;
				else children.first.chromosomes[i][j] = 0;
			}
		}

	for(int i=0; i<N; i++)
		for(int j=0; j<chromosomeSize; j++)
		{
			double mut = rnd.random();
			if(mut < mutationRate)
			{
				if(children.second.chromosomes[i][j] == 0)
					children.second.chromosomes[i][j] = 1;
				else children.second.chromosomes[i][j] = 0;
			}
		}

	for(int i=0; i<N; i++)
	{
		children.first.solution.push_back(from_binary(children.first.chromosomes[i], range, gray));
		children.second.solution.push_back(from_binary(children.second.chromosomes[i], range, gray));
	}

	children.first.fitness = fitnessCheck(children.first.solution);
	children.second.fitness = fitnessCheck(children.second.solution);
}

void mutate(pair<Individual,Individual>& children)
{
	int flipBit = rnd.random(chromosomeSize);			// select random "gene" to mutate for child 1

	for(int i=0; i<N; i++)
	{
		if(children.first.chromosomes[i][flipBit] == 0)
			children.first.chromosomes[i][flipBit] = 1;
		else children.first.chromosomes[i][flipBit] = 0;
	}

	flipBit = rnd.random(chromosomeSize);				// select new random "gene" to mutate for child 2

	for(int i=0; i<N; i++)
	{
		if(children.second.chromosomes[i][flipBit] == 0)
			children.second.chromosomes[i][flipBit] = 1;
		else children.second.chromosomes[i][flipBit] = 0;
	}

	for(int i=0; i<N; i++)
	{
		children.first.solution.push_back(from_binary(children.first.chromosomes[i], range, gray));
		children.second.solution.push_back(from_binary(children.second.chromosomes[i], range, gray));
	}

	children.first.fitness = fitnessCheck(children.first.solution);
	children.second.fitness = fitnessCheck(children.second.solution);
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

// Return the fittest individual
Individual fittest()
{
	double best;

	for (int i=0; i <populationSize; i++)
	{
		if(i == 0)
		{
			bestI = i;
			best = population[i].fitness;
		}

		else if(population[i].fitness < best)
		{
			best = population[i].fitness;
			bestI = i;
		}
	}

	return population[bestI];
}

// Wrapper function to choose the right replacement strategy at runtime
void replacement(pair<Individual,Individual> children)
{
	switch(replacementStrategy)
	{
		case 'W':
			replaceWeakest(children);
			break;

		case 'T':
			replaceWithTournament(children);
			break;

		case 'R':
			replaceRandom(children);
			break;

		case 'O':
			replaceOldest(children);
			break;
	}
}

// Replace weakest of 2 individuals
void replaceWithTournament(pair<Individual,Individual> children)
{
	short child;
	Individual toReplace = tournamentSelectionWeaker();

	if(children.first.fitness < children.second.fitness)
	{
		child = 0;
		if(toReplace.fitness > children.first.fitness)
		{
			children.first.number = toReplace.number;
			population[toReplace.number] = children.first;
		}
	}

	else
	{
		child = 1;
		if (toReplace.fitness > children.second.fitness)
		{
			children.second.number = toReplace.number;
			population[toReplace.number] = children.second;
		}
	}

	toReplace = tournamentSelectionWeaker();

	if(child == 0)
	{
		if(children.second.fitness < toReplace.fitness)
		{
			children.second.number = toReplace.number;
			population[toReplace.number] = children.second;
		}
	}
	else															// if we just replaced the weakest with child 2
	{
		if(children.first.fitness < toReplace.fitness)				// if child 1 has lower fitness than the weakest fitness
		{
			children.first.number = toReplace.number;
			population[toReplace.number] = children.first;					// replace the weakest with child 1
		}
	}
}

// Replace weakest
void replaceWeakest(pair<Individual,Individual> children)
{
	short child;
	Individual weakest = weakestLink();

	if(children.first.fitness < children.second.fitness)			// if child 1 has lower fitness
	{
		child = 0;
		if(weakest.fitness > children.first.fitness)				// if child 1 fitness is lower than the weakest fitness
		{
			children.first.number = weakest.number;
			population[weakestI] = children.first;					// replace the weakest with child 1
		}
	}
	else
	{
		child = 1;
		if(weakest.fitness > children.second.fitness)				// if child 2 has lower fitness than the weakest fitness
		{
			children.second.number = weakest.number;
			population[weakestI] = children.second;					// replace the weakest with child 2
		}
	}

	weakest = weakestLink();

	if(child == 0)													// if we just replaced the weakest with child 1
	{
		if(children.second.fitness < weakest.fitness)				// if child 2 has lower fitness than the weakest fitness
		{
			children.second.number = weakest.number;
			population[weakestI] = children.second;					// replace the weakest with child 2
		}
	}
	else															// if we just replaced the weakest with child 2
	{
		if(children.first.fitness < weakest.fitness)				// if child 1 has lower fitness than the weakest fitness
		{
			children.first.number = weakest.number;
			population[weakestI] = children.first;					// replace the weakest with child 1
		}
	}
}

void replaceRandom(pair<Individual,Individual> children)
{
	int rand1 = rnd.random(populationSize);
	int rand2 = rnd.random(populationSize);

	children.first.number = population[rand1].number;
	population[rand1] = children.first;

	children.second.number = population[rand2].number;
	population[rand2] = children.second;
}

void replaceOldest(pair<Individual,Individual> children)
{
	Individual toReplace = oldest();

	children.first.number = toReplace.number;
	population[toReplace.number] = children.first;

	toReplace = oldest();
	children.second.number = toReplace.number;
	population[toReplace.number] = children.second;
}

Individual oldest()
{
	int oldestI = 0;
	int maxAge;
	for(int i=0; i<populationSize; i++)
	{
		if(i == 0)
		{
			maxAge = population[i].age;
			oldestI = i;
		}
		else if(maxAge < population[i].age)
		{
			maxAge = population[i].age;
			oldestI = i;
		}
	}

	return population[oldestI];
}

double evaluate()
{
	double best;

	for(int i=0; i<populationSize; i++)
	{
		if(i == 0)
		{
			best = population[i].fitness;
			bestI = i;
		}
		else if(population[i].fitness < best)
		{
			best = population[i].fitness;
			bestI = i;
		}
	}

	return best;
}

double fitnessCheck(vector<double> solution)
{
	switch(fitnessFunction)
	{
		case 1:
			return f1(solution);
			break;
		case 2:
			return f2(solution);
			break;
		case 3:
			return f3(solution);
			break;
		default:
			return 0.0;
			break;
	}
}

void updateAge(int cycles)
{
	for(int i=0; i<populationSize; i++)
	{
		population[i].age = cycles;
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

// Shekel’s Foxholes
// N = 2
// −65.536 ≤ xi ≤ 65.535
// Optimum is f (x) ≈ 0.99804 at (32, 32)
double f3(const vector<double>& xs)
{
	static double f5_arr[2][25] =
	{
		{
			-32.0, -16.0, 0.0, 16.0, 32.0,
			-32.0, -16.0, 0.0, 16.0, 32.0,
			-32.0, -16.0, 0.0, 16.0, 32.0,
			-32.0, -16.0, 0.0, 16.0, 32.0,
			-32.0, -16.0, 0.0, 16.0, 32.0
		},

		{
			-32.0, -32.0, -32.0, -32.0, -32.0,
			-16.0, -16.0, -16.0, -16.0, -16.0,
			0.0, 	0.0, 	0.0, 	0.0, 	0.0,
		   16.0, 16.0, 16.0, 16.0, 16.0,
		   32.0, 32.0, 32.0, 32.0, 32.0
		}
	};

	double x = xs[0];
	double y = xs[1];
	double sum = 0.0;

	for(int i=0; i<=24; ++i)
	{
		double diff1 = x - f5_arr[0][i];
		double diff2 = y - f5_arr[1][i];
		double subsum=(diff1 * diff1 * diff1 * diff1 * diff1 * diff1) + (diff2 * diff2 * diff2 * diff2 * diff2 * diff2);
		subsum += double(i + 1);
		subsum = double(1) / subsum;
		sum += subsum;
	}

	return 500.0 - (1.0 / (0.002 + sum) );
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
