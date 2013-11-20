#include <iostream>
#include <vector>
#include <time.h>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include "mtrandom.h"

using namespace std;

struct Individual
{
	int age;
	int number;
	double fitness;
	bool tabu;
	vector<int> chromosomes[10];
	vector<double> solution;
};

void display();
void usage();
void readFile(string file);
void initialize(vector<Individual>& population);
pair<Individual, Individual> parentSelection();
Individual tabuSelection();
Individual randomSelection();
Individual tournamentSelection();
Individual tournamentSelectionWeaker();
pair<Individual,Individual> crowdSelection();
pair<Individual,Individual> crossover(pair<Individual,Individual> parents);
pair<Individual,Individual> fixedPointCrossover(pair<Individual,Individual> parents);
pair<Individual,Individual> randomPointCrossover(pair<Individual,Individual> parents);
pair<Individual,Individual> twoPointCrossover(pair<Individual,Individual> parents);
void mutate(pair<Individual,Individual>& children);
void fixedMutation(pair<Individual,Individual>& children);
void randomMutation(pair<Individual,Individual>& children);
void replacement(pair<Individual,Individual> children);
void replaceWithTournament(pair<Individual,Individual> children);
void replaceWeakest(pair<Individual,Individual> children);
void replaceRandom(pair<Individual,Individual> children);
void replaceRandom2(pair<Individual,Individual> children);
void replaceOldest(pair<Individual,Individual> children);
void replaceFirstWeaker(pair<Individual,Individual> children);
Individual oldest();
Individual weakestLink();
Individual fittest();
double evaluate();
double fitnessCheck(vector<double> solution);
void updateAge(int cycles);

double f1(const vector<double>& xs);
double f2(const vector<double>& xs);
double f3(const vector<double>& xs);
double rana(const vector<double>& xs);
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
int N;
bool gray;
char replacementStrategy;
int iterations = 0;
double mutationRate;
char selectionStrategy;
char crossoverStrategy;
int maxIterations;
char mutation;
double currentBest;
double lastBest;
int numberOfIterations;

int main(int argc, char *argv[])
{
	cout << "Evolutionary computing!" << endl;
	if(argc != 2)
	{
		usage();
		exit(1);
	}
	else
	{
		rnd.seed_random(time(NULL));
		file = argv[1];
		readFile(file);
	}

		initialize(population);
		
		do
		{
			parents = parentSelection();
			children = crossover(parents);
			mutate(children);
			if(useAge) updateAge(iterations);
			replacement(children);
			cout << std::fixed;
	//		cout << setprecision(10) << "Best solution: " << evaluate() << " x: " << population[bestI].solution[0] << " y:" << population[bestI].solution[1] << " iteration#: " << iterations << "\r";
			currentBest = evaluate();

			if(iterations == 0)
			{
				lastBest = currentBest;
				numberOfIterations = 0;
			}
			
			else if(currentBest < lastBest)
			{
				lastBest = currentBest;
				numberOfIterations = iterations;
			}

			cout << setprecision(9) << "#Iterations: " << iterations << " Best solution found: " << evaluate() << "\r";
			iterations++;
		}while(!done && iterations < maxIterations);

		cout << endl;
		cout << setprecision(9) << " #Iterations: " << numberOfIterations << " Best solution found: " << lastBest << endl;

//	cout << endl;
//	display();

    return 0;
}

void display()
{
	for (int i = 0; i < populationSize; i++)
	{
		cout << "#" << i << " Fitness: " << population[i].fitness << endl;
		cout << "#" << i << " tabu: " << population[i].tabu << endl;
	}
}

void usage()
{
	cout << "Please provide the config file you want to run as an argument" << endl;
	cout << "The config file has 9 parameters which are as follows (in this order):" << endl;
	cout << "Problem number (1-4)" << endl;
	cout << "Population size (whatever you want)" << endl;
	cout << "Selection strategy (R for random selection) (T for tournament selection) (C for crowd selection)" << endl;
	cout << "Use gray coded binary represantation or not (0 for false, 1 for true)" << endl;
	cout << "Replacement strategy to use (W replace weakest) (T replace with tournament) (R replace random) (O replace oldest)" << endl;
	cout << "Mutation type (F for fixed one bit mutation) (R for random mutation using mutation rate)" << endl;
	cout << "Mutation rate (any number from 0 to 1) (only used if R is provided for mutation type)" << endl;
	cout << "Crossover strategy (F for fixed 1 point crossover) (R for random 1 point crossover) (T for random 2 point crossover)" << endl;
	cout << "Maximum number of iterations (whatever you want)" << endl;
}

void readFile(string file)
{
	fstream myfile(file.c_str(), ios_base::in);		// create the filestream

	string selection;
	string replacement;
	string crossover;

    if(myfile)
    {
    	myfile >> fitnessFunction >> populationSize >> selectionStrategy >> gray >> replacementStrategy >> mutation >> mutationRate >> crossoverStrategy >> maxIterations;
		switch(fitnessFunction)						// problem to solve
		{
			case 1:									// problem 1
				range = make_pair(-5.12 , 5.11);	// create the correct range for problem 1
				N = 1;								// set number of chromosomes
				break;
			case 2:									// problem 2
				range = make_pair(-2.048, 2.047);	// create the correct range for problem 2
				N = 2;								// set number of chromosomes
				break;
			case 3:									// problem 3
				range = make_pair(-65.536, 65.535);	// create the correct range for problem 3
				N = 2;								// set number of chromosomes
				break;
			case 4:									// problem 4
				range = make_pair(-512.0, 511.0);	// create the correct range for problem 4
				N = 10;								// set number of chromosomes
				break;
			default:
				break;
		}
		switch(selectionStrategy)
		{
			case 'R':
				selection = "Random selection";
				break;
			case 'T':
				selection = "Tournament selection";
				break;
			case 'C':
				selection = "Crowd selection";
				break;
			default:
				selection = "Tournament selection";
				break;
		}
		switch(replacementStrategy)					// W, T, R, or O
		{
			case 'W':
				replacement = "Replace weakest";
				break;
			case 'O':								// if replacement strategy is replace oldest
				useAge = true;						// set useAge to true
				replacement = "Replace oldest";
				break;
			case 'T':
				replacement = "Replace with tournament selection";
				break;
			case 'R':
				replacement = "Replace random";
				break;
			case 'F':
				replacement = "Replace first weaker";
				break;
			default:
				replacement = "Replace with tournament selection";
				break;
		}
		switch(crossoverStrategy)
		{
			case 'F':
				crossover = "Fixed 1 point crossover";
				break;
			case 'R':
				crossover = "Random 1 point crossover";
				break;
			case 'T':
				crossover = "Random 2 point crossover";
				break;
			default:
				crossover = "Random 2 point crossover";
		}
    }
    else
    {
		cout << "Couldn't read from file" << endl;
		usage();
		exit(1);									// exit the program
	}

	cout << "Problem: " << fitnessFunction << endl;
	cout << "Population size: " << populationSize << endl;
	cout << "Selection strategy: " << selection << endl;
	cout << "Gray coded: " << gray << endl;
	cout << "Replacement strategy: " << replacement <<  endl;
	cout << "Mutation type: " << mutation << endl;
	cout << "Mutation rate: " << mutationRate << endl;
	cout << "Crossover strategy: " << crossover << endl;
	cout << "Max iterations: " << maxIterations << endl;
}

// Initialize the population
void initialize(vector<Individual>& population)
{
	double best;						// best fitness value
	double weakest;						// weakest fitness value

	for(int i=0; i<populationSize; i++)	// initialize population
	{
		Individual individual;			// create a new individual
//		individual.chromosomes[] = new vector<int> chromosomes[N];
		individual.number = i;			// assign the new individual a number
		individual.tabu = false;

		for(int j=0; j<N; j++)			// loop through the chromosomes array
		{
			vector<int> chromosome;															// create a new chromosome vector
			double tmprnd = rnd.random(range.first, range.second);							// get a new random number in the correct range
			individual.chromosomes[j] = chromosome;											// put the chromosome vector into the array
			individual.chromosomes[j] = to_binary(tmprnd, range, chromosomeSize, gray);		// put the binary encoded string into the chromosome vector
			individual.solution.push_back(tmprnd);											// push the random number
		}

		individual.fitness = fitnessCheck(individual.solution);								// get the fitness value for the individual

		if(i == 0)									// if it's the first individual
		{
			weakestI = i;							// weakest index is the first number
			weakest = individual.fitness;			// weakest fitness is the first fitness value
			best = individual.fitness;				// best fitness is the first fitness value
			bestI = i;								// best index is the first number
		}

		else if(individual.fitness > weakest)		// if current individual's fitness is greater than weakest
		{
			weakest = individual.fitness;			// set the weakest to current individual's fitness
			weakestI = i;							// set weakest index to i
		}

		else if(individual.fitness < best)			// if current individual's fitness is less than best
		{
			best = individual.fitness;				// set the best to current individual's fitness
			bestI = i;								// set best index to i
		}

		population.push_back(individual);			// push the new individual into the population
	}
}

// Wrapper function to select the right parent selection strategy
pair<Individual, Individual> parentSelection()
{
	Individual parent1, parent2;
	pair<Individual,Individual> parents;

	switch(selectionStrategy)							// decide on what selection function to use
	{
		case 'R':										// selectionStrategy R means parents are selected randomly

			parent1 = randomSelection();				// select parent 1 randomly
			parent2 = randomSelection();				// select parent 2 randomly

			while(parent1.number == parent2.number)		// while we have the same parents
				parent2 = randomSelection();			// select a new random parent 2

			parent1.tabu = true;
			parent2.tabu = true;

			parents = make_pair(parent1, parent2);		// make the pair
			return parents;								// return parents pair

		case 'T':										// selectionStrategy T means parents are selected with tournament selection

			parent1 = tournamentSelection();			// select parent 1 with tournament selection
			parent2 = tournamentSelection();			// select parent 2 with tournament selection

			while(parent1.number == parent2.number)		// while we have the same parents
				parent2 = tournamentSelection();		// select a new parent 2 with tournament selection

			parents = make_pair(parent1, parent2);		// make the pair
			return parents;								// return parents pair

		case 'C':										// selectionStrategy C means parents are selected with crowd selection

			return crowdSelection();					// return a pair of two parents with crowd selection

		case 'A':

			parent1 = tabuSelection();
			parent2 = tabuSelection();

			parents = make_pair(parent1, parent2);
			return parents;

		default:										// default is tournament selection

			parent1 = tournamentSelection();			// select parent 1 with tournament selection
			parent2 = tournamentSelection();			// select parent 2 with tournament selection

			while(parent1.number == parent2.number)		// while we have the same parents
				parent2 = tournamentSelection();		// select a new parent 2 with tournament selection

			parents = make_pair(parent1, parent2);		// make the pair
			return parents;								// return parents pair
	}
}

Individual tabuSelection()
{
	Individual individual;

	for(int i=0; i<populationSize; i++)
	{
		if(!population[i].tabu)
		{
			individual = population[i];
			individual.tabu = true;
			break;
		}
	}

//	cout << "Individual to return tabu: " << individual.tabu << endl;

	return individual;
}

// Select a random individual from the population
Individual randomSelection()
{
	Individual individual;									// create a new individual
	individual = population[rnd.random(populationSize)];	// get a random individual from the population

	return individual;										// return the individual
}

// Select an individual with tournament selection
Individual tournamentSelection()
{
	int rand1 = rnd.random(populationSize);							// get a new random number
	int rand2 = rnd.random(populationSize);							// get a new random number

	while(rand1 == rand2) rand2 = rnd.random(populationSize);		// while the numbers are the same we get a new random number for rand2

	Individual candidate1 = population[rand1];						// get candidate 1 from population
	Individual candidate2 = population[rand2];						// get candidate 2 from population

	if(candidate1.fitness < candidate2.fitness) return candidate1;	// if candidate 1 has lower fitness value we return candidate 1

	else return candidate2;											// else, candidate 2 has lower (or equal) and we return candidate 2
}

// Select an individual to replace with tournament selection
Individual tournamentSelectionWeaker()
{
	int rand1 = rnd.random(populationSize);							// get a new random number
	int rand2 = rnd.random(populationSize);							// get a new random number

	while(rand1 == rand2) rand2 = rnd.random(populationSize);		// while the numbers are the same we get a new random number for rand2

	Individual candidate1 = population[rand1];						// get candidate 1 from population
	Individual candidate2 = population[rand2];						// get candidate 2 from population

	if(candidate1.fitness > candidate2.fitness) return candidate1;	// if candidate 1 has higher fitness value we return candidate 1

	else return candidate2;											// else, candidate 2 has higher (or equal) and we return candidate 2
}

pair<Individual,Individual> crowdSelection()
{
	vector<Individual> matingPool(10);								// create a vector of size 10 for the mating pool
	int randParent1, closest;										
	Individual parent1, candidate;									// create two individuals to become parents

	pair<Individual,Individual> parents;							// create a pair to return

	for(int i=0; i<10; i++)											// loop through the mating pool
	{
		matingPool[i] = population[rnd.random(populationSize)];		// intialize mating pool with random individuals from the population
	}

	randParent1 = rnd.random(10);									// get a random number
	parent1 = matingPool[randParent1];								// select parent 1 randomly from mating pool

	double distance, bestDistance;

	for(int i=0; i<matingPool.size(); i++)							// loop through the mating pool
	{
		if(parent1.number == matingPool[i].number)					// if we have two instances of the same individual
			continue;												// skip this step in the loop

		candidate = matingPool[i];									// candidate to mate will be the current individual in mating pool

		// calculate distance between parent 1 and the mating candidate
		if(N == 2)
		distance = sqrt( pow(parent1.solution[0] - candidate.solution[0], 2) + pow(parent1.solution[1] - candidate.solution[1], 2) );
/*
		else if(N == 10)
		distance = sqrt( pow(parent1.solution[0] - candidate.solution[0], 2) + pow(parent1.solution[1] - candidate.solution[1], 2) +
						 pow(parent1.solution[2] - candidate.solution[2], 2) + pow(parent1.solution[3] - candidate.solution[3], 2) +
						 pow(parent1.solution[4] - candidate.solution[4], 2) + pow(parent1.solution[5] - candidate.solution[5], 2) +
						 pow(parent1.solution[6] - candidate.solution[6], 2) + pow(parent1.solution[7] - candidate.solution[7], 2) +
						 pow(parent1.solution[8] - candidate.solution[8], 2) + pow(parent1.solution[9] - candidate.solution[9], 2) );
*/
		// if 10 dimensions calculate the minkowski distance between parent 1 and the mating candidate
		else if(N == 10)
		distance = pow( pow(parent1.solution[0] - candidate.solution[0], 10) + pow(parent1.solution[1] - candidate.solution[1], 10) +
						 pow(parent1.solution[2] - candidate.solution[2], 10) + pow(parent1.solution[3] - candidate.solution[3], 10) +
						 pow(parent1.solution[4] - candidate.solution[4], 10) + pow(parent1.solution[5] - candidate.solution[5], 10) +
						 pow(parent1.solution[6] - candidate.solution[6], 10) + pow(parent1.solution[7] - candidate.solution[7], 10) +
						 pow(parent1.solution[8] - candidate.solution[8], 10) + pow(parent1.solution[9] - candidate.solution[9], 10), 1/10);

		if(i == 0)									// if it's the first step in the loop
		{
			bestDistance = distance;				// best distance set to current distance
			closest = i;							// index to closest set to current i
		}
		else if(distance < bestDistance)			// if distance is less than the best distance
		{
			bestDistance = distance;				// set best distance to current distance
			closest = i;							// index to closest set to current i
		}
	}

	parents = make_pair(parent1, matingPool[closest]);	// make the pair with parent 1 and the closest candidate
	
	return parents;										// return the pair
}

// Wrapper function to select correct crossover strategy
pair<Individual,Individual> crossover(pair<Individual,Individual> parents)
{
	switch(crossoverStrategy)
	{
		case 'F':									// fixed point crossover
			return fixedPointCrossover(parents);
		case 'R':									// random point crossover
			return randomPointCrossover(parents);
        case 'T':									// two point crossover
            return twoPointCrossover(parents);
		default:									// if no argument was provided
			return twoPointCrossover(parents);		// default is two point crossover
	}
}

// Create and return a pair of children with fixed point crossover
pair<Individual,Individual> fixedPointCrossover(pair<Individual,Individual> parents)
{
	pair<Individual, Individual> children;		// create an empty pair of individuals
	Individual child1, child2;					// create two new individuals

	for(int i=0; i<N; i++)									// loop through each chromosome
	{
		for(int j=0; j<chromosomeSize/2; j++)				// loop from beginning to middle of each chromosome
		{
			child1.chromosomes[i].push_back(parents.first.chromosomes[i][j]);	// give child 1 first five genes of first parent
			child2.chromosomes[i].push_back(parents.second.chromosomes[i][j]);	// give child 2 first five genes of second parent
		}

	//for(int i=0; i<N; i++)									// loop through each chromosome
		for(int j=chromosomeSize/2; j<chromosomeSize; j++)	// loop from middle to end of each chromosome
		{
			child1.chromosomes[i].push_back(parents.second.chromosomes[i][j]);	// give child 1 last five genes of second parent
			child2.chromosomes[i].push_back(parents.first.chromosomes[i][j]);	// give child 2 last five genes of first parent
		}
	}

	child1.tabu = false;
	child2.tabu = false;

	children = make_pair(child1, child2);		// make the pair

	return children;				// return the pair
}

pair<Individual,Individual> randomPointCrossover(pair<Individual,Individual> parents)
{
	pair<Individual, Individual> children;			// create an empty pair of individuals
	Individual child1, child2;						// create two new individuals

	for(int i=0; i<N; i++)									// loop through each chromosome
	{
		int point = rnd.random(0, chromosomeSize);			// get a random number to split chromosome on

		for(int j=0; j<point; j++)									// loop from beginning to the random point of each chromosome
		{
			child1.chromosomes[i].push_back(parents.first.chromosomes[i][j]);	// give child 1 first x number of genes from first parent
			child2.chromosomes[i].push_back(parents.second.chromosomes[i][j]);	// give child 2 first x number of genes from second parent
		}
		for(int j=point; j<chromosomeSize; j++)						// loop from the random point to end of each chromosome
		{
			child1.chromosomes[i].push_back(parents.first.chromosomes[i][j]);	// give child 1 last x number of genes from first parent
			child2.chromosomes[i].push_back(parents.second.chromosomes[i][j]);	// give child 2 last x number of genes from second parent
		}
	}

	child1.tabu = false;
	child2.tabu = false;

	children = make_pair(child1, child2);		// make the pair

	return children;			// return the pair
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

	child1.tabu = false;
	child2.tabu = false;

	children = make_pair(child1, child2);

	return children;

}

void mutate(pair<Individual,Individual>& children)
{
	switch(mutation)
	{
		case 'R':
			randomMutation(children);
			break;
		case 'F':
			fixedMutation(children);
			break;
		default:
			randomMutation(children);
			break;
	}
}

void randomMutation(pair<Individual,Individual>& children)
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

void fixedMutation(pair<Individual,Individual>& children)
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
//			replaceRandom(children);
			replaceRandom2(children);
			break;

		case 'O':
			replaceOldest(children);
			break;

		case 'F':
			replaceFirstWeaker(children);
			break;

		default:
			replaceWithTournament(children);
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

void replaceRandom2(pair<Individual,Individual> children)
{
	Individual toReplace1, toReplace2;

	toReplace1 = population[rnd.random(populationSize)];

	if(children.first.fitness < toReplace1.fitness)
	{
		children.first.number = toReplace1.number;
		population[toReplace1.number] = children.first;
	}

	toReplace2 = population[rnd.random(populationSize)];

	if(children.second.fitness < toReplace2.fitness)
	{
		children.second.number = toReplace2.number;
		population[toReplace2.number] = children.second;
	}
}

void replaceOldest(pair<Individual,Individual> children)
{
	Individual toReplace = oldest();

	if(children.first.fitness < toReplace.fitness)
	{
		children.first.number = toReplace.number;
		population[toReplace.number] = children.first;	
	}

	toReplace = oldest();
	if(children.second.fitness < toReplace.fitness)
	{
		children.second.number = toReplace.number;
		population[toReplace.number] = children.second;
	}
}

void replaceFirstWeaker(pair<Individual,Individual> children)
{
	for (int i=0; i<populationSize; i++)
	{
		if(children.first.fitness < population[i].fitness)
		{
			children.first.number = population[i].number;
			population[i] = children.first;
		}
	}

	for (int i=0; i<populationSize; i++)
	{
		if(children.second.fitness < population[i].fitness)
		{
			children.second.number = population[i].number;
			population[i] = children.second;
		}
	}
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
		case 4:
			return rana(solution);
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

// Rana’s F2 function
// Minimize f ( ̄) = big complicated thing...
// N = 10
// −512.0 ≤ xi ≤ 511.0
// Optimum value f ( ̄) ≈ −511
double rana(const vector<double>& xs)
{
	const double RANA_WEIGHTS[]= {0.3489, 0.1848, 0.3790, 0.4386,
								  0.9542, 0.1430, 0.7849, 0.3689,
								  0.9767, 0.8163};

	const double RANA_SUM_WEIGHTS=5.3953;
	double sum=0.0;

	for(unsigned int i=0; i<xs.size(); ++i)
	{
		double x1 = xs[i];
		double x2 = xs[(i+1) % xs.size()];
		sum += RANA_WEIGHTS[i] * (x1 * sin(sqrt(fabs(x2 + 1.0 - x1))) *
		cos(sqrt(fabs(x1 + x2 + 1.0))) + (x2 + 1.0) *
		cos(sqrt(fabs(x2 + 1.0 - x1))) *
		sin(sqrt(fabs(x1 + x2 + 1.0))));
	}

	return sum / RANA_SUM_WEIGHTS;
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
