#ifndef PARALLEL_MINIMUM_CUT_SQUAREROOTCUT_HPP
#define PARALLEL_MINIMUM_CUT_SQUAREROOTCUT_HPP

#include "mpi.h"
#include "AdjacencyListGraph.hpp"
#include "SequentialSquareRootCut.hpp"
#include "InputIterator.hpp"
#include "IteratedSparseSampling.hpp"

/**
 * Implements the sqrt(n) `sparse' minimum cut algorithm. This is the top level class
 * that abstracts away the algorithm implementation as well as MPI details.
 */
class SquareRootCut {
	MPI_Comm communicator_;
	int p_, rank_;
	MPI_Datatype mpi_edge_t_;
	double base_case_multiplier_;
	static const int odd_color_ = std::numeric_limits<int>::max();
	void initializeDatatype();
	typedef graph_slice<long> GraphSlice;
	/**
	 * HC minimum group size. A power of 2
	 */

public:
	static const unsigned group_size_ = 2;

	enum Variant { LOW_CONCURRENCY, HIGH_CONCURRENCY };
	struct Result {
		AdjacencyListGraph::Weight weight;
		Variant variant;
		unsigned trials;
		double cuttingTime;
	};

	/**
	 * The communicator ownership is exclusive to SquareRootCut until all members have performed
	 * a runMaster/runWorker.
	 */
	SquareRootCut(MPI_Comm comm, double base_case_multiplier = 2) : communicator_(comm) {
		MPI_Comm_size(communicator_, &p_);
		MPI_Comm_rank(communicator_, &rank_);
		base_case_multiplier_ = base_case_multiplier;
	}

	bool master() const {
		return rank_ == 0;
	}

	unsigned processors() const {
		return (unsigned) p_;
	}

	/**
	 * \return Are the groups singular?
	 */
	bool lowConcurrency(unsigned vertex_count, unsigned edge_count, double success_probability) const;

	Result runMaster(InputIterator & input, double success_probability, uint32_t seed);

	void runWorker(InputIterator & input, double success_probability, uint32_t seed);


	/**
	 * Its too late...
	 * @param n
	 * @param seed
	 * @return
	 */
	Result runGIGAMaster(unsigned long n, double success_probability, uint32_t seed);

	void runGIGAWorker(unsigned long n, double success_probability, uint32_t seed);


	/**
	 * \param n number of vertices
	 * \param m number of edges
	 * \param success_probability Minimum success probability
	 */
	unsigned numberOfTrials(unsigned n, unsigned m, double success_probability) const;

protected:

	/**
	 * \param success_probability Minimum success probability
	 */
	double cPrime(double success_probability) const;

	unsigned intermediate_size(unsigned n, unsigned m) const;

	/**
	 * \param graph
	 * \param success_probability Minimum success probability
	 * \param seed                PRNG seed. Workers will be seeded with {seed + process_index}
	 * \return 2-tuple of (minimum cut weight, number of trials per processor)
	 */
	Result runLowConcurrencyMaster(InputIterator & input, double success_probability, uint32_t seed);

	void runLowConcurrencyWorker();
	/**
	 * Let all nodes arrange into groups. The groups will have a 'local' master that will guide the shrinking.
	 */
	Result participateInGroup(InputIterator & input, double success_probability, uint32_t seed);

	Result participateInGIGAGroup(unsigned long n, double success_probability, uint32_t seed);

	Result runConcurrentMaster(InputIterator & input, double success_probability, uint32_t seed);

	void runConcurrentWorker(InputIterator & input, double success_probability, uint32_t seed);
};


#endif //PARALLEL_MINIMUM_CUT_SQUAREROOTCUT_HPP
