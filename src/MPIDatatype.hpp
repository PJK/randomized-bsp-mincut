#ifndef PARALLEL_MINIMUM_CUT_MPIDATATYPE_HPP
#define PARALLEL_MINIMUM_CUT_MPIDATATYPE_HPP

#include "AdjacencyListGraph.hpp"
#include "mpi.h"

/**
 * Trait container to allow external implementation for any type
 */
template <typename ElementType>
struct MPIDatatype {};

template <>
struct MPIDatatype<int> {
	static MPI_Datatype constructType() {
		return MPI_INT;
	}
};

template <>
struct MPIDatatype<AdjacencyListGraph::Edge> {
	// TODO global registry -- free datatype
	static MPI_Datatype edge_type;
	static bool initialized;

	static MPI_Datatype constructType() {
		if (!initialized) {
			int blocklengths[3] = {1, 1, 1};

			// This leaks abstraction :(
			MPI_Datatype types[3] = { MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED_LONG };

			MPI_Aint offsets[3];

			offsets[0] = offsetof(AdjacencyListGraph::Edge, from);
			offsets[1] = offsetof(AdjacencyListGraph::Edge, to);
			offsets[2] = offsetof(AdjacencyListGraph::Edge, weight);

			MPI_Type_create_struct(3, blocklengths, offsets, types, &edge_type);
			MPI_Type_commit(&edge_type);

			initialized = true;
		}
		return edge_type;
	};
};

#endif //PARALLEL_MINIMUM_CUT_MPIDATATYPE_HPP
