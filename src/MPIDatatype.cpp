#include "MPIDatatype.hpp"

bool MPIDatatype<AdjacencyListGraph::Edge>::initialized = false;
MPI_Datatype MPIDatatype<AdjacencyListGraph::Edge>::edge_type;
