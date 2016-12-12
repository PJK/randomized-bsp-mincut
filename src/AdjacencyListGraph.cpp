#include "AdjacencyListGraph.hpp"
#include <unordered_set>

void AdjacencyListGraph::addEdge(unsigned from, unsigned to, Weight weight)
{
	if (&edges_ != parent_edges_) {
		edges_.assign(parent_edges_->begin(), parent_edges_->end());
		parent_edges_ = &edges_;
	}
	assert(from < vertex_count_);
	assert(to < vertex_count_);

	std::tie(from, to) = normalize(from, to);
	edges_.push_back({ from, to, weight });
}

void AdjacencyListGraph::addEdge(Edge e)
{
	addEdge(e.from, e.to, e.weight);
}

unsigned AdjacencyListGraph::maxVertexID() const
{
	unsigned max = { 0 };
	for (auto const & edge : *parent_edges_)
		max = std::max(max, std::max(edge.from, edge.to));

	return max;
}

/**
 * @return The graph with singleton vertices removed and edges renamed so that it is a connected graph on [0, vertex_count)
 */
std::unique_ptr<AdjacencyListGraph> AdjacencyListGraph::compact() const {
	std::unique_ptr<AdjacencyListGraph> result(new AdjacencyListGraph(vertex_count()));
	// Mapping old_id -> new_id
	std::vector<unsigned> mapping(maxVertexID() + 1, std::numeric_limits<unsigned>::max());

	// Construct the mapping
	unsigned next_id { 0 };
	for (auto const & edge : *parent_edges_) {
		if (mapping.at(edge.from) == std::numeric_limits<unsigned>::max())
			mapping.at(edge.from) = next_id++;
		if (mapping.at(edge.to) == std::numeric_limits<unsigned>::max())
			mapping.at(edge.to) = next_id++;
	}

	for (auto const & edge : *parent_edges_)
		result->addEdge(mapping.at(edge.from), mapping.at(edge.to), edge.weight);

	return result;
}

void AdjacencyListGraph::weaklyContractEdge(unsigned from, unsigned to)
{
	// std::cout << "Contracting " << from << " -- " << to << std::endl;
	// Since we are in a weak contraction phase, we need to check for the actual state of partitions
	unsigned actual_from = { disjoint_sets_.find(from) },
			actual_to = { disjoint_sets_.find(to) };

	// The edge does not exist. This is OK since this might be due to a previous contraction from the same sample
	if (actual_from == actual_to)
		return;

	disjoint_sets_.unify(actual_from, actual_to);
	assert(disjoint_sets_.find(actual_from) == disjoint_sets_.find(actual_to));

	vertex_count_--;
}

/**
 * Contraction with loop removal and edge merging
 *
 * Beware: This does not check for existence of the edge! It will happily merge any vertices you pass to it!
 */
void AdjacencyListGraph::contractEdge(unsigned from, unsigned to)
{
	weaklyContractEdge(from, to);
	//finalizeContractionPhase();
	// TODO assert canonical form
}

//
//namespace std
//{
//	template<>
//	struct hash<AdjacencyListGraph::Edge>
//	{
//		typedef AdjacencyListGraph::Edge argument_type;
//		typedef std::size_t result_type;
//
//		result_type operator()(argument_type const& e) const
//		{
//			// https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
//			//return (e.from + e.to) * (e.from + e.to + 1) / 2 + e.to;
//			return 32000 * e.from + e.to;
//		}
//	};
//}

struct EdgeHasher {
	typedef AdjacencyListGraph::Edge argument_type;
	typedef std::size_t result_type;
	unsigned vertex_count_;
	// Highly potent magic - Hasher is supposed to be stateless. Luckily, this is C++
	mutable std::vector<unsigned> mapping;
	sitmo::prng_engine * random_;

	EdgeHasher(unsigned vertex_count, unsigned max_vertex_id, sitmo::prng_engine * random) : vertex_count_(vertex_count),
																							 mapping(max_vertex_id + 1, std::numeric_limits<unsigned>::max()),
																							 random_(random)
	{}

	result_type operator()(argument_type const& e) const
	{
		// https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
		//return (e.from + e.to) * (e.from + e.to + 1) / 2 + e.to;

		// Retarted icc doesn't properly handle initializer list
		AdjacencyListGraph::Edge eval = e;
		std::tie(eval.from, eval.to) = std::make_tuple(relabel(e.from), relabel(e.to));
		return eval.from xor eval.to;
	}

	unsigned relabel(unsigned vtex) const
	{
		auto rep = mapping.at(vtex);
		if (rep == std::numeric_limits<unsigned>::max()) {
			unsigned label = random_->operator()();
			mapping.at(vtex) = label;
			return label;
		} else return rep;
	}
};

/**
 * Remove loop and merge edges after a series of weaklyContractEdge() calls
 */
void AdjacencyListGraph::finalizeContractionPhase(sitmo::prng_engine * random)
{
	std::unordered_set<Edge, EdgeHasher> uniq_edges(1, EdgeHasher(vertex_count_, maxVertexID(), random));
	// Asymptotically pessimistic, constantly optimistic :)
	uniq_edges.reserve(unsigned(pow(vertex_count_, 2.0/3) / 4));

	for (Edge edge : *parent_edges_) {
		std::tie(edge.to, edge.from) = normalize(disjoint_sets_.find(edge.to),  disjoint_sets_.find(edge.from));

		if (edge.from == edge.to)
			continue;

		std::unordered_set<Edge, EdgeHasher>::iterator reference = uniq_edges.find(edge);
		if (reference == uniq_edges.end())
			uniq_edges.insert(edge);
		else {
			// Hack around the constness of the reference (The iterator is constant so that the elements cannot change,
			// but this fugly hack is OK since our hash function is not dependent on the value of weight. <3
			const_cast<Weight &>(reference->weight) = reference->weight + edge.weight;
		}
	}

	edges_.resize(uniq_edges.size());
	edges_.assign(uniq_edges.begin(), uniq_edges.end());
	parent_edges_ = &edges_;
}

AdjacencyListGraph::EdgeList const & AdjacencyListGraph::edges() const
{
	return *parent_edges_;
}


std::ostream & operator<< (std::ostream & out, AdjacencyListGraph const & graph)
{
	for (auto & edge : *graph.parent_edges_)
		out << edge << std::endl;
	return out;
}

std::tuple<unsigned, unsigned> AdjacencyListGraph::normalize(unsigned v1, unsigned v2) const
{
	return std::tuple<unsigned, unsigned>(std::min(v1, v2), std::max(v1, v2));
}
