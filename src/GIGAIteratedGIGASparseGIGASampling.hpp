#ifndef PARALLEL_MINIMUM_CUT_XXXXITERATEDSPARSESAMPLING_HPP
#define PARALLEL_MINIMUM_CUT_XXXXITERATEDSPARSESAMPLING_HPP

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

#include "sum_tree.hpp"
#include "prng_engine.hpp"
#include "DisjointSets.hpp"
#include "sorting/SamplingSorter.hpp"
#include "MPIDatatype.hpp"
#include "recursive-contract/graph_slice.hpp"
#include "utils.hpp"

/**
 * Implements ISS, the matrix reduction, and runs the local trials.
 *
 * This class performs the logic of a node within a group that is specified by the communicator.
 */
// TODO move this to a CPP file
class GIGAIteratedGIGASparseGIGASampling {
	MPI_Comm communicator_;
	int rank_, color_, group_size_;
	unsigned target_size_;
	std::vector<AdjacencyListGraph::Edge> edges_slice_;
	const int root_rank_ = 0;
	const float epsilon_ = 0.1f;
	int32_t seed_with_offset_;
	sitmo::prng_engine random_engine_;
	MPI_Datatype mpi_edge_t_;
	unsigned vertex_count_;
	unsigned long n_;

public:

	GIGAIteratedGIGASparseGIGASampling(MPI_Comm communicator, int color, int group_size, unsigned long n, int32_t seed_with_offset, unsigned int target_size) :
			communicator_(communicator),
			color_(color),
			group_size_(group_size),
			target_size_(target_size),
			seed_with_offset_(seed_with_offset),
			random_engine_(seed_with_offset),
			vertex_count_(n),
			n_(n)
	{
		MPI_Comm_rank(communicator_, &rank_);
		mpi_edge_t_ = MPIDatatype<AdjacencyListGraph::Edge>::constructType();
	}

	~GIGAIteratedGIGASparseGIGASampling() {
	}

	int rank() const {
		return rank_;
	}

	bool master() const {
		return rank_ == root_rank_;
	}

	/**
	 * \return Is this the last node in its group?
	 */
	bool last() const {
		return rank_ == group_size_ - 1;
	}

	/**
	 * Read a ~ 1/(group size) slice of edges
	 */
	void loadSlice() {
		unsigned long total_edges = n_ * (n_ - 1) / 2;
		unsigned slice_portion = total_edges / group_size_;
		unsigned slice_from = slice_portion * rank();
		// The last node takes any leftover edges
		unsigned slice_to = last() ? total_edges : slice_portion * (rank() + 1);

		// If it works, it aint stupid, stupid. Enumerate edges.....
		unsigned long edge_ctr = 0;
		for (unsigned long i = 0; i < n_; i++) {
			for (unsigned long j = i + 1; j < n_; j++) {
				if (edge_ctr >= slice_from && edge_ctr < slice_to) {
					edges_slice_.push_back({i, j, 100});
				}
				edge_ctr++;
			}
		}

#ifndef NDEBUG
		std::cout << "Rank " << rank() << " ---> " << edges_slice_ << std::endl;
#endif
		assert(edges_slice_.size() == slice_to - slice_from);
	}

	/**
	 *
	 * \return Has the desired number of vertices been reached?
	 */
	bool samplingTrial() {
		// Compute weight sums
		AdjacencyListGraph::Weight slice_weight = std::accumulate(
				edges_slice_.begin(),
				edges_slice_.end(),
				0,
				[](AdjacencyListGraph::Weight state, AdjacencyListGraph::Edge & edge) { return state + edge.weight; }
		);

		std::vector<AdjacencyListGraph::Weight> weights;
		if (master()) {
			weights.resize(group_size_);
		}

#ifndef NDEBUG
		int size;
		MPI_Type_size(MPI_UNSIGNED_LONG, &size);
		assert(sizeof(AdjacencyListGraph::Weight) == size);
#endif

		MPI_Gather(&slice_weight, 1, MPI_UNSIGNED_LONG, weights.data(), 1, MPI_UNSIGNED_LONG, root_rank_, communicator_);

		// TODO check the sum against the input in debug mode for small inputs

		if (master()) {
			initiateSampling(weights);
		} else {
			acceptSamplingRequest();
			receiveAndApplyMapping();
		}

		return vertex_count_ == target_size_;
	}

	/**
	 * Maps edge endpoints after contraction
	 */
	void receiveAndApplyMapping() {
		std::vector<unsigned> vertex_map(vertex_count_);
		MPI_Bcast(vertex_map.data(), vertex_map.size(), MPI_UNSIGNED, root_rank_, communicator_);

		applyMapping(vertex_map);
		MPI_Bcast(&vertex_count_, 1, MPI_UNSIGNED, root_rank_, communicator_);
	}

	/**
	 * Apply the map to all endpoints, dropping loops
	 * @param vertex_map
	 */
	void applyMapping(const std::vector<unsigned> & vertex_map) {
		assert(vertex_map.size() == vertex_count_);

		std::vector<AdjacencyListGraph::Edge> updated_edges;

		for (auto edge : edges_slice_) {
			edge.from = vertex_map.at(edge.from);
			edge.to = vertex_map.at(edge.to);
			if (edge.from != edge.to) {
				updated_edges.push_back(edge);
			}
		}

		edges_slice_.swap(updated_edges);
	}

	/**
	 * \param edges array of {edge_count >= 0} edges
	 * \param[out] vertices_map preallocated map of {vertex_count >= 0} vertices. Will be filled with partitions label from [0, vertex_count)
	 * \param components_count the desired number of connected components
	 * \param[out] resulting_vertex_count how many vertices remain
	 * \param true if the described prefix exists
	 * I don't trust this code, it has just been copied over -- My past self has written it
	 */
	bool prefixConnectedComponents(const std::vector<AdjacencyListGraph::Edge> & edges,
								   std::vector<unsigned> & vertex_map,
								   unsigned components_count,
								   unsigned & resulting_vertex_count)
	{
		DisjointSets<unsigned> dsets(vertex_map.size());

		if (components_count == 0 || vertex_map.size() == 0) {
			return true;
		}

		unsigned components_active = vertex_map.size();
		bool found = false;

		size_t i = 0;
		for (; i < edges.size() && components_active > components_count; i++) {
			int v1_set = dsets.find(edges.at(i).from),
					v2_set = dsets.find(edges.at(i).to);
			if (v1_set != v2_set) {
				components_active--;
				dsets.unify(v1_set, v2_set);
			}
		}

		if (components_active == components_count) {
			found = true;
		}

		// Also relabel the components to be in [0, new_vertex_count)!
		const long mapping_undefined = -1l;
		std::vector<long> component_labels(vertex_map.size(), mapping_undefined);
		size_t next_label = 0;
		for (unsigned j = 0; j < vertex_map.size(); j++) {
			if (component_labels.at(dsets.find(j)) == mapping_undefined) {
				component_labels.at(dsets.find(j)) = next_label++;
			}

			vertex_map.at(j) = component_labels.at(dsets.find(j));
		}

		resulting_vertex_count = components_active;
		return found;
	}


	bool initiateSampling(std::vector<AdjacencyListGraph::Weight> const & weights) {
		/**
		 * Compute selection per processor
		 */
		unsigned number_of_edges_to_sample = (unsigned) std::pow((float) n_, 1 + epsilon_ / 2);
		// [i] = how many i should sample
		// We use ints to allow using this for the displacements
		std::vector<int> edges_per_processor(group_size_, 0);
		// {number_of_edges_to_sample} processor indices in the order of sampling
		std::vector<size_t> edges_symbolic_order;

		sum_tree<AdjacencyListGraph::Weight> index(weights.data(), group_size_);
		AdjacencyListGraph::Weight sum = index.root();

		std::uniform_int_distribution<AdjacencyListGraph::Weight> uniform_int(1, sum);

		for (size_t i = 0; i < number_of_edges_to_sample; i++) {
			size_t selection = index.lower_bound(uniform_int(random_engine_));

			edges_per_processor[selection]++;
			edges_symbolic_order.push_back(selection);
		}

		/**
		 * Scatter sampling requests
		 */
		int edges_to_sample_locally;
		MPI_Scatter(edges_per_processor.data(), 1, MPI_INT, &edges_to_sample_locally, 1, MPI_INT, 0, communicator_);

		/**
		 * Take part in sampling
		 */
		std::vector<AdjacencyListGraph::Edge> samples = sample(edges_to_sample_locally);

		/**
		 * Gather samples
		 */
		// Allocate space
		std::vector<AdjacencyListGraph::Edge> global_samples(number_of_edges_to_sample);
		// Calculate displacement vector
		std::vector<int> relative_displacements(edges_per_processor);
		relative_displacements.insert(relative_displacements.begin(), 0); // Start at offset zero
		relative_displacements.pop_back(); // Last element not needed

		std::vector<int> displacements;
		std::partial_sum(relative_displacements.begin(), relative_displacements.end(), std::back_inserter(displacements));

		assert(master());
		MPI_Gatherv(
				samples.data(),
				edges_to_sample_locally,
				mpi_edge_t_,
				global_samples.data(),
				edges_per_processor.data(),
				displacements.data(),
				mpi_edge_t_,
				root_rank_,
				communicator_
		);

		assert(global_samples.size() == number_of_edges_to_sample);

		/**
		 * Reorder according to the symbolic order
		 */
		std::vector<AdjacencyListGraph::Edge> ordered_global_samples;
		ordered_global_samples.reserve(global_samples.size());
		std::vector<int> edges_consumed_per_processor(group_size_, 0);

		for (auto processor_selection : edges_symbolic_order) {
			ordered_global_samples.push_back(
					global_samples.at(displacements.at(processor_selection) + edges_consumed_per_processor.at(processor_selection))
			);
			edges_consumed_per_processor.at(processor_selection)++;
		}

		/**
		 * Incremental prefix scan
		 */
		std::vector<unsigned> vertex_map(vertex_count_);
		unsigned resulting_vertex_count;
		prefixConnectedComponents(
				global_samples,
				vertex_map,
				target_size_,
				resulting_vertex_count
		);

		/**
		 * Broadcast mapping
		 */
		MPI_Bcast(vertex_map.data(), vertex_map.size(), MPI_UNSIGNED, root_rank_, communicator_);
		applyMapping(vertex_map);

		/**
		 * Broadcast updated graph size (needed for correct allocations in the next trial)
		 */
		MPI_Bcast(&resulting_vertex_count, 1, MPI_UNSIGNED, root_rank_, communicator_);
		vertex_count_ = resulting_vertex_count;

		// Yay we are done
		return vertex_count_ == target_size_;
	}

	void acceptSamplingRequest() {
		int edges_to_sample_locally;
		MPI_Scatter(nullptr, 1, MPI_INT, &edges_to_sample_locally, 1, MPI_INT, 0, communicator_);

//		std::cout << rank_ << " of " << color_ << " will sample " << edges_to_sample_locally << std::endl;

		std::vector<AdjacencyListGraph::Edge> samples = sample(edges_to_sample_locally);

		MPI_Gatherv(
				samples.data(),
				edges_to_sample_locally,
				mpi_edge_t_,
				nullptr,
				nullptr,
				nullptr,
				mpi_edge_t_,
				root_rank_,
				communicator_
		);
	}

	std::vector<AdjacencyListGraph::Edge> sample(unsigned edge_count) {
		/**
		 * Preprocessing
		 */
		std::vector<AdjacencyListGraph::Weight> edge_weights;
		std::transform(
				edges_slice_.begin(),
				edges_slice_.end(),
				back_inserter(edge_weights),
				[](AdjacencyListGraph::Edge edge) { return edge.weight; }
		);
		std::vector<AdjacencyListGraph::Edge> edges;

		// <3 <3 <3 sum_tree cannot store <= sequences. The client code should handle this, obviously
		if (edges_slice_.size() == 1) {
			edges.resize(edge_count, *edges_slice_.begin());
			return edges;
		}

		assert(edge_weights.size() == edges_slice_.size());
		sum_tree<AdjacencyListGraph::Weight> index(edge_weights.data(), edge_weights.size());

		edges.reserve(edge_count);

		// TODO possibly keep the computed weight
		std::uniform_int_distribution<AdjacencyListGraph::Weight> uniform_int(1, std::accumulate(edge_weights.begin(), edge_weights.end(), 0));

		/**
		 * Sampling
		 */
		for (size_t i = 0; i < edge_count; i++) {
			edges.push_back(
					edges_slice_.at(index.lower_bound(uniform_int(random_engine_)))
			);
		}

		return edges; // NRVO
	}

	/**
	 * Shrinks the graph
	 */
	void shrink() {
		while (!samplingTrial()) {}
	}

	/**
	 * Reduce results across all nodes
	 * @return input for recursive contract
	 */
	graph_slice<long> reduce() {
		// To enable correct sorting
		std::for_each(
				edges_slice_.begin(),
				edges_slice_.end(),
				[](AdjacencyListGraph::Edge & e) { e.normalize(); }
		);

		SamplingSorter<AdjacencyListGraph::Edge> sorter(communicator_, std::move(edges_slice_), seed_with_offset_);
		std::vector<AdjacencyListGraph::Edge> sorted_slice = sorter.sort();

		// Fun part: reduce edges
		// Note that MPI_Scan is not directly applicable since it produces the same number of outputs, and there
		// is no simple way to have 'carry' and still mark all edges but one for each pair of vertices void
		std::vector<AdjacencyListGraph::Edge> locally_reduced_slice;
		if (! sorted_slice.empty()) {
			auto sorted_edges_iterator = sorted_slice.begin();
			locally_reduced_slice.push_back(*(sorted_edges_iterator++));
			while (sorted_edges_iterator != sorted_slice.end()) {
				if (locally_reduced_slice.back() == *sorted_edges_iterator) {
					locally_reduced_slice.back().weight += sorted_edges_iterator->weight;
				} else {
					locally_reduced_slice.push_back(*sorted_edges_iterator);
				}
				++sorted_edges_iterator;
			}
		}

		/*
		 * Now it gets even funnier. We need to merge the edges on the boundaries of slices.
		 * Slice can be either:
		 *  - empty
		 *  - singular
		 *  - have two or more distinct edges (regular)
		 *
		 * There can be a segment of singular processor of length up to p:
		 * ... || (a, b) || (a, b) || .... || (a, b) || (c, d) || ...
		 * Moreover, there may be any number of empty processors at the end.
		 *
		 * How we deal with this:
		 * Every processor 'offers' it's first edge to the preceding one. These edges are all-gathered.
		 * Processors that don't have an edge to offer will supply a dummy value.
		 *
		 * If a processor was offered a dummy edge, it ignores it.
		 *
		 * All processors that offered an edge except the first one remote that edge from their slice.
		 *
		 * Now we consider processors that were offered a real edge. Since it may be that there are
		 * several following processors that offer the same value, the leftmost processor will collect
		 * all the edges that should be combined with the first edge it receives
		 *
		 * Note that there might be gaps with processors owning no edges after this process.
		 */

		/** Boundary edges offred by each processor */
		std::vector<AdjacencyListGraph::Edge> boundary_edges(group_size_);

		/** Dummy value to signal empty slice */
		AdjacencyListGraph::Edge dummy_edge = {
				std::numeric_limits<unsigned >::max(),
				std::numeric_limits<unsigned >::max(),
				std::numeric_limits<AdjacencyListGraph::Weight>::max()
		};

		if (locally_reduced_slice.empty()) {
			locally_reduced_slice.push_back(dummy_edge);
		}

		MPI_Allgather(
				&locally_reduced_slice.front(),
				1,
				mpi_edge_t_,
				boundary_edges.data(),
				1,
				mpi_edge_t_,
				communicator_
		);

		// Give up the edge we've offered
		if (rank_ != 0) {
			locally_reduced_slice.erase(locally_reduced_slice.begin());
		}

		// If we are not the last processor and we've been offered a real edge
		if ((rank_ < (group_size_ - 1)) && (boundary_edges.at(rank_ + 1) != dummy_edge)) {
			// ... and we are the leftmost processor for this streak of edges
			if (boundary_edges.at(rank_) != boundary_edges.at(rank_ + 1)) {
				// ... then collect all the edges ...
				size_t next_edge_to_collect = rank_ + 1;

				// ... if the offered edge is different, add it explicitly ...
				if (locally_reduced_slice.back() != boundary_edges.at(next_edge_to_collect)) {
					locally_reduced_slice.push_back(boundary_edges.at(next_edge_to_collect++)); // and move on to the next edges
				}

				// ... and take all matching edges
				while ((next_edge_to_collect < size_t(group_size_)) && (boundary_edges.at(next_edge_to_collect) == locally_reduced_slice.back())) {
					locally_reduced_slice.back().weight += boundary_edges.at(next_edge_to_collect++).weight;
				}
			}
		}

		// TODO scope the locals so that we can throw away the sorted list


		/*
		 * Matrix construction time.
		 *
		 * Lukas hath said "for a matrix n' rows and a subset of currently p' processors:
		 * the matrix is distributed row-wise such that each processor holds ceil(n/p) consecutive rows of the matrix."
		 * "note also that the matrix needs to be padded with zeroes, such that p' divides the number of rows."
		 *
		 * A minor problem is that it is easy to determine where to send the edges based on the first component (they
		 * are sorted by it and thus we can compute non-uniform all-to-all displacements, but this will only initialize
		 * one half of the entries, the corresponding symmetric entries will be missing. We just add a transpose manually.
		 *
		 * TODO is there a nicer way to do it already built in?
		 */

		{
			int rows_per_processor = (int) std::ceil(float(vertex_count_) / group_size_);
			int row_col_size = rows_per_processor * group_size_;

			/** p_i will receive edges_per_processor[i] edges to place in its slice of rows */
			std::vector<int> edges_per_processor(group_size_, 0);

			// Duplicate edges to get both triangles
			{
				size_t original_size = locally_reduced_slice.size();
				for (size_t i { 0 }; i < original_size; i++) {
					AdjacencyListGraph::Edge e = locally_reduced_slice.at(i);
					std::swap(e.from, e.to);
					locally_reduced_slice.push_back(e);
				}
			}

			std::sort(locally_reduced_slice.begin(), locally_reduced_slice.end());

			for (auto const edge : locally_reduced_slice) {
				edges_per_processor.at(edge.from / rows_per_processor)++;
			}

			/** Receive displacements */
			std::vector<int> edges_to_receive(group_size_);

			/* Exchange segment sizes */
			MPI_Alltoall(
					edges_per_processor.data(),
					1,
					MPI_INT,
					edges_to_receive.data(),
					1,
					MPI_INT,
					communicator_
			);

			/** Compute offsets */
			std::vector<int> arriving_edges_offsets = MPIUtils::prefix_offsets(edges_to_receive);
			std::vector<int> edges_groups_offsets = MPIUtils::prefix_offsets(edges_per_processor);
			int number_of_edges_to_receive = std::accumulate(edges_to_receive.begin(), edges_to_receive.end(), 0);
			/** Target buffer */
			std::vector<AdjacencyListGraph::Edge> incoming_edges_buffer(number_of_edges_to_receive);

			MPI_Alltoallv(
					locally_reduced_slice.data(),
					edges_per_processor.data(),
					edges_groups_offsets.data(),
					mpi_edge_t_,
					incoming_edges_buffer.data(),
					edges_to_receive.data(),
					arriving_edges_offsets.data(),
					mpi_edge_t_,
					communicator_
			);

			// C++ it like it's 1998
			long * rows_slice = new long[row_col_size * rows_per_processor];
			std::fill(rows_slice, rows_slice + row_col_size * rows_per_processor, 0l);

			/** Global index of our first row */
			size_t row_offset = rows_per_processor * rank_;

			for (auto edge : incoming_edges_buffer) {
				rows_slice[(edge.from - row_offset) * row_col_size + edge.to] = edge.weight;
			}

			DebugUtils::print(rank_, [&](std::ostream & out) {
				out << "About to construct slice with: vertex count: " << vertex_count_
					<< " rows_per_processor: " << rows_per_processor
					<< " row_col_size: " << row_col_size;
			});

			graph_slice<long> graph_slice(
					vertex_count_,
					rows_per_processor,
					rank_,
					row_col_size,
					rows_slice
			);

			assert (graph_slice.invariant());

			return graph_slice;
		}
	}

	unsigned vertexCount() {
		return vertex_count_;
	}
};

#endif //PARALLEL_MINIMUM_CUT_ITERATEDSPARSESAMPLING_HPP
