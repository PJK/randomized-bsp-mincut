#include "utils.hpp"
#include <numeric>
#include <algorithm>

std::vector<int> MPIUtils::prefix_offsets(const std::vector<int> & sizes) {
	std::vector<int> offsets;
	offsets.push_back(0);
	std::partial_sum(
			sizes.begin(),
			--sizes.end(),
			std::back_inserter(offsets)
	);
	return offsets; // NRVO
}