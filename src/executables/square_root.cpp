#include "../SquareRootCut.hpp"
#include "input/InputIterator.hpp"
#include "../utils.hpp"

int main(int argc, char* argv[])
{
	if (argc != 4) {
		std::cout << "Usage: square_root PROBABILITY INPUT_FILE SEED" << std::endl;
		return 1;
	}

	float success_probability { std::stof(argv[1], nullptr) };
	uint32_t seed = { (uint32_t) std::stoi(argv[3]) };

	MPI_Init(&argc, &argv);

	SquareRootCut cutter(MPI_COMM_WORLD);

	InputIterator input(argv[2]);

	if (cutter.master()) {
		SquareRootCut::Result res = cutter.runMaster(input, success_probability, seed);
		std::cout << argv[2] << ","
				  << seed << ","
				  << cutter.processors() << ","
				  << input.vertexCount() << ","
				  << input.edgeCount() << ","
				  << res.cuttingTime << ","
				  << res.trials << ","
				  << (res.variant == SquareRootCut::Variant::HIGH_CONCURRENCY ? "high" : "low") << ","
				  << res.weight << std::endl;
	} else {
		cutter.runWorker(input, success_probability, seed);
	}

	MPI_Finalize();
}
