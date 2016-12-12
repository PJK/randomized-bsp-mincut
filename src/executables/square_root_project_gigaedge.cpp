#include "../SquareRootCut.hpp"
#include "input/InputIterator.hpp"
#include "../utils.hpp"

int main(int argc, char* argv[])
{
	if (argc != 4) {
		std::cout << "Usage: square_root prob N SEED" << std::endl;
		return 1;
	}

	float success_probability { std::stof(argv[1], nullptr) };
	unsigned long n = { std::stoul(argv[2]) };
	uint32_t seed = { (uint32_t) std::stoi(argv[3]) };

	MPI_Init(&argc, &argv);

	SquareRootCut cutter(MPI_COMM_WORLD);


	if (cutter.master()) {
		SquareRootCut::Result res = cutter.runGIGAMaster(n, success_probability, seed);
		std::cout << argv[2] << ","
				  << seed << ","
				  << cutter.processors() << ","
				  << n << ","
				  << n * (n - 1) / 2 << ","
				  << res.cuttingTime << ","
				  << res.trials << ","
				  << (res.variant == SquareRootCut::Variant::HIGH_CONCURRENCY ? "high" : "low") << ","
				  << res.weight << std::endl;
	} else {
		cutter.runGIGAWorker(n,success_probability, seed);
	}

	MPI_Finalize();
}
