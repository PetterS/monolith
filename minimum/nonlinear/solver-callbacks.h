// Petter Strandmark 2014.
#ifndef MINIMUM_NONLINEAR_SOLVER_CALLBACKS_H
#define MINIMUM_NONLINEAR_SOLVER_CALLBACKS_H

#include <fstream>

#include <minimum/nonlinear/solver.h>

namespace minimum {
namespace nonlinear {

// Saves the current point to a file at
// every iteration.
class FileCallback {
   public:
	FileCallback(std::ofstream& file_) : file(file_) {}

	bool operator()(const CallbackInformation& information) const {
		for (int i = 0; i < information.x->size(); ++i) {
			file << (*information.x)[i] << " ";
		}
		file << std::endl;
		return true;
	}

   private:
	std::ofstream& file;
};
}  // namespace nonlinear
}  // namespace minimum

#endif
