// Petter Strandmark 2013.
#ifndef MINIMUM_NONLINEAR_FUNCTION_SERIALIZER_H
#define MINIMUM_NONLINEAR_FUNCTION_SERIALIZER_H

#include <iostream>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/nonlinear/function.h>

namespace minimum {
namespace nonlinear {

class Serialize {
	friend MINIMUM_NONLINEAR_API std::ostream& operator<<(std::ostream& out,
	                                                      const Serialize& serializer);
	friend MINIMUM_NONLINEAR_API std::istream& operator>>(std::istream& in,
	                                                      const Serialize& serializer);

   public:
	Serialize(const Function& function)
	    : readonly_function(&function),
	      writable_function(nullptr),
	      user_space(nullptr),
	      factory(nullptr) {}

	Serialize(Function* function,
	          std::vector<double>* input_user_space,
	          const TermFactory& input_factory)
	    : readonly_function(nullptr),
	      writable_function(function),
	      user_space(input_user_space),
	      factory(&input_factory) {}

   private:
	const Function* readonly_function;
	Function* writable_function;
	std::vector<double>* user_space;
	const TermFactory* factory;
};

MINIMUM_NONLINEAR_API std::ostream& operator<<(std::ostream& out, const Serialize& serializer);
MINIMUM_NONLINEAR_API std::istream& operator>>(std::istream& in, Serialize& serializer);
}  // namespace nonlinear
}  // namespace minimum

#endif
