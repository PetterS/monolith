#pragma once

#include <utility>

namespace minimum {
namespace core {

//
// at_scope_exit( statement; ) executes statement at the end
// of the current scope.
//
template <typename F>
class ScopeGuard {
   public:
	ScopeGuard(F&& f) : f(std::forward<F>(f)) {}

	ScopeGuard(ScopeGuard&& guard) : f(std::move(guard.f)), active(guard.active) {
		guard.dismiss();
	}

	~ScopeGuard() {
		if (active) {
			f();
		}
	}

	ScopeGuard(const ScopeGuard&) = delete;
	ScopeGuard& operator=(const ScopeGuard&) = delete;

	void dismiss() { active = false; }

   private:
	F f;
	bool active = true;
};

template <typename F>
[[nodiscard]] ScopeGuard<F> make_scope_guard(F&& f) {
	return std::move(ScopeGuard<F>(std::forward<F>(f)));
};
}  // namespace core
}  // namespace minimum

#define MINIMUM_JOIN_PP_SYMBOLS_HELPER(arg1, arg2) arg1##arg2
#define MINIMUM_JOIN_PP_SYMBOLS(arg1, arg2) MINIMUM_JOIN_PP_SYMBOLS_HELPER(arg1, arg2)
#define at_scope_exit(code)                                          \
	auto MINIMUM_JOIN_PP_SYMBOLS(spii_scope_exit_guard_, __LINE__) = \
	    ::minimum::core::make_scope_guard([&]() { code; })
