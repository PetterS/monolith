#pragma once

#ifdef _WIN32
#ifdef minimum_algorithms_EXPORTS
#define MINIMUM_ALGORITHMS_API __declspec(dllexport)
#else
#define MINIMUM_ALGORITHMS_API __declspec(dllimport)
#endif
#else
#define MINIMUM_ALGORITHMS_API
#endif
