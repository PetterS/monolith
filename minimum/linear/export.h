#pragma once

#ifdef _WIN32
#ifdef minimum_linear_EXPORTS
#define MINIMUM_LINEAR_API __declspec(dllexport)
#else
#define MINIMUM_LINEAR_API __declspec(dllimport)
#endif
#else
#define MINIMUM_LINEAR_API
#endif
