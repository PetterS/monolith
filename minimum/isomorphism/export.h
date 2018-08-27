#pragma once

#ifdef _WIN32
#ifdef minimum_isomorphism_EXPORTS
#define MINIMUM_ISOMORPHISM_API __declspec(dllexport)
#else
#define MINIMUM_ISOMORPHISM_API __declspec(dllimport)
#endif
#else
#define MINIMUM_ISOMORPHISM_API
#endif
