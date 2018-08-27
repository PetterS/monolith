#pragma once

#ifdef _WIN32
#ifdef minimum_core_EXPORTS
#define MINIMUM_CORE_API __declspec(dllexport)
#else
#define MINIMUM_CORE_API __declspec(dllimport)
#endif
#else
#define MINIMUM_CORE_API
#endif
