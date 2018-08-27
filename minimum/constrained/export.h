#pragma once

#ifdef _WIN32
#ifdef minimum_constrained_EXPORTS
#define MINIMUM_CONSTRAINED_API __declspec(dllexport)
#define MINIMUM_CONSTRAINED_API_EXTERN_TEMPLATE
#else
#define MINIMUM_CONSTRAINED_API __declspec(dllimport)
#define MINIMUM_CONSTRAINED_API_EXTERN_TEMPLATE extern
#endif
#else
#define MINIMUM_CONSTRAINED_API
#define MINIMUM_CONSTRAINED_API_EXTERN_TEMPLATE
#endif
