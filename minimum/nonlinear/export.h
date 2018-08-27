#pragma once

#ifdef _WIN32
#ifdef minimum_nonlinear_EXPORTS
#define MINIMUM_NONLINEAR_API __declspec(dllexport)
#define MINIMUM_NONLINEAR_API_EXTERN_TEMPLATE
#else
#define MINIMUM_NONLINEAR_API __declspec(dllimport)
#define MINIMUM_NONLINEAR_API_EXTERN_TEMPLATE extern
#endif
#else
#define MINIMUM_NONLINEAR_API
#define MINIMUM_NONLINEAR_API_EXTERN_TEMPLATE
#endif
