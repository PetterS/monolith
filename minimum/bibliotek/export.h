#pragma once

#ifdef _WIN32
#ifdef minimum_bibliotek_EXPORTS
#define MINIMUM_BIBLIOTEK_API __declspec(dllexport)
#else
#define MINIMUM_BIBLIOTEK_API __declspec(dllimport)
#endif
#else
#define MINIMUM_BIBLIOTEK_API
#endif
