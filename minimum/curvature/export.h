#pragma once

#ifdef _WIN32
#ifdef minimum_curvature_EXPORTS
#define RESEARCH_CURVATURE_API __declspec(dllexport)
#else
#define RESEARCH_CURVATURE_API __declspec(dllimport)
#endif
#else
#define RESEARCH_CURVATURE_API
#endif
