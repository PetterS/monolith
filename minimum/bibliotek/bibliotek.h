#pragma once

#include <minimum/bibliotek/export.h>

extern "C" {

MINIMUM_BIBLIOTEK_API int add_function(int a, int b);
MINIMUM_BIBLIOTEK_API int sum_function(int* d, int n);

MINIMUM_BIBLIOTEK_API int minimum_linear_test(int* data, int n);

MINIMUM_BIBLIOTEK_API void* minimum_linear_ip(const char* proto_data, int size);
MINIMUM_BIBLIOTEK_API void minimum_linear_ip_destroy(void* ip);
MINIMUM_BIBLIOTEK_API void minimum_linear_ip_save_mps(void* ip, const char* file_name);

MINIMUM_BIBLIOTEK_API void* minimum_linear_solver_default();
MINIMUM_BIBLIOTEK_API void* minimum_linear_solver_glpk();
MINIMUM_BIBLIOTEK_API void* minimum_linear_solver_minisat();
MINIMUM_BIBLIOTEK_API void* minimum_linear_solver_constraint();
MINIMUM_BIBLIOTEK_API void minimum_linear_solver_destroy(void* solver);

MINIMUM_BIBLIOTEK_API void* minimum_linear_solutions(void* solver, void* ip);
MINIMUM_BIBLIOTEK_API int minimum_linear_solutions_get(void* solutions,
                                                       void* ip_,
                                                       double* solution,
                                                       int len);
MINIMUM_BIBLIOTEK_API void minimum_linear_solutions_destroy(void* solutions);

MINIMUM_BIBLIOTEK_API extern const char* error_string;
}