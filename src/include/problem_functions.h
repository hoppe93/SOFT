#ifndef _PROBLEM_FUNCTIONS_H
#define _PROBLEM_FUNCTIONS_H

#include "rkf45.h"
#include "IO_data.h"
#include "arguments.h"

typedef struct {
  solution_data* (*output)(solution_data*);
  ode_solution* (*init)(vector*, initial_data*);
  vector* (*equation)(double, vector*); 
}problem;

problem* use_problem(arguments*);

solution_data* output_GCM(solution_data*);
solution_data* output_no_GCM(solution_data*);

ode_solution* init_GCM(vector*, initial_data*);
ode_solution* init_no_GCM(vector*, initial_data*); 

#endif/*_PROBLEM_FUNCTIONS_H*/
