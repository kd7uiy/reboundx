#ifndef _LIBrebxf_H
#define _LIBrebxf_H
#ifndef M_PI 
#define M_PI 3.1415926535879323846
#endif
#include "rebound.h"

struct rebxf_params {
	int allocatedN;

	double *tau_a;
	double *tau_e;
	double *tau_inc;
	double *tau_pomega;

	double e_damping_p; // p parameter from Deck & Batygin (2015) for how e-damping
	// is coupled to a-damping at order e^2
	// p = 1 : e-damping at const angular momentum.  p = 0 : no contribution to a-damping
	// equal to p/3 with p defined as in Goldreich & Schlichting 2014
};

double* rebxf_get_tau_a(struct reb_simulation* sim);
void rebxf_set_tau_a(struct reb_simulation* sim, double* tau_a);

double* rebxf_get_tau_e(struct reb_simulation* sim);
void rebxf_set_tau_e(struct reb_simulation* sim, double* tau_e);

double* rebxf_get_tau_inc(struct reb_simulation* sim);
void rebxf_set_tau_inc(struct reb_simulation* sim, double* tau_inc);

double* rebxf_get_tau_pomega(struct reb_simulation* sim);
void rebxf_set_tau_pomega(struct reb_simulation* sim, double* tau_pomega);

struct rebxf_params* rebxf_addxf(struct reb_simulation* sim);

void rebxf_forces(struct reb_simulation* const sim);

void rebxf_modify_elements(struct reb_simulation* const sim);

#endif