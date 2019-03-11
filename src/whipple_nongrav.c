/** * @file central_force.c
 * @brief   Whipple non-gravitational force
 * @author  Ben Pearson <whereisroadster.com@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2019 Ben Pearson
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Whipple Non-gravity$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 Ben Pearson
 * Implementation Paper    *In progress*
 * Based on                A Review of Comets and Non Gravitational Forces http://adsabs.harvard.edu/full/1994IAUS..160..241Y
 * ======================= ===============================================
 * 
 * Adds a force that acts against gravity, using the formula outlined in A Review of Comets and Non Gravitational Forces, which is as follows:
 * A1*g(r)*r^+A2*g(r)*T^, where g(r)=alpha*(r/r0)^-m*(1+(r/r0)^n)^-k. There is assumed to be only a single source for all of these effects.
 *
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * radiation_source             No          Radiation source, if included. If not included, defaults to 0.
 * Aradial       (double)       No          Corresponds to A1, radial force. Defaults to 0. (AU/(ephemeris day)^2).
 * Atransverse   (double)       No          Corresponds to A2, translocation force. Defaults to 0. (AU/(ephemeris day)^2)
 * r0            (double)       No          Heliocentric distance inside which the bulk of solar insolation goes to sublimation. Defaults to 2.808 for ice. (au)
 * alpha         (double)       No          Normalization factor. Defaults to 0.111262
 * m             (double)       No          Defaults to 2.15
 * n             (double)       No          Defaults to 5.093
 * k             (double)       No          Defaults to 4.6142
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reboundx.h"

static double rebx_get_param_double_default_val(struct reb_particle const *particles, const char* param, double defaultVal) {
	double *ret = rebx_get_param_check(particles, param, REBX_TYPE_DOUBLE);
	if (ret == NULL) {
		return defaultVal;
	}
	return *ret;
}

static void rebx_calculate_whipple_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int source_index, const int N){
    const struct reb_particle source = sim->particles[source_index];
#pragma omp parallel for
    for (int i=0; i<N; i++){
        double ar = rebx_get_param_double_default_val(&particles[i], "Aradial", 0);
		double at = rebx_get_param_double_default_val(&particles[i], "Atransverse", 0);
		if (ar != 0 && at != 0){
			double r0 = rebx_get_param_double_default_val(&particles[i], "r0", 2.808);
			double alpha = rebx_get_param_double_default_val(&particles[i], "alpha", 0.111262);
			double m = rebx_get_param_double_default_val(&particles[i], "m", 2.15);
			double n = rebx_get_param_double_default_val(&particles[i], "n", 5.093);
			double k = rebx_get_param_double_default_val(&particles[i], "k", 4.6142);

			const struct reb_particle p = particles[i];
			const double dx = p.x - source.x;
			const double dy = p.y - source.y;
			const double dz = p.z - source.z;

			const double r = sqrt(dx*dx + dy*dy + dz*dz);
			const double g = alpha * pow(r / r0,-m) * pow(1 + pow(r/ r0,n),-k);

			particles[i].ax += ar * g * dx + at * g * p.vx;
			particles[i].ay += ar * g * dy + at * g * p.vx;
			particles[i].az += ar * g * dz + at * g * p.vz;
		}
    }
}

void rebx_whipple(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N){
	int source_found=0;
    for (int i=0; i<N; i++){
        if (rebx_get_param_check(&particles[i], "radiation_source", REBX_TYPE_INT) != NULL){
            source_found = 1;
            rebx_calculate_whipple_force(sim, particles, i, N);
        }
    }
    if (!source_found){
        rebx_calculate_whipple_force(sim, particles, 0, N);    // default source to index 0 if "radiation_source" not found on any particle
    }
}