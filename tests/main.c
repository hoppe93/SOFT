/* Main test program */

#include <stdio.h>
#include <string.h>
#include <magnetic_field.h>
#include <settings.h>
#include <sycamera.h>
#include <tools.h>
#include <util.h>

#include "test.h"

/* Initialization helpers (optional to call) */
void init_magnetic_circular(void) {
	magnetic_init();
	magnetic_handler *mh = magnetic_handler_select("circular");
	struct general_settings *mset = malloc(sizeof(struct general_settings));
	mset->name = "circular";
	mset->n = 4;
	mset->setting = malloc(sizeof(char*)*mset->n);
	mset->value = malloc(sizeof(char*)*mset->n);
	mset->setting[0] = setname("B0");
	mset->value[0] = setname("6");
	mset->setting[1] = setname("major_radius");
	mset->value[1] = setname("0.68");
	mset->setting[2] = setname("safety_factor");
	mset->value[2] = setname("1");
	mset->setting[3] = setname("minor_radius");
	mset->value[3] = setname("0.22");

	mh->init(mset);
	mh->init_run();
	mh->init_particle();
}
void init_sycamera(void) {
	const unsigned int NSETTINGS = 6;

	struct general_settings *ts = malloc(sizeof(struct general_settings));
	ts->name = setname("sycamera");
	ts->setting = malloc(sizeof(char*)*NSETTINGS);
	ts->value = malloc(sizeof(char*)*NSETTINGS);

	ts->setting[0] = setname("aperture");
	ts->value[0] = setname("0.006");
	ts->setting[1] = setname("direction");
	ts->value[1] = setname("-1,-4,0");
	ts->setting[2] = setname("radiation");
	ts->value[2] = setname("synchrotron");
	ts->setting[3] = setname("cone");
	ts->value[3] = setname("delta");
	ts->setting[4] = setname("vision_angle");
	ts->value[4] = setname("0.8");
	ts->setting[5] = setname("position");
	ts->value[5] = setname("0,-1.069,-0.2265");
	ts->setting[6] = setname("product");
	ts->value[6] = setname("image");
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		fprintf(stderr, "ERROR: No test specified.\n");
		return 1;
	} else if (argc > 2) {
		fprintf(stderr, "ERROR: Only one test can be run at a time.\n");
		return 1;
	}

	if (!strcmp(argv[1], "bsspec")) {
		return test_bss();
	} else if (!strcmp(argv[1], "bsdist")) {
		return test_bsdist();
	} else if (!strcmp(argv[1], "distfunc")) {
		return test_distfunc();
	} else if (!strcmp(argv[1], "hyperbola")) {
		return test_hyp();
	} else {
		fprintf(stderr, "ERROR: Unrecognized test: '%s'.\n", argv[1]);
	}

	return 0;
}
