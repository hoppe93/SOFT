/* Manages sycamera output processors */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "sycout.h"
#include "util.h"

const int NUMBER_OF_SYCOUTS=6;
sycout_type *all_sycouts=NULL;
int *sycout_selected=NULL, sycout_nselected=0;

void sycout_init_handler(void) {
    all_sycouts = malloc(sizeof(sycout_type) * NUMBER_OF_SYCOUTS);

    all_sycouts[0].name = setname("image");
    all_sycouts[0].deinit_run = sycout_image_deinit_run;
    all_sycouts[0].init = sycout_image_init;
    all_sycouts[0].init_run = sycout_image_init_run;
    all_sycouts[0].init_particle = sycout_image_init_particle;
    all_sycouts[0].step = sycout_image_step;
    all_sycouts[0].write = sycout_image_write;

    all_sycouts[1].name = setname("topview");
    all_sycouts[1].deinit_run = sycout_topview_deinit_run;
    all_sycouts[1].init = sycout_topview_init;
    all_sycouts[1].init_run = sycout_topview_init_run;
    all_sycouts[1].init_particle = sycout_topview_init_particle;
    all_sycouts[1].step = sycout_topview_step;
    all_sycouts[1].write = sycout_topview_write;

	all_sycouts[2].name = setname("space3d");
	all_sycouts[2].deinit_run = sycout_space3d_deinit_run;
    all_sycouts[2].init = sycout_space3d_init;
    all_sycouts[2].init_run = sycout_space3d_init_run;
    all_sycouts[2].init_particle = sycout_space3d_init_particle;
    all_sycouts[2].step = sycout_space3d_step;
    all_sycouts[2].write = sycout_space3d_write;

	all_sycouts[3].name = setname("spectrometer");
	all_sycouts[3].deinit_run = sycout_spectrometer_deinit_run;
    all_sycouts[3].init = sycout_spectrometer_init;
    all_sycouts[3].init_run = sycout_spectrometer_init_run;
    all_sycouts[3].init_particle = sycout_spectrometer_init_particle;
    all_sycouts[3].step = sycout_spectrometer_step;
    all_sycouts[3].write = sycout_spectrometer_write;

	all_sycouts[4].name = setname("green");
	all_sycouts[4].deinit_run = sycout_green_deinit_run;
    all_sycouts[4].init = sycout_green_init;
    all_sycouts[4].init_run = sycout_green_init_run;
    all_sycouts[4].init_particle = sycout_green_init_particle;
    all_sycouts[4].step = sycout_green_step;
    all_sycouts[4].write = sycout_green_write;

    all_sycouts[5].name = setname("polimage");
    all_sycouts[5].deinit_run = sycout_polimage_deinit_run;
    all_sycouts[5].init = sycout_polimage_init;
    all_sycouts[5].init_run = sycout_polimage_init_run;
    all_sycouts[5].init_particle = sycout_polimage_init_particle;
    all_sycouts[5].step = sycout_polimage_step;
    all_sycouts[5].write = sycout_polimage_write;
}
void sycout_prepare_run(void) {
    int i;
    for (i = 0; i < sycout_nselected; i++) {
        all_sycouts[sycout_selected[i]].init_run();
    }
}
void sycout_init_particle(particle *p) {
    int i;
    for (i = 0; i < sycout_nselected; i++) {
        all_sycouts[sycout_selected[i]].init_particle(p);
    }
}
void sycout_deinit_run(void) {
    int i;
    for (i = 0; i < sycout_nselected; i++) {
        all_sycouts[sycout_selected[i]].deinit_run();
    }
}
void sycout_step(struct sycout_data *data) {
    int i;
    for (i = 0; i < sycout_nselected; i++) {
        all_sycouts[sycout_selected[i]].step(data);
    }
}
void sycout_select(const char *name, struct general_settings *set, int nsets) {
    int i, j;
    for (i = 0; i < NUMBER_OF_SYCOUTS; i++) {
        if (!strcmp(all_sycouts[i].name, name)) {
            sycout_selected = realloc(sycout_selected, sizeof(int)*(sycout_nselected+1));
            if (sycout_selected == NULL) {
                fprintf(stderr, "Memory error! Unable to select sycout '%s'\n", name);
                exit(-1);
            }

            sycout_selected[sycout_nselected] = i;
            sycout_nselected++;

            for (j = 0; j < nsets; j++) {
                if (!strcmp((set+j)->name, name)) {
                    (all_sycouts+i)->init(set+j);
                    return;
                }
            }

            fprintf(stderr, "WARNING: No settings provided for sycout map!\n");
            return;
        }
    }

    fprintf(stderr, "Unrecognized sycout: %s\n", name);
    exit(-1);
}

void sycout_output_all(int mpi_rank, int nprocesses) {
    int i;
    for (i = 0; i < sycout_nselected; i++) {
        all_sycouts[sycout_selected[i]].write(mpi_rank, nprocesses);
    }
}
