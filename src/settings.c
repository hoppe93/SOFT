/* Load settings */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "domain.h"
#include "global.h"
#include "particles.h"
#include "settings.h"
#include "util.h"

char DEBUG_OUTPUT=0;
FILE *settings_file;
char settings_eof=0;
settings_token *tkn;
int settings_line=1, settings_charpos=0;
char settings_token_names[10][10]={
	"invalid",
	"name",
	"value",
	"magnetic",
	"tool",
	"particle",
	"particles",
	"sycout",
	"{",
	"}"
};

particle *default_particle;

/* Get the next (meaningful) character from input file */
char settings_getc(void) {
	char c, q;
	if (settings_eof) return 0;

	do {
		c = (char)fgetc(settings_file);
		settings_charpos++;
		if (c == '\n') {
			settings_line++;
			settings_charpos=0;
		}

		if (c==EOF) {settings_eof=1; return 0;}
		else if (c=='#') {
			do {q = settings_getc();} while (q!=0 && q != '\n');
			return q;
		//} else if (c!=' '&&c!='\t') return c;
		} else return c;
	} while (c != EOF);

	settings_eof = 1;
	return 0;
}
void settings_rewind(void) {
	fseek(settings_file, -1, SEEK_CUR);
	settings_charpos--;
}
enum settings_token_type token(void) {
	char buffer[1024], c;
	int length=0, val=0;

	do {
		c = settings_getc();
		switch (c) {
			case '\n':
			case ';':
				if (length > 0)
					goto END_OF_LOOP;
				else break;
			case '{':
			case '}':
			case '=':
				if (length > 0) {
					settings_rewind();
					goto END_OF_LOOP;
				} else if (c == '{') return TKN_BLOCK_START;
				else if (c=='}') return TKN_BLOCK_END;
				else val = 1;		/* Beginning with '=' means it's a value */
				break;
			case 0:	/* = End of file */
				return TKN_EOF;
			default:
				if ((c==' '||c=='\t')&&length>0) {
					goto END_OF_LOOP;
				} else if (c!=' '&&c!='\t') {
					buffer[length++] = c;
				}
				break;
		}
	} while (c != 0);

END_OF_LOOP:
	buffer[length] = 0;

	if (val) {
		tkn->type = TKN_VALUE;
		tkn->val = malloc(length+1);
		strcpy(tkn->val, buffer);
	} else if (!strcmp(buffer, "magnetic")) {
		tkn->type = TKN_MAGNETIC;
	} else if (!strcmp(buffer, "tool")) {
		tkn->type = TKN_TOOL;
	} else if (!strcmp(buffer, "particle")) {
		tkn->type = TKN_PARTICLE;
	} else if (!strcmp(buffer, "particles")) {
		tkn->type = TKN_PARTICLES;
	} else if (!strcmp(buffer, "sycout")) {
		tkn->type = TKN_SYCOUT;
	} else {
		tkn->type = TKN_NAME;
		tkn->val = malloc(length+1);
		strcpy(tkn->val, buffer);
	}

	return tkn->type;
}
void expect(enum settings_token_type type) {
	if (token() != type) {
		fprintf(stderr, "ERROR: line %d: %d: Expected %s!\n", settings_line, settings_charpos, settings_token_names[type]);
		exit(-1);
	}
}

void settings_interpret_value(particle *p, char *name, char *val) {
	if (!strcmp(name, "charge"))
		p->charge = atof(val)*CHARGE;
	else if (!strcmp(name, "mass"))
		p->mass = atof(val)*AMU_TO_KG;
	else if (!strcmp(name, "r0"))
		p->r0 = atodpn(val, 3, p->r0);
	else if (!strcmp(name, "t0"))
		p->t0 = atof(val);
	else if (!strcmp(name, "tend"))
		p->tend = atof(val);
	else if (!strcmp(name, "v0"))
		p->v0 = atodpn(val, 3, p->v0);
	else if (!strcmp(name, "vpar"))
		p->vpar = atof(val);
	else if (!strcmp(name, "vperp"))
		p->vperp = atof(val);
	else {
		fprintf(stderr, "ERROR: Unexpected particle setting: '%s'!\n", name);
		exit(-1);
	}

	free(tkn->val);
}
void settings_init_particle(particle *p) {
	p->t0 		= default_particle->t0;
	p->tend 	= default_particle->tend;
	p->mass		= default_particle->mass;
	p->charge	= default_particle->charge;
	p->vpar		= default_particle->vpar;
	p->vperp	= default_particle->vperp;
	p->gc_position = 1;

	if (default_particle->v0 != NULL) {
		p->v0 = malloc(sizeof(double)*3);
		p->v0[0] = default_particle->v0[0];
		p->v0[1] = default_particle->v0[1];
		p->v0[2] = default_particle->v0[2];
	} else p->v0 = malloc(sizeof(double)*3);

	if (default_particle->r0 != NULL) {
		p->r0 = malloc(sizeof(double)*3);
		p->r0[0] = default_particle->r0[0];
		p->r0[1] = default_particle->r0[1];
		p->r0[2] = default_particle->r0[2];
	} else p->r0 = malloc(sizeof(double)*3);
}

void settings_interpret_object(struct general_settings *obj) {
	expect(TKN_NAME);

	/* Add new tool setting */
	//set->sycouts = realloc(set->sycouts, sizeof(struct sycout_settings)*(set->nsycouts+1));
	obj->name = tkn->val;

	/* Expect { */
	expect(TKN_BLOCK_START);

	//int n = set->nsycouts++;
	obj->setting = NULL;
	obj->value = NULL;
	obj->n = 0;

	/* Read all values */
	while (!settings_eof && token() != TKN_BLOCK_END) {
		if (tkn->type != TKN_NAME) {
			fprintf(stderr, "Expected setting name!\n");
			exit(-1);
		}

		char *name = tkn->val;
		expect(TKN_VALUE);

		obj->setting = realloc(obj->setting, sizeof(char*)*(obj->n+1));
		obj->value   = realloc(obj->value,   sizeof(char*)*(obj->n+1));

		obj->setting[obj->n] = malloc(strlen(name)+1);
		obj->value[obj->n] = malloc(strlen(tkn->val)+1);
		strcpy(obj->setting[obj->n], name);
		strcpy(obj->value[obj->n], tkn->val);

		free(tkn->val);
		free(name);

		obj->n++;
	}
}

settings *settings_interpret() {
	settings *set = malloc(sizeof(settings));
	set->magfield = NULL;
	set->magnetic = NULL;
	set->nmagnetic = 0;
	set->tolerance = 0;
	set->maxtimestep = 0;
	set->equation = NULL;
	set->tool = NULL;
	set->tools = NULL;
	set->ntools = 0;
	set->sycouts = NULL;
	set->nsycouts = 0;
	set->particles = NULL;
	set->nparts = 0;
	set->warnoncollision = 0;
	//set->no_outer_wall = 0;
	set->threads = omp_get_max_threads();
	set->interptimestep = 0;
	set->nodrifts = 0;

	while (!settings_eof) {
		switch (token()) {
			case TKN_NAME: {
				char *name = tkn->val;
				expect(TKN_VALUE);

				if (!strcmp(name, "debug")) {
					DEBUG_OUTPUT = atoi(tkn->val);
				} else if (!strcmp(name, "domain_has_outer_wall")) {
					if (!strcmp(tkn->val, "no")) {
						printf("Removing outer wall of domain\n");
						domain_remove_outer_wall();
					} else if (strcmp(tkn->val, "yes")) {
						fprintf(stderr, "Unrecognized choice for 'domain_has_outer_wall'!\n");
						exit(-1);
					}
				} else if (!strcmp(name, "interptimestep")) {
					set->interptimestep = atof(tkn->val);
				} else if (!strcmp(name, "magnetic_field")) {
					set->magfield = tkn->val;
				} else if (!strcmp(name, "nodrifts")) {
					if (!strcmp(tkn->val, "no")) {
						set->nodrifts = 0;
					} else if (!strcmp(tkn->val, "yes")) {
						printf("Dropping drift terms from the guiding-center equations of motion\n");
						set->nodrifts = 1;
					}
				} else if (!strcmp(name, "threads")) {
					set->threads = atoi(tkn->val);
					if (set->threads <= 0) {
						printf("ERROR: Bad number of threads selected: %d\n", set->threads);
						exit(EXIT_FAILURE);
					}
					if (set->threads != omp_get_max_threads() && set->threads != 1)
						printf("WARNING: Optimal number of threads to use is %d\n", omp_get_max_threads());

				} else if (!strcmp(name, "tolerance")) {
					set->tolerance = atof(tkn->val);
				} else if (!strcmp(name, "maxtimestep")) {
					set->maxtimestep = atof(tkn->val);
				} else if (!strcmp(name, "useequation")) {
					set->equation = tkn->val;
				} else if (!strcmp(name, "usetool")) {
					set->tool = tkn->val;
				} else if (!strcmp(name, "warnoncollision")) {
					if (!strcmp(tkn->val, "yes"))
						set->warnoncollision = 1;
					else if (!strcmp(tkn->val, "no"))
						set->warnoncollision = 0;
					else
						printf("WARNING: Unrecognized value for 'warnoncollision': '%s'. Valid values are 'yes' or 'no'.\n", tkn->val);
				} else
					/* Else, we assume it's a particle setting */
					settings_interpret_value(default_particle, name, tkn->val);

				free(name);
			} break;
			case TKN_MAGNETIC: {
				set->magnetic = realloc(set->magnetic, sizeof(struct general_settings)*(set->nmagnetic+1));
				int n = set->nmagnetic++;
				settings_interpret_object(set->magnetic+n);
			} break;
			case TKN_SYCOUT: {
				/* Add new sycout setting */
				set->sycouts = realloc(set->sycouts, sizeof(struct general_settings)*(set->nsycouts+1));
				int n = set->nsycouts++;
				settings_interpret_object(set->sycouts+n);
			} break;
			case TKN_TOOL: {
				/* Add new tool setting */
				set->tools = realloc(set->tools, sizeof(struct general_settings)*(set->ntools+1));
				int n = set->ntools++;
				settings_interpret_object(set->tools+n);
			} break;
			case TKN_PARTICLE: {
				set->particles = realloc(set->particles, sizeof(particle)*(set->nparts+1));
				settings_init_particle(set->particles+set->nparts);

				expect(TKN_BLOCK_START);
				while (!settings_eof && token() != TKN_BLOCK_END) {
					if (tkn->type != TKN_NAME) {
						fprintf(stderr, "Expected setting name!\n");
						exit(-1);
					}

					char *name = tkn->val;
					expect(TKN_VALUE);
					settings_interpret_value(set->particles+set->nparts, name, tkn->val);

					free(name);
				}
				set->nparts++;
			} break;
			case TKN_PARTICLES: {
				set->particlespec = malloc(sizeof(struct particlespec));
				struct particlespec *p = set->particlespec;

				p->inputtype = PARTICLES_NONE;
				p->gentype = PARTICLES_GT_QUEUE;
				p->gc_position = 1;
				p->mass = 0.000548579909 * AMU_TO_KG;
				p->charge = -CHARGE;

				double db[3];

				expect(TKN_BLOCK_START);
				while (!settings_eof && token() != TKN_BLOCK_END) {
					if (tkn->type != TKN_NAME) {
						fprintf(stderr, "Expected setting name!\n");
						exit(-1);
					}

					char *name = tkn->val;
					expect(TKN_VALUE);

					if (!strcmp(name, "charge"))
						p->charge = atof(tkn->val)*CHARGE;
					else if (!strcmp(name, "cospitch")) {
						double *val = atodpn(tkn->val, 3, db);
						p->cospitch0 = val[0];
						p->cospitch1 = val[1];
						p->cospitchn = (int)val[2];
						p->inputtype |= PARTICLES_COSPITCH;
					} else if (!strcmp(name, "mass"))
						p->mass = atof(tkn->val)*AMU_TO_KG;
					else if (!strcmp(name, "gc_position")) {
						if (!strcmp(tkn->val, "no"))
							p->gc_position = 0;
						else if (strcmp(tkn->val, "yes")) {
							fprintf(stderr, "ERROR: Unrecognized choice for 'gc_position': %s\n", tkn->val);
							exit(-1);
						}
					} else if (!strcmp(name, "r")) {
						double *val = atodpn(tkn->val, 3, db);
						p->r0 = val[0];
						p->r1 = val[1];
						p->rn = (int)val[2];
						p->rcoord = PARTICLES_RC_MAJOR;
					} else if (!strcmp(name, "rdyn")) {
						double *val = atodpn(tkn->val, 2, db);
						p->rdynmax = val[0];
						p->rn = (int)val[1];
						p->rcoord = PARTICLES_RC_DYNAMIC;
					} else if (!strcmp(name, "t")) {
						double *val = atodpn(tkn->val, 2, db);
						p->t0 = val[0];
						p->tend = val[1];
					} else if (!strcmp(name, "p")) {
						double *val = atodpn(tkn->val, 3, db);
						p->p0 = val[0];
						p->p1 = val[1];
						p->pn = (int)val[2];
						p->inputtype |= PARTICLES_P;
					} else if (!strcmp(name, "pitch")) {
						double *val = atodpn(tkn->val, 3, db);
						p->pitch0 = val[0];
						p->pitch1 = val[1];
						p->pitchn = (int)val[2];
						p->inputtype |= PARTICLES_PITCH;
					} else if (!strcmp(name, "ppar")) {
						double *val = atodpn(tkn->val, 3, db);
						p->ppar0 = val[0];
						p->ppar1 = val[1];
						p->pparn = (int)val[2];
						p->inputtype |= PARTICLES_PPAR;
					} else if (!strcmp(name, "pperp")) {
						double *val = atodpn(tkn->val, 3, db);
						p->pperp0 = val[0];
						p->pperp1 = val[1];
						p->pperpn = (int)val[2];
						p->inputtype |= PARTICLES_PPERP;
					} else if (!strcmp(name, "generation")) {
						if (!strcmp(tkn->val, "queue")) {
							p->gentype = PARTICLES_GT_QUEUE;
						} else if (!strcmp(tkn->val, "even"))
							p->gentype = PARTICLES_GT_EVEN;
					} else {
						fprintf(stderr, "ERROR: Unexpected particle setting: '%s'!\n", name);
						exit(-1);
					}

					free(tkn->val);
					free(name);
				}

				if (p->mass <= 0.0) {
					fprintf(stderr, "ERROR: Invalid mass set: %e kg\n", p->mass);
					exit(-1);
				}
				if (p->charge == 0.0) {
					fprintf(stderr, "ERROR: The charge must be different from zero.\n");
					exit(-1);
				}
			} break;
			case TKN_EOF: break;
			default: fprintf(stderr, "ERROR: line %d: %d: Unexpected token at beginning of statement!\n", settings_line, settings_charpos); exit(-1);
		}
	}

	return set;
}

void settings_init(char *filename) {
	/* Open settings file */
	settings_file = fopen(filename, "r");

	if (!settings_file) {
		fprintf(stderr, "ERROR: Unable to open file '%s'!\n", filename);
		perror("ERROR");
		exit(-1);
	}

	/* Allocate token and particle */
	tkn = malloc(sizeof(settings_token));
	default_particle = malloc(sizeof(particle));

	/* Setup default particle */
	default_particle->t0 = 0;
	default_particle->tend = -1.0;
	default_particle->mass = 0.000548579909 * AMU_TO_KG;
	default_particle->charge = -CHARGE;
	default_particle->vpar = 0;
	default_particle->vperp = 0;
	default_particle->v0 = NULL;
	default_particle->r0 = NULL;
}

settings *load_settings(char *filename) {
	/* Initialize settings interpreter */
	settings_init(filename);

	/* Create settings object */
	settings *set = settings_interpret();

	/* Clean up */
	fclose(settings_file);

	/* Verify all required settings are set */
	char *e=NULL;
	if (set->magfield == NULL) e="No magnetic field set!";
	//else if (set->domain == NULL) e = "No domain set!";
	else if (set->equation == NULL) e = "No equation set!";
	else if (set->tool == NULL) e = "No tool set!";

	if (e != NULL) {
		fprintf(stderr, "ERROR: %s\n", e);
		exit(-1);
	}
	return set;
}

void settings_token_test(char *file) {
	settings_init(file);

	FILE *out;
	out = fopen("settings.txt", "w");

	if (!out) {
		fprintf(stderr, "ERROR: Unable to create file settings.txt!\n");
		perror("ERROR");
		exit(-1);
	}

	while (!settings_eof) {
		switch (token()) {
			case TKN_NAME: fprintf(out, "%s", tkn->val); break;
			case TKN_VALUE: fprintf(out, "=%s\n", tkn->val); break;
			case TKN_MAGNETIC: fprintf(out, "magnetic "); break;
			case TKN_TOOL: fprintf(out, "tool "); break;
			case TKN_PARTICLE: fprintf(out, "particle "); break;
			case TKN_PARTICLES: fprintf(out, "particles "); break;
			case TKN_BLOCK_START: fprintf(out, "{\n"); break;
			case TKN_BLOCK_END: fprintf(out, "}\n"); break;
			case TKN_EOF: break;
			default: fprintf(stderr, "Invalid token!\n"); exit(-1);
		}
	}

	fclose(out);
}
