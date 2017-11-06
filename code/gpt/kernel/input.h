/* input.h - Header file for (GDF) input functions */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

double *gptinputdoublearray(struct gdfdata *ds, char *name, int *len) ;
double *gptinputdoublearraypoints(struct gdfdata *ds, char *name, int points) ;
double *gptinputdoublearraypointsnull(struct gdfdata *ds, char *name, int points) ;
struct gdfdata *gptinputgetgroup(struct gdfmem *gm, char *filename, char *name) ;
