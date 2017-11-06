/* elem.h - element header file */

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "gdf.h"

#include "gps.h"
#include "compat.h"
#include "axis.h"
#include "init.h"
#include "sim.h"
#include "memman.h"
#include "misc.h"
#include "odeman.h"
#include "output.h"
#include "input.h"
#include "parman.h"
#include "parse.h"
#include "poisson.h"
#include "scatter.h"
#include "dists.h"

#define X  (par->r[0])
#define Y  (par->r[1])
#define Z  (par->r[2])
#define BX (par->B[0])
#define BY (par->B[1])
#define BZ (par->B[2])
#define EX (par->E[0])
#define EY (par->E[1])
#define EZ (par->E[2])
