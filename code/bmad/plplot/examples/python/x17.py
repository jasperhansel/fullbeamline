#!/usr/bin/env python

# Copyright 2002 Gary Bishop
# Copyright 2004-2016 Alan W. Irwin

# This file is part of PLplot.

# PLplot is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; version 2 of the License.

# PLplot is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.

# You should have received a copy of the GNU Library General Public License
# along with the file PLplot; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

# Plots a simple stripchart with four pens.

# Append to effective python path so that can find plplot modules.
from plplot_python_start import *

import sys
import plplot as w
from numpy import *

# Parse and process command line arguments
w.plparseopts(sys.argv, w.PL_PARSE_FULL)

# Initialize plplot
w.plinit()

from time import sleep

def main(w):
    nsteps = 1000

# If db is used the plot is much more smooth. However, because of the
#   async X behaviour, one does not have a real-time scripcharter.

#    w.plsetopt("db", "")
#    w.plsetopt("np", "")

# User sets up plot completely except for window and data
# Eventually settings in place when strip chart is created will be
# remembered so that multiple strip charts can be used simultaneously.
#

# Specify some reasonable defaults for ymin and ymax
# The plot will grow automatically if needed (but not shrink)

    ymin = -0.1
    ymax = 0.1

# Specify initial tmin and tmax -- this determines length of window.
# Also specify maximum jump in t
# This can accomodate adaptive timesteps

    tmin = 0.
    tmax = 10.
    tjump = 0.3	# percentage of plot to jump

# Axes options same as w.plbox.
# Only automatic tick generation and label placement allowed
# Eventually I ll make this fancier

    colbox = 1
    collab = 3
    styline = [2, 3, 4, 5]
    colline = [2, 3, 4, 5]

    legline = ["sum", "sin", "sin*noi", "sin+noi"]

    xlab = 0.
    ylab = 0.25	# legend position

    autoy = 1	# autoscale y
    acc = 1	# don t scrip, accumulate

    w.pladv(0)
    w.plvsta()

# Register our error variables with PLplot
# From here on, we're handling all errors here

    #w.plsError(&pl_errcode, errmsg)

    id1 = w.plstripc("bcnst", "bcnstv",
                   tmin, tmax, tjump, ymin, ymax,
                   xlab, ylab,
                   autoy, acc,
                   colbox, collab,
                   colline, styline, legline,
                   "t", "", "Strip chart demo")

# Let plplot handle errors from here on

    #w.plsError(NULL, NULL)

    autoy = 0	# autoscale y
    acc = 1	# accumulate

# This is to represent a loop over time
# Let's try a random walk process

    y1 = y2 = y3 = y4 = 0.0
    dt = 0.1

    for n in range(nsteps):
        sleep(0.01)
	t = n * dt
	noise = w.plrandd() - 0.5
	y1 = y1 + noise
	y2 = sin(t*pi/18.)
	y3 = y2 * noise
	y4 = y2 + noise/3.

        # There is no need for all pens to have the same number of
        # points or beeing equally time spaced.
		
        if n%2:	
	    w.plstripa(id1, 0, t, y1)
	if n%3:
	    w.plstripa(id1, 1, t, y2)
	if n%4:
	    w.plstripa(id1, 2, t, y3)
	if n%5:
	    w.plstripa(id1, 3, t, y4)

    # Destroy strip chart and it's memory

    w.plstripd(id1)

    # Restore defaults
    # No defaults changed so nothing to restore

main(w)
w.plend()
