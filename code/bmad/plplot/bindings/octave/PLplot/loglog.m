## Copyright (C) 1998-2003 Joao Cardoso.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## This file is part of plplot_octave.

## usage: loglog (x, y)
##        loglog (x1, y1, x2, y2, ...)
##        loglog (x, y, fmt)
##
## Make a 2D plot of y versus x using a log scale for the x axis.
##
## See the help message for the plot command for a description of how
## the arguments are interpreted.

function loglog (varargin)

  __plt__ ("lgxy", varargin{:});

endfunction
