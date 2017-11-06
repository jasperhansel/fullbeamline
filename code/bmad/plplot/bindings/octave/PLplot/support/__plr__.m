## Copyright (C) 1996 John W. Eaton
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.

## Author: jwe
## Modified: jc

function __plr__ (theta, rho, fmt)

  if (nargin == 1)
    [nr, nc] = size (theta);
    if (nr == 1)
      theta = theta.';
      tmp = nr;
      nr = nc;
      nc = tmp;
    endif
    theta_i = imag (theta);
    if (any (theta_i))
      rho = theta_i;
      theta = real (theta);
    else
      rho = theta;
      theta = (1:nr)';
    endif
  endif

  if (nargin <= 2)
    if (any (imag (theta)))
      theta = real (theta);
    endif
    if (any (imag (rho)))
      rho = real (rho);
    endif
    if (isscalar (theta))
      if (isscalar (rho))
        x = rho * cos (theta);
        y = rho * sin (theta);
        __plt2ss__ (x, y, fmt);
      endif
    elseif (isvector (theta))
      if (isvector (rho))
        if (length (theta) != length (rho))
          error ("polar: vector lengths must match");
        endif
        if (rows (rho) == 1)
          rho = rho.';
        endif
        if (rows (theta) == 1)
          theta = theta.';
        endif
        x = rho .* cos (theta);
        y = rho .* sin (theta);
        __plt2vv__ (x, y, fmt);
      elseif (ismatrix (rho))
        [t_nr, t_nc] = size (theta);
        if (t_nr == 1)
          theta = theta.';
          tmp = t_nr;
          t_nr = t_nc;
          t_nc = tmp;
        endif
        [r_nr, r_nc] = size (rho);
        if (t_nr != r_nr)
          rho = rho.';
          tmp = r_nr;
          r_nr = r_nc;
          r_nc = tmp;
        endif
        if (t_nr != r_nr)
          error ("polar: vector and matrix sizes must match");
        endif
        x = diag (cos (theta)) * rho;
        y = diag (sin (theta)) * rho;
        __plt2vm__ (x, y, fmt);
      endif
    elseif (ismatrix (theta))
      if (isvector (rho))
        [r_nr, r_nc] = size (rho);
        if (r_nr == 1)
          rho = rho.';
          tmp = r_nr;
          r_nr = r_nc;
          r_nc = tmp;
        endif
        [t_nr, t_nc] = size (theta);
        if (r_nr != t_nr)
          theta = rho.';
          tmp = t_nr;
          t_nr = t_nc;
          t_nc = tmp;
        endif
        if (r_nr != t_nr)
          error ("polar: vector and matrix sizes must match");
        endif
        diag_r = diag (r);
        x = diag_r * cos (theta);
        y = diag_r * sin (theta);
        __plt2mv__ (x, y, fmt);
      elseif (ismatrix (rho))
        if (size (rho) != size (theta))
          error ("polar: matrix dimensions must match");
        endif
        x = rho .* cos (theta);
        y = rho .* sin (theta);
        __plt2mm__ (x, y, fmt);
      endif
    endif
  else
    usage ("__plr__ (x [, y])");
  endif

endfunction
