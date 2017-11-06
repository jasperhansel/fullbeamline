## Copyright (C) 2003 John W. Eaton
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
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, write to the Free
## Software Foundation, 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

## -*- texinfo -*-
## @deftypefn {Function File} {} struct_contains (@var{expr}, @var{name})
## This function has been deprecated.  Use isfield instead.
## @end deftypefn

## Author: jwe

## Changed by Rafael Laboissiere for use with the octave-plplot package
## on Sat May  8 16:42:53 UTC 2004

function retval = struct_contains (varargin)

  retval = isstruct (varargin{1}) && isfield (varargin{:});

endfunction
