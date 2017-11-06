# Copyright 2002 Gary Bishop
# Copyright 2002-2016 Alan W. Irwin
# Copyright 2008 Andrew Ross
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

# Wrap raw python interface to C API, plplotc, with this user-friendly version
# which implements some useful variations of the argument lists.

from plplotc import *
import types
import numpy

# Redefine plcont to have the user-friendly interface
# Allowable syntaxes:

# plcont( z, [kx, lx, ky, ly], clev, [pltr, [pltr_data] or [xg, yg, [wrap]]])
# N.B. Brackets represent options here and not python lists!

# All unbracketed arguments within brackets must all be present or all be
# missing.  Furthermore, z must be a 2D array, kx, lx, ky, ly must all be
# integers, clev must be a 1D array, pltr can be a function reference or
# string, pltr_data is an optional arbitrary data object, xg and yg are
# optional 1D or 2D arrays and wrap (which only works if xg and yg
# are specified) is 0, 1, or 2.

# If pltr is a string it must be either "pltr0", "pltr1", or "pltr2" to
# refer to those built-in transformation functions.  Alternatively, the
# function names pltr0, pltr1, or pltr2 may be specified to refer to
# the built-in transformation functions or an arbitrary name for a
# user-defined transformation function may be specified.  Such functions
# must have x, y, and optional pltr_data arguments and return arbitrarily
# transformed x' and y' in a tuple.  The built-in pltr's such as pltr1 and
# pltr2 use pltr_data = tuple(xg, yg), and for this oft-used case (and any
# other user-defined pltr which uses a tuple of two arrays for pltr_data),
# we also provide optional xg and yg arguments separately as an alternative
# to the tuple method of providing these data. Note, that pltr_data cannot
# be in the argument list if xg and yg are there, and vice versa. Also note
# that the built-in pltr0 and some user-defined transformation functions
# ignore the auxiliary pltr_data (or the alternative xg and yg) in which
# case neither pltr_data nor xg and yg need to be specified.

_plcont = plcont
def plcont(z, *args):
    z = numpy.asarray(z)
    if len(z.shape) != 2:
	raise ValueError, "Expected 2D z array"

    if len(args) > 4 and type(args[0]) == types.IntType:
	for i in range(1,4):
	    if type(args[i]) != types.IntType:
		raise ValueError, "Expected 4 ints for kx,lx,ky,ly"

	else:
	    # these 4 args are the kx, lx, ky, ly ints
	    ifdefault_range = 0
	    kx,lx,ky,ly = args[0:4]
	    args = args[4:]
    else:
	ifdefault_range = 1

    if len(args) > 0:
	clev = numpy.asarray(args[0])
	if len(clev.shape) !=1:
	    raise ValueError, "Expected 1D clev array"
	args = args[1:]
    else:
	raise ValueError, "Missing clev argument"

    if len(args) > 0 and ( \
    type(args[0]) == types.StringType or \
    type(args[0]) == types.FunctionType or \
    type(args[0]) == types.BuiltinFunctionType):
	pltr = args[0]
	# Handle the string names for the callbacks though specifying the
	# built-in function name directly (without the surrounding quotes)
	# or specifying any user-defined transformation function
	# (following above rules) works fine too.
	if type(pltr) == types.StringType:
	    if pltr == "pltr0":
		pltr = pltr0
	    elif pltr == "pltr1":
		pltr = pltr1
	    elif pltr == "pltr2":
		pltr = pltr2
	    else:
		raise ValueError, "pltr string is unrecognized"

	args = args[1:]
	# Handle pltr_data or separate xg, yg, [wrap]
	if len(args) == 0:
	    # Default pltr_data
	    pltr_data = None
	elif len(args) == 1:
	    #Must be pltr_data
	    pltr_data = args[0]
	    args = args[1:]
	elif len(args) >= 2:
	    xg = numpy.asarray(args[0])
	    if len(xg.shape) < 1 or len(xg.shape) > 2:
		raise ValueError, "xg must be 1D or 2D array"
	    yg = numpy.asarray(args[1])
	    if len(yg.shape) != len(xg.shape):
		raise ValueError, "yg must have same number of dimensions as xg"
	    args = args[2:]
	    # wrap only relevant if xg and yg specified.
	    if len(args) > 0:
	     if type(args[0]) == types.IntType:
	      wrap = args[0]
	      args = args[1:]
	      if len(xg.shape) == 2 and len(yg.shape) == 2 and  \
	      z.shape == xg.shape and z.shape == yg.shape:
		# handle wrap
		if wrap == 1:
		    z = numpy.resize(z, (z.shape[0]+1, z.shape[1]))
		    xg = numpy.resize(xg, (xg.shape[0]+1, xg.shape[1]))
		    yg = numpy.resize(yg, (yg.shape[0]+1, yg.shape[1]))
		elif wrap == 2:
		    z = numpy.transpose(numpy.resize( \
		    numpy.transpose(z), (z.shape[1]+1, z.shape[0])))
		    xg = numpy.transpose(numpy.resize( \
		    numpy.transpose(xg), (xg.shape[1]+1, xg.shape[0])))
		    yg = numpy.transpose(numpy.resize( \
		    numpy.transpose(yg), (yg.shape[1]+1, yg.shape[0])))
		elif wrap != 0:
		    raise ValueError, "Invalid wrap specifier, must be 0, 1 or 2."
	      elif wrap != 0:
		  raise ValueError, "Non-zero wrap specified and xg and yg are not 2D arrays"
	     else:
		 raise ValueError, "Specified wrap is not an integer"
	    pltr_data = (xg, yg)
    else:
	# default is identity transformation
	pltr = pltr0
	pltr_data = None
    if len(args) > 0:
	raise ValueError, "Too many arguments for plcont"
    if ifdefault_range:
	# Default is to take full range (still using fortran convention
	# for indices which is embedded in the PLplot library API)
	kx = 1
	lx = z.shape[0]
	ky = 1
	ly = z.shape[1]
    _plcont(z, kx, lx, ky, ly, clev, pltr, pltr_data)
plcont.__doc__ = _plcont.__doc__

# Redefine plvect to have the user-friendly interface
# Allowable syntaxes:

# plvect( u, v, scaling, [pltr, [pltr_data] or [xg, yg, [wrap]]])
_plvect = plvect
def plvect(u, v, *args):
    u = numpy.asarray(u)
    v = numpy.asarray(v)

    if len(u.shape) != 2:
        raise ValueError, "Expected 2D u array"
    if len(v.shape) != 2:
        raise ValueError, "Expected 2D v array"
    if (u.shape[0] != v.shape[0]) or (u.shape[1] != v.shape[1]) :
        raise ValueError, "Expected u and v arrays to be the same dimensions"

    if len(args) > 0 and (type(args[0]) == types.FloatType or type(args[0]) == numpy.float64) :
        scaling = args[0]
        args = args[1:]
    else:
        raise ValueError, "Missing scaling argument"

    if len(args) > 0 and ( \
    type(args[0]) == types.StringType or \
    type(args[0]) == types.FunctionType or \
    type(args[0]) == types.BuiltinFunctionType):
        pltr = args[0]
        # Handle the string names for the callbacks though specifying the
        # built-in function name directly (without the surrounding quotes)
        # or specifying any user-defined transformation function
        # (following above rules) works fine too.
        if type(pltr) == types.StringType:
            if pltr == "pltr0":
                pltr = pltr0
            elif pltr == "pltr1":
                pltr = pltr1
            elif pltr == "pltr2":
                pltr = pltr2
            else:
                raise ValueError, "pltr string is unrecognized"

        args = args[1:]
        # Handle pltr_data or separate xg, yg, [wrap]
        if len(args) == 0:
            # Default pltr_data
            pltr_data = None
        elif len(args) == 1:
            #Must be pltr_data
            pltr_data = args[0]
            args = args[1:]
        elif len(args) >= 2:
            xg = numpy.asarray(args[0])
            if len(xg.shape) < 1 or len(xg.shape) > 2:
                raise ValueError, "xg must be 1D or 2D array"
            yg = numpy.asarray(args[1])
            if len(yg.shape) != len(xg.shape):
                raise ValueError, "yg must have same number of dimensions as xg"
            args = args[2:]
            # wrap only relevant if xg and yg specified.
            if len(args) > 0:
             if type(args[0]) == types.IntType:
              wrap = args[0]
              args = args[1:]
              if len(xg.shape) == 2 and len(yg.shape) == 2 and  \
              u.shape == xg.shape and u.shape == yg.shape:
                # handle wrap
                if wrap == 1:
                    u = numpy.resize(u, (u.shape[0]+1, u.shape[1]))
                    v = numpy.resize(v, (v.shape[0]+1, v.shape[1]))
                    xg = numpy.resize(xg, (xg.shape[0]+1, xg.shape[1]))
                    yg = numpy.resize(yg, (yg.shape[0]+1, yg.shape[1]))
                elif wrap == 2:
                    u = numpy.transpose(numpy.resize( \
                    numpy.transpose(u), (u.shape[1]+1, u.shape[0])))
                    v = numpy.transpose(numpy.resize( \
                    numpy.transpose(v), (v.shape[1]+1, v.shape[0])))
                    xg = numpy.transpose(numpy.resize( \
                    numpy.transpose(xg), (xg.shape[1]+1, xg.shape[0])))
                    yg = numpy.transpose(numpy.resize( \
                    numpy.transpose(yg), (yg.shape[1]+1, yg.shape[0])))
                elif wrap != 0:
                    raise ValueError, "Invalid wrap specifier, must be 0, 1 or 2."
              elif wrap != 0:
                  raise ValueError, "Non-zero wrap specified and xg and yg are not 2D arrays"
             else:
                 raise ValueError, "Specified wrap is not an integer"
            pltr_data = (xg, yg)
    else:
        # default is identity transformation
        pltr = pltr0
        pltr_data = None
    if len(args) > 0:
        raise ValueError, "Too many arguments for plvect"
    _plvect(u, v, scaling, pltr, pltr_data)
plvect.__doc__ = _plvect.__doc__

# Redefine plimagefr to have the user-friendly interface
# Allowable syntaxes:

# plimagefr( img, xmin, xmax, ymin, ymax, zmin, zmax, valuemin, valuemax, [pltr, [pltr_data] or [xg, yg, [wrap]]])
_plimagefr = plimagefr
def plimagefr(img, *args):
    img = numpy.asarray(img)

    if len(img.shape) != 2:
        raise ValueError, "Expected 2D img array"

    if len(args) >= 8 :
        for i in range(8) :
            if (type(args[i]) != types.FloatType and \
                type(args[i]) != numpy.float64 and \
                type(args[i]) != types.IntType) :
                raise ValueError, "Expected 8 numbers for xmin, xmax, ymin, ymax, zmin, zmax, valuemin, valuemax"
        else:
            # These 8 args are xmin, xmax, ymin, ymax, zmin, zmax, valuemin, valuemax
            xmin, xmax, ymin, ymax, zmin, zmax, valuemin, valuemax = args[0:8]
            args = args[8:]
    else:
        raise ValueError, "Expected 8 numbers for xmin, xmax, ymin, ymax, zmin, zmax, valuemin, valuemax"

    if len(args) > 0 and ( \
    type(args[0]) == types.StringType or \
    type(args[0]) == types.FunctionType or \
    type(args[0]) == types.BuiltinFunctionType):
        pltr = args[0]
        # Handle the string names for the callbacks though specifying the
        # built-in function name directly (without the surrounding quotes)
        # or specifying any user-defined transformation function
        # (following above rules) works fine too.
        if type(pltr) == types.StringType:
            if pltr == "pltr0":
                pltr = pltr0
            elif pltr == "pltr1":
                pltr = pltr1
            elif pltr == "pltr2":
                pltr = pltr2
            else:
                raise ValueError, "pltr string is unrecognized"

        args = args[1:]
        # Handle pltr_data or separate xg, yg, [wrap]
        if len(args) == 0:
            # Default pltr_data
            pltr_data = None
        elif len(args) == 1:
            #Must be pltr_data
            pltr_data = args[0]
            args = args[1:]
        elif len(args) >= 2:
            xg = numpy.asarray(args[0])
            if len(xg.shape) < 1 or len(xg.shape) > 2:
                raise ValueError, "xg must be 1D or 2D array"
            yg = numpy.asarray(args[1])
            if len(yg.shape) != len(xg.shape):
                raise ValueError, "yg must have same number of dimensions as xg"
            args = args[2:]
            # wrap only relevant if xg and yg specified.
            if len(args) > 0:
             if type(args[0]) == types.IntType:
              wrap = args[0]
              args = args[1:]
              if len(xg.shape) == 2 and len(yg.shape) == 2 and  \
              img.shape[0] == xg.shape[0]-1 and img.shape[1] == xg.shape[1]-1:
                # handle wrap
                if wrap == 1:
                    img = numpy.resize(img, (img.shape[0]+1, u.shape[1]))
                    xg = numpy.resize(xg, (xg.shape[0]+1, xg.shape[1]))
                    yg = numpy.resize(yg, (yg.shape[0]+1, yg.shape[1]))
                elif wrap == 2:
                    img = numpy.transpose(numpy.resize( \
                    numpy.transpose(img), (img.shape[1]+1, img.shape[0])))
                    xg = numpy.transpose(numpy.resize( \
                    numpy.transpose(xg), (xg.shape[1]+1, xg.shape[0])))
                    yg = numpy.transpose(numpy.resize( \
                    numpy.transpose(yg), (yg.shape[1]+1, yg.shape[0])))
                elif wrap != 0:
                    raise ValueError, "Invalid wrap specifier, must be 0, 1 or 2."
              elif wrap != 0:
                  raise ValueError, "Non-zero wrap specified and xg and yg are not 2D arrays"
             else:
                 raise ValueError, "Specified wrap is not an integer"
            pltr_data = (xg, yg)
    else:
        # default is NULL callback arguments.
        pltr = None
        pltr_data = None
    if len(args) > 0:
        raise ValueError, "Too many arguments for plimagefr"
    _plimagefr(img, xmin, xmax, ymin, ymax, zmin, zmax, valuemin, valuemax, pltr, pltr_data)
plimagefr.__doc__ = _plimagefr.__doc__

# Redefine plshades to have the user-friendly interface
# Allowable syntaxes:

# plshades(z,  [xmin, xmax, ymin, ymax,] clev, \
# fill_width, [cont_color, cont_width,], rect, \
# [pltr, [pltr_data] or [xg, yg, [wrap]]])

_plshades = plshades
def plshades(z, *args):
    z = numpy.asarray(z)
    if len(z.shape) != 2:
	raise ValueError, "Expected 2D z array"

    if len(args) > 4 and \
    (type(args[0]) == types.FloatType or type(args[0]) == numpy.float64 or type(args[0]) == types.IntType) and \
    (type(args[1]) == types.FloatType or type(args[1]) == numpy.float64 or type(args[1]) == types.IntType) and \
    (type(args[2]) == types.FloatType or type(args[2]) == numpy.float64 or type(args[2]) == types.IntType) and \
    (type(args[3]) == types.FloatType or type(args[3]) == numpy.float64 or type(args[3]) == types.IntType):
	# These 4 args are xmin, xmax, ymin, ymax
	xmin, xmax, ymin, ymax = args[0:4]
	args = args[4:]
    else:
	# These values are ignored if pltr and pltr_data are defined in any case.
	# So pick some convenient defaults that work for the pltr0, None case
	xmin = -1.
	xmax = 1.
	ymin = -1.
	ymax = 1.

    # clev must be present.
    if len(args) > 0:
	clev = numpy.asarray(args[0])
	if len(clev.shape) !=1:
	    raise ValueError, "Expected 1D clev array"
	args = args[1:]
    else:
	raise ValueError, "Missing clev argument"

    # fill_width must be present
    if len(args) > 0 and (type(args[0]) == types.FloatType or type(args[0]) == numpy.float64):
	fill_width = args[0]
	args = args[1:]
    else:
	raise ValueError, "fill_width argument must be present and of types.FloatType or numpy.float64 type"

    # cont_color and cont_width are optional.
    if len(args) > 2 and \
    type(args[0]) == types.IntType and \
    (type(args[1]) == types.FloatType or type(args[1]) == numpy.float64):
	# These 2 args are
	cont_color, cont_width = args[0:2]
	args = args[2:]
    else:
	# Turn off contouring.
	cont_color, cont_width = (0,0.)

    # rect must be present.
    if len(args) > 0 and type(args[0]) == types.IntType:
	rect = args[0]
	args = args[1:]
    else:
	raise ValueError, "Missing rect argument"

    if len(args) > 0 and ( \
    type(args[0]) == types.NoneType or \
    type(args[0]) == types.StringType or \
    type(args[0]) == types.FunctionType or \
    type(args[0]) == types.BuiltinFunctionType):
	pltr = args[0]
	# Handle the string names for the callbacks though specifying the
	# built-in function name directly (without the surrounding quotes)
	# or specifying any user-defined transformation function
	# (following above rules) works fine too.
	if type(pltr) == types.StringType:
	    if pltr == "pltr0":
		pltr = pltr0
	    elif pltr == "pltr1":
		pltr = pltr1
	    elif pltr == "pltr2":
		pltr = pltr2
	    else:
		raise ValueError, "pltr string is unrecognized"

	args = args[1:]
	# Handle pltr_data or separate xg, yg, [wrap]
	if len(args) == 0:
	    # Default pltr_data
	    pltr_data = None
	elif len(args) == 1:
	    #Must be pltr_data
	    pltr_data = args[0]
	    args = args[1:]
	elif len(args) >= 2:
	    xg = numpy.asarray(args[0])
	    if len(xg.shape) < 1 or len(xg.shape) > 2:
		raise ValueError, "xg must be 1D or 2D array"
	    yg = numpy.asarray(args[1])
	    if len(yg.shape) != len(xg.shape):
		raise ValueError, "yg must have same number of dimensions as xg"
	    args = args[2:]
	    # wrap only relevant if xg and yg specified.
	    if len(args) > 0:
	     if type(args[0]) == types.IntType:
	      wrap = args[0]
	      args = args[1:]
	      if len(xg.shape) == 2 and len(yg.shape) == 2 and  \
	      z.shape == xg.shape and z.shape == yg.shape:
		# handle wrap
		if wrap == 1:
		    z = numpy.resize(z, (z.shape[0]+1, z.shape[1]))
		    xg = numpy.resize(xg, (xg.shape[0]+1, xg.shape[1]))
		    yg = numpy.resize(yg, (yg.shape[0]+1, yg.shape[1]))
		elif wrap == 2:
		    z = numpy.transpose(numpy.resize( \
		    numpy.transpose(z), (z.shape[1]+1, z.shape[0])))
		    xg = numpy.transpose(numpy.resize( \
		    numpy.transpose(xg), (xg.shape[1]+1, xg.shape[0])))
		    yg = numpy.transpose(numpy.resize( \
		    numpy.transpose(yg), (yg.shape[1]+1, yg.shape[0])))
		elif wrap != 0:
		    raise ValueError, "Invalid wrap specifier, must be 0, 1 or 2."
	      elif wrap != 0:
		  raise ValueError, "Non-zero wrap specified and xg and yg are not 2D arrays"
	     else:
		 raise ValueError, "Specified wrap is not an integer"
	    pltr_data = (xg, yg)
    else:
	# default is identity transformation
	pltr = pltr0
	pltr_data = None
    if len(args) > 0:
	raise ValueError, "Too many arguments for plshades"

    _plshades(z, xmin, xmax, ymin, ymax, clev, \
    fill_width, cont_color, cont_width, rect, pltr, pltr_data)
plshades.__doc__ = _plshades.__doc__

# Redefine plshade to have the user-friendly interface
# Allowable syntaxes:

# _plshade(z, [xmin, xmax, ymin, ymax,] \
# shade_min, shade_max, sh_cmap, sh_color, sh_width, \
# [min_color, min_width, max_color, max_width,] rect, \
# [pltr, [pltr_data] or [xg, yg, [wrap]]])

_plshade = plshade
def plshade(z, *args):
    z = numpy.asarray(z)
    if len(z.shape) != 2:
	raise ValueError, "Expected 2D z array"

    # Extra check on shade_min = float on end is absolutely necessary
    # to unambiguously figure out where we are in the argument list.
    if len(args) > 9 and \
    (type(args[0]) == types.FloatType or type(args[0]) == numpy.float64 or type(args[0]) == types.IntType) and \
    (type(args[1]) == types.FloatType or type(args[1]) == numpy.float64 or type(args[1]) == types.IntType) and \
    (type(args[2]) == types.FloatType or type(args[2]) == numpy.float64 or type(args[2]) == types.IntType) and \
    (type(args[3]) == types.FloatType or type(args[3]) == numpy.float64 or type(args[3]) == types.IntType) and \
    (type(args[4]) == types.FloatType or type(args[4]) == numpy.float64) :
	# These 4 args are xmin, xmax, ymin, ymax
	xmin, xmax, ymin, ymax = args[0:4]
	args = args[4:]
    else:
	# These values are ignored if pltr and pltr_data are defined in any case.
	# So pick some convenient defaults that work for the pltr0, None case
	xmin = -1.
	xmax = 1.
	ymin = -1.
	ymax = 1.

    # shade_min, shade_max, sh_cmap, sh_color, sh_width, must be present.
    # sh_color can be either integer or float.
    if len(args) > 5 and \
    (type(args[0]) == types.FloatType or type(args[0]) == numpy.float64) and \
    (type(args[1]) == types.FloatType or type(args[1]) == numpy.float64) and \
    type(args[2]) == types.IntType and \
    (type(args[3]) == types.FloatType or type(args[3]) == numpy.float64 or type(args[3]) == types.IntType) and \
    (type(args[4]) == types.FloatType or type(args[4]) == numpy.float64):
	shade_min, shade_max, sh_cmap, sh_color, sh_width = args[0:5]
	args = args[5:]
    else:
	raise ValueError, \
	"shade_min, shade_max, sh_cmap, sh_color, sh_width, must be present with sh_cmap of types.IntType type and the rest of types.FloatType or numpy.float64 type"

    # min_color, min_width, max_color, max_width are optional.
    if len(args) > 4 and \
    type(args[0]) == types.IntType and \
    (type(args[1]) == types.FloatType or type(args[1]) == numpy.float64) and \
    type(args[2]) == types.IntType and \
    (type(args[3]) == types.FloatType or type(args[3]) == numpy.float64):
	# These 4 args are
	min_color, min_width, max_color, max_width = args[0:4]
	args = args[4:]
    else:
	# Turn off boundary colouring
	min_color, min_width, max_color, max_width = (0,0.,0,0.)

    # rect must be present.
    if len(args) > 0 and type(args[0]) == types.IntType:
	rect = args[0]
	args = args[1:]
    else:
	raise ValueError, "Missing rect argument"

    if len(args) > 0 and ( \
    type(args[0]) == types.NoneType or \
    type(args[0]) == types.StringType or \
    type(args[0]) == types.FunctionType or \
    type(args[0]) == types.BuiltinFunctionType):
	pltr = args[0]
	# Handle the string names for the callbacks though specifying the
	# built-in function name directly (without the surrounding quotes)
	# or specifying any user-defined transformation function
	# (following above rules) works fine too.
	if type(pltr) == types.StringType:
	    if pltr == "pltr0":
		pltr = pltr0
	    elif pltr == "pltr1":
		pltr = pltr1
	    elif pltr == "pltr2":
		pltr = pltr2
	    else:
		raise ValueError, "pltr string is unrecognized"

	args = args[1:]
	# Handle pltr_data or separate xg, yg, [wrap]
	if len(args) == 0:
	    # Default pltr_data
	    pltr_data = None
	elif len(args) == 1:
	    #Must be pltr_data
	    pltr_data = args[0]
	    args = args[1:]
	elif len(args) >= 2:
	    xg = numpy.asarray(args[0])
	    if len(xg.shape) < 1 or len(xg.shape) > 2:
		raise ValueError, "xg must be 1D or 2D array"
	    yg = numpy.asarray(args[1])
	    if len(yg.shape) != len(xg.shape):
		raise ValueError, "yg must have same number of dimensions as xg"
	    args = args[2:]
	    # wrap only relevant if xg and yg specified.
	    if len(args) > 0:
	     if type(args[0]) == types.IntType:
	      wrap = args[0]
	      args = args[1:]
	      if len(xg.shape) == 2 and len(yg.shape) == 2 and  \
	      z.shape == xg.shape and z.shape == yg.shape:
		# handle wrap
		if wrap == 1:
		    z = numpy.resize(z, (z.shape[0]+1, z.shape[1]))
		    xg = numpy.resize(xg, (xg.shape[0]+1, xg.shape[1]))
		    yg = numpy.resize(yg, (yg.shape[0]+1, yg.shape[1]))
		elif wrap == 2:
		    z = numpy.transpose(numpy.resize( \
		    numpy.transpose(z), (z.shape[1]+1, z.shape[0])))
		    xg = numpy.transpose(numpy.resize( \
		    numpy.transpose(xg), (xg.shape[1]+1, xg.shape[0])))
		    yg = numpy.transpose(numpy.resize( \
		    numpy.transpose(yg), (yg.shape[1]+1, yg.shape[0])))
		elif wrap != 0:
		    raise ValueError, "Invalid wrap specifier, must be 0, 1 or 2."
	      elif wrap != 0:
		  raise ValueError, "Non-zero wrap specified and xg and yg are not 2D arrays"
	     else:
		 raise ValueError, "Specified wrap is not an integer"
	    pltr_data = (xg, yg)
    else:
	# default is identity transformation
	pltr = pltr0
	pltr_data = None
    if len(args) > 0:
	raise ValueError, "Too many arguments for plshade"

    _plshade(z,  xmin, xmax, ymin, ymax, \
    shade_min, shade_max, sh_cmap, sh_color, sh_width, \
    min_color, min_width, max_color, max_width, rect, pltr, pltr_data)
plshade.__doc__ = _plshade.__doc__

# Redefine plscmap1l to have the user-friendly interface
# Allowable syntaxes:

# plscmap1l(itype, pos, coord1, coord2, coord3[, alt_hue_path])

_plscmap1l = plscmap1l
def plscmap1l(itype, pos, coord1, coord2, coord3, *args):

    pos = numpy.asarray(pos)
    if len(pos.shape) != 1:
	raise ValueError, "Expected 1D pos array"
		
    if len(args) == 0:
	# Default alt_hue_path
	alt_hue_path = numpy.zeros(pos.shape[0]-1,dtype="int")
    elif len(args) == 1:
	alt_hue_path = numpy.asarray(args[0])
    else:
	raise ValueError, "Too many arguments to plscmap1l"
    _plscmap1l(itype, pos, coord1, coord2, coord3, alt_hue_path)
plscmap1l.__doc__ = _plscmap1l.__doc__
