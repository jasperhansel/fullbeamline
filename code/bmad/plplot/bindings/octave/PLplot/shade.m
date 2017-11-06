## Copyright (C) 1998-2003  Joao Cardoso.
## Copyright (C) 2004  Rafael Laboissiere
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2 of the License, or (at your
## option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## This file is part of plplot_octave.

## shade(x, y, z [, levels [, contour] )
## shade(z [, levels [, contour] )
##
## plot shade of matrix z versus vectors x,y.
## level can be a scalar speficying the number of shade levels,
## or a vector, specifying where the shade level must be
## if contour exists, each shade level will be contoured with 'countour' color

function shade(x, y, z, levels, cont )

  global __pl
  strm = __pl_init;

  unwind_protect

  old_empty_list_elements_ok = warning("query", "Octave:empty-list-elements");
  warning("off","Octave:empty-list-elements");

  if (nargin == 1 && ismatrix(x))
    levels = 20;
    cont = 0;
    z = x; y = 1:rows(z); x = 1:columns(z);
  elseif (nargin == 2 && ismatrix(x))
    cont = 0; levels = y;
    z = x; y = 1:rows(z); x = 1:columns(z);
  elseif (nargin == 3 && isscalar(z))
    levels = y; cont = z;
    z = x; y = 1:rows(z); x = 1:columns(z);
  elseif (nargin == 3)
    levels = 20;
    cont = 0;
  elseif (nargin == 4)
    cont = 0;
  endif

  if (isscalar(levels) && levels == 1)	# segmentation violation!
    levels = 2;
  endif

  if (rows(x) > 1 && columns(x) > 1 && rows(y) > 1 && columns(y) > 1)
    xymat = 1;
  else
    xymat = 0;
  endif

  ## plot color and pen width of boundary of shade region

  max_color = 0; max_width = 0;

  if (cont)
    cont_color = cont; cont_width = 1;
  else
    cont_color = 0; cont_width = 0;
  endif

  if (ishold == 0)
    xmm = xm = min(min(x)); xM = max(max(x));
    ymm = ym = min(min(y)); yM = max(max(y));
    zm = min(min(z)); zM = max(max(z));

    if (__pl.axis_st(strm))
      xm = __pl.axis(strm,1); xM = __pl.axis(strm,2);
      if (xymat == 0)
        ix = find(x >= xm & x <= xM);
        x=x(ix); z=z(:,ix);
        xmm = min(x);
      endif

      if (length(__pl.axis(strm,:)) >= 4)	
	ym = __pl.axis(strm,3); yM = __pl.axis(strm,4);
        if (xymat == 0)
	  iy = find(y >= ym & y <= yM);
	  y=y(iy); z=z(iy,:);
	  ymm = min(y);
        endif
      else
	__pl.axis(strm,3) = ym; __pl.axis(strm,4) = yM;
      endif
      if (length(__pl.axis(strm,:)) == 6)
	zm = __pl.axis(strm,5); zM = __pl.axis(strm,6);
      else
	__pl.axis(strm,5) = zm; __pl.axis(strm,6) = zM;
      endif
    else	# make axis() return current axis
      __pl.axis(strm,1) = xm; __pl.axis(strm,2) = xM;
      __pl.axis(strm,3) = ym; __pl.axis(strm,4) = yM;
      __pl.axis(strm,5) = zm; __pl.axis(strm,6) = zM;		
    endif

    __pl.plcol(strm) = 1;
    plcol0(15);pllsty(1);
    __pl_plenv(xm, xM, ym, yM, 0, -2);
  else
    if (columns(__pl.axis(strm,:)) != 6)
      error("You must contour/shade plot something before entering hold mode");
    endif
    xmm = xm = __pl.axis(strm,1); xM = __pl.axis(strm,2);
    ymm = ym = __pl.axis(strm,3); yM = __pl.axis(strm,4);
    zm = __pl.axis(strm,5); zM = __pl.axis(strm,6);
    if (xymat == 0)
      ix = find(x >= xm & x <= xM);
      iy = find(y >= ym & y <= yM);
      z = z( iy, ix );
      x = x( ix );
      y = y( iy );
    endif
  endif

  maxx = max(max(x)); maxy = max(max(y)); minx = min(min(x)); miny = min(min(y));
  if (columns(x)>1 && rows(x) == 1)
    x = x';
  endif
  if (columns(y)>1 && rows(y) == 1)
    y = y';
  endif

  if (!isscalar(levels))
    n = length(levels)-1;
    clevel = levels;
    cclevel = clevel;
  else
    n = levels;
    clevel = linspace(zm, zM, levels+1);
    cclevel = linspace(zm, zM, levels);
  endif

  __pl.type(strm) = -2;
  __pl.lab_str = "";
  __pl.plcol(strm) = 1;
  __pl.pllsty(strm) = 1;	
  __pl.lab_pos(strm) = 1;

  plpsty(0);
  if (0) ## plshades() is slower than several calls to plshade() !? and plshades() sometimes fails ?!
    for i = 1:n
      plshade1(z', 0, minx, maxx, miny, maxy,
	      clevel(i), clevel(i+1),
	      1, (i-1) / (n-1), 1,
	      cont_color, cont_width, max_color, max_width, 1, x, y);
    endfor
  else
    if (columns(x) == 1 && columns(y) == 1)
      plshades1(z', minx, maxx, miny, maxy,
	     clevel', 1, cont_color, cont_width, 1, x, y);
    else
      plshades2(z', minx, maxx, miny, maxy,
	     clevel', 1, cont_color, cont_width, 1, x', y');
    endif
  endif

  for i = 1:n
    __pl.lab_str = [__pl.lab_str; sprintf("%#+.2G", cclevel(i))];
    __pl.lab_col(strm,__pl.lab_pos(strm)) = __pl.plcol(strm);
    __pl.lab_lsty(strm,__pl.lab_pos(strm)) = __pl.pllsty(strm);
    __pl.lab_pos(strm) = __pl.lab_pos(strm) + 1;				
    __pl.plcol(strm) = rem(__pl.plcol(strm), 15)+1;
    if  (__pl.line_style(strm))
      __pl.pllsty(strm) = rem(__pl.pllsty(strm), 8)+1;
    endif
  endfor

  if (__pl.grid(strm))		# this has to be done after shading
    plcol0(15);
    plbox("bcnsgt",0,0,"bcnsgtv",0,0)
  else
    plcol0(15);
    plbox("bcnst",0,0,"bcnstv",0,0)
  endif

  if (__pl.legend(strm))
    __pl_draw_legend;
  endif

  plcol0(15);
  pllab(tdeblank(__pl.xlabel(strm,:)), tdeblank(__pl.ylabel(strm,:)), tdeblank(__pl.tlabel(strm,:)));

  plflush;
  __pl.items(strm) = 1; # for now!

  unwind_protect_cleanup

  warning(old_empty_list_elements_ok.state, "Octave:empty-list-elements");

  end_unwind_protect


endfunction
