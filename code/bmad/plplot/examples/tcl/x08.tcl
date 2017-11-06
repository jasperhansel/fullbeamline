#----------------------------------------------------------------------------
# Copyright (C) 2004  Alan W. Irwin
#
# This file is part of PLplot.
#
# PLplot is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library General Public License as published
# by the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# PLplot is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with PLplot; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#----------------------------------------------------------------------------

# Does a series of 3-d plots for a given data set, with different
# viewing options in each plot.

# Routine for restoring colour map1 to default.
# See static void plcmap1_def(void) in plctrl.c for reference.
proc restore_cmap1_8 {w} {
   # For center control points, pick black or white, whichever is closer to bg
   # Be careful to pick just short of top or bottom else hue info is lost
   $w cmd plgcolbg rbg gbg bbg
   set vertex [expr {($rbg + $gbg + $bbg)/(3.*255.)}]
   if {$vertex < 0.5} {
      set vertex 0.01
      set midpt 0.10
   } else {
      set vertex 0.99
      set midpt 0.90
   }
   # Independent variable of control points.
   matrix i f 6 = {0., 0.44, 0.50, 0.50, 0.56, 1.}
   # Hue for control points.  Blue-violet to red
   matrix h f 6 = {260., 260., 260., 0., 0., 0.}
   # Lightness ranging from medium to vertex to medium
   # double quotes required rather than braces so that variables
   # get evaluated for initialization.
   matrix l f 6  = "0.5, $midpt, $vertex, $vertex, $midpt, 0.5"
   # Saturation is complete for default
   matrix s f 6 = {1., 1., 1., 1., 1., 1.}
   # Integer flag array is zero (no interpolation along far-side of colour
   # figure)
   matrix rev i 6 = {0, 0, 0, 0, 0, 0}
   # Default number of cmap1 colours
   $w cmd plscmap1n 128
   # Interpolate between control points to set up default cmap1.
   $w cmd plscmap1l 0 i h l s rev
}

# Routine for initializing color map 1 in HLS space.
# Basic grayscale variation from half-dark (which makes more interesting
# looking plot compared to dark) to light.
proc cmap1_init_8 {w gray} {
   # Independent variable of control points.
   matrix i f 2 = {0., 1.}
   if {$gray == 1} {
      # Hue for control points.  Doesn't matter since saturation is zero.
      matrix h f 2 = {0., 0.}
      # Lightness ranging from half-dark (for interest) to light.
      matrix l f 2 = {0.5, 1.}
      # Gray scale has zero saturation.
      matrix s f 2 = {0., 0.}
   } else {
      # Hue ranges from blue (240 deg) to red (0 or 360 deg)
      matrix h f 2 = {240., 0.}
      # Lightness and saturation are constant (values taken from C example).
      matrix l f 2 = {0.6, 0.6}
      matrix s f 2 = {0.8, 0.8}
   }

   # Integer flag array is zero (no interpolation along far-side of colour
   # figure
   matrix rev i 2 = {0, 0}
   # Number of cmap1 colours is 256 in this case.
   $w cmd plscmap1n 256
   # Interpolate between control points to set up default cmap1.
   $w cmd plscmap1l 0 i h l s rev
}

proc x08 {{w loopback}} {

    # these should be defined elsewhere.
    set rosen 0

    matrix alt f 2 = {60.0, 40.0}
    matrix az  f 2 = {30.0, -30.0}

    set xpts 35
    set ypts 45
    set n_col 256
    set two_pi [expr {2.0 * $::PLPLOT::PL_PI} ]

    matrix x f $xpts
    matrix y f $ypts
    matrix z f $xpts $ypts

    set dx [expr {2.0/($xpts-1)}]
    set dy [expr {2.0/($ypts-1)}]
    for {set i 0} {$i < $xpts} {incr i} {
	x $i = [expr {-1.0 + $i * $dx}]
	if {$rosen == 1} {
	   x $i = [expr {1.5* [x $i]}]
	}
    }

    for {set i 0} {$i < $ypts} {incr i} {
	y $i = [expr {-1.0 + $i *$dy}]
	if {$rosen == 1} {
	   y $i = [expr {0.5 + [y $i]}]
	}
    }

    for {set i 0} {$i < $xpts} {incr i} {
	set xx [x $i]
	for {set j 0} {$j < $ypts} {incr j} {
	    set yy [y $j]
	    if {$rosen == 1} {
	      z $i $j = [expr {(pow(1. - $xx,2) + \
		100 * pow($yy - pow($xx,2),2))}]
	      set zz [z $i $j]
	      if {$zz > 0.} {
		z $i $j = [expr {log($zz)}]
	      } else {
		z $i $j = [expr {-5.0}]
	      }
	    } else {
	    set r [expr {sqrt( $xx * $xx + $yy * $yy )} ]

	    z $i $j = [expr {exp(-$r * $r) * cos( $two_pi * $r )} ]
	    }
	}
    }

    # parameters of ellipse (in x, y index coordinates) that limits the data.
    # x0, y0 correspond to the exact floating point centre of the index
    # range.
    set x0 [expr {0.5 * ( $xpts - 1 )}]
    set a  [expr {0.9 * $x0}]
    set y0 [expr {0.5 * ( $ypts - 1 )}]
    set b  [expr {0.7 * $y0}]

    matrix zlimited  f $xpts $ypts
    matrix indexymin i $xpts
    matrix indexymax i $xpts

    set indexxmin 0
    set indexxmax $xpts
    for {set i $indexxmin} {$i < $indexxmax} {incr i} {
        indexymin $i = 0
        indexymax $i = $ypts
        for {set j $indexymin} {$j < $indexymax} {incr j} {
            zlimited $i $j = [expr {1.0e300}]
        }
    }

    for {set i $indexxmin} {$i < $indexxmax} {incr i} {
        set square_root [expr {sqrt(1. - min( 1., pow( ( $i - $x0 ) / $a, 2. ) ) )}]
        # Add 0.5 to find nearest integer and therefore preserve symmetry
        # with regard to lower and upper bound of y range.
        indexymin $i = [expr {max( 0, int( 0.5 + $y0 - $b * $square_root ) )}]
        # indexymax calculated with the convention that it is 1
        # greater than highest valid index.
        indexymax $i = [expr {min( $ypts, 1 + int( 0.5 + $y0 + $b * $square_root ) )}]

        for { set j [indexymin $i]} {$j < [indexymax $i]} {incr j} {
            zlimited $i $j = [z $i $j]
        }
    }

    set zmin [z min [ expr $xpts * $ypts]]
    set zmax [z max [ expr $xpts * $ypts]]

    set nlev 10
    matrix clev f $nlev
    set step [expr {($zmax-$zmin)/($nlev+1)}]
    for {set i 0} {$i < $nlev} {incr i} {
	clev $i = [expr {$zmin + ($i+1) * $step}]
    }

    $w cmd pllightsource 1. 1. 1.
    for {set k 0} {$k < 2} {incr k} {
	for {set ifshade 0} {$ifshade < 5} {incr ifshade} {
	    $w cmd pladv 0
	    $w cmd plvpor 0.0 1.0 0.0 0.9
	    $w cmd plwind -1.0 1.0 -0.9 1.1
	    $w cmd plcol0 3
	    set title [format "#frPLplot Example 8 - Alt=%.0f, Az=%.0f" \
			   [alt $k] [az $k]]
	    $w cmd plmtex "t" 1.0 0.5 0.5 $title
	    $w cmd plcol0 1
	    if {$rosen == 1} {
	       $w cmd plw3d 1.0 1.0 1.0 -1.5 1.5 -0.5 1.5 $zmin $zmax [alt $k] [az $k]
	    } else {
	       $w cmd plw3d 1.0 1.0 1.0 -1.0 1.0 -1.0 1.0 $zmin $zmax [alt $k] [az $k]
	    }
	    $w cmd plbox3 "bnstu" "x axis" 0.0 0 \
	                  "bnstu" "y axis" 0.0 0 \
	    		  "bcdmnstuv" "z axis" 0.0 0
  	    $w cmd plcol0 2
	    # diffuse light surface plot
	    if {$ifshade == 0} {
	       cmap1_init_8 $w 1
	       $w cmd plsurf3d x y z 0
	    # magnitude colored plot
	    } elseif {$ifshade == 1} {
	       cmap1_init_8 $w 0
	
	       # Note: the [expr] command is essential here - plsurf3d doesn't accept
	       # the hexadecimal form apparently!
	       $w cmd plsurf3d x y z [expr {$::PLPLOT::MAG_COLOR}]
	    # magnitude colored plot with faceted squares
	    } elseif {$ifshade == 2} {
	       cmap1_init_8 $w 0
	       $w cmd plsurf3d x y z [expr {$::PLPLOT::MAG_COLOR | $::PLPLOT::FACETED}]
	    # magnitude colored plot with contours.
	    } elseif {$ifshade == 3}  {
	       cmap1_init_8 $w 0
	       $w cmd plsurf3d x y z \
		 [expr {$::PLPLOT::MAG_COLOR | $::PLPLOT::SURF_CONT | $::PLPLOT::BASE_CONT}] clev
	    } else {
	       cmap1_init_8 $w 0
	       $w cmd plsurf3dl x y z \
		 [expr {$::PLPLOT::MAG_COLOR | $::PLPLOT::SURF_CONT | $::PLPLOT::BASE_CONT}] clev \
		 $indexxmin indexymin indexymax
	    }
        }
    }

# Restore defaults
    # $w cmd plcol0 1
    restore_cmap1_8 $w
}
