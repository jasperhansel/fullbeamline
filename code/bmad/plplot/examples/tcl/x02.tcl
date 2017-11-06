proc x02 {{w loopback}} {

    x02_demo1 $w
    x02_demo2 $w

# Restore defaults
    $w cmd plfont 1
    $w cmd plssub 1 1
    $w cmd pleop
    # $w cmd plcol0 1
}

# Demonstrates multiple windows and default color map 0 palette.
# Divides screen into 16 regions.

proc x02_demo1 {w} {

    $w cmd plbop

    $w cmd plssub 4 4

    x02_draw_windows $w 16

    $w cmd pleop
}

# Demonstrates multiple windows, user-modified color map 0 palette, and
# HLS -> RGB translation.  Divides screen into 100 regions.

proc x02_demo2 {w} {

    $w cmd plbop

    $w cmd plssub 10 10

# Set up cmap0
# Use 100 custom colors in addition to base 16
    set npts 100
    set offset 16
    set ntot [expr $npts + $offset]
    matrix r f $ntot
    matrix g f $ntot
    matrix b f $ntot

# Min & max lightness values
    set lmin 0.15
    set lmax 0.85

# The faster way to allocate also crashes, so stick with the slow way until we
# can resolve this.  Set the following to 1 to see the bug.. :)
    set see_the_bug 0
    if !$see_the_bug {
        $w cmd plscmap0n $ntot
    }

    for {set i 0} {$i < $npts} {incr i} {
        set i1 [expr $i + $offset]

    # Bounds on HLS, from plhlsrgb() commentary --
    #	hue		[0., 360.]	degrees
    #	lightness	[0., 1.]	magnitude
    #	saturation	[0., 1.]	magnitude

    # Vary hue uniformly from left to right
        set h [expr (360. / 10. ) * ( $i % 10 )]
    # Vary lightness uniformly from top to bottom, between min & max
        set l [expr $lmin + ($lmax - $lmin) * ($i / 10) / 9.]
    # Use max saturation
        set s 1.0

        $w cmd plhlsrgb $h $l $s r1 g1 b1
	# puts [format "%3d %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f" $i1 $h $l $s $r1 $g1 $b1]
        if $see_the_bug {
            r $i1 = [expr int($r1 * 255.001)]
            g $i1 = [expr int($g1 * 255.001)]
            b $i1 = [expr int($b1 * 255.001)]
        } else {
            set r2 [expr int($r1 * 255.001)]
            set g2 [expr int($g1 * 255.001)]
            set b2 [expr int($b1 * 255.001)]
  	    # puts [format "%3d %3d %3d %3d" $i1 $r2 $g2 $b2]
            $w cmd plscol0 $i1 $r2 $g2 $b2
        }
    }

    if $see_the_bug {
       # Load default cmap0 colors into our custom set
       for {set i 0} {$i < $offset} {incr i} {
	  plgcol0 $i r1 g1 b1
	  r $i = [expr int($r1)]
	  g $i = [expr int($g1)]
	  b $i = [expr int($b1)]
       }
       # temporary
       for {set i 0} {$i < $ntot} {incr i} {
	  set r1 [expr [r $i]]
	  set g1 [expr [g $i]]
	  set b1 [expr [b $i]]
  	  puts [format "%3d %3d %3d %3d" $i $r1 $g1 $b1]
       }
       # The following call currently segfaults.
       $w cmd plscmap0 r g b
    }

    x02_draw_windows $w 100 $offset

    $w cmd pleop
}


# Draws a set of numbered boxes with colors according to cmap0 entry.

proc x02_draw_windows { w nw {cmap0_offset 0} } {

    $w cmd plschr 0.0 3.5
    $w cmd plfont 4

    for {set i 0} {$i < $nw} {incr i} {
	$w cmd plcol0 [expr $i + $cmap0_offset]
	$w cmd pladv 0
	set vmin 0.1
	set vmax 0.9
	for {set j 0} {$j <= 2} {incr j} {
	    $w cmd plwidth [expr $j+1]
	    $w cmd plvpor $vmin $vmax $vmin $vmax
	    $w cmd plwind 0.0 1.0 0.0 1.0
	    $w cmd plbox "bc" 0.0 0 "bc" 0.0 0
	    set vmin [expr $vmin+0.1]
	    set vmax [expr $vmax-0.1]
	}
	$w cmd plwidth 1
	$w cmd plptex 0.5 0.5 1.0 0.0 0.5 $i
    }
}
