(*
   Illustrates backdrop plotting of world, US maps.
   Contributed by Wesley Ebisuzaki.

   Updated for OCaml by Hezekiah Carty

Copyright 2007, 2008  Hezekiah M. Carty

This file is part of PLplot.

PLplot is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

PLplot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with PLplot.  If not, see <http://www.gnu.org/licenses/>.
*)

open Plplot
open Printf

(*--------------------------------------------------------------------------*\
 * mapform19
 *
 * Defines specific coordinate transformation for example 19.
 * Not to be confused with mapform in src/plmap.c.
 * x[], y[] are the coordinates to be plotted.
\*--------------------------------------------------------------------------*)

let pi = atan 1.0 *. 4.0

(* The mapform callback gets a single pair of x and y coordinates and their
   index.
   It should return the properly transformed coordinates a an (x,y) tuple. *)
let mapform19 x y =
  let radius = 90.0 -. y in
  let xp = radius *. cos (x *. pi /. 180.0) in
  let yp = radius *. sin (x *. pi /. 180.0) in
  (xp, yp)

(* "Normalize" longitude values so that they always fall between -180.0 and
   180.0 *)
let normalize_longitude lon =
  if lon >= -180.0 && lon <= 180.0 then
    lon
  else
    let times = floor ((abs_float lon +. 180.0) /. 360.0) in
    if lon < 0.0 then
      lon +. 360.0 *. times
    else
      lon -. 360.0 *. times

(* A custom axis labeling function for longitudes and latitudes. *)
let geolocation_labeler axis value =
  let label_val, direction_label =
    match axis with
    | PL_Y_AXIS ->
        let label_val = value in
        label_val,
          if label_val > 0.0 then
            " N"
          else if label_val < 0.0 then
            " S"
          else
            "Eq"
    | PL_X_AXIS ->
        let label_val = normalize_longitude value in
        label_val,
          if label_val > 0.0 then
            " E"
          else if label_val < 0.0 then
            " W"
          else
            ""
    | PL_Z_AXIS -> invalid_arg "Invalid axis - only X or Y are supported"
  in
  if axis = PL_Y_AXIS && label_val = 0.0 then
    (* A special case for the equator *)
    sprintf "%s" direction_label
  else
    sprintf "%.0f%s" (abs_float label_val) direction_label

(*--------------------------------------------------------------------------*\
 * main
 *
 * Shows two views of the world map.
\*--------------------------------------------------------------------------*)

let () =
  (* Longitude (x) and latitude (y) *)
  let miny = -70.0 in
  let maxy = 80.0 in

  (* Parse and process command line arguments *)
  plparseopts Sys.argv [PL_PARSE_FULL];

  plinit();

  (* Cartesian plots *)
  (* Most of world *)

  let minx = -170.0 in
  let maxx = minx +. 360.0 in

  (* Setup a custom latitude and longitude-based scaling function. *)
  plslabelfunc geolocation_labeler;

  plcol0 1;
  plenv minx maxx miny maxy 1 70;
  (* No transform function is passed to plmap.  Since we have not set one yet,
     it defaults to using an identity function (xp = x, yp = y) *)
  plmap "usaglobe" minx maxx miny maxy;

  (* The Americas *)
  let minx = 190.0 in
  let maxx = 340.0 in

  plcol0 1;
  plenv minx maxx miny maxy 1 70;
  (* Again, we have not set a transform.  Everything remains in a Cartesian
     projection. *)
  plmap "usaglobe" minx maxx miny maxy;

  (* Clear the labeling function *)
  plunset_labelfunc ();

  (* Polar, Northern hemisphere *)
  let minx = 0.0 in
  let maxx = 360.0 in

  plenv (-75.) 75. (-75.) 75. 1 (-1);
  (* Now we set mapform19 as our callback.  It will be used for every following
     plmap and plmeridians call until set to something else. *)
  plset_mapform mapform19;
  plmap "globe" minx maxx miny maxy;

  pllsty 2;
  (* This call to plmeridians is also after the set_mapform call, so it uses
     the same projection as the plmap call above. *)
  plmeridians 10.0 10.0 0.0 360.0 (-10.0) 80.0;

  plunset_mapform ();

  (* Polar, Northern hemisphere, this time with a PLplot-wide transform *)
  let minx = 0.0 in
  let maxx = 360.0 in

  plstransform mapform19;

  pllsty 1;
  plenv (-75.) 75. (-75.) 75. 1 (-1);
  (* No need to set the map transform here as the global transform will be
     used. *)
  plmap "globe" minx maxx miny maxy;

  pllsty 2;
  plmeridians 10.0 10.0 0.0 360.0 (-10.0) 80.0;

  (* Show Baltimore, MD on the map *)
  plcol0 2;
  plssym 0.0 2.0;
  plpoin [|-76.6125|] [|39.2902778|] 18;
  plssym 0.0 1.0;
  plptex ~-.76.6125 43.0 0.0 0.0 0.0 "Baltimore, MD";

  (* For OCaml, this is how the global transform is cleared *)
  plunset_transform ();

  plend ();
  ()
