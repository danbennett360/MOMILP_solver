## 
The code in this directory is used to compare the solutions produced
by  polyscip and MOP.   

The input file file format is 
   [ d1, d2, ..., dn] var=value var=value ...

command:
   **comp** file1 file2 _options_

   Note: the filenames must be given as the first parameters.

   **options** 

   * -noDistance : do not show the distance between the two surfaces.  Note this will not be calculated if the two sets of points are the same
   * -noFull : do not show the differences between non-whatever variables.   
   * -noDiff : do not show the difference between points in file 1 and file 2

---
## 7/2/18 
* Continued to change compare.c
** Added checking of whatever you call the other variables generated
** reorganized the code to allow calling of only some tests.
** added command line flags

## 7/25/18
* Set difference was going backwards, I think.
* Add extra rows to the problems being solved so that a negative multiple of the extreme directions can be used to create "affine" combination of new point. If a negative multiple is used, it shows the point sits outside the polyhedron.

## 7/25/18 - Entry 2
* Modified use of the negative multiple of the extreme directions to only be useable by passing the flag "-Ext". It is not helpful in most cases -- it should be used if an error "Unable to solve/get ..." is encountered.

## 7/26/18 
* Completely reworked the optimization routine that generates the distance. I took advantage of CPLEX's ability to optimize over a quadratic objective, so the output of the problem is no longer the minimum multiplier in the convex combination, but a true distance. The output of the problem is now the minimum distance from a single point in the second set of points to the polytope formed from the first set of points. As a result, by taking the maximum of these minimum distances we produce the Hausdorff distance (https://en.wikipedia.org/wiki/Hausdorff_distance) from the polytope constructed from the first set of points to the second set of points (not its polytope, though). THIS is a valid metric and therefore acceptable to report. 
* These changes nullified the "-Ext" flag. It is gone.
* Added a "-RelDist" flag that, when called, computes the magnitude of the postition vector associated with each point in set 2 (i.e., the distance of the point from the origin), takes the minimum of these magnitudes, and reports the Hausdorff distance as a measure relative to this magnitude, i.e., RelDist = HausdorffDist/min(magnitudes).
