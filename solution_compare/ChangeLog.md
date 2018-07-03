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
