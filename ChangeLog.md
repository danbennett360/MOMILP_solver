## July 2, 2018
 * Added Epsilon() to several of the classes in an attempt to allow the user to chagne epsilon at run time.
 ** I'm not sure if this is completely correct, but it seems to work.
 
 ## July 3, 2018
 * Added ability to read *.mop files.
 ** I cheated here. I generate the first problem by reading the original *.mop file, then for each additional problem I just copied the file used to create the previous problem to a new file with the first occurence of an "N" row moved to the bottom of the "ROWS" section and read in the modified file. It works on the one example I tried -- I verified that the solution is exactly the same using the *.mop file and the 3 *.lp files.
 
 ## July 15, 2018
 * Note that some of these changes were from a previous iteration, but I made so many that they need documented.
 * Added several command line flags: (i) -showprogress (T or F) which allows for printing a summary of the results so far periodically -- very useful for debugging, helps us to know if things are going wrong (for example, if the number of saved solutions never increases, or increases by the same amount every time a summary is printed), (ii) -progressvalue (positive integer) which specifies that the summary should be printed every '-progressvalue' iterations, (iii) -debug (0, 1, 2) which allows the user to turn debugging on/off without recompiling -- setting 1 turns on a check for repeated simplices and a check for negative normals, setting 2 is same as setting 1 AND turns on the standard out printing associated with the DEBUG variable, (iv) -reldist (T or F) which allows for computing distance between points and simplices relatively -- this should cause less need for 'tuning' epsilon between problems and between dimensions.
 * Decreased size of epsilon when checking for 'shadowing' of adjacent simplices.
 * Moved check for positive normal vector components to after the vector has been normalized and added that these checks have epsilon tolerance. This has two important implications: (i) -normalize is no longer a valid command line flag because normalization is necessary for correct performance, (ii) not using epsilon tolerance when checking the sign of normal components was causing a vector such as (-2e-16, 34, 5, 15) to be seen as problematic, i.e., having a negative component even though the 'negative' component is essentially zero (this was the cause of the 'There is a simplex on the stack that is not oriented correctly! Exiting!' error).
 * Added a global epsilon so that it is easy to set the epsilon value for the simplex class.
 * Added a check for dominated points in the computation of the extreme points that form the initial simplex. If any of these points are repeated or dominate one another it means that there is a pair of objectives that are not conflicting, i.e., they always produce the same solutions, and thus both objectives from the pair need not be considered together.
