## July 2, 2018
 * Added Epsilon() to several of the classes in an attempt to allow the user to chagne epsilon at run time.
 ** I'm not sure if this is completely correct, but it seems to work.
 
 ## July 3, 2018
 * Added ability to read *.mop files.
 ** I cheated here. I generate the first problem by reading the original *.mop file, then for each additional problem I just copied the file used to create the previous problem to a new file with the first occurence of an "N" row moved to the bottom of the "ROWS" section and read in the modified file. It works on the one example I tried -- I verified that the solution is exactly the same using the *.mop file and the 3 *.lp files.
