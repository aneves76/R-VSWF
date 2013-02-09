factorial <- function(n) {
# If the library hasn't been loaded yet, load it
if (!is.loaded(symbol.For('factorial_R_wrapper'))) {
    dyn.load('factorial.so')
    }
# Call the subroutine, giving n as a Fortran integer, and giving
# result as an array with space for 1 integer.
# A list of parameters values after the function is called
# is returned, assigned to the same names as are given before the
# = signs in the arguments.
# NOTE:  factorial_R_wrapper has its R in lower-case now.  This
# is because Fortran is case insensitive, so usually names are
# all converted to lower-case.
returned_data = .Fortran('factorial_R_wrapper', 
                         n=as.integer(n), 
                         result=integer(1))
# Return the value of the result parameter 
return(returned_data$result)
}
