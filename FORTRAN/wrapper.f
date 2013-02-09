       SUBROUTINE factorial_R_wrapper(n, answer)
       INTEGER n, answer, factorial
       EXTERNAL factorial
       answer = factorial(n)
       END SUBROUTINE
