      INTEGER FUNCTION factorial(n)
      INTEGER n, i
      factorial = 1
      DO i = 2,n
          factorial = factorial * i
      END DO
      END FUNCTION
