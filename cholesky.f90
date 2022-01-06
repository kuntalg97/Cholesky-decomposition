!     Cholesky Decomposition
!     Kuntal Ghosh
!     September 2020

!---------------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND (15,307)
      INTEGER :: i,j,k,ios,n
      REAL (KIND=dp), ALLOCATABLE :: x(:,:),L(:,:),LT(:,:),mulmat(:,:)
      !INTEGER, ALLOCATABLE :: x(:,:) 
      REAL (KIND=dp) :: sum

      OPEN (1, file = 'matrix.dat', status = 'old')
    
      n = 0
      linecount: DO 
            READ (1,*,IOSTAT=ios)
            if (ios/=0) EXIT linecount
            n = n + 1
      END DO linecount
      REWIND(1)

      ALLOCATE (x(n,n),L(n,n),LT(n,n),mulmat(n,n))

      DO i = 1,n
            READ (1,*) (x(i,j),j=1,n)
      END DO

      WRITE (*,*) "Original matrix -->"
      DO i = 1,n          !Use when x(n,n) has real elements
            WRITE (*, fmt = '(5F16.6)') (x(i,j),j=1,n) 
      END DO

      !DO i = 1,n         !Use when x(n,n) has integer elements
      !      DO j = 1,n   !Run j = 1,i to print only lower triangular face
      !          WRITE (*, fmt = '(i3,$)') x(i,j)
      !      END DO
      !      WRITE (*,*)
      !END DO

      DO i = 1,n
            DO j = 1,n
                  L(i,j) = 0._dp
            END DO
      END DO      

      DO j = 1,n
            DO i = j,n
                  IF (j.eq.i) THEN
                          sum = 0._dp
                          DO k = 1,j-1
                               sum = sum + L(i,k)**2
                          END DO
                          L(i,j) = dsqrt(x(i,j) - sum)
                  ELSE
                          sum = 0._dp
                          DO k = 1,j-1
                               sum = sum + L(i,k)*L(j,k)
                          END DO
                          L(i,j) = (x(i,j) - sum)/L(j,j)
                  END IF
            END DO
      END DO      

      WRITE (*,*)
      WRITE (*,*) "Lower triangular matrix -->"
      DO i = 1,n  
            WRITE (*, fmt = '(5F16.6)') (L(i,j),j=1,n) 
      END DO

      LT = transpose(L)

      WRITE (*,*) "Tranpose of L matrix -->"
      DO i = 1,n  
            WRITE (*, fmt = '(5F16.6)') (LT(i,j),j=1,n) 
      END DO

      mulmat = matmul(L,LT)

      WRITE (*,*) "L*LT matrix -->"
      DO i = 1,n  
            WRITE (*, fmt = '(5F16.6)') (mulmat(i,j),j=1,n) 
      END DO

      DEALLOCATE (x,L,LT,mulmat)

      STOP
      END
