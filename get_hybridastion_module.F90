MODULE MyModule
    CONTAINS

    SUBROUTINE initialize_mat(a, b, steps, num_orbitals)
        IMPLICIT NONE
        COMPLEX(8), INTENT(OUT) :: a(:,:,:,:)
        COMPLEX(8), INTENT(OUT) :: b(:,:,:,:)
        INTEGER, INTENT(IN) :: steps, num_orbitals
        INTEGER :: r, i, j, spin
      
        do r = 1, steps
            do spin = 1, 2
                do i = 1, num_orbitals 
                    do j = 1, num_orbitals 
                      a(r, spin, i, j) = CMPLX(REAL(i) + REAL(j), REAL(i) - REAL(j), 8) ! Create a complex number
                      if (i == j) THEN
                        b(r, spin, i, j) = CMPLX(r, 0, 8)
                      else
                        b(r, spin, i, j) = CMPLX(0, r * r)
                      end if
                      write(*, *) a(r, spin, i, j), b(r, spin, i, j), i, j 
                    end do
                end do
                print *, "the spin is " , spin, "the steps is ", steps
                WRITE(*, '(A)', ADVANCE='YES') 
            end do    
        end do
    END SUBROUTINE initialize_mat

    SUBROUTINE multiple_mat(c,a, b, steps)
        IMPLICIT NONE
        COMPLEX(8), INTENT(OUT) :: c(:,:,:,:)
        COMPLEX(8), INTENT(IN) :: a(:,:,:,:)
        COMPLEX(8), INTENT(IN) :: b(:,:,:,:)
        INTEGER, INTENT(IN) :: steps
        INTEGER :: r, i, j, spin
      
        do r = 1, steps
            do spin = 1, 2
                c(r, spin, :,:) = MATMUL(a(r, spin, :,:), b(r, spin, :,:))
            end do
        end do
    END SUBROUTINE multiple_mat

    SUBROUTINE PrintMatrix(c, steps, num_orbitals)
        IMPLICIT NONE
        COMPLEX(8), INTENT(IN) :: c(:,:,:,:)
        INTEGER, INTENT(IN) :: steps, num_orbitals
        INTEGER :: r, i, j, spin
      
        do r = 1, steps
            do spin = 1, 2
                do i = 1, num_orbitals
                  do j = 1, num_orbitals
                    write(*, *) c(r, spin, i, j), i, j
                  end do
                end do
                WRITE(*, '(A)', ADVANCE='YES')     
            end do
        end do
    END SUBROUTINE PrintMatrix

    SUBROUTINE ReadDataFromFile(filename, A, num_orbitals)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(50) :: filename_complete
        INTEGER, INTENT(IN) :: num_orbitals
        complex(8), INTENT(OUT) :: A(:,:,:,:)
        real(8) :: energy, real_part, imag_part
        INTEGER :: i, j, r = 1, ios
        LOGICAL :: end_of_file
      
        DO i = 1, num_orbitals
            DO j = 1, num_orbitals
              r = 1
              !Create the file name as a string
              WRITE(filename_complete, '(A,I0,A,I0,A)') filename, i-1, '_', j-1, '.dat'
              print *, filename_complete
                OPEN(10, FILE=filename_complete, STATUS='OLD', ACTION='READ')
                DO WHILE (.TRUE.)
                    READ(10, *, IOSTAT=ios) energy, real_part, imag_part
                    ! Check for end of file
                    IF (ios /= 0) THEN
                        end_of_file = .TRUE.
                        EXIT
                    END IF
                    A(r, 1, i, j) = CMPLX(real_part, imag_part, 8)
                    A(r, 2, i, j) = CMPLX(real_part, imag_part, 8)
                    !PRINT *, "Energy:", energy, "Real:", real_part, "Imag:", imag_part, "r: ", r,&
                    !   "A(r, 1, i, j): ", A(r, 1, i, j)
                    r = r + 1
                END DO
            
              ! Close the file
              CLOSE(10)
            END DO
        END Do     
      END SUBROUTINE ReadDataFromFile  
  END MODULE MyModule