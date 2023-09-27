MODULE MyModule
    CONTAINS

    SUBROUTINE get_lesser_hybridisation(gf_retarded, gf_lesser, se_lesser, lesser_hybridisation, num_orbitals, steps)
        IMPLICIT NONE
        INTEGER :: i, j, r, spin 
        INTEGER, INTENT(IN) :: steps, num_orbitals
        COMPLEX*16, INTENT(IN) :: gf_retarded(:,:,:,:), gf_lesser(:,:,:,:), se_lesser(:,:,:,:)
        COMPLEX*16, INTENT(OUT) :: lesser_hybridisation(:,:,:,:)
        COMPLEX*16, ALLOCATABLE :: gf_adv(:,:)
        ALLOCATE(gf_adv(num_orbitals, num_orbitals))
        
        do r = 1, steps
            do spin = 1, 2
                ! Query the optimal workspace size for each iteration
                call ComputeAdjoint(gf_retarded(r,spin, :,:), gf_adv, num_orbitals)

                if (r == 1 .and. spin == 1) then
                    do i =1, num_orbitals
                        do j = 1, num_orbitals
                            write(*, *) gf_retarded(r,spin, i,j), gf_adv(i,j), i, j
                        end do
                    end do
                end if
            end do
        end do

    END SUBROUTINE get_lesser_hybridisation


    SUBROUTINE ReadDataFromFile(filename, A, num_orbitals)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(50) :: filename_complete
        INTEGER, INTENT(IN) :: num_orbitals
        COMPLEX*16, INTENT(OUT) :: A(:,:,:,:)
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
                    if (r == 1 .and. filename == "gf_retarded_") then
                        PRINT *, "Energy:", energy, "Real:", real_part, "Imag:", imag_part, "r: ", r,&
                            "A(r, 1, i, j): ", A(r, 1, i, j)
                    end if
 
                    r = r + 1
                END DO
            
              ! Close the file
              CLOSE(10)
            END DO
        END Do     
      END SUBROUTINE ReadDataFromFile  



    
    SUBROUTINE ComputeAdjoint(A, AdjA, N)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        COMPLEX*16, INTENT(IN) :: A(N, N)
        COMPLEX*16, INTENT(OUT) :: AdjA(N, N)
        INTEGER :: i, j
    
        ! Compute the adjoint of the complex matrix A
        DO i = 1, N
            DO j = 1, N
                AdjA(j, i) = CONJG(A(i, j))
            END DO
        END DO
    END SUBROUTINE ComputeAdjoint


  END MODULE MyModule


  MODULE MytestModule
    private

    public :: initialize_and_invert_matrix
    CONTAINS
    subroutine initialize_and_invert_matrix(M, A)
        implicit none
        integer, intent(in) :: M
        complex*16, intent(inout), dimension(M, M) :: A
        complex*16, allocatable, dimension(:) :: WORK
        integer, allocatable, dimension(:) :: IPIV
        integer i, j, info, error
    
        allocate(WORK(M), IPIV(M), stat=error)
        if (error /= 0) then
          print *, "error: not enough memory"
          return
        end if
    
        ! Definition of the test matrix A
        do i = 1, M
          do j = 1, M
            if (j == i) then
              A(i, j) = (1, 1)
            else
              A(i, j) = (2, 0)
            end if
          end do
        end do
    
        call ZGETRF(M, M, A, M, IPIV, info)
        if (info == 0) then
          write(*,*) "LU factorization succeeded"
        else
          write(*,*) "LU factorization failed"
        end if
    
        call ZGETRI(M, A, M, IPIV, WORK, M, info)
        if (info == 0) then
          write(*,*) "Matrix inversion succeeded"
        else
          write(*,*) "Matrix inversion failed"
        end if
    
        ! Print the inverted matrix
        do i = 1, M
          do j = 1, M
            write(*,*) A(i, j)
          end do
        end do
    
        deallocate(IPIV, WORK, stat=error)
        if (error /= 0) then
          print *, "error: failed to release memory"
        end if
      end subroutine initialize_and_invert_matrix
    
    subroutine initialize_and_invert(m, num_orbitals)
        implicit none
        real(8), intent(out) :: m(:,:)
        integer, intent(in) :: num_orbitals
        integer i, j
        integer :: info
        integer, allocatable :: ipiv(:)
        real(8), allocatable :: work(:)
        integer :: LWORK
        real(8) :: real_LWORK
    
        ! Allocate memory for workspace and pivot arrays
        allocate(ipiv(num_orbitals))
        
        ! Initialize the matrix
        do i = 1, num_orbitals
            do j = 1, num_orbitals
                m(i, j) = real(i + j, 8)
                write(*, *) m(i, j) 
            end do
        end do
    
        ! Query the optimal workspace size for LAPACK matrix inversion
        LWORK = -1
        call dgetrf(num_orbitals, num_orbitals, m, num_orbitals, ipiv, info)
        print *, "here"
        ! Check for errors during factorization
        if (info /= 0) then
            print *, "Error: LU factorization failed. INFO =", info
            return
        end if
        
        ! Query the workspace size
        real_LWORK = 0.0  ! Initialize to zero
        print *, "trying to invert"
        call dgetri(num_orbitals, m, num_orbitals, ipiv, work, LWORK, info)
        print *, "tried the first inversion"
        if (info /= 0) then
            ! If LWORK is too small, requery with the correct size
            LWORK = int(real_LWORK) + 1
            allocate(work(LWORK))
        else
            ! LWORK is already correct, allocate workspace
            allocate(work(LWORK))
        end if
        print *, "here"
        ! Perform matrix inversion using LAPACK
        call dgetri(num_orbitals, m, num_orbitals, ipiv, work, LWORK, info)
    
        ! Deallocate workspace and pivot array
        deallocate(work)
        deallocate(ipiv)
        
        ! Check for errors during matrix inversion
        if (info /= 0) then
            print *, "Error: Matrix inversion failed. INFO =", info
        end if
    
        do i = 1, num_orbitals
            do j = 1, num_orbitals
                write(*, *) m(i, j) 
            end do
        end do
    end subroutine initialize_and_invert
    
    



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


  END MODULE MytestModule