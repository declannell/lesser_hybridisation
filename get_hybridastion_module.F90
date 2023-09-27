MODULE MyModule
    contains

    SUBROUTINE get_lesser_hybridisation(gf_retarded, gf_lesser, se_lesser, lesser_hybridisation, num_orbitals, steps)
        implicit none
        integer :: r, spin 
        integer, intent(in) :: steps, num_orbitals
        complex*16, intent(in) :: gf_lesser(:,:,:,:), se_lesser(:,:,:,:)
        complex*16, intent(inout) :: lesser_hybridisation(:,:,:,:), gf_retarded(:,:,:,:)
        complex*16, allocatable :: gf_adv(:,:)
        allocate(gf_adv(num_orbitals, num_orbitals))
        
        do r = 1, steps
            do spin = 1, 2
                !this is calculating the adjoint of gf_retarded
                call ComputeAdjoint(gf_retarded(r,spin, :,:), gf_adv, num_orbitals)
                
                !this overwrites gf_adv and gf_retarded as their inverse
                call invert_matrix(num_orbitals, gf_adv)
                call invert_matrix(num_orbitals, gf_retarded(r, spin, :,:))
                
                !calculation of the lesser hybridiation
                lesser_hybridisation(r, spin, :, :) = matmul(gf_retarded(r, spin, :, :),&
                     matmul(gf_lesser(r, spin, :, :) , gf_adv)) - se_lesser(r, spin, :, :)
            end do
        end do

    END SUBROUTINE get_lesser_hybridisation


    SUBROUTINE ReadDataFromFile(filename, A, num_orbitals, energy)
        implicit none
        character(len=*), intent(in) :: filename
        character(50) :: filename_complete
        integer, intent(in) :: num_orbitals
        complex*16, intent(out) :: A(:,:,:,:)
        real (8), intent(out) :: energy(:)
        real(8) real_part, imag_part
        integer :: i, j, r = 1, ios
        logical :: end_of_file
      
        do i = 1, num_orbitals
            do j = 1, num_orbitals
              r = 1
              !Create the file name as a string
              write(filename_complete, '(A,I0,A,I0,A)') filename, i-1, '_', j-1, '.dat'
              print *, "Reading: ", filename_complete
                open(10, FILE=filename_complete, STATUS='OLD', ACTION='READ')
                do while (.TRUE.)
                    read(10, *, IOSTAT=ios) energy(r), real_part, imag_part
                    ! Check for end of file
                    if (ios /= 0) then
                        end_of_file = .TRUE.
                        exit
                    end if
                    A(r, 1, i, j) = CMPLX(real_part, imag_part, 8)
                    A(r, 2, i, j) = CMPLX(real_part, imag_part, 8) 
                    r = r + 1
                end do
            
              ! Close the file
              close(10)
            end do
        end do     
      END SUBROUTINE ReadDataFromFile  

      SUBROUTINE invert_matrix(M, A)
        implicit none
        integer, intent(in) :: M
        complex*16, intent(inout), dimension(M, M) :: A
        complex*16, allocatable, dimension(:) :: WORK
        integer, allocatable, dimension(:) :: IPIV
        integer info, error
    
        allocate(WORK(M), IPIV(M), stat=error)
        if (error /= 0) then
          print *, "error: not enough memory"
          return
        end if
    
        call ZGETRF(M, M, A, M, IPIV, info)
        if (info /= 0) then
          write(*,*) "LU factorization failed"
        end if
    
        call ZGETRI(M, A, M, IPIV, WORK, M, info)
        if (info /= 0) then
          write(*,*) "Matrix inversion failed"
        end if
    
        deallocate(IPIV, WORK, stat=error)
    END SUBROUTINE invert_matrix

    

    SUBROUTINE ComputeAdjoint(A, AdjA, N)
        implicit none
        integer, intent(in) :: N
        complex*16, intent(in) :: A(N, N)
        complex*16, intent(out) :: AdjA(N, N)
        integer :: i, j
    
        ! Compute the adjoint of the complex matrix A
        do i = 1, N
            do j = 1, N
                AdjA(j, i) = CONJG(A(i, j))
            end do
        end do
    END SUBROUTINE ComputeAdjoint


    SUBROUTINE WriteToFile(filename, object, num_orbitals, steps, energy)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(in) :: num_orbitals, steps
        complex*16, intent(in) :: object(:,:,:,:)
        real (8), intent(in) ::energy(:)
        integer :: i, j, r
        character(80) :: filename_complete
    
        do i = 1, num_orbitals
            do j = 1, num_orbitals
                ! Generate the complete filename including path
                write(filename_complete, '(A,"_",I0,"_",I0,"_fortran.dat")') filename, i-1, j-1

                print *, "Writing to file: ", filename_complete
    
                ! Open the file for writing (creates a new file)
                open(10, FILE=filename_complete, STATUS='REPLACE', ACTION='WRITE')
    
                ! Write data to the file
                do r = 1, steps
                    write(10, *) energy(r), REAL(object(r, 1, i, j)), AIMAG(object(r, 2, i, j))
                end do
    
                ! Close the file
                CLOSE(10)
            end do
        end do
    END SUBROUTINE
  END MODULE MyModule
