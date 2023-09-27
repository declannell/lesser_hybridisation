  PROGRAM ComplexMatrixExample
    use MyModule
    use MytestModule

        implicit none
        integer,parameter::M = 2
        complex*16,allocatable,dimension(:,:)::A
        allocate(A(M,M))
        call initialize_and_invert_matrix(M, A)
        deallocate(A)
        
    !IMPLICIT NONE
    !integer, parameter :: dp = selected_real_kind(15, 307)
    !integer, parameter :: steps = 200     ! Update with the desired number of steps    
    !INTEGER, PARAMETER :: num_orbitals = 2
!
    !COMPLEX(8), ALLOCATABLE :: gf_retarded(:, :, :,:), gf_lesser(:, :, :, :), se_lesser(:, :, :, :), lesser_hybridisation(:,:,:,:)
    !real(8), allocatable :: m(:,:)
!
!
    !! Allocate memory for the complex matrix
    !ALLOCATE(gf_retarded(steps, 2, num_orbitals, num_orbitals))
    !ALLOCATE(gf_lesser(steps, 2, num_orbitals, num_orbitals))
    !ALLOCATE(se_lesser(steps, 2, num_orbitals, num_orbitals))
    !ALLOCATE(lesser_hybridisation(steps, 2, num_orbitals, num_orbitals))
    !allocate(m(num_orbitals, num_orbitals))
!
!
    !call initialize_and_invert(m, num_orbitals)
    !!call initialize_mat(a, b, steps, num_orbitals)
    !!call multiple_mat(c, a, b, steps)
    !!call PrintMatrix(c, steps, num_orbitals)
    !call ReadDataFromFile("gf_retarded_", gf_retarded, num_orbitals)
    !call ReadDataFromFile("gf_lesser_", gf_lesser, num_orbitals)
    !call ReadDataFromFile("se_lesser_", se_lesser, num_orbitals)
!
    !call get_lesser_hybridisation(gf_retarded, gf_lesser, se_lesser, lesser_hybridisation, num_orbitals, steps)
    !
    !DEALLOCATE(gf_retarded, gf_lesser, se_lesser, lesser_hybridisation)
  END PROGRAM ComplexMatrixExample
  





!program main
!    implicit none
!
!    !TYPE :: complex_num
!    !    REAL, DIMENSION(2) :: data
!    !END TYPE complex_num
!
!    !parameter means const. real is double
!    integer, parameter :: dp = selected_real_kind(15, 307)
!    integer, parameter :: num_orbitals = 2 ! Update with the desired number of orbitals
!    integer, parameter :: steps = 200     ! Update with the desired number of steps
!    integer :: i 
!    integer :: j
!
!    real(dp), dimension(steps) :: energy
!    !real, dimension(2)
!    complex:: x,y, j1
!    x = (1, 1)
!    y = (1, -1)
!    j1 = (0, 1)
!
!
!    COMPLEX(8), ALLOCATABLE :: A(:,:) ! Define a complex matrix
!      
!        ! Allocate memory for the complex matrix
!    ALLOCATE(A(num_orbitals, num_orbitals))
!      
!    do i = 1, num_orbitals
!        do j = 1, num_orbitals
!          A(i, j) = CMPLX(REAL(i) + REAL(j), REAL(i) - REAL(j), 8) ! Create a complex number
!          write(*, *) A(i, j)
!        end do
!    end do
!        ! You can set values in A as needed
!      
!        ! Perform operations using LAPACK (e.g., matrix factorization)
!        ! For example, you can use zgetrf for LU factorization
!      
!        ! Deallocate the matrix when done
!    DEALLOCATE(A)
!        
!
!      
!
!    write(*,*)x+y
!    !5 format(2x,2f9.5)
!    write(*,*) conjg(x),j1*x*y
!    !6 format(2x,2f9.5,2x,2f9.5)
!    write(*,*) real(y)
!    !7 format (2x,f9.5)
!
!    print *, num_orbitals, steps
!    ! Declare a variable of type complex_num and initialize it
!    ! Define the complex_num type
!
!    CALL InitializeComplexVariable()
!  ! Declare a subroutine to demonstrate initializing the complex_num variable
!  CONTAINS
!
!  SUBROUTINE InitializeComplexVariable()
!    TYPE(complex_num) :: z
!
!    ! Initialize the complex variable
!    z%data(1) = 3.0  ! Real part
!    z%data(2) = 4.0  ! Imaginary part
!
!    ! Print the complex number
!    PRINT *, "Complex Number:", z%data(1), " + ", z%data(2), "i"
!  END SUBROUTINE InitializeComplexVariable
!  
!END PROGRAM main
!
!! Call the subroutine to demonstrate variable initialization
!
!
!    ! Read data
!    !call read_object(steps, num_orbitals, gf_retarded, "gf_retarded", energy)
!    !call read_object(steps, num_orbitals, gf_lesser, "gf_lesser", energy)
!    !call read_object(steps, num_orbitals, se_lesser, "se_lesser", energy)
!!
!    !! Calculate hybridisation
!    !call get_hybridisation(steps, num_orbitals, gf_retarded, gf_lesser, se_lesser, hybridisation_lesser)
!!
!    !! Print results to files
!    !call print_to_file(steps, num_orbitals, hybridisation_lesser, "hybridisation_lesser", energy)
!
!!contains
!!
!!    !subroutine read_object(steps, num_orbitals, object, filename, energy)
!!        !integer, intent(in) :: steps, num_orbitals
!!        !complex(dp), intent(out) :: object(steps, num_orbitals, num_orbitals)
!!        !real(dp), intent(out) :: energy(steps)
!!        !character(len=*), intent(in) :: filename
!!        !complex(dp) :: j1
!!        !integer :: i, j, line_count
!!        !character(256) :: line
!!        !real(dp) :: num1, num2, num3
!!        !character(256) :: var
!!        !complex(dp), dimension(steps) :: object_temp
!!
!!        !j1 = sqrt(-1.0_dp)
!!
!!        !do i = 1, num_orbitals
!            !do j = 1, num_orbitals
!                !var = trim(filename) // '_' // trim(itoa(i)) // '_' // trim(itoa(j)) // '.dat'
!                !print *, "reading file ", var
!                !line_count = 0
!                !open(unit=10, file=var, status='old')
!                !do
!                    !read(10, *, iostat=ios) num1, num2, num3
!                    !if (ios /= 0) exit
!                    !line_count = line_count + 1
!                    !object_temp(line_count) = num2 + num3 * j1
!                    !energy(line_count) = num1
!                !end do
!                !close(10)
!                !do r = 1, steps
!                    !object(r, i, j) = object_temp(r)
!                !end do
!            !end do
!!        !end do
!!    !end subroutine read_object
!!
!!    !subroutine get_hybridisation(steps, num_orbitals, gf_retarded, gf_lesser, se_lesser, hybridisation_lesser)
!!        !integer, intent(in) :: steps, num_orbitals
!!        !complex(dp), intent(inout) :: gf_retarded(steps, 2, num_orbitals, num_orbitals), gf_lesser(steps, 2, num_orbitals, num_orbitals)
!!        !complex(dp), intent(in) :: se_lesser(steps, 2, num_orbitals, num_orbitals)
!!        !complex(dp), intent(inout) :: hybridisation_lesser(steps, 2, num_orbitals, num_orbitals)
!!        !integer :: r, spin
!!        !complex(dp), dimension(num_orbitals, num_orbitals) :: gf_advanced, gf_advanced_inverse, gf_retarded_inverse
!!
!!        !do r = 1, steps
!            !do spin = 1, 2
!                !! Use the LAPACK routine DSYEV
!                !call DSYEV('V', 'U', num_orbitals, gf_retarded(r, spin, :, :), num_orbitals, energy, info)
!!
!                !gf_advanced = adjoint(gf_retarded(r, spin, :, :))
!                !gf_advanced_inverse = matmul(gf_advanced, 1.0_dp / det(gf_advanced))
!                !gf_retarded_inverse = matmul(gf_retarded(r, spin, :, :), 1.0_dp / det(gf_retarded(r, spin, :, :)))
!!
!                !hybridisation_lesser(r, spin, :, :) = matmul(matmul(gf_retarded_inverse, gf_lesser(r, spin, :, :)), gf_advanced_inverse) - hybridisation_lesser(r, spin, :, :)
!            !end do
!!        !end do
!!    !end subroutine get_hybridisation
!!
!!    !subroutine print_to_file(steps, num_orbitals, object, filename, energy)
!!        !integer, intent(in) :: steps, num_orbitals
!!        !complex(dp), intent(in) :: object(steps, num_orbitals, num_orbitals)
!!        !real(dp), intent(in) :: energy(steps)
!!        !character(len=*), intent(in) :: filename
!!        !character(256) :: var
!!        !integer :: i, j, r
!!        !open(unit=10, file=trim(filename) // '_00.dat', status='replace')
!!        !do i = 1, num_orbitals
!            !do j = 1, num_orbitals
!                !var = trim(filename) // '_' // trim(itoa(i)) // '_' // trim(itoa(j)) // '.dat'
!                !open(unit=10, file=var, status='replace')
!                !do r = 1, steps
!                    !write(10, *) energy(r), real(object(r, i, j)), aimag(object(r, i, j))
!                !end do
                !close(10)
            !end do
!        !end do
!    !end subroutine print_to_file
!

!