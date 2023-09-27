  PROGRAM ComplexMatrixExample
    use MyModule
       
    implicit none
    integer, parameter :: steps = 200     ! Update with the desired number of steps    
    integer, parameter :: num_orbitals = 2
    complex*16, allocatable :: gf_retarded(:, :, :,:), gf_lesser(:, :, :, :), se_lesser(:, :, :, :), lesser_hybridisation(:,:,:,:)
    real(8), dimension(:), allocatable :: energy

    allocate(gf_retarded(steps, 2, num_orbitals, num_orbitals), gf_lesser(steps, 2, num_orbitals, num_orbitals))
    allocate(se_lesser(steps, 2, num_orbitals, num_orbitals), lesser_hybridisation(steps, 2, num_orbitals, num_orbitals))
    allocate(energy(steps))

    !this is reading the texts which I generated using my multiple orbital sigma 2 code
    call ReadDataFromFile("gf_retarded_", gf_retarded, num_orbitals, energy)
    call ReadDataFromFile("gf_lesser_", gf_lesser, num_orbitals, energy)
    call ReadDataFromFile("se_lesser_", se_lesser, num_orbitals, energy)

    !this gets the lesser hybridisation
    call get_lesser_hybridisation(gf_retarded, gf_lesser, se_lesser, lesser_hybridisation, num_orbitals, steps)

    !this subroutine prints the lesser hybridisation to file
    call WriteToFile("lesser_hybridisation_", lesser_hybridisation, num_orbitals, steps, energy)
    deallocate(gf_retarded, gf_lesser, se_lesser, lesser_hybridisation, energy)
  END PROGRAM ComplexMatrixExample
  
