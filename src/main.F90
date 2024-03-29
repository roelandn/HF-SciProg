program HartreeFock

   ! Hartree-Fock implementation
   ! Roeland Neugarten, March 2024

   use molecular_structure
   use ao_basis
   use compute_integrals
   use diagonalization
   use define_system
     implicit none

     ! Variable containing the molecular structure
     type(molecular_structure_t) :: molecule
     ! Variable containing the atomic orbital basis
     type(basis_set_info_t) :: ao_basis

     integer  :: n_AO, n_occ
     integer  :: kappa, lambda, iteration
     real(8)  :: start, finish, start_step, finish_step
     real(8)  :: E_HF
     real(8), allocatable :: F(:,:),V(:,:),T(:,:),S(:,:), C(:,:), eps(:), D(:,:), hcore(:,:)

    real(8) :: Dn, Dnm1=0.D0, conv

     ! The following large array can be eliminated when Fock matrix contruction is implemented
     real(8), allocatable :: ao_integrals (:,:,:,:)
     integer, parameter :: max_iterations = 50
     
     call cpu_time(start)

     ! Definition of the molecule and GTOs
     call read_xyz(molecule, ao_basis, n_occ, conv)

     n_AO = ao_basis%nao

     print*, "N_occ = ", n_occ
     print*, "N_AO  = ", n_AO
     print*, ""
     print*, "Convergence requirement = ", conv

     ! Compute the overlap matrix
     allocate (S(n_AO,n_AO))
     call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)

     ! Compute the kinetic matrix
     allocate (T(n_AO,n_AO))
     call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)

     ! Compute the potential matrix
     allocate (V(n_AO,n_AO))
     call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)

     ! Compute the core Hamiltonian matrix (the potential is positive, we scale with -e = -1 to get to the potential energy matrix)
     allocate (hcore(n_AO,n_AO))
     allocate (F(n_AO, n_AO))
     hcore = T - V
     F = hcore

     allocate (C(n_AO,n_AO))
     allocate (eps(n_AO))

     allocate (D(n_AO,n_AO))

     allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))

    ! Compute all 2-electron integrals
     call generate_2int (ao_basis,ao_integrals)

     write(*,'(A5, 2A25, A15)') "i", "E", "dDn", "CPU time"
    do iteration = 1, max_iterations

      call cpu_time(start_step)

     ! Diagonalize the Fock matrix
     call solve_genev (F,S,C,eps)

     ! Form the density matrix
     D = 0
     do lambda = 1, n_ao
        do kappa = 1, n_ao
           D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))
       end do
     end do
     
    !Compute Fock matrix
    F = hcore
    do lambda = 1, n_ao
      do kappa = 1, n_ao
        F(kappa, lambda) = F(kappa, lambda) + 2.D0 * sum(D*ao_integrals(kappa, lambda, :, :))
        F(kappa, lambda) = F(kappa, lambda) - 1.D0 * sum(D*ao_integrals(kappa, :, :, lambda))  
      end do
    end do

    E_HF = sum((hcore+F)*D) 
  
    !Compute convergence with norm of [F,D] commutator
    Dn = sqrt(sum((matmul(F,D)-matmul(D,F))**2))

    if(abs(Dn-Dnm1).le.conv) then
      write(*,*) "Converged with dDn = ", Dn
      exit
    else
      call cpu_time(finish_step)
      write(*,'(I5,2F25.13,F15.6)') iteration, E_HF, Dn, finish_step-start_step
      Dnm1 = Dn
    end if
    
  end do

  call cpu_time(finish)
  print*, "Calculation time:           ", finish-start, "seconds"
  print*, "The Hartree-Fock energy:    ", E_HF
   end
