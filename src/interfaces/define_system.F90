module define_system
implicit none

!Define the basis set and molecule from user file mol.xyz

contains
  subroutine define_basis(ao_basis, charge, coords)
    use ao_basis

    real(8), allocatable :: charge(:), coords(:,:)
    integer :: i
    type(basis_set_info_t), intent(inout) :: ao_basis
    type(basis_func_info_t) :: gto
    do i=1, size(charge)
       if(charge(i).eq.1) then
           !Hydrogen atom                   l               coord                     exp
           call add_shell_to_basis(ao_basis,0,(/coords(1,i),coords(2,i),coords(3,i)/),0.1D0)
           call add_shell_to_basis(ao_basis,0,(/coords(1,i),coords(2,i),coords(3,i)/),1.0D0)
           call add_shell_to_basis(ao_basis,0,(/coords(1,i),coords(2,i),coords(3,i)/),3.0D0)
       else
           !Other elements
           ! s-orbitals                     l               coord                     exp
           call add_shell_to_basis(ao_basis,0,(/coords(1,i),coords(2,i),coords(3,i)/),0.1D0)
           call add_shell_to_basis(ao_basis,0,(/coords(1,i),coords(2,i),coords(3,i)/),0.35D0)
           call add_shell_to_basis(ao_basis,0,(/coords(1,i),coords(2,i),coords(3,i)/),1.0D0)
           call add_shell_to_basis(ao_basis,0,(/coords(1,i),coords(2,i),coords(3,i)/),3.0D0)
           call add_shell_to_basis(ao_basis,0,(/coords(1,i),coords(2,i),coords(3,i)/),10.0D0)
           ! p-orbitals                     l               coord                     exp
           call add_shell_to_basis(ao_basis,1,(/coords(1,i),coords(2,i),coords(3,i)/),0.2D0)
           call add_shell_to_basis(ao_basis,1,(/coords(1,i),coords(2,i),coords(3,i)/),1.0D0)
           call add_shell_to_basis(ao_basis,1,(/coords(1,i),coords(2,i),coords(3,i)/),5.0D0)
           ! d-orbitals                     l               coord                     exp
           call add_shell_to_basis(ao_basis,2,(/coords(1,i),coords(2,i),coords(3,i)/),1.0D0)
       end if
    end do

  end subroutine

  subroutine read_xyz(molecule, ao_basis, n_occ, conv)
     use ao_basis
     use molecular_structure
     implicit none
     integer :: i, len, io, loc, j
     integer, intent(out) :: n_occ
     type(basis_set_info_t), intent(inout) :: ao_basis
     real(8), intent(out) :: conv
     type(basis_func_info_t) :: gto

     type(molecular_structure_t), intent(inout) :: molecule
     
     character(2), allocatable :: atoms(:)
     character(2) :: table(10)
     real(8), allocatable :: charge(:), coords(:,:)

     table = (/"H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne"/)

     open(10, file="mol.xyz")

     !Determine the length of the xyz file first
     len = 0
     do
       read(10, *, iostat=io)
       if(io/=0) exit
       len = len + 1
     end do

     allocate(charge(len-2))
     allocate(atoms(len-2))
     allocate(coords(3,len-2))
     
     !Discard first line of xyz file and read convergence criterion
     rewind 10
     read(10,*)
     read(10,*) conv
     !Read in atom element and (x,y,z) coordinates
     do i=1, len-2
       read(10,*) atoms(i), coords(1,i), coords(2,i), coords(3,i)
       charge(i) = findloc(table, atoms(i), 1)
     end do
    
     print*, "System:"
     do i=1, len-2
        print*, atoms(i)
        print*, coords(:,i)
     end do

     if(mod(sum(charge),2.0) .eq. 0) then
        n_occ = sum(charge)/2
     else
        stop 'Please submit a closed-shell system.'
     end if


     call add_atoms_to_molecule(molecule, charge, coords)
     call define_basis(ao_basis, charge, coords)

  end subroutine

end module