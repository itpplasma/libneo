module libneo_collisions
   use, intrinsic :: iso_c_binding
   use libneo_species, only: species_t

   implicit none

   integer, parameter :: dp = kind(1.0d0)

   interface
      function gsl_sf_gamma(x) bind(c, name='gsl_sf_gamma')
         import :: c_double
         implicit none
         real(c_double), value :: x
         real(c_double) :: gsl_sf_gamma
      end function gsl_sf_gamma
   end interface

   interface
      function gsl_sf_gamma_inc_P(a, x) bind(c, name='gsl_sf_gamma_inc_P')
         import :: c_double
         implicit none
         real(c_double), value :: a, x
         real(c_double) :: gsl_sf_gamma_inc_P
      end function gsl_sf_gamma_inc_P
   end interface

contains

   subroutine fill_species_arr_coulomb_log(num_species, species_arr)

      implicit none
      integer, intent(in) :: num_species
      type(species_t), intent(inout) :: species_arr(num_species)

      integer :: i, j
      integer :: interaction_type

      if (.not. allocated(species_arr(1)%coulomb_log)) then
         do i = 1, num_species
            allocate (species_arr(i)%coulomb_log(num_species))
         end do
      end if

      do i = 1, num_species
         do j = i, num_species
            interaction_type = species_arr(i)%typ_int + species_arr(j)%typ_int
            call calc_coulomb_log(interaction_type, species_arr(i), species_arr(j), species_arr(i)%coulomb_log(j))
            if (i /= j) species_arr(j)%coulomb_log(i) = species_arr(i)%coulomb_log(j)
         end do
      end do

   end subroutine

   subroutine calc_coulomb_log(interaction_type, species_a, species_b, coulomb_log)
      ! determines the Coulomb logarithm for a given interaction type
      ! units are in cgs (temperature is in eV)
      ! interaction_type (integer): "ee" = 2, "ei" = 3, "ie" = 3, "ii" = 4
      ! species a is electrons for ee, ei and ie; species b is ions

      use math_constants, only: m_e, m_p
      use libneo_species, only: species_t

      implicit none

      integer, parameter :: EE = 2, EI = 3, II = 4 ! or IE = 3

      integer, intent(in) :: interaction_type
      type(species_t), intent(in) :: species_a, species_b
      real(dp), intent(out) :: coulomb_log

      select case (interaction_type)
      case (EE)
         coulomb_log = 23.5d0 - log(sqrt(species_a%dens)*species_a%temp**(-5.0d0/4.0d0)) - &
                       sqrt(10.0d0**(-5.0d0) + (log(species_a%temp) - 2.0d0)**2.0d0/16.0d0)
      case (EI) ! or IE
         if (species_b%temp*m_e/species_b%mass < species_a%temp &
             .and. species_a%temp < 10.0d0*species_b%charge_num**2.0d0) then
            coulomb_log = 23.0d0 - log(sqrt(species_a%dens)*species_b%charge_num* &
                                       species_a%temp**(-1.5d0))
         elseif (species_b%temp*m_e/species_b%mass < 10.0d0*species_b%charge_num**2.0d0 &
                 .and. 10.0d0*species_b%charge_num**2.0d0 < species_a%temp) then
            coulomb_log = 24.0d0 - log(sqrt(species_a%dens)*species_a%temp**(-1.0d0))
         elseif (species_a%temp < species_b%temp*m_e/species_b%mass) then
            coulomb_log = 16.0d0 - log(sqrt(species_b%dens)*species_b%temp**(-1.5d0) &
                                       *species_b%charge_num**2.0d0*species_b%mass/m_p)
         end if
      case (II)
         coulomb_log = 23.0d0 - log(species_a%charge_num*species_b%charge_num &
             *(species_a%mass + species_b%mass)/(species_a%mass*species_b%temp &
             + species_b%mass*species_a%temp) &
             *sqrt(species_a%dens*species_a%charge_num**2.0d0/species_a%temp &
             + species_b%dens*species_b%charge_num**2.0d0/species_b%temp))
      case default
         print *, "Unknown interaction type"
      end select

      if (coulomb_log < 0.0d0 .or. coulomb_log > 30.0d0) then
         stop "Coulomb logarithm has unexpected value"
      end if

   end subroutine

   subroutine calc_perp_coll_freq(vel, species_a, species_b, coulomb_log, coll_freq)

      ! determines the perpendicular collision frequency between two species a and b
      ! a is the test species, b is the background species
      ! units are in cgs (temperature is in eV)

      use math_constants, only: pi, ev_to_cgs, E
      use libneo_species, only: species_t

      implicit none

      real(dp), intent(in) :: vel, coulomb_log
      type(species_t), intent(in) :: species_a, species_b
      real(dp), intent(out) :: coll_freq
      real(dp) :: nue_0, z_ab
      real(dp) :: psi_of_x
      real(dp) :: p = 1.5d0

      nue_0 = 4.0d0*pi*species_a%charge_num**2.0d0*species_b%charge_num**2.0d0*E**4.0d0 &
              *coulomb_log*species_b%dens/(species_a%mass**2.0d0*vel**3.0d0)
      z_ab = species_b%mass*vel**2.0d0/(2.0d0*species_b%temp*ev_to_cgs)
      psi_of_x = lower_incomplete_gamma(p, z_ab)

      coll_freq = 2.0d0*nue_0*((1.0d0 - 1.0d0/(2.0d0*z_ab))*psi_of_x*2.0d0/sqrt(pi) &
                               + 2.0d0/sqrt(pi)*z_ab**0.5d0*exp(-z_ab))

   end subroutine

   subroutine calc_perp_coll_freq_slow_limit_ee(vel, dens_e, temp_b, coulomb_log, coll_freq)

      use math_constants, only: m_e, ev_to_cgs

      implicit none

      real(dp), intent(in) :: vel, dens_e, temp_b, coulomb_log
      real(dp), intent(out) :: coll_freq

      coll_freq = 5.8d0*10d0**(-6d0)*dens_e*coulomb_log*temp_b**-0.5d0*(m_e*vel**2d0*0.5d0/ev_to_cgs)**(-1.0d0)

   end subroutine

   subroutine calc_perp_coll_freq_fast_limit_ee(vel, dens_e, coulomb_log, coll_freq)

      use math_constants, only: m_e, ev_to_cgs

      implicit none
      real(dp), intent(in) :: vel, dens_e, coulomb_log
      real(dp), intent(out) :: coll_freq

      coll_freq = 7.7d0*10.0d0**(-6.0d0)*dens_e*coulomb_log*(m_e*vel**2.0d0*0.5d0/ev_to_cgs)**(-1.5d0)

   end subroutine

   function lower_incomplete_gamma(a, x) result(gamma)
      real(dp), intent(in) :: a, x
      real(dp) :: gamma

      gamma = gsl_sf_gamma_inc_P(a, x)*gsl_sf_gamma(a)
   end function lower_incomplete_gamma
end module
