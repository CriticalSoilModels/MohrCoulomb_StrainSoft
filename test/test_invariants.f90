program check_invariants
    use kind_precision_module, only: dp, i32
    use mod_MCSS_funcs, only: CalculateInvariants
implicit none

integer(kind = i32) :: IntGlo
real(kind = dp)     :: Sig(6), p, J, Lode, S3TA

Sig = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]

! Define the values for the invariants
call CalculateInvariants(IntGlo, Sig, p, J, Lode, S3TA)

print *, "sigma_bar: ", calc_abbo_sloan_sigma_bar(Sig)

print *, "Mean Stress: ", p
print *, "Dev. Stress: ", J
print *, "Lode Angle : ", Lode
! This is an intermediate value to get to the lode angle.
! Quick look seems to be used to calc K  
print *, "S3TA :", S3TA


contains

 function calc_abbo_sloan_sigma_bar(sigma) result (sigma_bar)
    ! Input
    real(kind = dp) :: sigma(6)
    ! Output
    real(kind = dp) :: sigma_bar

    ! Local variables
    real(kind = dp) :: mean_stress, dev_stress(6)
    
    ! Calc the mean stress
    mean_stress = sum(sigma(1:3)) / 3.0_dp

    ! Calc the deviatoric stress
    dev_stress(:) = sigma(:)

    dev_stress(1:3) = dev_stress(1:3) - mean_stress

    ! Calc sigma_bar
    sigma_bar = sqrt( 0.5 * sum(dev_stress(1:3)**2) + sum(dev_stress(4:6)**2) )

end function calc_abbo_sloan_sigma_bar

end program check_invariants
