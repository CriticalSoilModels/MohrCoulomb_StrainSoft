program check_yield_function
    use kind_precision_module, only: dp, i32
    use mod_MCSS_funcs, only: DetermineYieldFunctionValue
    
    implicit none
    
call DetermineYieldFunctionValue(IntGlo, Sig, c, phi, F)

print *, "Yield function value"
end program check_yield_function