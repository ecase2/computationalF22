! ECON899. Problem set 1.
! The code solves the neoclassical growth model with stochastic productivity.
! It is a modified version of the script "hw2ns_pablo.f95" from Dean Corbae's website. 
! Modifications are introduced by Emily Case, Hanna Han, Anna Lukianova (09/13/2022).

    
! ===========================================================================================
module parameters
implicit none
   REAL, PARAMETER  	:: b = 0.99, d = 0.025, a = 0.36
   REAL, PARAMETER  	:: klb = 0.01, inc = 0.025, kub = 75.
   INTEGER, PARAMETER 	:: length_grid_k = (kub-klb)/inc + 1
   REAL, PARAMETER      :: zg = 1.25, zb = 0.2
   REAL, PARAMETER      :: Pgg = 0.977, Pbg = 0.023, Pgb = 0.074, Pbb = 0.926 
   REAL , PARAMETER 	:: toler   = 1.e-4						! Numerical tolerance
end module parameters 

! ============================================================================================
module global
use  parameters
implicit none
   REAL     :: Kgrid(length_grid_k), value(length_grid_k, 2), g_k(length_grid_k, 2)
   REAL     :: vtmp(length_grid_k, length_grid_k, 2)
   REAL     :: value_new(length_grid_k, 2)
end module global 

! ============================================================================================

PROGRAM  HW2Stochastic
   use parameters
   use global 
   
   REAL 		        	:: total, etime, dist
   REAL, DIMENSION(2)  		:: elapsed

   call solution

   total=etime(elapsed)


	PRINT*,'--------------------------------------------------'
	PRINT*,'total time elpased =',total
	PRINT*,'--------------------------------------------------'


END PROGRAM HW2Stochastic

! ============================================================================================
subroutine solution
USE parameters
USE global

   IMPLICIT  NONE

   INTEGER :: iter, index_k, index_kp, index_z
   REAL    :: diff, k, kp, c, z, pr1, pr2
   
   INTEGER :: i = 1


   do while (i<=length_grid_k)   !do loop for assigning capital grid K
     Kgrid(i) = klb + (i-1)*inc
     !write(*,*) i, Kgrid(i)
     i = i + 1
   end do
  
    
   iter = 1
   diff = 1000.d0
   value(:, 1) = 0.*Kgrid		!Initial Value guess
   value(:, 2) = 0.*Kgrid
   
	do while (diff>= toler)
        do index_z = 1, 2
            if (index_z == 1) then 
                z  = zg 
                pr1 = Pgg 
                pr2 = Pbg 
            endif 
            if (index_z == 2) then 
                z  = zb 
                pr1 = Pgb 
                pr2 = Pbb 
            endif
            do index_k = 1, length_grid_k				! Capital grid
                k = Kgrid(index_k)
                            vtmp(index_k,:, index_z) = -1.0e-16

                do index_kp = 1, length_grid_k
                    kp = Kgrid(index_kp)
                    c = z*k**a+(1.-d)*k-kp

                    if (c>0.) then
                            vtmp(index_k, index_kp, index_z) = log(c)+b*(pr1*value(index_kp, 1)+pr2*value(index_kp, 2))
                    endif

                enddo

                value_new(index_k, index_z) = MAXVAL(vtmp(index_k,:, index_z))
                            g_k(index_k, index_z) 	   = Kgrid(MAXLOC(vtmp(index_k,:, index_z),1))
            enddo

            diff  = maxval(abs(value_new-value))/ABS(value_new(length_grid_k,1))
            value = value_new 


        !    print*, 'Iteration =',iter,'sup_norm =',diff
          !          iter = iter+1

        enddo
    enddo

	print *, ' '
	print *, 'Successfully converged with sup_norm ', diff
    	!print *, g_k
    
    !CALL vcDrawCurve@(d, Kgrid, g_k, length_grid_k)


    open (UNIT=1,FILE='fun.dat',STATUS='replace')
	do index_k = 1, length_grid_k
        	WRITE(UNIT=1,FMT=*) value(index_k, 1), value(index_k, 2), g_k(index_k, 1), g_k(index_k, 2)

        end do
	close (UNIT=1)

    
end subroutine 