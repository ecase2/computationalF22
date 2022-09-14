! ===========================================================================================
module parameters
implicit none
   REAL, PARAMETER  	:: b = 0.99, d = 0.025, a = 0.36
   REAL, PARAMETER  	:: klb = 0.01, inc = 0.025, kub = 75.0
   REAL, PARAMETER      :: zg = 1.25, zb = 0.2
   REAL, PARAMETER      :: gg = 0.977, gb = 0.023, bg = 0.074, bb = 0.926
   INTEGER, PARAMETER 	:: length_grid_k = (kub-klb)/inc + 1
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
   REAL    :: diff, k, kp, c, z, pone, ptwo
   
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
                ptwo = gg
                pone = gb
            endif 
            if (index_z == 2) then 
                z  = zb 
                pone = bb 
                ptwo = bg
            endif
            do index_k = 1, length_grid_k				! Capital grid
                k = Kgrid(index_k)
                            vtmp(index_k,:, index_z) = -1.0e-16

                do index_kp = 1, length_grid_k
                    kp = Kgrid(index_kp)
                    c = z*k**a+(1.-d)*k-kp

                    if (c>0.) then
                        vtmp(index_k, index_kp, index_z) = log(c)+b*value(index_kp, index_z)*pone + b*value(index_kp, index_z)*ptwo
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
    
  !  open (UNIT=1,FILE='polfun.dat',STATUS='replace')
!	do index_k = 1, length_grid_k
 !           WRITE(UNIT=1,FMT=*) g_k_g(index_k), g_k_b(index_k)
  !  end do
   ! close (UNIT=1)
    
end subroutine 