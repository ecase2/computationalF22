module parameters
implicit none
   REAL, PARAMETER  	:: b = 0.99, d = 0.025, a = 0.36
   REAL, PARAMETER  	:: klb = 0.01, inc = 0.025, kub = 75.0
   REAL, PARAMETER      :: zg = 1.25, zb = 0.2
   REAL, PARAMETER      :: gg = 0.977, gb = 0.023, bg = 0.074, bb = 0.926
   INTEGER, PARAMETER 	:: length_grid_k = (kub-klb)/inc + 1
   REAL , PARAMETER 	:: toler   = 1.e-4						! Numerical tolerance
end module

! ============================================================================================
module global
USE parameters
implicit none
   REAL 		        :: Kgrid(length_grid_k), valueg(length_grid_k), valueb(length_grid_k), g_k_g(length_grid_k), g_k_b(length_grid_k)
   REAL                 :: vtmpg(length_grid_k,length_grid_k), vtmpb(length_grid_k,length_grid_k), value_newb(length_grid_k), value_newg(length_grid_k)
end module

! ============================================================================================

PROGRAM  HW2Stochastic
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
   REAL    :: diff, k, kp, c, z, value, value_new, vtmp, g_k, pone, ptwo
   
   INTEGER :: i = 1


   do while (i<=length_grid_k)   !do loop for assigning capital grid K
     Kgrid(i) = klb + (i-1)*inc
     !write(*,*) i, Kgrid(i)
     i = i + 1
   end do


   iter = 1
   diff = 1000.d0
   valueb = 0.*Kgrid		!Initial Value guess
   valueg = 0.*Kgrid
   
	do while (diff>= toler)
        do index_z = 1, 2
            if (index_z == 1) then 
                z  = zg 
                ptwo = gg
                pone = gb
                vtmp = vtmpg
                value = valueg
                value_new = value_newg
                g_k = g_k_g
            endif 
            if (index_z == 2) then 
                z  = zb 
                pone = bb 
                ptwo = bg
                vtmp = vtmpb
                value = valueb
                value_new = value_newb
                g_k = g_k_b
            endif
            do index_k = 1, length_grid_k				! Capital grid
                k = Kgrid(index_k)
                            vtmp(index_k,:) = -1.0e-16

                do index_kp = 1, length_grid_k
                    kp = Kgrid(index_kp)
                    c = z*k**a+(1.-d)*k-kp

                    if (c>0.) then
                                    vtmp(index_k,index_kp) = log(c)+b*valueb(index_kp)*pone + b*valueg(index_kp)*ptwo
                        endif

                enddo

                value_new(index_k) = MAXVAL(vtmp(index_k,:))
                            g_k(index_k) 	   = Kgrid(MAXLOC(vtmp(index_k,:),1))
            enddo

            diff  = maxval(abs(value_new-value))/ABS(value_new(length_grid_k))
            if (index_z == 1) then 
                valueg = value_new
            endif
            if (index_z == 2) then 
                valueb = value_new
            endif 


            print*, 'Iteration =',iter,'sup_norm =',diff
                    iter = iter+1

        enddo
    enddo

	print *, ' '
	print *, 'Successfully converged with sup_norm ', diff
    	!print *, g_k
    
    !CALL vcDrawCurve@(d, Kgrid, g_k, length_grid_k)


        open (UNIT=1,FILE='valuefun',STATUS='replace')
	do index_k = 1, length_grid_k
        	WRITE(UNIT=1,FMT=*) valueb(index_k)
            WRITE(UNIT=1,FMT=*) valueg(index_k)
        end do
	close (UNIT=1)

end subroutine
