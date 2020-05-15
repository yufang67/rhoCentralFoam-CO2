MODULE Transprop
!
!
     USE def_variables     
     USE def_constants
     USE properties
!     USE interp_functions
!     USE non_linear_solvers
!      
     IMPLICIT NONE
!
     PRIVATE
     PUBLIC  CO2TRANS_PROP
!
!
     CONTAINS
!
!===============================================================================
!
SUBROUTINE CO2TRANS_PROP(cp_out, vis_out, dummy, dummy2, dummy3, dummy4, flagloc,dummy5,dummy6) 
!
!===============================================================================
!
!
! Input: u_in (the specific internal energy)
!        v_in (the specific volume v_in) 
! Output: cp_out  (heat capacity) 
!         vis_out (viscosity)
!
!===============================================================================
REAL(pr),INTENT(OUT)  :: cp_out, vis_out, dummy, dummy2, dummy3, dummy4
!INTEGER, INTENT(OUT)  ::
INTEGER,INTENT(IN)   :: flagloc
!
!Local Variables
INTEGER :: i, j, flag_TP, j_sat, i_R, i_L, flag_loca,Niter,exitflag
REAL(pr) :: v_min, v_max, v_sat, v_sat_log, delta,&
&            qual,press,temp,sound,T_guess,out2
REAL(pr) :: duL_dp,  duV_dp,  dvL_dp,  dvV_dp
REAL(pr) :: vL, uL, vV, uV, pp, ratio, du_dp_x, dv_dp_x, entrov, entrol
REAL(pr) :: x_check, v_check, f_out
REAL(pr) :: gamma_pg, e_pg, v_in, u_in
REAL(pr) :: resnorm,out3,T_gs, p_gs, c_gs


x_out   = 1_pr
a_out   = 1_pr
!res     = 0_pr
f_out   = 0_pr
!
!
SELECT CASE (flag_loca)
!
CASE( 0 )
STOP '** Locating the points in the physical domaion failed '
!
!###### LL #####
CASE( 1 )
     CALL Lin_int_Left_Low(T_out, p_out, c_out, u_in, v_in)

! sound speed correction
!
     IF ( (u_in >-200.0e3_pr) .AND. (v_in >1.9e-3_pr) ) THEN
        CALL sound_speed(T_out,v_in,c_out)        
     ENDIF 
!
     CALL entropy(T_out, v_in, entro)
     x_out   = 0_pr
     a_out   = 0_pr
!
!###### LH #####
CASE( 2 )
!
!   IF ( (u_in > e_F) .OR. ((u_in < e_G) .AND. (v_in>v_G)) ) THEN
    IF ( u_in > e_F) THEN
    CALL Lin_int_Left_High(T_gs, p_gs, c_gs, u_in, v_in)
!   
    CALL  New_Rap1D(1, T_out, out2, resnorm, Niter,&
     &              exitflag, u_in,T_gs,v_in, out3) 
    CALL  pressure(T_out,v_in,p_out)
    CALL  sound_speed(T_out,v_in,c_out)
!
   ELSE 
     CALL Lin_int_Left_High(T_out, p_out, c_out, u_in, v_in)
   ENDIF
! speed sound correction
   IF ( (u_in < e_G) .AND. (v_in >v_G) ) THEN
        CALL sound_speed(T_out,v_in,c_out)
   ENDIF
!   
   IF ( (p_out > p_cr) .AND. (T_out > T_cr) ) THEN
       x_out   = 0_pr
       a_out   = 0_pr
   ENDIF
   CALL entropy(T_out, v_in, entro)
!
!
!###### R ######
CASE( 3 )
     delta        = (x_mesh_max - x_mesh_min)/(MMM_R-1)
     x_mesh_R     = x_mesh_min + (/(i*delta, i=0,MMM_R-1)/)
     CALL Lin_int(T_out, p_out, c_out, u_in, v_in, NNN_R, MMM_R, x_mesh_R, y_mesh_R, &
&         spline_pmin, spline_Vsat, TTT_R, ppp_R, ccc_R) 

     CALL entropy(T_out, v_in, entro)
!## HT ## 
CASE( 4 )
     delta       = (x_mesh_max - x_mesh_min)/(MMM_HT-1)
     x_mesh_HT   =  x_mesh_min + (/(i*delta, i=0,MMM_HT-1)/)
     CALL Lin_int_Log10(T_out, p_out, c_out, u_in, v_in, NNN_HT, MMM_HT, x_mesh_HT, y_mesh_HT, &
&         spline_right_HT, spline_left_HT, TTT_HT, ppp_HT, ccc_HT)
!
!
   IF ( (p_out > p_cr) .AND. (T_out > T_cr) ) THEN
       x_out   = 0_pr
       a_out   = 0_pr
   ENDIF
   CALL entropy(T_out, v_in, entro)
!
!
!######## TP #########
CASE( 5 )
!@@@@@@@@@@ TPH @@@@@@@@@@
   IF (u_in .GE. e_tri_R) THEN  
!
         delta       = (x_mesh_max - x_mesh_min)/(MMM_TPH-1)
         x_mesh_TPH  =  x_mesh_min + (/(i*delta, i=0,MMM_TPH-1)/)
         CALL Lin_TP_Log10(T_out, p_out, c_out, x_out, u_in, v_in, NNN_TPH, MMM_TPH,&
&             x_mesh_TPH, y_mesh_TPH,spline_Vsat, spline_left_TPH,&
&             TTT_TPH, ppp_TPH, ccc_TPH, xxx_TPH)
!
!         CALL New_Rap1D(2, press, qual, res, Niter,&
!&                       exitflag, u_in, p_out, v_in, temp)
!
!         p_out = press
!         T_out = temp
!         x_out = qual

        CALL satprop (3, p_out, dummy, vV, vL, uV, uL)
!        CALL satderiv(3, p_out, duL_dp, duV_dp, dvL_dp, dvV_dp)
!!
        x_out   = (v_in - vL) / (vV - vL)
!        ratio   = ((uV - uL)/(vV - vL))   ! (J/kg)/(m3/kg)
!        du_dp_x = x_out * duV_dp + (1_pr - x_out) * duL_dp
!        dv_dp_x = x_out * dvV_dp + (1_pr - x_out) * dvL_dp
!!
!       c_out = SQRT((p_out  + ratio)/(du_dp_x - ratio * dv_dp_x)) * v_in ! (m/s)
!!         
        a_out = x_out*vV/v_in
!!
         IF (((a_out .GT. 1_pr) .AND. (x_out .GT. 8e-1_pr)) .OR. (x_out .GT. a_out)) THEN
            a_out = x_out
         ENDIF
         CALL entropy(T_out, vV, entrov)
         CALL entropy(T_out, vL, entrol)
         entro = x_out*entrov + (1.0-qual)*entrol
!@@@@@@@@@@@@@ TPL @@@@@@@@@@@@@
   ELSEIF (u_in .LE. e_cr) THEN
!
         delta       = (x_mesh_max - x_mesh_min)/(MMM_TPL-1)
         x_mesh_TPL  =  x_mesh_min + (/(i*delta, i=0,MMM_TPL-1)/)
         CALL Lin_TP_Log10(T_out, p_out, c_out, x_out, u_in, v_in, NNN_TPL, MMM_TPL,&
&             x_mesh_TPL,y_mesh_TPL,spline_right_TPL, spline_Lsat_LL,&
&             TTT_TPL, ppp_TPL, ccc_TPL, xxx_TPL)
!      
!      IF( (u_in < -199e3_pr) .OR.(v_in > 2.34e-3_pr) )THEN 
!         CALL New_Rap1D(2, press, qual, res, Niter,&
!&                       exitflag, u_in, p_out, v_in, temp)
!!
!         IF (res > 1_pr) THEN
!            print*, 'TPL interpolation guess prob', res 
!!           print*, 'u,v',u_in,v_in, 'p_iter',press, 'p_tb', p_out
!        ENDIF
!        p_out = press
!!       T_out = temp
!!       x_out = qual
         CALL satprop (3, p_out, dummy, vV, vL, uV, uL)
!        CALL satderiv(3, p_out, duL_dp, duV_dp, dvL_dp, dvV_dp)
!!
         x_out   = (v_in - vL) / (vV - vL)
!
!       ratio = ((uV - uL)/(vV - vL))   ! (J/kg)/(m3/kg)
!       du_dp_x = x_out * duV_dp + (1_pr - x_out) * duL_dp
!!      dv_dp_x = x_out * dvV_dp + (1_pr - x_out) * dvL_dp
!
!       c_out = SQRT((p_out  + ratio)/(du_dp_x - ratio * dv_dp_x)) * v_in ! (m/s)
!!         
         a_out = x_out*vV/v_in
!!
         IF (((a_out .GT. 1_pr) .AND. (x_out .GT. 8e-1_pr)) .OR. (x_out .GT. a_out)) THEN
            a_out = x_out
         ENDIF
         CALL entropy(T_out, vV, entrov)
         CALL entropy(T_out, vL, entrol)
         entro = x_out*entrov + (1.0-qual)*entrol
!      ENDIF
!@@@@@@@@@@@@@ TPM @@@@@@@@@@@@@@@@
   ELSE  
!
         delta       = (x_mesh_max - x_mesh_min)/(MMM_TPM-1)
         x_mesh_TPM  =  x_mesh_min + (/(i*delta, i=0,MMM_TPM-1)/)
         CALL Lin_TP_Log10(T_out, p_out, c_out, x_out, u_in, v_in, NNN_TPM, MMM_TPM,&
&             x_mesh_TPM,y_mesh_TPM,spline_right_TPM, spline_left_TPM,&
&             TTT_TPM, ppp_TPM, ccc_TPM, xxx_TPM)
!
!       IF ( (u_in > -183e3_pr) .OR. (v_in >2.3408e-3_pr) ) THEN    
!         CALL New_Rap1D(2, press, qual, res, Niter,&
!     &                exitflag, u_in, p_out, v_in, temp)
!
!         IF (res > 1_pr) THEN
!            print*, 'TPM interpolation guess prob', res
!            print*, 'u,v',u_in,v_in, 'p_iter',press, 'p_tb', p_out
!         ENDIF
!!        p_out = press
!         T_out = temp
!         x_out = qual
!
         CALL satprop (3, p_out, dummy, vV, vL, uV, uL)
!        CALL satderiv(3, p_out, duL_dp, duV_dp, dvL_dp, dvV_dp)
!
         x_out   = (v_in - vL) / (vV - vL)
!
!       ratio = ((uV - uL)/(vV - vL))   ! (J/kg)/(m3/kg)
!       du_dp_x = x_out * duV_dp + (1_pr - x_out) * duL_dp
!       dv_dp_x = x_out * dvV_dp + (1_pr - x_out) * dvL_dp
!
!       c_out = SQRT((p_out  + ratio)/(du_dp_x - ratio * dv_dp_x)) * v_in ! (m/s)
!!         
        a_out = x_out*vV/v_in
!
        IF (((a_out .GT. 1_pr) .AND. (x_out .GT. 8e-1_pr)) .OR. (x_out .GT. a_out)) THEN
            a_out = x_out
        ENDIF
        CALL entropy(T_out, vV, entrov)
        CALL entropy(T_out, vL, entrol)
        entro = x_out*entrov + (1.0-qual)*entrol
!    ENDIF    
   ENDIF
! ENDIF
!################## LP ###############  perfect gas OR  span-wagner (iterative)
CASE( 6 )
!        
        IF (u_in .LT. -123.74e3_pr) THEN
! In Solid phase, but we push it back to LP
         u_in  = -123.74e3_pr + 1e3_pr
         f_out = 2_pr
         Print*,'** out of range (SOLID) in Interp_table case(6)'
        ENDIF
! Using directly the EoS of span-wagner for small pressure region
!
!        CALL New_Rap1D(1,temp,out2,res,Niter,&
!                        exitflag,u_in,T_guess,v_in,out2)
!        IF (res > 10e-10) THEN
!          print*, "res IN Interp_table", res, "iter", Niter
!          STOP '** conversion from internal energy to temperature in small pressure  region failed in Interp_table case(6)'
!        ENDIF        
!
!        CALL pressure(temp, v_in, press)
!        CALL sound_speed(temp, v_in, sound)
!
!        T_out = temp
!        p_out = press
!        c_out = sound
!
!        IF ((p_out .GE. 0.5_pr) .OR. (p_out .LT. 0_pr)) THEN
!         STOP '** Problem in small pressure region in Interp_table case (6)' 
!        ENDIF
!--------------  perfect gas formulation--------------------------------------------- 
         gamma_pg = 1.313_pr
         e_pg     = 236.0294e3_pr
         p_out    = (gamma_pg-1_pr)*(u_in+e_pg)/v_in
         c_out    = sqrt(gamma_pg*p_out*v_in)
         T_out    = p_out*v_in/R
         x_out    = 1_pr
         a_out    = 1_pr
         IF ((p_out .GE. 0.5e6_pr) .OR. (p_out .LT. 1e-2_pr) .OR.& 
&            (c_out .LT. 1e-2_pr ) .OR. (T_out .LT. 1e-2_pr)) THEN
         STOP '** Problem Perfect gas formulation in Interp_table case (6)' 
         ENDIF
END SELECT
!
!
!######### CHECK #########
! Pressure
  IF ( (p_out /= p_out) .OR. (p_out < 0.0) ) THEN
     print*, 'LOOK_UP TABLE pressure negative or inf'
     print*, 'Input v= ', v_in, 'e= ',u_in
     STOP
  ENDIF
! Temperature
  IF ( (T_out /= T_out) .OR. (T_out < 0.0) ) THEN
     print*, 'LOOK_UP TABLE Temperature negative or inf'
     print*, 'Input v= ', v_in, 'e= ',u_in
     STOP
  ENDIF 
! Speed of sound
  IF ( (c_out /= c_out) .OR. (c_out < 0.0) ) THEN
     print*, 'LOOK_UP TABLE sound speed negative or inf'
     print*, 'Input v= ', v_in, 'e= ',u_in
     STOP
  ENDIF
! Quality and mass fraction
  IF ( (x_out /= x_out) .OR. (x_out < 0.0) .OR. &
&      (a_out /= a_out) .OR. (a_out < 0.0) ) THEN
     print*, 'LOOK_UP TABLE quanlity and mass fraction  negative or inf'
     print*, 'Input v= ', v_in, 'e= ',u_in
     STOP
  ENDIF
!
END SUBROUTINE CO2TRANS_PROP
!
!
!
END MODULE Transprop
