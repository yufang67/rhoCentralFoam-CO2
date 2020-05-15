extern "C" 
{
//void __properties_MOD_pressure(double* ii,double* jj, double* kk);
void __properties_MOD_inter_energy(double* temp,double* vol, double* energy);
void __grid_MOD_make_grid();
void __interp_table_MOD_co2bllt_equi(double* p_out, double* T_out,double* c_out,
                                     double* x_out, double* a_out,double* dummy,
                                     double* u_in,  double* v_in, int*    flag);
void __solver_eos_MOD_brentroots2(int* MODE, double* out_1,double* out_2,double* residue,
                                  int* Niter,double* GGG,double* lower,double* upper,double* in_1,double* in_2);

void __solver_eos_MOD_eos_1d(int* MODE,double* es,double* out_2,double* residue,int* Niter,int* flag,double* Qp5,double* eguess,double*vbc,double* in_2);
void __derivees_MOD_co2der(double* dp_dv_u,double* dp_du_v ,double* e0,double* v0,double* T0,double* p0,double* res2,double* res3,double* res4);
void __transprop_MOD_co2visco(double* vis,double* vv,double* TT,double* x_out,double* pp,int* flag);
void __transprop_MOD_cpco2(double* cp_out,double* v_in,double* vp_out,double* x_in,double* T_in,double* p_in,int* flag_loca);
void __transprop_MOD_cvco2(double* cv_out,double* v_in,double* vp_out,double* x_in,double* T_in,double* p_in,int* flag_loca);
void __transprop_MOD_dedrco2(double* dedr_out,double* v_in,double* vp_out,double* x_in,double* T_in,double* p_in,int* flag_loca);
void __transprop_MOD_co2conduc2phase(double* lambda_out,double* v_in,double* vp_out,double* x_in,double* T_in,double* p_in,int* flag_loca);
//void __properties_MOD_heat_cap_p(double* temp,double* vol, double* cp)
//void __properties_MOD_heat_cap_v(double* temp,double* vol, double* cv)
//void grid_construction_left_low_();
//void grid_construction_left_high_();
//void grid_construction_high_temperature_();
}
//extern double* __def_variables_MOD_ccc_ll;



