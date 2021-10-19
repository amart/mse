#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <pk13_3.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_report1 = new ofstream("mceval.dat");
  styr.allocate("styr");
  endyr.allocate("endyr");
  base_endyr.allocate("base_endyr");
  hcr_styr.allocate("hcr_styr");
  rcrage.allocate("rcrage");
  trmage.allocate("trmage");
  nbins1.allocate("nbins1");
  nbins2.allocate("nbins2");
  nbins3.allocate("nbins3");
  cattot.allocate(styr,endyr,"cattot");
  cattot_log_sd.allocate(styr,endyr,"cattot_log_sd");
  nyrs_fsh.allocate("nyrs_fsh");
  fshyrs.allocate(1,nyrs_fsh,"fshyrs");
  multN_fsh.allocate(1,nyrs_fsh,"multN_fsh");
  ac_yng_fsh.allocate(1,nyrs_fsh,"ac_yng_fsh");
  ac_old_fsh.allocate(1,nyrs_fsh,"ac_old_fsh");
  nyrslen_fsh.allocate("nyrslen_fsh");
  fshlenyrs.allocate(1,nyrslen_fsh,"fshlenyrs");
  multNlen_fsh.allocate(1,nyrslen_fsh,"multNlen_fsh");
  rwlk_sd.allocate(styr,endyr-1,"rwlk_sd");
  rwlk_sd_short.allocate(styr,base_endyr);
  catp.allocate(1,nyrs_fsh,rcrage,trmage,"catp");
  lenp.allocate(1,nyrslen_fsh,1,nbins1,"lenp");
  wt_fsh.allocate(styr,endyr,rcrage,trmage,"wt_fsh");
  nyrs_srv1_bs.allocate("nyrs_srv1_bs");
  srvyrs1_bs.allocate(1,nyrs_srv1_bs,"srvyrs1_bs");
  indxsurv1_bs.allocate(1,nyrs_srv1_bs,"indxsurv1_bs");
  indxsurv_log_sd1_bs.allocate(1,nyrs_srv1_bs,"indxsurv_log_sd1_bs");
  nyrs_srv1_ek.allocate("nyrs_srv1_ek");
  srvyrs1_ek.allocate(1,nyrs_srv1_ek,"srvyrs1_ek");
  indxsurv1_ek.allocate(1,nyrs_srv1_ek,"indxsurv1_ek");
  indxsurv_log_sd1_ek.allocate(1,nyrs_srv1_ek,"indxsurv_log_sd1_ek");
  nyrs_srv1_dy.allocate("nyrs_srv1_dy");
  srvyrs1_dy.allocate(1,nyrs_srv1_dy,"srvyrs1_dy");
  indxsurv1_dy.allocate(1,nyrs_srv1_dy,"indxsurv1_dy");
  indxsurv_log_sd1_dy.allocate(1,nyrs_srv1_dy,"indxsurv_log_sd1_dy");
  yrfrct_srv1.allocate(styr,endyr,"yrfrct_srv1");
  nyrsac_srv1.allocate("nyrsac_srv1");
  srv_acyrs1.allocate(1,nyrsac_srv1,"srv_acyrs1");
  multN_srv1.allocate(1,nyrsac_srv1,"multN_srv1");
  ac_yng_srv1.allocate(1,nyrsac_srv1,"ac_yng_srv1");
  ac_old_srv1.allocate(1,nyrsac_srv1,"ac_old_srv1");
  nyrslen_srv1.allocate("nyrslen_srv1");
  srv_lenyrs1.allocate(1,nyrslen_srv1,"srv_lenyrs1");
  multNlen_srv1.allocate(1,nyrslen_srv1,"multNlen_srv1");
  srvp1.allocate(1,nyrsac_srv1,rcrage,trmage,"srvp1");
  srvlenp1.allocate(1,nyrslen_srv1,1,nbins3,"srvlenp1");
  wt_srv1.allocate(styr,endyr,rcrage,trmage,"wt_srv1");
  nyrs_srv2.allocate("nyrs_srv2");
  srvyrs2.allocate(1,nyrs_srv2,"srvyrs2");
  indxsurv2.allocate(1,nyrs_srv2,"indxsurv2");
  indxsurv_log_sd2.allocate(1,nyrs_srv2,"indxsurv_log_sd2");
  yrfrct_srv2.allocate(styr,endyr,"yrfrct_srv2");
  nyrsac_srv2.allocate("nyrsac_srv2");
  srv_acyrs2.allocate(1,nyrsac_srv2,"srv_acyrs2");
  multN_srv2.allocate(1,nyrsac_srv2,"multN_srv2");
  ac_yng_srv2.allocate(1,nyrsac_srv2,"ac_yng_srv2");
  ac_old_srv2.allocate(1,nyrsac_srv2,"ac_old_srv2");
  nyrslen_srv2.allocate("nyrslen_srv2");
  srv_lenyrs2.allocate(1,nyrslen_srv2,"srv_lenyrs2");
  multNlen_srv2.allocate(1,nyrslen_srv2,"multNlen_srv2");
  srvp2.allocate(1,nyrsac_srv2,rcrage,trmage,"srvp2");
  srvlenp2.allocate(1,nyrslen_srv2,1,nbins2,"srvlenp2");
  wt_srv2.allocate(styr,endyr,rcrage,trmage,"wt_srv2");
  nyrs_srv3.allocate("nyrs_srv3");
  srvyrs3.allocate(1,nyrs_srv3,"srvyrs3");
  indxsurv3.allocate(1,nyrs_srv3,"indxsurv3");
  indxsurv_log_sd3.allocate(1,nyrs_srv3,"indxsurv_log_sd3");
  yrfrct_srv3.allocate(styr,endyr,"yrfrct_srv3");
  nyrsac_srv3.allocate("nyrsac_srv3");
  srv_acyrs3.allocate(1,nyrsac_srv3,"srv_acyrs3");
  multN_srv3.allocate(1,nyrsac_srv3,"multN_srv3");
  nyrslen_srv3.allocate("nyrslen_srv3");
  srv_lenyrs3.allocate(1,nyrslen_srv3,"srv_lenyrs3");
  multNlen_srv3.allocate(1,nyrslen_srv3,"multNlen_srv3");
  srvp3.allocate(1,nyrsac_srv3,rcrage,trmage,"srvp3");
  srvlenp3.allocate(1,nyrslen_srv3,1,nbins2,"srvlenp3");
  wt_srv3.allocate(styr,endyr,rcrage,trmage,"wt_srv3");
  age_trans.allocate(rcrage,trmage,rcrage,trmage,"age_trans");
  len_trans1.allocate(rcrage,trmage,1,nbins1,"len_trans1");
  len_trans2.allocate(rcrage,trmage,1,nbins2,"len_trans2");
  len_trans3.allocate(rcrage,trmage,1,nbins3,"len_trans3");
  wt_pop.allocate(styr,endyr,rcrage,trmage,"wt_pop");
  wt_spawn.allocate(styr,endyr,rcrage,trmage,"wt_spawn");
  mat_old.allocate(rcrage,trmage,"mat_old");
  mat.allocate(rcrage,trmage,"mat");
  wt_pop_proj.allocate(rcrage,trmage,"wt_pop_proj");
  wt_spawn_proj.allocate(rcrage,trmage,"wt_spawn_proj");
  wt_fsh_proj.allocate(rcrage,trmage,"wt_fsh_proj");
  wt_srv_proj.allocate(rcrage,trmage,"wt_srv_proj");
  Ftarget.allocate(endyr+1,endyr+5,"Ftarget");
  B40.allocate("B40");
  log_mean_recr_proj.allocate("log_mean_recr_proj");
  sigmasq_recr.allocate("sigmasq_recr");
#define MAX_BISECT_ITER  128            // maximum number of iterations for the bisection
#define BISECT_TOL       0.00001        // maximum tolerance for bisection
  if (endyr > base_endyr)
  {
    phase_future = 4;
  }
  else
  {
    phase_future = -4;
  }
  endyr_fsh_dev = base_endyr;
   // mess about with the random walk deviations so that the base_endyr+1 and beyond values don't matter in the likelihood calcs
  for (i = styr; i < (endyr_fsh_dev - 1); i++)
  {
      rwlk_sd_short(i) = rwlk_sd(i);
  }
  rwlk_sd_short(endyr_fsh_dev - 1) = rwlk_sd(endyr-1);
  Tier3_alpha = 0.05;
}

void model_parameters::initializationfunction(void)
{
  mean_log_initN.set_initial_value(0.0);
  mean_log_recruit.set_initial_value(0.0);
  mean_log_F.set_initial_value(-1.6);
  M.set_initial_value(0.30);
  log_q1_bs.set_initial_value(0.0);
  log_q1_ek.set_initial_value(0.0);
  log_q1_dy.set_initial_value(0.3);
  log_q2.set_initial_value(0.0);
  log_q3.set_initial_value(-1.6);
  slp1_fsh_dev.set_initial_value(0.0);
  inf1_fsh_dev.set_initial_value(0.0);
  slp2_fsh_dev.set_initial_value(0.0);
  inf2_fsh_dev.set_initial_value(0.0);
  log_slp1_fsh_mean.set_initial_value(1.0);
  inf1_fsh_mean.set_initial_value(4.0);
  log_slp2_fsh_mean.set_initial_value(1.0);
  inf2_fsh_mean.set_initial_value(8.0);
  log_slp1_fsh_mean_const.set_initial_value(1.0);
  inf1_fsh_mean_const.set_initial_value(4.0);
  log_slp2_fsh_mean_const.set_initial_value(1.0);
  inf2_fsh_mean_const.set_initial_value(8.0);
  log_slp2_srv1.set_initial_value(1.0);
  inf2_srv1.set_initial_value(9.0);
  log_slp1_srv2.set_initial_value(-0.8);
  inf1_srv2.set_initial_value(47.0);
  log_slp2_srv2.set_initial_value(-.14);
  inf2_srv2.set_initial_value(7.0);
  log_slp1_srv3.set_initial_value(0.0);
  inf1_srv3.set_initial_value(5.0);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  M.allocate(0.1,0.5,-1,"M");
  mean_log_initN.allocate(-15,15,-1,"mean_log_initN");
  dev_log_initN.allocate(rcrage+1,trmage,-15,15,-2,"dev_log_initN");
  initN.allocate(rcrage+1,trmage,"initN");
  #ifndef NO_AD_INITIALIZE
    initN.initialize();
  #endif
  mean_log_recruit.allocate(-15,15,1,"mean_log_recruit");
  dev_log_recruit.allocate(styr,endyr,-15,15,4,"dev_log_recruit");
  log_recr_proj.allocate(endyr+1,endyr+5,-5,5,10,"log_recr_proj");
  recruit_proj.allocate(endyr+1,endyr+5,"recruit_proj");
  N_proj.allocate(endyr+1,endyr+5,rcrage,trmage,"N_proj");
  #ifndef NO_AD_INITIALIZE
    N_proj.initialize();
  #endif
  F_proj.allocate(endyr+1,endyr+5,"F_proj");
  #ifndef NO_AD_INITIALIZE
    F_proj.initialize();
  #endif
  Z_proj.allocate(endyr+1,endyr+5,rcrage,trmage,"Z_proj");
  #ifndef NO_AD_INITIALIZE
    Z_proj.initialize();
  #endif
  C_proj.allocate(endyr+1,endyr+5,rcrage,trmage,"C_proj");
  #ifndef NO_AD_INITIALIZE
    C_proj.initialize();
  #endif
  Nsrv_proj.allocate(endyr+1,endyr+5,rcrage,trmage,"Nsrv_proj");
  #ifndef NO_AD_INITIALIZE
    Nsrv_proj.initialize();
  #endif
  slctfsh_proj.allocate(rcrage,trmage,"slctfsh_proj");
  #ifndef NO_AD_INITIALIZE
    slctfsh_proj.initialize();
  #endif
  Ecattot_proj.allocate(endyr+1,endyr+5,"Ecattot_proj");
  #ifndef NO_AD_INITIALIZE
    Ecattot_proj.initialize();
  #endif
  Esumbio_proj.allocate(endyr+1,endyr+5,"Esumbio_proj");
  Espawnbio_proj.allocate(endyr+1,endyr+5,"Espawnbio_proj");
  Esrv_proj.allocate(endyr+1,endyr+5,"Esrv_proj");
  Exrate_proj.allocate(endyr+1,endyr+5,"Exrate_proj");
  sbio.allocate("sbio");
  #ifndef NO_AD_INITIALIZE
  sbio.initialize();
  #endif
  log_slp1_fsh_mean.allocate(-5,5,4,"log_slp1_fsh_mean");
  inf1_fsh_mean.allocate(1,5,4,"inf1_fsh_mean");
  log_slp2_fsh_mean.allocate(-5,5,4,"log_slp2_fsh_mean");
  inf2_fsh_mean.allocate(7,20,4,"inf2_fsh_mean");
  slp1_fsh_dev.allocate(styr,endyr_fsh_dev,-5,5,5,"slp1_fsh_dev");
  inf1_fsh_dev.allocate(styr,endyr_fsh_dev,-5,5,5,"inf1_fsh_dev");
  slp2_fsh_dev.allocate(styr,endyr_fsh_dev,-5,5,5,"slp2_fsh_dev");
  inf2_fsh_dev.allocate(styr,endyr_fsh_dev,-5,5,5,"inf2_fsh_dev");
  slp1_fsh.allocate(styr,endyr_fsh_dev,"slp1_fsh");
  #ifndef NO_AD_INITIALIZE
    slp1_fsh.initialize();
  #endif
  inf1_fsh.allocate(styr,endyr_fsh_dev,"inf1_fsh");
  #ifndef NO_AD_INITIALIZE
    inf1_fsh.initialize();
  #endif
  slp2_fsh.allocate(styr,endyr_fsh_dev,"slp2_fsh");
  #ifndef NO_AD_INITIALIZE
    slp2_fsh.initialize();
  #endif
  inf2_fsh.allocate(styr,endyr_fsh_dev,"inf2_fsh");
  #ifndef NO_AD_INITIALIZE
    inf2_fsh.initialize();
  #endif
  log_slp1_fsh_mean_const.allocate(-5,5,phase_future,"log_slp1_fsh_mean_const");
  inf1_fsh_mean_const.allocate(1,5,phase_future,"inf1_fsh_mean_const");
  log_slp2_fsh_mean_const.allocate(-5,5,phase_future,"log_slp2_fsh_mean_const");
  inf2_fsh_mean_const.allocate(7,20,phase_future,"inf2_fsh_mean_const");
  slp1_fsh_const.allocate("slp1_fsh_const");
  #ifndef NO_AD_INITIALIZE
  slp1_fsh_const.initialize();
  #endif
  inf1_fsh_const.allocate("inf1_fsh_const");
  #ifndef NO_AD_INITIALIZE
  inf1_fsh_const.initialize();
  #endif
  slp2_fsh_const.allocate("slp2_fsh_const");
  #ifndef NO_AD_INITIALIZE
  slp2_fsh_const.initialize();
  #endif
  inf2_fsh_const.allocate("inf2_fsh_const");
  #ifndef NO_AD_INITIALIZE
  inf2_fsh_const.initialize();
  #endif
  log_slp2_srv1.allocate(-5,5,7,"log_slp2_srv1");
  inf2_srv1.allocate(5,20,7,"inf2_srv1");
  srv1_age1.allocate(0,2,7,"srv1_age1");
  log_slp1_srv2.allocate(-5,5,7,"log_slp1_srv2");
  inf1_srv2.allocate(1,50,-1,"inf1_srv2");
  log_slp2_srv2.allocate(-5,5,8,"log_slp2_srv2");
  inf2_srv2.allocate(5,10,8,"inf2_srv2");
  srv2_age1.allocate(0,2,7,"srv2_age1");
  log_slp1_srv3.allocate(-5,5,9,"log_slp1_srv3");
  inf1_srv3.allocate(1,20,9,"inf1_srv3");
  mean_log_F.allocate(-10,10,1,"mean_log_F");
  dev_log_F.allocate(styr,endyr,-10,10,2,"dev_log_F");
  F.allocate(styr,endyr,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  log_q1_bs.allocate(-10,10,5,"log_q1_bs");
  log_q1_ek.allocate(-10,10,5,"log_q1_ek");
  log_q1_dy.allocate(-10,10,5,"log_q1_dy");
  log_q2.allocate(-10,10,-1,"log_q2");
  log_q3.allocate(-10,10,6,"log_q3");
  q1_bs.allocate("q1_bs");
  #ifndef NO_AD_INITIALIZE
  q1_bs.initialize();
  #endif
  q1_ek.allocate("q1_ek");
  #ifndef NO_AD_INITIALIZE
  q1_ek.initialize();
  #endif
  q1_dy.allocate("q1_dy");
  #ifndef NO_AD_INITIALIZE
  q1_dy.initialize();
  #endif
  q2.allocate("q2");
  #ifndef NO_AD_INITIALIZE
  q2.initialize();
  #endif
  q3.allocate("q3");
  #ifndef NO_AD_INITIALIZE
  q3.initialize();
  #endif
  avgR.allocate("avgR");
  #ifndef NO_AD_INITIALIZE
  avgR.initialize();
  #endif
  avgR_CV.allocate("avgR_CV");
  #ifndef NO_AD_INITIALIZE
  avgR_CV.initialize();
  #endif
  SB0.allocate("SB0");
  #ifndef NO_AD_INITIALIZE
  SB0.initialize();
  #endif
  R0.allocate("R0");
  #ifndef NO_AD_INITIALIZE
  R0.initialize();
  #endif
  phi0.allocate("phi0");
  #ifndef NO_AD_INITIALIZE
  phi0.initialize();
  #endif
  SBcurr.allocate("SBcurr");
  #ifndef NO_AD_INITIALIZE
  SBcurr.initialize();
  #endif
  F100.allocate("F100");
  #ifndef NO_AD_INITIALIZE
  F100.initialize();
  #endif
  SB100.allocate("SB100");
  #ifndef NO_AD_INITIALIZE
  SB100.initialize();
  #endif
  SBtarget.allocate("SBtarget");
  #ifndef NO_AD_INITIALIZE
  SBtarget.initialize();
  #endif
  F40.allocate("F40");
  #ifndef NO_AD_INITIALIZE
  F40.initialize();
  #endif
  SB40.allocate("SB40");
  #ifndef NO_AD_INITIALIZE
  SB40.initialize();
  #endif
  F35.allocate("F35");
  #ifndef NO_AD_INITIALIZE
  F35.initialize();
  #endif
  SB35.allocate("SB35");
  #ifndef NO_AD_INITIALIZE
  SB35.initialize();
  #endif
  F20.allocate("F20");
  #ifndef NO_AD_INITIALIZE
  F20.initialize();
  #endif
  SB20.allocate("SB20");
  #ifndef NO_AD_INITIALIZE
  SB20.initialize();
  #endif
  F_ABC.allocate("F_ABC");
  #ifndef NO_AD_INITIALIZE
  F_ABC.initialize();
  #endif
  ABC.allocate("ABC");
  #ifndef NO_AD_INITIALIZE
  ABC.initialize();
  #endif
  F_OFL.allocate("F_OFL");
  #ifndef NO_AD_INITIALIZE
  F_OFL.initialize();
  #endif
  OFL.allocate("OFL");
  #ifndef NO_AD_INITIALIZE
  OFL.initialize();
  #endif
  N.allocate(styr,endyr+5,rcrage,trmage,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  endN.allocate(rcrage,trmage,"endN");
  all_recruits.allocate(styr,endyr+5,"all_recruits");
  #ifndef NO_AD_INITIALIZE
    all_recruits.initialize();
  #endif
  Z.allocate(styr,endyr,rcrage,trmage,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  C.allocate(styr,endyr,rcrage,trmage,"C");
  #ifndef NO_AD_INITIALIZE
    C.initialize();
  #endif
  Nsrv1.allocate(styr,endyr,rcrage,trmage,"Nsrv1");
  #ifndef NO_AD_INITIALIZE
    Nsrv1.initialize();
  #endif
  slctsrv1.allocate(rcrage,trmage,"slctsrv1");
  #ifndef NO_AD_INITIALIZE
    slctsrv1.initialize();
  #endif
  Nsrv2.allocate(styr,endyr,rcrage,trmage,"Nsrv2");
  #ifndef NO_AD_INITIALIZE
    Nsrv2.initialize();
  #endif
  slctsrv2.allocate(rcrage,trmage,"slctsrv2");
  #ifndef NO_AD_INITIALIZE
    slctsrv2.initialize();
  #endif
  Nsrv3.allocate(styr,endyr,rcrage,trmage,"Nsrv3");
  #ifndef NO_AD_INITIALIZE
    Nsrv3.initialize();
  #endif
  slctsrv3.allocate(rcrage,trmage,"slctsrv3");
  #ifndef NO_AD_INITIALIZE
    slctsrv3.initialize();
  #endif
  slctfsh_base.allocate(rcrage,trmage,"slctfsh_base");
  #ifndef NO_AD_INITIALIZE
    slctfsh_base.initialize();
  #endif
  slctfsh.allocate(styr,endyr,rcrage,trmage,"slctfsh");
  #ifndef NO_AD_INITIALIZE
    slctfsh.initialize();
  #endif
  slctfsh_const.allocate(rcrage,trmage,"slctfsh_const");
  #ifndef NO_AD_INITIALIZE
    slctfsh_const.initialize();
  #endif
  Eecocon.allocate(styr,endyr,"Eecocon");
  #ifndef NO_AD_INITIALIZE
    Eecocon.initialize();
  #endif
  Eec.allocate(styr,endyr,rcrage,trmage,"Eec");
  #ifndef NO_AD_INITIALIZE
    Eec.initialize();
  #endif
  Ecattot.allocate(styr,endyr,"Ecattot");
  #ifndef NO_AD_INITIALIZE
    Ecattot.initialize();
  #endif
  Ecatp.allocate(styr,endyr,rcrage,trmage,"Ecatp");
  #ifndef NO_AD_INITIALIZE
    Ecatp.initialize();
  #endif
  Elenp.allocate(styr,endyr,1,nbins1,"Elenp");
  #ifndef NO_AD_INITIALIZE
    Elenp.initialize();
  #endif
  Eindxsurv1_bs.allocate(styr,endyr,"Eindxsurv1_bs");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv1_bs.initialize();
  #endif
  Eindxsurv1_ek.allocate(styr,endyr,"Eindxsurv1_ek");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv1_ek.initialize();
  #endif
  Eindxsurv1_dy.allocate(styr,endyr,"Eindxsurv1_dy");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv1_dy.initialize();
  #endif
  Esrvp1.allocate(styr,endyr,rcrage,trmage,"Esrvp1");
  #ifndef NO_AD_INITIALIZE
    Esrvp1.initialize();
  #endif
  Esrvlenp1.allocate(styr,endyr,1,nbins3,"Esrvlenp1");
  #ifndef NO_AD_INITIALIZE
    Esrvlenp1.initialize();
  #endif
  Eindxsurv2.allocate(styr,endyr,"Eindxsurv2");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv2.initialize();
  #endif
  Esrvp2.allocate(styr,endyr,rcrage,trmage,"Esrvp2");
  #ifndef NO_AD_INITIALIZE
    Esrvp2.initialize();
  #endif
  Esrvlenp2.allocate(styr,endyr,1,nbins2,"Esrvlenp2");
  #ifndef NO_AD_INITIALIZE
    Esrvlenp2.initialize();
  #endif
  Eindxsurv3.allocate(styr,endyr,"Eindxsurv3");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv3.initialize();
  #endif
  Esrvp3.allocate(styr,endyr,rcrage,trmage,"Esrvp3");
  #ifndef NO_AD_INITIALIZE
    Esrvp3.initialize();
  #endif
  Esrvlenp3.allocate(styr,endyr,1,nbins2,"Esrvlenp3");
  #ifndef NO_AD_INITIALIZE
    Esrvlenp3.initialize();
  #endif
  loglik.allocate(1,22,"loglik");
  #ifndef NO_AD_INITIALIZE
    loglik.initialize();
  #endif
  llcatp.allocate(1,nyrs_fsh,"llcatp");
  #ifndef NO_AD_INITIALIZE
    llcatp.initialize();
  #endif
  lllenp.allocate(1,nyrslen_fsh,"lllenp");
  #ifndef NO_AD_INITIALIZE
    lllenp.initialize();
  #endif
  llsrvp1.allocate(1,nyrsac_srv1,"llsrvp1");
  #ifndef NO_AD_INITIALIZE
    llsrvp1.initialize();
  #endif
  llsrvlenp1.allocate(1,nyrslen_srv1,"llsrvlenp1");
  #ifndef NO_AD_INITIALIZE
    llsrvlenp1.initialize();
  #endif
  llsrvp2.allocate(1,nyrsac_srv2,"llsrvp2");
  #ifndef NO_AD_INITIALIZE
    llsrvp2.initialize();
  #endif
  llsrvlenp2.allocate(1,nyrslen_srv2,"llsrvlenp2");
  #ifndef NO_AD_INITIALIZE
    llsrvlenp2.initialize();
  #endif
  llsrvp3.allocate(1,nyrsac_srv3,"llsrvp3");
  #ifndef NO_AD_INITIALIZE
    llsrvp3.initialize();
  #endif
  llsrvlenp3.allocate(1,nyrslen_srv3,"llsrvlenp3");
  #ifndef NO_AD_INITIALIZE
    llsrvlenp3.initialize();
  #endif
  recruit.allocate(styr,endyr,"recruit");
  sd_avg_rec.allocate("sd_avg_rec");
  Espawnbio.allocate(styr,endyr+1,"Espawnbio");
  Esumbio.allocate(styr,endyr+1,"Esumbio");
  res_fish.allocate(1,nyrs_fsh,rcrage,2*trmage-rcrage+1,"res_fish");
  #ifndef NO_AD_INITIALIZE
    res_fish.initialize();
  #endif
  res_srv1.allocate(1,nyrsac_srv1,rcrage,2*trmage-rcrage+1,"res_srv1");
  #ifndef NO_AD_INITIALIZE
    res_srv1.initialize();
  #endif
  res_srv2.allocate(1,nyrsac_srv2,rcrage,2*trmage-rcrage+1,"res_srv2");
  #ifndef NO_AD_INITIALIZE
    res_srv2.initialize();
  #endif
  res_srv3.allocate(1,nyrsac_srv3,rcrage,2*trmage-rcrage+1,"res_srv3");
  #ifndef NO_AD_INITIALIZE
    res_srv3.initialize();
  #endif
  res_srv3len.allocate(1,nyrslen_srv3,1,2*nbins2,"res_srv3len");
  #ifndef NO_AD_INITIALIZE
    res_srv3len.initialize();
  #endif
  pearson_fish.allocate(1,nyrs_fsh,rcrage,trmage,"pearson_fish");
  #ifndef NO_AD_INITIALIZE
    pearson_fish.initialize();
  #endif
  pearson_srv1.allocate(1,nyrsac_srv1,rcrage,trmage,"pearson_srv1");
  #ifndef NO_AD_INITIALIZE
    pearson_srv1.initialize();
  #endif
  pearson_srv2.allocate(1,nyrsac_srv2,rcrage,trmage,"pearson_srv2");
  #ifndef NO_AD_INITIALIZE
    pearson_srv2.initialize();
  #endif
  pearson_srv3.allocate(1,nyrsac_srv3,rcrage,trmage,"pearson_srv3");
  #ifndef NO_AD_INITIALIZE
    pearson_srv3.initialize();
  #endif
  pearson_srv3len.allocate(1,nyrslen_srv3,1,nbins2,"pearson_srv3len");
  #ifndef NO_AD_INITIALIZE
    pearson_srv3len.initialize();
  #endif
  effN_fsh.allocate(1,nyrs_fsh,"effN_fsh");
  #ifndef NO_AD_INITIALIZE
    effN_fsh.initialize();
  #endif
  effN_srv1.allocate(1,nyrsac_srv1,"effN_srv1");
  #ifndef NO_AD_INITIALIZE
    effN_srv1.initialize();
  #endif
  effN_srv2.allocate(1,nyrsac_srv2,"effN_srv2");
  #ifndef NO_AD_INITIALIZE
    effN_srv2.initialize();
  #endif
  effN_srv3.allocate(1,nyrsac_srv3,"effN_srv3");
  #ifndef NO_AD_INITIALIZE
    effN_srv3.initialize();
  #endif
  RMSE_srv1_bs.allocate("RMSE_srv1_bs");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv1_bs.initialize();
  #endif
  RMSE_srv1_ek.allocate("RMSE_srv1_ek");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv1_ek.initialize();
  #endif
  RMSE_srv1_dy.allocate("RMSE_srv1_dy");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv1_dy.initialize();
  #endif
  RMSE_srv2.allocate("RMSE_srv2");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv2.initialize();
  #endif
  RMSE_srv3.allocate("RMSE_srv3");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv3.initialize();
  #endif
  var_prof.allocate("var_prof");
  objfun.allocate("objfun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  for (i=1;i<=nyrs_fsh;i++)
    {
  for (j=rcrage;j<=trmage;j++)
    {
    if(j<ac_yng_fsh(i))
      {
      catp(i,ac_yng_fsh(i)) += catp(i,j);
      catp(i,j) = 0;
      }
    if(j>ac_old_fsh(i))
      {
      catp(i,ac_old_fsh(i)) += catp(i,j);
      catp(i,j) = 0;
      }
    }}
  for (i=1;i<=nyrsac_srv1;i++)
    {
  for (j=rcrage;j<=trmage;j++)
    {
    if(j<ac_yng_srv1(i))
      {
      srvp1(i,ac_yng_srv1(i)) += srvp1(i,j);
      srvp1(i,j) = 0;
      }
    if(j>ac_old_srv1(i))
      {
      srvp1(i,ac_old_srv1(i)) += srvp1(i,j);
      srvp1(i,j) = 0;
      }
    }}
  for (i=1;i<=nyrsac_srv2;i++)
    {
  for (j=rcrage;j<=trmage;j++)
    {
    if(j<ac_yng_srv2(i))
      {
      srvp2(i,ac_yng_srv2(i)) += srvp2(i,j);
      srvp2(i,j) = 0;
      }
    if(j>ac_old_srv2(i))
      {
      srvp2(i,ac_old_srv2(i)) += srvp2(i,j);
      srvp2(i,j) = 0;
      }
    }}
  mcmc_iter = 0;
  y_frac_sp = 0.21;             // this matches what MWD uses; rounded value for 15 March
  o = 0.00001;
  var_prof.set_stepnumber(30);
  var_prof.set_stepsize(0.1);
}

void model_parameters::userfunction(void)
{
  objfun =0.0;
  ofstream& report1= *pad_report1;
  Convert_log_parameters();
  Selectivity();
  Mortality();
  Numbers_at_age();
  Catch_at_age();
  Expected_values();
  if(last_phase())
  {
  Projections();
  }
  Objective_function();
  MCMC_output();
}

void model_parameters::Convert_log_parameters(void)
{
  ofstream& report1= *pad_report1;
  for (j=rcrage+1;j<=trmage;j++)
  {
  initN(j) = mfexp(mean_log_recruit +  dev_log_recruit(styr) - M*double(j-rcrage)+dev_log_initN(j));
  }
  initN(trmage) /= (1.0 - mfexp(-M));
  recruit = mfexp(mean_log_recruit +  dev_log_recruit);
  F = mfexp(mean_log_F + dev_log_F);
  q1_bs = mfexp(log_q1_bs);
  q1_ek = mfexp(log_q1_ek);
  q1_dy = mfexp(log_q1_dy);
  q2 = mfexp(log_q2);
  q3 = mfexp(log_q3);
}

void model_parameters::Selectivity(void)
{
  ofstream& report1= *pad_report1;
   for (i=styr;i<=endyr_fsh_dev;i++)
   {
   slp1_fsh(i)=mfexp(log_slp1_fsh_mean+slp1_fsh_dev(i));
   inf1_fsh(i)=inf1_fsh_mean+inf1_fsh_dev(i);
   slp2_fsh(i)=mfexp(log_slp2_fsh_mean+slp2_fsh_dev(i));
   inf2_fsh(i)=inf2_fsh_mean+inf2_fsh_dev(i);
   for (j=rcrage;j<=trmage;j++)
   {
   slctfsh(i,j) = (1/(1+mfexp(-(slp1_fsh(i))*(double(j)-(inf1_fsh(i))))))*
                  (1-1/(1+mfexp(-(slp2_fsh(i))*(double(j)-(inf2_fsh(i))))));
   }
   slctfsh(i)=slctfsh(i)/slctfsh(i,6);
   }
   // ZTA - calculate constant future fishery selectivity
   if (endyr > base_endyr)
   {
       slctfsh_const.initialize();
       slp1_fsh_const=mfexp(log_slp1_fsh_mean_const);
       inf1_fsh_const=inf1_fsh_mean_const;
       slp2_fsh_const=mfexp(log_slp2_fsh_mean_const);
       inf2_fsh_const=inf2_fsh_mean_const;
       for (j = rcrage; j <= trmage; j++)
       {
           slctfsh_const(j) = (1.0/(1.0+mfexp(-(slp1_fsh_const)*(double(j)-(inf1_fsh_const)))))*
                              (1.0-(1.0/(1.0+mfexp(-(slp2_fsh_const)*(double(j)-(inf2_fsh_const))))));
       }
       slctfsh_const /= max(slctfsh_const);
       for (i = (base_endyr+1); i <= endyr; i++)
       {
           slctfsh(i) = slctfsh_const;
       }
   }
   for (j=rcrage;j<=trmage;j++)
   {
       slctfsh_base(j) = (1.0/(1.0+mfexp(-(mfexp(log_slp1_fsh_mean))*(double(j)-(inf1_fsh_mean)))))*
                         (1.0-1.0/(1.0+mfexp(-(mfexp(log_slp2_fsh_mean))*(double(j)-(inf2_fsh_mean)))));
   }
   slctfsh_base=slctfsh_base/max(slctfsh_base);
  for (j=rcrage;j<=trmage;j++)
    {
    slctsrv1(j) = (1-1/(1+mfexp(-mfexp(log_slp2_srv1)*(double(j)-inf2_srv1))));
    }
    slctsrv1=slctsrv1/slctsrv1(2);
	slctsrv1(rcrage)=srv1_age1;
  for (j=rcrage;j<=trmage;j++)
    {
    slctsrv2(j) = (1/(1+mfexp(-mfexp(log_slp1_srv2)*(double(j)-inf1_srv2))))
  *(1-1/(1+mfexp(-mfexp(log_slp2_srv2)*(double(j)-inf2_srv2))));
    }
    slctsrv2=slctsrv2/slctsrv2(7);
	slctsrv2(rcrage)=srv2_age1;
  for (j=rcrage;j<=trmage;j++)
    {
    slctsrv3(j) = (1/(1+mfexp(-mfexp(log_slp1_srv3)*(double(j)-inf1_srv3))));
    }
    slctsrv3=slctsrv3/slctsrv3(10);
}

void model_parameters::Mortality(void)
{
  ofstream& report1= *pad_report1;
  for (i=styr;i<=endyr;i++)
    {
  for (j=rcrage;j<=trmage;j++)
    {
    Z(i,j)=(F(i)*slctfsh(i,j))+M;
    }}
}

void model_parameters::Numbers_at_age(void)
{
  ofstream& report1= *pad_report1;
  N.initialize();
  N(styr)(rcrage+1,trmage)=initN;
  for (i=styr;i<=endyr;i++)
    {
    N(i,rcrage)=recruit(i);
    }
  for (i=styr;i<endyr;i++)
    {
  for (j=rcrage;j<trmage;j++)
    {
    N(i+1,j+1)=N(i,j)*mfexp(-Z(i,j));
    }
    N(i+1,trmage)+=N(i,trmage)*mfexp(-Z(i,trmage));
    }
  endN=N(endyr);
  sd_avg_rec = mean(recruit((hcr_styr+rcrage),(endyr-1)));
}

void model_parameters::Catch_at_age(void)
{
  ofstream& report1= *pad_report1;
  for (i=styr;i<=endyr;i++)
    {
  for (j=rcrage;j<=trmage;j++)
    {
    C(i,j)=N(i,j)*((F(i)*slctfsh(i,j))/Z(i,j))*(1-mfexp(-Z(i,j)));
    Eec(i,j)=N(i,j)*(M/Z(i,j))*(1-mfexp(-Z(i,j)));
    Nsrv1(i,j)=slctsrv1(j)*N(i,j)*mfexp(-yrfrct_srv1(i)*Z(i,j));
    Nsrv2(i,j)=slctsrv2(j)*N(i,j)*mfexp(-yrfrct_srv2(i)*Z(i,j));
    Nsrv3(i,j)=slctsrv3(j)*N(i,j)*mfexp(-yrfrct_srv3(i)*Z(i,j));
     }}
}

void model_parameters::Expected_values(void)
{
  ofstream& report1= *pad_report1;
  for (i=styr;i<=endyr;i++)
    {
    Ecattot(i) = 1000000*sum(elem_prod(C(i),wt_fsh(i)));
    Eecocon(i) = 1000000*sum(elem_prod(Eec(i),wt_pop(i)));
    Ecatp(i) = (C(i)/sum(C(i)))*age_trans;
    Elenp(i) = Ecatp(i) * len_trans1;
    Eindxsurv1_bs(i)= q1_bs*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv1(i)*Z(i))),slctsrv1),wt_srv1(i)));
    Eindxsurv1_ek(i)= q1_ek*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv1(i)*Z(i))),slctsrv1),wt_srv1(i)));
    Eindxsurv1_dy(i)= q1_dy*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv1(i)*Z(i))),slctsrv1),wt_srv1(i)));
    Esrvp1(i) = (Nsrv1(i)/sum(Nsrv1(i)))*age_trans;
    Esrvlenp1(i) = Esrvp1(i) * len_trans3;
    Eindxsurv2(i)= q2*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv2(i)*Z(i))),slctsrv2),wt_srv2(i)));
    Esrvp2(i) = (Nsrv2(i)/sum(Nsrv2(i)))*age_trans;
    Esrvlenp2(i) = Esrvp2(i) * len_trans2;
    Eindxsurv3(i)= q3*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv3(i)*Z(i))),slctsrv3),wt_srv3(i)));
    Esrvp3(i) = (Nsrv3(i)/sum(Nsrv3(i)))*age_trans;
    Esrvlenp3(i) = Esrvp3(i) * len_trans2;
    Esumbio(i)= N(i)(rcrage+2,trmage)*wt_pop(i)(rcrage+2,trmage);
    Espawnbio(i)= sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-0.21*Z(i))),wt_spawn(i)),0.5*mat));
    }
  for (i=1;i<=nyrs_fsh;i++)
    {
  for (j=rcrage;j<=trmage;j++)
    {
    if(j<ac_yng_fsh(i))
      {
      Ecatp(fshyrs(i),ac_yng_fsh(i)) += Ecatp(fshyrs(i),j);
      Ecatp(fshyrs(i),j) = 0;
      }
    if(j>ac_old_fsh(i))
      {
      Ecatp(fshyrs(i),ac_old_fsh(i)) += Ecatp(fshyrs(i),j);
      Ecatp(fshyrs(i),j) = 0;
      }
    }}
  for (i=1;i<=nyrsac_srv1;i++)
    {
  for (j=rcrage;j<=trmage;j++)
    {
    if(j<ac_yng_srv1(i))
      {
      Esrvp1(srv_acyrs1(i),ac_yng_srv1(i)) += Esrvp1(srv_acyrs1(i),j);
      Esrvp1(srv_acyrs1(i),j) = 0;
      }
    if(j>ac_old_srv1(i))
      {
      Esrvp1(srv_acyrs1(i),ac_old_srv1(i)) += Esrvp1(srv_acyrs1(i),j);
      Esrvp1(srv_acyrs1(i),j) = 0;
      }
    }}
  for (i=1;i<=nyrsac_srv2;i++)
    {
  for (j=rcrage;j<=trmage;j++)
    {
    if(j<ac_yng_srv2(i))
      {
      Esrvp2(srv_acyrs2(i),ac_yng_srv2(i)) += Esrvp2(srv_acyrs2(i),j);
      Esrvp2(srv_acyrs2(i),j) = 0.;
      }
    if(j>ac_old_srv2(i))
      {
      Esrvp2(srv_acyrs2(i),ac_old_srv2(i)) += Esrvp2(srv_acyrs2(i),j);
      Esrvp2(srv_acyrs2(i),j) = 0;
      }
    }}
}

void model_parameters::Projections(void)
{
  ofstream& report1= *pad_report1;
  dvariable sbio, sumbio;
  sbio = sumbio = 0.0;
  all_recruits.initialize();
  all_recruits = column(N,rcrage);
 for (i=endyr+1;i<=endyr+5;i++)
    {
    recruit_proj(i)=mfexp(log_recr_proj(i)+(sigmasq_recr/2));
    all_recruits(i) = mean(recruit((hcr_styr+rcrage),(endyr-1)));
    }
 for (i=endyr+1;i<=endyr+5;i++)
    {
    N_proj(i,rcrage)=recruit_proj(i);
    }
  for (j=rcrage;j<trmage;j++)
    {
    N_proj(endyr+1,j+1)=N(endyr,j)*mfexp(-Z(endyr,j));
    }
    N_proj(endyr+1,trmage)+=N(endyr,trmage)*mfexp(-Z(endyr,trmage));
  endyr_avg_slct=endyr-1;
  styr_avg_slct=endyr_avg_slct-4;
  for (j=rcrage;j<=trmage;j++)
    {
    slctfsh_proj(j) = 0;
   for (i=styr_avg_slct;i<=endyr_avg_slct;i++)
    {
    slctfsh_proj(j) += slctfsh(i,j);
    }
    }
   slctfsh_proj=slctfsh_proj/max(slctfsh_proj);
  for (i=endyr+1;i<=endyr+5;i++)
    {
    // ZTA - assuming that N() are the numbers-at-age at the beginning of the year,
    // ZTA - this updates N() for use in calculate_curr_bio_ref_points()
    N(i) = N_proj(i);
    // ZTA - calculate values for decision rule
    F100  = get_spr_rates(1.00, slctfsh_proj, (i-1));
    SB100 = SBcurr;
    F40   = get_spr_rates(0.40, slctfsh_proj, (i-1));
    SB40  = SBcurr;
    F35   = get_spr_rates(0.35, slctfsh_proj, (i-1));
    SB35  = SBcurr;
    F20   = get_spr_rates(0.20, slctfsh_proj, (i-1));
    SB20  = SBcurr;
    SBtarget = SB40 * (F35 / F40);
    // ZTA - initialize
    F_proj(i) = F40;
    // cout << "in Projections:\tyear " << i << "\tF40 " << F40 << "\tSBtarget " << SBtarget << endl;
   for (loop=1;loop<=10;loop++)
    {
    for (j=rcrage;j<=trmage;j++)
    {
    Z_proj(i,j)=(F_proj(i)*slctfsh_proj(j))+M;
    }
    sbio = sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-0.21*Z_proj(i))),wt_spawn_proj),0.5*mat));
  // MWD old version
    // F_proj(i)=Ftarget(i);
    // if (sbio < B40)
    // {
    // F_proj(i)=Ftarget(i)*(((sbio/B40)-0.05)/(1-0.05));
    // }
    // }
    if (sbio < (0.25 * SB20))
    {
        F_proj(i) = 0.0;
    }
    else if (sbio < SB20)
    {
        // F_proj(i) = 0.0;
        // one one-thousandth of the current biomass caught as bycatch in other fisheries
        sumbio = sum(elem_prod(N_proj(i)(rcrage+2,trmage),wt_pop_proj(rcrage+2,trmage)));
        F_proj(i) = solve_for_fishing_mortality((0.001*1000000.0*sumbio),i);
    }
    else if (sbio < SBtarget)
    {
        F_proj(i)=F40*(((sbio/SBtarget)-0.05)/(1.0-0.05));
    }
    else
    {
        F_proj(i)=F40;
    }
    }
  for (j=rcrage;j<=trmage;j++)
    {
    Z_proj(i,j)=(F_proj(i)*slctfsh_proj(j))+M;
    }
    sbio = sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-0.21*Z_proj(i))),wt_spawn_proj),0.5*mat));
    if (sbio < SB20)
    {
        sumbio = sum(elem_prod(N_proj(i)(rcrage+2,trmage),wt_pop_proj(rcrage+2,trmage)));
        cout << "in pk13_3::Projections:\tyear " << i << ", B40 " << B40 << ", SB40 " << SB40  << ", SB20 " << SB20 << ", Low SB " << sbio << ", sumbio " << sumbio << ", F_ABC " << F_proj(i) << endl;
    }
  if(i<endyr+5)
  {
  for (j=rcrage;j<trmage;j++)
   {
   N_proj(i+1,j+1)=N_proj(i,j)*mfexp(-Z_proj(i,j));
   }
   N_proj(i+1,trmage)+=N_proj(i,trmage)*mfexp(-Z_proj(i,trmage));
   N(i+1) = N_proj(i+1);
  }
  for (j=rcrage;j<=trmage;j++)
    {
    C_proj(i,j)=N_proj(i,j)*((F_proj(i)*slctfsh_proj(j))/Z_proj(i,j))*(1-mfexp(-Z_proj(i,j)));
    Nsrv_proj(i,j)=N_proj(i,j)*mfexp(-yrfrct_srv1(endyr)*Z_proj(i,j));
    }
    Ecattot_proj(i) = 1000000*sum(elem_prod(C_proj(i),wt_fsh_proj));
    Esumbio_proj(i)= N_proj(i)(rcrage+2,trmage)*wt_pop_proj(rcrage+2,trmage);
    Exrate_proj(i)=Ecattot_proj(i)/(1000000*Esumbio_proj(i));
    Espawnbio_proj(i)= sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-0.21*Z_proj(i))),wt_spawn_proj),0.5*mat));
    Esrv_proj(i)= q1_ek*sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-yrfrct_srv1(endyr)*Z_proj(i))),slctsrv1),wt_srv_proj));
    }
}

void model_parameters::Objective_function(void)
{
  ofstream& report1= *pad_report1;
  loglik(1) = -.5*norm2(elem_div((log(cattot)-log(Ecattot)),cattot_log_sd));
  for (i=1;i<=nyrs_fsh;i++)
    {
    llcatp(i) = 0;
  for (j=ac_yng_fsh(i);j<=ac_old_fsh(i);j++)
    {
      llcatp(i) += multN_fsh(i)*(catp(i,j)+o)*log((Ecatp(fshyrs(i),j)+o)/(catp(i,j)+o));
    res_fish(i,j)=catp(i,j);
    res_fish(i,trmage-rcrage+j+1)=Ecatp(fshyrs(i),j);
  if(multN_fsh(i)>0)
    {
	// pearson_fish(i,j)=(catp(i,j)-Ecatp(fshyrs(i),j))/sqrt((Ecatp(fshyrs(i),j)*(1.-Ecatp(fshyrs(i),j)))/multN_fsh(i));
    }
	}
  if(multN_fsh(i)>0)
    {
	effN_fsh(i) = sum(elem_prod(Ecatp(fshyrs(i)),(1-Ecatp(fshyrs(i)))))/sum(square(catp(i)-Ecatp(fshyrs(i))));
	}
	}
  loglik(2) = sum(llcatp);
  for (i=1;i<=nyrslen_fsh;i++)
    {
    lllenp(i) = 0;
  for (j=1;j<=nbins1;j++)
    {
      lllenp(i) += multNlen_fsh(i)*(lenp(i,j)+o)*log((Elenp(fshlenyrs(i),j)+o)/(lenp(i,j)+o));
    }}
  loglik(3) = sum(lllenp);
  loglik(4) = -.5*norm2(elem_div(
       (log(indxsurv1_bs)-log(Eindxsurv1_bs(srvyrs1_bs))+square(indxsurv_log_sd1_bs)/2.),indxsurv_log_sd1_bs));
  loglik(4) += -.5*norm2(elem_div(
       (log(indxsurv1_ek)-log(Eindxsurv1_ek(srvyrs1_ek))+square(indxsurv_log_sd1_ek)/2.),indxsurv_log_sd1_ek));
  loglik(4) += -.5*norm2(elem_div(
       (log(indxsurv1_dy)-log(Eindxsurv1_dy(srvyrs1_dy))+square(indxsurv_log_sd1_dy)/2.),indxsurv_log_sd1_dy));
  RMSE_srv1_bs= sqrt(norm2(log(indxsurv1_bs)-log(Eindxsurv1_bs(srvyrs1_bs))+square(indxsurv_log_sd1_bs)/2.)/nyrs_srv1_bs);
  RMSE_srv1_ek= sqrt(norm2(log(indxsurv1_ek)-log(Eindxsurv1_ek(srvyrs1_ek))+square(indxsurv_log_sd1_ek)/2.)/nyrs_srv1_ek);
  RMSE_srv1_dy= sqrt(norm2(log(indxsurv1_dy)-log(Eindxsurv1_dy(srvyrs1_dy))+square(indxsurv_log_sd1_dy)/2.)/nyrs_srv1_dy);
  for (i=1;i<=nyrsac_srv1;i++)
    {
    llsrvp1(i) = 0;
	for (j=ac_yng_srv1(i);j<=ac_old_srv1(i);j++)
    {
      llsrvp1(i) += multN_srv1(i)*(srvp1(i,j)+o)*log((Esrvp1(srv_acyrs1(i),j)+o)/(srvp1(i,j)+o));
	res_srv1(i,j)=srvp1(i,j);
    res_srv1(i,trmage-rcrage+j+1)=Esrvp1(srv_acyrs1(i),j);
	if(multN_srv1(i)>0)
    {
	// pearson_srv1(i,j)=(srvp1(i,j)-Esrvp1(srv_acyrs1(i),j))/sqrt((Esrvp1(srv_acyrs1(i),j)*(1.-Esrvp1(srv_acyrs1(i),j)))/multN_srv1(i));
    }
    }
  if(multN_srv1(i)>0)
    {
    effN_srv1(i) = sum(elem_prod(Esrvp1(srv_acyrs1(i)),(1-Esrvp1(srv_acyrs1(i)))))/sum(square(srvp1(i)-Esrvp1(srv_acyrs1(i))));
	}
	}
  loglik(5) = sum(llsrvp1);
  for (i=1;i<=nyrslen_srv1;i++)
    {
    llsrvlenp1(i) = 0;
  for (j=1;j<=nbins3;j++)
    {
      llsrvlenp1(i) += multNlen_srv1(i)*(srvlenp1(i,j)+o)*log((Esrvlenp1(srv_lenyrs1(i),j)+o)/(srvlenp1(i,j)+o));
    }}
   loglik(6) = sum(llsrvlenp1);
  loglik(7) = -.5*norm2(elem_div(
       (log(indxsurv2)-log(Eindxsurv2(srvyrs2))+square(indxsurv_log_sd2)/2.),indxsurv_log_sd2));
  RMSE_srv2= sqrt(norm2(log(indxsurv2)-log(Eindxsurv2(srvyrs2))+square(indxsurv_log_sd2)/2.)/nyrs_srv2);
  for (i=1;i<=nyrsac_srv2;i++)
    {
    llsrvp2(i) = 0;
  for (j=ac_yng_srv2(i);j<=ac_old_srv2(i);j++)
    {
      llsrvp2(i) += multN_srv2(i)*(srvp2(i,j)+o)*log((Esrvp2(srv_acyrs2(i),j)+o)/(srvp2(i,j)+o));
	  res_srv2(i,j)=srvp2(i,j);
      res_srv2(i,trmage-rcrage+j+1)=Esrvp2(srv_acyrs2(i),j);
	if(multN_srv2(i)>0)
    {
	// pearson_srv2(i,j)=(srvp2(i,j)-Esrvp2(srv_acyrs2(i),j))/sqrt((Esrvp2(srv_acyrs2(i),j)*(1.-Esrvp2(srv_acyrs2(i),j)))/multN_srv2(i));
    }
    }
  if(multN_srv2(i)>0)
    {
    effN_srv2(i) = sum(elem_prod(Esrvp2(srv_acyrs2(i)),(1-Esrvp2(srv_acyrs2(i)))))/sum(square(srvp2(i)-Esrvp2(srv_acyrs2(i))));
	}
	}
  loglik(8) = sum(llsrvp2);
  for (i=1;i<=nyrslen_srv2;i++)
    {
    llsrvlenp2(i) = 0;
  for (j=1;j<=nbins2;j++)
    {
      llsrvlenp2(i) += multNlen_srv2(i)*(srvlenp2(i,j)+o)*log((Esrvlenp2(srv_lenyrs2(i),j)+o)/(srvlenp2(i,j)+o));
    }}
   loglik(9) = sum(llsrvlenp2);
  loglik(10) = 0;
    loglik(11) = -.5*norm2(elem_div(
       (log(indxsurv3)-log(Eindxsurv3(srvyrs3))+square(indxsurv_log_sd3)/2.),indxsurv_log_sd3));
    RMSE_srv3= sqrt(norm2(log(indxsurv3)-log(Eindxsurv3(srvyrs3))+square(indxsurv_log_sd3)/2.)/nyrs_srv3);
	 //age composition
  for (i=1;i<=nyrsac_srv3;i++)
    {
    llsrvp3(i) = 0;
  for (j=rcrage;j<=trmage;j++)
    {
      llsrvp3(i) += multN_srv3(i)*(srvp3(i,j)+o)*log((Esrvp3(srv_acyrs3(i),j)+o)/(srvp3(i,j)+o));
	  res_srv3(i,j)=srvp3(i,j);
      res_srv3(i,trmage-rcrage+j+1)=Esrvp3(srv_acyrs3(i),j);
	if(multN_srv3(i)>0)
    {
	// pearson_srv3(i,j)=(srvp3(i,j)-Esrvp3(srv_acyrs3(i),j))/sqrt((Esrvp3(srv_acyrs3(i),j)*(1.-Esrvp3(srv_acyrs3(i),j)))/multN_srv3(i));
    }
    }
  if(multN_srv3(i)>0)
    {
    effN_srv3(i) = sum(elem_prod(Esrvp3(srv_acyrs3(i)),(1-Esrvp3(srv_acyrs3(i)))))/sum(square(srvp3(i)-Esrvp3(srv_acyrs3(i))));
	}
	}
   loglik(12) = sum(llsrvp3);
  for (i=1;i<=nyrslen_srv3;i++)
    {
    llsrvlenp3(i) = 0;
  for (j=1;j<=nbins2;j++)
    {
      llsrvlenp3(i) += multNlen_srv3(i)*(srvlenp3(i,j)+o)*log((Esrvlenp3(srv_lenyrs3(i),j)+o)/(srvlenp3(i,j)+o));
      res_srv3len(i,j)=srvlenp3(i,j);
      res_srv3len(i,nbins2+j)=Esrvlenp3(srv_lenyrs3(i),j);
	if(multNlen_srv3(i)>0)
    {
	// pearson_srv3len(i,j)=(srvlenp3(i,j)-Esrvlenp3(srv_lenyrs3(i),j))/sqrt((Esrvlenp3(srv_lenyrs3(i),j)*(1.-Esrvlenp3(srv_lenyrs3(i),j)))/multNlen_srv3(i));
    }
    }}
   loglik(13) = sum(llsrvlenp3);
   loglik(14) = 0;
   loglik(15)=0;
   loglik(16)=0;
   loglik(17) = 0;
  loglik(18)= 0;
  loglik(18) += -0.5*square(dev_log_recruit(styr)/1.0);
  loglik(18) += -0.5*norm2(dev_log_recruit(styr+1,styr+7)/1.0);
  loglik(18) += -0.5*norm2(dev_log_recruit(endyr-1,endyr)/1.0);
  loglik(19)  = -0.5*norm2(elem_div(first_difference(slp1_fsh_dev),rwlk_sd_short(styr,endyr_fsh_dev-1)));
  loglik(19) += -0.5*norm2(elem_div(first_difference(inf1_fsh_dev),4.0*rwlk_sd_short(styr,endyr_fsh_dev-1)));
  loglik(19) += -0.5*norm2(elem_div(first_difference(slp2_fsh_dev),rwlk_sd_short(styr,endyr_fsh_dev-1)));
  loglik(19) += -0.5*norm2(elem_div(first_difference(inf2_fsh_dev),4.0*rwlk_sd_short(styr,endyr_fsh_dev-1)));
  if(last_phase())
  {
  loglik(20) =  -(1/(2.0*sigmasq_recr))*norm2(log_recr_proj - log_mean_recr_proj);
  }
  else
  {
  loglik(20)=0;
  }
  loglik(21)= 0;
   loglik(22)= -(1/(2.0*(square(0.0244)+square(0.000001))))*square(log_q1_dy- log_q1_ek-.124);
  objfun = -sum(loglik);
    var_prof=Esumbio(endyr);
}

dvariable model_parameters::solve_for_fishing_mortality(dvariable catch_level, int naa_yr)
{
  ofstream& report1= *pad_report1;
    RETURN_ARRAYS_INCREMENT();
    // given a catch amount (in metric tonnes) for the following year, calculate the associate fishing mortality
    int a, i = 0;
    dvariable upperFy, lowerFy, testFy;
    dvariable Ftmp, function_value, exp_value, wt_a;
    // default value
    Ftmp = 1.0e-6;
    // initial values for Fy search
    upperFy = 2.0;
    lowerFy = 0.0;
    // find Fy here by BISECTION method
    // function_val should be 0 at Fy (see eqn B3, Punt 1995)
    testFy = (upperFy + lowerFy) / 2.0;
    function_value = 0.0;
    for (a = rcrage; a <= trmage; a++)
    {
        exp_value = M + (slctfsh_proj(a) * testFy);
        function_value += (wt_fsh_proj(a) * slctfsh_proj(a) * testFy * N_proj(naa_yr,a) * (1.0 - mfexp(-exp_value)) / exp_value);
    }
    function_value *= (1.0e9/1.0e3);    // convert to metric tonnes
    function_value -= catch_level;
    while (fabs(function_value) > BISECT_TOL && i < MAX_BISECT_ITER)
    {
        // printf("iteration %d function value %f Fy %f\n", i, function_value, testFy);
        if (function_value > 0.0)
        {
            upperFy = testFy;
        }
        else
        {
            lowerFy = testFy;
        }
        testFy = (upperFy + lowerFy) / 2.0;
        function_value = 0.0;
        for (a = rcrage; a <= trmage; a++)
        {
            exp_value = M + (slctfsh_proj(a) * testFy);
            function_value += (wt_fsh_proj(a) * slctfsh_proj(a) * testFy * N_proj(naa_yr,a) * (1.0 - mfexp(-exp_value)) / exp_value);
        }
        function_value *= (1.0e9/1.0e3);    // convert to metric tonnes
        function_value -= catch_level;
        i++;
    }
    Ftmp = testFy;
    RETURN_ARRAYS_DECREMENT();
    return(Ftmp);
}

void model_parameters::calculate_curr_bio_ref_points(dvar_vector sel, int curr_yr)
{
  ofstream& report1= *pad_report1;
    // calculate the biological reference points for curr_yr+1
    // NOTE:  this assumes that N(curr_yr+1,a) has been filled in and
    // is the numbers-at-age at the beginning of year (curr_yr+1)
    int a, i, next_yr;
    dvariable spbio_tmp;
    dvariable f_target, spbio_target, spbio_floor;
    dvariable curr_spbio, curr_totbio;
    dvar_vector Z_ABC(rcrage,trmage), Z_OFL(rcrage,trmage);
    F100.initialize();
    SB100.initialize();
    SBtarget.initialize();
    F40.initialize();
    SB40.initialize();
    F35.initialize();
    SB35.initialize();
    F20.initialize();
    SB20.initialize();
    F_ABC.initialize();
    ABC.initialize();
    F_OFL.initialize();
    OFL.initialize();
    Z_ABC.initialize();
    Z_OFL.initialize();
    next_yr = curr_yr + 1;
    F100  = get_spr_rates(1.00, sel, curr_yr);
    SB100 = SBcurr;
    F40   = get_spr_rates(0.40, sel, curr_yr);
    SB40  = SBcurr;
    F35   = get_spr_rates(0.35, sel, curr_yr);
    SB35  = SBcurr;
    F20   = get_spr_rates(0.20, sel, curr_yr);
    SB20  = SBcurr;
    SBtarget = SB40 * (F35 / F40);
    // initialize
    F_ABC = F40;
    // iterate 10 times to get stable values for F_ABC and F_OFL (number matches MWD's code above)
    for (i = 1; i <= 10; i++)
    {
        f_target     = F40;
        spbio_target = SBtarget;
        spbio_floor  = SB20;
        Esumbio(next_yr) = sum(elem_prod(wt_pop_proj((rcrage+2),trmage),N(next_yr)((rcrage+2),trmage)));
        Espawnbio(next_yr) = 0.0;
        for (a = rcrage; a <= trmage; a++)
        {
            Espawnbio(next_yr) += (0.5 * wt_spawn_proj(a) * N(next_yr,a) * mat(a) * mfexp(-0.21 * (M + (sel(a) * F_ABC))));
        }
        spbio_tmp = Espawnbio(next_yr);
        if (spbio_tmp < (0.25 * spbio_floor))
        {
            F_ABC = 0.0;
        }
        else if (spbio_tmp < spbio_target)
        {
            F_ABC = f_target * ((spbio_tmp / spbio_target) - Tier3_alpha) / (1.0 - Tier3_alpha);
        }
        else
        {
            F_ABC = f_target;
        }
        if (spbio_tmp <= spbio_floor)
        {
            cout << "in pk10_1::calculate_curr_bio_ref_points:\tSB <= SB20% in year " << next_yr << ":\tSB20 "<< spbio_floor << "\tSpawning biomass " << Espawnbio(next_yr) << "\tTotal biomass " << Esumbio(next_yr) << endl;
        }
        f_target     = F35;
        spbio_target = SB40;
        if (spbio_tmp < spbio_target)
        {
            F_OFL = f_target * ((spbio_tmp / spbio_target) - Tier3_alpha) / (1.0 - Tier3_alpha);
        }
        else
        {
            F_OFL = f_target;
        }
    }
    Z_ABC = (sel * F_ABC) + M;
    Z_OFL = (sel * F_OFL) + M;
    OFL = ABC = 0.0;
    for (a = rcrage; a <= trmage; a++)
    {
        ABC += (wt_fsh_proj(a) * N(next_yr,a) * (1.0 - mfexp(-Z_ABC(a))) * (sel(a) * F_ABC) / (Z_ABC(a)));
        OFL += (wt_fsh_proj(a) * N(next_yr,a) * (1.0 - mfexp(-Z_OFL(a))) * (sel(a) * F_OFL) / (Z_OFL(a)));
    }
    // convert from millions of metric tonnes to metric tonnes
    ABC *= 1000000.0;
    OFL *= 1000000.0;
}

dvariable model_parameters::get_spr_rates(dvariable spr_percent, dvar_vector sel, int curr_yr)
{
  ofstream& report1= *pad_report1;
    RETURN_ARRAYS_INCREMENT();
    dvariable df=1.0e-3;
    dvariable F1, F2, F3;
    dvariable yld1, yld2, yld3;
    dvariable dyld, dyldp;
    F1.initialize();
    F1 = 0.2; // starting point for Fspr...depends on species...maybe M
    // Newton Raphson stuff to go here
    for (int ii=1; ii<=10; ii++) // arbitrary fixed intervals
    {
        F2     = F1 + df;
        F3     = F1 - df;
        yld1   = -1000.0 * square(log(spr_percent/spr_ratio(F1, sel, curr_yr)));
        yld2   = -1000.0 * square(log(spr_percent/spr_ratio(F2, sel, curr_yr)));
        yld3   = -1000.0 * square(log(spr_percent/spr_ratio(F3, sel, curr_yr)));
        dyld   = (yld2 - yld3) / (2.0 * df);           // First derivative (to find the root of this)
        dyldp  = (yld3 - (2.0 * yld1) + yld2) / (df*df);  // Newton-Raph approximation 2nd deriv
        F1    -= (dyld / dyldp);
    }
    RETURN_ARRAYS_DECREMENT();
    return(F1);
}

dvariable model_parameters::spr_unfished(int curr_yr)
{
  ofstream& report1= *pad_report1;
    RETURN_ARRAYS_INCREMENT();
    int jj, yy, lb, ub;
    dvariable Ntmp, SBtmp;
    Ntmp.initialize();
    SBtmp.initialize();
    // avgR is the average number of age rcrage FEMALE recruits during the period (hcr_styr+rcrage) through curr_yr-1
    lb = hcr_styr + rcrage;
    ub = curr_yr - 1;
    avgR = mean(all_recruits(lb,ub));
    avgR_CV = std_dev(log(all_recruits(lb,ub)));
    // cout << "in spr_unfished: avgR = " << avgR << endl;
    Ntmp = 0.5 * avgR;
    SBtmp = 0.0;
    for (jj=rcrage; jj < trmage; jj++)
    {
        SBtmp += (Ntmp * mat(jj) * wt_spawn_proj(jj) * mfexp(-0.21 * M));
        Ntmp  *= (mfexp(-1.0 * M));
    }
    Ntmp /= (1.0 - mfexp(-1.0 * M));
    SBtmp += (Ntmp * mat(jj) * wt_spawn_proj(jj) * mfexp(-0.21 * M));
    // cout << "in spr_unfished: SBtmp = " << SBtmp << endl;
    RETURN_ARRAYS_DECREMENT();
    return(SBtmp);
}

dvariable model_parameters::spr_ratio(dvariable trial_F, dvar_vector sel, int curr_yr)
{
  ofstream& report1= *pad_report1;
    RETURN_ARRAYS_INCREMENT();
    /* uses following globals:
      wt(1,nages)       wt at age (female spawning)
      natmort(1,nages)  natural mortality
      p_mature(1,nages)  proportion mature (females)
      yrfrac   Fraction of year which defines peak spawning
      phizero  Spawning biomass per recruit with no fishing...
    */
    dvar_vector Ntmp(rcrage,trmage), srvtmp(rcrage,trmage);
    int jj=1;
    Ntmp.initialize();
    srvtmp.initialize();
    SBcurr = 0.0;
    for (jj = rcrage; jj <= trmage; jj++)
    {
        srvtmp(jj) = (sel(jj) * trial_F) + M;
    }
    SB0 = spr_unfished(curr_yr);
    phi0 = SB0 / avgR;  // SBPR
    // cout << "in spr_ratio: phi0 = " << phi0 << endl;
    jj = rcrage;
    Ntmp(jj) = 0.5 * avgR;
    SBcurr  += (Ntmp(jj) * mat(jj) * wt_spawn_proj(jj) * mfexp(-0.21 * srvtmp(jj)));
    for (jj=rcrage+1 ; jj < trmage; jj++)
    {
        Ntmp(jj) = Ntmp(jj-1) * mfexp(-1.0 * srvtmp(jj-1));
        SBcurr  += (Ntmp(jj) * mat(jj) * wt_spawn_proj(jj) * mfexp(-0.21 * srvtmp(jj)));
    }
    Ntmp(trmage) = (Ntmp(trmage-1) * mfexp(-1.0 * srvtmp(trmage-1)) / (1.0 - mfexp(-1.0 * srvtmp(trmage))));
    SBcurr += (Ntmp(trmage) * mat(trmage) * wt_spawn_proj(trmage) * mfexp(-0.21 * srvtmp(trmage)));
    RETURN_ARRAYS_DECREMENT();
    return(SBcurr / SB0);
}

void model_parameters::MCMC_output(void)
{
  ofstream& report1= *pad_report1;
  if(mceval_phase())
  {
  report1<<mean(recruit(1978,endyr-1));
  report1<<" ";
  report1<<endyr-2;
  report1<<" ";
  report1<<F(endyr-2);
  report1<<" ";
  report1<<Espawnbio(endyr-2);
  report1<<" ";
  report1<<recruit(endyr-2);
  report1<<" ";
  report1<<endyr-1;
  report1<<" ";
  report1<<F(endyr-1);
  report1<<" ";
  report1<<Espawnbio(endyr-1);
  report1<<" ";
  report1<<recruit(endyr-1);
  report1<<" ";
  report1<<endyr;
  report1<<" ";
  report1<<F(endyr);
  report1<<" ";
  report1<<Espawnbio(endyr);
  report1<<" ";
  report1<<recruit(endyr);
  report1<<" ";
  report1<<endyr+1;
  report1<<" ";
  report1<<F_proj(endyr+1);
  report1<<" ";
  report1<<Espawnbio_proj(endyr+1);
  report1<<" ";
  report1<<recruit_proj(endyr+1);
  report1<<" ";
  report1<<endyr+2;
  report1<<" ";
  report1<<F_proj(endyr+2);
  report1<<" ";
  report1<<Espawnbio_proj(endyr+2);
  report1<<" ";
  report1<<recruit_proj(endyr+2);
  report1<<" ";
  report1<<endyr+3;
  report1<<" ";
  report1<<F_proj(endyr+3);
  report1<<" ";
  report1<<Espawnbio_proj(endyr+3);
  report1<<" ";
  report1<<recruit_proj(endyr+3);
  report1<<" ";
  report1<<endyr+4;
  report1<<" ";
  report1<<F_proj(endyr+4);
  report1<<" ";
  report1<<Espawnbio_proj(endyr+4);
  report1<<" ";
  report1<<recruit_proj(endyr+4);
  report1<<" ";
  report1<<endyr+5;
  report1<<" ";
  report1<<F_proj(endyr+5);
  report1<<" ";
  report1<<Espawnbio_proj(endyr+5);
  report1<<" ";
  report1<<recruit_proj(endyr+5);
  report1 << endl;
  }
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1.e0, 1.e-1, 1.e-4, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{1000, 1000, 1000, 1000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  if (last_phase())
  {
      calculate_curr_bio_ref_points(slctfsh_proj,endyr);
  }
  report << "Objective function" << endl;
  report << objfun << endl;
  report << "Likelihood components" << endl;
  report << loglik << endl;
  report << " " << endl;
  report << "Natural mortality" << endl;
  report << M << endl;
  report << "Ecosystem comsumption" << endl;
  report << Eecocon << endl;
  report << "Initial age comp" << endl;
  report << initN << endl;
  report << "Recruits" << endl;
  report << recruit << endl;
  report << " " << endl;
  report << "Total catch" << endl;
  report << cattot << endl;
  report << "Expected total catch" << endl;
  report << Ecattot << endl;
  report << "Fishing mortalities" << endl;
  report << F << endl;
  report << "Base fishery selectivity" << endl;
  report << slctfsh_base << endl;
  report << "Selectivity means" << endl;
  report << log_slp1_fsh_mean << endl;
  report << inf1_fsh_mean << endl;
  report << log_slp2_fsh_mean << endl;
  report << inf2_fsh_mean << endl;
  report << "Selectivity deviances" << endl;
  report << slp1_fsh_dev << endl;
  report << inf1_fsh_dev << endl;
  report << slp2_fsh_dev << endl;
  report << inf2_fsh_dev << endl;
  report << "Selectivity vectors" << endl;
  report <<  slp1_fsh << endl;
  report <<  inf1_fsh << endl;
  report <<  slp2_fsh << endl;
  report <<  inf2_fsh << endl;
  report << "Fishery selectivity" << endl;
  report << slctfsh << endl;
  report << "Fishery age composition likelihoods" << endl;
  report << llcatp << endl;
  report << "Fishery  age composition" << endl;
  report << catp << endl;
  report << "Expected fishery age composition" << endl;
  report << Ecatp << endl;
  report << "Observed and expected age comp" << endl;
  report << res_fish << endl;
  report << "Pearson residuals age comp" << endl;
  report << pearson_fish << endl;
  report << "Input N" << endl;
  report << multN_fsh << endl;
  report << "Effective N age comp" << endl;
  report << effN_fsh << endl;
  report << "Fishery length composition likelihoods" << endl;
  report << lllenp << endl;
  report << "Fishery length composition" << endl;
  report << lenp << endl;
  report << "Expected length composition" << endl;
  report << Elenp << endl;
  report << " " << endl;
  report << "Survey 1 q" << endl;
  report << q1_bs << endl;
  report << q1_ek << endl;
  report << q1_dy << endl;
  report << "Selectivity parameters" << endl;
  report << log_slp2_srv1 << endl;
  report << inf2_srv1 << endl;
  report << "Survey 1 selectivity" << endl;
  report << slctsrv1 << endl;
  report << "Expected survey 1 index" << endl;
  report << Eindxsurv1_bs << endl;
  report << Eindxsurv1_ek << endl;
  report << Eindxsurv1_dy << endl;
  report << "RMSE" << endl;
  report << RMSE_srv1_bs << endl;
  report << RMSE_srv1_ek << endl;
  report << RMSE_srv1_dy << endl;
  report << "Survey 1 age composition likelihoods" << endl;
  report << llsrvp1 << endl;
  report << "Survey 1 age composition" << endl;
  report << srvp1 << endl;
  report << "Expected survey 1 age composition" << endl;
  report << Esrvp1 << endl;
  report << "Observed and expected age comp" << endl;
  report << res_srv1 << endl;
  report << "Pearson residuals age comp" << endl;
  report << pearson_srv1 << endl;
  report << "Input N" << endl;
  report << multN_srv1 << endl;
  report << "Effective N age comp" << endl;
  report << effN_srv1 << endl;
  report << "Survey 1 length composition likelihoods" << endl;
  report << llsrvlenp1 << endl;
  report << "Survey 1 length composition" << endl;
  report << srvlenp1 << endl;
  report << "Expected survey 1 length composition" << endl;
  report << Esrvlenp1 << endl;
  report << " " << endl;
  report << "Survey 2 q" << endl;
  report << q2 << endl;
  report << "Selectivity parameters" << endl;
  report << log_slp1_srv2 << endl;
  report << inf1_srv2 << endl;
  report << log_slp2_srv2 << endl;
  report << inf2_srv2 << endl;
  report << "Survey 2 selectivity" << endl;
  report << slctsrv2 << endl;
  report << "Expected survey 2 index" << endl;
  report << Eindxsurv2 << endl;
  report << "RMSE" << endl;
  report << RMSE_srv2 << endl;
  report << "Survey 2 age composition likelihoods" << endl;
  report << llsrvp2 << endl;
  report << "Survey 2 age composition" << endl;
  report << srvp2 << endl;
  report << "Expected survey 2 age composition" << endl;
  report << Esrvp2 << endl;
  report << "Observed and expected age comp" << endl;
  report << res_srv2 << endl;
  report << "Pearson residuals age comp" << endl;
  report << pearson_srv2 << endl;
  report << "Input N" << endl;
  report << multN_srv2 << endl;
  report << "Effective N age comp" << endl;
  report << effN_srv2 << endl;
  report << "Survey 2 length composition likelihoods" << endl;
  report << llsrvlenp2 << endl;
  report << "Survey 2 length composition" << endl;
  report << srvlenp2 << endl;
  report << "Expected survey 2 length composition" << endl;
  report << Esrvlenp2 << endl;
  report << " " << endl;
   report << "Survey 3 q" << endl;
  report << q3 << endl;
  report << "Selectivity parameters" << endl;
  report << log_slp1_srv3 << endl;
  report << inf1_srv3 << endl;
  report << "Survey 3 selectivity" << endl;
  report << slctsrv3 << endl;
  report << "Expected survey 3 index" << endl;
  report << Eindxsurv3 << endl;
  report << "RMSE" << endl;
  report << RMSE_srv3 << endl;
  report << "Survey 3 age composition likelihoods" << endl;
  report << llsrvp3 << endl;
  report << "Survey 3 age composition" << endl;
  report << srvp3 << endl;
  report << "Expected survey 3 age composition" << endl;
  report << Esrvp3 << endl;
  report << "Observed and expected age comp" << endl;
  report << res_srv3 << endl;
  report << "Pearson residuals age comp" << endl;
  report << pearson_srv3 << endl;
  report << "Input N" << endl;
  report << multN_srv3 << endl;
  report << "Effective N age comp" << endl;
  report << effN_srv3 << endl;
  report << "Survey 3 length composition likelihoods" << endl;
  report << llsrvlenp3 << endl;
  report << "Survey 3 length composition" << endl;
  report << srvlenp3 << endl;
  report << "Expected survey 3 length composition" << endl;
  report << Esrvlenp3 << endl;
  report << "Observed and expected length comp" << endl;
  report << res_srv3len << endl;
  report << "Pearson residuals length comp" << endl;
  report << pearson_srv3len << endl;
  report << "Expected summary (age 3+) biomass" << endl;
  report << Esumbio << endl;
  report << "Expected spawning biomass" << endl;
  report << Espawnbio << endl;
  report << "Numbers at age" << endl;
  report << N << endl;
  report << endl;
  report << "Projection output" << endl;
  report << "Recruits" << endl;
  report << recruit_proj << endl;
  report << "Log recruitment" << endl;
  report << log_recr_proj << endl;
  report << "Variances" << endl;
  report <<  sigmasq_recr<< endl;
  report << "Numbers at age" << endl;
  report << N_proj << endl;
  report << "Catch at age" << endl;
  report << C_proj << endl;
  report << "Survey numbers at age" << endl;
  report << Nsrv_proj << endl;
  report << "Weight at age" << endl;
  report << "Population" << endl;
  report <<     wt_pop_proj << endl;
  report << "Spawning" << endl;
  report <<     wt_spawn_proj << endl;
  report << "Fishery" << endl;
  report <<     wt_fsh_proj << endl;
  report << "Fishery selectivity" << endl;
  report <<    slctfsh_proj << endl;
  report << "Total catches & summary biomass" << endl;
  report << "Total catches" << endl;
  report <<     Ecattot_proj << endl;
  report << "Summary biomass" << endl;
  report <<     Esumbio_proj << endl;
  report << "Spawning biomass" << endl;
  report <<     Espawnbio_proj << endl;
  report << "Survey biomass" << endl;
  report <<     Esrv_proj << endl;
  report << "Ftarget B40" << endl;
  report <<     Ftarget << endl;
  report <<     B40 << endl;
  report << "Fishing mortality" << endl;
  report <<     F_proj << endl;
  report << endl;
  if (endyr > base_endyr)
  {
      report << "Future recruitment" << endl;
      report << N((base_endyr+1),rcrage);
      if (endyr > (base_endyr+1))
      {
          for (i = (base_endyr+2); i <= endyr; i++)
          {
              report << " " << N(i,rcrage);
          }
      }
      report << endl;
      report << endl;
  }
  if (last_phase())
  {
      report << endl;
      report << "Values of biological reference points for year " << (endyr+1) << endl;
      report << "xx\tFxx\tSBxx" << endl;
      report << "1.00\t" << F100 << "\t" << SB100 << endl;
      report << "0.40\t" << F40 << "\t" << SB40 << endl;
      report << "0.35\t" << F35 << "\t" << SB35 << endl;
      report << "0.20\t" << F20 << "\t" << SB20 << endl;
      report << endl;
      report << "SB0 = " << SB0 << ", avgR = " << avgR << ", avgR_CV = " << avgR_CV << ", phi0 = " << phi0 << endl;
      report << "SBtarget = " << SBtarget << ", which is SB" << (100.0 * (SBtarget / (SB40 / 0.4))) << endl;
      report << endl;
      report << "Decision rule parameter values for year " << (endyr+1) << endl;
      report << "Estimated " << (endyr+1) << " spawning biomass\t" << Espawnbio(endyr+1) << endl;
      report << "Estimated " << (endyr+1) << " age 3+ biomass\t" << Esumbio(endyr+1) << endl;
      report << "F_ABC\t" << F_ABC << endl;
      report << "ABC\t" << ABC << endl;
      report << "F_OFL\t" << F_OFL << endl;
      report << "OFL\t" << OFL << endl;
      report << endl;
  }
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_report1;
  pad_report1 = NULL;
}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
 arrmblsize = 3000000;
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
 gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
