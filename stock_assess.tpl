// GOA pollock age-structured model
// One fishery, three surveys: acoustic, bottom trawl, ADFG crab/groundfish
// Double logistic selectivity for fishery
// Random walks in selectivity
// Logistic or double logistic for surveys
// Five year projection
// Shelikof Strait EIT split into three catchability periods
// Biosonics MillerFreeman and OscarDyson
// Reference model for 2013
//
// 26 March 2014 - edited by ZTA for use in the GOA walleye pollock MSE
//

DATA_SECTION

  !!CLASS ofstream report1("mceval.dat")

  init_int styr                                  // Starting year for population model
  init_int endyr                                 // Ending year for population model
  init_int base_endyr                            // Ending year for base model
  init_int hcr_styr                              // Starting year recruitment calcs for the harvest control rule
  init_int rcrage                                // Recruitment age
  init_int trmage                                // Last modeled age
  init_int nbins1                                // Number of length bins in transitiom matrix 1
  init_int nbins2                                // Number of length bins in transitiom matrix 2
  init_int nbins3                                // Number of length bins in transitiom matrix 3

//Fishery
  init_vector cattot(styr,endyr)                 // Total catch in tons
  init_vector cattot_log_sd(styr,endyr)          // Total catch (cv) = sdev of log(cattot)

  init_int nyrs_fsh                              // Number of fishery age comps
  init_ivector fshyrs(1,nyrs_fsh)                // Years for the fishery age comps
  init_vector multN_fsh(1,nyrs_fsh)              // Multinomial sample size by year
  init_ivector ac_yng_fsh(1,nyrs_fsh)            // Accumulation of lower ages
  init_ivector ac_old_fsh(1,nyrs_fsh)            // Accumulation of upper ages

  init_int nyrslen_fsh                           // Number of fishery length comps
  init_ivector fshlenyrs(1,nyrslen_fsh)          // Years for the fishery length comps
  init_vector multNlen_fsh(1,nyrslen_fsh)        // Multinomial sample size by year

  init_vector rwlk_sd(styr,endyr-1)              // Random walk stdevs
  vector rwlk_sd_short(styr,base_endyr)
  init_matrix catp(1,nyrs_fsh,rcrage,trmage)     // Catch proportions at age
  init_matrix lenp(1,nyrslen_fsh,1,nbins1)       // Catch proportions at age

  init_matrix wt_fsh(styr,endyr,rcrage,trmage)   // Weight at age by year
                                                 // For matrices indices for rows, then col
//Survey 1 (Acoustic)
//Biosonics
  init_int nyrs_srv1_bs                             // Number of survey biomass estimates
  init_ivector srvyrs1_bs(1,nyrs_srv1_bs)           // Years in which surveys occured
  init_vector indxsurv1_bs(1,nyrs_srv1_bs)          // Survey index
  init_vector indxsurv_log_sd1_bs(1,nyrs_srv1_bs)   // Survey index (cv) = sdev of log(indxsurv)
//EK500
  init_int nyrs_srv1_ek                             // Number of survey biomass estimates
  init_ivector srvyrs1_ek(1,nyrs_srv1_ek)           // Years in which surveys occured
  init_vector indxsurv1_ek(1,nyrs_srv1_ek)           // Survey index
  init_vector indxsurv_log_sd1_ek(1,nyrs_srv1_ek)    // Survey index (cv) = sdev of log(indxsurv)
//Dyson
  init_int nyrs_srv1_dy                             // Number of survey biomass estimates
  init_ivector srvyrs1_dy(1,nyrs_srv1_dy)           // Years in which surveys occured
  init_vector indxsurv1_dy(1,nyrs_srv1_dy)           // Survey index
  init_vector indxsurv_log_sd1_dy(1,nyrs_srv1_dy)    // Survey index (cv) = sdev of log(indxsurv)



  init_vector yrfrct_srv1(styr,endyr)            // Fraction of year to midpoint of survey
  init_int nyrsac_srv1                           // Number of survey age comps
  init_ivector srv_acyrs1(1,nyrsac_srv1)         // Years for the survey age comp
  init_vector multN_srv1(1,nyrsac_srv1)          // Multinomial sample size by year
  init_ivector ac_yng_srv1(1,nyrsac_srv1)        // Accumulation of lower ages
  init_ivector ac_old_srv1(1,nyrsac_srv1)        // Accumulation of upper ages
  init_int nyrslen_srv1                          // Number of survey length comps
  init_ivector srv_lenyrs1(1,nyrslen_srv1)       // Years for the survey length comps
  init_vector multNlen_srv1(1,nyrslen_srv1)      // Multinomial sample size by year
  init_matrix srvp1(1,nyrsac_srv1,rcrage,trmage) // Survey proportions at age
  init_matrix srvlenp1(1,nyrslen_srv1,1,nbins3)  // Survey proportions at length
  init_matrix wt_srv1(styr,endyr,rcrage,trmage)  // Survey weights at age
                                                 // Note full dimensions for weight matrix
//Survey 2 (Bottom trawl)
  init_int nyrs_srv2                             // Number of surveys
  init_ivector srvyrs2(1,nyrs_srv2)              // Years in which surveys occured
  init_vector indxsurv2(1,nyrs_srv2)             // Survey index
  init_vector indxsurv_log_sd2(1,nyrs_srv2)      // Survey index (cv) = sdev of log(indxsurv)
  init_vector yrfrct_srv2(styr,endyr)            // Fraction of year to midpoint of survey
  init_int nyrsac_srv2                           // Number of survey age comps
  init_ivector srv_acyrs2(1,nyrsac_srv2)         // Years for the survey age comp
  init_vector multN_srv2(1,nyrsac_srv2)          // Multinomial sample size by year
  init_ivector ac_yng_srv2(1,nyrsac_srv2)        // Accumulation of lower ages
  init_ivector ac_old_srv2(1,nyrsac_srv2)        // Accumulation of upper ages
  init_int nyrslen_srv2                          // Number of survey length comps
  init_ivector srv_lenyrs2(1,nyrslen_srv2)       // Years for the survey length comps
  init_vector multNlen_srv2(1,nyrslen_srv2)      // Multinomial sample size by year
  init_matrix srvp2(1,nyrsac_srv2,rcrage,trmage) // Survey proportions at age
  init_matrix srvlenp2(1,nyrslen_srv2,1,nbins2)  // Survey proportions at length
  init_matrix wt_srv2(styr,endyr,rcrage,trmage)  // Survey weights at age

//Survey 3 (ADFG coastal survey)

  init_int nyrs_srv3                             // Number of survey biomass estimates
  init_ivector srvyrs3(1,nyrs_srv3)              // Years in which surveys occured
  init_vector indxsurv3(1,nyrs_srv3)             // Survey index
  init_vector indxsurv_log_sd3(1,nyrs_srv3)      // Survey index (cv) = sdev of log(indxsurv)
  init_vector yrfrct_srv3(styr,endyr)            // Fraction of year to midpoint of survey
  init_int nyrsac_srv3                           // Number of survey age comps
  init_ivector srv_acyrs3(1,nyrsac_srv3)         // Years for the survey age comps
  init_vector multN_srv3(1,nyrsac_srv3)          // Multinomial sample size by year
  init_int nyrslen_srv3                          // Number of survey length comps
  init_ivector srv_lenyrs3(1,nyrslen_srv3)       // Years for the survey length comps
  init_vector multNlen_srv3(1,nyrslen_srv3)      // Multinomial sample size by year
  init_matrix srvp3(1,nyrsac_srv3,rcrage,trmage) // Survey proportions at age
  init_matrix srvlenp3(1,nyrslen_srv3,1,nbins2)  // Survey proportions at length
  init_matrix wt_srv3(styr,endyr,rcrage,trmage)  // Survey weights at age

//Age error transition matrix
  init_matrix age_trans(rcrage,trmage,rcrage,trmage)
//Age to length transition matrix
  init_matrix len_trans1(rcrage,trmage,1,nbins1)
  init_matrix len_trans2(rcrage,trmage,1,nbins2)
  init_matrix len_trans3(rcrage,trmage,1,nbins3)

//Population vectors
  init_matrix wt_pop(styr,endyr,rcrage,trmage)   // Population weight at age
  init_matrix wt_spawn(styr,endyr,rcrage,trmage) // Population weight at age at spawning (April 15)
// Anne's maturity vector
  init_vector mat_old(rcrage,trmage)                 // Proportion mature
  init_vector mat(rcrage,trmage)                 // Proportion mature

//Projection parameters
  init_vector wt_pop_proj(rcrage,trmage)         // Projection population weight at age at start of year
  init_vector wt_spawn_proj(rcrage,trmage)       // Projection population weight at age at spawning (April 15)
  init_vector wt_fsh_proj(rcrage,trmage)         // Projection fishery weight at age
  init_vector wt_srv_proj(rcrage,trmage)         // Projection arbitrary survey weight at age
  init_vector Ftarget(endyr+1,endyr+5)
  init_number B40
// mean log recruitment
  init_number   log_mean_recr_proj
 // Variance for log recr, recruitment indices
  init_number sigmasq_recr

  int styr_avg_slct
  int endyr_avg_slct
  int endyr_fsh_dev

  int i                                          // Index for year
  int j                                          // Index for age
  int loop

  int phase_future

  int mcmc_iter

  number y_frac_sp                              // time of peak spawning

  number Tier3_alpha                            // alpha for the NPFMC Tier 3 decision rule, 0.05

  number o                                       // A small number

  !!#define MAX_BISECT_ITER  128            // maximum number of iterations for the bisection
  !!#define BISECT_TOL       0.00001        // maximum tolerance for bisection


 LOCAL_CALCS


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


 END_CALCS


INITIALIZATION_SECTION

  mean_log_initN      0.0   // Mean log initial age composition
  mean_log_recruit    0.0   // Mean log recruitment
  mean_log_F        -1.6   // Mean log fishing mortality

  M                  0.30  // Natural mortality
  log_q1_bs          0.0   // Survey 1 catchability
  log_q1_ek          0.0   // Survey 1 catchability
  log_q1_dy          0.3   // Survey 1 catchability
  log_q2             0.0   // Survey 2 catchability
  log_q3            -1.6   // Survey 3 catchability

  slp1_fsh_dev       0.0   // Selectivity deviance terms
  inf1_fsh_dev       0.0
  slp2_fsh_dev       0.0
  inf2_fsh_dev       0.0

//Starting values for selectivity curves
  log_slp1_fsh_mean  1.0
  inf1_fsh_mean      4.0
  log_slp2_fsh_mean  1.0
  inf2_fsh_mean      8.0

// for constant future fishery selectivity
  log_slp1_fsh_mean_const  1.0
  inf1_fsh_mean_const      4.0
  log_slp2_fsh_mean_const  1.0
  inf2_fsh_mean_const      8.0

  log_slp2_srv1      1.0
  inf2_srv1          9.0

//  log_slp1_srv2      0.0
//  inf1_srv2          4.0
//  log_slp2_srv2      0.0
//  inf2_srv2          8.0

  log_slp1_srv2     -0.8
  inf1_srv2          47.0
  log_slp2_srv2      -.14
  inf2_srv2          7.0


  log_slp1_srv3      0.0
  inf1_srv3          5.0

PARAMETER_SECTION

//Population parameters

  init_bounded_number M(0.1,0.5,-1)

  init_bounded_number mean_log_initN(-15,15,-1)
  init_bounded_dev_vector dev_log_initN(rcrage+1,trmage,-15,15,-2)
  vector initN(rcrage+1,trmage)

  init_bounded_number mean_log_recruit(-15,15,1)
  init_bounded_dev_vector dev_log_recruit(styr,endyr,-15,15,4)

//Forward projections
  init_bounded_vector log_recr_proj(endyr+1,endyr+5,-5,5,10)

  sdreport_vector recruit_proj(endyr+1,endyr+5)
  matrix N_proj(endyr+1,endyr+5,rcrage,trmage)
  vector F_proj(endyr+1,endyr+5)
  matrix Z_proj(endyr+1,endyr+5,rcrage,trmage)
  matrix C_proj(endyr+1,endyr+5,rcrage,trmage)
  matrix Nsrv_proj(endyr+1,endyr+5,rcrage,trmage)
  vector slctfsh_proj(rcrage,trmage)
  vector Ecattot_proj(endyr+1,endyr+5)
  sdreport_vector Esumbio_proj(endyr+1,endyr+5)
  sdreport_vector Espawnbio_proj(endyr+1,endyr+5)
  sdreport_vector Esrv_proj(endyr+1,endyr+5)
  sdreport_vector Exrate_proj(endyr+1,endyr+5)
  number sbio

//Selectivity parameters

//Fishery selectivity

    init_bounded_number log_slp1_fsh_mean(-5,5,4)
    init_bounded_number inf1_fsh_mean(1,5,4)
    init_bounded_number log_slp2_fsh_mean(-5,5,4)
    init_bounded_number inf2_fsh_mean(7,20,4)

    init_bounded_dev_vector slp1_fsh_dev(styr,endyr_fsh_dev,-5,5,5)
    init_bounded_dev_vector inf1_fsh_dev(styr,endyr_fsh_dev,-5,5,5)
    init_bounded_dev_vector slp2_fsh_dev(styr,endyr_fsh_dev,-5,5,5)
    init_bounded_dev_vector inf2_fsh_dev(styr,endyr_fsh_dev,-5,5,5)

    vector slp1_fsh(styr,endyr_fsh_dev)
    vector inf1_fsh(styr,endyr_fsh_dev)
    vector slp2_fsh(styr,endyr_fsh_dev)
    vector inf2_fsh(styr,endyr_fsh_dev)

// for constant future fishery selectivity

    init_bounded_number log_slp1_fsh_mean_const(-5,5,phase_future)
    init_bounded_number inf1_fsh_mean_const(1,5,phase_future)
    init_bounded_number log_slp2_fsh_mean_const(-5,5,phase_future)
    init_bounded_number inf2_fsh_mean_const(7,20,phase_future)

    number slp1_fsh_const
    number inf1_fsh_const
    number slp2_fsh_const
    number inf2_fsh_const

//Acoustic survey selectivity
//    init_bounded_number log_slp1_srv1(-10,5,7)
//    init_bounded_number inf1_srv1(1,20,7)
    init_bounded_number log_slp2_srv1(-5,5,7)
    init_bounded_number inf2_srv1(5,20,7)
	init_bounded_number srv1_age1(0,2,7)
//    init_bounded_number srv1_age2(0,2,7)
//    init_bounded_number srv1_age3(0,2,7)

//Trawl selectivity
    init_bounded_number log_slp1_srv2(-5,5,7)
//    init_bounded_number inf1_srv2(1,50,7)
    init_bounded_number inf1_srv2(1,50,-1)

    init_bounded_number log_slp2_srv2(-5,5,8)
    init_bounded_number inf2_srv2(5,10,8)
    init_bounded_number srv2_age1(0,2,7)

//ADFG selectivity
    init_bounded_number log_slp1_srv3(-5,5,9)
    init_bounded_number inf1_srv3(1,20,9)
//    init_bounded_number log_slp2_srv3(-5,5,3)
//    init_bounded_number inf2_srv3(3,20,3)

//Fishing mortality and survey catchablility

  init_bounded_number mean_log_F(-10,10,1)
  init_bounded_dev_vector dev_log_F(styr,endyr,-10,10,2)
  vector F(styr,endyr)

  init_bounded_number log_q1_bs(-10,10,5)
  init_bounded_number log_q1_ek(-10,10,5)
  init_bounded_number log_q1_dy(-10,10,5)
  init_bounded_number log_q2(-10,10,-1)
//  init_bounded_number log_q2(-10,10,5)
  init_bounded_number log_q3(-10,10,6)

  number q1_bs
  number q1_ek
  number q1_dy
  number q2
  number q3

  number avgR;
  number avgR_CV;
  number SB0;
  number R0;
  number phi0;
  number SBcurr;
  number F100;
  number SB100;
  number SBtarget;
  number F40;
  number SB40;
  number F35;
  number SB35;
  number F20;
  number SB20;
  number F_ABC;
  number ABC;
  number F_OFL;
  number OFL;

  //Dependent parameters

  matrix N(styr,endyr+5,rcrage,trmage)
  sdreport_vector endN(rcrage,trmage)
  vector all_recruits(styr,endyr+5)
  matrix Z(styr,endyr,rcrage,trmage)
  matrix C(styr,endyr,rcrage,trmage)
  matrix Nsrv1(styr,endyr,rcrage,trmage)
  vector slctsrv1(rcrage,trmage)
  matrix Nsrv2(styr,endyr,rcrage,trmage)
  vector slctsrv2(rcrage,trmage)
  matrix Nsrv3(styr,endyr,rcrage,trmage)
  vector slctsrv3(rcrage,trmage)

  vector slctfsh_base(rcrage,trmage)
  matrix slctfsh(styr,endyr,rcrage,trmage)
  vector slctfsh_const(rcrage,trmage)

  vector Eecocon(styr,endyr)
  matrix Eec(styr,endyr,rcrage,trmage)

  vector Ecattot(styr,endyr)
  matrix Ecatp(styr,endyr,rcrage,trmage)
  matrix Elenp(styr,endyr,1,nbins1)
  vector Eindxsurv1_bs(styr,endyr)
  vector Eindxsurv1_ek(styr,endyr)
  vector Eindxsurv1_dy(styr,endyr)
  matrix Esrvp1(styr,endyr,rcrage,trmage)
  matrix Esrvlenp1(styr,endyr,1,nbins3)

  vector Eindxsurv2(styr,endyr)
  matrix Esrvp2(styr,endyr,rcrage,trmage)
  matrix Esrvlenp2(styr,endyr,1,nbins2)

  vector Eindxsurv3(styr,endyr)
  matrix Esrvp3(styr,endyr,rcrage,trmage)
  matrix Esrvlenp3(styr,endyr,1,nbins2)

   vector loglik(1,22)

  vector llcatp(1,nyrs_fsh)
  vector lllenp(1,nyrslen_fsh)

  vector llsrvp1(1,nyrsac_srv1)
  vector llsrvlenp1(1,nyrslen_srv1)

  vector llsrvp2(1,nyrsac_srv2)
  vector llsrvlenp2(1,nyrslen_srv2)

  vector llsrvp3(1,nyrsac_srv3)
  vector llsrvlenp3(1,nyrslen_srv3)

  sdreport_vector recruit(styr,endyr)
  sdreport_number sd_avg_rec

  sdreport_vector Espawnbio(styr,endyr+1)
  sdreport_vector Esumbio(styr,endyr+1)

  //residual output matrices
  matrix res_fish(1,nyrs_fsh,rcrage,2*trmage-rcrage+1)
  matrix res_srv1(1,nyrsac_srv1,rcrage,2*trmage-rcrage+1)
  matrix res_srv2(1,nyrsac_srv2,rcrage,2*trmage-rcrage+1)
  matrix res_srv3(1,nyrsac_srv3,rcrage,2*trmage-rcrage+1)
  matrix res_srv3len(1,nyrslen_srv3,1,2*nbins2)

  matrix pearson_fish(1,nyrs_fsh,rcrage,trmage)
  matrix pearson_srv1(1,nyrsac_srv1,rcrage,trmage)
  matrix pearson_srv2(1,nyrsac_srv2,rcrage,trmage)
  matrix pearson_srv3(1,nyrsac_srv3,rcrage,trmage)
  matrix pearson_srv3len(1,nyrslen_srv3,1,nbins2)

  vector effN_fsh(1,nyrs_fsh)
  vector effN_srv1(1,nyrsac_srv1)
  vector effN_srv2(1,nyrsac_srv2)
  vector effN_srv3(1,nyrsac_srv3)

  number RMSE_srv1_bs
  number RMSE_srv1_ek
  number RMSE_srv1_dy
  number RMSE_srv2
  number RMSE_srv3

//Objective function

  likeprof_number var_prof
  objective_function_value objfun

PRELIMINARY_CALCS_SECTION

// Do upper and lower accumulations in age composition data

// Fishery
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

// Survey 1
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

// Survey 2
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
//  var_prof.set_stepnumber(10);
  var_prof.set_stepsize(0.1);


PROCEDURE_SECTION
//Use C++ syntax for the procedure section

//Calls to functions

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

FUNCTION Convert_log_parameters

// Assume close to equilibrium at f=0 at the start
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

FUNCTION Selectivity

// Fishery selectivity

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
//   slctfsh(i)=slctfsh(i)/max(slctfsh(i));
// The plan would be to check and adjust the max selected age as needed
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

// Mean selectivity for SPR calcs and projections
//   slctfsh(endyr,2)=0.0327824;
//   slctfsh(endyr,3)=0.127648;
//   slctfsh(endyr,4)=0.386186;
//   slctfsh(endyr,5)=0.731703;
//   slctfsh(endyr,6)=0.924972;
//   slctfsh(endyr,7)=0.98532;
//   slctfsh(endyr,8)=1;
//   slctfsh(endyr,9)=0.977554;
//   slctfsh(endyr,10)=0.426926;


//Survey 1 selectivity
  for (j=rcrage;j<=trmage;j++)
    {

    slctsrv1(j) = (1-1/(1+mfexp(-mfexp(log_slp2_srv1)*(double(j)-inf2_srv1))));
//    mfexp(-mfexp(log_slp2_srv1)*(double(j)-2));
    }

    slctsrv1=slctsrv1/slctsrv1(2);
	slctsrv1(rcrage)=srv1_age1;
//    slctsrv1(rcrage)=srv1_age2;
//    slctsrv1(rcrage+1)=srv1_age3;

//Survey 2 selectivity
  for (j=rcrage;j<=trmage;j++)
    {
    slctsrv2(j) = (1/(1+mfexp(-mfexp(log_slp1_srv2)*(double(j)-inf1_srv2))))
  *(1-1/(1+mfexp(-mfexp(log_slp2_srv2)*(double(j)-inf2_srv2))));
//    mfexp(-mfexp(log_slp2_srv2)*(double(j)-2));
    }
    slctsrv2=slctsrv2/slctsrv2(7);
	slctsrv2(rcrage)=srv2_age1;

//Survey 3 selectivity
  for (j=rcrage;j<=trmage;j++)
    {
    slctsrv3(j) = (1/(1+mfexp(-mfexp(log_slp1_srv3)*(double(j)-inf1_srv3))));
//*(1-1/(1+mfexp(-mfexp(log_slp2_srv3)*(double(j)-inf2_srv3))));
    }
    slctsrv3=slctsrv3/slctsrv3(10);

FUNCTION Mortality

  for (i=styr;i<=endyr;i++)
    {
  for (j=rcrage;j<=trmage;j++)
    {
    Z(i,j)=(F(i)*slctfsh(i,j))+M;
    }}


FUNCTION Numbers_at_age

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

FUNCTION Catch_at_age

// Catch at age and survey numbers at age

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

FUNCTION Expected_values

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

// 3+ biomass
    Esumbio(i)= N(i)(rcrage+2,trmage)*wt_pop(i)(rcrage+2,trmage);
// Alternative: 2+ biomass
//    Esumbio(i)= N(i)(rcrage+1,trmage)*wt_pop(i)(rcrage+1,trmage);
// Total biomass
//    Esumbio(i)= N(i)(rcrage,trmage)*wt_pop(i)(rcrage,trmage);
// Alternative:Total biomass at spawning (use for apportioning to management area)
//    Esumbio(i)= sum(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv1(i)*Z(i))),wt_srv1(i)));
// 2+ biomass at spawning (use for apportioning to management area)
//    Esumbio(i)= sum(elem_prod(elem_prod(N(i)(rcrage+1,trmage),mfexp(-yrfrct_srv1(i)*Z(i)(rcrage+1,trmage))),wt_srv1(i)(rcrage+1,trmage)));
// 3+ biomass at spawning
//    Esumbio(i)= sum(elem_prod(elem_prod(N(i)(rcrage+2,trmage),mfexp(-yrfrct_srv1(i)*Z(i)(rcrage+2,trmage))),wt_srv1(i)(rcrage+2,trmage)));

    Espawnbio(i)= sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-0.21*Z(i))),wt_spawn(i)),0.5*mat));
//For retrospective comparison, use below AND the pre-1999 weight at age (in dat file)
//    Espawnbio(i)= sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-0.21*Z(i))),wt_spawn(i)),0.5*mat_old));
    }


// Do upper and lower accumulation in expected age composition

// Fishery
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

// Survey 1

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

// Survey 2
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

FUNCTION Projections

  dvariable sbio, sumbio;

  sbio = sumbio = 0.0;

  all_recruits.initialize();
  all_recruits = column(N,rcrage);

//Recruitments

 for (i=endyr+1;i<=endyr+5;i++)
    {
//  recruitment with bias correction
    recruit_proj(i)=mfexp(log_recr_proj(i)+(sigmasq_recr/2));
// for MCMC projections to get the prob < B20
//   recruit_proj(i)=mfexp(log_recr_proj(i));
// or just use average recruitment after 1977
// note that endyr-1 is the last year for mean
//    recruit_proj(i)=mean(recruit(1978,endyr-1));
    all_recruits(i) = mean(recruit((hcr_styr+rcrage),(endyr-1)));
    }

 for (i=endyr+1;i<=endyr+5;i++)
    {
    N_proj(i,rcrage)=recruit_proj(i);
    }

//Initialize the age composition

//  Standard projection
  for (j=rcrage;j<trmage;j++)
    {
    N_proj(endyr+1,j+1)=N(endyr,j)*mfexp(-Z(endyr,j));
    }
    N_proj(endyr+1,trmage)+=N(endyr,trmage)*mfexp(-Z(endyr,trmage));

// Set 2007 year class to mean
// note that endyr-1 is the last year for mean
//   for (j=rcrage;j<trmage;j++)
//   {
//   N_proj(endyr+1,j+1)=N(endyr,j)*mfexp(-Z(endyr,j));
//   }
//   N_proj(endyr+1,rcrage+1)=mean(recruit(1979,endyr-1))*mfexp(-Z(endyr,rcrage));
//   N_proj(endyr+1,trmage)+=N(endyr,trmage)*mfexp(-Z(endyr,trmage));

// Averaging window for selectivity
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

//Forward projections

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

// have to get Z twice
//  Tuning loop to get the spawning biomass adjustment right
   for (loop=1;loop<=10;loop++)
    {

    for (j=rcrage;j<=trmage;j++)
    {
    Z_proj(i,j)=(F_proj(i)*slctfsh_proj(j))+M;
    }
    sbio = sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-0.21*Z_proj(i))),wt_spawn_proj),0.5*mat));

  // MWD old version
//  Set the fishing mortality rate

    // F_proj(i)=Ftarget(i);
    // if (sbio < B40)
    // {
    // F_proj(i)=Ftarget(i)*(((sbio/B40)-0.05)/(1-0.05));
// // SSL control rule
// //   F_proj(i)=Ftarget(i)*(((sbio/B40)-0.2)/(1-0.2));
    // }
    // }

//  Set the fishing mortality rate
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

// Total mortality
  for (j=rcrage;j<=trmage;j++)
    {
    Z_proj(i,j)=(F_proj(i)*slctfsh_proj(j))+M;
    }

    sbio = sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-0.21*Z_proj(i))),wt_spawn_proj),0.5*mat));
    if (sbio < SB20)
    {
        sumbio = sum(elem_prod(N_proj(i)(rcrage+2,trmage),wt_pop_proj(rcrage+2,trmage)));
        cout << "in stock_assess::Projections:\tyear " << i << ", B40 " << B40 << ", SB40 " << SB40  << ", SB20 " << SB20 << ", Low SB " << sbio << ", sumbio " << sumbio << ", F_ABC " << F_proj(i) << endl;
    }

//  Numbers at age

  if(i<endyr+5)
  {
  for (j=rcrage;j<trmage;j++)
   {
   N_proj(i+1,j+1)=N_proj(i,j)*mfexp(-Z_proj(i,j));
   }
   N_proj(i+1,trmage)+=N_proj(i,trmage)*mfexp(-Z_proj(i,trmage));
   N(i+1) = N_proj(i+1);
  }
// Catches

  for (j=rcrage;j<=trmage;j++)
    {
    C_proj(i,j)=N_proj(i,j)*((F_proj(i)*slctfsh_proj(j))/Z_proj(i,j))*(1-mfexp(-Z_proj(i,j)));
//    Nsrv_proj(i,j)=q1_ek*slctsrv1(j)*N_proj(i,j)*mfexp(-yrfrct_srv1(endyr)*Z_proj(i,j));
    Nsrv_proj(i,j)=N_proj(i,j)*mfexp(-yrfrct_srv1(endyr)*Z_proj(i,j));
    }

//  Total catches and biomass

    Ecattot_proj(i) = 1000000*sum(elem_prod(C_proj(i),wt_fsh_proj));
// 3+ biomass
    Esumbio_proj(i)= N_proj(i)(rcrage+2,trmage)*wt_pop_proj(rcrage+2,trmage);
// Alternative: 2+ biomass
//    Esumbio_proj(i)= N_proj(i)(rcrage+1,trmage)*wt_pop_proj(rcrage+1,trmage);
    Exrate_proj(i)=Ecattot_proj(i)/(1000000*Esumbio_proj(i));
    Espawnbio_proj(i)= sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-0.21*Z_proj(i))),wt_spawn_proj),0.5*mat));
    Esrv_proj(i)= q1_ek*sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-yrfrct_srv1(endyr)*Z_proj(i))),slctsrv1),wt_srv_proj));
//    Esrv_proj(i)= sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-yrfrct_srv1(endyr)*Z_proj(i))),slctsrv1),wt_srv_proj));
    }


FUNCTION Objective_function

// Fishery likelihoods

//Total catch
  loglik(1) = -.5*norm2(elem_div((log(cattot)-log(Ecattot)),cattot_log_sd));
//Age composition
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

//Length composition
  for (i=1;i<=nyrslen_fsh;i++)
    {
    lllenp(i) = 0;
  for (j=1;j<=nbins1;j++)
    {
      lllenp(i) += multNlen_fsh(i)*(lenp(i,j)+o)*log((Elenp(fshlenyrs(i),j)+o)/(lenp(i,j)+o));
    }}
  loglik(3) = sum(lllenp);


// Survey 1 likelihoods
//Total biomass
  loglik(4) = -.5*norm2(elem_div(
       (log(indxsurv1_bs)-log(Eindxsurv1_bs(srvyrs1_bs))+square(indxsurv_log_sd1_bs)/2.),indxsurv_log_sd1_bs));
  loglik(4) += -.5*norm2(elem_div(
       (log(indxsurv1_ek)-log(Eindxsurv1_ek(srvyrs1_ek))+square(indxsurv_log_sd1_ek)/2.),indxsurv_log_sd1_ek));
  loglik(4) += -.5*norm2(elem_div(
       (log(indxsurv1_dy)-log(Eindxsurv1_dy(srvyrs1_dy))+square(indxsurv_log_sd1_dy)/2.),indxsurv_log_sd1_dy));
  RMSE_srv1_bs= sqrt(norm2(log(indxsurv1_bs)-log(Eindxsurv1_bs(srvyrs1_bs))+square(indxsurv_log_sd1_bs)/2.)/nyrs_srv1_bs);
  RMSE_srv1_ek= sqrt(norm2(log(indxsurv1_ek)-log(Eindxsurv1_ek(srvyrs1_ek))+square(indxsurv_log_sd1_ek)/2.)/nyrs_srv1_ek);
  RMSE_srv1_dy= sqrt(norm2(log(indxsurv1_dy)-log(Eindxsurv1_dy(srvyrs1_dy))+square(indxsurv_log_sd1_dy)/2.)/nyrs_srv1_dy);

//Age composition
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
//length composition
  for (i=1;i<=nyrslen_srv1;i++)
    {
    llsrvlenp1(i) = 0;
  for (j=1;j<=nbins3;j++)
    {
      llsrvlenp1(i) += multNlen_srv1(i)*(srvlenp1(i,j)+o)*log((Esrvlenp1(srv_lenyrs1(i),j)+o)/(srvlenp1(i,j)+o));
    }}
   loglik(6) = sum(llsrvlenp1);

// Survey 2 likelihoods

//Total biomass
  loglik(7) = -.5*norm2(elem_div(
       (log(indxsurv2)-log(Eindxsurv2(srvyrs2))+square(indxsurv_log_sd2)/2.),indxsurv_log_sd2));
  RMSE_srv2= sqrt(norm2(log(indxsurv2)-log(Eindxsurv2(srvyrs2))+square(indxsurv_log_sd2)/2.)/nyrs_srv2);

//Age composition
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
//length composition

  for (i=1;i<=nyrslen_srv2;i++)
    {
    llsrvlenp2(i) = 0;
  for (j=1;j<=nbins2;j++)
    {
      llsrvlenp2(i) += multNlen_srv2(i)*(srvlenp2(i,j)+o)*log((Esrvlenp2(srv_lenyrs2(i),j)+o)/(srvlenp2(i,j)+o));
    }}
   loglik(9) = sum(llsrvlenp2);

  loglik(10) = 0;

// Survey 3 likelihoods

//Total biomass
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
//length composition
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

//Constraints on recruitment
  loglik(18)= 0;
//  loglik(18)+= -0.5*square((mean_log_initN-mean_log_recruit)/0.3);
  loglik(18) += -0.5*square(dev_log_recruit(styr)/1.0);
  loglik(18) += -0.5*norm2(dev_log_recruit(styr+1,styr+7)/1.0);
//  loglik(18) += -0.5*norm2(dev_log_initN/1.0);
  loglik(18) += -0.5*norm2(dev_log_recruit(endyr-1,endyr)/1.0);
//  loglik(18) += -0.5*square(dev_log_recruit(endyr)/1.0);
//  Stronger constraint on endyear recruitment dev
//  loglik(18) += -0.5*square(dev_log_recruit(endyr)/0.25);


//Normal process error on selectivity deviations
//Note rwlk_sd(styr,endyr-1)

  loglik(19)  = -0.5*norm2(elem_div(first_difference(slp1_fsh_dev),rwlk_sd_short(styr,endyr_fsh_dev-1)));
  loglik(19) += -0.5*norm2(elem_div(first_difference(inf1_fsh_dev),4.0*rwlk_sd_short(styr,endyr_fsh_dev-1)));
  loglik(19) += -0.5*norm2(elem_div(first_difference(slp2_fsh_dev),rwlk_sd_short(styr,endyr_fsh_dev-1)));
  loglik(19) += -0.5*norm2(elem_div(first_difference(inf2_fsh_dev),4.0*rwlk_sd_short(styr,endyr_fsh_dev-1)));

// Recruitment in projection mode

  if(last_phase())
  {
  loglik(20) =  -(1/(2.0*sigmasq_recr))*norm2(log_recr_proj - log_mean_recr_proj);
  }
  else
  {
  loglik(20)=0;
  }
  loglik(21)= 0;

//Likelihood for the vessel comparison experiment between the oscar dyson and the miller freeman
   loglik(22)= -(1/(2.0*(square(0.0244)+square(0.000001))))*square(log_q1_dy- log_q1_ek-.124);
//   loglik(22)= -(1/(2.0*square(0.001)))*square(log_q1_dy- log_q1_ek-.124);
//   loglik(22)= -(1/(2.0*square(0.1)))*square(log_q1_dy- log_q1_ek-.124);
//   loglik(22)= 0;

  objfun = -sum(loglik);

// Variable to do a likelihood profile over
    var_prof=Esumbio(endyr);
//  var_prof=mean(recruit(1979,endyr-2));
//    var_prof=log_q2;


FUNCTION dvariable solve_for_fishing_mortality(dvariable catch_level, int naa_yr)

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


FUNCTION void calculate_curr_bio_ref_points(dvar_vector sel, int curr_yr)

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
            cout << "in stock_assess::calculate_curr_bio_ref_points:\tSB <= SB20% in year " << next_yr << ":\tSB20 "<< spbio_floor << "\tSpawning biomass " << Espawnbio(next_yr) << "\tTotal biomass " << Esumbio(next_yr) << endl;
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

// Jim Ianelli's code for calculating Fxx%
// START

FUNCTION dvariable get_spr_rates(dvariable spr_percent, dvar_vector sel, int curr_yr)

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

FUNCTION dvariable spr_unfished(int curr_yr)

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

FUNCTION dvariable spr_ratio(dvariable trial_F, dvar_vector sel, int curr_yr)

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

// end of Jim Ianelli's code


FUNCTION MCMC_output
  if(mceval_phase())
  {
// Dumping data to mceval.dat
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


RUNTIME_SECTION

  convergence_criteria 1.e0, 1.e-1, 1.e-4, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7
  maximum_function_evaluations 1000, 1000, 1000, 1000

TOP_OF_MAIN_SECTION
 arrmblsize = 3000000;
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
 gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);


REPORT_SECTION

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
//  report << log_slp1_srv1 << endl;
//  report << inf1_srv1 << endl;
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

