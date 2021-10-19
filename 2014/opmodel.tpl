// Gulf of Alaska walleye pollock operating model
// temporal/seasonal components
// no spatial components
// ZTA
// Version 0.AF, January 2007
// using ADMB v7.1.1, with MinGW g++ 3.4.4 at AFSC and MS VC++ 6.0 at UW
//
// DATA_SECTION comes from MWD's 2011 GOA pollock stock assessment model, pk11_1.tpl
//
// 2008-02-14, Chapters 5 and 6, base scenario
//
// 2011-08-22, update with data through 2010
//
// 2011-12-02, update with data through 2011
//
// 2014-03-31, update with data through 2013 and change rcrage from 2 to 1

DATA_SECTION

// BEGIN:  MWD GOA pollock stock assessment data

!! ad_comm::change_datafile_name("pk13_3.dat");

  init_int styr                                  // Starting year for population model
  init_int endyr                                 // Ending year for population model
  init_int base_endyr                            // Ending year for base model
  init_int hcr_styr                              // Starting year for recruitment calculations for the harvest control rule
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
  init_vector indxsurv1_ek(1,nyrs_srv1_ek)          // Survey index
  init_vector indxsurv_log_sd1_ek(1,nyrs_srv1_ek)   // Survey index (cv) = sdev of log(indxsurv)
//Dyson
  init_int nyrs_srv1_dy                             // Number of survey biomass estimates
  init_ivector srvyrs1_dy(1,nyrs_srv1_dy)           // Years in which surveys occured
  init_vector indxsurv1_dy(1,nyrs_srv1_dy)          // Survey index
  init_vector indxsurv_log_sd1_dy(1,nyrs_srv1_dy)   // Survey index (cv) = sdev of log(indxsurv)

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
// Variance for log recr
  init_number sigmasq_recr

// END:  MWD GOA pollock stock assessment data

  int styr_fsh_sel_dev                          // updated first value for fsh sel devs

  int iii                                        // Index for year
  int jjj                                        // Index for age

  number o                                       // A small number


  // new section for operating model parameters
!! ad_comm::change_datafile_name("opmodel.dat");

  init_int st_age                           // starting age for age classes to keep track of
  init_int end_age                          // ending age for age classes to keep track of
  init_int om_hcr_styr                      // starting year in the OPERATING MODEL for recruitment calculations for the harvest control rule
  init_int om_rec_avg_styr                  // starting year for the OPERATING MODEL for calculating average recruitment and std dev (used for generating future recruitment)

  int nyears                                // number of years to keep track of
  int nages                                 // number of ages to keep track of

  init_int nsel_fsh                         // number of fishery selectivity curves
  init_vector fsh_frac(1,nsel_fsh)          // fraction of the year when fishing occurs
  init_int nsel_srv                         // number of survey selectivity curves
  init_vector srv_frac(1,nsel_srv)          // fraction of the year when survey occurs
  init_number sp_frac                       // fraction of the year when spawning occurs

  int nsel_tot                              // total number of selectivity curves
  int nyears_simple                         // number of years (starting at styr) of continuous fishing
  int nyears_complex                        // number of years (until endyr) of seasonal fishing

  init_int styr_fsh_sel                     // year that multiple fishing seasons start
  init_matrix frac_catch(styr_fsh_sel,endyr,1,nsel_fsh) // catch by year and season (fraction of year)

  init_int nyrs_fsh_paa_all                 // number of years of expanded fishery c-a-a data (2004 SAFE)
  init_ivector yrs_fsh_paa_all(1,nyrs_fsh_paa_all)      // years of expanded fishery c-a-a data
  init_vector multN_fsh_paa_all(1,nyrs_fsh_paa_all)     // effective sample size
  init_matrix fsh_paa_all(1,nyrs_fsh_paa_all,st_age,end_age)    // expanded fishery c-a-a data

  init_int nyrs_srv_1_paa_all               // number of years of expanded EIT survey p-a-a data (2004 SAFE)
  init_ivector yrs_srv_1_paa_all(1,nyrs_srv_1_paa_all)  // years of expanded EIT survey p-a-a data
  init_vector multN_srv_1_paa_all(1,nyrs_srv_1_paa_all) // effective sample size
  init_matrix srv_1_paa_all(1,nyrs_srv_1_paa_all,st_age,end_age)    // expanded EIT survey p-a-a data

  init_int nyrs_srv_2_paa_all               // number of years of expanded EIT survey p-a-a data (2004 SAFE)
  init_ivector yrs_srv_2_paa_all(1,nyrs_srv_2_paa_all)  // years of expanded EIT survey p-a-a data
  init_vector multN_srv_2_paa_all(1,nyrs_srv_2_paa_all) // effective sample size
  init_matrix srv_2_paa_all(1,nyrs_srv_2_paa_all,st_age,end_age)    // expanded EIT survey p-a-a data

  init_vector M(st_age,end_age)             // initial values for natural mortality at age

  int first_rec_year                        // first year that recruits will be calculated
  int last_rec_year                         // last year that recruits will be calculated

  init_int s_r_relat                        // Stock-Recruit relationship: 1 - B-H, 2 - R, 3 - S, 4 - avg R
  int h_phase                               // don't estimate h for avg R S-R relationship (R0 = avg R)

  init_int nsel_lengths                     // number of lengths in selectivity curves
  init_vector sel_lengths(1,nsel_lengths)   // midpoint of selectivity length bins
  init_int phase_fsh_sel                    // phase to turn on fishery selectivity curves
  init_ivector phase_q(1,nsel_srv)          // phase to turn on survey catchability estimation
  init_ivector phase_srv_sel_age1(1,nsel_srv)   // phase to turn on survey selectivity value for age-st_age fish
  init_ivector phase_srv_sel_a1(1,nsel_srv) // phase to turn on survey selectivity curves (ascending part)
  init_ivector phase_srv_sel_b1(1,nsel_srv) // phase to turn on survey selectivity curves (ascending part)
  init_ivector phase_srv_sel_a2(1,nsel_srv) // phase to turn on survey selectivity curves (descending part)
  init_ivector phase_srv_sel_b2(1,nsel_srv) // phase to turn on survey selectivity curves (descending part)
  init_matrix age_trans_all(st_age,end_age,st_age,end_age)  // ageing error matrix for expanded proportions-at-age data
  init_matrix age_len_trans(st_age,end_age,1,nsel_lengths)  // age-length transition matrix
  vector maa_old(st_age,end_age)            // old maturity-at-age matrix
  vector maa(st_age,end_age)                // maturity-at-age matrix

  init_vector waa_pop(st_age,end_age)       // average population weight-at-age
  init_vector waa_fsh(st_age,end_age)       // average fishery weight-at-age
  init_matrix waa_srv(1,nsel_srv,st_age,end_age)    // average weight-at-age for surveys

  init_int nyrs_proj                        // number of years to do projections
  init_matrix init_M_proj(endyr+1,endyr+nyrs_proj,st_age,end_age)   // M-at-age values to use in projections
  init_int nsubsample                       // number of MCMC iterations to skip
  int mcmc_iter                             // number of MCMC iterations

  init_int complex_flag                     // 0 - determinstic, simple model; 1 - recruitment error only; 2 - recruitment and survey obs error; 3 - full model
  init_int future_rec_flag                  // 0 - generate future recruitment randomly using avg level of recruitment and sigmaR; 1 - future recruitment is historical recruitment sampled with replacement
  init_int debug_flag                       // 1 - debug on; 0 - debug off
  init_int catch_flag                       // 0 - no catch in future projections; 1 - calculate catch
  init_int BRP_M_flag                       // 0 - use annual M-at-age for calculating BRPs; 1 - use average M-at-age since 1977+st_age for calculating BRPs


!! if (catch_flag == 3)
!! {
  // new section for future catches
!!     ad_comm::change_datafile_name("future_catch.dat");

  init_matrix proj_catch_array(1,100,endyr+1,endyr+nyrs_proj);   // values for annual catch in the future

!! }


  number sigma_f                            // sigmaF (sigmaR for equilibrium numbers in first year)
  number sigma_r                            // sigmaR (recruitment deviance)

  number Tier3_alpha                        // alpha for the NPFMC Tier 3 decision rule, 0.05

  number mult_factor                        // ABC when SB < SB20% is mult_factor * age_3_biomass, 0.001 * 0.001 => metric tonnes

  number conv_factor                        // 1e9

  number R0_mean                            // mean for prior on R0
  number R0_stddev                          // std dev for prior on R0

  !!#define MAX_IDF_LINES    4096           // maximum number of lines in pk13_proj.dat
  !!#define MAX_IDF_LINE_LEN 4096           // maximum line length of lines in pk13_proj.dat
  !!#define MAX_BISECT_ITER  128            // maximum number of iterations for the bisection
  !!#define BISECT_TOL       0.000001       // maximum tolerance for bisection


  !!cout << "Data check" << endl;
  !!cout << "nyrs_srv2\t" << nyrs_srv2 << endl;
  !!cout << "srv_frac\t" << srv_frac << endl;
  !!cout << "phase_q\t" << phase_q << endl;
  !!cout << "nyrs_proj\t" << nyrs_proj << endl;



 LOCAL_CALCS

  nsel_tot       = nsel_fsh + nsel_srv;
  nyears_simple  = styr_fsh_sel - styr;
  nyears_complex = endyr - styr_fsh_sel + 1;

  nyears         = endyr - styr + 1;
  nages          = end_age - st_age + 1;

  first_rec_year = styr;
  last_rec_year  = endyr;

  // styr_fsh_sel_dev = styr + 11;         // HARDCODED - skip the years that MWD has random walk SD as 0.001
  // rwlk_sd *= 3.0;                       // HARDCODED - test
  // rwlk_sd(styr_fsh_sel_dev) *= 2.0;     // HARDCODED - increase the SD for the first year used
  styr_fsh_sel_dev = styr;                  // stop messing about with the rwlk_sd vector

  if (s_r_relat == 1 || s_r_relat == 2 || s_r_relat == 3)
  {
      h_phase = 2;
  }
  else
  {
      h_phase = -2;
  }

  // map MWD maturity at age to (larger) maturity at age
  maa_old.initialize();
  maa.initialize();
  maa_old = maa = 0.0;
  if (trmage <= end_age)
  {
      maa_old(rcrage,trmage) = mat_old(rcrage,trmage);
      maa(rcrage,trmage)     = mat(rcrage,trmage);

      if (trmage < end_age)
      {
          for (int a = (trmage + 1); a <= end_age; a++)
          {
              maa_old(a) = 1.0;
              maa(a)     = 1.0;
          }
      }
  }

  sigma_f = 1.0;
  sigma_r = 1.0;

  Tier3_alpha = 0.05;

  mult_factor = 0.000001;

  conv_factor = 1000000000.;    // convert from kg to millions of metric tonnes

  R0_mean   = 1.0446e+009;
  R0_stddev = 0.1 * R0_mean;    // AEP-specified CV of 0.1

  o = 0.00001;

  mcmc_iter  = 0;


 END_CALCS


INITIALIZATION_SECTION

  // using a PIN file instead


PARAMETER_SECTION

  init_bounded_number log_R0(5.0,35.0,1);
  init_bounded_number log_h(-1.6,1.6,h_phase);
  init_bounded_number log_q(-10.0,10.0,-1);
  init_bounded_vector log_q_fsh(1,nsel_fsh,-10.0,10.0,-1);
  init_bounded_number log_q_srv_bs(-10.0,10.0,5);
  init_bounded_number log_q_srv_ek(-10.0,10.0,5);
  init_bounded_number_vector log_q_srv(1,nsel_srv,-10.0,10.0,phase_q);
  init_bounded_number log_initR(5.0,35.0,1);
  init_bounded_number mean_log_Fmort(-10.0,10.0,1);
  init_bounded_dev_vector log_Fmort_dev(styr,endyr,-10.0,10.0,2);
  init_bounded_dev_vector log_rec_dev(first_rec_year,last_rec_year,-15.0,15.0,3);
  init_bounded_dev_vector log_init_dev(st_age+1,end_age,-15.0,15.0,-2);

  // comparing MWD's code:  slope == b, inflection == a
  init_bounded_number log_fsh_sel_cont_a1(0.0,2.0,phase_fsh_sel);       // in MWD model, bounds are 1 and 5; starting point is 4
  init_bounded_number log_fsh_sel_cont_b1(-5.0,5.0,phase_fsh_sel);      // in MWD model, bounds are -5 and 5; starting point is 1
  init_bounded_number log_fsh_sel_cont_a2(1.0,3.0,phase_fsh_sel);       // in MWD model, bounds are 7 and 20; starting point is 7
  init_bounded_number log_fsh_sel_cont_b2(-5.0,5.0,phase_fsh_sel);      // in MWD model, bounds are -5 and 5; starting point is 1

  init_bounded_vector log_fsh_sel_a1(1,nsel_fsh,0.0,5.0,-phase_fsh_sel);
  init_bounded_vector log_fsh_sel_b1(1,nsel_fsh,-5.0,5.0,-phase_fsh_sel);
  init_bounded_vector log_fsh_sel_a2(1,nsel_fsh,1.0,5.0,-phase_fsh_sel);
  init_bounded_vector log_fsh_sel_b2(1,nsel_fsh,-5.0,5.0,-phase_fsh_sel);

  init_bounded_number_vector log_srv_sel_age1(1,nsel_srv,-15.0,1.0,phase_srv_sel_age1); // in MWD model, bounds are 0 and 2

  init_bounded_number_vector log_srv_sel_a1(1,nsel_srv,0.0,5.0,phase_srv_sel_a1);   // in MWD model, bounds are (X, 1, X, 1, X, 1) and (X, 50, X, 20, X, 100); starting points are (X, 7, X, 5, X, 5)
  init_bounded_number_vector log_srv_sel_b1(1,nsel_srv,-5.0,5.0,phase_srv_sel_b1);  // in MWD model, bounds are -5 and 5; starting points are (X, 0, X, 0, X, 0)
  init_bounded_number_vector log_srv_sel_a2(1,nsel_srv,1.0,5.0,phase_srv_sel_a2);   // in MWD model, bounds are (7, 5, X, X, X, 5) and (20, 10, X, X, X, 100); starting points are (9, 8, X, X, X, 8)
  init_bounded_number_vector log_srv_sel_b2(1,nsel_srv,-5.0,5.0,phase_srv_sel_b2);  // in MWD model, bounds are -5 and 5, except for srv 6 (-20,5); starting points are (1, 0, X, X, X, -10)

  // init_number_vector class is not defined in all versions of ADMB
  // init_vector log_srv_sel_a1(1,nsel_srv,-phase_fsh_sel);
  // init_vector log_srv_sel_b1(1,nsel_srv,-phase_fsh_sel);
  // init_vector log_srv_sel_a2(1,nsel_srv,-phase_fsh_sel);
  // init_vector log_srv_sel_b2(1,nsel_srv,-phase_fsh_sel);

  init_bounded_dev_vector fsh_sel_cont_a1_dev(styr_fsh_sel_dev,endyr,-5.0,5.0,phase_fsh_sel);
  init_bounded_dev_vector fsh_sel_cont_b1_dev(styr_fsh_sel_dev,endyr,-5.0,5.0,phase_fsh_sel);
  init_bounded_dev_vector fsh_sel_cont_a2_dev(styr_fsh_sel_dev,endyr,-5.0,5.0,phase_fsh_sel);
  init_bounded_dev_vector fsh_sel_cont_b2_dev(styr_fsh_sel_dev,endyr,-5.0,5.0,phase_fsh_sel);

  matrix N(styr,endyr+1,st_age,end_age);
  // matrix NL(styr,endyr+1,1,nsel_lengths);
  matrix C(styr,endyr,st_age,end_age);
  matrix F(styr,endyr,st_age,end_age);
  matrix Slen(1,nsel_tot+1,1,nsel_lengths);
  matrix Sage(1,nsel_tot+1,st_age,end_age);
  matrix Sage_fsh(styr,endyr,st_age,end_age);
  matrix Z(styr,endyr,st_age,end_age);
  matrix expZ(styr,endyr,st_age,end_age);
  matrix expZsp(styr,endyr,st_age,end_age);
  3darray expZfrac(1,nsel_srv,styr,endyr,st_age,end_age);

  vector initN(st_age,end_age);
  vector Fmort(styr,endyr);
  vector spawn_biomass(styr-1,endyr+1);
  vector total_biomass(styr-1,endyr+1);
  vector age_3_plus_biomass(styr-1,endyr+1);
  vector estsrvbio_bs(styr,endyr);
  vector estsrvbio_ek(styr,endyr);
  matrix estsrvbio(1,nsel_srv,styr,endyr);

  number fsh_sel_cont_a1;
  number fsh_sel_cont_b1;
  number fsh_sel_cont_a2;
  number fsh_sel_cont_b2;

  vector fsh_sel_a1(1,nsel_fsh);
  vector fsh_sel_b1(1,nsel_fsh);
  vector fsh_sel_a2(1,nsel_fsh);
  vector fsh_sel_b2(1,nsel_fsh);

  vector srv_sel_age1(1,nsel_srv);
  vector srv_sel_a1(1,nsel_srv);
  vector srv_sel_b1(1,nsel_srv);
  vector srv_sel_a2(1,nsel_srv);
  vector srv_sel_b2(1,nsel_srv);

  number R0;
  number h;
  number steepness;
  number q;
  vector q_fsh(1,nsel_fsh);
  number q_srv_bs;
  number q_srv_ek;
  vector q_srv(1,nsel_srv);
  number initR;

  // extra data not in MWD's model
  vector pred_catch(styr,endyr);
  matrix pred_fsh_paa_all(styr,endyr,st_age,end_age);
  matrix pred_fsh_paa(styr,endyr,rcrage,trmage);
  3darray pred_srv_paa_all(1,nsel_srv,styr,endyr,st_age,end_age);
  3darray pred_srv_paa(1,nsel_srv,styr,endyr,rcrage,trmage);

  matrix pred_fsh_pal(1,nyrslen_fsh,1,nbins1);
  matrix pred_srv_1_pal(1,nyrslen_srv1,1,nbins3);
  matrix pred_srv_2_pal(1,nyrslen_srv2,1,nbins2);
  matrix pred_srv_3_pal(1,nyrslen_srv3,1,nbins2);

  // projection parameters
  matrix N_proj(endyr+1,endyr+nyrs_proj+1,st_age,end_age);
  matrix C_proj(endyr+1,endyr+nyrs_proj,st_age,end_age);
  vector F_proj(endyr+1,endyr+nyrs_proj);
  matrix M_proj(endyr+1,endyr+nyrs_proj,st_age,end_age);
  matrix Z_proj(endyr+1,endyr+nyrs_proj,st_age,end_age);
  matrix expZ_proj(endyr+1,endyr+nyrs_proj,st_age,end_age);
  matrix expZsp_proj(endyr+1,endyr+nyrs_proj,st_age,end_age);
  3darray expZfrac_proj(1,nsel_srv,endyr+1,endyr+nyrs_proj,st_age,end_age);
  matrix fsh_paa_all_proj(endyr,endyr+nyrs_proj,st_age,end_age);
  3darray srv_paa_all_proj(1,nsel_srv,endyr+1,endyr+nyrs_proj,st_age,end_age);
  vector catch_proj(endyr+1,endyr+nyrs_proj);
  vector spawn_biomass_proj(endyr+1,endyr+nyrs_proj);
  vector total_biomass_proj(endyr+1,endyr+nyrs_proj);
  vector age_3_plus_biomass_proj(endyr+1,endyr+nyrs_proj);
  matrix estsrvbio_proj(1,nsel_srv,endyr+1,endyr+nyrs_proj);
  3darray estsrvN_proj(1,nsel_srv,endyr+1,endyr+nyrs_proj,st_age,end_age);
  vector estsrvbio_proj_2_sd(endyr+1,endyr+nyrs_proj);
  vector estsrvbio_proj_2_multN(endyr+1,endyr+nyrs_proj);

  vector estFOFL_proj(endyr+1,endyr+nyrs_proj);
  vector estFABC_proj(endyr+1,endyr+nyrs_proj);
  vector estABC_proj(endyr+1,endyr+nyrs_proj);
  matrix estABCaa_proj(endyr+1,endyr+nyrs_proj,st_age,end_age);

  number mcmc_avg_rec;
  number mcmc_avg_log_rec;
  number mcmc_std_dev_avg_log_rec;
  number mcmc_CV_avg_log_rec;
  number mcmc_avg_rec_years;
  number mcmc_max_rec;

  // initialize the random number generator used for the projections
  !!CLASS random_number_generator mcmc_rng(1);

  // parameters for decision rule
  vector F100(endyr,endyr+nyrs_proj);
  vector SB100(endyr,endyr+nyrs_proj);
  vector SBtarget(endyr,endyr+nyrs_proj);
  vector F40(endyr,endyr+nyrs_proj);
  vector SB40(endyr,endyr+nyrs_proj);
  vector F35(endyr,endyr+nyrs_proj);
  vector SB35(endyr,endyr+nyrs_proj);
  vector F20(endyr,endyr+nyrs_proj);
  vector SB20(endyr,endyr+nyrs_proj);
  vector Sage_fsh_avg(st_age,end_age);
  vector M_at_age(st_age,end_age);
  number F_ABC;
  number ABC;
  number F_OFL;
  number OFL;
  number phi0;
  number avgR;
  number avgR_CV;
  number SB0;
  number SBcurr;

  // calculated value for B0 over a moving window of 25 years of avg rec level
  // see 2005 GOA pollock SAFE, pg 67 (pg 27 in PDF)
  vector B0_mw(endyr+1,endyr+nyrs_proj);

  vector f(1,20);

  objective_function_value obj_fun;

  // sdreport_number sd_q_srv_bs;
  // sdreport_vector sd_q_srv(1,nsel_srv);
  sdreport_vector sd_recruits_1(styr,endyr+1);
  sdreport_number sd_avg_rec_1;
  sdreport_number sd_avg_rec_1_alt;
  sdreport_vector sd_recruits_2(styr,endyr+1);
  sdreport_number sd_avg_rec_2;
  sdreport_number sd_avg_rec_2_alt;
  sdreport_vector sd_total_biomass(styr,endyr);
  sdreport_vector sd_spawn_biomass(styr,endyr);
  sdreport_vector sd_age_3_plus_biomass(styr,endyr);
  // sdreport_vector sd_estsrvbio_bs(styr,endyr);
  // sdreport_matrix sd_estsrvbio(1,nsel_srv,styr,endyr);

  likeprof_number var_prof


PRELIMINARY_CALCS_SECTION   // from MWD's model

// Do upper and lower accumulations in age composition data

// Fishery
  for (iii=1;iii<=nyrs_fsh;iii++)
    {
  for (jjj=rcrage;jjj<=trmage;jjj++)
    {

    if(jjj<ac_yng_fsh(iii))
      {
      catp(iii,ac_yng_fsh(iii)) += catp(iii,jjj);
      catp(iii,jjj) = 0;
      }
    if(jjj>ac_old_fsh(iii))
      {
      catp(iii,ac_old_fsh(iii)) += catp(iii,jjj);
      catp(iii,jjj) = 0;
      }
    }}

// Survey 1
  for (iii=1;iii<=nyrsac_srv1;iii++)
    {
  for (jjj=rcrage;jjj<=trmage;jjj++)
    {
    if(jjj<ac_yng_srv1(iii))
      {
      srvp1(iii,ac_yng_srv1(iii)) += srvp1(iii,jjj);
      srvp1(iii,jjj) = 0;
      }
    if(jjj>ac_old_srv1(iii))
      {
      srvp1(iii,ac_old_srv1(iii)) += srvp1(iii,jjj);
      srvp1(iii,jjj) = 0;
      }
    }}

// Survey 2
  for (iii=1;iii<=nyrsac_srv2;iii++)
    {
  for (jjj=rcrage;jjj<=trmage;jjj++)
    {
    if(jjj<ac_yng_srv2(iii))
      {
      srvp2(iii,ac_yng_srv2(iii)) += srvp2(iii,jjj);
      srvp2(iii,jjj) = 0;
      }
    if(jjj>ac_old_srv2(iii))
      {
      srvp2(iii,ac_old_srv2(iii)) += srvp2(iii,jjj);
      srvp2(iii,jjj) = 0;
      }
    }}

  var_prof.set_stepnumber(30);
//  var_prof.set_stepnumber(10);
  var_prof.set_stepsize(0.1);


PROCEDURE_SECTION

    calculate_selectivity();

    calculate_mortality();

    calculate_n_a_a();

    calculate_c_a_a();

    calculate_length_comps();

    evaluate_objective_function();

    calculate_projections();


FUNCTION calculate_selectivity

    dvariable maxSel;
    int a,i,j;

    fsh_sel_cont_a1 = mfexp(log_fsh_sel_cont_a1);
    fsh_sel_cont_b1 = mfexp(log_fsh_sel_cont_b1);
    fsh_sel_cont_a2 = mfexp(log_fsh_sel_cont_a2);
    fsh_sel_cont_b2 = mfexp(log_fsh_sel_cont_b2);

    fsh_sel_a1      = mfexp(log_fsh_sel_a1);
    fsh_sel_b1      = mfexp(log_fsh_sel_b1);
    fsh_sel_a2      = mfexp(log_fsh_sel_a2);
    fsh_sel_b2      = mfexp(log_fsh_sel_b2);

    srv_sel_age1    = mfexp(log_srv_sel_age1);
    srv_sel_a1      = mfexp(log_srv_sel_a1);
    srv_sel_b1      = mfexp(log_srv_sel_b1);
    srv_sel_a2      = mfexp(log_srv_sel_a2);
    srv_sel_b2      = mfexp(log_srv_sel_b2);

    Slen.initialize();
    Sage.initialize();

    // if (debug_flag) cout << fsh_sel_cont_b2 << " " << fsh_sel_b2 << " " << srv_sel_b2 << endl;
    // if (debug_flag) cout << endl;

    // selectivity curves 1 to nsel_fsh are for the fisheries (seasons in later years)
    // selectivity curves nsel_fsh + 1 to nsel_fsh + nsel_srv are for the surveys
    // selectivity curve nsel_tot (nsel_fsh + nsel_srv) + 1 is for continuous fishing
    // for (l = 1; l <= nsel_lengths; l++)
    // {
    //     // all curves are double logistic, courtesy of Dorn and Methot 1990

    //     for (i = 1; i <= nsel_fsh; i++)
    //     {
    //         Slen(i,l) = (1.0 / (1.0 + mfexp(-fsh_sel_b1(i)*(sel_lengths(l) - fsh_sel_a1(i)))))*(1.0 - (1.0 / (1.0 + mfexp(-fsh_sel_b2(i)*(sel_lengths(l) - fsh_sel_a2(i))))));
    //     }

    //     for (i = 1; i <= nsel_srv; i++)
    //     {
    //         Slen(i+nsel_fsh,l) = (1.0 / (1.0 + mfexp(-srv_sel_b1(i)*(sel_lengths(l) - srv_sel_a1(i)))))*(1.0 - (1.0 / (1.0 + mfexp(-srv_sel_b2(i)*(sel_lengths(l) - srv_sel_a2(i))))));
    //     }

    //     Slen(nsel_tot+1,l) = (1.0 / (1.0 + mfexp(-fsh_sel_cont_b1*(sel_lengths(l) - fsh_sel_cont_a1))))*(1.0 - (1.0 / (1.0 + mfexp(-fsh_sel_cont_b2*(sel_lengths(l) - fsh_sel_cont_a2)))));
    // }

    // selectivity curves 1 to nsel_fsh are for the fisheries (seasons in later years)
    // selectivity curves nsel_fsh + 1 to nsel_fsh + nsel_srv are for the surveys
    // selectivity curve nsel_tot (nsel_fsh + nsel_srv) + 1 is for continuous fishing
    for (a = st_age; a <= end_age; a++)
    {
        // all curves are double logistic, courtesy of Dorn and Methot 1990

        for (i = 1; i <= nsel_fsh; i++)
        {
            Sage(i,a) = (1.0 / (1.0 + mfexp(-fsh_sel_b1(i)*(double(a) - fsh_sel_a1(i)))))*(1.0 - (1.0 / (1.0 + mfexp(-fsh_sel_b2(i)*(double(a) - fsh_sel_a2(i))))));
        }

        for (i = 1; i <= nsel_srv; i++)
        {
            Sage(i+nsel_fsh,a) = (1.0 / (1.0 + mfexp(-srv_sel_b1(i)*(double(a) - srv_sel_a1(i)))))*(1.0 - (1.0 / (1.0 + mfexp(-srv_sel_b2(i)*(double(a) - srv_sel_a2(i))))));
        }

        Sage(nsel_tot+1,a) = (1.0 / (1.0 + mfexp(-fsh_sel_cont_b1*(double(a) - fsh_sel_cont_a1))))*(1.0 - (1.0 / (1.0 + mfexp(-fsh_sel_cont_b2*(double(a) - fsh_sel_cont_a2)))));
    }

    // set selectivity for upper ages to selectivity at trmage
    // MAY CHANGE LATER WHEN DATA IS AVAILABLE TO AGE 15
    // for version 0.16 and above:  extended p-a-a data (ages 1 to 15) is
    // available for the fishery and surveys 1 and 2 (2004 GOA pollock SAFE)
    if (trmage < end_age)
    {
        for (j = (trmage + 1); j <= end_age; j++)
        {
            i = 3 + nsel_fsh;   // survey 3
            Sage(i,j) = Sage(i,trmage);
        }
    }

    // match MWD's model by having a separate value for age-st_age fish
    for (i = 1; i <= nsel_srv; i++)
    {
        // include them in the selectivity curve ONLY if they will be estimated
        if (phase_srv_sel_age1(i) > 0)
        {
            Sage(i+nsel_fsh,st_age) = srv_sel_age1(i);
        }
    }

    // rescale age selectivity curves
    for (i = 1; i <= nsel_tot+1; i++)
    {
        maxSel = max(Sage(i));
        Sage(i) /= maxSel;
    }

    // match MWD's model by having a separate value for age-st_age fish
    // for (i = 1; i <= nsel_srv; i++)
    // {
        // include them in the selectivity curve ONLY if they will be estimated
    //     if (phase_srv_sel_age1(i) > 0)
    //     {
    //         Sage(i+nsel_fsh,st_age) = srv_sel_age1(i);
    //     }
    // }

    // convert from selectivity-at-length to selectivity-at-age
    for (i = 1; i <= nsel_tot+1; i++)
    {
        Slen(i) = Sage(i) * age_len_trans;
    }

    // rescale length selectivity curves
    for (i = 1; i <= nsel_tot+1; i++)
    {
        maxSel = max(Slen(i));
        Slen(i) /= maxSel;
    }

    if (debug_flag) cout << "end of calc_sel" << endl;


FUNCTION calculate_mortality

    dvariable maxSel;
    int y,a,i;

    Sage_fsh.initialize();
    Fmort.initialize();

    // annual fishing mortality
    Fmort = mfexp(mean_log_Fmort + log_Fmort_dev);

    // covers the period for changes in fishery selectivity
    for (y = styr_fsh_sel_dev; y <= endyr; y++)
    {
        // annual fishery selectivity
        for (a = st_age; a <= end_age; a++)
        {
            Sage_fsh(y,a) = (1.0 / (1.0 + mfexp(-(fsh_sel_cont_b1 * mfexp(fsh_sel_cont_b1_dev(y))) * (double(a) - (fsh_sel_cont_a1 + fsh_sel_cont_a1_dev(y)))))) * (1.0 - (1.0 / (1.0 + mfexp(-(fsh_sel_cont_b2 * mfexp(fsh_sel_cont_b2_dev(y))) * (double(a) - (fsh_sel_cont_a2 + fsh_sel_cont_a2_dev(y)))))));
        }

        // normalize selectivity
        maxSel = max(Sage_fsh(y));
        Sage_fsh(y) /= maxSel;
    }

    // fishery selectivity is constant over the period where MWD has very low SD on random walk
    for (y = styr; y < styr_fsh_sel_dev; y++)
    {
        Sage_fsh(y) = Sage_fsh(styr_fsh_sel_dev);
    }

    for (y = styr; y <= endyr; y++)
    {
        // ages covered by selectivity
        for (a = st_age; a <= end_age; a++)
        {
            F(y,a) = Sage_fsh(y,a) * Fmort(y);
            Z(y,a) = M(a) + F(y,a);
        }
    }

    expZ   = mfexp(-Z);

    // survival with fishing until spawning
    expZsp = mfexp(-sp_frac * Z);

    for (y = styr; y <= endyr; y++)
    {
        // this assumes continuous fishing - NEED TO FIX
        expZfrac(1,y) = mfexp(-yrfrct_srv1(y) * Z(y));
        expZfrac(2,y) = mfexp(-yrfrct_srv2(y) * Z(y));
        expZfrac(3,y) = mfexp(-yrfrct_srv3(y) * Z(y));
    }

    if (debug_flag) cout << "end of calc_mort" << endl;


FUNCTION calculate_n_a_a

    dvariable totn, wt_a, wt_sp_a, spbio_tmp, Ntmp, srvbio_tmp;
    int y,a,i,j,j_max;
    int acc_age;

    R0       = mfexp(log_R0);
    h        = mfexp(log_h);
    q        = mfexp(log_q);
    q_fsh    = mfexp(log_q_fsh);
    q_srv_bs = mfexp(log_q_srv_bs);
    q_srv_ek = mfexp(log_q_srv_ek);
    q_srv    = mfexp(log_q_srv);
    initR    = mfexp(log_initR);

    // equilibrium age structure
    initN(st_age) = initR * mfexp(log_rec_dev(styr));
    for (a = (st_age+1); a <= end_age; a++)
    {
        initN(a) = initN(a-1) * mfexp(-M(a-1));
    }
    initN(end_age) /= (1.0 - mfexp(-M(end_age)));

    // put in first year deviances but not for the recruits (age 1)
    for (a = (st_age+1); a <= end_age; a++)
    {
        initN(a) *= mfexp(log_init_dev(a));
    }

    // calculate the equilibrium spawning and total biomass
    spawn_biomass.initialize();
    total_biomass.initialize();
    age_3_plus_biomass.initialize();

    for (a = (st_age+1); a <= end_age; a++)
    {
        if (a < rcrage || a > trmage)
        {
            if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
            {
                wt_a    = wt_pop(styr,trmage);
                wt_sp_a = wt_spawn(styr,trmage);
            }
            else
            {
                wt_a    = waa_pop(a);
                wt_sp_a = waa_srv(1,a);
            }

            spawn_biomass(styr-1) += (0.5 * initN(a) * maa(a) * wt_sp_a * mfexp(-sp_frac * M(a)));
            total_biomass(styr-1) += (initN(a) * wt_a);

            if (a >= 3)
            {
                age_3_plus_biomass(styr-1) += (initN(a) * wt_a);
            }
        }
        else
        {
            spawn_biomass(styr-1) += (0.5 * initN(a) * maa(a) * wt_spawn(styr,a) * mfexp(-sp_frac * M(a)));
            total_biomass(styr-1) += (initN(a) * wt_pop(styr,a));

            if (a >= 3)
            {
                age_3_plus_biomass(styr-1) += (initN(a) * wt_pop(styr,a));
            }
        }
    }

    if (debug_flag) cout << "init of calc_n_a_a" << endl;

    // cout << "in calculate_n_a_a: R0 = " << R0 << ", h = " << h << ", SpBio(styr-1) = " << spawn_biomass(styr-1) << ", s_r_relat = " << s_r_relat << endl;

    // calculate the equilibrium number of recruits (age 1 in following year)
    spbio_tmp = spawn_biomass(styr-1) / conv_factor;
    steepness = 0.0;

    phi0.initialize();
    Ntmp.initialize();
    Ntmp = 0.5;
    for (a = st_age ; a < end_age; a++)
    {
        if (a < rcrage)
        {
            wt_sp_a = waa_srv(1,a);
        }
        else if (a > trmage)
        {
            if (complex_flag == 0 || complex_flag == 1 || complex_flag == 2)
            {
                wt_sp_a = wt_spawn(styr,trmage);
            }
            else
            {
                wt_sp_a = waa_srv(1,a);
            }
        }
        else
        {
            wt_sp_a = wt_spawn(styr,a);
        }

        phi0 += (Ntmp * maa(a) * wt_sp_a * mfexp(-sp_frac * M(a)));
        Ntmp *= mfexp(-M(a));
    }
    Ntmp /= (1.0 - mfexp(-M(end_age)));
    phi0 += (Ntmp * maa(end_age) * wt_sp_a * mfexp(-sp_frac * M(end_age)));

    // if (last_phase()) cout << "phi0 in year " << styr << " is " << phi0 << endl;

    if (s_r_relat == 1)
    {
        // Beverton-Holt S-R relationship (gamma = -1)
        N(styr,st_age) = conv_factor * (4.0 * R0 * h * spbio_tmp) / ((phi0 * R0 * (1.0 - h)) + (((5.0 * h) - 1.0) * spbio_tmp));
        steepness = h;
    }
    else if (s_r_relat == 2)
    {
        // Ricker S-R relationship (gamma = 0)
        N(styr,st_age) = conv_factor * (1.0 / phi0) * spbio_tmp * mfexp(h * (1.0 - (spbio_tmp / (R0 * phi0))));
        steepness = mfexp(h) / (mfexp(h) + 4.0);
    }
    else if (s_r_relat == 3)
    {
        // Schaefer S-R relationship (gamma = 1)
        N(styr,st_age) = conv_factor * R0 * spbio_tmp * (1.0 - (h * spbio_tmp));
    }
    else
    {
        // average recruits S-R relationship
        N(styr,st_age) = initR;
    }

    N(styr)((st_age+1),end_age) = initN((st_age+1),end_age);

    if (debug_flag) cout << "begin of calc_n_a_a" << endl;

    estsrvbio_bs.initialize();
    estsrvbio_ek.initialize();
    estsrvbio.initialize();
    pred_srv_paa_all.initialize();
    pred_srv_paa.initialize();

    for (y = styr; y <= endyr; y++)
    {
        // apply recruitment errors
        N(y,st_age) *= mfexp(log_rec_dev(y));

        for (a = (st_age+1); a <= end_age; a++)
        {
            // this assumes continuous fishing
            N(y+1,a) = N(y,a-1) * expZ(y,a-1);
        }
        N(y+1,end_age) = (N(y,end_age) * expZ(y,end_age)) + (N(y,end_age-1) * expZ(y,end_age-1));

        for (a = (st_age+1); a <= end_age; a++)
        {
            // spawning takes place after fishing seasons A and B
            // and survey 1 take place; need to fix this
            if (a < rcrage || a > trmage)
            {
                if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
                {
                    wt_a    = wt_pop(y,trmage);
                    wt_sp_a = wt_spawn(y,trmage);
                }
                else
                {
                    wt_a    = waa_pop(a);
                    wt_sp_a = waa_srv(1,a);
                }

                spawn_biomass(y) += (0.5 * N(y,a) * maa(a) * wt_sp_a * expZsp(y,a));
                total_biomass(y) += (N(y,a) * wt_a);

                if (a >= 3)
                {
                    age_3_plus_biomass(y) += (N(y,a) * wt_a);
                }
            }
            else
            {
                spawn_biomass(y) += (0.5 * N(y,a) * maa(a) * wt_spawn(y,a) * expZsp(y,a));
                total_biomass(y) += (N(y,a) * wt_pop(y,a));

                if (a >= 3)
                {
                    age_3_plus_biomass(y) += (N(y,a) * wt_pop(y,a));
                }
            }
        }

        // calculate the number of recruits (age 1 in following year)
        spbio_tmp = spawn_biomass(y) / conv_factor;
        steepness = 0.0;

        phi0.initialize();
        Ntmp.initialize();
        Ntmp = 0.5;
        for (a = st_age ; a < end_age; a++)
        {
            if (a < rcrage)
            {
                wt_sp_a = waa_srv(1,a);
            }
            else if (a > trmage)
            {
                if (complex_flag == 0 || complex_flag == 1 || complex_flag == 2)
                {
                    wt_sp_a = wt_spawn(y,trmage);
                }
                else
                {
                    wt_sp_a = waa_srv(1,a);
                }
            }
            else
            {
                wt_sp_a = wt_spawn(y,a);
            }

            phi0 += (Ntmp * maa(a) * wt_sp_a * mfexp(-sp_frac * M(a)));
            Ntmp *= mfexp(-M(a));
        }
        Ntmp /= (1.0 - mfexp(-M(end_age)));
        phi0 += (Ntmp * maa(end_age) * wt_sp_a * mfexp(-sp_frac * M(end_age)));
        // if (last_phase()) cout << "phi0 in year " << y << " is " << phi0 << endl;

        if (s_r_relat == 1)
        {
            // Beverton-Holt S-R relationship (gamma = -1)
            N(y+1,st_age) = conv_factor * (4.0 * R0 * h * spbio_tmp) / ((phi0 * R0 * (1.0 - h)) + (((5.0 * h) - 1.0) * spbio_tmp));
            steepness = h;
        }
        else if (s_r_relat == 2)
        {
            // Ricker S-R relationship (gamma = 0)
            N(y+1,st_age) = conv_factor * (1.0 / phi0) * spbio_tmp * mfexp(h * (1.0 - (spbio_tmp / (R0 * phi0))));
            steepness = mfexp(h) / (mfexp(h) + 4.0);
        }
        else if (s_r_relat == 3)
        {
            // Schaefer S-R relationship (gamma = 1)
            N(y+1,st_age) = conv_factor * R0 * spbio_tmp * (1.0 - (h * spbio_tmp));
        }
        else
        {
            // average recruits S-R relationship
            N(y+1,st_age) = R0;
        }

        // calculate survey biomass, converting from kg to millions of metric tonnes
        for (a = st_age; a <= end_age; a++)
        {
            if (a < rcrage || a > trmage)
            {
                if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
                {
                    // separate q and same selectivity for survey 1 Biosonics, survey 1 EK500, and survey 1 Dyson
                    i = 1;
                    srvbio_tmp       = wt_srv1(y,trmage) * Sage(i+nsel_fsh,a) * N(y,a) * expZfrac(i,y,a) / conv_factor;
                    estsrvbio_bs(y) += (q_srv_bs * srvbio_tmp);
                    estsrvbio_ek(y) += (q_srv_ek * srvbio_tmp);
                    estsrvbio(i,y)  += (q_srv(i) * srvbio_tmp);

                    i = 2;
                    estsrvbio(i,y) += (q_srv(i) * wt_srv2(y,trmage) * Sage(i+nsel_fsh,a) * N(y,a) * expZfrac(i,y,a) / conv_factor);

                    i = 3;
                    estsrvbio(i,y) += (q_srv(i) * wt_srv3(y,trmage) * Sage(i+nsel_fsh,a) * N(y,a) * expZfrac(i,y,a) / conv_factor);
                }
                else
                {
                    // separate q and same selectivity for survey 1 BBiosonics, survey 1 EK500, and survey 1 Dyson
                    i = 1;
                    srvbio_tmp       = waa_srv(i,a) * Sage(i+nsel_fsh,a) * N(y,a) * expZfrac(i,y,a) / conv_factor;
                    estsrvbio_bs(y) += (q_srv_bs * srvbio_tmp);
                    estsrvbio_ek(y) += (q_srv_ek * srvbio_tmp);
                    estsrvbio(i,y)  += (q_srv(i) * srvbio_tmp);

                    for (i = 2; i <= nsel_srv; i++)
                    {
                        estsrvbio(i,y) += (q_srv(i) * waa_srv(i,a) * Sage(i+nsel_fsh,a) * N(y,a) * expZfrac(i,y,a) / conv_factor);
                    }
                }
            }
            else
            {
                // separate q and same selectivity for survey 1 Biosonics, survey 1 EK500, and survey 1 Dyson
                i = 1;
                srvbio_tmp       = wt_srv1(y,a) * Sage(i+nsel_fsh,a) * N(y,a) * expZfrac(i,y,a) / conv_factor;
                estsrvbio_bs(y) += (q_srv_bs * srvbio_tmp);
                estsrvbio_ek(y) += (q_srv_ek * srvbio_tmp);
                estsrvbio(i,y)  += (q_srv(i) * srvbio_tmp);

                i = 2;
                estsrvbio(i,y) += (q_srv(i) * wt_srv2(y,a) * Sage(i+nsel_fsh,a) * N(y,a) * expZfrac(i,y,a) / conv_factor);

                i = 3;
                estsrvbio(i,y) += (q_srv(i) * wt_srv3(y,a) * Sage(i+nsel_fsh,a) * N(y,a) * expZfrac(i,y,a) / conv_factor);
            }
        }

        // calculate survey proportions-at-age
        for (i = 1; i <= nsel_srv; i++)
        {
            // AEP wonders if elem_prod is doing what it is supposed to be doing
            totn = sum(elem_prod(elem_prod(Sage(i+nsel_fsh), N(y)(st_age,end_age)), expZfrac(i,y)(st_age,end_age)));
            if (totn > 0)
            {
                pred_srv_paa_all(i,y) = (elem_prod(elem_prod(Sage(i+nsel_fsh), N(y)(st_age,end_age)), expZfrac(i,y)(st_age,end_age))) / totn;

                // apply ageing error
                pred_srv_paa_all(i,y) = pred_srv_paa_all(i,y) * age_trans_all;

                // reset for this particular year
                pred_srv_paa(i,y) = 0.0;

                // cover the common ages by default
                pred_srv_paa(i,y)(rcrage,trmage) = pred_srv_paa_all(i,y)(rcrage,trmage);

                // accumulate older ages into MWD plus group
                if (end_age > trmage)
                {
                    for (j = (trmage + 1); j <= end_age; j++)
                    {
                        pred_srv_paa(i,y,trmage) += pred_srv_paa_all(i,y,j);
                    }
                }
                // now renormalize the proportions
                if (sum(pred_srv_paa(i,y)) > 0.0)
                {
                    pred_srv_paa(i,y) /= sum(pred_srv_paa(i,y));
                }

                // need to do proportions-at-length but length bins currently don't match

            }
            else
            {
                pred_srv_paa_all(i,y) = 0.0;
                pred_srv_paa(i,y) = 0.0;
            }
        }
    }

    // now accumulate proportions-at-age to their respective age sets
    for (i = 1; i <= nsel_srv; i++)
    {
        if (i == 1)
        {
            j_max = nyrsac_srv1;
            for (j = 1; j <= j_max; j++)
            {
                y = srv_acyrs1(j);

                // reset for this particular year
                pred_srv_paa(i,y) = 0.0;

                pred_srv_paa(i,y)(ac_yng_srv1(j),ac_old_srv1(j)) = pred_srv_paa_all(i,y)(ac_yng_srv1(j),ac_old_srv1(j));

                // accumulate, since st_age (1) < rcrage <= ac_yng_srvN(j)
                if (rcrage < ac_yng_srv1(j))
                {
                    for (a = rcrage; a < ac_yng_srv1(j); a++)
                    {
                        pred_srv_paa(i,y,ac_yng_srv1(j)) += pred_srv_paa_all(i,y,a);
                    }
                }

                // accumulate, since ac_old_srvN(j) <= trmage (10) < end_age (15)
                for (a = (ac_old_srv1(j) + 1); a <= end_age; a++)
                {
                    pred_srv_paa(i,y,ac_old_srv1(j)) += pred_srv_paa_all(i,y,a);
                }

                // now renormalize the proportions
                if (sum(pred_srv_paa(i,y)) > 0.0)
                {
                    pred_srv_paa(i,y) /= sum(pred_srv_paa(i,y));
                }
            }
        }
        else if (i == 2)
        {
            j_max = nyrsac_srv2;
            for (j = 1; j <= j_max; j++)
            {
                y = srv_acyrs2(j);

                // reset for this particular year
                pred_srv_paa(i,y) = 0.0;

                pred_srv_paa(i,y)(ac_yng_srv2(j),ac_old_srv2(j)) = pred_srv_paa_all(i,y)(ac_yng_srv2(j),ac_old_srv2(j));

                // accumulate, since st_age (1) < rcrage <= ac_yng_srvN(j)
                if (rcrage < ac_yng_srv2(j))
                {
                    for (a = rcrage; a < ac_yng_srv2(j); a++)
                    {
                        pred_srv_paa(i,y,ac_yng_srv2(j)) += pred_srv_paa_all(i,y,a);
                    }
                }

                // accumulate, since ac_old_srvN(j) <= trmage (10) < end_age (15)
                for (a = (ac_old_srv2(j) + 1); a <= end_age; a++)
                {
                    pred_srv_paa(i,y,ac_old_srv2(j)) += pred_srv_paa_all(i,y,a);
                }

                // now renormalize the proportions
                if (sum(pred_srv_paa(i,y)) > 0.0)
                {
                    pred_srv_paa(i,y) /= sum(pred_srv_paa(i,y));
                }
            }
        }
        else if (i == 3)
        {
            j_max = nyrsac_srv3;
            for (j = 1; j <= j_max; j++)
            {
                y = srv_acyrs3(j);

                // reset for this particular year
                pred_srv_paa(i,y) = 0.0;

                pred_srv_paa(i,y)(rcrage,trmage) = pred_srv_paa_all(i,y)(rcrage,trmage);

                // accumulate, since trmage (10) < end_age (15)
                for (a = (trmage + 1); a <= end_age; a++)
                {
                    pred_srv_paa(i,y,trmage) += pred_srv_paa_all(i,y,a);
                }

                // now renormalize the proportions
                if (sum(pred_srv_paa(i,y)) > 0.0)
                {
                    pred_srv_paa(i,y) /= sum(pred_srv_paa(i,y));
                }
            }
        }
    }

    // sdreport variables
    // sd_q_srv_bs           = q_srv_bs;
    // sd_q_srv              = q_srv;
    sd_recruits_1         = column(N,1);
    sd_avg_rec_1          = mean(sd_recruits_1((om_hcr_styr+st_age),(endyr-1)));
    sd_avg_rec_1_alt      = mean(sd_recruits_1((om_rec_avg_styr+st_age),(endyr-1)));
    sd_recruits_2         = column(N,2);
    sd_avg_rec_2          = mean(sd_recruits_2((om_hcr_styr+st_age+1),(endyr-1)));
    sd_avg_rec_2_alt      = mean(sd_recruits_2((om_rec_avg_styr+st_age+1),(endyr-1)));
    sd_total_biomass      = total_biomass(styr,endyr);
    sd_spawn_biomass      = spawn_biomass(styr,endyr);
    sd_age_3_plus_biomass = age_3_plus_biomass(styr,endyr);
    // sd_estsrvbio_bs       = estsrvbio_bs;
    // sd_estsrvbio          = estsrvbio;

    // set the recruitment in the next year to be the average recruitment since om_hcr_styr
    N(endyr+1,st_age) = sd_avg_rec_1;

    if (debug_flag) cout << "end of calc_n_a_a" << endl;


FUNCTION calculate_c_a_a

    dvariable est_catch_num, wt_a;
    int y,a,j;

    C.initialize();
    pred_fsh_paa_all.initialize();
    pred_fsh_paa.initialize();

    if (debug_flag) cout << "start of calc_c_a_a" << endl;

    for (y = styr; y <= endyr; y++)
    {
        // calculate catch-at-age numbers
        for (a = st_age; a <= end_age; a++)
        {
            C(y,a) = F(y,a) * N(y,a) * (1.0 - expZ(y,a)) / Z(y,a);
        }

        est_catch_num = sum(C(y));

        // calculate proportion-at-age
        for (a = st_age; a <= end_age; a++)
        {
            if (est_catch_num > 0 && C(y,a) > 0)
            {
                pred_fsh_paa_all(y,a) = C(y,a) / est_catch_num;
            }
            else
            {
                pred_fsh_paa_all(y,a) = 0.0;
            }
        }

        // apply ageing error
        pred_fsh_paa_all(y) = pred_fsh_paa_all(y) * age_trans_all;

        // reset for this particular year
        pred_fsh_paa(y) = 0.0;

        // cover the common ages by default
        pred_fsh_paa(y)(rcrage,trmage) = pred_fsh_paa_all(y)(rcrage,trmage);

        // accumulate older ages into MWD plus group
        if (end_age > trmage)
        {
            for (a = (trmage + 1); a <= end_age; a++)
            {
                pred_fsh_paa(y,trmage) += pred_fsh_paa_all(y,a);
            }
        }

        // now renormalize the proportions
        if (sum(pred_fsh_paa(y)) > 0.0)
        {
            pred_fsh_paa(y) /= sum(pred_fsh_paa(y));
        }

       // need to do proportions-at-length but length bins currently don't match

    }

    // now accumulate proportions-at-age to their respective age sets
    for (j = 1; j <= nyrs_fsh; j++)
    {
        y = fshyrs(j);

        // reset for this particular year
        pred_fsh_paa(y) = 0.0;

        pred_fsh_paa(y)(ac_yng_fsh(j),ac_old_fsh(j)) = pred_fsh_paa_all(y)(ac_yng_fsh(j),ac_old_fsh(j));

        // accumulate, since st_age (1) < rcrage <= ac_yng_fsh(j)
        if (rcrage < ac_yng_fsh(j))
        {
            for (a = rcrage; a < ac_yng_fsh(j); a++)
            {
                pred_fsh_paa(y,ac_yng_fsh(j)) += pred_fsh_paa_all(y,a);
            }
        }

        // accumulate, since ac_old_fsh(j) <= trmage (10) < end_age (15)
        if (end_age > ac_old_fsh(j))
        {
            for (a = (ac_old_fsh(j) + 1); a <= end_age; a++)
            {
                   pred_fsh_paa(y,ac_old_fsh(j)) += pred_fsh_paa_all(y,a);
            }
        }

        // now renormalize the remaining proportions
        if (sum(pred_fsh_paa(y)) > 0.0)
        {
            pred_fsh_paa(y) /= sum(pred_fsh_paa(y));
        }
    }

    pred_catch.initialize();

    for (y = styr; y <= endyr; y++)
    {
        // ages covered by MWD model/data are 2-10;
        // ages covered in this model are 1-15
        for (a = st_age; a <= end_age; a++)
        {
            if (a < rcrage || a > trmage)
            {
                if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
                {
                    wt_a = wt_fsh(y,trmage);
                }
                else
                {
                    wt_a = waa_fsh(a);
                }

                pred_catch(y) += (wt_a * C(y,a));
            }
            else
            {
                pred_catch(y) += (wt_fsh(y,a) * C(y,a));
            }
        }
    }

    // convert from kg to tonnes
    pred_catch /= 1000.0;

    if (debug_flag) cout << "end of calc_c_a_a" << endl;


FUNCTION calculate_length_comps

    int uu,yy;

    pred_fsh_pal.initialize();
    pred_srv_1_pal.initialize();
    pred_srv_2_pal.initialize();
    pred_srv_3_pal.initialize();

    // calculate fishery catch-at-length
    for (uu = 1; uu <= nyrslen_fsh; uu++)
    {
        yy = fshlenyrs(uu);
        // pred_fsh_pal(i) = trans(len_trans1) * pred_fsh_paa(y);
        pred_fsh_pal(uu) = pred_fsh_paa(yy) * len_trans1;
    }

    // calculate EIT survey length composition
    for (uu = 1; uu <= nyrslen_srv1; uu++)
    {
        yy = srv_lenyrs1(uu);
        // pred_srv_1_pal(i) = trans(len_trans3) * pred_srv_paa(1,y);
        pred_srv_1_pal(uu) = pred_srv_paa(1,yy) * len_trans3;
    }

    // calculate NMFS GOA Triennial/Biennial bottom trawl survey length composition
    for (uu = 1; uu <= nyrslen_srv2; uu++)
    {
        yy = srv_lenyrs2(uu);
        // pred_srv_2_pal(i) = trans(len_trans2) * pred_srv_paa(2,y);
        pred_srv_2_pal(uu) = pred_srv_paa(2,yy) * len_trans2;
    }

    // calculate ADF&G coastal survey length composition
    for (uu = 1; uu <= nyrslen_srv3; uu++)
    {
        yy = srv_lenyrs3(uu);
        // pred_srv_3_pal(i) = trans(len_trans2) * pred_srv_paa(3,y);
        pred_srv_3_pal(uu) = pred_srv_paa(3,yy) * len_trans2;
    }

    if (debug_flag) cout << "end of calc_length_comps" << endl;


FUNCTION void calculate_bio_ref_points(dvar_vector sel, int curr_yr, dvar_vector Mvec)

    // calculate the biological reference points for curr_yr+1
    // NOTE:  this assumes that N(curr_yr+1,a) has been filled in and
    // is the numbers-at-age at the beginning of year (curr_yr+1)

    int a, i, next_yr;
    dvariable wt_a, wt_sp_a, spbio_tmp;
    dvariable f_target, spbio_target, spbio_floor, curr_spbio;
    dvar_vector Z_ABC(st_age,end_age), Z_OFL(st_age,end_age);

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

    F100(curr_yr)  = get_spr_rates(1.00, sel, curr_yr, Mvec);
    SB100(curr_yr) = SBcurr;

    F40(curr_yr)   = get_spr_rates(0.40, sel, curr_yr, Mvec);
    SB40(curr_yr)  = SBcurr;

    F35(curr_yr)   = get_spr_rates(0.35, sel, curr_yr, Mvec);
    SB35(curr_yr)  = SBcurr;

    F20(curr_yr)   = get_spr_rates(0.20, sel, curr_yr, Mvec);
    SB20(curr_yr)  = SBcurr;

    SBtarget(curr_yr) = SB40(curr_yr) * (F35(curr_yr) / F40(curr_yr));

    // initialize
    F_ABC = F40(curr_yr);

    // iterate 10 times to get stable values for F_ABC and F_OFL
    for (i = 1; i <= 10; i++)
    {
        f_target     = F40(curr_yr);
        spbio_target = SBtarget(curr_yr);
        spbio_floor  = SB20(curr_yr);

        spawn_biomass(next_yr) = total_biomass(next_yr) = age_3_plus_biomass(next_yr) = 0.0;
        for (a = (st_age+1); a <= end_age; a++)
        {
            // spawning takes place after fishing seasons A and B
            // and survey 1 take place; need to fix this
            if (a < rcrage || a > trmage)
            {
                if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
                {
                    wt_a    = wt_pop_proj(trmage);
                    wt_sp_a = wt_spawn_proj(trmage);
                }
                else
                {
                    wt_a    = waa_pop(a);
                    wt_sp_a = waa_srv(1,a);
                }

                spawn_biomass(next_yr) += (0.5 * N(next_yr,a) * maa(a) * wt_sp_a * mfexp(-sp_frac * (Mvec(a) + (sel(a) * F_ABC))));
                total_biomass(next_yr) += (N(next_yr,a) * wt_a);

                if (a >= 3)
                {
                    age_3_plus_biomass(next_yr) += (N(next_yr,a) * wt_a);
                }
            }
            else
            {
                spawn_biomass(next_yr) += (0.5 * N(next_yr,a) * maa(a) * wt_spawn_proj(a) * mfexp(-sp_frac * (Mvec(a) + (sel(a) * F_ABC))));
                total_biomass(next_yr) += (N(next_yr,a) * wt_pop_proj(a));

                if (a >= 3)
                {
                    age_3_plus_biomass(next_yr) += (N(next_yr,a) * wt_pop_proj(a));
                }
            }
        }

        spbio_tmp = spawn_biomass(next_yr) / conv_factor;

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
            cout << "calculate_bio_ref_points:\tSB <= SB20% in year " << next_yr << ":\tSB20 "<< spbio_floor << "\tSpawning biomass " << spawn_biomass(next_yr) << "\tTotal biomass " << total_biomass(next_yr) << endl;
        }

        f_target     = F35(curr_yr);
        spbio_target = SB40(curr_yr);

        if (spbio_tmp < spbio_target)
        {
            F_OFL = f_target * ((spbio_tmp / spbio_target) - Tier3_alpha) / (1.0 - Tier3_alpha);
        }
        else
        {
            F_OFL = f_target;
        }
    }

    Z_ABC = (sel * F_ABC) + Mvec;
    Z_OFL = (sel * F_OFL) + Mvec;

    OFL = ABC = 0.0;
    for (a = st_age; a <= end_age; a++)
    {
        if (a < rcrage || a > trmage)
        {
            if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
            {
                wt_a    = wt_fsh_proj(trmage);
            }
            else
            {
                wt_a    = waa_fsh(a);
            }
        }
        else
        {
            wt_a = wt_fsh_proj(a);
        }

        ABC += (wt_a * N(next_yr,a) * (1.0 - mfexp(-Z_ABC(a))) * (sel(a) * F_ABC) / (Z_ABC(a)));
        OFL += (wt_a * N(next_yr,a) * (1.0 - mfexp(-Z_OFL(a))) * (sel(a) * F_OFL) / (Z_OFL(a)));
    }

    // convert from kg to mt
    ABC /= 1000.0;
    OFL /= 1000.0;

    if (debug_flag) cout << "end of calculate_bio_ref_points" << endl;


FUNCTION evaluate_objective_function

    int y,i,j,k;
    bool extYear = false;

    f.initialize();

    f = 0.0;

    // likelihood function for catch totals
    for (y = styr; y <= endyr; y++)
    {
        if (cattot(y) > 0 && pred_catch(y) > 0 && cattot_log_sd(y) > 0)
        {
            f(1) += (0.5 * square((log(cattot(y) + o) - log(pred_catch(y) + o)) / cattot_log_sd(y)));
        }
    }

    // likelihood function for fishery catch-at-age
    for (i = 1; i <= nyrs_fsh; i++)
    {
        y = fshyrs(i);

        // check if there is extended p-a-a data for this year
        extYear = false;
        for (k = 1; k <= nyrs_fsh_paa_all && (!extYear); k++)
        {
            extYear = (y == yrs_fsh_paa_all(k));
        }

        if (!extYear)
        {
            for (j = ac_yng_fsh(i); j <= ac_old_fsh(i); j++)
            {
                f(2) -= (multN_fsh(i) * (catp(i,j) + o) * (log(pred_fsh_paa(y,j) + o) - log(catp(i,j) + o)));
            }
        }
    }

    // likelihood function for fishery catch-at-length
    for (i = 1; i <= nyrslen_fsh; i++)
    {
        for (j = 1; j <= nbins1; j++)
        {
            f(3) -= (multNlen_fsh(i) * (lenp(i,j) + o) * (log(pred_fsh_pal(i,j) + o) - log(lenp(i,j) + o)));
        }
    }

    // likelihood function for survey biomass estimates
    for (i = 1; i <= nyrs_srv1_bs; i++)     // BS survey (EIT survey)
    {
        y = srvyrs1_bs(i);
        f(5) += (0.5 * square((log(indxsurv1_bs(i)) - log(estsrvbio_bs(y)) + (square(indxsurv_log_sd1_bs(i)) / 2.0)) / indxsurv_log_sd1_bs(i)));
    }

    for (i = 1; i <= nyrs_srv1_ek; i++)     // EK500 survey (EIT survey)
    {
        y = srvyrs1_ek(i);
        f(5) += (0.5 * square((log(indxsurv1_ek(i)) - log(estsrvbio_ek(y)) + (square(indxsurv_log_sd1_ek(i)) / 2.0)) / indxsurv_log_sd1_ek(i)));
    }

    for (i = 1; i <= nyrs_srv1_dy; i++)     // Dyson survey (EIT survey)
    {
        y = srvyrs1_dy(i);
        f(5) += (0.5 * square((log(indxsurv1_dy(i)) - log(estsrvbio(1,y)) + (square(indxsurv_log_sd1_dy(i)) / 2.0)) / indxsurv_log_sd1_dy(i)));
    }

    for (i = 1; i <= nyrs_srv2; i++)        // Triennial/Biennial bottom trawl survey
    {
        y = srvyrs2(i);
        f(5) += (0.5 * square((log(indxsurv2(i)) - log(estsrvbio(2,y)) + (square(indxsurv_log_sd2(i)) / 2.0)) / indxsurv_log_sd2(i)));
    }

    for (i = 1; i <= nyrs_srv3; i++)        // ADF&G coastal survey
    {
        y = srvyrs3(i);
        f(5) += (0.5 * square((log(indxsurv3(i)) - log(estsrvbio(3,y)) + (square(indxsurv_log_sd3(i)) / 2.0)) / indxsurv_log_sd3(i)));
    }

    // likelihood function for survey prop-at-age
    for (i = 1; i <= nyrsac_srv1; i++)      // EIT survey
    {
        y = srv_acyrs1(i);

        // check if there is extended p-a-a data for this year
        extYear = false;
        for (k = 1; k <= nyrs_srv_1_paa_all && (!extYear); k++)
        {
            extYear = (y == yrs_srv_1_paa_all(k));
        }

        if (!extYear)
        {
            for (j = ac_yng_srv1(i); j <= ac_old_srv1(i); j++)
            {
                f(6) -= (multN_srv1(i) * (srvp1(i,j) + o) * (log(pred_srv_paa(1,y,j) + o) - log(srvp1(i,j) + o)));
            }
        }
    }

    for (i = 1; i <= nyrsac_srv2; i++)      // NMFS GOA Triennial/Biennial bottom trawl survey
    {
        y = srv_acyrs2(i);

        // check if there is extended p-a-a data for this year
        extYear = false;
        for (k = 1; k <= nyrs_srv_2_paa_all && (!extYear); k++)
        {
            extYear = (y == yrs_srv_2_paa_all(k));
        }

        if (!extYear)
        {
            for (j = ac_yng_srv2(i); j <= ac_old_srv2(i); j++)
            {
                f(6) -= (multN_srv2(i) * (srvp2(i,j) + o) * (log(pred_srv_paa(2,y,j) + o) - log(srvp2(i,j) + o)));
            }
        }
    }

    for (i = 1; i <= nyrsac_srv3; i++)      // ADF&G coastal survey
    {
        y = srv_acyrs3(i);
        for (j = rcrage; j <= trmage; j++)
        {
            f(6) -= (multN_srv3(i) * (srvp3(i,j) + o) * (log(pred_srv_paa(3,y,j) + o) - log(srvp3(i,j) + o)));
        }
    }

    // likelihood function for survey prop-at-length
    for (i = 1; i <= nyrslen_srv1; i++)     // EIT survey
    {
        for (j = 1; j <= nbins3; j++)
        {
            f(7) -= (multNlen_srv1(i) * (srvlenp1(i,j) + o) * (log(pred_srv_1_pal(i,j) + o) - log(srvlenp1(i,j) + o)));
        }
    }

    for (i = 1; i <= nyrslen_srv2; i++)     // NMFS GOA Triennial/Biennial bottom trawl survey
    {
        for (j = 1; j <= nbins2; j++)
        {
            f(7) -= (multNlen_srv2(i) * (srvlenp2(i,j) + o) * (log(pred_srv_2_pal(i,j) + o) - log(srvlenp2(i,j) + o)));
        }
    }

    for (i = 1; i <= nyrslen_srv3; i++)     // ADF&G coastal survey
    {
        for (j = 1; j <= nbins2; j++)
        {
            f(7) -= (multNlen_srv3(i) * (srvlenp3(i,j) + o) * (log(pred_srv_3_pal(i,j) + o) - log(srvlenp3(i,j) - o)));
        }
    }

    // likelihood function for expanded p-a-a data from 2004 GOA walleye pollock SAFE (pg 69 and 73)
    for (i = 1; i <= nyrs_fsh_paa_all; i++)     // extended fishery p-a-a data
    {
        y = yrs_fsh_paa_all(i);
        for (j = st_age; j <= end_age; j++)
        {
            f(8) -= (multN_fsh_paa_all(i) * (fsh_paa_all(i,j) + o) * (log(pred_fsh_paa_all(y,j) + o) - log(fsh_paa_all(i,j) + o)));
        }
    }

    for (i = 1; i <= nyrs_srv_1_paa_all; i++)   // extended EIT survey p-a-a data
    {
        y = yrs_srv_1_paa_all(i);
        for (j = st_age; j <= end_age; j++)
        {
            f(8) -= (multN_srv_1_paa_all(i) * (srv_1_paa_all(i,j) + o) * (log(pred_srv_paa_all(1,y,j) + o) - log(srv_1_paa_all(i,j) + o)));
        }
    }

    for (i = 1; i <= nyrs_srv_2_paa_all; i++)   // extended Bottom trawl survey p-a-a data
    {
        y = yrs_srv_2_paa_all(i);
        for (j = st_age; j <= end_age; j++)
        {
            f(8) -= (multN_srv_2_paa_all(i) * (srv_2_paa_all(i,j) + o) * (log(pred_srv_paa_all(2,y,j) + o) - log(srv_2_paa_all(i,j) + o)));
        }
    }

    // likelihood function for recruitment residuals
    // f(9)  = (0.5 / square(sigma_r)) * norm2(log_rec_dev + square(sigma_r)/2.0);
    // match pk13_3.tpl
    f(9) = (0.5 * square(log_rec_dev(first_rec_year)/1.0)) + (0.5 * norm2(log_rec_dev(first_rec_year+1,first_rec_year+7)/1.0)) + (0.5 * norm2(log_rec_dev(last_rec_year-1,last_rec_year)/1.0));

    // likelihood function for residuals in the first year
    f(10) = (0.5 / square(sigma_f)) * norm2(log_init_dev);

    // random walk for changing fishery selectivity
    f(12)  = (0.5 * norm2(elem_div(first_difference(fsh_sel_cont_a1_dev), 4.0 * rwlk_sd(styr_fsh_sel_dev,endyr-1))));
    f(12) += (0.5 * norm2(elem_div(first_difference(fsh_sel_cont_b1_dev), 1.0 * rwlk_sd(styr_fsh_sel_dev,endyr-1))));
    f(12) += (0.5 * norm2(elem_div(first_difference(fsh_sel_cont_a2_dev), 4.0 * rwlk_sd(styr_fsh_sel_dev,endyr-1))));
    f(12) += (0.5 * norm2(elem_div(first_difference(fsh_sel_cont_b2_dev), 1.0 * rwlk_sd(styr_fsh_sel_dev,endyr-1))));

    // likelihood component for scaling srv1 Dyson q to Miller Freeman EK500 q (see pk10_1.tpl, line1243)
    f(17) = (1.0 / (2.0 * (square(0.0244) + square(0.000001)))) * square(log_q_srv(1) - log_q_srv_ek - 0.124);

    // put a prior on R0 for the S-R relationships
    if (s_r_relat == 1 || s_r_relat == 2 || s_r_relat == 3)
    {
        f(18) = (0.5 / square(R0_stddev)) * square((1.0e9 * R0) - R0_mean);
    }

    // likelihood function to ensure that N(1964,1) and N(1965,1) are close
    f(20) = 1000.0 * square(log(N(styr,st_age)) - log(N(styr+1,st_age)));

    obj_fun = sum(f);

    var_prof=mean(sd_recruits_1((om_hcr_styr+st_age),endyr-1));

    if (debug_flag) cout << "end of calc_obj_fun" << endl;


REPORT_SECTION

    int y,a,i,j;

    if (debug_flag) cout << "begin of report_section" << endl;

    if (last_phase())
    {
        // calculate average fishery selectivity from 1992 through endyr-1
        // NOTE:  MWD uses (endyr-5) through (endyr-1) in pk10_1.tpl
        Sage_fsh_avg.initialize();
        Sage_fsh_avg = 0.0;
        for (a = st_age; a <= end_age; a++)
        {
            // for (y = 1992; y < endyr; y++)
            for (y = (endyr-5); y < endyr; y++)
            {
                Sage_fsh_avg(a) += Sage_fsh(y,a);
            }
        }
        Sage_fsh_avg /= max(Sage_fsh_avg);

        calculate_bio_ref_points(Sage_fsh_avg,endyr,M);
    }

    report << "R0" << endl;
    report << R0 << endl;
    report << "h" << endl;
    report << h << endl;
    report << "steepness" << endl;
    report << steepness << endl;
    report << "initR" << endl;
    report << initR << endl;
    report << "M" << endl;
    report << M << endl;
    report << "initial age structure" << endl;
    report << (initN(st_age,end_age) / initN(st_age)) << endl;
    report << "q survey" << endl;
    report << q_srv_bs << "\t" << q_srv_ek << "\t" << q_srv << endl;
    report << endl;

    report << "year \t catch \t est catch" << endl;
    for (y = styr; y <= endyr; y++)
    {
        report << y << "\t" << cattot(y) << "\t" << pred_catch(y) << endl;
    }
    report << endl;

    report << "survey 1 BS biomass estimates" << endl;
    for (i = 1; i <= nyrs_srv1_bs; i++)     // BS survey (EIT survey)
    {
        y = srvyrs1_bs(i);
        report << y << "\t" << indxsurv1_bs(i) << "\t" << indxsurv_log_sd1_bs(i) << "\t" << estsrvbio_bs(y) << endl;
    }
    report << endl;

    report << "survey 1 EK biomass estimates" << endl;
    for (i = 1; i <= nyrs_srv1_ek; i++)     // EK500 survey (EIT survey)
    {
        y = srvyrs1_ek(i);
        report << y << "\t" << indxsurv1_ek(i) << "\t" << indxsurv_log_sd1_ek(i) << "\t" << estsrvbio_ek(y) << endl;
    }
    report << endl;

    report << "survey 1 Dyson biomass estimates" << endl;
    for (i = 1; i <= nyrs_srv1_dy; i++)     // Dyson survey (EIT survey)
    {
        y = srvyrs1_dy(i);
        report << y << "\t" << indxsurv1_dy(i) << "\t" << indxsurv_log_sd1_dy(i) << "\t" << estsrvbio(1,y) << endl;
    }
    report << endl;

    report << "survey 2 biomass estimates" << endl;
    for (i = 1; i <= nyrs_srv2; i++)        // Triennial/Biennial bottom trawl survey
    {
        y = srvyrs2(i);
        report << y << "\t" << indxsurv2(i) << "\t" << indxsurv_log_sd2(i) << "\t" << estsrvbio(2,y) << endl;
    }
    report << endl;

    report << "survey 3 biomass estimates" << endl;
    for (i = 1; i <= nyrs_srv3; i++)        // ADF&G coastal survey
    {
        y = srvyrs3(i);
        report << y << "\t" << indxsurv3(i) << "\t" << indxsurv_log_sd3(i) << "\t" << estsrvbio(3,y) << endl;
    }
    report << endl;

    // report << "Selectivity at length" << endl;
    // report << Slen << endl;
    // report << endl;
    report << "Selectivity at age" << endl;
    report << "fishery A\t" << Sage(1) << endl;
    report << "fishery B\t" << Sage(2) << endl;
    report << "fishery C\t" << Sage(3) << endl;
    report << "fishery D\t" << Sage(4) << endl;
    report << "survey 1\t" << Sage(5) << endl;
    report << "survey 2\t" << Sage(6) << endl;
    report << "survey 3\t" << Sage(7) << endl;
    report << "fishery\t" << Sage(8) << endl;
    report << endl;

    report << "Fishery selectivity deviances" << endl;
    report << "a1\t" << fsh_sel_cont_a1_dev << endl;
    report << "b1\t" << fsh_sel_cont_b1_dev << endl;
    report << "a2\t" << fsh_sel_cont_a2_dev << endl;
    report << "b2\t" << fsh_sel_cont_b2_dev << endl;
    report << endl;

    report << "Fishing mortality" << endl;
    report << Fmort << endl;
    report << endl;

    report << "Year\tTotBio\tSpBio\tAge 3+ Bio\tNumbers-at-age" << endl;
    for (y = styr; y <= (endyr + 1); y++)
    {
        report << y << "\t" << (total_biomass(y) / conv_factor) << "\t" << (spawn_biomass(y) / conv_factor) << "\t" << (age_3_plus_biomass(y) / conv_factor) << "\t" << N(y) << endl;
    }
    report << endl;

    report << "Fishery c-a-a, actual" << endl;
    for (i = 1; i <= nyrs_fsh; i++)
    {
        y = fshyrs(i);
        // if (multN_fsh(i) > 0)
        // {
            report << y << "\t" << catp(i)(ac_yng_fsh(i),ac_old_fsh(i)) << endl;
        // }

    }
    report << endl;
    report << "Fishery c-a-a, estimated" << endl;
    for (i = 1; i <= nyrs_fsh; i++)
    {
        y = fshyrs(i);
        // if (multN_fsh(i) > 0)
        // {
            report << y << "\t" << pred_fsh_paa(y)(ac_yng_fsh(i),ac_old_fsh(i)) << endl;
        // }
    }
    report << endl;

    report << "Survey 1 age composition, actual" << endl;
    for (i = 1; i <= nyrsac_srv1; i++)      // EIT survey
    {
        y = srv_acyrs1(i);
        // if (multN_srv1(i) > 0)
        // {
            report << y << "\t" << srvp1(i)(ac_yng_srv1(i),ac_old_srv1(i)) << endl;
        // }
    }
    report << endl;
    report << "Survey 1 age composition, estimated" << endl;
    for (i = 1; i <= nyrsac_srv1; i++)      // EIT survey
    {
        y = srv_acyrs1(i);
        // if (multN_srv1(i) > 0)
        // {
            report << y << "\t" << pred_srv_paa(1,y)(ac_yng_srv1(i),ac_old_srv1(i)) << endl;
        // }
    }
    report << endl;

    report << "Survey 2 age composition, actual" << endl;
    for (i = 1; i <= nyrsac_srv2; i++)      // NMFS GOA Triennial/Biennial bottom trawl survey
    {
        y = srv_acyrs2(i);
        // if (multN_srv2(i) > 0)
        // {
            report << y << "\t" << srvp2(i)(ac_yng_srv2(i),ac_old_srv2(i)) << endl;
        // }
    }
    report << endl;
    report << "Survey 2 age composition, estimated" << endl;
    for (i = 1; i <= nyrsac_srv2; i++)      // NMFS GOA Triennial/Biennial bottom trawl survey
    {
        y = srv_acyrs2(i);
        // if (multN_srv2(i) > 0)
        // {
            report << y << "\t" << pred_srv_paa(2,y)(ac_yng_srv2(i),ac_old_srv2(i)) << endl;
        // }
    }
    report << endl;

    report << "Survey 3 age composition, actual" << endl;
    for (i = 1; i <= nyrsac_srv3; i++)      // ADF&G coastal survey
    {
        y = srv_acyrs3(i);
        // if (multN_srv3(i) > 0)
        // {
            report << y << "\t" << srvp3(i) << endl;
        // }
    }
    report << endl;
    report << "Survey 3 age composition, estimated" << endl;
    for (i = 1; i <= nyrsac_srv3; i++)      // ADF&G coastal survey
    {
        y = srv_acyrs3(i);
        // if (multN_srv3(i) > 0)
        // {
            report << y << "\t" << pred_srv_paa(3,y) << endl;
        // }
    }
    report << endl;

    report << "Fishery c-a-l, actual" << endl;
    for (i = 1; i <= nyrslen_fsh; i++)
    {
        y = fshlenyrs(i);
        // if (multNlen_fsh(i) > 0)
        // {
            report << y << "\t" << lenp(i) << endl;
        // }

    }
    report << endl;
    report << "Fishery c-a-l, estimated" << endl;
    for (i = 1; i <= nyrslen_fsh; i++)
    {
        y = fshlenyrs(i);
        // if (multNlen_fsh(i) > 0)
        // {
            report << y << "\t" << pred_fsh_pal(i) << endl;
        // }
    }
    report << endl;

    report << "Survey 1 length composition, actual" << endl;
    for (i = 1; i <= nyrslen_srv1; i++)      // EIT survey
    {
        y = srv_lenyrs1(i);
        // if (multNlen_srv1(i) > 0)
        // {
            report << y << "\t" << srvlenp1(i) << endl;
        // }
    }
    report << endl;
    report << "Survey 1 length composition, estimated" << endl;
    for (i = 1; i <= nyrslen_srv1; i++)      // EIT survey
    {
        y = srv_lenyrs1(i);
        // if (multNlen_srv1(i) > 0)
        // {
            report << y << "\t" << pred_srv_1_pal(i) << endl;
        // }
    }
    report << endl;

    report << "Survey 2 length composition, actual" << endl;
    for (i = 1; i <= nyrslen_srv2; i++)      // NMFS GOA Triennial/Biennial bottom trawl survey
    {
        y = srv_lenyrs2(i);
        // if (multNlen_srv2(i) > 0)
        // {
            report << y << "\t" << srvlenp2(i) << endl;
        // }
    }
    report << endl;
    report << "Survey 2 length composition, estimated" << endl;
    for (i = 1; i <= nyrslen_srv2; i++)      // Triennial/Biennial bottom trawl survey
    {
        y = srv_lenyrs2(i);
        // if (multNlen_srv2(i) > 0)
        // {
            report << y << "\t" << pred_srv_2_pal(i) << endl;
        // }
    }
    report << endl;

    report << "Survey 3 length composition, actual" << endl;
    for (i = 1; i <= nyrslen_srv3; i++)      // ADF&G coastal survey
    {
        y = srv_lenyrs3(i);
        // if (multNlen_srv3(i) > 0)
        // {
            report << y << "\t" << srvlenp3(i) << endl;
        // }
    }
    report << endl;
    report << "Survey 3 length composition, estimated" << endl;
    for (i = 1; i <= nyrslen_srv3; i++)      // ADF&G coastal survey
    {
        y = srv_lenyrs3(i);
        // if (multNlen_srv3(i) > 0)
        // {
            report << y << "\t" << pred_srv_3_pal(i) << endl;
        // }
    }
    report << endl;

    report << "Fishery extended catch-at-age, actual" << endl;
    for (i = 1; i <= nyrs_fsh_paa_all; i++)     // extended fishery p-a-a data
    {
        y = yrs_fsh_paa_all(i);
        report << y << "\t" << fsh_paa_all(i) << endl;
    }
    report << endl;
    report << "Fishery extended catch-at-age, estimated" << endl;
    for (i = 1; i <= nyrs_fsh_paa_all; i++)     // extended fishery p-a-a data
    {
        y = yrs_fsh_paa_all(i);
        report << y << "\t" << pred_fsh_paa_all(y) << endl;
    }
    report << endl;

    report << "Survey 1 extended age composition, actual" << endl;
    for (i = 1; i <= nyrs_srv_1_paa_all; i++)   // extended EIT survey p-a-a data
    {
        y = yrs_srv_1_paa_all(i);
        report << y << "\t" << srv_1_paa_all(i) << endl;
    }
    report << endl;
    report << "Survey 1 extended age composition, estimated" << endl;
    for (i = 1; i <= nyrs_srv_1_paa_all; i++)   // extended EIT survey p-a-a data
    {
        y = yrs_srv_1_paa_all(i);
        report << y << "\t" << pred_srv_paa_all(1,y) << endl;
    }
    report << endl;

    report << "Survey 2 extended age composition, actual" << endl;
    for (i = 1; i <= nyrs_srv_2_paa_all; i++)   // extended Bottom trawl survey p-a-a data
    {
        y = yrs_srv_2_paa_all(i);
        report << y << "\t" << srv_2_paa_all(i) << endl;
    }
    report << endl;
    report << "Survey 2 extended age composition, estimated" << endl;
    for (i = 1; i <= nyrs_srv_2_paa_all; i++)   // extended Bottom trawl survey p-a-a data
    {
        y = yrs_srv_2_paa_all(i);
        report << y << "\t" << pred_srv_paa_all(2,y) << endl;
    }
    report << endl;

    report << "Fishery selectivity by year" << endl;
    for (i = styr; i <= endyr; i++)
    {
        report << i << "\t" << Sage_fsh(i) << endl;
    }
    report << endl;

    report << "Estimated survey biomass by year" << endl;
    for (i = styr; i <= endyr; i++)
    {
        report << i << "\t" << estsrvbio_bs(i) << "\t" << estsrvbio_ek(i);
        for (j = 1; j <= nsel_srv; j++)
        {
            report << "\t" << estsrvbio(j,i);
        }
        report << endl;
    }
    report << endl;

    // print out some stuff about the decision rule using the average fishery selectivity
    if (last_phase())
    {
        report << "Decision rule parameter values" << endl;
        report << "xx\tFxx\tSBxx" << endl;

        report << "1.00\t" << F100(endyr) << "\t" << SB100(endyr) << endl;

        report << "0.40\t" << F40(endyr) << "\t" << SB40(endyr) << endl;

        report << "0.35\t" << F35(endyr) << "\t" << SB35(endyr) << endl;

        report << "0.20\t" << F20(endyr) << "\t" << SB20(endyr) << endl;

        report << endl;
        avgR_CV = std_dev(log(sd_recruits_1((om_hcr_styr+st_age),(endyr-1))));
        report << "SB0 = " << SB0 << ", avgR = " << avgR << ", avgR_CV = " << avgR_CV << ", phi0 = " << phi0 << endl;
        report << "SBtarget = " << SBtarget(endyr) << ", which is SB" << (100.0 * (SBtarget(endyr) / (SB40(endyr) / 0.4))) << endl;
        report << endl;

        report << "Decision rule parameter values for year " << (endyr+1) << endl;
        report << "Projected fishing selectivity (avg of prev 5 years)\t" << Sage_fsh_avg << endl;
        report << "Estimated " << (endyr+1) << " spawning biomass\t" << (spawn_biomass(endyr+1) / conv_factor) << endl;
        report << "Estimated " << (endyr+1) << " total biomass\t" << (total_biomass(endyr+1) / conv_factor) << endl;
        report << "Estimated " << (endyr+1) << " age 3+ biomass\t" << (age_3_plus_biomass(endyr+1) / conv_factor) << endl;
        report << "F_ABC\t" << F_ABC << endl;
        report << "ABC\t" << ABC << endl;
        report << "F_OFL\t" << F_OFL << endl;
        report << "OFL\t" << OFL << endl;
        report << endl;
    }

    report << "objective function" << endl;
    report << obj_fun << endl;
    report << "objective function components" << endl;
    report << f << endl;
    report << endl;

    if (debug_flag) cout << "end of report_section" << endl;


FUNCTION calculate_projections

    int uu,yy,y,a;
    dvariable TAC, prev_spawn_biomass;
    dvector actual_paa(st_age,end_age);
    dvector rand_paa(st_age,end_age);

    // initialize all of the projection parameters
    N_proj.initialize();
    C_proj.initialize();
    F_proj.initialize();
    M_proj.initialize();
    Z_proj.initialize();
    expZ_proj.initialize();
    expZsp_proj.initialize();
    expZfrac_proj.initialize();
    fsh_paa_all_proj.initialize();
    srv_paa_all_proj.initialize();
    catch_proj.initialize();
    spawn_biomass_proj.initialize();
    total_biomass_proj.initialize();
    age_3_plus_biomass_proj.initialize();
    estsrvbio_proj.initialize();
    estsrvN_proj.initialize();
    estsrvbio_proj_2_sd.initialize();
    estsrvbio_proj_2_multN.initialize();

    F100.initialize();
    SB100.initialize();
    SBtarget.initialize();
    F40.initialize();
    SB40.initialize();
    F35.initialize();
    SB35.initialize();
    F20.initialize();
    SB20.initialize();

    Sage_fsh_avg.initialize();

    estFABC_proj.initialize();
    estABC_proj.initialize();
    estABCaa_proj.initialize();

    B0_mw.initialize();

    mcmc_avg_rec.initialize();
    mcmc_avg_log_rec.initialize();
    mcmc_avg_rec_years.initialize();
    mcmc_std_dev_avg_log_rec.initialize();
    mcmc_CV_avg_log_rec.initialize();
    mcmc_max_rec.initialize();

    avgR_CV.initialize();

    if (mceval_phase())
    {
        ofstream pfile("opmodel.prj", ios::app | ios::ate);
        ofstream mcmcfile("opm_mcmc.dat", ios::out | ios::app);

        mcmc_iter++;

        // M-at-age for the projections
        M_proj = init_M_proj;

        // calculate average fishery selectivity from 1992 through endyr-1
        // NOTE:  MWD uses (endyr-5) through (endyr-1) in pk10_1.tpl
        Sage_fsh_avg.initialize();
        Sage_fsh_avg = 0.0;
        for (a = st_age; a <= end_age; a++)
        {
            // for (y = 1992; y < endyr; y++)
            for (y = (endyr-5); y < endyr; y++)
            {
                Sage_fsh_avg(a) += Sage_fsh(y,a);
            }
        }
        Sage_fsh_avg /= max(Sage_fsh_avg);

        calculate_bio_ref_points(Sage_fsh_avg,endyr,M);

        // calculate values used to generate future recruitment
        mcmc_avg_rec_years = (endyr-1) - (om_rec_avg_styr+st_age) + 1;
        mcmc_avg_rec = mean(sd_recruits_1((om_rec_avg_styr+st_age),(endyr-1)));
        mcmc_avg_log_rec = mean(log(sd_recruits_1((om_rec_avg_styr+st_age),(endyr-1))));
        avgR_CV = mcmc_CV_avg_log_rec = std_dev(log(sd_recruits_1((om_rec_avg_styr+st_age),(endyr-1))));
        mcmc_max_rec = max(sd_recruits_1(styr,(endyr-1)));

        // output the parameter values for this MCMC iterations
        if (mcmcfile.is_open())
        {
            mcmcfile << mcmc_iter << "\t" << obj_fun <<  "\t"  << log_R0 << "\t" << log_h << "\t" << log_q << "\t" << log_q_fsh << "\t" << log_q_srv_bs << "\t" << log_q_srv_ek << "\t" << log_q_srv << "\t" << log_initR << "\t" << Fmort << "\t" << log_rec_dev << "\t" << log_init_dev << "\t" << log_fsh_sel_cont_a1 << "\t" << log_fsh_sel_cont_b1 << "\t" << log_fsh_sel_cont_a2 << "\t" << log_fsh_sel_cont_b2 << "\t" << log_srv_sel_age1 << "\t" << log_srv_sel_a1 << "\t" << log_srv_sel_b1 << "\t" << log_srv_sel_a2 << "\t" << log_srv_sel_b2 << "\t" << sd_avg_rec_1 << "\t" << sd_recruits_1 << "\t" << sd_avg_rec_2 << "\t" << sd_recruits_2 << "\t" << sd_total_biomass << "\t" << sd_spawn_biomass << "\t" << sd_age_3_plus_biomass << "\t" << SB0 << "\t" << avgR << "\t" << avgR_CV << "\t" << (total_biomass(endyr+1) / conv_factor) << "\t" << (spawn_biomass(endyr+1) / conv_factor) << "\t" << F_ABC << "\t" << ABC << "\t" << F_OFL << "\t" << OFL << endl;
        }

        // after the burn-in period
        // if (mcmc_iter > 1000)
        // sample from all MCMC iterations
        // do MLE first
        if (mcmc_iter == 1 || (mcmc_iter % nsubsample) == 0)
        {
            random_number_generator tmp_rng(get_iter_seed(mcmc_iter));
            mcmc_rng = tmp_rng;

            // calculate values used to generate future recruitment
            mcmc_avg_rec_years = (endyr-1) - (om_rec_avg_styr + st_age) + 1;
            mcmc_avg_rec = mean(sd_recruits_1((om_rec_avg_styr+st_age),(endyr-1)));
            mcmc_avg_log_rec = mean(log(sd_recruits_1((om_rec_avg_styr+st_age),(endyr-1))));
            mcmc_std_dev_avg_log_rec = 0.0;
            for (uu = (om_rec_avg_styr+st_age) ; uu < endyr; uu++)
            {
                mcmc_std_dev_avg_log_rec += square(log(sd_recruits_1(uu)) - mcmc_avg_log_rec);
            }
            mcmc_CV_avg_log_rec = sqrt(mcmc_std_dev_avg_log_rec / (mcmc_avg_rec_years - 1.0));
            mcmc_max_rec = max(sd_recruits_1(styr,(endyr-1)));

            cout << "in calculate_projections: MCMC avg rec values (" << om_rec_avg_styr << " on):" << "\t" << mcmc_iter << "\tavg age-" << st_age << " rec = " << mcmc_avg_rec << "\tavg log rec = " << mcmc_avg_log_rec << "\tCV = " << mcmc_CV_avg_log_rec << "\tmax rec = " << mcmc_max_rec << endl;


            // calculated above
            // calculate average fishery selectivity from 1992 through endyr-1
            // NOTE:  MWD uses (endyr-5) through (endyr-1) in pk10_1.tpl
            // Sage_fsh_avg.initialize();
            // Sage_fsh_avg = 0.0;
            // for (a = st_age; a <= end_age; a++)
            // {
                // for (y = 1992; y < endyr; y++)
                // for (y = (endyr-5); y < endyr; y++)
                // {
                    // Sage_fsh_avg(a) += Sage_fsh(y,a);
                // }
            // }
            // Sage_fsh_avg /= max(Sage_fsh_avg);


            // initialize the first year of the fishery proportions-at-age
            actual_paa = value(C(endyr)(st_age,end_age) / sum(C(endyr)(st_age,end_age))) * age_trans_all;
            if (complex_flag == 0 || complex_flag == 1)
            {
                rand_paa = actual_paa;
            }
            else
            {
                rand_paa = get_multinomial_sample(actual_paa, st_age, end_age, multN_fsh(nyrs_fsh));
            }
            fsh_paa_all_proj(endyr) = rand_paa;


            // initialize numbers-at-age in first projection year; N(endyr,st_age) = R0
            N_proj = 0.0;
            N_proj(endyr+1)(st_age,end_age) = N(endyr+1)(st_age,end_age);


            // read in the vector of future catches for this MCMC iteration
            if (catch_flag == 3)
            {
                if (mcmc_iter == 1)
                {
                    catch_proj = proj_catch_array(1);
                }
                else
                {
                    catch_proj = proj_catch_array((mcmc_iter / nsubsample));
                }
            }

            if (debug_flag) cout << "in calculate_projections: starting year loop for MCMC vector " << mcmc_iter << endl;

            // *** THE LOOP FOR THE PROJECTION PERIOD ***
            // loop over all nyrs_proj into the future
            for (yy = (endyr + 1); yy <= (endyr + nyrs_proj); yy++)
            {
                // can put some variability on this later
                // M_proj(yy)(st_age,end_age) = M(st_age,end_age);
                if (BRP_M_flag == 1)
                {
                    // calculate the average M-at-age for 1977+st_age through year yy
                    for (a = st_age; a <= end_age; a++)
                    {
                        M_at_age(a) = (((endyr - (1977 + st_age) + 1) * M(a)) + sum(column(M_proj,a)(endyr+1,yy))) / ((endyr - (1977 + st_age) + 1) + (yy - (endyr + 1) + 1));
                    }
                }
                else
                {
                    // use the current M-at-age
                    M_at_age = M_proj(yy);
                }

                // estimated TAC for previous year is this year's catch
                // could add some lognormally-distributed error to this later
                if (catch_flag == 1 || catch_flag == 2)
                {
                    compose_input_data_file(yy-1);
                }

                if (yy == (endyr + 1))
                {
                    prev_spawn_biomass = spawn_biomass(endyr);
                }
                else
                {
                    prev_spawn_biomass = spawn_biomass_proj(yy-1);
                }

                if (debug_flag) cout << "in calculate_projections: " << idf.countNodes() << " vs. " << idf_proj.countNodes() << endl;

                // run estimation model if ABC will be calculated by estimation model
                if (catch_flag == 1 || catch_flag == 2)
                {
                    // write input data to file for estimator model
                    write_idf_to_estimator_file();

                    // run estimator model
                    run_estimator_model(mcmc_iter,yy);
                }

                // get results of estimator model
                TAC = parse_estimator_results(mcmc_iter,yy);

                // if (TAC >= prev_spawn_biomass)
                // {
                    // cout << "ERROR in calculate_projections: TAC (" << TAC << ") higher than spawning biomass (" << prev_spawn_biomass << ")" << endl;
                    // TAC = 1.0;
                // }

                if (catch_flag == 1)
                {
                    // estimation model was run, but want no catch applied to "true" population
                    catch_proj(yy) = 1.0;
                }
                else if (catch_flag == 4)
                {
                    // set the TAC to the 'true' ABC
                    calculate_F_ABC(yy,M_at_age);

                    catch_proj(yy) = estABC_proj(yy);
                }
                else
                {
                    // can put lognormally-distributed implementation error here
                    catch_proj(yy) = max(1.0,value(TAC));
                }

                F_proj(yy) = calculate_F_proj(yy,M_proj(yy));

                calculate_pop_params(yy);

                calculate_F_ABC(yy,M_at_age);

                write_curryr_params_to_file(mcmc_iter,yy);

                cout << "in calculate_projections: simulation " << mcmc_iter << " year " << yy << ": est catch = " << catch_proj(yy) << ", est F = " << F_proj(yy) << ", M-at-age = " << M_at_age << endl;
            }
        }

        // close files
        mcmcfile.close();
        pfile.close();
    }

    if (debug_flag) cout << "end of calculate_projections" << endl;


FUNCTION int get_iter_seed(int iter)

    int ii;
    char s[80];

    int seed = -1;
    ifstream fseed("seeds.txt", ios::in);

    // check iter bounds
    if (iter < 1 || iter > 10000)
    {
        iter = 5050;
    }

    // make sure that the file is open (not available in all C/C++ compilers)
    // assure(fseed, "seeds.txt");

    // get line iter
    fseed.seekg(0, ios::beg);
    if (iter > 1)
    {
        for (ii = 1; ii < iter; ii++)
        {
            fseed.getline(s,80,'\n');
        }
    }
    fseed.getline(s,80);
    seed = atoi(s);

    // close file
    fseed.close();

    // cout << "seed " << iter << " is " << seed << ", s = " << s << endl;

    if (debug_flag) cout << "end of get_iter_seed" << endl;

    return(seed);


FUNCTION void compose_input_data_file(int curryr)

    char tmpstr[MAX_IDF_LINE_LEN];
    char catstr[MAX_IDF_LINE_LEN];
    double srv5;
    int ii,uu;

    // use MWD's input data file (pk13_3.dat) as the estimator output file
    if (curryr == endyr)
    {
        // read the input file in once ONLY
        if (!readIDFfile)
        {
            // read in pk13_3.dat data file
            int nline = 0;
            int len_s;
            ifstream ifs("pk13_3.dat", ios::in);
            char s[MAX_IDF_LINE_LEN];

            // go to the beginning of the file
            ifs.seekg(0);

            if (debug_flag) cout << "in compose_input_data_file: reading in DAT file..." << endl;

            while(ifs.getline(s,MAX_IDF_LINE_LEN,'\n'))
            {
                // if the line ends in a control character, remove it
                len_s = strlen(s);
                if (iscntrl(s[len_s-1])) s[len_s-1] = '\0';

                if (debug_flag) cout << "in compose_input_data_file: read in line \'" << s << "\'" << endl;

                idf.addNodeEnd(adstring(s));
                nline++;

                if (debug_flag) cout << "in compose_input_data_file: line " << nline << ", strlen = " << strlen(s) << ", \'" << s << "\'" << endl;
            }

            readIDFfile = true;
        }

        // copy the original input data file object - DOES THIS WORK HERE?
        idf_proj.copyStringLinkList(idf);
    }
    else
    {
        // set pk13_proj.dat to pk13_proj.dat + this year's data
        // idf_proj = idf_proj + new stuff;

        int rep1, rep2, rep3, rep4, rep5, rep6, rep7, rep8, rep9, repA, repB, repC, repD, repE, repF;
        int app1, app2, app3, app4, app5, app6, app7, app8, app9, appA, appB, appC, appD, appE, appF, appG, appH;
        int appI, appJ, appK, appL, appM, appN, appO, appP, appQ, appR, appS, appT, appU, appV, appW, appX, appY;
        int abf1, abf2, abf3, abf4, abf5, abf6, abf7, abf8, abf9, abfA, abfB, abfC, abfD, abfE, abfF, abfG, abfH;
        int fshpaa, srvpaa, srvyrs;
        dvector vecpaa(st_age,end_age);

        // NOTE (v 0.27):  not doing any length data yet, just FISHERY and SURVEY 1
        // projecting fishery catch & CV, catch p-a-a, w-a-a,
        // survey 1 (Dyson survey) indices & CV, p-a-a, w-a-a,
        // survey 2 (bottom trawl) indices & CV, p-a-a, w-a-a BIENNIALLY,
        // survey 3 (ADF&G coastal survey) indices & CV, p-a-a, w-a-a

        // REPLACE: find places to replace strings after this string
        rep1 = idf_proj.findNodeNum("#Starting year ending year");
        rep2 = idf_proj.findNodeNum("#Number of fishery age comps");
        rep3 = idf_proj.findNodeNum("#Transition standard deviations for random walk (endyr - styr) year for change");
        rep4 = idf_proj.findNodeNum("#Number of Dyson surveys");
        rep5 = idf_proj.findNodeNum("#Survey 1 number of age comps");
        rep6 = idf_proj.findNodeNum("#Number of survey 2");
        rep7 = idf_proj.findNodeNum("#Survey 2 number of age comps");
        rep8 = idf_proj.findNodeNum("#Number of survey 3");
        rep9 = idf_proj.findNodeNum("#Survey 3 number of age comps");
        // repA = idf_proj.findNodeNum(adstring("#Number of survey 5"));
        // repB = idf_proj.findNodeNum(adstring(""));
        // repC = idf_proj.findNodeNum(adstring(""));
        // repD = idf_proj.findNodeNum(adstring(""));
        // repE = idf_proj.findNodeNum(adstring(""));
        // repF = idf_proj.findNodeNum(adstring(""));

        // cout << "in compose_input_data_file: REPLACE\t" << rep1 << "\t" << rep2 << "\t" << rep3 << "\t" << rep4 << "\t" << rep5 << "\t" << rep6 << "\t" << rep7 << "\t" << rep8 << "\t" << rep9 << endl;

        sprintf(tmpstr,"%d\t%d", int(styr), curryr);
        idf_proj.replaceNodeString(rep1+1,tmpstr);

        // annual data
        fshpaa = int(nyrs_fsh) + (curryr - endyr);
        sprintf(tmpstr,"%d", fshpaa);
        idf_proj.replaceNodeString(rep2+1,tmpstr);

        sprintf(tmpstr,"%f", rwlk_sd(styr));
        for (ii = styr+1; ii <= (endyr - 2); ii++)
        {
            sprintf(catstr,"\t%f", rwlk_sd(ii));
            strcat(tmpstr,catstr);
        }
        for (ii = 1; ii <= (curryr - endyr); ii++)
        {
            sprintf(catstr,"\t%f", rwlk_sd(endyr-2));
            strcat(tmpstr,catstr);
        }
        sprintf(catstr,"\t%f", rwlk_sd(endyr-1));
        strcat(tmpstr,catstr);
        idf_proj.replaceNodeString(rep3+2,tmpstr);

        // survey 1 - annual data
        srvyrs = int(nyrs_srv1_dy) + (curryr - endyr);
        sprintf(tmpstr,"%d", srvyrs);
        idf_proj.replaceNodeString(rep4+1,tmpstr);
        srvpaa = int(nyrsac_srv1) + (curryr - endyr);
        sprintf(tmpstr,"%d", srvpaa);
        idf_proj.replaceNodeString(rep5+1,tmpstr);

        // survey 2 - biennial data (in odd years)
        if (curryr % 2 == 1)
        {
            srvyrs = int(nyrs_srv2) + ((curryr - endyr + 1) / 2);
            sprintf(tmpstr,"%d", srvyrs);
            idf_proj.replaceNodeString(rep6+1,tmpstr);
            srvpaa = int(nyrsac_srv2) + ((curryr - endyr + 1) / 2);
            sprintf(tmpstr,"%d", srvpaa);
            idf_proj.replaceNodeString(rep7+1,tmpstr);
        }

        // survey 3 - annual data
        srvyrs = int(nyrs_srv3) + (curryr - endyr);
        sprintf(tmpstr,"%d", srvyrs);
        idf_proj.replaceNodeString(rep8+1,tmpstr);
        // survey 3 age comp data (in even years)
        if (curryr % 2 == 0)
        {
            srvpaa = int(nyrsac_srv3) + ((curryr - endyr + 1) / 2);
            sprintf(tmpstr,"%d", srvpaa);
            idf_proj.replaceNodeString(rep9+1,tmpstr);
        }

        // APPEND: find places to append onto strings after this string
        app1 = idf_proj.findNodeNum("#Fishery annual catches (tons) all years");
        app2 = idf_proj.findNodeNum("#Fishery annual catch CV");
        app3 = idf_proj.findNodeNum("#Fishery years with age comp");
        app4 = idf_proj.findNodeNum("#Fishery multinomial sample sizes");
        app5 = idf_proj.findNodeNum("#Fishery lower and upper accumulation ages");
        app6 = idf_proj.findNodeNum("#Transition standard deviations for random walk (endyr - styr) year for change");
        app7 = idf_proj.findNodeNum("#Years in which survey 1 occured");
        app8 = idf_proj.findNodeNum("#Survey indices for survey 1");
        app9 = idf_proj.findNodeNum("#Survey 1 CV");
        appA = idf_proj.findNodeNum("#Fraction of year to midpoint of survey 1 (have to include all yrs to get expected values for non-survey yrs)");
        appB = idf_proj.findNodeNum("#Survey 1 years for age comps");
        appC = idf_proj.findNodeNum("#Survey 1 multinomial sample sizes");
        appD = idf_proj.findNodeNum("#Survey 1 lower and upper accumulation ages");
        appE = idf_proj.findNodeNum("#Years in which survey 2 occured");
        appF = idf_proj.findNodeNum("#Survey 2 indices");
        appG = idf_proj.findNodeNum("#Survey 2 CV");
        appH = idf_proj.findNodeNum("#Fraction of year to midpoint of survey 2 (have to include all yrs to get expected values for non-survey yrs)");
        appI = idf_proj.findNodeNum("#Survey 2 years for age comps");
        appJ = idf_proj.findNodeNum("#Survey 2 multinomial sample sizes");
        appK = idf_proj.findNodeNum("#Survey 2 lower and upper accumulation ages");
        // appL = idf_proj.findNodeNum(adstring("#Fraction of year to midpoint of survey 3 (have to include all yrs to get expected values for non-survey yrs)"));
        appM = idf_proj.findNodeNum("#Years in which survey 3 occured");
        appN = idf_proj.findNodeNum("#Survey 3 indices");
        appO = idf_proj.findNodeNum("#Survey 3 CV");
        appP = idf_proj.findNodeNum("#Fraction of year to midpoint of survey 3 (have to include all yrs to get expected values for non-survey yrs)");
        appQ = idf_proj.findNodeNum("#Survey 3 years for the age comps");
        appR = idf_proj.findNodeNum("#Survey 3 multinomial samples size");
        // appS = idf_proj.findNodeNum(adstring("#Years in which survey 5 occurred (plus 1 year)"));
        // appT = idf_proj.findNodeNum(adstring("#Survey 5 indices"));
        // appU = idf_proj.findNodeNum(adstring("#Survey 5 CV"));
        // appV = idf_proj.findNodeNum(adstring(""));
        // appW = idf_proj.findNodeNum(adstring(""));
        // appX = idf_proj.findNodeNum(adstring(""));
        // appY = idf_proj.findNodeNum(adstring("#Fraction of year to midpoint of survey 6 (have to include all yrs to get expected values for non-survey yrs)"));

        // cout << "in compose_input_data_file: APPEND\t" << app1 << "\t" << app2 << "\t" << app3 << "\t" << app4 << "\t" << app5 << "\t" << app6 << "\t" << app7 << "\t" << app8 << "\t" << app9 << endl;
        // cout << "in compose_input_data_file: APPEND\t" << appA << "\t" << appB << "\t" << appC << "\t" << appD << "\t" << appE << "\t" << appF << "\t" << appG << "\t" << appH << "\t" << appI << endl;
        // cout << "in compose_input_data_file: APPEND\t" << appJ << "\t" << appK << "\t" << appM << "\t" << appN << "\t" << appO << "\t" << appP << "\t" << appQ << "\t" << appR << endl;

        sprintf(tmpstr, "\t%d", curryr);
        idf_proj.appendNodeString(app1+1,tmpstr);
        sprintf(tmpstr, "\t%f", value(catch_proj(curryr)));
        idf_proj.appendNodeString(app1+2,tmpstr);

        sprintf(tmpstr, "\t%f", cattot_log_sd(endyr));
        idf_proj.appendNodeString(app2+1,tmpstr);

        sprintf(tmpstr, "\t%d", (curryr-1));
        idf_proj.appendNodeString(app3+1,tmpstr);

        sprintf(tmpstr, "\t%f", double(multN_fsh(nyrs_fsh)));
        idf_proj.appendNodeString(app4+1,tmpstr);

        sprintf(tmpstr, "\t%d", int(ac_yng_fsh(nyrs_fsh)));
        idf_proj.appendNodeString(app5+1,tmpstr);
        sprintf(tmpstr, "\t%d", int(ac_old_fsh(nyrs_fsh)));
        idf_proj.appendNodeString(app5+2,tmpstr);

        sprintf(tmpstr, "\t%d", (curryr-1));
        idf_proj.appendNodeString(app6+1,tmpstr);

        sprintf(tmpstr, "\t%d", curryr);
        idf_proj.appendNodeString(app7+1,tmpstr);

        sprintf(tmpstr, "\t%f", value(estsrvbio_proj(1,curryr)));
        idf_proj.appendNodeString(app8+1,tmpstr);

        sprintf(tmpstr, "\t%f", double(indxsurv_log_sd1_dy(nyrs_srv1_dy)));
        idf_proj.appendNodeString(app9+1,tmpstr);

        sprintf(tmpstr, "\t%d", curryr);
        idf_proj.appendNodeString(appA+1,tmpstr);
        sprintf(tmpstr, "\t%f", yrfrct_srv1(endyr));
        idf_proj.appendNodeString(appA+2,tmpstr);

        sprintf(tmpstr, "\t%d", curryr);
        idf_proj.appendNodeString(appB+1,tmpstr);

        sprintf(tmpstr, "\t%f", double(multN_srv1(nyrsac_srv1)));
        idf_proj.appendNodeString(appC+1,tmpstr);

        sprintf(tmpstr, "\t%d", int(ac_yng_srv1(nyrsac_srv1)));
        idf_proj.appendNodeString(appD+1,tmpstr);
        sprintf(tmpstr, "\t%d", int(ac_old_srv1(nyrsac_srv1)));
        idf_proj.appendNodeString(appD+2,tmpstr);

        if (curryr % 2 == 1)
        {
            sprintf(tmpstr, "\t%d", curryr);
            idf_proj.appendNodeString(appE+1,tmpstr);

            sprintf(tmpstr, "\t%f", value(estsrvbio_proj(2,curryr)));
            idf_proj.appendNodeString(appF+1,tmpstr);

            sprintf(tmpstr, "\t%f", value(estsrvbio_proj_2_sd(curryr)));
            idf_proj.appendNodeString(appG+1,tmpstr);

            sprintf(tmpstr, "\t%d", curryr);
            idf_proj.appendNodeString(appI+1,tmpstr);

            sprintf(tmpstr, "\t%f", value(estsrvbio_proj_2_multN(curryr)));
            idf_proj.appendNodeString(appJ+1,tmpstr);

            sprintf(tmpstr, "\t%d", int(ac_yng_srv2(nyrsac_srv2)));
            idf_proj.appendNodeString(appK+1,tmpstr);
            sprintf(tmpstr, "\t%d", int(ac_old_srv2(nyrsac_srv2)));
            idf_proj.appendNodeString(appK+2,tmpstr);
        }

        sprintf(tmpstr, "\t%d", curryr);
        idf_proj.appendNodeString(appH+1,tmpstr);
        sprintf(tmpstr, "\t%f", yrfrct_srv2(endyr));
        idf_proj.appendNodeString(appH+2,tmpstr);


        sprintf(tmpstr, "\t%d", curryr);
        idf_proj.appendNodeString(appM+1,tmpstr);

        sprintf(tmpstr, "\t%f", value(estsrvbio_proj(3,curryr)));
        idf_proj.appendNodeString(appN+1,tmpstr);

        sprintf(tmpstr, "\t%f", double(indxsurv_log_sd3(nyrs_srv3)));
        idf_proj.appendNodeString(appO+1,tmpstr);

        sprintf(tmpstr, "\t%d", curryr);
        idf_proj.appendNodeString(appP+1,tmpstr);
        sprintf(tmpstr, "\t%f", yrfrct_srv3(endyr));
        idf_proj.appendNodeString(appP+2,tmpstr);

        if (curryr % 2 == 0)
        {
            sprintf(tmpstr, "\t%d", curryr);
            idf_proj.appendNodeString(appQ+1,tmpstr);

            sprintf(tmpstr, "\t%f", double(multN_srv3(nyrsac_srv3)));
            idf_proj.appendNodeString(appR+1,tmpstr);
        }


        // ADD LINE: find places to add new strings before this string, usually prop-at-age data
        abf1 = idf_proj.findNodeNum("#Fishery catch at length (proportions) (years w/length comps X number of length bins) Length transition matrix 1");
        abf2 = idf_proj.findNodeNum("#Survey 1 data (acoustic)");
        abf3 = idf_proj.findNodeNum("#Survey 1 numbers at length (proportions) Length transition matrix 3");
        abf4 = idf_proj.findNodeNum("#Survey 2 data (Bottom trawl)");
        abf5 = idf_proj.findNodeNum("#Survey 2 numbers at length (proportions) Length transition matrix 2");
        // abf6 = idf_proj.findNodeNum(adstring(""));
        // abf7 = idf_proj.findNodeNum(adstring("#Survey 3 data (Egg production)"));
        // abf8 = idf_proj.findNodeNum(adstring(""));
        // abf9 = idf_proj.findNodeNum(adstring(""));
        abfA = idf_proj.findNodeNum("#Survey 3 data (ADFG coastal)");
        abfB = idf_proj.findNodeNum("#Survey 3 numbers at length (proportions) Length transition matrix 2");
        // abfC = idf_proj.findNodeNum(adstring(""));
        // abfD = idf_proj.findNodeNum(adstring("#Survey 5 data (Mckelvey age 1 index)"));
        // abfE = idf_proj.findNodeNum(adstring(""));
        abfF = idf_proj.findNodeNum("#Age error transition (10 X 10)");
        abfG = idf_proj.findNodeNum("#Population weight at age (kg) at spawning (Use Sheklikof strait EIT survey estimates)");
        abfH = idf_proj.findNodeNum("#Old maturity at age");

        // cout << "in compose_input_data_file: ADD BEFORE\t" << abf1 << "\t" << abf2 << "\t" << abf3 << "\t" << abf4 << "\t" << abf5 << "\t" << abfA << "\t" << abfB << "\t" << abfF << "\t" << abfG << "\t" << abfH << endl;

        // add them in reverse order, so that the line numbers are correct

        // this assumes that rcrage = 1 and trmage = 10

        // isn't there a better way to do this?
        sprintf(tmpstr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", wt_spawn_proj(1), wt_spawn_proj(2), wt_spawn_proj(3), wt_spawn_proj(4), wt_spawn_proj(5), wt_spawn_proj(6), wt_spawn_proj(7), wt_spawn_proj(8), wt_spawn_proj(9), wt_spawn_proj(10));
        idf_proj.addNodeBefore(abfH-1,tmpstr);
        sprintf(tmpstr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", wt_pop_proj(1), wt_pop_proj(2), wt_pop_proj(3), wt_pop_proj(4), wt_pop_proj(5), wt_pop_proj(6), wt_pop_proj(7), wt_pop_proj(8), wt_pop_proj(9), wt_pop_proj(10));
        idf_proj.addNodeBefore(abfG-1,tmpstr);


        sprintf(tmpstr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", wt_srv3(endyr,1), wt_srv3(endyr,2), wt_srv3(endyr,3), wt_srv3(endyr,4), wt_srv3(endyr,5), wt_srv3(endyr,6), wt_srv3(endyr,7), wt_srv3(endyr,8), wt_srv3(endyr,9), wt_srv3(endyr,10));
        idf_proj.addNodeBefore(abfF-1,tmpstr);

        if (curryr % 2 == 0)
        {
            vecpaa = value(srv_paa_all_proj(3,curryr)(st_age,end_age) / sum(srv_paa_all_proj(3,curryr)(rcrage,end_age)));
            if (end_age > trmage)
            {
                for (uu = (trmage + 1); uu <= end_age; uu++)
                {
                    vecpaa(trmage) += vecpaa(uu);
                }
            }
            sprintf(tmpstr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", vecpaa(1), vecpaa(2), vecpaa(3), vecpaa(4), vecpaa(5), vecpaa(6), vecpaa(7), vecpaa(8), vecpaa(9), vecpaa(10));
            idf_proj.addNodeBefore(abfB-1,tmpstr);
        }


        sprintf(tmpstr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", wt_srv2(endyr,1), wt_srv2(endyr,2), wt_srv2(endyr,3), wt_srv2(endyr,4), wt_srv2(endyr,5), wt_srv2(endyr,6), wt_srv2(endyr,7), wt_srv2(endyr,8), wt_srv2(endyr,9), wt_srv2(endyr,10));
        idf_proj.addNodeBefore(abfA-1,tmpstr);

        if (curryr % 2 == 1)
        {
            vecpaa = value(srv_paa_all_proj(2,curryr)(st_age,end_age) / sum(srv_paa_all_proj(2,curryr)(rcrage,end_age)));
            if (end_age > trmage)
            {
                for (uu = (trmage + 1); uu <= end_age; uu++)
                {
                    vecpaa(trmage) += vecpaa(uu);
                }
            }
            // this assumes that rcrage = 1 and trmage = 10
            sprintf(tmpstr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", vecpaa(1), vecpaa(2), vecpaa(3), vecpaa(4), vecpaa(5), vecpaa(6), vecpaa(7), vecpaa(8), vecpaa(9), vecpaa(10));
            idf_proj.addNodeBefore(abf5-1,tmpstr);
        }


        sprintf(tmpstr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", wt_srv1(endyr,1), wt_srv1(endyr,2), wt_srv1(endyr,3), wt_srv1(endyr,4), wt_srv1(endyr,5), wt_srv1(endyr,6), wt_srv1(endyr,7), wt_srv1(endyr,8), wt_srv1(endyr,9), wt_srv1(endyr,10));
        idf_proj.addNodeBefore(abf4-1,tmpstr);

        vecpaa = value(srv_paa_all_proj(1,curryr)(st_age,end_age) / sum(srv_paa_all_proj(1,curryr)(rcrage,end_age)));
        if (end_age > trmage)
        {
            for (uu = (trmage + 1); uu <= end_age; uu++)
            {
                vecpaa(trmage) += vecpaa(uu);
            }
        }
        // this assumes that rcrage = 1 and trmage = 10
        sprintf(tmpstr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", vecpaa(1), vecpaa(2), vecpaa(3), vecpaa(4), vecpaa(5), vecpaa(6), vecpaa(7), vecpaa(8), vecpaa(9), vecpaa(10));
        idf_proj.addNodeBefore(abf3-1,tmpstr);


        sprintf(tmpstr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", wt_fsh_proj(1), wt_fsh_proj(2), wt_fsh_proj(3), wt_fsh_proj(4), wt_fsh_proj(5), wt_fsh_proj(6), wt_fsh_proj(7), wt_fsh_proj(8), wt_fsh_proj(9), wt_fsh_proj(10));
        idf_proj.addNodeBefore(abf2-1,tmpstr);

        vecpaa = value(fsh_paa_all_proj(curryr-1)(st_age,end_age) / sum(fsh_paa_all_proj(curryr-1)(rcrage,end_age)));
        if (end_age > trmage)
        {
            for (uu = (trmage + 1); uu <= end_age; uu++)
            {
                vecpaa(trmage) += vecpaa(uu);
            }
        }
        // this assumes that rcrage = 1 and trmage = 10
        sprintf(tmpstr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", vecpaa(1), vecpaa(2), vecpaa(3), vecpaa(4), vecpaa(5), vecpaa(6), vecpaa(7), vecpaa(8), vecpaa(9), vecpaa(10));
        idf_proj.addNodeBefore(abf1-1,tmpstr);
    }

    if (debug_flag) cout << "end of compose_input_data_file" << endl;


FUNCTION dvariable calculate_F_proj(int curryr,dvar_vector Mvec)

    RETURN_ARRAYS_INCREMENT();

    dvariable Fproj = 0.01;
    int a, ii = 0;
    dvariable lowerF, upperF, testF, catch_diff;
    dvariable est_catch = 0.0;
    dvariable estZ = 0.0;
    dvariable wt_a = 0.0;

    // range of F values to search
    upperF = 2.0;
    lowerF = 0.0;

    testF = (lowerF + upperF) / 2.0;

    // calculate estimated catch using testF
    for (a = st_age; a <= end_age; a++)
    {
        estZ = (Mvec(a) + (Sage_fsh_avg(a) * testF));

        if (a < rcrage || a > trmage)
        {
            if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
            {
                wt_a = wt_fsh_proj(trmage);
            }
            else
            {
                wt_a = waa_fsh(a);
            }
        }
        else
        {
            wt_a = wt_fsh_proj(a);
        }

        // use current year's naa (naa at start of year)
        est_catch += (wt_a * (Sage_fsh_avg(a) * testF) * N_proj(curryr,a) * (1.0 - mfexp(-estZ)) / estZ);
    }
    est_catch /= 1000.0;    // convert to metric tonnes
    catch_diff = est_catch - catch_proj(curryr);

    // solve for F by bisection

    // use curryr's naa to get curryr catch total

    while (fabs(catch_diff) > BISECT_TOL && ii < MAX_BISECT_ITER)
    {
        // cout << "calculate_F_proj: est_catch = " << est_catch << " testF = " << testF << " " << i << endl;

        if (catch_diff > 0.0)
        {
            upperF = testF;
        }
        else
        {
            lowerF = testF;
        }

        testF = (lowerF + upperF) / 2.0;

        est_catch = 0.0;
        for (a = st_age; a <= end_age; a++)
        {
            estZ = Mvec(a) + (Sage_fsh_avg(a) * testF);

            if (a < rcrage || a > trmage)
            {
                if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
                {
                    wt_a = wt_fsh_proj(trmage);
                }
                else
                {
                    wt_a = waa_fsh(a);
                }
            }
            else
            {
                wt_a = wt_fsh_proj(a);
            }

            // use current year's naa (naa at start of year)
            est_catch += (wt_a * (Sage_fsh_avg(a) * testF) * N_proj(curryr,a) * (1.0 - mfexp(-estZ)) / estZ);
        }

        est_catch /= 1000.0;    // convert to metric tonnes
        catch_diff = est_catch - catch_proj(curryr);

        ii++;
    }

    Fproj = testF;

    if (debug_flag) cout << "end of calculate_F_proj" << endl;

    RETURN_ARRAYS_DECREMENT();

    return(Fproj);


FUNCTION void calculate_pop_params(int curryr)

    int a, i, multN, sample_idx;
    dvector actual_paa(st_age,end_age);
    dvector rand_paa(st_age,end_age);
    dvariable wt_a, wt_sp_a, wt_srv_a, sd_srv, Ntmp, spbio_tmp;

    // first the Ms and Zs

    // ages covered by selectivity
    for (a = st_age; a <= end_age; a++)
    {
        Z_proj(curryr,a) = M_proj(curryr,a) + (Sage_fsh_avg(a) * F_proj(curryr));
    }

    expZ_proj = mfexp(-Z_proj);

    // survival with fishing until spawning
    expZsp_proj = mfexp(-sp_frac * Z_proj);

    // this assumes continuous fishing - NEED TO FIX
    expZfrac_proj(1,curryr) = mfexp(-yrfrct_srv1(endyr) * Z_proj(curryr));
    expZfrac_proj(2,curryr) = mfexp(-yrfrct_srv2(endyr) * Z_proj(curryr));
    expZfrac_proj(3,curryr) = mfexp(-yrfrct_srv3(endyr) * Z_proj(curryr));

    // calculate the naa at the beginning of next year (given this F [catch]) and related parameters
    for (a = (st_age+1); a <= end_age; a++)
    {
        // this assumes continuous fishing
        N_proj(curryr+1,a) = N_proj(curryr,a-1) * expZ_proj(curryr,a-1);
    }
    N_proj(curryr+1,end_age) = (N_proj(curryr,end_age) * expZ_proj(curryr,end_age)) + (N_proj(curryr,end_age-1) * expZ_proj(curryr,end_age-1));

    total_biomass_proj(curryr) = spawn_biomass_proj(curryr) = age_3_plus_biomass_proj(curryr) = 0.0;
    for (a = (st_age+1); a <= end_age; a++)
    {
        // spawning takes place after fishing seasons A and B
        // and survey 1 take place; need to fix this
        if (a < rcrage || a > trmage)
        {
            if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
            {
                wt_a    = wt_pop_proj(trmage);
                wt_sp_a = wt_spawn_proj(trmage);
            }
            else
            {
                wt_a    = waa_pop(a);
                wt_sp_a = waa_srv(1,a);
            }

            spawn_biomass_proj(curryr) += (0.5 * N_proj(curryr,a) * maa(a) * wt_sp_a * expZsp_proj(curryr,a));
            total_biomass_proj(curryr) += (N_proj(curryr,a) * wt_a);

            if (a >= 3)
            {
                age_3_plus_biomass_proj(curryr) += (N_proj(curryr,a) * wt_a);
            }
        }
        else
        {
            spawn_biomass_proj(curryr) += (0.5 * N_proj(curryr,a) * maa(a) * wt_spawn_proj(a) * expZsp_proj(curryr,a));
            total_biomass_proj(curryr) += (N_proj(curryr,a) * wt_pop_proj(a));

            if (a >= 3)
            {
                age_3_plus_biomass_proj(curryr) += (N_proj(curryr,a) * wt_pop_proj(a));
            }
        }
    }

    // calculate the number of recruits (age 1 in following year)
    spbio_tmp = (spawn_biomass_proj(curryr) / conv_factor);
    steepness = 0.0;

    phi0.initialize();
    Ntmp.initialize();
    Ntmp = 0.5;
    for (a = st_age ; a < end_age; a++)
    {
        if (a < rcrage)
        {
            wt_sp_a = waa_srv(1,a);
        }
        else if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
        {
            wt_sp_a = wt_spawn_proj(trmage);
        }
        else
        {
            wt_sp_a = wt_spawn_proj(a);
        }

        phi0 += (Ntmp * maa(a) * wt_sp_a * mfexp(-sp_frac * M_proj(curryr,a)));
        Ntmp *= mfexp(-M_proj(curryr,a));
    }
    Ntmp /= (1.0 - mfexp(-M_proj(curryr,end_age)));
    phi0 += (Ntmp * maa(end_age) * wt_sp_a * mfexp(-sp_frac * M_proj(curryr,end_age)));
    cout << "phi0 in year " << curryr << " is " << phi0 << endl;

    if (s_r_relat == 1)
    {
        // Beverton-Holt S-R relationship (gamma = -1)
        N_proj(curryr+1,st_age) = conv_factor * (4.0 * R0 * h * spbio_tmp) / ((phi0 * R0 * (1.0 - h)) + (((5.0 * h) - 1.0) * spbio_tmp));
        steepness = h;
    }
    else if (s_r_relat == 2)
    {
        // Ricker S-R relationship (gamma = 0)
        N_proj(curryr+1,st_age) = conv_factor * (1.0 / phi0) * spbio_tmp * mfexp(h * (1.0 - (spbio_tmp / (R0 * phi0))));
        steepness = mfexp(h) / (mfexp(h) + 4.0);
    }
    else if (s_r_relat == 3)
    {
        // Schaefer S-R relationship (gamma = 1)
        N_proj(curryr+1,st_age) = conv_factor * R0 * spbio_tmp * (1.0 - (h * spbio_tmp));
    }
    else
    {
        // average recruits S-R relationship
        N_proj(curryr+1,st_age) = mcmc_avg_rec;
    }

    // put in some lognormally-distributed process error
    // ERROR:  N_proj(curryr+1) *= something   vs.  N_proj(curryr+1,st_age) *= something
    if (complex_flag == 0)
    {
        // do nothing here; don't apply the process error
    }
    else
    {
        N_proj(curryr+1,st_age) *= (mfexp((mcmc_CV_avg_log_rec * randn(mcmc_rng)) - (0.5 * mcmc_CV_avg_log_rec * mcmc_CV_avg_log_rec)));
    }

    // generated recruits should not be greater than the historical maximum
    // N_proj(curryr+1,st_age) = min(N_proj(curryr+1,st_age),mcmc_max_rec);
    if (N_proj(curryr+1,st_age) > mcmc_max_rec)
    {
        N_proj(curryr+1,st_age) = mcmc_max_rec;
    }


    // calculate survey biomass, converting from kg to millions of metric tonnes
    for (a = st_age; a <= end_age; a++)
    {
        if (a < rcrage || a > trmage)
        {
            if (complex_flag == 0 || complex_flag == 1 || complex_flag == 2)
            {
                i = 1;
                if (a < rcrage)
                {
                    wt_srv_a = waa_srv(i,a);
                }
                else if (a > trmage)
                {
                    wt_srv_a = wt_srv1(endyr,trmage);
                }
                estsrvN_proj(i,curryr,a) = Sage(i+nsel_fsh,a) * N_proj(curryr,a) * expZfrac_proj(i,curryr,a);
                estsrvbio_proj(i,curryr) += (q_srv(i) * wt_srv_a * estsrvN_proj(i,curryr,a) / conv_factor);

                i = 2;
                if (a < rcrage)
                {
                    wt_srv_a = waa_srv(i,a);
                }
                else if (a > trmage)
                {
                    wt_srv_a = wt_srv2(endyr,trmage);
                }
                estsrvN_proj(i,curryr,a) = Sage(i+nsel_fsh,a) * N_proj(curryr,a) * expZfrac_proj(i,curryr,a);
                estsrvbio_proj(i,curryr) += (q_srv(i) * wt_srv_a * estsrvN_proj(i,curryr,a) / conv_factor);

                i = 3;
                if (a < rcrage)
                {
                    wt_srv_a = waa_srv(i,a);
                }
                else if (a > trmage)
                {
                    wt_srv_a = wt_srv3(endyr,trmage);
                }
                estsrvN_proj(i,curryr,a) = Sage(i+nsel_fsh,a) * N_proj(curryr,a) * expZfrac_proj(i,curryr,a);
                estsrvbio_proj(i,curryr) += (q_srv(i) * wt_srv_a * estsrvN_proj(i,curryr,a) / conv_factor);
            }
            else
            {
                for (i = 1; i <= nsel_srv; i++)
                {
                    estsrvN_proj(i,curryr,a) = Sage(i+nsel_fsh,a) * N_proj(curryr,a) * expZfrac_proj(i,curryr,a);
                    estsrvbio_proj(i,curryr) += (q_srv(i) * waa_srv(i,a) * estsrvN_proj(i,curryr,a) / conv_factor);
                }
            }
        }
        else
        {
            i = 1;
            estsrvN_proj(i,curryr,a) = Sage(i+nsel_fsh,a) * N_proj(curryr,a) * expZfrac_proj(i,curryr,a);
            estsrvbio_proj(i,curryr) += (q_srv(i) * wt_srv1(endyr,a) * estsrvN_proj(i,curryr,a) / conv_factor);

            i = 2;
            estsrvN_proj(i,curryr,a) = Sage(i+nsel_fsh,a) * N_proj(curryr,a) * expZfrac_proj(i,curryr,a);
            estsrvbio_proj(i,curryr) += (q_srv(i) * wt_srv2(endyr,a) * estsrvN_proj(i,curryr,a) / conv_factor);

            i = 3;
            estsrvN_proj(i,curryr,a) = Sage(i+nsel_fsh,a) * N_proj(curryr,a) * expZfrac_proj(i,curryr,a);
            estsrvbio_proj(i,curryr) += (q_srv(i) * wt_srv3(endyr,a) * estsrvN_proj(i,curryr,a) / conv_factor);
        }
    }

    // NEED TO PUT SOME OBSERVATION ERROR IN THE SURVEYS, but with what variance?
    if (complex_flag == 0 || complex_flag == 1)
    {
        // don't apply any observation error

        // pick random values for survey 2 CV and MN sample size
        sample_idx = get_uniform_sample_int(1,nyrs_srv2);
        sd_srv = indxsurv_log_sd2(sample_idx);
        estsrvbio_proj_2_sd(curryr) = sd_srv;
        sample_idx = get_uniform_sample_int(1,nyrsac_srv2);
        multN = multN_srv2(sample_idx);
        estsrvbio_proj_2_multN(curryr) = multN;
    }
    else
    {
        i = 1;
        sd_srv = indxsurv_log_sd1_dy(nyrs_srv1_dy);
        estsrvbio_proj(i,curryr) *= (mfexp((sd_srv * randn(mcmc_rng)) - (0.5 * sd_srv * sd_srv)));

        i = 2;
        // sd_srv = indxsurv_log_sd2(nyrs_srv2);
        sample_idx = get_uniform_sample_int(1,nyrs_srv2);
        sd_srv = indxsurv_log_sd2(sample_idx);
        estsrvbio_proj_2_sd(curryr) = sd_srv;
        estsrvbio_proj(i,curryr) *= (mfexp((sd_srv * randn(mcmc_rng)) - (0.5 * sd_srv * sd_srv)));

        i = 3;
        sd_srv = indxsurv_log_sd3(nyrs_srv3);
        estsrvbio_proj(i,curryr) *= (mfexp((sd_srv * randn(mcmc_rng)) - (0.5 * sd_srv * sd_srv)));
    }

    // calculate catch-at-age numbers
    for (a = st_age; a <= end_age; a++)
    {
        C_proj(curryr,a) = (Sage_fsh_avg(a) * F_proj(curryr)) * N_proj(curryr,a) * (1.0 - expZ_proj(curryr,a)) / Z_proj(curryr,a);
    }

    // more randomized stuff

    // what does fill_multinomial() actually do, since it returns a vector of integers?
    // NOTE: might have to write my own fill_multinomial routine
    actual_paa = value(C_proj(curryr)(st_age,end_age) / sum(C_proj(curryr)(st_age,end_age))) * age_trans_all;
    // rand_paa.fill_multinomial(rng,actual_paa);  <- this doesn't do what the ADMB documentation says it does
    if (complex_flag == 0 || complex_flag == 1)
    {
        rand_paa = actual_paa;
    }
    else
    {
        rand_paa = get_multinomial_sample(actual_paa, st_age, end_age, multN_fsh(nyrs_fsh));
    }
    fsh_paa_all_proj(curryr) = rand_paa;

    for (i = 1; i <= nsel_srv; i++)
    {
        actual_paa = value(estsrvN_proj(i,curryr)(st_age,end_age) / sum(estsrvN_proj(i,curryr)(st_age,end_age))) * age_trans_all;
        // rand_paa.fill_multinomial(rng,actual_paa);    <- this doesn't do what the ADMB documentation says it does

        if (i == 1)
        {
            multN = multN_srv1(nyrsac_srv1);
        }
        else if (i == 2)
        {
            // multN = multN_srv2(nyrsac_srv2);
            sample_idx = get_uniform_sample_int(1,nyrsac_srv2);
            multN = multN_srv2(sample_idx);
            estsrvbio_proj_2_multN(curryr) = multN;
        }
        else if (i == 3)
        {
            multN = multN_srv3(nyrsac_srv3);
        }
        else
        {
            multN = 60;
        }

        if (complex_flag == 0 || complex_flag == 1)
        {
            rand_paa = actual_paa;
        }
        else
        {
            rand_paa = get_multinomial_sample(actual_paa, st_age, end_age, multN);
        }
        srv_paa_all_proj(i,curryr) = rand_paa;
    }

    if (debug_flag) cout << "end of calculate_pop_params" << endl;


FUNCTION void calculate_F_ABC(int curryr,dvar_vector Mvec)

    int a, i;
    dvariable wt_a, wt_sp_a, spbio_tmp;
    dvariable curr_spbio, curr_age3_bio, f_target, spbio_target, sb_floor;
    dvar_vector Z_test(st_age,end_age), expZ_test(st_age,end_age), expZsp_test(st_age,end_age), Z_ABC(st_age,end_age);


    // calculate all of the biological reference points
    // F50(curryr)  = get_spr_rates(0.50, Sage_fsh_avg, curryr);
    // SB50(curryr) = SBcurr;
    F100(curryr)  = get_spr_rates(1.00, Sage_fsh_avg, curryr, Mvec);
    SB100(curryr) = SBcurr;
    F40(curryr)   = get_spr_rates(0.40, Sage_fsh_avg, curryr, Mvec);
    SB40(curryr)  = SBcurr;
    F35(curryr)   = get_spr_rates(0.35, Sage_fsh_avg, curryr, Mvec);
    SB35(curryr)  = SBcurr;
    F20(curryr)   = get_spr_rates(0.20, Sage_fsh_avg, curryr, Mvec);
    SB20(curryr)  = SBcurr;

    SBtarget(curryr) = SB40(curryr) * (F35(curryr) / F40(curryr));

    // cout << "in calculate_pop_params: F47 = " << F47(curryr) << ", SB47 = " << SB47(curryr) << ", F40 = " << F40(curryr) << ", SB40 = " << SB40(curryr) << ", F35 = " << F35(curryr) << ", SB35 = " << SB35(curryr) << endl;

    // calculate the F_ABC and ABC for the 'true' population
    // NOTE:  MWD calls the variable B40, but in the file pk1.dat, it is actually B47 (0.265)  <<<------
    // in the 2005 GOA pollock SAFE pg 68 (28 in PDF), the average recruit level is 755 million
    // and F40 = 0.276 >= F_ABC, B40 = 0.224, B0 = 0.559, F35 = F_OFL = 0.326
    Z_test.initialize();
    expZ_test.initialize();
    expZsp_test.initialize();
    Z_ABC.initialize();

    // initialize F_ABC with the fishing mortality F40
    estFABC_proj(curryr) = F40(curryr);

    // iterate 5 times to get stable values for F_ABC, ABC, and F_OFL
    for (i = 1; i <= 10; i++)
    {
        // calculate total mortality using the current F_ABC
        for (a = st_age; a <= end_age; a++)
        {
            Z_test(a) = Mvec(a) + (Sage_fsh_avg(a) * estFABC_proj(curryr));
        }

        expZ_test = mfexp(-Z_test);

        // survival with fishing until spawning
        expZsp_test = mfexp(-sp_frac * Z_test);

        // initialize current spawning biomass and age 3+ biomass
        curr_spbio = curr_age3_bio = 0.0;
        for (a = (st_age+1); a <= end_age; a++)
        {
            // spawning takes place after fishing seasons A and B
            // and survey 1 take place; need to fix this
            if (a < rcrage || a > trmage)
            {
                if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
                {
                    wt_a    = wt_pop_proj(trmage);
                    wt_sp_a = wt_spawn_proj(trmage);
                }
                else
                {
                    wt_a    = waa_pop(a);
                    wt_sp_a = waa_srv(1,a);
                }

                curr_spbio += (0.5 * N_proj(curryr,a) * maa(a) * wt_sp_a * expZsp_test(a));

                if (a >= 3)
                {
                    curr_age3_bio += (N_proj(curryr,a) * wt_a);
                }
            }
            else
            {
                curr_spbio += (0.5 * N_proj(curryr,a) * maa(a) * wt_spawn_proj(a) * expZsp_test(a));

                if (a >= 3)
                {
                    curr_age3_bio += (N_proj(curryr,a) * wt_pop_proj(a));
                }
            }
        }


        // calculate the F_ABC for the 'true' population
        f_target = F40(curryr);
        spbio_target = SBtarget(curryr);
        sb_floor = SB20(curryr);

        spbio_tmp = (curr_spbio / conv_factor);

        if (spbio_tmp < (0.25 * sb_floor))
        {
            cout << "low SB in year " << curryr << ":\tSB " << curr_spbio << "\tAge 3+ " << curr_age3_bio;
            estFABC_proj(curryr) = 0.0;
            cout << "\tFABC " << estFABC_proj(curryr) << endl;
        }
        else if (spbio_tmp < sb_floor)
        {
            // reduce fishing mortality so that the catch is 0.001 * age 3+ biomass
            cout << "low SB in year " << curryr << ":\tSB " << curr_spbio << "\tAge 3+ " << curr_age3_bio;
            estFABC_proj(curryr) = solve_for_fishing_mortality((mult_factor * curr_age3_bio),curryr);
            cout << "\tFABC " << estFABC_proj(curryr) << endl;
        }
        else if (spbio_tmp < spbio_target)
        {
            estFABC_proj(curryr) = f_target * ((spbio_tmp / spbio_target) - Tier3_alpha) / (1.0 - Tier3_alpha);
        }
        else
        {
            estFABC_proj(curryr) = f_target;
        }


        // calculate the ABC for the 'true' population
        estABC_proj(curryr) = 0.0;
        for (a = st_age; a <= end_age; a++)
        {
            Z_ABC(a) = Mvec(a) + (Sage_fsh_avg(a) * estFABC_proj(curryr));

            estABCaa_proj(curryr,a) = N_proj(curryr,a) * (1.0 - mfexp(-Z_ABC(a))) * (Sage_fsh_avg(a) * estFABC_proj(curryr)) / (Z_ABC(a));

            if (a < rcrage || a > trmage)
            {
                if (a > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
                {
                    wt_a = wt_fsh_proj(trmage);
                }
                else
                {
                    wt_a = waa_fsh(a);
                }
            }
            else
            {
                wt_a = wt_fsh_proj(a);
            }

            estABC_proj(curryr) += (estABCaa_proj(curryr,a) * wt_a);
        }
        estABC_proj(curryr) /= (1000.0);


        // calculate the F_OFL for the 'true' population
        // calculate overfishing limit fishing mortality using OFL decision rule
        f_target = F35(curryr);
        spbio_target = SB40(curryr);
        sb_floor = SB20(curryr);

        if (spbio_tmp < (0.25 * sb_floor))
        {
            estFOFL_proj(curryr) = 0.0;
        }
        else if (spbio_tmp < spbio_target)
        {
            estFOFL_proj(curryr) = f_target * ((spbio_tmp / spbio_target) - Tier3_alpha) / (1.0 - Tier3_alpha);
        }
        else
        {
            estFOFL_proj(curryr) = f_target;
        }
    }

    if (debug_flag) cout << "end of calculate_F_ABC" << endl;


FUNCTION dvariable solve_for_fishing_mortality(dvariable catch_level, int naa_yr)

    RETURN_ARRAYS_INCREMENT();

    // given a catch amount for the following year, calculate the associate fishing mortality
    int aa, ii = 0;
    dvariable upperFy, lowerFy, testFy;
    dvariable Ftmp, function_value, exp_value, wt_a;

    // default value
    Ftmp = 1.0e-3;

    // initial values for Fy search
    upperFy = 2.0;
    lowerFy = 0.0;

    // find Fy here by BISECTION method
    // function_val should be 0 at Fy (see eqn B3, Punt 1995)
    testFy = (upperFy + lowerFy) / 2.0;

    function_value = 0.0;
    for (aa = st_age; aa <= end_age; aa++)
    {
        if (aa < rcrage || aa > trmage)
        {
            if (aa > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
            {
                wt_a = wt_fsh_proj(trmage);
            }
            else
            {
                wt_a = waa_fsh(aa);
            }
        }
        else
        {
            wt_a = wt_fsh_proj(aa);
        }

        if (naa_yr == endyr)
        {
            exp_value = M(aa) + (Sage_fsh_avg(aa) * testFy);
            function_value += (wt_a * Sage_fsh_avg(aa) * testFy * N(naa_yr,aa) * (1.0 - exp(-exp_value)) / exp_value);
        }
        else
        {
            exp_value = M_proj(naa_yr,aa) + (Sage_fsh_avg(aa) * testFy);
            function_value += (wt_a * Sage_fsh_avg(aa) * testFy * N_proj(naa_yr,aa) * (1.0 - exp(-exp_value)) / exp_value);
        }
    }
    function_value /= 1000.0;       // conversion to metric tonnes
    function_value -= catch_level;

    while (fabs(function_value) > BISECT_TOL && ii < MAX_BISECT_ITER)
    {
        // cout << "iteration " << ii << " function value " << function_value << " Fy " << testFy << endl;

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
        for (aa = st_age; aa <= end_age; aa++)
        {
            if (aa < rcrage || aa > trmage)
            {
                if (aa > trmage && (complex_flag == 0 || complex_flag == 1 || complex_flag == 2))
                {
                    wt_a = wt_fsh_proj(trmage);
                }
                else
                {
                    wt_a = waa_fsh(aa);
                }
            }
            else
            {
                wt_a = wt_fsh_proj(aa);
            }

            if (naa_yr == endyr)
            {
                exp_value = M(aa) + (Sage_fsh_avg(aa) * testFy);
                function_value += (wt_a * Sage_fsh_avg(aa) * testFy * N(naa_yr,aa) * (1.0 - exp(-exp_value)) / exp_value);
            }
            else
            {
                exp_value = M_proj(naa_yr,aa) + (Sage_fsh_avg(aa) * testFy);
                function_value += (wt_a * Sage_fsh_avg(aa) * testFy * N_proj(naa_yr,aa) * (1.0 - exp(-exp_value)) / exp_value);
            }
        }
        function_value /= 1000.0;       // conversion to metric tonnes
        function_value -= catch_level;

        ii++;
    }

    Ftmp = testFy;

    if (debug_flag) cout << "end of solve_for_fishing_mortality" << endl;

    RETURN_ARRAYS_DECREMENT();

    return(Ftmp);


// Jim Ianelli's code for calculating Fxx%
// START

FUNCTION dvariable get_spr_rates(dvariable spr_percent, dvar_vector sel, int curryr, dvar_vector Mvec)

    RETURN_ARRAYS_INCREMENT();

    double df=1.e-3;
    dvariable F1, F2, F3;
    dvariable yld1, yld2, yld3;
    dvariable dyld, dyldp;

    F1.initialize();
    F1 = 0.3; // starting point for Fspr...depends on species...maybe M

    // Newton Raphson stuff to go here
    for (int ii=1; ii<=10; ii++) // arbitrary fixed intervals
    {
        F2     = F1 + df;
        F3     = F1 - df;
        yld1   = -1000. * square(log(spr_percent/spr_ratio(F1, sel, curryr, Mvec)));
        yld2   = -1000. * square(log(spr_percent/spr_ratio(F2, sel, curryr, Mvec)));
        yld3   = -1000. * square(log(spr_percent/spr_ratio(F3, sel, curryr, Mvec)));
        dyld   = (yld2 - yld3) / (2. * df);           // First derivative (to find the root of this)
        dyldp  = (yld3 - (2. * yld1) + yld2) / (df*df);  // Newton-Raph approximation 2nd deriv
        F1    -= (dyld / dyldp);
    }

    if (debug_flag) cout << "end of get_spr_rates" << endl;

    RETURN_ARRAYS_DECREMENT();

    return(F1);

FUNCTION dvariable spr_unfished(int curryr, dvar_vector Mvec)

    RETURN_ARRAYS_INCREMENT();

    int jj, yy, nCVyrs;
    dvariable Ntmp, SBtmp, wt_tmp;

    SBtmp.initialize();
    avgR.initialize();

    // avgR is the average number of age 1 recruits during the period 1978 through curryr-1
    avgR = 0.0;
    if (curryr == endyr)
    {
        for (yy = (om_hcr_styr+st_age); yy < endyr; yy++)
        {
            avgR += N(yy,st_age);
        }
        avgR /= double((endyr-1) - (om_hcr_styr+st_age) + 1);
    }
    else if (curryr == (endyr+1))
    {
        for (yy = (om_hcr_styr+st_age); yy <= endyr; yy++)
        {
            avgR += N(yy,st_age);
        }
        avgR /= double((curryr-1) - (om_hcr_styr+st_age) + 1);
    }
    else
    {
        for (yy = (om_hcr_styr+st_age); yy <= endyr; yy++)
        {
            avgR += N(yy,st_age);
        }
        for (yy = (endyr+1); yy < curryr; yy++)
        {
            avgR += N_proj(yy,st_age);
        }
        avgR /= double((curryr-1) - (om_hcr_styr+st_age) + 1);
    }
    avgR /= conv_factor;

    // cout << "in spr_unfished: avgR = " << avgR << endl;

    Ntmp = 0.5 * avgR;
    SBtmp = 0.0;
    for (jj=st_age; jj < end_age; jj++)
    {
        if (jj < rcrage)
        {
            wt_tmp = waa_srv(1,jj);
        }
        else if (jj > trmage)
        {
            wt_tmp = wt_spawn_proj(trmage);
        }
        else
        {
            wt_tmp = wt_spawn_proj(jj);
        }

        SBtmp += (Ntmp * maa(jj) * wt_tmp * mfexp(-sp_frac * Mvec(jj)));
        Ntmp  *= (mfexp(-Mvec(jj)));
    }

    // jj == end_age here
    Ntmp /= (1.0 - mfexp(-Mvec(jj)));
    SBtmp += (Ntmp * maa(jj) * wt_tmp * mfexp(-sp_frac * Mvec(jj)));

    // cout << "in spr_unfished: SBtmp = " << SBtmp << endl;

    if (debug_flag) cout << "end of spr_unfished" << endl;

    RETURN_ARRAYS_DECREMENT();

    return(SBtmp);

FUNCTION dvariable spr_ratio(dvariable trial_F, dvar_vector sel, int curryr, dvar_vector Mvec)

    RETURN_ARRAYS_INCREMENT();

    /* uses following globals:
      wt(1,nages)       wt at age (female spawning)
      natmort(1,nages)  natural mortality
      p_mature(1,nages)  proportion mature (females)
      yrfrac   Fraction of year which defines peak spawning
      phizero  Spawning biomass per recruit with no fishing...
    */

    dvar_vector Ntmp(st_age,end_age);
    dvar_vector srvtmp(st_age,end_age);
    dvar_vector Ftmp(st_age,end_age);
    dvariable wt_tmp;
    int jj=1;

    SBcurr.initialize();
    Ntmp.initialize();
    srvtmp.initialize();

    for (jj = st_age; jj <= end_age; jj++)
    {
        Ftmp(jj)   = sel(jj) * trial_F;
        srvtmp(jj) = Ftmp(jj) + Mvec(jj);
    }

    SB0 = spr_unfished(curryr, Mvec);
    phi0 = SB0 / avgR;  // SBPR
    // cout << "in spr_ratio: phi0 = " << phi0 << endl;

    jj = st_age;
    if (jj < rcrage)
    {
        wt_tmp = waa_srv(1,jj);
    }
    else if (jj > trmage)
    {
        wt_tmp = wt_spawn_proj(trmage);
    }
    else
    {
        wt_tmp = wt_spawn_proj(jj);
    }

    Ntmp(jj) = 0.5 * avgR;
    SBcurr += (Ntmp(jj) * maa(jj) * wt_tmp * mfexp(-sp_frac * srvtmp(jj)));

    for (jj=(st_age+1) ; jj < end_age; jj++)
    {
        if (jj < rcrage)
        {
            wt_tmp = waa_srv(1,jj);
        }
        else if (jj > trmage)
        {
            wt_tmp = wt_spawn_proj(trmage);
        }
        else
        {
            wt_tmp = wt_spawn_proj(jj);
        }

        Ntmp(jj) = Ntmp(jj-1) * mfexp(-srvtmp(jj-1));
        SBcurr  += (Ntmp(jj) * maa(jj) * wt_tmp * mfexp(-sp_frac * srvtmp(jj)));
    }
    Ntmp(end_age) = (Ntmp(end_age-1) * mfexp(-srvtmp(end_age-1)) / (1.0 - mfexp(-srvtmp(end_age))));
    SBcurr += (Ntmp(end_age) * maa(end_age) * wt_tmp * mfexp(-sp_frac * srvtmp(end_age)));

    if (debug_flag) cout << "end of spr_ratio" << endl;

    RETURN_ARRAYS_DECREMENT();

    return(SBcurr/SB0);

// end of Jim Ianelli's code

FUNCTION dvector get_multinomial_sample(dvector obs_p, int st_idx, int end_idx, int nsamples)

    RETURN_ARRAYS_INCREMENT();

    int i, j, nbin, ndraws;
    dvector est_p(st_idx,end_idx), add_p(st_idx,end_idx);
    double prop, unif_val;

    // initialize everything
    est_p.initialize();
    add_p.initialize();
    ndraws = 0;
    prop = 0.0;
    unif_val = 0.0;

    // number of draws used to fill the estimated proportion array
    ndraws = max((8 * (end_idx - st_idx + 1)),nsamples);

    // set up the additive proportion array
    add_p = 0.0;
    for (i = st_idx; i <= end_idx; i++)
    {
        prop += obs_p[i];
        add_p[i] = prop;
    }

    // fill estimated proportions array
    est_p = 0.0;
    for (j = 0; j < ndraws; j++)
    {
        // draw a random Uniform[0,1] random value
        unif_val = randu(mcmc_rng);

        // find bin where it goes
        nbin = -1;
        for (i = st_idx+1; i <= end_idx && nbin == -1; i++)
        {
            if (add_p[i] > unif_val && unif_val >= add_p[i-1])
            {
                nbin = i;
            }
        }

        // check if this value goes in the first bin (bin st_idx)
        if (nbin == -1 && add_p[st_idx] >= unif_val)
        {
            nbin = st_idx;
        }

        // check if this value goes in the last bin (bin end_idx)
        if (nbin == -1 && unif_val >= add_p[end_idx])
        {
            nbin = end_idx;
        }

        // if bin still hasn't been found then there is an error somewhere
        if (nbin == -1)
        {
            printf("ERROR: no bin found for %f %d %d %d\n", unif_val, st_idx, end_idx, ndraws);
            nbin = st_idx;
        }

        // increment count in this bin
        est_p[nbin] += 1;
    }

    // now normalize the estimated proportions
    est_p /= ((double)ndraws);

    if (debug_flag) cout << "end of get_multinomial_sample" << endl;

    RETURN_ARRAYS_DECREMENT();

    return(est_p);


FUNCTION int get_uniform_sample_int(int st_idx, int end_idx)

    int isample = -1;
    int ibegin, iend, idiff;

    if (st_idx < end_idx)
    {
        ibegin = st_idx;
        iend   = end_idx;
    }
    else
    {
        // if the range of integers is backwards, then reverse them
        ibegin = end_idx;
        iend   = st_idx;
    }
    idiff = iend - ibegin + 1;

    // get a random uniform number and generate a random integer within the range
    // the min() is in case the randu() returns 1
    isample = min(ibegin + (int)floor(randu(mcmc_rng) * ((double)idiff)),iend);

    if (debug_flag) cout << "end of get_uniform_sample_int" << endl;

    // cout << "get_uniform_sample_int: isample = " << isample << endl;

    return(isample);


FUNCTION void write_curryr_params_to_file(int niter, int curryr)

    ofstream tdf("opm_proj.dat", ios::out | ios::app);

    tdf << "OPM: " << niter << "\t" << curryr << "\t" << catch_proj(curryr) << "\t" << F_proj(curryr) << "\t" << estFOFL_proj(curryr) << "\t" << F40(curryr) << "\t" << SBtarget(curryr) << "\t" << estABC_proj(curryr) << "\t" << estFABC_proj(curryr) << "\t" << (total_biomass_proj(curryr) / conv_factor)  << "\t" << (spawn_biomass_proj(curryr) / conv_factor) << "\t" << (age_3_plus_biomass_proj(curryr) / conv_factor) << "\t" << N_proj(curryr) << "\t" << Z_proj(curryr) << endl;

    // close file
    tdf.close();

    if (curryr == (endyr + nyrs_proj))
    {
        int yy;

        ofstream trf("opm_proj.rep", ios::out | ios::app);

        // output parameters to compare with pk13_3.rep file contents
        trf << "Simulation\t" << niter << "\tYear\t" << curryr << endl;
        trf << "Recruits" << endl;
        trf << "Hist\t" << column(N,st_age)(styr,endyr) << "\tProj\t" << column(N_proj,st_age)(endyr+1,curryr) << endl;
        trf << "Spawning biomass" << endl;
        trf << "Hist\t" << spawn_biomass(styr,endyr) << "\tProj\t" << spawn_biomass_proj(endyr+1,curryr) << endl;
        trf << "Age 3+ biomass" << endl;
        trf << "Hist\t" << age_3_plus_biomass(styr,endyr) << "\tProj\t" << age_3_plus_biomass_proj(endyr+1,curryr) << endl;
        trf << "Survey biomass" << endl;
        trf << "1\tHist\t" << estsrvbio_bs(styr,endyr) << endl;
        trf << "1\tHist\t" << estsrvbio_ek(styr,endyr) << endl;
        trf << "1\tHist\t" << estsrvbio(1)(styr,endyr) << "\tProj\t" << estsrvbio_proj(1)(endyr+1,curryr) << endl;
        trf << "2\tHist\t" << estsrvbio(2)(styr,endyr) << "\tProj\t" << estsrvbio_proj(2)(endyr+1,curryr) << endl;
        trf << "3\tHist\t" << estsrvbio(3)(styr,endyr) << "\tProj\t" << estsrvbio_proj(3)(endyr+1,curryr) << endl;
        trf << "Fishing mortality" << endl;
        trf << "Hist\t" << Fmort(styr,endyr) << "\tProj\t" << F_proj(endyr+1,curryr) << endl;
        trf << "Fishery selectivity" << endl;
        for (yy = styr; yy <= endyr; yy++)
        {
            trf << yy << "\t" << Sage_fsh(yy) << endl;
        }
        trf << "Avg\t" << Sage_fsh_avg << endl;
        trf << "Init N proj\t" << N_proj(endyr+1)(st_age,end_age) << endl;
        trf << "Q" << endl;
        trf << q_srv_bs << "\t" << q_srv_ek << "\t" << q_srv << endl;
        trf << "Survey selectivity" << endl;
        trf << "1\t" << Sage(1+nsel_fsh) << endl;
        trf << "2\t" << Sage(2+nsel_fsh) << endl;
        trf << "3\t" << Sage(3+nsel_fsh) << endl;
        trf << "Catch applied\t" << catch_proj(endyr+1,curryr) << endl;
        trf << "F applied\t" << F_proj(endyr+1,curryr) << endl;
        trf << "Biological reference points" << endl;
        trf << "SBtarget\t" << SBtarget(endyr+1,curryr) << endl;
        trf << "SB40\t" << SB40(endyr+1,curryr) << endl;
        trf << "F40\t" << F40(endyr+1,curryr) << endl;
        trf << "SB20\t" << SB20(endyr+1,curryr) << endl;
        trf << "F ABC\t" << estFABC_proj(endyr+1,curryr) << endl;
        trf << "ABC\t" << estABC_proj(endyr+1,curryr) << endl;
        trf << "F OFL\t" << estFOFL_proj(endyr+1,curryr) << endl;
        trf << endl;

        // close file
        trf.close();
    }

    if (debug_flag) cout << "end of write_curryr_params_to_file" << endl;


FUNCTION write_idf_to_estimator_file

    int i;
    int nLines = idf_proj.countNodes();
    adstring s;

    ofstream ofs("pk13_proj.dat", ios::out);

    // cout << "in write_idf_to_estimator_file: nLines = " << nLines << endl;

    if (nLines > 0)
    {
        // this is kind of crude
        // is there an easier way (with the cursor node perhaps?)
        // to get the next node's string?
        for (i = 1; i <= nLines; i++)
        {
            s = idf_proj.getNodeString(i);
            ofs << s << endl;

            // cout << "write_idf_to_estimator_file: string " << i << " is " << s << endl;
        }
    }

    // close file
    ofs.close();

    if (debug_flag) cout << "end of write_idf_to_estimator_file" << endl;


FUNCTION void run_estimator_model(int niter, int curryr)

    int retval = 0;
    char s[MAX_IDF_LINE_LEN];

    sprintf(s,"./runest.bat %d %d",niter,curryr);

    retval = system(s);

    if (retval != 0)
    {
        cout << "run_estimator_model: Error running estimator model, error " << retval << endl;
    }

    if (debug_flag) cout << "end of run_estimator_model" << endl;


FUNCTION dvariable parse_estimator_results(int niter, int curryr)

    RETURN_ARRAYS_INCREMENT();

    dvariable TAC = 0.0, tot_bio = 0.0, sp_bio = 0.0, est_F = 0.0, est_Fxx = 0.0, est_SBxx = 0.0;;
    char s[MAX_IDF_LINE_LEN], old_sp_bio[MAX_IDF_LINE_LEN], future_rec_age2[MAX_IDF_LINE_LEN];
    char *stok = NULL;
    bool foundIt = false;

    strcpy(s,"");
    strcpy(old_sp_bio,"");
    strcpy(future_rec_age2,"0");

    // if estimation model has been run, then parse the results
    if (catch_flag == 1 || catch_flag == 2)
    {
        ifstream erf("pk13_3.rep", ios::in);

        // go to the beginning of the file
        erf.seekg(0, ios::beg);

        // look for the header for the estimated historical spawning biomass
        foundIt = false;
        while((!foundIt) && erf.getline(s,MAX_IDF_LINE_LEN,'\n'))
        {
            foundIt = (strcmp(s, "Expected spawning biomass") == 0);
        }

        // if we found the estimated historical spawning biomass, then parse it
        if (foundIt)
        {
            // get the next line with all of the data
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');
            // cout << "in parse_estimator_results: s = " << s << endl;

            if (s != NULL)
            {
                cout << "in parse_estimator_results: estimated historical spawning biomass is " << s << endl;
                strcpy(old_sp_bio, s);
            }
        }

        // look for the header for the projected TAC
        foundIt = false;
        while((!foundIt) && erf.getline(s,MAX_IDF_LINE_LEN,'\n'))
        {
            foundIt = (strcmp(s, "Total catches") == 0);
        }

        // if we found the projected TAC, then parse it
        if (foundIt)
        {
            // get the next line with all of the data
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');
            // cout << "in parse_estimator_results: s = " << s << endl;

            // the first character in a printed array in ADMB is a space,
            // so skip it and get the first token
            stok = strtok(s+1," ");
            if (stok != NULL)
            {
                cout << "in parse_estimator_results: TAC token is " << stok << endl;
                TAC = double(atof(stok));
            }
        }

        // look for the header for the projected total biomass
        foundIt = false;
        while((!foundIt) && erf.getline(s,MAX_IDF_LINE_LEN,'\n'))
        {
            // cout << "in parse_estimator_results: s = " << s << endl;
            foundIt = (strcmp(s, "Summary biomass") == 0);
        }

        // if we found the projected total biomass, then parse it
        if (foundIt)
        {
            // get the next line with all of the data
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');
            // cout << "in parse_estimator_results: s = " << s << endl;

            // the first character in a printed array in ADMB is a space,
            // so skip it and get the first token
            stok = strtok(s+1," ");
            if (stok != NULL)
            {
                cout << "in parse_estimator_results: tot_bio token is " << stok << endl;
                tot_bio = double(atof(stok));
            }
        }

        // look for the header for the projected spawning biomass
        foundIt = false;
        while((!foundIt) && erf.getline(s,MAX_IDF_LINE_LEN,'\n'))
        {
            // cout << "in parse_estimator_results: s = " << s << endl;
            foundIt = (strcmp(s, "Spawning biomass") == 0);
        }

        // if we found the projected spawning biomass, then parse it
        if (foundIt)
        {
            // get the next line with all of the data
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');
            // cout << "in parse_estimator_results: s = " << s << endl;

            // the first character in a printed array in ADMB is a space,
            // so skip it and get the first token
            stok = strtok(s+1," ");
            if (stok != NULL)
            {
                cout << "in parse_estimator_results: sp_bio token is " << stok << endl;
                sp_bio = double(atof(stok));
            }
        }

        // look for the header for the projected fishing mortality
        foundIt = false;
        while((!foundIt) && erf.getline(s,MAX_IDF_LINE_LEN,'\n'))
        {
            // cout << "in parse_estimator_results: s = " << s << endl;
            foundIt = (strcmp(s, "Fishing mortality") == 0);
        }

        // if we found the projected fishing mortality, then parse it
        if (foundIt)
        {
            // get the next line with all of the data
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');
            // cout << "in parse_estimator_results: s = " << s << endl;

            // the first character in a printed array in ADMB is a space,
            // so skip it and get the first token
            stok = strtok(s+1," ");
            if (stok != NULL)
            {
                cout << "in parse_estimator_results: est_F token is " << stok << endl;
                est_F = double(atof(stok));
            }
        }

        // look for the header for the estimated future age 2 recruits
        foundIt = false;
        while((curryr > (endyr+1)) && (!foundIt) && erf.getline(s,MAX_IDF_LINE_LEN,'\n'))
        {
            foundIt = (strcmp(s, "Future recruitment") == 0);
        }

        // if we found the estimated future age 2 recruits, then parse it (them)
        if (foundIt)
        {
            // get the next line with all of the data
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');
            // cout << "in parse_estimator_results: s = " << s << endl;

            if (s != NULL)
            {
                cout << "in parse_estimator_results: estimated future age 2 recruits is " << s << endl;
                strcpy(future_rec_age2, s);
            }
        }

        // look for the header for the decision rule values
        foundIt = false;
        while((!foundIt) && erf.getline(s,MAX_IDF_LINE_LEN,'\n'))
        {
            // cout << "in parse_estimator_results: s = " << s << endl;
            foundIt = (strcmp(s, "Decision rule values") == 0);
        }

        // if we found the projected fishing mortality, then parse it
        if (foundIt)
        {
            // skip three lines
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');

            // get the next line with all of the data for 0.47
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');
            // cout << "in parse_estimator_results: s = " << s << endl;

            // get the first token
            stok = strtok(s," \t");
            if (stok != NULL)
            {
                cout << "in parse_estimator_results: xx% token is " << stok << endl;
            }

            // get the second token
            stok = strtok(NULL," \t");
            if (stok != NULL)
            {
                cout << "in parse_estimator_results: Fxx% token is " << stok << endl;
                est_Fxx = double(atof(stok));
            }

            // get the third token
            stok = strtok(NULL," \t");
            if (stok != NULL)
            {
                cout << "in parse_estimator_results: SBxx% token is " << stok << endl;
                est_SBxx = double(atof(stok));
            }

            // get the next line with all of the data for 0.40
            erf.getline(s,MAX_IDF_LINE_LEN,'\n');
            // cout << "in parse_estimator_results: s = " << s << endl;

            // get the first token
            stok = strtok(s," \t");
            if (stok != NULL)
            {
                cout << "in parse_estimator_results: xx% token is " << stok << endl;
            }

            // get the second token
            stok = strtok(NULL," \t");
            if (stok != NULL)
            {
                cout << "in parse_estimator_results: Fxx% token is " << stok << endl;
                est_Fxx = double(atof(stok));
            }
        }

        // close file
        erf.close();
    }
    else if (catch_flag == 3)
    {
        TAC = catch_proj(curryr);
    }
    else if (catch_flag == 4)
    {
        TAC = 1.0;
    }

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    // open log file
    ofstream tdf("opm_proj.dat", ios::out | ios::app);

    tdf << asctime(timeinfo);
    tdf << "EST: " << niter << "\t" << curryr << "\t" << est_Fxx << "\t" << est_SBxx << "\t" << TAC << "\t" << est_F << "\t" << tot_bio << "\t" << sp_bio  << "\t" << old_sp_bio << "\t" << future_rec_age2 << endl;

    // close file
    tdf.close();

    if (debug_flag) cout << "end of parse_estimator_results" << endl;

    RETURN_ARRAYS_DECREMENT();

    return(TAC);


TOP_OF_MAIN_SECTION

  gradient_structure::set_NUM_DEPENDENT_VARIABLES(2000);
  gradient_structure::set_MAX_NVAR_OFFSET(2000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(400000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(4000000);
  arrmblsize = 54000000; // use instead of gradient_structure::set_ARRAY_MEMBLOCK_SIZE


RUNTIME_SECTION

    convergence_criteria 1.e-1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-5, 1.e-6, 1.e-7
    maximum_function_evaluations 1000, 1000, 1000, 2000, 2000, 2000, 5000, 5000


GLOBALS_SECTION

  #include <stdio.h>
  #include <stdlib.h>
  #include <iostream>
  #include <string.h>
//  #include <process.h>
  #include <errno.h>
  #include <time.h>

  #include <admodel.h>

    // NOTE: some of the code is a bit archaic and more C-like because the
    // Borland C++ 5.2 compiler/linker SUCKS and doesn't support std namespace classes
    // perhaps it's an include file issue.  or something else.

    // using namespace std;

    // the class components and method names are hopefully self-explanatory
    class stringNode
    {
        private:
            adstring *as;
            stringNode* prevNode;
            stringNode* nextNode;

        public:
            stringNode();
            void setPrevNode(stringNode *prev);
            stringNode* getPrevNode(void);
            void setNextNode(stringNode *next);
            stringNode* getNextNode(void);
            void setNodeString(const char *s);
            adstring getNodeString(void);
            ~stringNode();

        friend class stringLinkList;
    };

    stringNode::stringNode()
    {
        // constructor
        as       = NULL;
        prevNode = NULL;
        nextNode = NULL;
    }

    stringNode::~stringNode()
    {
        // destructor

        delete as;

        as       = NULL;
        prevNode = NULL;
        nextNode = NULL;
    }

    void stringNode::setPrevNode(stringNode *ptr)
    {
        prevNode = ptr;
    }

    stringNode* stringNode::getPrevNode(void)
    {
        return prevNode;
    }

    void stringNode::setNextNode(stringNode *ptr)
    {
        nextNode = ptr;
    }

    stringNode* stringNode::getNextNode(void)
    {
        return nextNode;
    }

    void stringNode::setNodeString(const char *s)
    {
        if (as != NULL)
        {
            delete as;
            as = NULL;
        }

        // strncpy(as,s,strlen(s));
        as = new adstring(s);
    }

    adstring stringNode::getNodeString(void)
    {
        return (*as);
    }


    // the class components and method names are hopefully self-explanatory
    class stringLinkList
    {
        private:
            stringNode *head;
            stringNode *cursor;
            int nNodes;
        public:
            stringLinkList();
            void addNodeBegin(const char *s);
            void addNodeEnd(adstring s);
            void addNodeBefore(int num, const char *s);
            void addNodeAfter(int num, const char *s);
            int findNodeNum(const char *s);
            void appendNodeString(int num, const char *sS);
            void replaceNodeString(int num, const char *rS);
            void replaceNodeString(const char *oS, const char *rS);
            stringNode* getNode(int num);
            stringNode* getNode(const char *s);
            adstring getNodeString(int num);
            void deleteNode(int num);
            void deleteNode(const char *s);
            int countNodes();
            stringLinkList(const stringLinkList& sLL);
            void copyStringLinkList(const stringLinkList& sLL);
            ~stringLinkList();
    };

    stringLinkList::stringLinkList()
    {
        // constructor
        head   = NULL;
        cursor = NULL;
        nNodes = 0;
    }

    stringLinkList::~stringLinkList()
    {
        // destructor
        stringNode *tmp;

        cursor = NULL;
        nNodes = 0;

        while (head != NULL)
        {
            tmp = head;
            head = head->nextNode;
            delete tmp;
        }
    }

    void stringLinkList::addNodeBegin(const char *s)
    {
        stringNode *newSN;

        newSN           = new stringNode;
        newSN->setNodeString(s);
        newSN->nextNode = head;
        head->prevNode  = newSN;

        head            = newSN;

        nNodes++;
    }

    void stringLinkList::addNodeEnd(adstring s)
    {
        if (head == NULL)
        {
            // empty list
            head           = new stringNode;
            head->setNodeString(s);
            head->prevNode = NULL;
            head->nextNode = NULL;
        }
        else
        {
            // search for end of list
            stringNode *holdSN, *newSN;

            holdSN = head;
            while (holdSN->nextNode != NULL)
            {
                holdSN = holdSN->nextNode;
            }

            newSN            = new stringNode;
            newSN->setNodeString(s);
            newSN->nextNode  = NULL;
            newSN->prevNode  = holdSN;
            holdSN->nextNode = newSN;
        }

        nNodes++;
    }

    void stringLinkList::addNodeBefore(int n, const char *s)
    {
        stringNode *holdSN, *afterSN, *newSN;
        int i;

        holdSN = head;
        for (i = 1; i < n; i++)
        {
            if (holdSN == NULL)
            {
                cout << "stringLinkList::addNodeBefore: Bad link number " << n << ", nNodes = " << nNodes << endl;
                return;
            }

            holdSN = holdSN->nextNode;
        }

        afterSN = holdSN->nextNode;

        newSN             = new stringNode;
        newSN->setNodeString(s);
        newSN->nextNode   = holdSN->nextNode;
        newSN->prevNode   = holdSN;

        holdSN->nextNode  = newSN;
        afterSN->prevNode = newSN;

        nNodes++;
    }

    void stringLinkList::addNodeAfter(int n, const char *s)
    {
        stringNode *holdSN, *afterSN, *newSN;
        int i;

        holdSN = head;
        for (i = 1; i < n; i++)
        {
            if (holdSN == NULL)
            {
                cout << "stringLinkList::addNodeAfter: Bad link number " << n << ", nNodes = " << nNodes << endl;
                return;
            }

            holdSN = holdSN->nextNode;
        }

        afterSN = holdSN->nextNode;

        newSN             = new stringNode;
        newSN->setNodeString(s);
        newSN->nextNode   = holdSN->nextNode;
        newSN->prevNode   = holdSN;

        holdSN->nextNode  = newSN;
        afterSN->prevNode = newSN;

        nNodes++;
    }

    int stringLinkList::findNodeNum(const char *s)
    {
        stringNode *holdSN;
        bool foundIt = false;
        int nNode = 0;
        int len_s = strlen(s);

        holdSN = head;
        while (holdSN != NULL && !foundIt)
        {
            nNode++;

            foundIt = (strncmp(holdSN->getNodeString(),s,len_s) == 0);

            holdSN = holdSN->nextNode;
        }

        return nNode;
    }

    void stringLinkList::appendNodeString(int n, const char *appS)
    {
        stringNode *currNode;
        adstring currS;

        currNode = getNode(n);

        if (currNode != NULL)
        {
            // what's wrong with this?  is this the problem?
            // pointer vs. data
            currS = currNode->getNodeString();
            // strncat(currS,appS,strlen(appS));
            currS = currS + adstring(appS);
            currNode->setNodeString(currS);
        }
        else
        {
            cout << "stringLinkList::appendNodeString: did not find node " << n << endl;
        }
    }

    void stringLinkList::replaceNodeString(int n, const char *replS)
    {
        stringNode *currNode;

        currNode = getNode(n);

        if (currNode != NULL)
        {
            currNode->setNodeString(replS);
        }
        else
        {
            cout << "stringLinkList::replaceNodeString: did not find node " << n << endl;
        }
    }

    void stringLinkList::replaceNodeString(const char *oldS, const char *replS)
    {
        stringNode *currNode;

        currNode = getNode(oldS);

        if (currNode != NULL)
        {
            currNode->setNodeString(replS);
        }
        else
        {
            cout << "stringLinkList::replaceNodeString: did not find nodeString \'" << oldS << "\'" << endl;
        }
    }

    stringNode* stringLinkList::getNode(int n)
    {
        if (n <= 0 || n > nNodes)
        {
            cout << "stringLinkList::getNode: Bad link number " << n << ", nNodes = " << nNodes << endl;
            return NULL;
        }
        else
        {
            stringNode *holdSN;
            int i;

            holdSN = head;

            // don't go to the next link if the node to get is the first one
            if (n > 1)
            {
                for (i = 1; i < n; i++)
                {
                    holdSN = holdSN->nextNode;

                    if (holdSN == NULL)
                    {
                        cout << "stringLinkList::getNode: Bad link number " << n << ", nNodes = " << nNodes << endl;
                        return NULL;
                    }
                }
            }

            return holdSN;
        }
    }

    stringNode* stringLinkList::getNode(const char *s)
    {
        if (&s == NULL)
        {
            cout << "stringLinkList::getNode: Bad link string \'" << s << "\'" << endl;
            return NULL;
        }
        else
        {
            stringNode *holdSN;
            int i;
            bool foundIt = false;
            int len_s = strlen(s);

            holdSN = head;

            // if the first node has what we're looking for, then don't search
            if (strncmp(holdSN->getNodeString(),s,len_s) == 0)
            {
                // in ADMB 6.2.1, CLASS adstring doesn't have an operator for "!="
            }
            else
            {
                for (i = 2; i <= nNodes && !foundIt; i++)
                {
                    holdSN = holdSN->nextNode;

                    if (holdSN == NULL)
                    {
                        cout << "stringLinkList::getNode: Bad link string \'" << s << "\', nNodes = " << nNodes << endl;
                        return NULL;
                    }

                    foundIt = (strncmp(holdSN->getNodeString(),s,len_s) == 0);
                }
            }

            return holdSN;
        }
    }

    adstring stringLinkList::getNodeString(int n)
    {
        if (n <= 0 || n > nNodes)
        {
            cout << "stringLinkList::getNode: Bad link number " << n << ", nNodes = " << nNodes << endl;
            return ("");
        }
        else
        {
            stringNode *holdSN;
            int i;

            holdSN = head;

            // don't go to the next link if the node to get is the first one
            if (n > 1)
            {
                for (i = 1; i < n; i++)
                {
                    holdSN = holdSN->nextNode;

                    if (holdSN == NULL)
                    {
                        cout << "stringLinkList::getNode: Bad link number " << n << ", nNodes = " << nNodes << endl;
                        return ("");
                    }
                }
            }

            return (holdSN->getNodeString());
        }
    }

    void stringLinkList::deleteNode(int n)
    {
        if (n <= 0 || n > nNodes)
        {
            cout << "stringLinkList::getNode: Bad link number " << n << ", nNodes = " << nNodes << endl;
        }
        else
        {
            stringNode *holdSN, *nextSN, *prevSN;
            int i;

            holdSN = head;

            // don't go to the next link if the node to delete is the first one
            if (n > 1)
            {
                for (i = 1; i < n; i++)
                {
                    holdSN = holdSN->nextNode;

                    if (holdSN == NULL)
                    {
                        cout << "stringLinkList::getNode: Bad link number " << n << ", nNodes = " << nNodes << endl;
                        return;
                    }
                }
            }

            prevSN = holdSN->prevNode;
            nextSN = holdSN->nextNode;
            prevSN->nextNode = nextSN;
            nextSN->prevNode = prevSN;
            delete holdSN;

            nNodes--;
        }
    }

    void stringLinkList::deleteNode(const char *s)
    {
        stringNode *holdSN, *nextSN, *prevSN;
        int len_s = strlen(s);

        holdSN = head;

        // delete the first node
        if (strncmp(holdSN->getNodeString(),s,len_s) == 0)
        {
            head = head->nextNode;
            head->prevNode = NULL;
            delete holdSN;
            nNodes--;
            return;
        }
        else
        {
            prevSN = holdSN;

            while (holdSN != NULL)
            {
                if (strncmp(holdSN->getNodeString(),s,len_s) == 0)
                {
                    nextSN = holdSN->nextNode;
                    prevSN = holdSN->prevNode;
                    prevSN->nextNode = nextSN;
                    nextSN->prevNode = prevSN;
                    delete holdSN;
                    nNodes--;
                    return;
                }

                prevSN = holdSN;
                holdSN = holdSN->nextNode;
            }
        }

        nNodes--;
    }

    int stringLinkList::countNodes()
    {
        nNodes = 0;

        if (head != NULL)
        {
            stringNode *holdSN;

            holdSN = head;

            while (holdSN != NULL)
            {
                nNodes++;
                holdSN = holdSN->nextNode;
            }
        }

        return nNodes;
    }

    stringLinkList::stringLinkList(const stringLinkList& sLL)
    {
        // copy-constructor for explicit data copy, not address copy
        // what does the compiler default copy-constructor do?

        // initialize new object
        head   = NULL;
        cursor = NULL;
        nNodes = 0;

        if (sLL.head != NULL && sLL.nNodes > 0)
        {
            stringNode *holdSN;
            int i;

            holdSN = sLL.head;

            // if there is stuff, then make copies of it
            for (i = 1; i <= sLL.nNodes; i++)
            {
                addNodeEnd(holdSN->getNodeString());
                holdSN = holdSN->nextNode;
            }
        }

        cout << "stringLinkList copy-constructor: copied " << sLL.nNodes << " nodes to " << nNodes << " nodes" << endl;
    }

    void stringLinkList::copyStringLinkList(const stringLinkList& sLL)
    {
        // just like the copy-constructor, but does it work differently?

        stringNode *holdSN, *tmpSN;
        int delNodes = 0;

        // if head != NULL, then do what is in the object destructor
        if (head != NULL)
        {
            holdSN = head;
            while (holdSN != NULL)
            {
                tmpSN = holdSN;
                holdSN = holdSN->nextNode;
                delete tmpSN;
                delNodes++;
            }
            cout << "stringLinkList::copyStringLinkList: deleted " << delNodes << " nodes" << endl;
        }

        // initialize object
        head   = NULL;
        cursor = NULL;
        nNodes = 0;

        if (sLL.head != NULL && sLL.nNodes > 0)
        {
            stringNode *holdSN;
            int i;

            holdSN = sLL.head;

            // if there is stuff, then make copies of it
            for (i = 1; i <= sLL.nNodes; i++)
            {
                addNodeEnd(holdSN->getNodeString());
                holdSN = holdSN->nextNode;
            }
        }

        cout << "stringLinkList::copyStringLinkList: copied " << sLL.nNodes << " nodes to " << nNodes << " nodes" << endl;
    }

    stringLinkList idf;
    stringLinkList idf_proj;
    bool readIDFfile = false;
    time_t rawtime;
    struct tm *timeinfo;


