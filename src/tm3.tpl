 //      Filename:  tm3.tpl
 //
 //   Description:  Atka mackerel tagging model
 //
 //       Version:  1.0
 //       Created:  1/8/2014 12:04:39 PM
 //      Revision:  none
 //      Compiler:  gcc
 //
 //        Author:  Kim Rand
 //  Organization:  Resource Ecology and Fisheries Management
 //
 //=====================================================================================
 //
 DATA_SECTION
  // |---------------------|
  // | SIMULATION CONTROLS |
  // |---------------------|
  int simflag;
  int rseed
 LOC_CALCS
    simflag = 0;
    rseed   = 0;
    int opt,on;
    if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
    {
      simflag = 1;
      rseed   = atoi(ad_comm::argv[on+1]);
    }
 END_CALCS  
  // |------------------------|
  // | DATA AND CONTROL FILES |
  // |------------------------|
  init_adstring datafile;
  !! ad_comm::change_datafile_name(datafile); 
  init_number M_in;
  init_number M_cv;
  init_int M_phase;
  number lambda_M;
  !! lambda_M = 1/(2.*M_cv*M_cv);
  init_number RF_cv;
  init_int ph_RF;
  number lambda_RF
  !! lambda_RF = 1/(2.*RF_cv*RF_cv);
  init_number MeanWt;
  init_int n_events
  // Tag releases
  init_number T;
  init_vector C(1,n_events);
  init_vector obs_tags(1,n_events);
  init_ivector repindex(1,n_events);
  init_vector ndays(1,n_events);  // Days elapsed since tagging
  ivector RF_idx(1,n_events);
  !! RF_idx.initialize(); RF_idx(1)=1;
  !! for (int i=2;i<=n_events;i++) if (ndays(i)!=0) RF_idx(i) = RF_idx(i-1)+1; else RF_idx(i) = RF_idx(i-1) ;
  int n_RF        // number of rec factors
  !! n_RF = RF_idx(n_events);
  !! cout <<RF_idx<<endl;
  init_number two // # double tag recovered with two tags
  init_number one // # double tag recovered with only one tag
  init_number prop_double_tagged;
  init_vector surv_data(1,2);
  // init_int nrep_obs;
  init_matrix rep_obs(1,n_events,1,2);
  init_number p_male_rel
  init_number sample_rel
  vector tag_rel(1,2)
  !! tag_rel(1) = p_male_rel * T; tag_rel(2) = T - tag_rel(1);
  init_vector p_male_catch(1,n_events);
  init_vector sample_catch(1,n_events);
  matrix obs_catch_sex(1,2,1,n_events);
  !! obs_catch_sex(1) = elem_prod(p_male_catch , C); obs_catch_sex(2) = C - obs_catch_sex(1);
  //  Tag recovery male, female, unsexed per event
  init_vector p_male_rec(1,n_events);
  init_vector sample_rec(1,n_events);
  matrix obs_tags_sex(1,2,1,n_events);
  !! obs_tags_sex(1) = elem_prod(p_male_rec , obs_tags); 
  !! obs_tags_sex(2) = obs_tags - obs_tags_sex(1);
  // !! cout <<obs_tags_sex<<endl;exit(1);

  init_int icheck;
  number surv_init;
  !! if(icheck!=12345) {cout <<" data file problem "<<endl; exit(1);}
  !! surv_init=surv_data(1)/sum(surv_data);
  int i
  int iday
  int mcflag
  !! mcflag=0;

INITIALIZATION_SECTION
  lnNinit 6.
  surv surv_init

PARAMETER_SECTION
  init_number lnNinit(1);
  init_number lnM(M_phase);
  init_bounded_number p_male_rel_hat(0,1.,2);
  init_bounded_number ploss(0,1);
  init_bounded_number surv(0,1,2);
  init_bounded_vector repr(1,n_events,0.05,1,1);
  init_matrix RF(1,2,1,n_RF,ph_RF)
  number prob_no_tag;
  matrix N(1,2,0,n_events);
  matrix Tags(1,2,0,n_events);
  matrix pred_tags_sex(1,2,1,n_events);
  vector Ntot_Peterson(1,n_events);
  number ER;
  number ffpen
  number movement_penalty
  vector fcomp(1,6);
  number M;
  sdreport_vector Ninit(1,2);
  sdreport_number Biomass;
  // Reporting rates for commercial and charter
  sdreport_vector RepRate(1,2);
  objective_function_value f
PRELIMINARY_CALCS_SECTION
  if( simflag )
  {
    if(!global_parfile)
    {
      cerr << "Must have a gsmac.pin file to use the -sim command line option"<<endl;
      ad_exit(1);
    }
    cout<<"|--------------------------------------------|"        <<endl;
    cout<<"|*** RUNNING SIMULATION WITH RSEED = "<<rseed<<" ***|"<<endl;
    cout<<"|--------------------------------------------|"        <<endl;
    
    simulation_model();
    //exit(1);
  }

PROCEDURE_SECTION
  initialize_params();

  do_dynamics();
  get_likelihoods();

  if(mceval_phase()|sd_phase())
    set_output();


FUNCTION initialize_params
  N.initialize();
  Tags.initialize();
  prob_no_tag = (1.-prop_double_tagged) * ploss + prop_double_tagged * ploss*ploss;
  Ninit(1)    = p_male_rel_hat      * mfexp(lnNinit);
  Ninit(2)    = (1.-p_male_rel_hat) * mfexp(lnNinit);
  N(1,0)      = 1.e6*Ninit(1);
  N(2,0)      = 1.e6*Ninit(2);
  M           = M_in * mfexp(lnM);
  Tags(1,0)   = tag_rel(1)*(1-prob_no_tag)*surv;
  Tags(2,0)   = tag_rel(2)*(1-prob_no_tag)*surv;
 
FUNCTION do_dynamics
// do_sims equals one for estimation, zero for catch simulations
  dvar_vector ERtmp(1,2);      
  for (i=1; i<=n_events; i++)
  {
    N(1,i)               = N(1,i-1)*mfexp(-ndays(i)*M/365.)*mfexp(RF(1,RF_idx(i)));
    N(2,i)               = N(2,i-1)*mfexp(-ndays(i)*M/365.)*mfexp(RF(2,RF_idx(i)));
    ERtmp(1)             = obs_catch_sex(1,i)/N(1,i);
    ERtmp(2)             = obs_catch_sex(2,i)/N(2,i);
    dvariable sr         = 1.-ERtmp(1);
    if(sr< 0.01) { ffpen = 0.; dvariable  kludge; kludge   = posfun(sr,0.01,ffpen); ERtmp    = 1.-kludge; cout <<" Exceeded population with the catch... "<< kludge<<" "<<C(i)<<" "<<N(1,i)<<endl; fcomp(6)+= 1e6*ffpen; }
    N(1,i)               -= obs_catch_sex(1,i);                            // subtract off catch from remaining stock
    N(2,i)               -= obs_catch_sex(2,i);                            // subtract off catch from remaining stock
    pred_tags_sex(1,i)   = Tags(1,i-1)*ERtmp(1)*mfexp(-ndays(i)*M/365.);  // predicted tags a function of ER and previous tag numbers
    Tags(1,i)            = Tags(1,i-1) - pred_tags_sex(1,i);           // subtract off catch from remaining stock
    pred_tags_sex(2,i)   = Tags(2,i-1)*ERtmp(2)*mfexp(-ndays(i)*M/365.);  // predicted tags a function of ER and previous tag numbers
    Tags(2,i)            = Tags(2,i-1) - pred_tags_sex(2,i);           // subtract off catch from remaining stock
    pred_tags_sex(1,i) *= repr(i);
    pred_tags_sex(2,i) *= repr(i);
  }

FUNCTION void get_likelihoods()
  fcomp.initialize();
  dvariable ptmp;
  // if (active(RF))
  {
    ptmp          = N(1,0)/(N(1,0)+N(2,0));
    fcomp(3)     -= sample_rel * (p_male_rel+1e-3) * log(ptmp+1e-3);
    fcomp(3)     -= sample_rel * (1-p_male_rel+1e-3) * log(1-ptmp+1e-3);
  }
  for (i=1; i<=n_events; i++)
  {
    // cout<< pred_tags(i,tloc,loc)<<" "<<obs_tags(i,tloc,loc)<<" "<<log(pred_tags(i,tloc,loc));
    fcomp(1)     += pred_tags_sex(1,i)-(obs_tags_sex(1,i)*log(pred_tags_sex(1,i)));
    fcomp(1)     += pred_tags_sex(2,i)-(obs_tags_sex(2,i)*log(pred_tags_sex(2,i)));
    fcomp(2)     -= (rep_obs(i,1)*log(repr(i))) +
                    (rep_obs(i,2)*(log(1.-repr(i))));
    // if (active(RF))                
    {
      ptmp          = pred_tags_sex(1,i)/(pred_tags_sex(1,i)+pred_tags_sex(2,i));
      fcomp(3)     -= sample_rec(i) * (p_male_rec(i)+1e-3) * log(ptmp+1e-3);
      fcomp(3)     -= sample_rec(i) * (1-p_male_rec(i)+1e-3) * log(1-ptmp+1e-3);

      ptmp          = N(1,i)/(N(1,i)+N(2,i));
      fcomp(3)     -= sample_catch(i) * (p_male_catch(i)+1e-3) * log(ptmp+1e-3);
      fcomp(3)     -= sample_catch(i) * (1-p_male_catch(i)+1e-3) * log(1-ptmp+1e-3);
    }
  }
  fcomp(2) -= (one*log(2.*ploss*(1.-ploss)) + two*(log((1.-ploss)*(1.-ploss))));
  fcomp(2) -= (surv_data(1)*log(surv))+ (surv_data(2)*log(1.-surv));
// Jim's addition to objective function that penalizes movements being different into and out of the area
// option to change penalty in different phases... this can be useful to detect local minima
  if (active(lnM))
    fcomp(3)  = lambda_M*square(lnM);
  if (active(RF))
  {
    fcomp(4)  = lambda_RF*norm2(RF(1));
    fcomp(4) += lambda_RF*norm2(RF(2));
  }
  f=sum(fcomp);
  // cout <<fcomp<<endl;

REPORT_SECTION
  Biomass = sum(Ninit)*MeanWt;
  int i,ii;
  report<<"NLL "<<endl<<fcomp<<" "<<f<<endl;
  report<< "Ninit"<<endl;
  report<< Ninit  <<endl;
  report<< "Biomass"<<endl;
  report<< Biomass  <<endl;
  report<< "Obs_Tags"<<endl;
  report <<obs_tags<<endl;
  report<< "Pred_Tags "<<endl;
  report <<pred_tags_sex<<endl;
  report<< "ploss"<<endl;
  report<<ploss<<endl;
  report<< "N"<<endl;
  report<<N<<endl;
  report<<"prob_no_tag "<<prob_no_tag<<endl;
  report<<"surv "<<surv<<endl;
  report<<" rep_rate "<<repr<<endl;
  report<<"DaysElapsed pred_tags_male obs_tags_male pred_tags_fem obs_tags_fem RF_male RF_fem SexRatio"<<endl;
  double dcnt=0.;
  for (int i=1;i<=n_events;i++)
  {
    dcnt+=ndays(i);
    report<<dcnt<<" "
          << pred_tags_sex(1,i) <<" "
          << obs_tags_sex(1,i)  <<" "
          << pred_tags_sex(2,i) <<" "
          << obs_tags_sex(2,i)  <<" "
          << mfexp(RF(1,RF_idx(i)))     <<" "
          << mfexp(RF(2,RF_idx(i)))     <<" "
          << N(1,i)/(N(1,i)+N(2,i)) <<" "
          // << N(1,i)<<" "
          // << N(2,i)<<" "
          // << (N(1,0)+N(2,0)) - (N(1,i)+N(2,i))<<" "
          // << MeanWt*(N(1,i)+N(2,i))<<" "
          << endl;
  }
  
FUNCTION set_output
  {
    dvector nRepRate(1,2);
    dvariable MeanN=0;
    nRepRate.initialize();
    for (int i=1;i<=n_events;i++)
    {
      RepRate(repindex(i)) += repr(i);
      nRepRate(repindex(i))++;
      MeanN += N(1,i);
      MeanN += N(2,i);
    }
    MeanN /= double(n_events);
    RepRate = elem_div(RepRate,nRepRate);
    Biomass = sum(Ninit)*MeanWt;
    ER      = sum(C)/MeanN;
  }
  if(mceval_phase())
  {
    if (mcflag)
    {
      mcout<<Biomass<<" "<< Ninit<<" "<<ploss<<" "<<RepRate<<" "<<M<<" "<<ER<<" "<<surv<<" "<<f<<" "<<fcomp<<endl;
    }
    else
    {
      mcout<<"N  p_loss  Rep_Rate_F  Ref_Rate_SF M ER  Survival  ObjFun  ObjF_Tags ObjF_Repr ObjF_M x x x"<<endl;
      mcflag=1;
    }
  }

  /**

   * @brief Simulation model
   * @details Uses many of the same routines as the assessment
   * model, over-writes the observed data in memory with simulated 
   * data.
   * 
   */
FUNCTION simulation_model
  // random number generator
  random_number_generator rng(rseed);
  do_dynamics();

  for (i=1; i<=n_events; i++)
  {
    N(1,i)               = N(1,i-1)*mfexp(-ndays(i)*M/365.)*mfexp(RF(1,RF_idx(i)));
    N(2,i)               = N(2,i-1)*mfexp(-ndays(i)*M/365.)*mfexp(RF(2,RF_idx(i)));
    dvariable ERtmp;      
    ERtmp                = C(i)/(N(1,i)+N(2,i));
    obs_catch_sex(1,i)   = value(ERtmp*N(1,i)) ;                            // subtract off catch from remaining stock
    obs_catch_sex(2,i)   = value(ERtmp*N(2,i)) ;                            // subtract off catch from remaining stock
  }

  obs_tags_sex = value(pred_tags_sex);
  /* dvector drec_dev(syr+1,nyr);
  drec_dev.fill_randn(rng);
  rec_dev = exp(logSigmaR) * drec_dev;
  */
GLOBALS_SECTION
  #include <admodel.h>
  /**
  \def ECHO(object)
  Prints name and value of \a object on echoinput %ofstream file.
  */
   #undef ECHO
   #define ECHO(object) echoinput << #object << "\n" << object << endl;
   // #define ECHO(object,text) echoinput << object << "\t" << text << endl;
 
   /**
   \def check(object)
   Prints name and value of \a object on checkfile %ofstream output file.
   */
   #define check(object) checkfile << #object << "\n" << object << endl;
   // Open output files using ofstream
   ofstream echoinput("echoinput.rep");
   ofstream checkfile("checkfile.rep");
   ofstream mcout("mcout.rep");
