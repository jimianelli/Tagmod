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
  init_number MeanWt;
  init_int n_events
  // Tag releases
  init_number T;
  init_vector C(1,n_events);
  init_vector obs_tags(1,n_events);
  init_ivector repindex(1,n_events);
  init_vector ndays(1,n_events);
  init_number two // # double tag recovered with two tags
  init_number one // # double tag recovered with only one tag
  init_number prop_double_tagged;
  init_vector surv_data(1,2);
  // init_int nrep_obs;
  init_matrix rep_obs(1,n_events,1,2);
  init_number p_male_rel
  vector tag_rel(1,2)
  !! tag_rel(1) = p_male_rel * T; tag_rel(2) = T - tag_rel(1);
  init_vector p_male_catch(1,n_events);
  matrix Catch_sex(1,2,1,n_events);
  !! Catch_sex(1) = elem_prod(p_male_catch , C); Catch_sex(2) = C - Catch_sex(1);
 //  Tag recovery male, female, unsexed per event
  init_vector p_male_rec(1,n_events);
  matrix obs_tags_sex(1,2,1,n_events);
  !! obs_tags_sex(1) = elem_prod(p_male_rec , obs_tags); obs_tags_sex(2) = obs_tags - obs_tags_sex(1);
  // !! cout <<obs_tags_sex<<endl;exit(1);

  init_int ph_RF;
  init_int icheck;

 LOC_CALCS
  // count number of recruitment factors needed (one each fishing event)
 END_CALCS 
  number surv_init;
  number initpenal;
 LOC_CALCS
// command line arguments for controlling code (override data file inputs)
   adstring infile;
   int on;
   int nopt=0;
   char * end;
   initpenal=0;
   if( (on=option_match(ad_comm::argc,ad_comm::argv,"-initPen",nopt)) >-1)
   {
     initpenal=strtod(ad_comm::argv[on+1], &end);
   }

   cout<<initpenal<<endl;

   if(icheck!=12345)cout <<" data file problem "<<endl;
   if(icheck!=12345)exit(1);
   if( (on=option_match(ad_comm::argc,ad_comm::argv,"-ind",nopt)) >-1)
   {
     infile=(ad_comm::argv[on+1]);
   }
   surv_init=surv_data(1)/sum(surv_data);

 END_CALCS
 !!CLASS ofstream mout("mcout.rep");
 !!CLASS adstring tmpname(infile);
  int i
  int tloc
  int loc
  int kk
  int iday
  int mcflag
  !! mcflag=0;

INITIALIZATION_SECTION
  // Ninit 20
  lnNinit 6
  // pmove 0.001
  surv surv_init

PARAMETER_SECTION
  init_number lnNinit(1);
  init_number lnM(M_phase);
  init_bounded_number ploss(0,1);
  init_bounded_number surv(0,1,2);
  init_bounded_vector repr(1,n_events,0.05,1,1);
  init_matrix RF(1,2,1,n_events,ph_RF)
  number prob_no_tag;
  matrix N(1,2,0,n_events);
  matrix Tags(1,2,0,n_events);
  matrix pred_tags(1,2,1,n_events);
  vector Ntot_Peterson(1,n_events);
  number ER;
  number ffpen
  number movement_penalty
  vector fcomp(1,6);
  sdreport_number M;
  sdreport_vector Ninit(1,2);
  sdreport_number Biomass;
  // Reporting rates for commercial and charter
  sdreport_vector RepRate(1,2);
  objective_function_value f

PROCEDURE_SECTION
  N.initialize();
  Tags.initialize();
  fcomp.initialize();
  prob_no_tag = (1.-prop_double_tagged) * ploss + prop_double_tagged * ploss*ploss;
  Ninit(1)    = p_male_rel * mfexp(lnNinit);
  Ninit(2)    = (1.-p_male_rel) * mfexp(lnNinit);
  N(1,0)      = 1.e6*Ninit(1);
  N(2,0)      = 1.e6*Ninit(2);
  M           = M_in * mfexp(lnM);
  Tags(1,0)   = tag_rel(1)*(1-prob_no_tag)*surv;
  Tags(2,0)   = tag_rel(2)*(1-prob_no_tag)*surv;
  
  do_dynamics();
  get_likelihoods();

  if(mceval_phase()|sd_phase())
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
      mout<<Biomass<<" "<< Ninit<<" "<<ploss<<" "<<RepRate<<" "<<M<<" "<<ER<<" "<<surv<<" "<<f<<" "<<fcomp<<endl;
    }
    else
    {
      mout<<"N  p_loss  Rep_Rate_F  Ref_Rate_SF M ER  Survival  ObjFun  ObjF_Tags ObjF_Repr ObjF_M x x x"<<endl;
      mcflag=1;
    }
  }

FUNCTION do_dynamics
// do_sims equals one for estimation, zero for catch simulations
  for (i=1; i<=n_events; i++)
  {
    N(1,i)               = N(1,i-1)*mfexp(-ndays(i)*M/365.)*mfexp(RF(1,i));
    N(2,i)               = N(2,i-1)*mfexp(-ndays(i)*M/365.)*mfexp(RF(2,i));
    dvariable ERtmp;      
    ERtmp                = C(i)/(N(1,i)+N(2,i));
    dvariable sr         = 1.-ERtmp;
    if(sr< 0.01) { ffpen = 0.; dvariable  kludge; kludge   = posfun(sr,0.01,ffpen); ERtmp    = 1.-kludge; cout <<" Exceeded population with the catch... "<< kludge<<" "<<C(i)<<" "<<N(i)<<endl; fcomp(6)+= 1e6*ffpen; }
    N(1,i)               -= Catch_sex(1,i);                                     // subtract off catch from remaining stock
    N(2,i)               -= Catch_sex(2,i);                                     // subtract off catch from remaining stock
    pred_tags(1,i)       = Tags(1,i-1)*ERtmp*mfexp(-ndays(i)*M/365.);   // predicted tags a function of ER and previous tag numbers
    Tags(1,i)            = Tags(1,i-1) - pred_tags(1,i);               // subtract off catch from remaining stock
    pred_tags(2,i)       = Tags(2,i-1)*ERtmp*mfexp(-ndays(i)*M/365.);   // predicted tags a function of ER and previous tag numbers
    Tags(2,i)            = Tags(2,i-1) - pred_tags(2,i);               // subtract off catch from remaining stock
  }

FUNCTION void get_likelihoods()
  int i,tloc,kk,irf;
  for (i=1; i<=n_events; i++)
  {
    pred_tags(1,i) *= repr(i);
    pred_tags(2,i) *= repr(i);
    // cout<< pred_tags(i,tloc,loc)<<" "<<obs_tags(i,tloc,loc)<<" "<<log(pred_tags(i,tloc,loc));
    fcomp(1)     += pred_tags(1,i)-(obs_tags_sex(1,i)*log(pred_tags(1,i)));
    fcomp(1)     += pred_tags(2,i)-(obs_tags_sex(2,i)*log(pred_tags(2,i)));
    fcomp(2)     -= (rep_obs(i,1)*log(repr(i))) +
                    (rep_obs(i,2)*(log(1.-repr(i))));
  }
  fcomp(2) -= (one*log(2.*ploss*(1.-ploss)) + two*(log((1.-ploss)*(1.-ploss))));
  fcomp(2) -= (surv_data(1)*log(surv))+ (surv_data(2)*log(1.-surv));
  /*
  for (kk=1; kk<=2; kk++)
  {
    for (int j=1; j<=nrep_obs(kk); j++)
    {
      fcomp(2) -= (rep_obs(kk,j,1)*log(repr(kk))) +
                  (rep_obs(kk,j,2)*(log(1.-repr(kk))));
    }
  }
  */
// Jim's addition to objective function that penalizes movements being different into and out of the area
// option to change penalty in different phases... this can be useful to detect local minima
  fcomp(3)  = lambda_M*square(lnM);
  fcomp(4)  = norm2(RF(1));
  fcomp(4) += norm2(RF(2));

  f=sum(fcomp);

FUNCTION void get_Peterson()
// Peterson estimates by fishing event
  for (i=1; i<=n_events; i++)
  {
    int itmp=repindex(i);
    if(repindex(i)>itmp)
      itmp=repindex(i);
    dvariable rep=repr(itmp);
    if((obs_tags(i))>0)  
      Ntot_Peterson(i) = (T * (1.-prob_no_tag) *surv * C(i) * rep + 1) / (1+obs_tags(i)) -1 ;
    // T*(1-surv)*C*RepRate/(sumObs)
      // Ntot_Peterson(i)= sum(T*1-prob_no_tag*surv)*sum(elem_div(C(i),RF(i)))*rep/sum(obs_tags(i));
      // Ntot_Peterson(i) = (T * (1.-prob_no_tag) *surv * C(i) * rep ) / (obs_tags(i))  ;
    else  
      Ntot_Peterson(i)=0;
  }

 /* FUNCTION sim_catches
    Nsims.initialize();
    get_Peterson();
    if( num_events >0)
    {
         Nsims(0)=1e6*Ninit;
         do_dynamics(0,nsim_days,sim_days,Csims,Nsims);  // to do the catch simulations
    }
    */ 


REPORT_SECTION
// report section should be cleaned up!
  int i,ii;
  report<<"F "<<endl<<f<<" "<<fcomp<<endl;
  report<< "Ninit"<<endl;
  report<< Ninit  <<endl;
  report<< "Biomass"<<endl;
  report<< Biomass  <<endl;
  report<< "Obs_Tags"<<endl;
  report <<obs_tags<<endl;
  report<< "Pred_Tags "<<endl;
  report <<pred_tags<<endl;
   // for (ii=1; ii<=n_events; ii++)  report<<obs_tags(ii,1)<<" "<<obs_tags(ii,2)<<endl;
  report<< "ploss"<<endl;
  report<<ploss<<endl;
  // report<< "ER"<<endl;
  // report<<sum(C)/sum(N(0))<<endl;
  report<< "N"<<endl;
  report<<N<<endl;
  report<<"prob_no_tag "<<prob_no_tag<<endl;
  report<<"surv "<<surv<<endl;
  report<<" rep_rate "<<repr<<endl;
  report<<"DaysElapsed"<<endl;
  double dcnt=0.;
  for (int i=1;i<=n_events;i++)
  {
    dcnt+=ndays(i);
    report<<dcnt<<endl;
  }
