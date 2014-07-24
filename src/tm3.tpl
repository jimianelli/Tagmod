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
  init_int n_events
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
  init_int icheck;
  // init_int num_events;
  // init_int nsim_days;
  // init_matrix sim_dat(1,num_events,1,4);
  int num_RFs;
 LOC_CALCS
  // count number of recruitment factors needed (one each fishing event)
  /*
   num_RFs=0;
   for (int ix=1; ix<=n_events; ix++)
     for (int ii=1; ii<=2; ii++) if(C(ix,ii)>0)
       num_RFs++;
   nsim_days=nsim_days+1;
  */
 END_CALCS 
  // matrix Csims(0,nsim_days,1,2)
  // vector sim_days(1,nsim_days)
  number surv_init;
  number initpenal;
  // int ph_RF;
 LOC_CALCS
// set up variables needed for catch simulations
  /*
   Csims.initialize();
   sim_days=1;
   for (ix=1; ix<=num_events; ix++)
   {
     int kloc=sim_dat(ix,1);
     int day1=sim_dat(ix,3);
     int day2=sim_dat(ix,4);
     for (int ii=day1; ii<=day2; ii++)
     {
       Csims(ii,kloc)=sim_dat(ix,2)/(day2-day1+1);
     }
   }
  */
// command line arguments for controlling code (override data file inputs)
   adstring infile;
   int on;
   int nopt=0;
   char * end;
   if( (on=option_match(ad_comm::argc,ad_comm::argv,"-prms",nopt)) >-1)
   {
     // movement_switch=strtod(ad_comm::argv[on+1], &end);
     // movement_penalty_in=strtod(ad_comm::argv[on+2], &end);
     // Nratio_penal=strtod(ad_comm::argv[on+3], &end);
     // est_RF=strtod(ad_comm::argv[on+4], &end);

     // cout<<" options "<<movement_switch<<" "<<movement_penalty_in<<" "<<Nratio_penal<<" "<<est_RF<<endl;
   }
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
  // RFpar 1
  surv surv_init

PARAMETER_SECTION
  init_number lnNinit(1);
  init_number lnM(-3);
  init_bounded_number ploss(0,1);
  init_bounded_number surv(0,1,2);
  init_bounded_vector repr(1,n_events,0.05,1,1);
  number prob_no_tag;
  vector N(0,n_events);
  vector Tags(0,n_events);
  vector pred_tags(1,n_events);
  vector Ntot_Peterson(1,n_events);
  number ER;
  number ffpen
  number movement_penalty
  vector fcomp(1,6);
  sdreport_number M;
  sdreport_number Ninit;
  // Reporting rates for commercial and charter
  sdreport_vector RepRate(1,2);
  objective_function_value f

PROCEDURE_SECTION
  N.initialize();
  Tags.initialize();
  fcomp.initialize();
  prob_no_tag = (1.-prop_double_tagged) * ploss + prop_double_tagged * ploss*ploss;
  Ninit       = mfexp(lnNinit);
  N(0)        = 1.e6*Ninit;
  M           = M_in * mfexp(lnM);
  Tags(0)     = T*(1-prob_no_tag)*surv;

  do_dynamics();
  get_likelihoods();

  if(mceval_phase()|sd_phase())
  {
    dvector nRepRate(1,2);
    nRepRate.initialize();
    for (int i=1;i<=n_events;i++)
    {
      RepRate(repindex(i)) += repr(i);
      nRepRate(repindex(i))++;
    }
    RepRate = elem_div(RepRate,nRepRate);
  }
  if(mceval_phase())
  {
    if (mcflag)
    {
      mout<<Ninit<<" "<<ploss<<" "<<RepRate<<" "<<M<<" "<<ER<<" "<<surv<<" "<<f<<" "<<fcomp<<endl;
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
    N(i) = N(i-1)*mfexp(-ndays(i)*M/365.);
    dvariable kc;  
    kc           = C(i);
    ER           = kc/N(i);
    dvariable sr = 1.-ER;
    if(sr< 0.01) { ffpen    = 0.; dvariable  kludge; kludge   = posfun(sr,0.01,ffpen); kc       = (1.-kludge)*N(i); ER       = 1.-kludge; cout <<" Exceeded population with the catch... "<< kludge<<" "<<C(i)<<" "<<N(i)<<endl; fcomp(6)+= 1e6*ffpen; }
    N(i)        -= kc;                                     // subtract off catch from remaining stock
    pred_tags(i) = Tags(i-1)*ER*mfexp(-ndays(i)*M/365.);   // predicted tags a function of ER and previous tag numbers
    Tags(i)      = Tags(i-1) - pred_tags(i);               // subtract off catch from remaining stock
  }

FUNCTION void get_likelihoods()
  int i,tloc,kk,irf;
  for (i=1; i<=n_events; i++)
  {
    pred_tags(i) *= repr(i);
    // cout<< pred_tags(i,tloc,loc)<<" "<<obs_tags(i,tloc,loc)<<" "<<log(pred_tags(i,tloc,loc));
    fcomp(1)     += pred_tags(i)-(obs_tags(i)*log(pred_tags(i)));
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
  fcomp(3)  = 12.5*square(lnM);
  fcomp(3) += 12.5*square(lnM);

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
  report<< "Obs_Tags Pred_Tags "<<endl;
  for (ii=1; ii<=n_events; ii++)  report<<obs_tags(ii)<<" "<<pred_tags(ii)<<endl;
   // for (ii=1; ii<=n_events; ii++)  report<<obs_tags(ii,1)<<" "<<obs_tags(ii,2)<<endl;
  report<< "ploss"<<endl;
  report<<ploss<<endl;
  report<< "ER"<<endl;
  // report<<sum(C)/sum(N(0))<<endl;
  report<< "N"<<endl;
  report<<N<<endl;
  report<<" prob no tag "<<prob_no_tag<<endl;
  report<<" surv "<<surv<<endl;
  report<<" rep rate "<<repr<<endl;
  // report<<" prob move "<<pmove<<endl;
  for (i=0; i<=n_events; i++)
  {
     report<<i<<" "<<Tags(i)<<" "<<Tags(i)<<endl;
  }
  /*
  if (last_phase())
  {
    ofstream ofs("tag.out",ios::app);
    ofs<<fcomp(6)<<" "<<*(pad_tmpname)<<" "<<movement_switch<<" "<<movement_penalty<<" "<<Nratio_penal<<" "<<est_RF<<" "<<f<<" "<<Ninit<<" "<<sum(Ninit)<<" "<<
        pmove<<" "<<N(0,1)/sum(N(0))<<" "<<N(n_events,1)/sum(N(n_events))<<" "<<Ntot_Peterson/1e6<<endl;


// column(RF,1)<<" "<<column(RF,2)<<" "<<Ntot_Peterson/1e6<<endl;

//    ofs<<*(pad_tmpname)<<" "<<movement_switch<<" "<<movement_penalty<<" "<<Nratio_penal<<" "<<est_RF<<" "<<f<<" "<<Ninit<<" "<<sum(Ninit)<<" "<<
//        pmove<<" "<<N(0,1)/sum(N(0))<<" "<<N(n_events,1)/sum(N(n_events))<<" "<<repr<<" "<<surv<<" "<<ploss<<" "<<sum(C)/sum(N(0))<<" "<<
//        Ntot_Peterson/1e6<<endl;

   }
   report<<" RF est "<<RF<<endl<<" RF par "<<RFpar<<endl;
   if( num_events >0)
   {
     report<<"simulations with specified commercial catch at i "<<endl;
     report<<"day catch_inside catch_outside number_inside number_outside "<<endl; 
     for(int isim=0; isim<=nsim_days-1; isim++)  report<<isim<<" "<<Csims(isim)<<" "<<Nsims(isim)<<endl;
   }


  */
