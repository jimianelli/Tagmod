 // ===========================
 // Atka mackerel tagging model
 // Updates to Dec/2007
 //
 // ===========================
 DATA_SECTION
  init_int movement_switch;
  init_int movement_penalty_in;
  init_int Nratio_penal;
  init_int est_RF;
  init_number M;
  init_int n_events
  init_vector T(1,2)
  init_matrix C(1,n_events,1,2)
  init_3darray obs_tags(1,n_events,1,2,1,2)
  init_imatrix repindex(1,n_events,1,2)
  init_ivector ndays(1,n_events)
  init_number two
  init_number one
  init_vector prop_double_tagged(1,2);
  init_vector surv_data(1,2);
  init_ivector nrep_obs(1,2);
  init_3darray rep_obs(1,2,1,nrep_obs,1,2);
  init_matrix RF_in(1,n_events,1,2)
  init_int icheck;
  init_int num_events;
  init_int nsim_days;
  init_matrix sim_dat(1,num_events,1,4);
  int num_RFs;
 LOC_CALCS
// count number of recruitment factors needed (one each fishing event)
   num_RFs=0;
   for (int ix=1; ix<=n_events; ix++)
     for (int ii=1; ii<=2; ii++) if(C(ix,ii)>0)num_RFs++;
   nsim_days=nsim_days+1;
 END_CALCS 
  matrix Csims(0,nsim_days,1,2)
  vector sim_days(1,nsim_days)
  number surv_init;
  number initpenal;
  int ph_RF;
 LOC_CALCS
// set up variables needed for catch simulations
   Csims.initialize();
   sim_days=1;
   for (int ix=1; ix<=num_events; ix++)
   {
     int kloc=sim_dat(ix,1);
     int day1=sim_dat(ix,3);
     int day2=sim_dat(ix,4);
     for (int ii=day1; ii<=day2; ii++)
       Csims(ii,kloc)=sim_dat(ix,2)/(day2-day1+1);
   }
// command line arguments for controlling code (override data file inputs)
   adstring infile;
   int on;
   int nopt=0;
   char * end;
   if( (on=option_match(ad_comm::argc,ad_comm::argv,"-prms",nopt)) >-1)
   {
     movement_switch=strtod(ad_comm::argv[on+1], &end);
     movement_penalty_in=strtod(ad_comm::argv[on+2], &end);
     Nratio_penal=strtod(ad_comm::argv[on+3], &end);
     est_RF=strtod(ad_comm::argv[on+4], &end);

     cout<<" options "<<movement_switch<<" "<<movement_penalty_in<<" "<<Nratio_penal<<" "<<est_RF<<endl;
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
   ph_RF=-1;
   if(est_RF==1)ph_RF=1;

 END_CALCS
 !!CLASS ofstream mout("mcout.rep");
 !!CLASS adstring tmpname(infile);

INITIALIZATION_SECTION
  Ninit 20
  pmove 0.001
  surv surv_init
  RFpar 1

PARAMETER_SECTION
// added RF parameters
// if not estimated they are fixed at values in data file
  init_bounded_vector Ninit(1,2,0,3200,1);
  init_bounded_vector pmove(1,2,0.,1,2)
  init_bounded_number ploss(0,1,1);
  init_bounded_number surv(0,1,3);
  init_bounded_vector repr(1,2,0.05,1,1);
  init_bounded_vector RFpar(1,num_RFs,0.5,15,ph_RF);

  vector prob_no_tag(1,2);
  matrix RF(1,n_events,1,2);
  matrix N(0,n_events,1,2);
  matrix Nsims(0,nsim_days,1,2);
  3darray Tags(0,n_events,1,2,1,2);
  3darray pred_tags(1,n_events,1,2,1,2);
  sdreport_vector Ntot_Peterson(1,n_events);
  number reprate
  number ER
  number ffpen
  number movement_penalty
  vector fcomp(1,6);
  likeprof_number tot;
  objective_function_value f;

PRELIMINARY_CALCS_SECTION
   int ii,jj;
   RF=1;
   for (ii=1; ii<=n_events; ii++)
     for (jj=1; jj<=2; jj++) 
       if(RF_in(ii,jj)==0)
         RF_in(ii,jj)=1;
         
   RF=RF_in;

PROCEDURE_SECTION
  int time, tloc, loc, kk,iday;
  N.initialize();
  Tags.initialize();
  fcomp=0.;
  prob_no_tag=(1.-prop_double_tagged)*(ploss)+prop_double_tagged*
         (ploss*ploss);

  if(active(RFpar))
  {
    int irf=1;
    for (time=1; time<=n_events; time++)
    {
      for (loc=1; loc<=2; loc++)
      {
        if(C(time,loc)>0)RF(time,loc)=RFpar(irf++);  // get the RFpar for each fishing event
      }
    } 
  }

  N(0) = 1.e6*Ninit;
  tot  = sum(N(0));
  for (tloc=1; tloc<=2; tloc++) 
    Tags(0,tloc,tloc) = T(tloc)*(1-prob_no_tag(tloc))*surv;

  do_dynamics(1,n_events,ndays,C,N);

  get_likelihoods();

  if(mceval_phase())
   mout<<Ninit<<" "<<sum(Ninit)<<" "<<mfexp(-1*pmove)<<" "<<ploss<<" "<<repr<<" "<<surv<<" "<<RFpar<<" "<<f<<" "<<fcomp<<endl;

  if(last_phase())
  { 
    Nsims.initialize();
    get_Peterson();
    if( num_events >0)
    {
         Nsims(0)=1e6*Ninit;
         do_dynamics(0,nsim_days,sim_days,Csims,Nsims);  // to do the catch simulations
    }
  }



FUNCTION void do_dynamics(CONST int& do_sims, CONST int& nsteps, CONST dvector days, CONST dmatrix Cat, dvar_matrix Ntmp)
// do_sims equals one for estimation, zero for catch simulations
  int time,loc,kk;
  for (time=1; time<=nsteps; time++)
  {
    Ntmp(time-1)*=mfexp(-1*days(time)*M/365.);
    if(do_sims>0) Tags(time-1)*=mfexp(-1*days(time)*M/365.);
    if(movement_switch==1) do_movement1(do_sims, time, days(time), Ntmp);
    else if (movement_switch==2) do_movement2(do_sims, time, days(time), Ntmp);

    for (loc=1; loc<=2; loc++)
    {
      if(loc==1)kk=2;
      if(loc==2)kk=1;
      dvariable kC;  
      if(do_sims>0) kC=Cat(time,loc)/RF(time,loc);
      if(do_sims==0) kC=Cat(time,loc);
      ER=kC/Ntmp(time,loc);
      dvariable sr=1.-ER;
      if(sr< 0.01)
      {
         ffpen=0.;
         dvariable kludge=posfun(sr,0.01,ffpen);
         fcomp(6)+=1e6*ffpen;
         kC=(1.-kludge)*Ntmp(time,loc);
         ER=1.-kludge;
         cout <<" kludge "<< kludge<<" "<<Cat(time,loc)<<" "<<Ntmp(time,loc)<<endl;
      }
      Ntmp(time,loc)-=kC;
      if (do_sims>0)
      {
        pred_tags(time,loc,loc)=Tags(time,loc,loc)*ER;
        pred_tags(time,kk,loc)=Tags(time,kk,loc)*ER;
        Tags(time,loc,loc)-=pred_tags(time,loc,loc);
        Tags(time,kk,loc)-=pred_tags(time,kk,loc);
      }
    }
   }

FUNCTION void do_movement1(CONST int& do_sims, CONST int& time, CONST int& daytmp, dvar_matrix Ntmp)
// fish movement occurs once prior to each fishing event 
  int loc,tloc,kk;
    for (loc=1; loc<=2; loc++)
    {
      if(loc==1)kk=2;
      if(loc==2)kk=1;
      Ntmp(time,loc)+=Ntmp(time-1,loc)*mfexp(-1.*daytmp*pmove(loc));
      Ntmp(time,kk)+=Ntmp(time-1,loc)*(1.-mfexp(-1.*daytmp*pmove(loc)));
      if(do_sims>0)
      {
        for (tloc=1; tloc<=2; tloc++)
        {
          Tags(time,tloc,loc)+=Tags(time-1,tloc,loc)*
            mfexp(-1.*daytmp*pmove(loc));
          Tags(time,tloc,kk)+=Tags(time-1,tloc,loc)*
           (1.- mfexp(-1.*daytmp*pmove(loc)));
        }
      }
    }

FUNCTION void do_movement2(CONST int& do_sims, CONST int& time, CONST int& daytmp, dvar_matrix Ntmp)
// fish movement is daily
  int loc,tloc,kk,iday;  
    dvariable pstay1=mfexp(-1*pmove(1));
    dvariable pstay2=mfexp(-1*pmove(2));
    dvariable p1=Ntmp(time-1,1);
    dvariable p2=Ntmp(time-1,2);
    dvariable p1tmp;
    for (iday=1; iday<=daytmp; iday++)
      {
         p1tmp= p1*pstay1 + p2*(1-pstay2);
         p2= p2*pstay2 + p1*(1-pstay1);
         p1=p1tmp;
      }
    Ntmp(time,1)=p1;
    Ntmp(time,2)=p2;
    if(do_sims>0)
    {
      for (tloc=1; tloc<=2; tloc++)
      {
        if(tloc==1)kk=2;
        if(tloc==2)kk=1;
        dvariable pstay1 =mfexp(-1*pmove(tloc));
        dvariable pstay2 =mfexp(-1*pmove(kk));
        dvariable p1     =Tags(time-1,tloc,tloc);
        dvariable p2     =Tags(time-1,tloc,kk);
        for (iday=1; iday<=daytmp; iday++)
        {
           p1tmp= p1*pstay1 + p2*(1-pstay2);
           p2= p2*pstay2 + p1*(1-pstay1);
           p1=p1tmp;
        }
        Tags(time,tloc,tloc)=p1;
        Tags(time,tloc,kk)=p2;
      }
    }

FUNCTION void get_likelihoods()
  int time,tloc,loc,kk,irf;
  for (time=1; time<=n_events; time++)
  {
    for (loc=1; loc<=2; loc++)
    {
      if(C(time,loc)> 0.) 
      {
         for (tloc=1; tloc<=2; tloc++)
         {
            reprate                   = repr(repindex(time,loc));
            pred_tags(time,tloc,loc) *= reprate;
            fcomp(1)                 += pred_tags(time,tloc,loc) -
                                        (obs_tags(time,tloc,loc)*log(pred_tags(time,tloc,loc)));
         }
      }
    }
  }
  fcomp(2) -= (one*log(2.*ploss*(1.-ploss)) + two*(log((1.-ploss)*(1.-ploss))));
  fcomp(2) -= (surv_data(1)*log(surv))+ (surv_data(2)*log(1.-surv));
  for (kk  =1; kk<=2; kk++)
  {
    for (int j=1; j<=nrep_obs(kk); j++)
    {
      fcomp(2) -= (rep_obs(kk,j,1)*log(repr(kk))) +
                  (rep_obs(kk,j,2)*(log(1.-repr(kk))));
    }
  }
// Jim's addition to objective function that penalizes movements being different into and out of the area
// option to change penalty in different phases... this can be useful to detect local minima
  if (!last_phase()) 
    movement_penalty=initpenal;
  else 
    movement_penalty=movement_penalty_in;

  fcomp(3) += movement_penalty*square(log(pmove(1))-log(pmove(2)));

// adhoc priors for RFpars and population proportion penalty
  if(active(RFpar))
     fcomp(4)+= 0.5*sum(square(log(elem_div(RF,RF_in))/0.4));
  if(Nratio_penal > 0) 
  {
     dvariable r1=N(0,1)/sum(N(0));
     dvariable r2=N(n_events,1)/sum(N(n_events));
     fcomp(5)+=square(log(r1/r2)  );
  }
  f=sum(fcomp);


FUNCTION void get_Peterson()
// Peterson estimates by fishing event
  int time, loc;
  for (time=1; time<=n_events; time++)
  {
    int itmp=repindex(time,1);
    if(repindex(time,2)>itmp)
      itmp=repindex(time,2);
    dvariable rep = repr(itmp);
    if(sum(obs_tags(time))>0)  
      Ntot_Peterson(time) = sum(elem_prod(T,1-prob_no_tag)*surv)*sum(elem_div(C(time),RF(time)))*rep/sum(obs_tags(time));
    else  
      Ntot_Peterson(time) = 0;
  }


REPORT_SECTION
// report section should be cleaned up!
  int time,ii;
  report<<"F "<<endl<<f<<" "<<fcomp<<endl;
  report<<"switches "<<movement_switch<<" "<<movement_penalty<<" "<<Nratio_penal<<endl;

  report<< "Ninit"<<endl;
  report<< Ninit<<endl;
  report<< "Pred_Tags"<<endl;
  for (ii=1; ii<=n_events; ii++)  report<<pred_tags(ii,1)<<" "<<pred_tags(ii,2)<<endl;
  report<< "Obs_tags"<<endl;
   for (ii=1; ii<=n_events; ii++)  report<<obs_tags(ii,1)<<" "<<obs_tags(ii,2)<<endl;
  report<< "ploss"<<endl;
  report<<ploss<<endl;
  report<< "ER"<<endl;
  report<<sum(C)/sum(N(0))<<endl;
  report<< "N"<<endl;
  report<<N<<endl;
  report<<" prob no tag "<<prob_no_tag<<endl;
  report<<" surv "<<surv<<endl;
  report<<" rep rate "<<repr<<endl;
  report<<" prob move "<<pmove<<endl;
  for (time=0; time<=n_events; time++)
  {
     report<<time<<" "<<Tags(time,1)<<" "<<Tags(time,2)<<endl;
  }
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
     report<<"simulations with specified commercial catch at time "<<endl;
     report<<"day catch_inside catch_outside number_inside number_outside "<<endl; 
     for(int isim=0; isim<=nsim_days-1; isim++)  report<<isim<<" "<<Csims(isim)<<" "<<Nsims(isim)<<endl;
   }

