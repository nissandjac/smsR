// Create a file to run the TMB version of sprat
#include <TMB.hpp>
#include <iostream>


template <class Type>
vector<Type> cumsum(vector<Type> x) {
  int n = x.size();
  vector<Type> ans(n);
  ans[0] = x[0];
  for (int i = 1; i < n; i++) ans[i] = x[i] + ans[i-1];
  return ans;
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data input
  DATA_ARRAY(weca); // Weight in the beginning of the year
  DATA_ARRAY(west); // Weight in catch
  DATA_ARRAY(Surveyobs); // Surveys
  DATA_ARRAY(Catchobs); // Catch observations
//   // //
// // // // Age
  DATA_VECTOR(age); // ages
  DATA_INTEGER(nage); // Number of age groups
  DATA_INTEGER(nseason);
  DATA_INTEGER(nyears); // Number of years
  DATA_INTEGER(nsurvey); // Number of surveys
  DATA_INTEGER(recseason); // Season where recruitment occurs
  DATA_INTEGER(Fminage);
  DATA_INTEGER(Fmaxage); // Max age to fish
  DATA_INTEGER(useEffort); // use effort or estimate F
  DATA_INTEGER(useBlocks); // Use blocks for species selectivity?
  DATA_IVECTOR(CminageSeason); // Minimum age with fishing mortality per season
  DATA_IVECTOR(Qminage); // Minium ages in surveys
  DATA_IVECTOR(Qmaxage); // Maximum age in survey
  DATA_IVECTOR(Qlastage);
  DATA_IVECTOR(bidx); // Indexes for blocks of fishing mortality
  DATA_INTEGER(endFseason); // Season in last year where fishing ends
  DATA_VECTOR(isFseason);
  DATA_VECTOR(surveyEnd);
  DATA_VECTOR(surveyStart);
  DATA_VECTOR(surveySeason);
  DATA_ARRAY(power);
  DATA_ARRAY(nocatch); // Matrix sized (year x season) determines wehter F>0
  DATA_IVECTOR(Qidx);
  DATA_IARRAY(Qidx_CV); // survey catchability matrix
//   // // // Time varying arrays
//   //DATA_ARRAY(F0);
  DATA_ARRAY(M); // Natural mortality
  DATA_ARRAY(Mat); // Maturity
  DATA_ARRAY(effort); // Effort input
// //
  DATA_IARRAY(catchCV); // Ages included in catch variation
//
// // // Recruitment
  DATA_INTEGER(recmodel); // Which recruitment model?
  //DATA_SCALAR(alpha); // Alpha parameter for recruitment
  //DATA_SCALAR(beta); // Beta parameter for recruitment
  //DATA_SCALAR(logSDrec); // Recruitment deviations
  DATA_VECTOR(nllfactor); // weight of likelihood functions
  DATA_IVECTOR(estCV); // Determine which SDs are calculated or estimated
   // Parameters
  PARAMETER_VECTOR(logRin); // Recruitment
  PARAMETER_VECTOR(logNinit); // Initial distribution
  PARAMETER_VECTOR(Fyear);
  PARAMETER_ARRAY(Fseason);
  PARAMETER_ARRAY(logFage);
  PARAMETER_VECTOR(SDsurvey); // Survey variation
  PARAMETER_ARRAY(SDcatch); // Catch SD's
// //
  PARAMETER_VECTOR(logQ);
  PARAMETER(pin);
  PARAMETER(logalpha);
  PARAMETER(logbeta);
  PARAMETER(logSDrec);

//
array<Type>Catch(nage,nyears, nseason);
array<Type>CatchN(nage,nyears, nseason);
array<Type>Zsave(nage,nyears,nseason);
array<Type>Nsave(nage,nyears+1,nseason);
array<Type>survey(nage, nyears,nsurvey);
array<Type>Qsurv(nage, nsurvey);
array<Type>SDS(nage, nsurvey);
array<Type>SDC(nage, nseason);
array<Type>F0(nage, nyears, nseason);
array<Type>log_exp_pattern(nage, nseason);
array<Type>Fquarter(nage,nyears,nseason);
array<Type>Fsel(nage,nyears,nseason);
array<Type>p(nage,nsurvey);
array<Type>Fagein(nage, nyears);
array<Type>CVcatch(nage, nyears);
//array<Type>SDR_catch(nage, nyears);


vector<Type>SSB(nyears+1);
vector<Type>Rsave(nyears+1);
vector<Type>logRec(nyears+1);
vector<Type>term_logN_next(nage); // Numbers at age in the future year

vector<Type>Rin(nyears);
vector<Type>Catchtot(nyears);
// vector<Type>powerin(nsurvey);
vector<Type>Q(logQ.size());


Catch.setZero();
CatchN.setZero();
Catchtot.setZero();
Rsave.setZero();
SSB.setZero();
Zsave.setZero();
Nsave.setZero();
survey.setZero();
SDC.setZero();
SDS.setZero();
Qsurv.setZero();
Fsel.setZero();
log_exp_pattern.setZero();

// Retransform and set up for model
for(int time=0;time<nyears;time++){
      Rin(time)=exp(logRin(time));
}

Type alpha = exp(logalpha);
Type beta = exp(logbeta);
Type SDrec = exp(logSDrec);


for(int i=0;i<(logQ.size());i++){ //
  Q(i) = exp(logQ(i));
}

// Age contribution of fishing mortality

Fagein.setZero();
//
for(int i=0;i<(nage);i++){ //
  for(int j=0;j<(nyears);j++){

  if(i >= Fminage){

    if(i < Fmaxage){
      Fagein(i,j) = exp(logFage(i-Fminage, bidx(j)));
    }else{
      Fagein(i,j) = exp(logFage(Fmaxage-Fminage,bidx(j)));
      }
    }
   }
}

REPORT(Fagein)
//
// // Seasonal contribution of fishing mortality
Fquarter.setZero();
// //
if(useBlocks == 0){
  for(int j=0;j<(nyears);j++){
    for(int i=0;i<nage;i++){ // Loop over other ages (recruits excluded)
      for(int qrts=0;qrts<nseason;qrts++){ // Loop over other ages

        if(isFseason(qrts) == 1){
          if(i < CminageSeason(qrts)){
          Fquarter(i,j,qrts) = Fseason(i,qrts);

        }else{
          Fquarter(i,j,qrts) = Fseason(0,qrts);
        }
        }else{

          if(i <= CminageSeason(qrts)){
          Fquarter(i,j,qrts) = Type(1);
          }else{
          Fquarter(i,j,qrts) = Type(1)/nseason; // Where does this come from?
          }
        }
      }
    }
  }
}else{
  for(int j=0;j<(nyears);j++){
    for(int i=0;i<nage;i++){ // Loop over other ages (recruits excluded)
      for(int qrts=0;qrts<nseason;qrts++){ // Loop over other ages

        if(isFseason(qrts) == 1){
          if(i <= CminageSeason(qrts)){
          Fquarter(i,j,qrts) = Type(0.0);
          }else{
          Fquarter(i,j,qrts) = Fseason(qrts,bidx(j));
          }
        }else{

          if(i <= CminageSeason(qrts)){
          Fquarter(i,j,qrts) = Type(1);
          }else{
          Fquarter(i,j,qrts) = Type(1)/nseason; // Where does this come from?
          }
        }
      }
    }
  }
}
//
REPORT(Fquarter)
// REPORT(Fagein)
// REPORT(Fyear)
//
// // // // // // Annual contribution of fishing mortality
// F0.setZero();
// //
// //
if(useEffort == 1){
  for(int time=0;time<(nyears);time++){ // Loop over years excluding last one
    for(int i=0;i<nage;i++){ // Loop over other ages (recruits excluded)
      for(int qrts=0;qrts<nseason;qrts++){ // Loop over other seasons

          //  if(qrts != (nseason-1)){
              if(nocatch(time, qrts)>0){
              F0(i,time,qrts) = Fquarter(i,time,qrts)*Fagein(i,time)*effort(time, qrts);
              }


              if(time == (nyears-1)){
                  if(Fquarter(i, time, qrts) > 0 && Fagein(i,time)> 0){

                  log_exp_pattern(i,qrts) = log(Fquarter(i,time,qrts)*Fagein(i,time));
                }else{
                  log_exp_pattern(i, qrts) = -1000;
                }
              }

             Zsave(i,time,qrts) = F0(i,time,qrts)+M(i,time,qrts);
     }
    }
  }
}else{
  for(int time=0;time<(nyears);time++){ // Loop over years excluding last one
    for(int i=0;i<nage;i++){ // Loop over other ages (recruits excluded)
      for(int qrts=0;qrts<nseason;qrts++){ // Loop over other seasons

          if(nocatch(time, qrts)>0){
            F0(i,time,qrts) = Fquarter(i,time,qrts)*Fagein(i,time)*Fyear(time);
          }

          if(time == (nyears-1)){
                if(Fquarter(i, time, qrts) > 0 && Fagein(i,time)> 0){

                log_exp_pattern(i,qrts) = log(Fquarter(i,time,qrts)*Fagein(i,time));
              }else{
                log_exp_pattern(i, qrts) = -1000;
              }
            }
          Zsave(i,time,qrts) = F0(i,time,qrts)+M(i,time,qrts);

     }
    }
  }
}
// // //
// // //
// //
// //
for(int i=0;i<nage;i++){ // Loop over other ages //
  for(int k=0;k<(nsurvey);k++){
    //
     if(power(i,k) == 1){
       p(i,k) = pin;
     }else{
       p(i,k) = Type(1.0);
     }
    }
  }
//
REPORT(p)
// //
// //
// // //
if(nsurvey>1){
  for(int i=0;i<nage;i++){ // Loop over other ages (recruits excluded)
    for(int k=0;k<(nsurvey);k++){ // Loop over surveys
        if(i >= Qminage(k) && i <= Qmaxage(k)){

            if(i < Qlastage(k)){
              Qsurv(i,k) = Q(Qidx(k)+i-Qminage(k));
            }else{
              Qsurv(i,k) = Q(Qidx(k)+Qlastage(k)-Qminage(k));
            }
        }
        if(i > Qmaxage(k)){
          Qsurv(i,k) = Type(0.0);
        }
      }

    }
}else{
    for(int i=0;i<nage;i++){
      if(i >= Qminage(0) && i <= Qmaxage(0)){
        Qsurv(i,0) = Q(Qidx(0)+i-Qminage(0));
    }
    if(i > Qmaxage(0)){
      Qsurv(i,0) = Type(0.0);
    }
  }
}
// // // //
// // // // //
// // // // // //
if(nsurvey>1){
  for(int k=0;k<(nsurvey);k++){ // Loop over surveys
    for(int i=0;i<nage;i++){ // Loop over other ages
        if(i >= Qminage(k) && i <= Qmaxage(k)){
          SDS(i,k) = SDsurvey(Qidx_CV(i,k));
        }
      }

    }
}else{
  for(int i=0;i<nage;i++){ // Loop over other ages
      if(i >= Qminage(0) && i <= Qmaxage(0)){
          SDS(i,0) = SDsurvey(Qidx_CV(i,0));
        }
      }

    }
// //
// // // // Set up at Nat age in first Quarter
Nsave(0,0,0) = 0; // Nothing here in the fist year

for(int i=1;i<(nage);i++){ // Get initial ages from age 1
  Nsave(i,0,0) = exp(logNinit(i-1));
}
//
for(int time=0;time<(nyears);time++){ // Start time loop
  Rsave(time) = Rin(time);
}
// // // // //
for(int time=0;time<(nyears);time++){ // Start time loop
  for(int qrts=0; qrts<(nseason);qrts++){
      if(qrts == 0){ // Spawning stock biomass is from season 1
        for(int i=0;i<nage;i++){ // Loop over other ages
             SSB(time) += Nsave(i,time,0)*west(i,time,0)*Mat(i,time,0); // Fix SSB
          }
      }
      if(qrts == (recseason-1)){ // Recruitment season
         if(recmodel == 1){ // Different recruitment models
           Rsave(time) = alpha * SSB(time);

             if(SSB(time) >= beta){
                Rsave(time) = alpha*beta;
             }
            Nsave(0,time,2) = Rsave(time);
        }
          if(recmodel == 2){
            Rsave(time) = Rin(time);
            Nsave(0,time,qrts) = Rsave(time);
            logRec(time) = log(Rin(time));
          }
       }
        if(qrts < (nseason-1)){
         for(int i=0;i<(nage);i++){ // Loop over other ages
         Nsave(i,time,qrts+1) = Nsave(i,time,qrts)*exp(-Zsave(i,time,qrts));
        }

        for(int i=0;i<nage;i++){ // Loop over other ages (starts at age 1)

           if(F0(i, time, qrts)>0){
           Catch(i,time,qrts)= (F0(i,time,qrts)/(Zsave(i,time,qrts)))*(1-exp(-Zsave(i,time,qrts)))*Nsave(i,time,qrts)*weca(i,time,qrts);// Calculate the catch in kg
           CatchN(i,time,qrts)= (F0(i,time,qrts)/(Zsave(i,time,qrts)))*(1-exp(-Zsave(i,time,qrts)))*Nsave(i,time,qrts);// Calculate the catch in #s
         }

           for(int k=0;k<nsurvey;k++){

             if(qrts == (surveySeason(k)-1)){
              if(surveyEnd(k) == 0){
               survey(i,time,k) = Nsave(i,time,qrts)*Qsurv(i,k);
             }else{
               Type Ntmp = Qsurv(i,k)*Nsave(i,time,qrts)*(exp(-Zsave(i,time,qrts)*surveyStart(k)));
               survey(i,time,k) = Ntmp*(1-exp(-Zsave(i,time,qrts)*(surveyEnd(k)-surveyStart(k))))/(Zsave(i,time,qrts)*(surveyEnd(k)-surveyStart(k)));
             }
            }
           }
        }

      }else{ // Last quarter

        for(int i=0;i<(nage-1);i++){ // Loop over other ages
        Nsave(i+1,time+1,0) = Nsave(i,time,qrts)*exp(-Zsave(i,time,qrts));
       }
       // Plus group
       Nsave(nage-1,time+1,0) = Nsave(nage-2,time,qrts)*exp(-Zsave(nage-2,time,qrts))+Nsave(nage-1,time,qrts)*exp(-Zsave(nage-1,time,qrts));

       for(int i=0;i<nage;i++){ // Loop over other ages

            Catch(i,time,qrts) = (F0(i,time,qrts)/(Zsave(i,time,qrts)))*(1-exp(-Zsave(i,time,qrts)))*Nsave(i,time,qrts)*weca(i,time,qrts);// Calculate the catch in kg
            CatchN(i,time,qrts)= (F0(i,time,qrts)/(Zsave(i,time,qrts)))*(1-exp(-Zsave(i,time,qrts)))*Nsave(i,time,qrts);// Calculate the catch in #s


          // Calculate survey abundance
          for(int k=0;k<nsurvey;k++){

            if(qrts == (surveySeason(k)-1)){
             if(surveyEnd(k) == 0){
              survey(i,time,k) =  Nsave(i,time,qrts);//*Qsurv(i,k);
            }else{
              Type Ntmp = Nsave(i,time,qrts)*(exp(-Zsave(i,time,qrts)*surveyStart(k)));
              survey(i,time,k) = Ntmp*(1-exp(-Zsave(i,time,qrts)*(surveyEnd(k)-surveyStart(k))))/(Zsave(i,time,qrts)*(surveyEnd(k)-surveyStart(k)));//*Qsurv(i,k)*;
            }
           }
          }
         }
       }
 // //
 } // End Quarter loop
} // End time loop
// // //
// // //
// // // Calculate SSB and recruitment in the new year
// //
for(int i=0;i<nage;i++){ // Loop over other ages
     SSB(nyears) += Nsave(i,nyears,0)*west(i,nyears-1,0)*Mat(i,nyears-1,0); //
     term_logN_next(i) = log(Nsave(i, nyears,0));
}
// // // // //
// // // // // Stock recruitment
vector<Type> SRpred(nyears+1);

for(int time=0;time<(nyears+1);time++){ // Loop over years

 SRpred(time) = alpha + log(SSB(time));
       if(SSB(time) > beta){
          SRpred(time) = alpha+log(beta);
     }

     if(time == nyears){

       if(recseason == 1){
       term_logN_next(0) = SRpred(time);
     }else{
       term_logN_next(0) = -999;
     }
   }
  }
// // // // // // //
Rsave(nyears) = exp(SRpred(nyears));
logRec(nyears) = SRpred(nyears);
// //
// //
// // // Run catch residuals for preliminary calcs
array<Type> resid_catch(nage,nyears, nseason); // Save residuals for SDR calculation
// // For catch variability calculation
// //
for(int time=0;time<(nyears);time++){ // No catches in last year
  for(int i=0;i<nage;i++){ // Loop over other ages
    for(int qrts=0;qrts<nseason;qrts++){ // Loop over seasons

      if(Catchobs(i,time,qrts)> 0 && CatchN(i, time, qrts) > 0){ // Log likelihood
          resid_catch(i,time,qrts) = log(Catchobs(i,time,qrts))-log(CatchN(i,time,qrts));
          }else{
            resid_catch(i, time, qrts) = Type(-99.);
      }
    }
  }
}
//
// // Calculate catch CV internally
int ncatch = catchCV.rows();
int astart = 0;
int aend = 1;
// // //
array<Type> no(ncatch,nseason);
array<Type> sumx(ncatch,nseason);
array<Type> sumx2(ncatch,nseason);

REPORT(ncatch)

for(int k=0;k<(ncatch);k++){ // Loop over number of catch CVs
//
   for(int qrts=0;qrts<nseason;qrts++){ // Loop over other ages

    astart = catchCV(k,qrts);
    if(k == (ncatch-1)){
      aend = nage;
    }else{
      aend = catchCV(k+1, qrts);
    }

     for(int time=0;time<(nyears);time++){ // No catches in last year


       for(int i=astart;i<aend;i++){
          if((Catchobs(i,time,qrts)> 0 && CatchN(i, time, qrts) > 0)){ // Log likelihood
            sumx(k,qrts) += log(CatchN(i,time,qrts))-log(Catchobs(i,time,qrts)); //resid_catch(i,time,qrts)*resid_catch(i,time,qrts);
            sumx2(k,qrts) += pow(log(CatchN(i,time,qrts))-log(Catchobs(i,time,qrts)),2); //resid_catch(i,time,qrts)*resid_catch(i,time,qrts);
            no(k,qrts)+=1;
          }
        }
      }
     }
   }
// //
// // // Now assign SDR to each age
array<Type>SDR_catch2(nage, nseason);
SDR_catch2.setZero();

if(estCV(1) == 2){ // Calculate the catch CV
  for(int k=0;k<(ncatch);k++){ // Loop over number of catch CVs
    for(int qrts=0;qrts<nseason;qrts++){ // Loop over other ages
      astart = catchCV(k,qrts);

      if(k == (ncatch-1)){
        aend = nage;
      }else{
        aend = catchCV(k+1, qrts);
      }

      for(int i=astart;i<aend;i++){
          if(no(k, qrts)>0){ // Only calculate if there are any observations
          SDR_catch2(i,qrts) = sqrt(sumx2(k,qrts)/no(k,qrts));
          //SDR_catch2(i, qrts) = (no(k,qrts)*sumx2(k,qrts)-sumx(k,qrts)*sumx(k,qrts))/(no(k,qrts)*no(k,qrts));
          }
        }
      }
    }
}
//
if(estCV(1) == 0){ // Estimate
  // Fix CV of catches
  for(int k=0;k<(ncatch);k++){ // Loop over number of catch CVs
    for(int qrts=0;qrts<(nseason);qrts++){

      astart = catchCV(k,qrts);

      if(k == (ncatch-1)){
        aend = nage;
      }else{
        aend = catchCV(k+1, qrts);
      }

      for(int i=astart;i<aend;i++){
          SDR_catch2(i,qrts) = SDcatch(k);
          }

      // if(i < catchCV(0,qrts)){
      // SDR_catch2(i,qrts) = Type(0.0);
      // }
      // if( (i >= catchCV(0,qrts)) && (i < catchCV(1,qrts)) ){
      //  SDR_catch2(i,qrts) = SDcatch(i,qrts);
      // }
      // if(i >= (catchCV(1,qrts))){
      // SDR_catch2(i,qrts) = SDcatch(catchCV(1,qrts),qrts);
      // }
    }
  }
}
// //
REPORT(SDR_catch2)
// Add survey residuals
array<Type> resid_survey(nage,nyears,nsurvey); // Save residuals for SDR calculation
Type sumsurv = 0;

for(int time=0;time<nyears;time++){ // Loop over other ages
  for(int i=0;i<nage;i++){ // Loop over other ages
      for(int k=0;k<nsurvey;k++){ // Loop over surveys

        if(Surveyobs(i, time,k) > 0){ // Non existent values have a -1 flag

          resid_survey(i,time,k) = log(Surveyobs(i,time,k))-(log(survey(i, time,k))*p(i,k)+log(Qsurv(i,k)));
          sumsurv += resid_survey(i,time,k);
        }else{
           resid_survey(i,time,k) = Type(-99.);
         }
      }
  }
}


Type nllC = 0.0; // log likelihood for Catch
//
 for(int time=0;time<(nyears);time++){ // No catches in last year
   for(int i=0;i<nage;i++){ // Loop over other ages
     for(int qrts=0;qrts<nseason;qrts++){ // Loop over seasons
       if(Catchobs(i,time,qrts)> 0 && CatchN(i, time, qrts) > 0){ // Log likelihood

       nllC += -dnorm(log(CatchN(i, time, qrts)),log(Catchobs(i, time, qrts)), sqrt(SDR_catch2(i,qrts)), true);
       Catchtot(time) += Catch(i,time, qrts);

     }
    }
   }
 }


Type nllsurv = Type(0.0); // log likelihood for survey observations
array<Type> Surveyout(nage,nyears,nsurvey); // Save residuals for SDR calculation
Surveyout.setZero();

for(int time=0;time<nyears;time++){ // Loop over other ages
  for(int i=0;i<nage;i++){ // Loop over other ages
      for(int k=0;k<nsurvey;k++){ // Loop over surveys


        if(Surveyobs(i, time,k) > 0){ // Non existent values have a -1 flag
        // Export survey numbers
        Surveyout(i,time,k) = exp(log(survey(i, time,k))*p(i,k)+log(Qsurv(i,k)));

        //nllsurv += -dnorm(pow(log(survey(i, time, qrts,k),1)),log(Surveyobs(i, time, qrts,k)), SDS(i,k), true);
        nllsurv += -dnorm(log(survey(i, time,k))*p(i,k)+log(Qsurv(i,k)),log(Surveyobs(i, time,k)), SDS(i,k), true);

        //nllsurv += -dnorm(log(survey(i, time,k)),log(Surveyobs(i, time,k)), SDS(i,k), true);

        }else{
         survey(i,time,k) = Type(-1.0);
       }
      }
    }
 }




// // // //
// // // // // Penalty function for recruitment errors
Type resid_x = 0;
vector<Type> resid_x_export(nyears);
//
// // Calculate SDrec
//
for(int i=0;i<nyears;i++){ // Loop over years
  resid_x += log(Rsave(i))-log(SRpred(i));
  resid_x_export(i) = log(Rsave(i))-log(SRpred(i));
}
//

// //

if(estCV(2) == 2){// Calculate the standard deviation of recruitment
  SDrec = (nyears*pow(resid_x, 2)-pow(resid_x,2))/sqrt(nyears);
}
REPORT(SDrec)
// // //
Type prec = Type(0.0);
vector<Type>ptest(nyears);
//
// // //
for(int time=0;time<(nyears);time++){ // Loop over other years
      prec += -dnorm(SRpred(time),log(Rsave(time)), SDrec, true);
      ptest(time) = -dnorm(SRpred(time),log(Rsave(time)), SDrec, true);
}
// // // // // // //
// // // // // prec += pCV;
// // // // // // //
// // // // // // //
Type ans = 0.0;
//
ans = nllsurv*nllfactor(0)+nllC*nllfactor(1)+prec*nllfactor(2);
// // //
vector<Type> ansvec(3);
ansvec(0) = nllsurv;
ansvec(1) = nllC;
ansvec(2) = prec;
//

REPORT(ansvec)
REPORT(SRpred)
REPORT(ptest)
REPORT(SSB)
REPORT(F0)
REPORT(Catch)
REPORT(CatchN)
REPORT(Rsave)
REPORT(Nsave)
REPORT(Zsave)
REPORT(Fyear)
REPORT(M)
REPORT(Qsurv)
REPORT(term_logN_next)
REPORT(log_exp_pattern)
REPORT(Surveyout)
REPORT(resid_catch)
REPORT(resid_survey)
REPORT(SDS)
REPORT(p)
// REPORT(survey)
// REPORT(ans)
// REPORT(alpha)
// REPORT(Fage)
// REPORT(Fquarter)
// REPORT(Qsurv)
// REPORT(Catchtot)
// REPORT(Surveyobs)
// REPORT(nll)
// //
// //
ADREPORT(SSB)
ADREPORT(SDR_catch2)
ADREPORT(ansvec)
ADREPORT(SDrec)
ADREPORT(F0)
ADREPORT(Catch)
ADREPORT(CatchN)
ADREPORT(Qsurv)
ADREPORT(Surveyout)
ADREPORT(SDS)
ADREPORT(SDR_catch2)
ADREPORT(Rsave)
ADREPORT(resid_catch)
ADREPORT(resid_survey)
ADREPORT(log_exp_pattern)
ADREPORT(term_logN_next)
ADREPORT(logRec)
ADREPORT(Catchtot)
ADREPORT(alpha)
ADREPORT(logbeta)
ADREPORT(SDrec)
ADREPORT(SRpred)
// //
// Type ans = 0.0;


  return ans;
}