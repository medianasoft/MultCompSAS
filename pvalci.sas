/*

DESCRIPTION

The PvalCI macro computes lower limits of one-sided simultaneous confidence 
intervals for the following p-value-based multiple testing procedures:

Weighted Bonferroni procedure
Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L., Tamhane, A.C., 
Hsu, J.C. (2009). Multiple testing methodology. Multiple Testing Problems in 
Pharmaceutical Statistics. Dmitrienko, A., Tamhane, A.C., Bretz, F. 
(editors). Chapman and Hall/CRC Press, New York. 

Weighted Holm procedure
Holm, S. (1979). A simple sequentially rejective multiple test procedure. 
Scandinavian Journal of Statistics. 6, 65-70.

Fixed-sequence procedure
Westfall, P. H., Krishen, A. (2001). Optimally weighted, fixed sequence, 
and gatekeeping multiple testing procedures. Journal of Statistical Planning 
and Inference. 99, 25-40.

Weighted fallback procedure
Wiens, B. (2003). A fixed-sequence Bonferroni procedure for testing multiple 
endpoints. Pharmaceutical Statistics. 2, 211-215.
Wiens, B., Dmitrienko, A. (2005). The fallback procedure for evaluating a 
single family of hypotheses. Journal of Biopharmaceutical Statistics. 15, 929-942.

GENERAL REFERENCE

Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L., Tamhane, A.C., 
Hsu, J.C. (2009). Multiple testing methodology. Multiple Testing Problems in 
Pharmaceutical Statistics. Dmitrienko, A., Tamhane, A.C., Bretz, F. 
(editors). Chapman and Hall/CRC Press, New York. 

AUTHOR

Alex Dmitrienko

SOURCE

This macro and related macros can be downloaded from https://github.com/medianasoft/MultCompSAS 

ARGUMENTS
 
IN = Name of the data set with raw p-values (RAW_P variable),
     treatment differences (EST variable), standard errors (SE) 
     and weights (WEIGHT variable). 
     The hypotheses are assumed to be ordered from first to last.

COVPROB = Simultaneous coverage probability.

OUT = Name of the data set with lower simultaneous confidence limits.
*/

%macro pvalci(in,covprob,out);
proc iml;

    * Compute the weighted Bonferroni p-value for an intersection hypothesis;
    start bonf(p,u);
        bonfp=min(p[loc(u)]/u[loc(u)]);
        return(bonfp);
    finish bonf;

    * Read in raw p-values, treatment differences, standard errors and weights;
    use &in;        
    read all var {raw_p} into p; 
    read all var {est} into est; 
    read all var {se} into se; 
    read all var {weight} into weight;  
    w=t(weight);

    * One-sided familywise error rate;
    alpha=1-&covprob;

    * Step 1: Compute adjusted p-values;

        nhyps=nrow(p);
        nints=2**nhyps-1;
        h=j(nints,nhyps,0);

        * Variables for the Bonferroni procedure;
        bonferroniv=j(nints,nhyps,0);
        bonferronip=j(nints,nhyps,0); 

        * Variables for the Holm procedure;
        holmv=j(nints,nhyps,0);
        holmp=j(nints,nhyps,0); 

        * Variables for the fixed-sequence procedure;
        fixedseqp=j(nints,nhyps,0); 

        * Variables for the fallback procedure;
        fallbackv=j(nints,nhyps,0);
        fallbackp=j(nints,nhyps,0); 

        * Construct a closed family;
        do i=1 to nhyps;
            do j=0 to nints-1;
                k=floor(j/2**(nhyps-i));
                if k/2=floor(k/2) then h[j+1,i]=1;
            end;
        end;        

        * Compute local p-values for intersection hypotheses;
        do i=1 to nints;

            * Compute hypothesis weights for the fixed-sequence procedure;
            minindex=0;
            do j=1 to nhyps;
                if h[i,j]=1 & minindex=0 then minindex=j;
            end;

            * Compute hypothesis weights for the fallback procedure;
            do j=1 to nhyps;
                if (h[i,j]=1) then do;
                    if (j=1) then sum=w[1];
                    if (j>1) then do;
                        continue=1;
                        sum=w[j];
                        do k=j-1 to 1 by -1;
                            if (continue=1 && h[i,k]=0) then sum=sum+w[k];
                            if (h[i,k]=1) then continue=0;
                        end;
                    end;
                    fallbackv[i,j]=sum;      
                end;
            end;

            * Bonferroni procedure;
            bonferroniv[i,]=w#h[i,];
            bonferronip[i,]=h[i,]*bonf(p,bonferroniv[i,]);
            * Holm procedure;
            holmv[i,]=w#h[i,]/sum(w#h[i,]);
            holmp[i,]=h[i,]*bonf(p,holmv[i,]);
            * Fixed-sequence procedure;
            fixedseqp[i,]=h[i,]*p[minindex];
            * Fallback procedure;
            fallbackp[i,]=h[i,]*bonf(p,fallbackv[i,]);
        end;

        * Vectors of adjusted p-values;
        bonf=j(nhyps,1,0);    
        holm=j(nhyps,1,0);    
        fixedseq=j(nhyps,1,0);    
        fallback=j(nhyps,1,0);    

        * Compute adjusted p-values;    
        do i=1 to nhyps; 
            bonf[i]=max(bonferronip[,i]); 
            holm[i]=max(holmp[,i]); 
            fixedseq[i]=max(fixedseqp[,i]); 
            fallback[i]=max(fallbackp[,i]); 
        end;  

    * Step 2: Compute lower limits of one-sided simultaneous confidence intervals;   

        * Vectors of confidence limits;
        bonfci=j(nhyps,1,0);    
        holmci=j(nhyps,1,0);    
        fixedseqci=j(nhyps,1,0);    
        fallbackci=j(nhyps,1,0);    

        * Rejection/acceptance of null hypotheses;
        holmsi=(holm<=alpha);
        fixedseqsi=(fixedseq<=alpha);
        fallbacksi=(fallback<=alpha);

        * Bonferroni procedure;
        do i=1 to nhyps;
            bonfci[i]=est[i]-se[i]*probit(1-alpha*weight[i]);
        end;

        * Holm procedure;
        if sum(holmsi)=nhyps then do;
            do i=1 to nhyps;
                holmci[i]=0 <> est[i]-se[i]*probit(1-alpha*weight[i]);
            end;
        end;
        if sum(holmsi)<nhyps then do;
            do i=1 to nhyps; 
                if holmsi[i]=1 then holmci[i]=0;
                if holmsi[i]=0 then do;
                    adjalpha=alpha*w[i]/sum(w[loc(holmsi=0)]);               
                    holmci[i]=est[i]-se[i]*probit(1-adjalpha);
                end;
            end; 
        end;    

        * Fixed-sequence procedure;
        if sum(fixedseqsi)=0 then do;
            fixedseqci[1]=est[1]-se[1]*probit(1-alpha);             
            do i=2 to nhyps; 
                fixedseqci[i]=.;
            end;
        end;
        if sum(fixedseqsi)=nhyps then do;
            temp=est-se*probit(1-alpha);    
            fixedseqci=min(temp);
        end;
        if 0<sum(fixedseqsi) & sum(fixedseqsi)<nhyps then do;
            fixedseqci[1]=0;
            do i=2 to nhyps; 
                if fixedseqsi[i]=1 then fixedseqci[i]=0;
                if fixedseqsi[i]=0 & fixedseqsi[i-1]=1 then fixedseqci[i]=est[i]-se[i]*probit(1-alpha);
                if fixedseqsi[i]=0 & fixedseqsi[i-1]=0 then fixedseqci[i]=.;            
            end; 
        end;    

        * Fallback procedure;
        if sum(fallbacksi)=nhyps then do;
            do i=1 to nhyps;
                fallbackci[i]=0 <> (est[i]-se[i]*probit(1-alpha*w[i]));
            end;
        end;
        if sum(fallbacksi)<nhyps then do;
            accept=1-t(fallbacksi);
            * Find the intersection hypothesis corresponding to the set of accepted hypotheses;
            do i=1 to nints; 
                if h[i,]=accept then accrow=i;
            end;         
            do i=1 to nhyps; 
                if fallbacksi[i]=1 then do;
                    first=1;
                    * Identify all subsets of the set of accepted hypotheses;
                    do j=1 to nints;
                        if h[j,]<=accept then do;                        
                            alphaj=alpha*sum(fallbackv[j,]);
                            adjalpha=(alpha-alphaj)*w[i]/sum(w[loc(h[j,]=0)]);
                            if adjalpha>0 then templimit=0 <> (est[i]-se[i]*probit(1-adjalpha));
                            if adjalpha=0 then templimit=0;
                            if first=1 then do;
                                minlimit=templimit;
                                first=0;
                            end;
                            minlimit=minlimit >< templimit;                        
                        end;
                    end;
                    fallbackci[i]=minlimit;
                end;
                if fallbacksi[i]=0 then do;           
                    adjalpha=alpha*fallbackv[accrow,i]; 
                    fallbackci[i]=est[i]-se[i]*probit(1-adjalpha);
                end;
            end; 
        end;    

    * Organize the results;
    ci=j(nhyps,8,0);
    ci[,1]=est;
    ci[,2]=se;
    ci[,3]=weight;
    ci[,4]=est-se#probit(1-alpha);
    ci[,5]=bonfci;
    ci[,6]=holmci;
    ci[,7]=fixedseqci;
    ci[,8]=fallbackci;

    * Create output data set;
    create ci from ci[colname={est se weight univariate bonferroni holm fixedseq fallback}];
    append from ci;
    quit;

data &out;
    set ci;
    test=_n_;
    label test="Test"
          est="Estimate"
          se="Standard error"
          weight="Weight"
          univariate="Univariate"
          bonferroni="Bonferroni"
          holm="Holm"
          fixedseq="Fixed-sequence"
          fallback="Fallback";
    run;
%mend pvalci;
