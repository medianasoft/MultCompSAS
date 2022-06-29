/********************************************************

Key Multiplicity Issues in Clinical Trials (Part I)
Online training course by Alex Dmitrienko
Downloaded from https://github.com/medianasoft/MultCompSAS

********************************************************/


/********************************************************

Macros

********************************************************/

/*

DESCRIPTION

The PvalProc macro computes adjusted p-values for the following p-valued-based 
multiple testing procedures:

Weighted Bonferroni procedure
Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L., Tamhane, A.C., 
Hsu, J.C. (2009). Multiple testing methodology. Multiple Testing Problems in 
Pharmaceutical Statistics. Dmitrienko, A., Tamhane, A.C., Bretz, F. 
(editors). Chapman and Hall/CRC Press, New York. 

Weighted Holm procedure
Holm, S. (1979). A simple sequentially rejective multiple test procedure. 
Scandinavian Journal of Statistics. 6, 65-70.

Weighted Hommel procedure
Hommel, G. (1988). A stagewise rejective multiple test procedure based on a 
modified Bonferroni test. Biometrika. 75, 383-386.

Weighted Hochberg procedure
Tamhane, A.C., Liu, L. (2008). On weighted Hochberg procedures. Biometrika. 
95, 279-294.

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

ARGUMENTS
 
IN =  Name of the data set with weights (WEIGHT variable) and
      raw p-values (RAW_P variable). The hypotheses
      are assumed to be ordered from first to last.

OUT = Name of the data set with adjusted p-values.

*/

%macro pvalproc(in,out);
proc iml;    

    * Compute the weighted Bonferroni p-value for an intersection hypothesis;
    start bonf(p,u);
        bonfp=min(p[loc(u)]/u[loc(u)]);
        return(bonfp);
    finish bonf;

    * Compute the weighted Simes p-value for an intersection hypothesis;
    start simes(p,w);
        pp=t(p[loc(w)]);
        ww=t(w[loc(w)]);
        comb=pp//ww; 
        temp=comb;   
        comb[,rank(pp)]=temp;    
        simesp=min(comb[1,]/cusum(comb[2,]));
        return(simesp);
    finish simes;  

    * Compute the weighted incomplete Simes p-value for an intersection hypothesis;
    start incsimes(p,w);
        pp=t(p[loc(w)]);
        ww=t(w[loc(w)]);
        k=ncol(pp);
        if k>1 then do;
            comb=pp//ww; 
            temp=comb;   
            comb[,rank(pp)]=temp;  
            modw=ww;
            modw[1]=0;
            modw[2:k]=comb[2,1:k-1];
            incsimesp=min((1-cusum(modw))#comb[1,]/comb[2,]);
        end;
        if k=1 then incsimesp=pp;
        return(incsimesp);
    finish incsimes; 

    * Read in weights and raw p-values;
    use &in;
    read all var {weight raw_p} into data;        
    data=t(data);
    w=data[1,];
    p=data[2,];    
    nhyps=ncol(p); 
    nints=2**nhyps-1;
    h=j(nints,nhyps,0);

    * Variables for the Bonferroni procedure;
    bonferroniv=j(nints,nhyps,0);
    bonferronip=j(nints,nhyps,0); 

    * Variables for the Holm procedure;
    holmv=j(nints,nhyps,0);
    holmp=j(nints,nhyps,0); 

    * Variables for the Hommel procedure;
    hommelv=j(nints,nhyps,0); 
    hommelp=j(nints,nhyps,0); 

    * Variables for the Hochberg procedure;
    hochbergv=j(nints,nhyps,0); 
    hochbergp=j(nints,nhyps,0); 

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
        * Hommel procedure;
        hommelv[i,]=w#h[i,]/sum(w#h[i,]);
        hommelp[i,]=h[i,]*simes(p,hommelv[i,]);
        * Hochberg procedure;
        hochbergv[i,]=w#h[i,]/sum(w#h[i,]);
        hochbergp[i,]=h[i,]*incsimes(p,hochbergv[i,]);
        * Fixed-sequence procedure;
        fixedseqp[i,]=h[i,]*p[minindex];
        * Fallback procedure;
        fallbackp[i,]=h[i,]*bonf(p,fallbackv[i,]);
    end;

    * Compute adjusted p-values;    
    adjp=j(nhyps,7,0);
    do i=1 to nhyps; 
        adjp[i,1]=p[i]; 
        adjp[i,2]=max(bonferronip[,i]); 
        adjp[i,3]=max(holmp[,i]); 
        adjp[i,4]=max(hommelp[,i]); 
        adjp[i,5]=max(hochbergp[,i]); 
        adjp[i,6]=max(fixedseqp[,i]); 
        adjp[i,7]=max(fallbackp[,i]); 
    end;        

    * Create output data set;
    create adjp from adjp[colname={raw bonferroni holm hommel hochberg fixedseq fallback}];
    append from adjp;        
    quit;

data &out;
    set adjp;
    test=_n_;
    label test="Test"
          raw="Raw"
          bonferroni="Bonferroni"
          holm="Holm"
          hommel="Hommel"
          hochberg="Hochberg"
          fixedseq="Fixed-sequence"
          fallback="Fallback";
    run;
%mend pvalproc;

/*

DESCRIPTION

The ChainSer macro computes adjusted p-values for the serial chain procedure defined in

Millen, B., Dmitrienko, A. (2010). Chain procedures: A class of flexible closed 
testing procedures with clinical trial applications. Statistics in Biopharmaceutical Research. 

GENERAL REFERENCE

Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L., Tamhane, A.C., 
Hsu, J.C. (2009). Multiple testing methodology. Multiple Testing Problems in 
Pharmaceutical Statistics. Dmitrienko, A., Tamhane, A.C., Bretz, F. 
(editors). Chapman and Hall/CRC Press, New York. 

AUTHOR

Alex Dmitrienko

ARGUMENTS
 
IN  = Name of the data set with raw p-values and procedure parameters:
         RAW_P: raw p-value
         WEIGHT: hypothesis weight
         G1-Gn: transition parameters (n is the number of null hypotheses)
         The null hypotheses are assumed to be ordered from first to last.

OUT = Name of the data set with adjusted p-values.

*/

%macro chainser(in,out);

proc iml;    

    * Compute the weighted Bonferroni p-value for an intersection hypothesis;
    start bonf(p,u);
        bonfp=min(p[loc(u)]/u[loc(u)]);
        return(bonfp);
    finish bonf;

    * Read in weights and raw p-values;
    use &in;
    read all var _NUM_ into data;   
    nhyps=nrow(data); 
    * Raw p-values;
    p=t(data[,1]);
    * Hypothesis weights;
    w=t(data[,2]);    
    * Transition parameters;
    g=data[,3:ncol(data)];
 
    nints=2**nhyps-1;
    h=j(nints,nhyps,0);
    s=j(nints,nhyps,0);
    v=j(nints,nhyps,0);
    localp=j(nints,nhyps,0);

    * Construct a closed family;
    do i=1 to nhyps;
        do j=0 to nints-1;
            k=floor(j/2**(nhyps-i));
            if k/2=floor(k/2) then h[j+1,i]=1;
        end;
    end;            

    * Compute local p-values for intersection hypotheses;
    do i=1 to nints;
        s[i,1]=w[1];
        v[i,1]=h[i,1]*s[i,1];
        do j=2 to nhyps;
            s[i,j]=w[j];
            do k=1 to (j-1);
                s[i,j]=s[i,j]+(1-h[i,k])*g[k,j]*s[i,k];
            end;
            v[i,j]=h[i,j]*s[i,j];
        end;
        * Local p-value;
        localp[i,]=h[i,]*bonf(p,v[i,]);
    end;

    * Compute adjusted p-values;    
    adjp=j(nhyps,2,0);
    do i=1 to nhyps; 
        adjp[i,1]=p[i]; 
        adjp[i,2]=max(localp[,i]); 
    end;      

    * Create output data set;
    create adjp from adjp[colname={raw chain}];
    append from adjp;        
    quit;

* Output data set;
data &out;
    set adjp;
    test=_n_;
    label test="Test"
          raw="Raw"
          chain="Chain";
    run;

%mend chainser;

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

/*

DESCRIPTION

The ParProc macro computes adjusted p-values for the following parametric 
multiple testing procedures in one-sided testing problems with a balanced 
one-way layout and equally weighted null hypotheses:

Single-step Dunnett procedure
Dunnett, C.W. (1955). A multiple comparison procedure for comparing several 
treatments with a control. Journal of  the American Statistical Association. 
50, 1096-1121. 

Step-down Dunnett procedure
Naik, U.D. (1975). Some selection rules for comparing p processes with a standard. 
Communications in Statistics. Series A. 4, 519-535.
Marcus, R., Peritz, E., Gabriel, K.R. (1976). On closed testing procedures with 
special reference to ordered analysis of variance. Biometrika. 63, 655-660. 
Dunnett, C.W., Tamhane, A.C. (1991). Step-down multiple tests for comparing 
treatments with a control in unbalanced one-way layouts. Statistics in Medicine. 
10, 939-947.


GENERAL REFERENCE

Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L., Tamhane, A.C., 
Hsu, J.C. (2009). Multiple testing methodology. Multiple Testing Problems in 
Pharmaceutical Statistics. Dmitrienko, A., Tamhane, A.C., Bretz, F. 
(editors). Chapman and Hall/CRC Press, New York. 

AUTHOR

Alex Dmitrienko

ARGUMENTS
 
IN =  Name of the data set with test statistics (T variable). 

N =   Number of patients per treatment arm (assumed to be constant across 
      the treatment arms).

OUT = Name of the data set with adjusted p-values.

*/

%macro parproc(in,n,out);
proc iml;    

    * Read in test statistics;
    use &in;
    read all var {t} into stat; 
    stat=t(stat); 
    m=ncol(stat); 
    n=&n;
    nu=(m+1)*(n-1);

    * Adjusted p-values;    
    adjp=j(m,3,0);

    * Univariate p-values;
    do i=1 to m;
        adjp[i,1]=1-probt(stat[i],2*(n-1));
    end;

    * Single-step Dunnett procedure;
    do i=1 to m;
        adjp[i,2]=1-probmc("DUNNETT1",stat[i],.,nu,m);
    end;

    * Step-down Dunnett procedure;
    sorttest=j(m,1,0);    
    temp=j(m,1,0);
    r=m+1-rank(stat);
    sorttest[r]=stat;    
    do i=1 to m;
        if i=1 then do;
            temp[1]=1-probmc("DUNNETT1",sorttest[1],.,nu,m);
            max=temp[1];
        end;
        if i>1 then do;
            temp[i]=max <> (1-probmc("DUNNETT1",sorttest[i],.,nu,m-i+1));
            max=max <> temp[i];
        end;
    end;
    adjp[,3]=temp[r];

    * Create output data set;
    create adjp from adjp[colname={raw dunnett stepdunnett}];
    append from adjp;        
    quit;

data &out;
    set adjp;
    test=_n_;
    label test="Test"
          raw="Raw"
          dunnett="Single-step Dunnett"
          stepdunnett="Step-down Dunnett";
    run;

%mend parproc;

/*

DESCRIPTION

The ParCI macro computes lower limits of one-sided simultaneous confidence 
intervals for the following parametric multiple testing procedures in 
one-sided testing problems with a balanced one-way layout and equally 
weighted null hypotheses:

Single-step Dunnett procedure
Step-down Dunnett procedure
Bofinger, E. (1987). Step-down procedures for comparison with a control. 
Australian Journal of Statistics. 29, 348-364.
Stefansson, G., Kim, W.-C., Hsu, J.C. (1988). On confidence sets in multiple 
comparisons. Statistical Decision Theory and Related Topics IV. Gupta, S.S., 
Berger, J.O. (editors). Academic Press, New York, 89-104.

GENERAL REFERENCE

Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L., Tamhane, A.C., 
Hsu, J.C. (2009). Multiple testing methodology. Multiple Testing Problems in 
Pharmaceutical Statistics. Dmitrienko, A., Tamhane, A.C., Bretz, F. 
(editors). Chapman and Hall/CRC Press, New York. 

AUTHOR

Alex Dmitrienko

ARGUMENTS
 
IN =  Name of the data set with test statistics (T variable), means (EST
      variable) and common standard error (SE variable). 

N =   Number of patients per treatment arm (assumed to be constant across 
      the treatment arms).

COVPROB = Simultaneous coverage probability.

OUT = Name of the data set with lower limits of one-sided simultaneous 
      confidence intervals.

*/

%macro parci(in,n,covprob,out);
proc iml;    

    * Read in test statistics, means and common standard deviation;
    use &in;
    read all var {t} into stat; 
    read all var {est} into est; 
    read all var {se} into temp; 
    stat=t(stat); 
    se=temp[1];
    m=ncol(stat); 
    n=&n;
    nu=(m+1)*(n-1);

    * Lower limits of one-sided simultaneous confidence intervals;
    ci=j(m,5,0);
    ci[,1]=est;
    ci[,2]=se;

    est=t(est); 

    * Lower limits of one-sided univariate confidence intervals;
    do i=1 to m;
        ci[i,3]=est[i]-tinv(&covprob,2*(n-1))*se;
    end;

    * Single-step Dunnett procedure;
    c=probmc("DUNNETT1",.,&covprob,nu,m);
    do i=1 to m;
        ci[i,4]=est[i]-c*se;
    end;

    * Step-down Dunnett procedure;

    * Step 1: Compute adjusted p-values;
    
        sorttest=j(m,1,0);    
        temp=j(m,1,0);
        r=m+1-rank(stat);
        sorttest[r]=stat;    
        do i=1 to m;
            if i=1 then do;
                temp[1]=1-probmc("DUNNETT1",sorttest[1],.,nu,m);
                max=temp[1];
            end;
            if i>1 then do;
                temp[i]=max <> (1-probmc("DUNNETT1",sorttest[i],.,nu,m-i+1));
                max=max <> temp[i];
            end;
        end;
        adjp=temp[r];
        reject=(adjp<=1-&covprob);

    * Step 2: Compute lower limits of one-sided simultaneous confidence intervals;
    
        nreject=sum(reject);
        if (nreject=m) then do;
            c=tinv(&covprob,nu);
            do i=1 to m;
                ci[i,5]=0 <> (est[i]-c*se); 
            end;
        end;
        if (nreject<m) then do;
            do i=1 to m; 
                if reject[i]=1 then ci[i,5]=0;
                if reject[i]=0 then do;
                    c=probmc("DUNNETT1",.,&covprob,nu,m-nreject);               
                    ci[i,5]=est[i]-c*se;
                end;
            end; 
        end;    

    * Create output data set;
    create ci from ci[colname={est se univariate dunnett stepdunnett}];
    append from ci;        
    quit;

data &out;
    set ci;
    test=_n_;
    label test="Test"
          est="Mean"
          se="Standard error"
          univariate="Univariate"
          dunnett="Single-step Dunnett"
          stepdunnett="Step-down Dunnett";
    run;

%mend parci;


options nodate nocenter ls=250;

/********************************************************

Module C

Computation of adjusted p-values for stepwise procedures 
with a data-driven testing sequence

********************************************************/

* Example 4: Type 2 diabetes trial;

data ex4;
   input raw_p;
   weight=1/3;
   datalines;
   0.0111
   0.0065
   0.0293
run;

proc multtest pdata=ex4 holm;
    title1 "Example 4: Type 2 diabetes trial";
    title2 "MULTTEST procedure";
    run;

%pvalproc(in=ex4,out=adjp);
proc print data=adjp noobs label;
    title1 "Example 4: Type 2 diabetes trial";
    title2 "PVALPROC macro";
    format raw holm 6.4;
    var test raw holm;
    run;

* Example 1: Prostate cancer trial;

data ex1; 
    input raw_p weight;
    datalines; 
    0.0102 0.8 
    0.0181 0.2 
run;

%pvalproc(in=ex1,out=adjp);
proc print data=adjp noobs label;
    title1 "Example 1: Prostate cancer trial";
    title2 "PVALPROC macro";
    format raw holm 6.4;
    var test raw holm;
    run;

/********************************************************

Module D

Computation of adjusted p-values for stepwise procedures 
with a pre-specified testing sequence

********************************************************/

* Example 4: Type 2 diabetes trial;

data ex4;
    input raw_p weight;
    datalines;
    0.0061 0.3333
    0.0233 0.3333
    0.0098 0.3334
run;

%pvalproc(in=ex4,out=adjp);
proc print data=adjp noobs label;
    title1 "Example 4: Type 2 diabetes trial";
    title2 "PVALPROC macro";
    format raw fallback 6.4;
    var test raw fallback;
    run;

data ex4chain;
    input raw_p weight g1-g3;
    datalines;
    0.0061 0.3333 0.0 0.5 0.5
    0.0233 0.3333 0.0 0.0 1.0
    0.0098 0.3334 0.0 0.0 0.0
run;

%chainser(in=ex4chain,out=adjp);
proc print data=adjp noobs label;
    title1 "Example 4: Type 2 diabetes trial";
    title2 "CHAINSER macro";
    format raw chain 6.4;
    var test raw chain;
    run;

/********************************************************

Module E

Computation of adjusted p-values for parametric procedures

********************************************************/

* Example 4: Type 2 diabetes trial;

data ex4;
    input t;
    datalines;
    2.64
    1.93
    2.31
run;

%parproc(in=ex4,n=90,out=adjp);
proc print data=adjp noobs label;
    title1 "Example 4: Type 2 diabetes trial";
    title2 "PARPROC macro";
    format raw dunnett 6.4;
    var test raw dunnett;
    run;

/********************************************************

Module F

Computation of lower limits of one-sided simultaneous confidence 
intervals for nonparametric and parametric procedures

********************************************************/

* Example 4: Type 2 diabetes trial;

data ex4;
   input raw_p est sd;
   weight=1/3;
   se=sd*sqrt(2/90);
   datalines;
   0.0045 0.63 1.6
   0.0276 0.46 1.6
   0.0110 0.55 1.6
run;

%pvalci(in=ex4,covprob=0.975,out=adjci);
proc print data=adjci noobs label;
   title1 "Example 4: Type 2 diabetes trial";
   title2 "PVALCI macro";
   format univariate holm 6.3;
   var test univariate holm;
   run;

* Example 4: Type 2 diabetes trial;

data ex4;
   input t est sd;
   se=sd*sqrt(2/90);
   datalines;
   2.64 0.63 1.6
   1.93 0.46 1.6
   2.31 0.55 1.6
run;

%parci(in=ex4,n=90,covprob=0.975,out=adjci);
proc print data=adjci noobs label;
   title1 "Example 4: Type 2 diabetes trial";
   title2 "PARCI macro";
   format univariate dunnett 6.3;
   var test univariate dunnett;
   run;
