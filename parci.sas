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

SOURCE

This macro and related macros can be downloaded from https://github.com/medianasoft/MultCompSAS 

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
