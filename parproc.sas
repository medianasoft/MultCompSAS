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

SOURCE

This macro and related macros can be downloaded from https://github.com/medianasoft/MultCompSAS 

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
