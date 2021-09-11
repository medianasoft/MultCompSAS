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

AUTHOR

Alex Dmitrienko

SOURCE

This macro and related macros can be downloaded from https://github.com/medianasoft/MultCompSAS 

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

