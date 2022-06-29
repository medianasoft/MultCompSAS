/********************************************************

Key Multiplicity Issues in Clinical Trials (Part II)
Online training course by Alex Dmitrienko
Downloaded from https://github.com/medianasoft/MultCompSAS

********************************************************/


/********************************************************

Macros

********************************************************/

/*

DESCRIPTION

The MixGate macro implements gatekeeping procedures based on standard 
and modified mixture methods


GENERAL REFERENCE

Brechenmacher, T., Dmitrienko, A. (2017). Multiplicity Adjustment Methods. 
Analysis of Clinical Trials Using SAS (Second Edition). Dmitrienko, A., Koch, G. 
(editors). SAS Press. 

AUTHOR

Thomas Brechenmacher

ARGUMENTS

Indata=[dataset name] specifies the name of the dataset containing the
gatekeeping design of interest, i.e., the individual null hypotheses, 
hierarchical families and logical restrictions.

Method=[Standard, Modif] is used to specify whether the standard or 
modified mixture method is used to construct the gatekeeping 
procedure.

Test=[Bonf, Holm, Hochberg, Hommel] is the component procedure 
used in the gatekeeping procedure.

Gamma=[0<=numerical<1, _NULL_] is the truncation parameter for 
each family (except the last family where the truncation parameter is 
automatically set to 1). Familes are separated by #.

Alpha=[$0<numerical<1] is the global FWER (only used for power 
calculation).

rawp=[dataset name, _NULL_] specifies the name of the dataset 
containing the raw p-values used to compute power (only used for 
power calculation).

adjpout=[dataset name, Adpj] gives the name of the dataset to be 
created with adjusted p-values calculated based on the raw p-values 
contained in dataset study.

powout=[dataset name, Indpow] gives the name of the dataset to be 
created with marginal power for each individual hypothesis (only used for 
power calculation).

pvalout=[dataset name,_NULL_] gives the name of the dataset to be 
created with adjusted p-values calculated based on the raw p-values 
contained in dataset rawp (only used for power calculation).

***************************************************/
    
%macro MixGate(indata=,method=STANDARD,test=Bonf,gamma=_NULL_,alpha=0.025,rawp=_NULL_,
               adjpout=adjp,powout=indpow,pvalout=_NULL_);
    
/****************************** Macro variables ******************************/
%let test2=%upcase(&test);
%let gamma=%upcase(&gamma);
%let pvalout=%upcase(&pvalout);
%let rawp=%upcase(&rawp);
%let method=%upcase(&method);
%if %sysfunc(exist(&indata)) %then %do;    
    data _null_;
        set &indata END=EOF;
            i+1;
            call symput(compress("hyp"||put(i,best12.)),hyp);
            if EOF then do;
            call symput('n_fam',compress(put(family,best12.)));
            call symput('n_hyp',compress(put(i,best12.)));
            end;
    run;
%end;
/**************************************************************************/
    
/***************************** proc IML starts ****************************/
proc IML;

    /*Call dataset Indata to extract family indices*/
    EDIT &indata;
    READ all var {family} into family;

    /*Call dataset Indata to extract logical restriction info*/
    EDIT &indata;
    READ all var {parallel serial} into logical;
        
    /*************** Extract basic information ***************/
    /*Number of hypotheses*/
    n_hyp=nrow(family);
    CALL SYMPUT('n_hyp',left(char(n_hyp)));
    /*Number of families*/
    n_fam=max(family);
    CALL SYMPUT('n_fam',left(char(n_fam)));
    /*Number of null hypotheses in each family*/
    n_Htest=j(1,n_fam,0);
    do i=1 to n_hyp;
        n_Htest[1,family[i,1]]=n_Htest[1,family[i,1]]+1;
    end;
    /*Transpose the Family vector*/
    family=family`;
    /*Number of intersection hypotheses*/
    n_int=2**n_hyp-1;
    /*Serial logical restrictions*/
    parallel=logical[,1];
    serial=logical[,2];
    /*Same logical restrictions data used to check for logical restrictions*/
    serial2=serial;
    parallel2=parallel;
    /********************************************************/

    /**** Raw p-values and number of simulations if any *****/
    %if &rawp eq _NULL_ %then %do;
        %let n_sim=0;
        EDIT &indata;
        READ all var {rawp} into rawp;
        rawp=rawp`;
    %end;
    %else %do;
        EDIT &rawp;
        READ all var _ALL_ into rawp;
        n_sim=nrow(rawp);
        CALL SYMPUT('n_sim',left(char(n_sim)));

        /*Matrix to contain individual power*/    
        indpow=j(1,n_hyp,0);
    %end;
    /********************************************************/
        
    /*********** Generate the 2**n-1 Intersections **********/
    inter=j(n_int,n_hyp,0);
    n=j(n_int,n_fam,0);
    m=j(n_int,n_fam,0);
    inter2=j(n_int,n_hyp,1);
    sercheck=j(n_int,n_hyp,1);
    parcheck=j(n_int,n_hyp,1);
    %if &method=MODIF %then %do;
        nstar=j(n_int,n_fam,0);
        intstar=j(n_int,n_hyp,1);
        sercheck2=j(n_int,n_hyp,1);
    %end;
    do i=1 to n_hyp;
        /*Intersections without logical restrictions*/
        do j=0 to n_int-1;
            k=floor(j/2**(n_hyp-i));
            if k/2=floor(k/2) then inter[j+1,i]=1;
        end;
        sercheck[,i]=inter[,i];
        /** Apply logical restrictions **/
        /*Serial logical restriction*/
        do while (index(serial2[i,1],"1")>0);
            sercheck[,i]=sercheck[,i]#(inter[,index(serial2[i,1],"1")]=0);
            %if &method=MODIF %then %do;
                sercheck2[,i]=sercheck2[,i]#(inter[,index(serial2[i,1],"1")]=0);
            %end;
            serial2[i,1]=concat(tranwrd(substr(serial2[i,1],1,
            index(serial2[i,1],"1")),'1','0'),substr(serial2[i,1],
            index(serial2[i,1],"1")+1));
        end;
        /*parallel logical restriction*/
        if index(parallel2[i,1],"1")>0 then do;
            do while (index(parallel2[i,1],"1")>0);
                parcheck[,i]=parcheck[,i]#(inter[,index(parallel2[i,1],"1")]=1);
                parallel2[i,1]=concat(tranwrd(substr(parallel2[i,1],1,
                index(parallel2[i,1],"1")),'1','0'),substr(parallel2[i,1],
                index(parallel2[i,1],"1")+1));
            end;
        end;
        else parcheck[,i]=0;
        parcheck[,i]=(parcheck[,i]=0);

        inter2[,i]=sercheck[,i]#parcheck[,i];
        %if &method=MODIF %then %do;
            intstar[,i]=sercheck2[,i]#parcheck[,i];
        %end;

        /*Determine n and m*/
        n[,family[,i]]=n[,family[,i]]+inter[,i];
        m[,family[,i]]=m[,family[,i]]+inter2[,i];
        %if &method=MODIF %then %do;
            nstar[,family[,i]]=nstar[,family[,i]]+intstar[,i];
        %end;
    end;  
    free sercheck parcheck serial2 parallel2 %if &method=MODIF %then %do;
    sercheck2 %end;;
    /*Division by 0 not feasible so replace all 0s by missing in m or nstar*/
    %if &method=STANDARD %then %do;
        m=choose(m=0,.,m);
    %end;
    %else %do;
        nstar=choose(nstar=0,.,nstar);
    %end;
    /***************************************************************/

    /*************** Remove duplicate intersections ****************/    
    /*Set-up ntemp to remove duplicate intersections*/
    %if &method=STANDARD %then %do;
        ntemp=n;
        do i=2 to n_fam;
            ntemp[,i]=n[,i]#(ntemp[,i-1]^=n_Htest[,i-1])+n_Htest[,i]
            #(ntemp[,i-1]=n_Htest[,i-1]);
        end;
        inter3=inter2||ntemp||n||m||(1:n_int)`;
    %end;
    %else %do;
        ntemp=nstar;
        do i=2 to n_fam;
            ntemp[,i]=nstar[,i]#(ntemp[,i-1]^=n_Htest[,i-1])+n_Htest[,i]
            #(ntemp[,i-1]=n_Htest[,i-1]);
        end;    
        inter3=inter2||ntemp||nstar||m||(1:n_int)`;
    %end;

    CALL sort(inter3,(1:n_hyp+3*n_fam+1),(1:n_hyp+n_fam));
    interuni=uniqueby(inter3,(1:n_hyp+n_fam-1),1:nrow(inter3));
    inter4=inter3[interuni,1:n_hyp+3*n_fam];
    %if &method=STANDARD %then %do;
        n=inter4[,n_hyp+n_fam+1:n_hyp+2*n_fam];
    %end;
    %else %do;
        nstar=inter4[,n_hyp+n_fam+1:n_hyp+2*n_fam];
    %end;
    m=inter4[,n_hyp+2*n_fam+1:n_hyp+3*n_fam];
    inter2=inter4[,1:n_hyp];
    inter2=choose(inter2=0,.,inter2);
    n_intuni=nrow(inter2);

    free inter3 inter4 ntemp;
    /***************************************************************/

    /***************************** Gamma ***************************/
    /*Define the gamma vector*/
    %if &test2 = BONF %then %do;
        gamma=j(1,n_fam-1,0)||{1};
    %end;
    %else %do;
        /*extract gamma vector*/
        pre_gam=concat(tranwrd("&gamma","#"," ")," 1");
        call symput('gamma2',pre_gam);
        free pre_gam;
        gamma= {&gamma2};
    %end;
    gamnorep=gamma;
    /***************************************************************/

    /* Generate Bonferroni mixing function's critical values: family weights */
    b=j(n_intuni,n_fam,.);
    test=j(n_intuni,1,1);

    %if &method=STANDARD %then %do;
        b[,1]=(n[,1]>0);
            
        do i=2 to n_fam;
            /*Test whether all previous families have no hypothesis*/
            test=test#(n[,i-1]=0);
            /*If family i-1 has no hypothesis then b3=1 otherwise b3=the usual 
            weight*/
            b1=1-n[,i-1];
            b2=(n_Htest[1,i-1]-n[,i-1])#(1-gamma[,i-1])/n_Htest[1,i-1];
            b3=b1<>b2;
            /*Apply the formula to calculate b*/
            b[,i]=test#(n[,i]^=0)<>b[,i-1]#b3;    
        end;
    %end;
    %else %do;
        b[,1]=(m[,1]>0);      
        do i=2 to n_fam;
            /*Test whether all previous families have no hypothesis*/
            test=test#(m[,i-1]=0);
            /*If family i-1 has no hypothesis then b3=1 otherwise b3=the usual 
            weight*/
            b1=1-m[,i-1];
            b2=(nstar[,i-1]-m[,i-1])#(1-gamma[,i-1])/nstar[,i-1];
            b3=b1<>b2;
            /*Apply the formula to calculate b*/
            b[,i]=test#(m[,i]>0)<>b[,i-1]#b3;    
        end;
    %end;
    free test b1 b2 b3;
    /*Division by 0 not feasible so replace all 0s by missing in b;*/
    b=choose(b=0,.,b);
    /*Division by 0 not feasible so replace all 0s by missing in m*/
    m=choose(m=0,.,m);
   /***************************************************************/

    /******* Prepare variables for the gatekeeping procedure ******/
    /*Repeat variables for all hypotheses*/
    do fam=1 to n_fam;
        if fam=1 then do;
            gamma2=repeat(gamma[,fam],1,n_Htest[,fam]);
            mrep=repeat(m[,fam],1,n_Htest[,fam]);
            %if &method=STANDARD %then %do;
                n_Htest2=repeat(n_Htest[,fam],1,n_Htest[,fam]);
            %end;
            %else %do;
                nstarep=repeat(nstar[,fam],1,n_Htest[,fam]);
            %end;
        end;
        else do;
            gamma2=gamma2||repeat(gamma[,fam],1,n_Htest[,fam]);
            mrep=mrep||repeat(m[,fam],1,n_Htest[,fam]);
            %if &method=STANDARD %then %do;
                n_Htest2=n_Htest2||repeat(n_Htest[,fam],1,n_Htest[,fam]);
            %end;
            %else %do;
                nstarep=nstarep||repeat(nstar[,fam],1,n_Htest[,fam]);
            %end;
        end;
    end;
    gamma=gamma2;
    m=mrep;
    %if &method=STANDARD %then %do;
        n_Htest=n_Htest2;
    %end;
    %else %do;
        nstar=nstarep;
    %end;
    /*Transpose raw p-values*/
    t_rawp=rawp`;
   /***************************************************************/
        
    /******************* Gatekeeping procedure ********************/
    %if &n_sim ne 0 %then %do; do simulation=1 to &n_sim;%end;
    %else %do;do simulation=1 to 1;%end;
         
        /*p-values ordering for each family*/
        temp=family`||t_rawp[,simulation]||(1:n_hyp)`;
        call sort(temp,{1 2});
        ord_p=temp[,3]`;
        free temp;

        interraw=rawp[simulation,]#inter2;
        interord=j(n_intuni,n_hyp);
        interadj=j(n_intuni,n_hyp);
        interp=j(n_intuni,n_hyp);
        interm=j(1,n_hyp);
        min_p=j(n_intuni,n_fam);
        final_p=j(n_intuni,1);
        temp=j(n_intuni,n_fam,.);
        min_p_b=j(n_intuni,n_fam);
        partemp=j(1,n_hyp,1);
            
        /*raw-pvalue order*/
        interord[,]=interraw[,ord_p];
        free interraw;
        prek=(interord>0);

        /*define k to be used in the truncated test formula*/
        do i=1 to n_hyp;    
            if i^=1 then do;
                if family[1,i] = family[1,i-1] then prek[,i]=prek[,i]+prek[,i-1];
            end;
        end;

        /*replace all zeros by 1*/
        k=prek+(prek=0);
        free prek;
     
        /*truncated test except for last family*/
        %if &method=STANDARD %then %do;
            %if &test2=HOMMEL %then %str(interadj=interord/
            (k#gamma/m + (1-gamma)/n_Htest););
            %else %if &test2=HOCHBERG %then %str(interadj=interord/
            (gamma/(m-k+1)+(1-gamma)/n_Htest););
            %else %if (&test2=HOLM or &test2=BONF) %then %str(interadj=interord/
            (gamma/m+(1-gamma)/n_Htest););
        %end;
        %else %do;
            %if &test2=HOMMEL %then %str(interadj=interord/
            (k#gamma/m+(1-gamma)/nstar););
            %else %if &test2=HOCHBERG %then %str(interadj=interord/
            (gamma/(m-k+1)+(1-gamma)/nstar););
            %else %if &test2=HOLM %then %str(interadj=interord/
            (gamma/m+(1-gamma)/nstar););
        %end;
        /*replace all missing by 1*/
        interadj=choose(interadj=.,1,interadj);
        /*Take the minimum of the calculated p-values in each family*/
        do i=1 to n_hyp;
            if i=1 then min_p[,1]=interadj[,1];
            else do;
                if family[1,i] = family[1,i-1] then min_p[,family[1,i]]=
                min_p[,family[1,i]]><interadj[,i];
                else min_p[,family[1,i]]=interadj[,i];
            end;
        end;
        /*divide by bi*/
        min_p_b=min_p/b;
            
        /*Replace missing values by 1*/
        temp2=(min_p_b=temp);
        min_p_b=min_p_b<>temp2;

        /*Compute intersection adjusted p-values*/
        do fam=1 to n_fam;
            final_p[,1]=final_p[,1]><min_p_b[,fam];
        end;
        free interord k interadj min_p min_p_b temp2 temp;    
            
        /*Compute the max over all intersections*/
        interp=final_p[,1]#inter2;
        free final_p;
        do intermax=1 to n_hyp;
            interm[,intermax]=max(interp[1:n_intuni,intermax]);        
        end;

        free interp;
                
        /*Add calculation for the final step (non-consonance for Hommel test)*/
        serial3=serial;
        parallel3=parallel;

        do i=1 to n_hyp;
            /*Serial logical restriction*/
            do while (index(serial3[i,1],"1")>0);
                interm[,i]=interm[,i]<>interm[,index(serial3[i,1],"1")];
                serial3[i,1]=concat(tranwrd(substr(serial3[i,1],1,
                index(serial3[i,1],"1")),'1','0'),substr(serial3[i,1],
                index(serial3[i,1],"1")+1));
            end;
            /*parallel logical restriction*/
            if index(parallel3[i,1],"1")>0 then do;
                do while (index(parallel3[i,1],"1")>0);
                    partemp[,i]=partemp[,i]><interm[,index(parallel3[i,1],"1")];
                    parallel3[i,1]=concat(tranwrd(substr(parallel3[i,1],1,
                    index(parallel3[i,1],"1")),'1','0'),substr(parallel3[i,1],
                    index(parallel3[i,1],"1")+1));
                end;
                interm[,i]=interm[,i]<>partemp[,i];
            end;
        end;
        free partemp serial3 parallel3;
        /********************************************************/
            
        /******************* calculate power ********************/
        %if &n_sim ne 0 %then %do;
            pow=(interm<=&alpha);
            /*Power for individual hypotheses*/
            indpow=indpow+pow;
        %end;
        /********************************************************/

        /**************** Output final datasets *****************/
        %else %do;
            allpval=gamnorep||interm;
            create adjp from allpval;
            append from allpval;
        %end;

        %if &pvalout ne _NULL_ and &n_sim ne 0 %then %do;
            if simulation=1 then &pvalout=interm;
            else &pvalout=&pvalout//interm;
        %end;
        free interm;
    end;

    %if &n_sim ne 0 %then %do;
        indpow=gamnorep||(indpow*100/&n_sim);
        create indpow from indpow;
        append from indpow;

        %if &pvalout ne _NULL_ %then %do;
            index=(1:&n_sim)`;
            call sort(index,1);
            &pvalout=index||repeat(gamnorep,&n_sim)||&pvalout;
            create &pvalout from &pvalout;
            append from &pvalout;
            free index;
        %end;
    %end;

quit;

/************************ proc IML ends ****************************/
/** Results when simulation **/
%if &n_sim ne 0 %then %do;
    /** Power for individual hypotheses **/
    data &powout (keep=gamma1-gamma%eval(&n_fam-1) hyp1-hyp%eval(&n_hyp));
        set indpow;    
            array gamma {%eval(&n_fam-1)};
            array hyp {&n_hyp};
            array col{%eval(&n_fam+&n_hyp)};
            %do i=1 %to &n_hyp;        
                hyp{&i}=col{&i+&n_fam};    
            %end;
            %do i=1 %to &n_fam-1; 
                gamma&i=col&i;
            %end;;
            label %do i=1 %to &n_hyp; 
                hyp&i="Individual power: &&hyp&i"
            %end;;    
    run;
        
%end;
/** Results when no simulation **/
%else %do;
    data &adjpout (keep=gamma1-gamma%eval(&n_fam-1) adj_p1-adj_p&n_hyp);
        set adjp;    
            array gamma {%eval(&n_fam-1)};
            array adj_p {&n_hyp};
            array col{%eval(&n_fam+&n_hyp)};
            do i=1 to &n_hyp;        
                adj_p{i}=col{i+&n_fam};    
            end;
            format adj_p1-adj_p&n_hyp PVALUE6.4;

            %do i=1 %to &n_fam-1; 
                gamma&i=col&i;
            %end;;
            label %do i=1 %to &n_hyp; 
            adj_p&i="Adjusted p-value: &&hyp&i"
            %end;;    
    run;
%end;
/** Output adjusted p-values for all simulations **/
%if &pvalout ne _NULL_ %then %do;
    data &pvalout (keep=simnum gamma1-gamma%eval(&n_fam-1) adj_p1-adj_p&n_hyp);
        length simnum 8;
        set &pvalout;                  
            array gamma {%eval(&n_fam-1)};
            array adj_p {&n_hyp};
            array col{%eval(&n_fam+&n_hyp+1)};
            simnum=col{1};
            do i=1 to &n_hyp;        
                adj_p{i}=col{1+i+&n_fam};    
            end;
            format adj_p1-adj_p&n_hyp PVALUE6.4;

            %do i=1 %to &n_fam-1; 
                gamma&i=col%eval(1+&i);
            %end;;
            label simnum="Simulation number" %do i=1 %to &n_hyp; 
                adj_p&i="Adjusted p-value: &&hyp&i"
            %end;;    
    run;
%end;
%mend MixGate;

/********************************************************

Module C

Computation of adjusted p-values for the Hochberg-based 
parallel gatekeeping procedure

********************************************************/

* Example 2: Type 2 diabetes mellitus trial;

data ex2;
    input hyp $ family parallel $ serial $ rawp;
    datalines;
    H1 1 000 000 0.0082 
    H2 1 000 000 0.0174 
    H3 2 110 000 0.0202 
run;

%MixGate(indata=ex2, method=Standard, 
         test=Hochberg, gamma=0.7, 
         adjpout=adjp);

proc print data=adjp noobs label;
    var adj_p1-adj_p3;
run;

/********************************************************

Module D

Computation of adjusted p-values for the Hochberg-based 
general gatekeeping procedure

********************************************************/

* Example 3: Schizophrenia trial;

data ex3;
    input hyp $ family parallel $ serial $ rawp;
    datalines;
    H1 1 0000 0000 0.0101 
    H2 1 0000 0000 0.0233 
    H3 2 0000 1000 0.0022 
    H4 2 0000 0100 0.0167 
run;

%MixGate(indata=ex3, method=Standard, 
         test=Hochberg, gamma=0.7, 
         adjpout=adjp);

proc print data=adjp noobs label;
    var adj_p1-adj_p4;
run;
