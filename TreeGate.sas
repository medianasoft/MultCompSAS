/*=============================================================================================================
---------------------------------------------------------------------------------------------------------------
Program Name:         TreeGate_v1.0.SAS
Main macro:			  TreeGate
sub macros:			  Powfunc, Report
Purpose:              1/ Compute adjusted p-values using the tree-structured gatekeeping procedure (Hommel, 
						 Hochberg or Holm-based).
					  2/ Compute power for prespecified power functions given a dataset of raw p-values.
					  3/ Compute familywise error rate given a dataset of raw p-values.
					  Please refer to the user manual for a detailed description of the macro.
Requirements:		  Dataset study must be created before running the macro.

%macro TreeGate(test=HOMMEL,exhaust=NO,gamma=_ALL_,gamby=0.1,gaminf=0,gamsup=0.9,powout=POWER,rawp=_NULL_,pvalout=_NULL_);

Macro options:        test    [Default=HOMMEL]: 	Multiplicity comparison procedure to be used with the gatekeeping
							   						procedure. A truncation version of the test is used in each family
													except the last one where the regular test is used. Available 
													values are Hommel, Hochberg and Holm.
					  gamma   [Default=_ALL_] 		Truncation fraction (gamma) value for each family. Values are separated  
							   						by #. Ex: gamma=0.5#0.9 for a gatekeeping with 3 familes where
													gamma1=0.5 and gamma2=0.9 (gamma3 is automatically set to 1).
													For gamma=_ALL_, results are calculated for gamma values from GAMINF
													to GAMSUP and incremented by GAMBY.
					  gaminf  [Default=0] 			Minimum value for gamma when GAMMA=_ALL_.
													Ex: for GAMINF=0.5, results are calculated for values of gamma from 0.5 to 0.9
					  gamsup  [Default=0.9] 		Maximum value for gamma when GAMMA=_ALL_.
													Ex: for GAMSUP=0.5, results are calculated for values of gamma from 0 to 0.5
					  gamby   [Default=0.1] 		Increment value for gamma when GAMMA=_ALL_.
													Ex: for GAMBY=0.1, results are calculated for values of gamma 0.0, 0.1, 
													0.2,E0.9
					  exhaust [Default=No]			If EXHAUST=Yes, the alpha-exhaustive gatekeeping procedure is performed.
					  alpha	  [Default=0.025]		Global familywise error rate
					  powout  [Default=POWER]		Ouput type for simulations. POWOUT=power requests power calculation and 
													POWOUT=fwer requests global familywise error rate calculation.
					  rawp   [Default=_NULL_]		Specify the dataset containing the raw p-values used for
							   						simulation. Only applicable for simulations. 
													When RAWP=_null_, no simulation is performed.
					  pvalout [Default=_NULL_]		Specify the dataset to contain all adjusted p-values
							   						calculated based on the raw p-values contained in RAWP.
													Only applicable for simulations. When PVALOUT=_null_, 
													the adjusted p-values are not saved. Note that the 
													running time of the macro is increased when adjusted
													p-values are requested.


Author:      	      Thomas Brechenmacher
OS/SAS Vnum:          Windows XP / 9.1
Date:                 09AUG2011

---------------------------------------------------------------------------------------------------------------
Note: 				  Macro SIMUL is given as an example and generate raw p-values though simulation for the
					  case of m families with n null hypotheses to be tested in each family. Macro options for 
					  macro simul are:
				      n_Htest			     		Number of null hypotheses to be tested in each family.
                      n_fam				     		Number of families.
                      n_subj				  		Number of subjects in each group (only balanced designs).
                      eff_size				  		Effect sizes. Effect sizes for each family are separated by # 
													and ordered by group within a family.
					  sigma					  		Correlation matrix.
					  n_sim					  		Number of simulations to be performed.
					  seed     [Default=5963]		Seed to be used for the random vector generating the test
							    					statistics. If seed=0 then the system time clock is used.
					  out							Name of the SAS dataset to contain the resulting raw p-values
---------------------------------------------------------------------------------------------------------------

---------------------------------------------------------------------------------------------------------------
 Modifications:
 Author/Date:         Reason:
 ----------------     -------

=============================================================================================================*/

/******************************* Simulation: create matrix of raw p-values *******************************/
%macro simul(n_Htest=,n_fam=,n_subj=,eff_size=,sigma=,n_sim=,seed=5963,out=);
	proc IML;
	/*Shape the correlation matrix and effect-size vector*/
	sigma=tranwrd("&sigma","#",",");
	mu=tranwrd("&eff_size","#"," ");
	call symput('sigma2',sigma);
	call symput('mu',mu);
	free sigma mu;

	/*Define the effect-size vector*/
	mu= {&mu}`;

	/*Define the covariance matrix of the group mean vector*/
	sigma= (1/&n_subj)*block(%do b=1 %to &n_Htest;{&sigma2},%end;{&sigma2}); 

	/*Define the contrast matrix*/
	contrast=j(&n_Htest*&n_fam,(&n_Htest+1)*&n_fam,0);
	do block=1 to &n_fam;
		do ligne=1 to &n_Htest;
			contrast[&n_Htest*(block-1)+ligne,block+ligne*&n_fam]=1;
			contrast[&n_Htest*(block-1)+ligne,block]=-1;
		end;
	end;	

	/*var(Z)=(sqrt(n_sim/2))*contrast*sigma*transpose(contrast)*/
	varZ=(&n_subj/2)*contrast*sigma*contrast`;

	/*generate the test statistic vector. test statistic=(Y2-Y1)/sqrt(2var/n_sim)*/
	CALL vnormal(Z,(sqrt(&n_subj/2))*mu,varZ,&n_sim %if &seed ne 0 %then ,&seed;);	
    /*p-values*/
	rawp=1-probnorm(Z);

	create &out from rawp;
		append from rawp;

	quit;
%mend;
/*****************************************************************************************************/

/********************************* Extract power definitions if any **********************************/
%macro powfunc();
	/*Organize data*/
	data def;
		set POWfunc;
			 defnraw+1;
	run;
	proc sort data=def;
		by powlabel;
	run;
	data def(drop=lastdef);
		length lastdef $40;
		set def end=EOF;
			lastdef=lag(powlabel);	
			if powlabel ne lastdef then defnum+1;
			if eof then call symput('n_def',compress(put(defnum,best12.)));
	run;
	proc sort data=def nodupkey;
		by powlabel defnum;
	run;
	/*Create macro variables for definition labels*/
	proc sort data=def;
		by defnraw;
	run;
	data def;
		set def(drop=defnum defnraw);
			defnum+1;
			call symput(compress("def"||put(defnum,best12.)),trim(left(powlabel)));	
	run;
	proc sort data=def;
		by powlabel;
	run;

	/*Prepare data for proc IML and define definition labels for proc report*/
	proc sort data=POWfunc out=pdefIML;
		by powlabel;
	run;
	data pdefIML(drop=powlabel defnum F1-F&n_fam H1-H&n_hyp i footnote label1-label&n_hyp)
		footnote (keep=powlabel defIML footnote);
		length defIML weight CF1-CF&n_fam FNUM1-FNUM&n_fam CH1-CH&n_hyp HNUM1-HNUM&n_hyp 8 footnote $1000 label1-label&n_hyp $50;
		merge pdefIML(in=a) def;
		by powlabel;
			if a;
			defIML=defnum;			
			array F{&n_fam};
			array CF{&n_fam};
			array Fnum{&n_fam};
			array H{&n_hyp};
			array CH{&n_hyp};
			array Hnum{&n_hyp};
			array label{&n_hyp};
			%do i=1 %to &n_hyp;
				label&i="&&hyp&i";
			%end;
			if weight=1 then footnote="Reject ";
			else footnote=compress(put(weight,best12.))||"*P(reject ";

			do i=1 to dim(F);
				if index(F{i},"-")>0 then do;
					Fnum{i}=input(tranwrd(F{i},"-",""),best12.);
					CF{i}=2;
					
					if i ne dim(F) then 
					footnote=compbl(footnote)||"exactly "||compress(put(Fnum{i},best12.))||" hyp. in F"||compress(put(i,best12.)||",");
					else footnote=compbl(footnote)||"exactly "||compress(put(Fnum{i},best12.))||" hyp. in F"||compress(put(i,best12.));
				end;
				else do;
					Fnum{i}=input(F{i},best12.);
					CF{i}=1;
					if Fnum{i} ne 0 then do;
						if i ne dim(F) then footnote=compbl(footnote)||"at least "||compress(put(Fnum{i},best12.))||" hyp. in F"||compress(put(i,best12.)||",");
						else footnote=compbl(footnote)||"at least "||compress(put(Fnum{i},best12.))||" hyp. in F"||compress(put(i,best12.));
					end;
				end;
			end;
			first=0;
			do i=1 to dim(H);
				if index(H{i},"-")>0 then do;
					Hnum{i}=input(tranwrd(H{i},"-",""),best12.);
					CH{i}=2;
					if first=0 and index(footnote,"hyp.")>0 then footnote=trim(left(trim(right(compbl(footnote)))))||",";
					if i ne dim(H) then footnote=compbl(footnote)||"not "||trim(label{i})||",";
					else footnote=compbl(footnote)||"not "||trim(left(right(label{i})));
					first=1;
				end;
				else do;
					Hnum{i}=input(H{i},best12.);
					if Hnum{i} eq 1 then do;
						 CH{i}=2;
						if first=0 and index(footnote,"hyp.")>0 then footnote=trim(left(trim(right(compbl(footnote)))))||",";
						if i ne dim(H) then footnote=compbl(footnote)||trim(label{i})||",";
						else footnote=compbl(footnote)||trim(label{i});
						first=1;
					end;
					else do;
						CH{i}=1;
					end;
				end;
			end;
			if substr(footnote,length(footnote),1)="," then footnote=substr(footnote,1,length(footnote)-1);
			if weight ne 1 then footnote=trim(footnote)||")";
	run;
	
	proc sort data=footnote;
		by defIML;
	run;
	data footnote(keep=powlabel footnote);
		length lastfoot $1000;
		set footnote;
		by defIML;
			i+1;
			lastdef=lag(defIML);	
			retain lastfoot "";
			if first.defIML then test="";
			if defIML ne lastdef then do;
				footnote=footnote;
				lastfoot=footnote;
			end;
			else do;
				footnote=trim(lastfoot)||"#+"||trim(footnote);
				lastfoot=footnote;
			end;
			if last.defIML;
	run;
%mend;
/*****************************************************************************************************/

/********************************* Report gatekeeping results ****************************************/
%macro report();
	OPTIONS center
	nonumber
    nodate
    nomprint
	nofmterr
    nomerror
    nomlogic
	notes
	source
    formdlim='' 
	mergenoby=warn
	formchar='=_'
	;
	options	ORIENTATION=LANDSCAPE 	LEFTMARGIN=1 in	RIGHTMARGIN=1 in
						 			TOPMARGIN=1.25 in	BOTTOMMARGIN=1 in
			PAPERSIZE=A4 FONT="Courier New" 8 ;
	options ls=145 ps=41;

	%if &exhaust=YES %then %str(title1 "&test-based gatekeeping procedure (alpha-exhaustive)";);
	%else %str(title1 "&test-based gatekeeping procedure";);
	/** Results when simulation **/
	%if &n_sim ne 0 %then %do;
		/*Apply format x.y;*/
		data _null_;
			/*Decimal number for power*/
			numlen=1/&n_sim;
			charpval=put(numlen,best12.);
			charperc=put(numlen*100,best12.);
			/*Calculate the overall length of the float: x*/
			if index(charpval,"E")=0 then lengpval=put(length(compress(charpval)),best12.);
			else lengpval=put(input(scan(compress(charpval),2,"E-"),best12.)+length(scan(compress(charpval),1,"E-"))+1,best12.);
			if index(charperc,"E")=0 then lengper1=put(length(compress(charperc)),best12.);
			else lengpval=put(input(scan(compress(charperc),2,"E-"),best12.)+length(scan(compress(charperc),1,"E-"))+1,best12.);
			/*Calculate the decimal numbers: y*/
			if index(charpval,".")>0 then decipval=put(length(scan(charpval,2,".")),best12.);
			else if index(charpval,"E")>0 then decipval=scan(compress(charpval),2,"E-");
			else decipval="0";	
			if index(charperc,".")>0 then deciperc=put(length(scan(charperc,2,".")),best12.);
			else if index(charperc,"E")>0 then deciperc=scan(compress(charperc),2,"E-"); 
			else deciperc="0";
			/*Special cases when n<=10*/
			if compress(lengper1) in ("2","3") and compress(deciperc)="0" then lengper1="1";
			/*the integer part in a percentage can either be 1, 2 or 3 characters long*/
			lengper2=compress(put(input(lengper1,best12.)+1,best12.));
			lengper3=compress(put(input(lengper1,best12.)+2,best12.));
			/*Compute the local macro variables with the overall length and decimal numbers*/
			formpval=compress(lengpval||"."||decipval);
			if compress(lengper1)="12" then do; /*Results are reported under the format best12, length=12 is a 
												  special case for percentage*/
				formper1=compress(lengper1||"."||"10");
				formper2=compress(lengper1||"."||"9");
				formper3=compress(lengper1||"."||"8");
			end;
			else if compress(lengper1)="11" then do; /*Results are reported under the format best12, length=11 
													   is a special case for percentage*/
				formper1=compress(lengper1||"."||"9");
				formper2=compress(lengper1||"."||"9");
				formper3=compress(lengper1||"."||"8");
			end;
			else do;
				formper1=compress(lengper1||"."||deciperc); /*general expression of the macro variables*/
				formper2=compress(lengper2||"."||deciperc);
				formper3=compress(lengper3||"."||deciperc);
			end;
			formpval=compress(lengpval||"."||decipval);

			call symput('formpval',formpval);
			call symput('formper1',formper1);
			call symput('formper2',formper2);
			call symput('formper3',formper3);
		run;
		%if &powout=POWER %then %do;
			/** Power for power functions **/
			title2 "Power results for %sysfunc(compress(&n_sim)) simulations";
			%if %sysfunc(exist(POWfunc)) %then %do;
				data power (keep=gamma1-gamma%eval(&n_fam-1) power1-power%eval(&n_def))
					 powpres (keep=gamma1-gamma%eval(&n_fam-1) powerc1-powerc%eval(&n_def));
				set powerIML;	
					array gamma {%eval(&n_fam-1)};
					array power {&n_def};
					array powerc {&n_def} $20;
					array col{%eval(&n_fam+&n_def)};
						%do i=1 %to &n_def;		
							power{&i}=col{&i+&n_fam};	
							if length(scan(compress(put(col{&i+&n_fam},best12.)),1,"."))=1 then do;
								powerc{&i}=put(col{&i+&n_fam},&formper1);
								format power&i &formper1;
							end;
							else if length(scan(compress(put(col{&i+&n_fam},best12.)),1,"."))=2 then do;
								powerc{&i}=put(col{&i+&n_fam},&formper2);
								format power&i &formper2;
							end;
							else if length(scan(compress(put(col{&i+&n_fam},best12.)),1,"."))=3 then do; 
								powerc{&i}=put(col{&i+&n_fam},&formper3);
								format power&i &formper3;
							end;

						%end;
					%do i=1 %to &n_fam-1; 
						gamma&i=col&i;
				    %end;;
					label %do i=1 %to &n_def; 
						power&i="Power: &&def&i"
				    %end;;	
				run;
				/*Mark highest power in each family*/
				%if &gamma=_ALL_ %then %do;
					proc sql;
					   create table powpres2 as select *,
					   %do i=1 %to &n_def-1;max(powerc&i) as max&i,%end;
					   max(powerc&n_def) as max&n_def
					   from powpres;
					quit;

					data powpres2;
						set powpres2;
							%do i=1 %to &n_def;
								if powerc&i=max&i then powerc&i=compress(powerc&i||"*");
							%end;
					run;
				%end;
				%else %do;
					data powpres2;
						set powpres;
					run;
				%end;
				/*Report power for each value of gamma*/
				title3 "Definition of power functions";
				proc report data=footnote headline headskip nowindows spacing=3 split='#' missing ps=39;
					column ('____#  ' powlabel footnote);

					
						define powlabel	/	order order=data	"Power function"	  width=20 left flow;
						define footnote / 	"Definition"         width=90 left flow;

						break after powlabel / skip;

				run;	

				title3 "Power results: Power functions";
				proc report data=powpres2 headline headskip nowindows spacing=3 split='#' missing ps=39;
					column ('____#  ' gamma1-gamma%eval(&n_fam-1) 
									  ("Power (%)#__" powerc1-powerc&n_def));

					%do i=1 %to %eval(&n_fam-2);
						define gamma&i	/	order order=data	"Gamma &i"	  width=7 left;
					%end;
					define gamma%eval(&n_fam-1)	/	id order order=data	"Gamma %eval(&n_fam-1)"	  width=7 left;
					%do i=1 %to &n_def;
						define powerc&i	/						"&&def&i" 	  width=12 left;
					%end;

					break after gamma1 / skip;
					%if &gamma=_ALL_ %then %do;
						compute after _page_;
							line @10 '*: Highest power for each definition';;
						endcomp;
					%end;

				run;	
			%end;
		
			/** Power for individual hypotheses **/
			data indpower (keep=gamma1-gamma%eval(&n_fam-1) hyp1-hyp%eval(&n_hyp))
					 indpres (keep=gamma1-gamma%eval(&n_fam-1) hypc1-hypc%eval(&n_hyp));
				set indIML;	
					array gamma {%eval(&n_fam-1)};
					array hyp {&n_hyp};
					array hypc {&n_hyp} $20;
					array col{%eval(&n_fam+&n_hyp)};
						%do i=1 %to &n_hyp;		
							hyp{&i}=col{&i+&n_fam};	
							if length(scan(compress(put(col{&i+&n_fam},best12.)),1,"."))=1 then do;
								hypc{&i}=put(col{&i+&n_fam},&formper1);
								format hyp&i &formper1;
							end;
							else if length(scan(compress(put(col{&i+&n_fam},best12.)),1,"."))=2 then do;
								hypc{&i}=put(col{&i+&n_fam},&formper2);
								format hyp&i &formper2;
							end;
							else if length(scan(compress(put(col{&i+&n_fam},best12.)),1,"."))=3 then do; 
								hypc{&i}=put(col{&i+&n_fam},&formper3);
								format hyp&i &formper3;
							end;

						%end;
					%do i=1 %to &n_fam-1; 
						gamma&i=col&i;
				    %end;;
					label %do i=1 %to &n_hyp; 
						hyp&i="Individual power: &&hyp&i"
				    %end;;	
			run;
			%if &gamma=_ALL_ %then %do;
				proc sql;
				   create table indpres2 as select *,
				   %do i=1 %to &n_hyp-1;max(hypc&i) as max&i,%end;
				   max(hypc&n_hyp) as max&n_hyp
				   from indpres;
				quit;

				data indpres2;
					set indpres2;
						%do i=1 %to &n_hyp;
							if hypc&i=max&i then hypc&i=compress(hypc&i||"*");
						%end;
				run;
			%end;
			%else %do;
				data indpres2;
					set indpres;
				run;
			%end;
			/*Report power for each value of gamma*/				
			title3 "Power results: Individual null hypotheses";
			proc report data=indpres2 headline headskip nowindows spacing=3 split='#' missing ps=39;
				column ('____#  ' gamma1-gamma%eval(&n_fam-1) 
								  ("Power (%)#__" hypc1-hypc&n_hyp));

				%do i=1 %to %eval(&n_fam-2);
					define gamma&i	/	order order=data	"Gamma &i"	  width=7 left;
				%end;
				define gamma%eval(&n_fam-1)	/	id order order=data	"Gamma %eval(&n_fam-1)"	  width=7 left;
				%do i=1 %to &n_hyp;
					define hypc&i	/						"&&hyp&i" 	  width=10 left;
				%end;

				break after gamma1 / skip;
				%if &gamma=_ALL_ %then %do;
					compute after _page_;
						line @10 '*: Highest power for each null hypothesis';;
					endcomp;
				%end;
			run;
		%end;
		/** FWER **/
		%if &powout=FWER %then %do;
			data FWER (keep=gamma1-gamma%eval(&n_fam-1) FWER)
				 FWERpres (keep=gamma1-gamma%eval(&n_fam-1) FWERc);
				length gamma1-gamma%eval(&n_fam-1) FWER 8 FWERc $20;
				set FWERIML;	
					FWER=col%eval(&n_fam+1);
					FWERc=compress(put(fwer,&formpval));
					format fwer &formpval;

					%do i=1 %to &n_fam-1; 
						gamma&i=col&i;
				    %end;;
					label FWER="FWER";
			run;
			title2 "Familywise error rate (FWER) for %sysfunc(compress(&n_sim)) simulations";
			title3;
			proc report data=FWERpres headline headskip nowindows spacing=3 split='#' missing ps=39;
				column ('____#  ' gamma1-gamma%eval(&n_fam-1) FWERc);

				%do i=1 %to %eval(&n_fam-1);
					define gamma&i	/	order order=data	"Gamma &i"	  width=7 left;
				%end; 
				define FWERc	/							"FWER"  		width=12 left;	

				break after gamma1 / skip;
			run;
		%end;
	%end;
	/** Results when no simulation **/
	%else %do;
		data adjp (keep=gamma1-gamma%eval(&n_fam-1) adj_p1-adj_p&n_hyp)
				 adj_pres (keep=gamma1-gamma%eval(&n_fam-1) adj_pc1-adj_pc&n_hyp);
			set adj_pIML;	
				array gamma {%eval(&n_fam-1)};
				array adj_p {&n_hyp};
				array adj_pc {&n_hyp} $20;
				array col{%eval(&n_fam+&n_hyp)};
					do i=1 to &n_hyp;		
						adj_p{i}=col{i+&n_fam};	
						adj_pc{i}=put(col{i+&n_fam},PVALUE6.4);
					end;
				format adj_p1-adj_p&n_hyp PVALUE6.4;

				%do i=1 %to &n_fam-1; 
					gamma&i=col&i;
			    %end;;
				label %do i=1 %to &n_hyp; 
					adj_p&i="Adjusted p-value: &&hyp&i"
			    %end;;	
		run;
		title2 "Adjusted p-values";
		title3;
		proc report data=adj_pres headline headskip nowindows spacing=3 split='#' missing ps=39;
			column ('____#  ' gamma1-gamma%eval(&n_fam-1) 
							  ("Adjusted p-value corresponding to#__" adj_pc1-adj_pc&n_hyp));

			%do i=1 %to %eval(&n_fam-2);
				define gamma&i	/	order order=data	"Gamma &i"	  width=7 left;
			%end;
			define gamma%eval(&n_fam-1)	/	id order order=data	"Gamma %eval(&n_fam-1)"	  width=7 left;
			%do i=1 %to &n_hyp;
				define adj_pc&i	/						"&&hyp&i" 	  width=10 left;
			%end;

			break after gamma1 / skip;
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
%mend;
/*****************************************************************************************************/

/******************************** Implement TreeGate procedure ************************************/

%macro TreeGate(test=Hommel,gamma=_ALL_,gaminf=0,gamsup=0.9,gamby=0.1,exhaust=NO,alpha=0.025,powout=POWER,rawp=_NULL_,pvalout=_NULL_);
	
	/********************************************* Macro variables *****************************************/
	%let test2=%upcase(&test);
	%let gamma=%upcase(&gamma);
	%let powout=%upcase(&powout);
	%let pvalout=%upcase(&pvalout);
	%let rawp=%upcase(&rawp);
	%let exhaust=%upcase(&exhaust);
	%if %sysfunc(exist(study)) %then %do;	
		data _null_;
				set study END=EOF;
				i+1;
				call symput(compress("hyp"||put(i,best12.)),hyp);
				if EOF then do;
					call symput('n_fam',compress(put(family,best12.)));
					call symput('n_hyp',compress(put(i,best12.)));
				end;
		run;
	%end;
	/*******************************************************************************************************/

	/******************************* Check existence of required information *******************************/	
	%if %sysfunc(exist(study)) %then %do;
		%let dsid = %sysfunc(open(study));
		%let chk1 = %sysfunc(varnum(&dsid, hyp));
		%let chk2 = %sysfunc(varnum(&dsid, family));
		%let chk3 = %sysfunc(varnum(&dsid, parallel));
		%let chk4 = %sysfunc(varnum(&dsid, serial));
		%let chk5 = %sysfunc(varnum(&dsid, rawp));
		%let rc = %sysfunc(close(&dsid));
	%end;
	%if %sysfunc(exist(POWfunc)) and &powout=POWER and &rawp ne _NULL_ %then %do;
		%let dsid = %sysfunc(open(POWfunc));
		%let chk6 = %sysfunc(varnum(&dsid, powlabel));
		%let chk7 = %sysfunc(varnum(&dsid, weight));
		%let rc = %sysfunc(close(&dsid));
	%end;
	%else %do;
		%let chk6=1;
		%let chk7=1;
	%end;

	%if not %sysfunc(exist(study)) %then %PUT ERROR: Dataset STUDY needs to be defined. Please refer to the user manual.;
	%else %if &powout=FWER and not %sysfunc(exist(FWERdef)) %then %PUT ERROR: Dataset FWERDEF needs to be defined. Please refer to the user manual.;
	%else %if &test2 ne HOMMEL and &test2 ne HOCHBERG and &test2 ne HOLM %then %PUT ERROR: Macro option TEST can only take the values HOMMEL, HOCHBERG or HOLM.;
	%else %if &rawp ne _NULL_ and &powout ne POWER and &powout ne FWER %then %PUT ERROR: Please set macro option POWOUT to either POWER or FWER.; 
	%else %if &exhaust ne YES and &exhaust ne NO %then %PUT ERROR: Macro option EXHAUST can only take the values YES or NO.;
	%else %if &chk1=0 or &chk2=0 or &chk3=0 or &chk4=0 or (&chk5=0 and &rawp=_NULL_) %then %do; 
		%if &chk1=0 %then %PUT ERROR: Variable HYP needs to be created in dataset STUDY;
		%if &chk2=0 %then %PUT ERROR: Variable FAMILY needs to be created in dataset STUDY;
		%if &chk3=0 %then %PUT ERROR: Variable PARALLEL needs to be created in dataset STUDY;
		%if &chk4=0 %then %PUT ERROR: Variable SERIAL needs to be created in dataset STUDY;
		%if &chk5=0 and &rawp=_NULL_ %then %PUT ERROR: Variable rawp needs to be created in dataset STUDY;
	%end;
	%else %if &chk6=0 or &chk7=0 %then %do;
		%if &chk6=0 %then %PUT ERROR: Variable POWLABEL needs to be created in dataset POWfunc;
		%if &chk7=0 %then %PUT ERROR: Variable WEIGHT needs to be created in dataset POWfunc;
	%end;
	%else %if &powout=FWER and not %sysfunc(exist(FWERdef)) %then %PUT ERROR: Dataset FWERDEF needs to be defined. Please refer to the user manual.;

	/*******************************************************************************************************/
	%else %do;
		/*********************************** Extract power definitions if any **********************************/
		/*Power functions*/
		%if %sysfunc(exist(POWfunc)) and &powout=POWER and &rawp ne _NULL_ %then %do;
			%powfunc;
		%end;
		
		/*******************************************************************************************************/
		
		/******************************************* proc IML starts *******************************************/
		proc IML;

		/*Call dataset study to extract family indices*/
		EDIT study;
		READ all var {family} into family;

		/*Call dataset study to extract logical restriction info*/
		EDIT study;
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
		/*Number of values for gamma*/
		%if &gamma eq _ALL_ %then %do;
			n_gam=floor((%sysevalf((&gamsup-&gaminf+&gamby)/&gamby)))**(n_fam-1);	
		%end;
		%else n_gam=1;;
		/********************************************************/

		/**** Raw p-values and number of simulations if any *****/
		%if &rawp eq _NULL_ %then %do;
			%let n_sim=0;
			EDIT study;
			READ all var {rawp} into rawp;
			rawp=rawp`;
		%end;
		%else %do;
			EDIT &rawp;
			READ all var _ALL_ into rawp;
			n_sim=nrow(rawp);
			CALL SYMPUT('n_sim',left(char(n_sim)));

			/*Matrix to contain individual power and FWER*/	
			indpow=j(n_gam,n_hyp,0);
			FWER=j(n_gam,1,0);
		%end;
		/********************************************************/

		/*************** Power definitions if any ***************/
		%if %sysfunc(exist(POWfunc)) and &powout=POWER and &rawp ne _NULL_ %then %do;		
			EDIT pdefIML;
				READ all var _ALL_ into POWfunc;
			powerDEF=j(n_gam,nrow(POWfunc),0);
			/*Number of power functions*/
			alldef=POWfunc[,{1 2}];
			n_def=POWfunc[nrow(POWfunc),1];

			/*Matrix to contain power*/
			power=j(n_gam,n_def,0);
		%end;
		/********************************************************/	

		/*********** FWER definition if appropriate *************/
		%if %sysfunc(exist(FWERdef)) and &powout=FWER and &rawp ne _NULL_ %then %do;
			EDIT FWERdef;
					READ all var _ALL_ into FWERdef;
		%end;
		/********************************************************/
		
	    /*********** Generate the 2**n-1 Intersections **********/
	    inter=j(n_int,n_hyp,0);
		n=j(n_int,n_fam,0);
		m=j(n_int,n_fam,0);
		inter2=j(n_int,n_hyp,1);
		sercheck=j(n_int,n_hyp,1);
		parcheck=j(n_int,n_hyp,1);
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
				serial2[i,1]=concat(tranwrd(substr(serial2[i,1],1,index(serial2[i,1],"1")),'1','0'),substr(serial2[i,1],index(serial2[i,1],"1")+1));
			end;
			/*parallel logical restriction*/
			if index(parallel2[i,1],"1")>0 then do;
				do while (index(parallel2[i,1],"1")>0);
					parcheck[,i]=parcheck[,i]#(inter[,index(parallel2[i,1],"1")]=1);
					parallel2[i,1]=concat(tranwrd(substr(parallel2[i,1],1,index(parallel2[i,1],"1")),'1','0'),substr(parallel2[i,1],index(parallel2[i,1],"1")+1));
				end;
			end;
			else parcheck[,i]=0;
			parcheck[,i]=(parcheck[,i]=0);

			inter2[,i]=sercheck[,i]#parcheck[,i];

			/*Determine n and m*/
			n[,family[,i]]=n[,family[,i]]+inter[,i];
			m[,family[,i]]=m[,family[,i]]+inter2[,i];
	    end;  
		free sercheck parcheck serial2 parallel2;
		/*Division by 0 not feasible so replace all 0s by missing in m*/
		m=choose(m=0,.,m);
		/********************************************************/

		/*************** Remove duplicate intersections ****************/	
		/*Set-up ntemp to remove duplicate intersections*/
		ntemp=n;
		do i=2 to n_fam;
			ntemp[,i]=n[,i]#(ntemp[,i-1]^=n_Htest[,i-1])+n_Htest[,i]#(ntemp[,i-1]=n_Htest[,i-1]);
		end;

		inter3=inter2||ntemp||n||m||(1:n_int)`;
		CALL sort(inter3,(1:n_hyp+3*n_fam+1),(1:n_hyp+n_fam));
		interuni=uniqueby(inter3,(1:n_hyp+n_fam-1),1:nrow(inter3));
		inter4=inter3[interuni,1:n_hyp+3*n_fam];
		n=inter4[,n_hyp+n_fam+1:n_hyp+2*n_fam];
		m=inter4[,n_hyp+2*n_fam+1:n_hyp+3*n_fam];
		inter2=inter4[,1:n_hyp];
		inter2=choose(inter2=0,.,inter2);
		n_intuni=nrow(inter2);

		free inter3 inter4 ntemp;
		/***************************************************************/

		/********************* Gamma matrix *********************/
		/*Define the gamma vector*/
		%if &gamma ^= _ALL_ %then %do;
			/*extract gamma vector*/
			pre_gam=concat(tranwrd("&gamma","#"," ")," 1");
			call symput('gamma2',pre_gam);
			free pre_gam;
			gamma= {&gamma2};
			gamnorep=gamma;
		%end;
		%else %do;
			gamma=j(n_intuni*n_gam,n_fam);
			gamnorep=j(n_gam,n_fam);
			i=1;
			%do ngam=1 %to %eval(&n_fam-1);
			do gam&ngam=&gaminf to &gamsup by &gamby;	
			%end;	
			%do ngam=1 %to %eval(&n_fam-1);
				gamma[(i-1)*n_intuni+1:i*n_intuni,&ngam]=gam&ngam;
				gamnorep[i,&ngam]=gam&ngam;
			%end;
			i=i+1;
			
			%do ngam=1 %to %eval(&n_fam-1);
				end;
			%end;
		%end;
		/********************************************************/

		/***** Generate Bonferroni mixing function's critical values: b denominators *****/
		b=j(n_intuni*n_gam,n_fam,.);
		test=j(n_intuni*n_gam,1,1);

		/*Repeat n and m for all gammas*/
		nrep=repeat(n,n_gam,1);
		mrep=repeat(m,n_gam,1);

		b[,1]=(nrep[,1]>0);
		
		do i=2 to n_fam;
			/*Test whether all previous families have no hypothesis*/
			test=test#(nrep[,i-1]=0);
			/*If family i-1 has no hypothesis then b3=1 otherwise b3=the usual b denominator*/
			b1=1-nrep[,i-1];
			b2=(n_Htest[1,i-1]-nrep[,i-1])#(1-gamma[,i-1])/n_Htest[1,i-1];
			b3=b1<>b2;
			/*Apply the formula to calculate b*/
			b[,i]=test#(nrep[,i]^=0)<>b[,i-1]#b3;	
		end;
		free test b1 b2 b3;
		/*Division by 0 not feasible so replace all 0s by missing in b;*/
		b=choose(b=0,.,b);
		/********************************************************************************/

		/******************** Alpha exhaustive gatekeeping procedure *********************/
		%if &exhaust=YES %then %do;
			extmp=j(n_intuni*n_gam,1,0);
			extmp2=j(n_intuni*n_gam,1,1);
			do fam=1 to n_fam;
				extmp=extmp+(mrep[,fam]>0);
			end;
			%if &gamma ^= _ALL_ %then %str(gamma=(extmp>1)#repeat(gamma,n_intuni)+(extmp=1)#extmp2;);
			%else %str(gamma=(extmp>1)#gamma+(extmp=1)#extmp2;);
			free extmp extmp2;
		%end;
		/********************************************************************************/

		/*************** Prepare variables for the gatekeeping procedure ****************/
		/*Repeat variables for all hypotheses*/
		do fam=1 to n_fam;
			if fam=1 then do;
				gamma2=repeat(gamma[,fam],1,n_Htest[,fam]);
				mrep2=repeat(mrep[,fam],1,n_Htest[,fam]);
				n_Htest2=repeat(n_Htest[,fam],1,n_Htest[,fam]);
			end;
			else do;
				gamma2=gamma2||repeat(gamma[,fam],1,n_Htest[,fam]);
				mrep2=mrep2||repeat(mrep[,fam],1,n_Htest[,fam]);
				n_Htest2=n_Htest2||repeat(n_Htest[,fam],1,n_Htest[,fam]);
			end;
			%if exhaust=Yes %then %do;
				extmp=extmp+(mrep2[,fam]>0);
			%end;
		end;
		gamma=gamma2;
		mrep=mrep2;
		n_Htest=n_Htest2;
		/*Transpose raw p-values*/
		t_rawp=rawp`;
		/*******************************************************************************/
		
		/*************** Gatekeeping procedure ******************/
		%if &n_sim ne 0 %then %do; do simulation=1 to &n_sim;%end;
		%else %do;do simulation=1 to 1;%end;
		 
			/*p-values ordering for each family*/
			temp=family`||t_rawp[,simulation]||(1:n_hyp)`;
			call sort(temp,{1 2});
			ord_p=temp[,3]`;
			free temp;

			interraw=rawp[simulation,]#inter2;
	      	interord=j(n_intuni,n_hyp);
			interadj=j(n_intuni*n_gam,n_hyp);
			interp=j(n_intuni*n_gam,n_hyp);
			interm=j(n_gam,n_hyp);
			min_p=j(n_intuni*n_gam,n_fam);
			final_p=j(n_intuni*n_gam,1);
			temp=j(n_intuni*n_gam,n_fam,.);
			min_p_b=j(n_intuni*n_gam,n_fam);
			partemp=j(n_gam,n_hyp,1);
			
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

			/*Repeat k and interord for all gammas*/
			interordr=repeat(interord,n_gam,1);
			krep=repeat(k,n_gam,1);
			free k;
	 
			/*truncated test except for last family*/
			%if &test2=HOMMEL %then %str(interadj=interordr/(krep#gamma/mrep + (1-gamma)/n_Htest););
			%else %if &test2=HOCHBERG %then %str(interadj=interordr/(gamma/(mrep-krep+1) + (1-gamma)/n_Htest););
			%else %if &test2=HOLM %then %str(interadj=interordr/(gamma/mrep + (1-gamma)/n_Htest););
			/*replace all missing by 1*/
			interadj=choose(interadj=.,1,interadj);
			/*Take the minimum of the calculated p-values in each family*/
			do i=1 to n_hyp;
				if i=1 then min_p[,1]=interadj[,1];
				else do;
					if family[1,i] = family[1,i-1]
					then min_p[,family[1,i]]=min_p[,family[1,i]]><interadj[,i];
					else min_p[,family[1,i]]=interadj[,i];
				end;
			end;
			/*divide by bi*/
			min_p_b=min_p/b;
			
			/*Replace missing values by 1*/
			temp2=(min_p_b=temp);
			min_p_b=min_p_b<>temp2;

			/*Compute interaction adjusted p-values*/
			do fam=1 to n_fam;
				final_p[,1]=final_p[,1]><min_p_b[,fam];
			end;
			free interordr krep interadj min_p min_p_b temp2 temp;	
			
		    /*Compute the max over all intersections*/
			interp=final_p[,1]#repeat(inter2,n_gam,1);
			free final_p;
			do intermax=1 to n_hyp;
				do gamval=1 to n_gam;
					interm[gamval,intermax]=max(interp[(gamval-1)*n_intuni+1:gamval*n_intuni,intermax]);		
				end;
			end;

			free interp;
				
			/*Add calculation for the final step (non-consonance for Hommel test)*/
			serial3=serial;
			parallel3=parallel;

			do i=1 to n_hyp;
				/*Serial logical restriction*/
				do while (index(serial3[i,1],"1")>0);
					interm[,i]=interm[,i]<>interm[,index(serial3[i,1],"1")];
					serial3[i,1]=concat(tranwrd(substr(serial3[i,1],1,index(serial3[i,1],"1")),'1','0'),substr(serial3[i,1],index(serial3[i,1],"1")+1));
				end;
				/*parallel logical restriction*/
				if index(parallel3[i,1],"1")>0 then do;
					do while (index(parallel3[i,1],"1")>0);
						partemp[,i]=partemp[,i]><interm[,index(parallel3[i,1],"1")];
						parallel3[i,1]=concat(tranwrd(substr(parallel3[i,1],1,index(parallel3[i,1],"1")),'1','0'),substr(parallel3[i,1],index(parallel3[i,1],"1")+1));
					end;
					interm[,i]=interm[,i]<>partemp[,i];
				end;
			end;
			free partemp serial3 parallel3;
			/********************************************************/
			
			/******* calculate power and FWER for each gamma ********/
			%if &n_sim ne 0 %then %do;

				pow=(interm<=&alpha);
				%if &powout=POWER %then %do;
					%if %sysfunc(exist(POWfunc)) %then %do;
						POWfunc2=j(n_gam,nrow(POWfunc),1);
						powfam=j(n_gam,n_fam,0);
						powfam2=POWfunc[,3:2+2*n_fam];
						powhyp=POWfunc[,2+2*n_fam+1:2+2*n_fam+2*n_hyp];
						do i=1 to n_hyp;
							powfam[,family[,i]]=powfam[,family[,i]]+pow[,i];
						end;

						
						do def=1 to nrow(POWfunc);
							do f=1 to n_fam;
								POWfunc2[,def]=POWfunc2[,def]#((powfam[,f]>=powfam2[def,f+n_fam])#(powfam2[def,f]=1)+(powfam[,f]=powfam2[def,f+n_fam])#(powfam2[def,f]=2));
							end;
							do i=1 to n_hyp;
								POWfunc2[,def]=POWfunc2[,def]#((powhyp[def,n_hyp+i]=pow[,i])#(powhyp[def,i]=2)+(powhyp[def,i]=1));
							end;
						end;
						powerDEF=powerDEF+POWfunc2;
						free powfam powfam2 POWfunc2 powhyp;

					%end;
				%end;


				/*Power for individual hypotheses*/
				indpow=indpow+pow;
				/*FWER*/
				%if &powout=FWER %then %do;
					preFWER=j(n_gam,1,0);
					do i=1 to n_hyp;
						preFWER=preFWER+(pow[,i]#(FWERdef[,i]=0));
					end;
					FWER=FWER+(preFWER>0);
					free preFWER;
				%end;

			%end;
			/********************************************************/

			/**************** Output final datasets *****************/
			%else %do;
				allpval=gamnorep||interm;
				create adj_pIML from allpval;
					append from allpval;
			%end;

			%if &pvalout ne _NULL_ and &n_sim ne 0 %then %do;
				if simulation=1 then &pvalout=interm;
				else &pvalout=&pvalout//interm;
			%end;
			free interm;
		end;

		%if &n_sim ne 0 %then %do;
			%if &powout=POWER %then %do;
				%if %sysfunc(exist(POWfunc)) %then %do;
					do def=1 to ncol(powerDEF);
						power[,alldef[def,1]]=power[,alldef[def,1]]+alldef[def,2]#powerDEF[,def];
					end;
					power=gamnorep||(power*100/&n_sim);	
					create POWERIML from power;
						append from power;
				%end;
				indpow=gamnorep||(indpow*100/&n_sim);
				create indIML from indpow;
					append from indpow;
			%end;
			%else %if &powout=FWER %then %do;
				FWER=gamnorep||(FWER/&n_sim);
				create FWERIML from FWER;
					append from FWER;
			%end;
			%if &pvalout ne _NULL_ %then %do;
				index=repeat((1:&n_sim)`,n_gam,1);
				call sort(index,1);
				&pvalout=index||repeat(gamnorep,&n_sim)||&pvalout;
				create &pvalout from &pvalout;
					append from &pvalout;
				free index;
			%end;
		%end;

		quit;

		/******************************************* proc IML ends *******************************************/

		/******************************************* Report results ******************************************/
		%report;
		proc datasets lib=work memtype=data;
			save study %if %sysfunc(exist(POWfunc)) %then POWfunc;
				 %if &powout=POWER and &n_sim ne 0 and %sysfunc(exist(POWfunc)) %then power;
				 %if &powout=POWER and &n_sim ne 0 %then indpower; %if &powout=FWER and &n_sim ne 0 %then FWER; 
				 %if %sysfunc(exist(FWERdef)) %then FWERdef; %if &n_sim=0 %then adjp; %if &n_sim ne 0 and &pvalout ne _NULL_ %then &pvalout;
				 %if &n_sim ne 0 %then &rawp;;
		run;
		quit;
	%end;
		/*****************************************************************************************************/
%mend;
