#define _USE_MATH_DEFINES
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_vaccine_oney_maps(NumericMatrix in_Comp, List parm, List siaparm, NumericMatrix beta_full,
	NumericVector pop_full, NumericVector cov1, double cov2, int t_start)
{
	/* This function runs measles transmission, ageing, and vaccination at a specific age group and timestep
	during calendar years for simulation. */

	// =================================================
	// Set-up (index, parameters)
	// =================================================
	int i_M = 0;    // Maternal immune
	int i_S = 1;    // Susceptible
	int i_I = 2;    // Infectious
	int i_R = 3;    // Recovered
	int i_V1S = 4;  // Vaccinated susceptible - 1 dose
	int i_V1I = 5;  // Vaccinated infectious  - 1 dose
	int i_V1R = 6;  // Vaccinated recovered   - 1 dose
	int i_V2S = 7;  // Vaccinated susceptible - 2 dose
	int i_V2I = 8;  // Vaccinated infectious  - 2 dose
	int i_V2R = 9;  // Vaccinated recovered   - 2 dose
	int i_V3S = 10; // Vaccinated susceptible - 3 dose
	int i_V3I = 11; // Vaccinated infectious  - 3 dose
	int i_V3R = 12; // Vaccinated recovered   - 3 dose
	int i_V1F = 13; // Vaccinated population  - 1 dose (counter)

	List outp;                           // output list, including compartment distribution and cases
	NumericMatrix trans_Comp(254, 14);   // compartments after including transmission, 254 age groups, 14 states
	NumericMatrix out_Comp(254, 14);     // compartments after including transmission and ageing, 254 age groups, 14 states
	NumericVector betta(254);            // contact rate reported by a contactor of age a, 254 age groups
	NumericVector cyc(254);              // case prevalence (adjusted for seasonality), 254 age groups
	NumericVector newinfect(254);        // new infections/cases, 254 age groups
	NumericVector newinfect0d(254);      // new zero-dose infections/cases, 254 age groups
	NumericVector newdoseRI1(254);       // MCV1 doses administrated through routine immunisation, 254 age groups
	NumericVector newdoseRI2(254);       // MCV2 doses administrated through routine immunisation, 254 age groups
	NumericVector newdoseSIAc1(254);     // SIA doses administrated through campaign with strategy C to zero-dose, 254 age groups
	NumericVector newdoseSIAc2(254);     // SIA doses administrated through campaign with strategy C to already-vaccinated, 254 age groups
	NumericVector newdoseSIAde1(254);    // SIA doses administrated through strategies D & E to zero-dose, 254 age groups
	NumericVector newdoseSIAde2(254);    // SIA doses administrated through strategies D & E to already-vaccinated, 254 age groups
	NumericVector newdoseSIAf1(254);     // SIA doses administrated through campaign with strategy F to zero-dose, 254 age groups
	NumericVector newdoseSIAf2(254);     // SIA doses administrated through campaign with strategy F to already-vaccinated, 254 age groups
	
	double tcycle = 0.0;                 // seasonality
	double lambda = 0.0;                 // force of infection
	double pop_fert_SR = 0.0;            // population of S and R at fertility age
	double pop_fert_R = 0.0;             // population of R at fertility age
	double prp_R = 0.0;                  // proportion of being born with maternal immunity
	double n_eff = 0.0;                  // population with effective protection of MCV1
	double n_v1 = 0.0;                   // population who receive MCV1
	double p_eff = 0.0;                  // proportion of effective vaccine protection among those receive MCV1
	double ve2 = 0.0;                    // 2nd dose vaccine efficay conditioned on MCV1

	double tstep = as<double>(parm["tstep"]);                // timesteps per year
	double gamma = as<double>(parm["gamma"]);                // recovery rate per timestep
	double amp = as<double>(parm["amp"]);	                 // amplification for seasonality
	NumericVector ve1 = as<NumericVector>(parm["ve1"]);      // vaccine efficacy of the first dose, 254 age groups
	double ve2plus = as<double>(parm["ve2plus"]);            // vaccine protection for two doses
	double age_w = 52/tstep;                                 // weekly ageing rate per timestep
	double age_y = 1/tstep;  				                 // annualy ageing rate per timestep
	double wane = (12/6)/tstep;  		                     // waning rate maternal immunity per timestep (duration of 6 months)
  
    int sia_rounds = as<int>(siaparm["sia_rounds"]);                  // number of SIA rounds in a single year
	IntegerVector alla0 = as<IntegerVector>(siaparm["a0"]);           // starting target age groups
	IntegerVector alla1 = as<IntegerVector>(siaparm["a1"]);           // ending target age groups
	NumericVector allsiacov = as<NumericVector>(siaparm["siacov"]);   // SIA coverage among total population
	IntegerVector allpoptype = as<IntegerVector>(siaparm["poptype"]); // type of targetd population based on strategy: 3-C, 4-D, 5-E, 6-F


    for (int t = t_start; t <= (t_start+tstep); ++t)
	{
		// =================================================
        // Transmission
        // =================================================
        tcycle = 1.0 + amp*sin(2.0*M_PI*t/tstep);     //Seasonality in Beta, M_PI = 3.14159265358979323846
        cyc = tcycle*(in_Comp(_,i_I) + in_Comp(_,i_V1I) + in_Comp(_,i_V2I) + in_Comp(_,i_V3I) + 1.0e-9);

		for (int a = 0; a < 254; ++a)
		{
			betta = beta_full(a,_);
            lambda = 1.0 - exp(-sum(betta*cyc));
            newinfect[a] += lambda*(in_Comp(a,i_S) + in_Comp(a,i_V1S) + in_Comp(a,i_V2S) + in_Comp(a,i_V3S));
			newinfect0d[a] += lambda*(in_Comp(a,i_S));
            //Rcout << "Age group = " << a+1 << "\n" // print out FOI to check
            //      << "FOI = " << lambda << "\n";

            trans_Comp(a,i_M)   = in_Comp(a,i_M);
            trans_Comp(a,i_S)   = in_Comp(a,i_S)   - lambda*in_Comp(a,i_S);
            trans_Comp(a,i_I)   = in_Comp(a,i_I)   + lambda*in_Comp(a,i_S)   - gamma*in_Comp(a,i_I);
            trans_Comp(a,i_R)   = in_Comp(a,i_R)   + gamma *in_Comp(a,i_I);
            trans_Comp(a,i_V1S) = in_Comp(a,i_V1S) - lambda*in_Comp(a,i_V1S);
            trans_Comp(a,i_V1I) = in_Comp(a,i_V1I) + lambda*in_Comp(a,i_V1S) - gamma*in_Comp(a,i_V1I);
            trans_Comp(a,i_V1R) = in_Comp(a,i_V1R) + gamma *in_Comp(a,i_V1I);
            trans_Comp(a,i_V2S) = in_Comp(a,i_V2S) - lambda*in_Comp(a,i_V2S);
            trans_Comp(a,i_V2I) = in_Comp(a,i_V2I) + lambda*in_Comp(a,i_V2S) - gamma*in_Comp(a,i_V2I);
            trans_Comp(a,i_V2R) = in_Comp(a,i_V2R) + gamma *in_Comp(a,i_V2I);
            trans_Comp(a,i_V3S) = in_Comp(a,i_V3S) - lambda*in_Comp(a,i_V3S);
            trans_Comp(a,i_V3I) = in_Comp(a,i_V3I) + lambda*in_Comp(a,i_V3S) - gamma*in_Comp(a,i_V3I);
            trans_Comp(a,i_V3R) = in_Comp(a,i_V3R) + gamma *in_Comp(a,i_V3I);
            trans_Comp(a,i_V1F) = in_Comp(a,i_V1F);
	    }
        //Rcout << "\ntrnansmission cycle completed\n";


        // =================================================
        // Ageing and routine vaccination (MCV1, MCV2)
        // =================================================
        //NumericMatrix out_Comp(254, 14);         // output compartments after including transmission and ageing, 254 age groups, 14 states

        // ageing process: 3-100 years old, no vaccination
        for (int a = 253; a > 155; --a)
		{
			out_Comp(a,i_M)   = trans_Comp(a,i_M)   - age_y*trans_Comp(a,i_M)   + age_y*trans_Comp(a-1,i_M)   - wane*trans_Comp(a,i_M);
            out_Comp(a,i_S)   = trans_Comp(a,i_S)   - age_y*trans_Comp(a,i_S)   + age_y*trans_Comp(a-1,i_S)   + wane*trans_Comp(a,i_M);
            out_Comp(a,i_I)   = trans_Comp(a,i_I)   - age_y*trans_Comp(a,i_I)   + age_y*trans_Comp(a-1,i_I)  ;
            out_Comp(a,i_R)   = trans_Comp(a,i_R)   - age_y*trans_Comp(a,i_R)   + age_y*trans_Comp(a-1,i_R)  ;
            out_Comp(a,i_V1S) = trans_Comp(a,i_V1S) - age_y*trans_Comp(a,i_V1S) + age_y*trans_Comp(a-1,i_V1S);
            out_Comp(a,i_V1I) = trans_Comp(a,i_V1I) - age_y*trans_Comp(a,i_V1I) + age_y*trans_Comp(a-1,i_V1I);
            out_Comp(a,i_V1R) = trans_Comp(a,i_V1R) - age_y*trans_Comp(a,i_V1R) + age_y*trans_Comp(a-1,i_V1R);
            out_Comp(a,i_V2S) = trans_Comp(a,i_V2S) - age_y*trans_Comp(a,i_V2S) + age_y*trans_Comp(a-1,i_V2S);
            out_Comp(a,i_V2I) = trans_Comp(a,i_V2I) - age_y*trans_Comp(a,i_V2I) + age_y*trans_Comp(a-1,i_V2I);
            out_Comp(a,i_V2R) = trans_Comp(a,i_V2R) - age_y*trans_Comp(a,i_V2R) + age_y*trans_Comp(a-1,i_V2R);
            out_Comp(a,i_V3S) = trans_Comp(a,i_V3S) - age_y*trans_Comp(a,i_V3S) + age_y*trans_Comp(a-1,i_V3S);
            out_Comp(a,i_V3I) = trans_Comp(a,i_V3I) - age_y*trans_Comp(a,i_V3I) + age_y*trans_Comp(a-1,i_V3I);
            out_Comp(a,i_V3R) = trans_Comp(a,i_V3R) - age_y*trans_Comp(a,i_V3R) + age_y*trans_Comp(a-1,i_V3R);
            out_Comp(a,i_V1F) = trans_Comp(a,i_V1F) - age_y*trans_Comp(a,i_V1F) + age_y*trans_Comp(a-1,i_V1F);
            //Rcout << a+1 << " ";
        }
        //Rcout << "\nYearly age groups done\n";

        // ageing process and routine vaccination: 2 weeks - 2 years old
        for (int a = 155; a > 0; --a)
		{
			if (a == 71)
			{
				// adjust MCV2 coverage as it applies to already-vaccinated population
				double p_allpop = sum(trans_Comp(71-1,_)) - trans_Comp(71-1,i_V1F);                                                    //proportion of total population
				double p_vaced = 1.0 - ((trans_Comp(71-1,i_M) + trans_Comp(71-1,i_S) + trans_Comp(71-1,i_I) + trans_Comp(71-1,i_R))/p_allpop); //proportion of already-vaccianted population
				double adjcov2 = 0.0;  // adjusted MCV2 coverage
				double surcov1 = cov1[71-1] ;  // MCV1 coverage and surplus MCV2 coverage delivered as MCV1
				if (cov2 > 0.0) 
				{
					if (cov2 <= p_vaced) 
					{
						adjcov2 = cov2/p_vaced;
					}
					else
					{
						adjcov2 = 1.0; 
						if (p_allpop > p_vaced) {surcov1 = cov1[71-1] + (cov2-p_vaced)/(p_allpop-p_vaced);} // surplus MCV2 doses given to zero-dose population
						if (surcov1 > 1.0){surcov1 = 1.0;}
					}
				//if (t == (t_start + tstep/2)) {Rcout << "MCV2 coverage: " << cov2 << " -> " << adjcov2 << ", MCV1 = " << surcov1 << ", % vaccinated: " << p_vaced << "\n";}
				}

				// MCV doses administrated				
				newdoseRI2[71] += age_w*cov2*p_allpop;  // include surplus MCV2 doses 
				
				// ageing process and routine vaccination: 72 weeks old
				out_Comp(71,i_M)   = trans_Comp(71,i_M)
								   - age_w*trans_Comp(71,i_M)
								   + age_w*trans_Comp(71-1,i_M)*(1.0-surcov1) // under 'perfect timeliness' assumption, cov1 = 0 for all ages other than 39 weeks old
								   - wane *trans_Comp(71,i_M);

				out_Comp(71,i_S)   = trans_Comp(71,i_S)
								   - age_w*trans_Comp(71,i_S)
								   + age_w*trans_Comp(71-1,i_S)*(1.0-surcov1)
								   + wane *trans_Comp(71,i_M);

				out_Comp(71,i_I)   = trans_Comp(71,i_I)
								   - age_w*trans_Comp(71,i_I)
								   + age_w*trans_Comp(71-1,i_I)*(1.0-surcov1);

				out_Comp(71,i_R)   = trans_Comp(71,i_R)
								   - age_w*trans_Comp(71,i_R)
								   + age_w*trans_Comp(71-1,i_R)*(1.0-surcov1);

				out_Comp(71,i_V1S) = trans_Comp(71,i_V1S)
								   - age_w*trans_Comp(71,i_V1S)
								   + age_w*trans_Comp(71-1,i_V1S)*(1.0-adjcov2)
								   + age_w*(trans_Comp(71-1,i_M)+trans_Comp(71-1,i_S))*surcov1*(1.0-ve1[71-1]);

				out_Comp(71,i_V1I) = trans_Comp(71,i_V1I)
								   - age_w*trans_Comp(71,i_V1I)
								   + age_w*trans_Comp(71-1,i_V1I)*(1.0-adjcov2)
								   + age_w*trans_Comp(71-1,i_I)*surcov1;

				out_Comp(71,i_V1R) = trans_Comp(71,i_V1R)
								   - age_w*trans_Comp(71,i_V1R)
								   + age_w*trans_Comp(71-1,i_V1R)*(1.0-adjcov2)
								   + age_w*(trans_Comp(71-1,i_M)+trans_Comp(71-1,i_S))*surcov1*ve1[71-1]
								   + age_w*trans_Comp(71-1,i_R)*surcov1;

				out_Comp(71,i_V1F) = trans_Comp(71,i_V1F)
								   - age_w*trans_Comp(71,i_V1F)
								   + age_w*trans_Comp(71-1,i_V1F)*(1.0-adjcov2)
								   + age_w*(trans_Comp(71-1,i_M)+trans_Comp(71-1,i_S))*surcov1*ve1[71-1];

				n_eff = out_Comp(71,i_V1F);  // number of children effectively protected by MCV1
				n_v1 = out_Comp(71,i_V1S) + out_Comp(71,i_V1I) + out_Comp(71,i_V1R);  // number of children received MCV1
				if (n_v1 > 0.0)
				{
					p_eff = n_eff/n_v1;  // proportion of effective protection among children received MCV1
					if (p_eff > ve2plus) {ve2 = 0.0;}  // vaccine efficacy of 2nd dose conditioned on 1st dose
					else {ve2 = (ve2plus - p_eff)/(1.0 - p_eff);}
				}
				else
				{
					ve2 = ve2plus;
				}
				//Rcout << "\n2nd dose efficacy: " << ve2 << "\n ";

				out_Comp(71,i_V2S) = trans_Comp(71,i_V2S)
								   - age_w*trans_Comp(71,i_V2S)
								   + age_w*trans_Comp(71-1,i_V2S)*(1.0-adjcov2)
								   + age_w*trans_Comp(71-1,i_V1S)*adjcov2*(1.0-ve2);

				out_Comp(71,i_V2I) = trans_Comp(71,i_V2I)
								   - age_w*trans_Comp(71,i_V2I)
								   + age_w*trans_Comp(71-1,i_V2I)*(1.0-adjcov2)
								   + age_w*trans_Comp(71-1,i_V1I)*adjcov2;

				out_Comp(71,i_V2R) = trans_Comp(71,i_V2R)
								   - age_w*trans_Comp(71,i_V2R)
							   	   + age_w*trans_Comp(71-1,i_V2R)*(1.0-adjcov2)
								   + age_w*trans_Comp(71-1,i_V1S)*adjcov2*ve2
								   + age_w*trans_Comp(71-1,i_V1R)*adjcov2;

				double ve3 = ve2;  // assume 3+ dose efficacy (0.0 or equal to ve2)
				
				out_Comp(71,i_V3S) = trans_Comp(71,i_V3S)
								   - age_w*trans_Comp(71,i_V3S)
								   + age_w*trans_Comp(71-1,i_V3S)*(1.0-adjcov2*ve3)   // (1.0-adjcov2 + adjcov2*(1.0-ve3))
								   + age_w*trans_Comp(71-1,i_V2S)*adjcov2*(1.0-ve3);

				out_Comp(71,i_V3I) = trans_Comp(71,i_V3I)
								   - age_w*trans_Comp(71,i_V3I)
								   + age_w*trans_Comp(71-1,i_V3I)
								   + age_w*trans_Comp(71-1,i_V2I)*adjcov2;

				out_Comp(71,i_V3R) = trans_Comp(71,i_V3R)
								   - age_w*trans_Comp(71,i_V3R)
								   + age_w*trans_Comp(71-1,i_V3R)
								   + age_w*trans_Comp(71-1,i_V2R)*adjcov2
								   + age_w*trans_Comp(71-1,i_V2S)*adjcov2*ve3
								   + age_w*trans_Comp(71-1,i_V3S)*adjcov2*ve3;
			}
			else
			{
				// MCV doses administrated				
				newdoseRI1[a] += age_w*cov1[a-1]*(trans_Comp(a-1,i_M)
										        + trans_Comp(a-1,i_S)
										        + trans_Comp(a-1,i_I)
										        + trans_Comp(a-1,i_R));
				
                //Rcout << a+1 << " doseRI1:" << newdoseRI1[a] << "\n ";
				
				out_Comp(a,i_M)   = trans_Comp(a,i_M)
			                      - age_w*trans_Comp(a,i_M)
                			      + age_w*trans_Comp(a-1,i_M)*(1.0-cov1[a-1])
           	    			      - wane*trans_Comp(a,i_M);

           	    out_Comp(a,i_S)   = trans_Comp(a,i_S)
                                  - age_w*trans_Comp(a,i_S)
                                  + age_w*trans_Comp(a-1,i_S)*(1.0-cov1[a-1])
           	    			      + wane*trans_Comp(a,i_M);

           	    out_Comp(a,i_I)   = trans_Comp(a,i_I)
                                  - age_w*trans_Comp(a,i_I)
           	    			      + age_w*trans_Comp(a-1,i_I)*(1.0-cov1[a-1]);

           	    out_Comp(a,i_R)   = trans_Comp(a,i_R)
                                  - age_w*trans_Comp(a,i_R)
                			      + age_w*trans_Comp(a-1,i_R)*(1.0-cov1[a-1]);

           	    out_Comp(a,i_V1S) = trans_Comp(a,i_V1S)
                                  - age_w*trans_Comp(a,i_V1S)
                	              + age_w*trans_Comp(a-1,i_V1S)
                	              + age_w*(trans_Comp(a-1,i_M)+trans_Comp(a-1,i_S))*cov1[a-1]*(1.0-ve1[a-1]);

           	    out_Comp(a,i_V1I) = trans_Comp(a,i_V1I)
                                  - age_w*trans_Comp(a,i_V1I)
                	              + age_w*trans_Comp(a-1,i_V1I)
                	              + age_w*trans_Comp(a-1,i_I)*cov1[a-1];

           	    out_Comp(a,i_V1R) = trans_Comp(a,i_V1R)
                                  - age_w*trans_Comp(a,i_V1R)
                	              + age_w*trans_Comp(a-1,i_V1R)
                	              + age_w*(trans_Comp(a-1,i_M)+trans_Comp(a-1,i_S))*cov1[a-1]*ve1[a-1]
           	    	              + age_w*trans_Comp(a-1,i_R)*cov1[a-1];

           	    out_Comp(a,i_V1F) = trans_Comp(a,i_V1F)
                                  - age_w*trans_Comp(a,i_V1F)
                                  + age_w*trans_Comp(a-1,i_V1F)
                	              + age_w*(trans_Comp(a-1,i_M)+trans_Comp(a-1,i_S))*cov1[a-1]*ve1[a-1];

           	    out_Comp(a,i_V2S) = trans_Comp(a,i_V2S) - age_w*trans_Comp(a,i_V2S) + age_w*trans_Comp(a-1,i_V2S);
                out_Comp(a,i_V2I) = trans_Comp(a,i_V2I) - age_w*trans_Comp(a,i_V2I) + age_w*trans_Comp(a-1,i_V2I);
                out_Comp(a,i_V2R) = trans_Comp(a,i_V2R) - age_w*trans_Comp(a,i_V2R) + age_w*trans_Comp(a-1,i_V2R);
                out_Comp(a,i_V3S) = trans_Comp(a,i_V3S) - age_w*trans_Comp(a,i_V3S) + age_w*trans_Comp(a-1,i_V3S);
                out_Comp(a,i_V3I) = trans_Comp(a,i_V3I) - age_w*trans_Comp(a,i_V3I) + age_w*trans_Comp(a-1,i_V3I);
                out_Comp(a,i_V3R) = trans_Comp(a,i_V3R) - age_w*trans_Comp(a,i_V3R) + age_w*trans_Comp(a-1,i_V3R);
			}
		}

		// ageing process: 1 week at birth
		pop_fert_SR = 0.0;
		pop_fert_R = 0.0;
		for (int a = 171; a < 189; ++a)  // fertility age: 18-35 years old
		{
			pop_fert_SR += (1.0 - out_Comp(a,i_I) - out_Comp(a,i_V1I) - out_Comp(a,i_V2I) - out_Comp(a,i_V3I))*pop_full[a];
			pop_fert_R += (out_Comp(a,i_R) + out_Comp(a,i_V1R) + out_Comp(a,i_V2R) + out_Comp(a,i_V3R))*pop_full[a];
		}
		prp_R = pop_fert_R/pop_fert_SR;  // Proportion of being born with maternal immunity
		
		//Rcout << "Proportion of immune = " << prp_R*100 << "%%\n";

		out_Comp(0,i_M)   = trans_Comp(0,i_M)   - age_w*trans_Comp(0,i_M)   + age_w*prp_R;
		out_Comp(0,i_S)   = trans_Comp(0,i_S)   - age_w*trans_Comp(0,i_S)   + age_w*(1.0-prp_R);
		out_Comp(0,i_I)   = trans_Comp(0,i_I)   - age_w*trans_Comp(0,i_I)  ;
		out_Comp(0,i_R)   = trans_Comp(0,i_R)   - age_w*trans_Comp(0,i_R)  ;
		out_Comp(0,i_V1S) = trans_Comp(0,i_V1S) - age_w*trans_Comp(0,i_V1S);
		out_Comp(0,i_V1I) = trans_Comp(0,i_V1I) - age_w*trans_Comp(0,i_V1I);
		out_Comp(0,i_V1R) = trans_Comp(0,i_V1R) - age_w*trans_Comp(0,i_V1R);
		out_Comp(0,i_V2S) = trans_Comp(0,i_V2S) - age_w*trans_Comp(0,i_V2S);
		out_Comp(0,i_V2I) = trans_Comp(0,i_V2I) - age_w*trans_Comp(0,i_V2I);
		out_Comp(0,i_V2R) = trans_Comp(0,i_V2R) - age_w*trans_Comp(0,i_V2R);
		out_Comp(0,i_V3S) = trans_Comp(0,i_V3S) - age_w*trans_Comp(0,i_V3S);
		out_Comp(0,i_V3I) = trans_Comp(0,i_V3I) - age_w*trans_Comp(0,i_V3I);
		out_Comp(0,i_V3R) = trans_Comp(0,i_V3R) - age_w*trans_Comp(0,i_V3R);
		out_Comp(0,i_V1F) = trans_Comp(0,i_V1F) - age_w*trans_Comp(0,i_V1F);


		// =================================================
		// SIA in the middle of a year
		// =================================================
		// spread SIA runs into different timesteps after reaching the middle of a year
		if ((sia_rounds > 0) && (t >= (t_start + tstep/2))) 
		{
			if ((t - (t_start + tstep/2)) < sia_rounds) 
			{
				int ird = t - (t_start + tstep/2);  // timesteps counting from the beginning of the year
				int a0 = alla0[ird];                // starting target age group for a specific SIA round
				int a1 = alla1[ird];                // ending target age group for a specific SIA round
				double siacov = allsiacov[ird];     // SIA coverage among total population for a specific SIA round
				int poptype = allpoptype[ird];      // type of targetd population: 0 - genereal, 1 - zero-dose, 2 - alreday-vaccinated		
				in_Comp = clone(out_Comp);          // temporary compartments after including transmission and routine vaccination
				
				for (int a = (a0-1); a < a1; ++a)
				{
					double siadose   = siacov*pop_full[a];                                                            // number of total SIA doses
					double p_allpop  = sum(in_Comp(a,_)) - in_Comp(a,i_V1F);                                          // proportion of total population
				    double p_0dose   = (in_Comp(a,i_M) + in_Comp(a,i_S) + in_Comp(a,i_I) + in_Comp(a,i_R))/p_allpop ; // proportion of already-vaccianted population
					double pop_0dose = p_0dose*pop_full[a];                                                           // number of zero-dose population
					double pop_vaced = (1.0 - p_0dose)*pop_full[a];                                                   // number of already-vaccinated population
					
					double siacov1 = 0.0, siacov2 = 0.0;  // SIA coverages for zero-dose and already-vaccinated populations
					
					if (pop_0dose < 0.5) {pop_0dose = 0.5;}
					if (pop_vaced < 0.5) {pop_vaced = 0.5;}	
					// Rcout << "timestep = " << t << ", age group = " << a0 << "-" << a1 << "\n ";
					// Rcout << "populations: " << pop_0dose << " (zero-dose) " << pop_vaced << " (vaccinated)\n";
		
					if (poptype == 3 || poptype == 6)  // Strategy C (regular) or F (one-time off)
					{						
						// SIA distribution based on the assumption of 7.7% never reached
						double pr_0dose = pop_0dose/(pop_0dose + pop_vaced);  // proportion of zero-dose population
						
						if ((siacov < 0.923) && (pr_0dose > 0.077))
						{ // doses given randomly to the population except for the never-reached 
							siacov1 = siadose*((pr_0dose-0.077)/(1-0.077))/pop_0dose;
							siacov2 = siadose*((1-pr_0dose)/(1-0.077))/pop_vaced;						
						}
						else 
						{ // doses first given to already-vaccinated and then to zero-dose populations				
							if (siadose > pop_vaced)
							{
								siacov1 = (siadose-pop_vaced)/pop_0dose; 
								siacov2 = 1.0;
							} 
							else
							{
								siacov1 = 0.0; 
								siacov2 = siadose/pop_vaced;
							}
						}
					}
					else
					{
						if (poptype == 4)  // Strategy D (HTR & MOV, first-dose)
						{
							if (siadose <= pop_0dose)
							{
								siacov1 = siadose/pop_0dose;
								siacov2 = 0.0;
							} else
							{
								siacov1 = 1.0;
								siacov2 = (siadose-pop_0dose)/pop_vaced;
							}
						} 
						else              // Strategy E (HTR & MOV, second-dose)
						{
							if (siadose <= pop_vaced)
							{
								siacov1 = 0.0;
								siacov2 = siadose/pop_vaced;
							} else 
							{
								siacov1 = (siadose-pop_vaced)/pop_0dose;
								siacov2 = 1.0;
							}
						}
					}
					//Rcout << "SIA strategy " << poptype << "\n";
					//Rcout << "SIA coverages: " << siacov1 << " (dose 1) " << siacov2 << " (dose 2)\n";
		
					if (siacov1 > 1.0) {siacov1 = 1.0;}
					if (siacov2 > 1.0) {siacov2 = 1.0;}
	
					double n_eff = 0.0, n_v1 = 0.0, p_eff = 0.0, ve2 = 0.0;  // initialise for later calculations

					// MCV doses administrated
				    if (poptype == 3)
				    {
				    	newdoseSIAc1[a] += siacov1*(in_Comp(a,i_M)
				    					 		  + in_Comp(a,i_S)
				    							  + in_Comp(a,i_I)
				    							  + in_Comp(a,i_R));
				    	newdoseSIAc2[a] += siacov2*(in_Comp(a,i_V1S)
				    							  + in_Comp(a,i_V1I)
				    							  + in_Comp(a,i_V1R)
				    							  + in_Comp(a,i_V2S)
				    					 		  + in_Comp(a,i_V2I)
				    							  + in_Comp(a,i_V2R)
												  + in_Comp(a,i_V3S)
												  + in_Comp(a,i_V3I)
												  + in_Comp(a,i_V3R));
						/*if (a == 52){
							Rcout << "Number of SIA dose C, one year old: " << newdoseSIAc1[a] + newdoseSIAc2[a]<< "\n";
						}*/	
				    }			
				    else
				    {   
				        if (poptype == 6) 
						{
							newdoseSIAf1[a] += siacov1*(in_Comp(a,i_M)
													  + in_Comp(a,i_S)
													  + in_Comp(a,i_I)
													  + in_Comp(a,i_R));
				    		newdoseSIAf2[a] += siacov2*(in_Comp(a,i_V1S)
													  + in_Comp(a,i_V1I)
													  + in_Comp(a,i_V1R)
													  + in_Comp(a,i_V2S)
													  + in_Comp(a,i_V2I)
													  + in_Comp(a,i_V2R)
													  + in_Comp(a,i_V3S)
													  + in_Comp(a,i_V3I)
													  + in_Comp(a,i_V3R));
						} 
						else 
						{ // poptype == 4 or poptype == 5
							newdoseSIAde1[a] += siacov1*(in_Comp(a,i_M)  
													   + in_Comp(a,i_S)
													   + in_Comp(a,i_I)
													   + in_Comp(a,i_R));
				    		newdoseSIAde2[a] += siacov2*(in_Comp(a,i_V1S) 
													   + in_Comp(a,i_V1I)
													   + in_Comp(a,i_V1R)
													   + in_Comp(a,i_V2S)
													   + in_Comp(a,i_V2I)
													   + in_Comp(a,i_V2R)
													   + in_Comp(a,i_V3S)
													   + in_Comp(a,i_V3I)
													   + in_Comp(a,i_V3R));
						}
				    }
					
					out_Comp(a,i_M)     = in_Comp(a,i_M)
										- in_Comp(a,i_M)*siacov1;
	
					out_Comp(a,i_S)     = in_Comp(a,i_S)
										- in_Comp(a,i_S)*siacov1;
	
					out_Comp(a,i_I)     = in_Comp(a,i_I)
										- in_Comp(a,i_I)*siacov1;
	
					out_Comp(a,i_R)     = in_Comp(a,i_R)
										- in_Comp(a,i_R)*siacov1;
	
					out_Comp(a,i_V1S)   = in_Comp(a,i_V1S)
										- in_Comp(a,i_V1S)*siacov2
										+ (in_Comp(a,i_M)+in_Comp(a,i_S))*siacov1*(1.0-ve1[a]);
	
					out_Comp(a,i_V1I)   = in_Comp(a,i_V1I)
										- in_Comp(a,i_V1I)*siacov2
										+ in_Comp(a,i_I)*siacov1;
	
					out_Comp(a,i_V1R)   = in_Comp(a,i_V1R)
										- in_Comp(a,i_V1R)*siacov2
										+ in_Comp(a,i_R)*siacov1
										+ (in_Comp(a,i_M)+in_Comp(a,i_S))*siacov1*ve1[a];
	
					out_Comp(a,i_V1F)   = in_Comp(a,i_V1F)
										- in_Comp(a,i_V1F)*siacov2
										+ (in_Comp(a,i_M)+in_Comp(a,i_S))*siacov1*ve1[a];
	
					n_eff = out_Comp(a,i_V1F);  // number of children effectively protected by MCV1
					n_v1 = out_Comp(a,i_V1S) + out_Comp(a,i_V1I) + out_Comp(a,i_V1R);  // number of children received MCV1
					if (n_v1 > 0.0)
					{
						p_eff = n_eff/n_v1;  // proportion of effective protection among children received MCV1
						if (p_eff > ve2plus) {ve2 = 0.0;}  // vaccine efficacy of 2nd dose conditioned on 1st dose
						else {ve2 = (ve2plus - p_eff)/(1.0 - p_eff);}
					}
					else
					{
						ve2 = ve2plus;
					}
	
					//Rcout << "\n2nd dose efficacy: " << ve2 << "\n ";
	
					out_Comp(a,i_V2S)   = in_Comp(a,i_V2S)
										- in_Comp(a,i_V2S)*siacov2
										+ in_Comp(a,i_V1S)*siacov2*(1.0-ve2);
	
					out_Comp(a,i_V2I)   = in_Comp(a,i_V2I)
										- in_Comp(a,i_V2I)*siacov2
										+ in_Comp(a,i_V1I)*siacov2;
	
					out_Comp(a,i_V2R)   = in_Comp(a,i_V2R)
										- in_Comp(a,i_V2R)*siacov2
										+ in_Comp(a,i_V1S)*siacov2*ve2
										+ in_Comp(a,i_V1R)*siacov2;
	
					double ve3 = ve2;  // assume 3+ dose efficacy equal to ve2
					out_Comp(a,i_V3S)   = in_Comp(a,i_V3S)
					                    - in_Comp(a,i_V3S)*siacov2
										+ in_Comp(a,i_V3S)*siacov2*(1.0-ve3)
										+ in_Comp(a,i_V2S)*siacov2*(1.0-ve3); 
	
					out_Comp(a,i_V3I)   = in_Comp(a,i_V3I)
										+ in_Comp(a,i_V2I)*siacov2;     // - in_Comp(a,i_V3I)*siacov2 + in_Comp(a,i_V3I)*siacov2
	
					out_Comp(a,i_V3R)   = in_Comp(a,i_V3R)
										+ in_Comp(a,i_V2R)*siacov2
										+ in_Comp(a,i_V2S)*siacov2*ve3
										+ in_Comp(a,i_V3S)*siacov2*ve3; // - in_Comp(a,i_V3R)*siacov2 + in_Comp(a,i_V3R)*siacov2	
					// Rcout << a+1 << " ";
				}
			}
		}
		in_Comp = clone(out_Comp);
    }
	
	outp["cases"] = newinfect;
	outp["cases0d"] = newinfect0d;
	outp["doseRI1"] = newdoseRI1;
	outp["doseRI2"] = newdoseRI2;
	outp["doseSIAc1"] = newdoseSIAc1;
	outp["doseSIAc2"] = newdoseSIAc2;
	outp["doseSIAde1"] = newdoseSIAde1;
	outp["doseSIAde2"] = newdoseSIAde2;
	outp["doseSIAf1"] = newdoseSIAf1;
	outp["doseSIAf2"] = newdoseSIAf2;
	outp["out_Comp"] = in_Comp;
	
	return outp;
}
