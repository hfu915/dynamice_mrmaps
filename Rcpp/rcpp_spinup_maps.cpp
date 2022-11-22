#define _USE_MATH_DEFINES
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List rcpp_spinup_maps(NumericMatrix in_Comp, List parm, NumericMatrix beta_full, NumericVector pop_full, int spinup_yr) {
	
	/* This function runs measles transmission and ageing for each age group during spin-up period. */

    // ================================================= 
    // Set-up (index, parameters)
    // =================================================
    int i_M = 0;    // Maternal Immune
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
    int i_V1F = 13; //  Vaccinated population  - 1 dose (counter) 
    
    NumericMatrix trans_Comp(254, 14);   // output compartments after including transmission, 254 age groups, 14 states 
	NumericVector betta(254);            // contact rate reported by a contactor of age a, 254 age groups
    NumericVector cyc(254);              // case prevalence (adjusted for seasonality), 254 age groups
	double tcycle = 0.0;                 // seasonality 
    double lambda = 0.0;                 // force of infection
    double pop_fert_SR = 0.0;            // population of S and R at fertility age
    double pop_fert_R = 0.0;             // population of R at fertility age
    double prp_R = 0.0;                  // proportion of being born with maternal immunity
    
    double tstep = as<double>(parm["tstep"]);
    double gamma = as<double>(parm["gamma"]);
    double amp = as<double>(parm["amp"]);	
    double age_w = 52/tstep;               // ageing - weekly
    double age_y = 1/tstep;  			   // ageing - annualy
    double wane = (12/6)/tstep;  		   // waning maternal immunity (6 months) - per timestep
	
    List outp;                                // output list, including compartments and cases
	NumericMatrix outp_incd(254, spinup_yr);  // output number of new infections at end of each year, 254 age groups
	NumericMatrix outp_susc(254, spinup_yr);  // output number of prevalent number of suscpetible population at end of each year, 254 age groups
	NumericMatrix outp_reco(254, spinup_yr);  // output number of prevalent number of recovered population at end of each year, 254 age groups
	
    for (int yr = 0; yr < spinup_yr; ++yr){
		
		NumericVector newinfect(254);      // new infections/cases, 254 age groups

		for (int t = tstep*yr + 1; t <= tstep*(yr+1); ++t) {
			
			// ================================================= 
			// Transmission
			// =================================================  
			tcycle = 1.0 + amp*sin(2.0*M_PI*t/tstep);     //Seasonality in Beta, M_PI = 3.14159265358979323846 
			cyc = tcycle*(in_Comp(_,i_I) + in_Comp(_,i_V1I) + in_Comp(_,i_V2I) + in_Comp(_,i_V3I) + 1.0e-9); 
			
			for (int a = 0; a < 254; ++a){
					
				betta = beta_full(a,_);
				lambda = 1.0 - exp(-sum(betta*cyc)); 
	
				newinfect[a] += lambda*(in_Comp(a,i_S) + in_Comp(a,i_V1S) + in_Comp(a,i_V2S) + in_Comp(a,i_V3S));
								
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
								
				/*// check whether there are negative values
				for (int cp = 0; cp <= i_V3R; ++cp){
					if (trans_Comp(a,cp) < 0.0){	
						Rcout << "time step = " << t << ", age group = " << a+1 << "\n";
						Rcout << "compartment = " << cp+1 << ", value = " << trans_Comp(a,cp) << "\n";
					}
				}*/
			}
			//Rcout << "\ntrnansmission cycle completed\n";
			
			
			// ================================================= 
			// Ageing
			// ================================================= 
			NumericMatrix out_Comp(254, 14);         // output compartments after including transmission and ageing, 254 age groups, 14 states
			
			// ageing process: 3-100 years old
			for (int a = 253; a > 155; --a){
															
				out_Comp(a,i_M)   = trans_Comp(a,i_M)   + age_y*trans_Comp(a-1,i_M)   - age_y*trans_Comp(a,i_M)   - wane*trans_Comp(a,i_M);     
				out_Comp(a,i_S)   = trans_Comp(a,i_S)   + age_y*trans_Comp(a-1,i_S)   - age_y*trans_Comp(a,i_S)   + wane*trans_Comp(a,i_M);     
				out_Comp(a,i_I)   = trans_Comp(a,i_I)   + age_y*trans_Comp(a-1,i_I)   - age_y*trans_Comp(a,i_I)  ; 
				out_Comp(a,i_R)   = trans_Comp(a,i_R)   + age_y*trans_Comp(a-1,i_R)   - age_y*trans_Comp(a,i_R)  ; 
				out_Comp(a,i_V1S) = trans_Comp(a,i_V1S) + age_y*trans_Comp(a-1,i_V1S) - age_y*trans_Comp(a,i_V1S); 
				out_Comp(a,i_V1I) = trans_Comp(a,i_V1I) + age_y*trans_Comp(a-1,i_V1I) - age_y*trans_Comp(a,i_V1I); 
				out_Comp(a,i_V1R) = trans_Comp(a,i_V1R) + age_y*trans_Comp(a-1,i_V1R) - age_y*trans_Comp(a,i_V1R); 
				out_Comp(a,i_V2S) = trans_Comp(a,i_V2S) + age_y*trans_Comp(a-1,i_V2S) - age_y*trans_Comp(a,i_V2S);	
				out_Comp(a,i_V2I) = trans_Comp(a,i_V2I) + age_y*trans_Comp(a-1,i_V2I) - age_y*trans_Comp(a,i_V2I); 
				out_Comp(a,i_V2R) = trans_Comp(a,i_V2R) + age_y*trans_Comp(a-1,i_V2R) - age_y*trans_Comp(a,i_V2R);	
				out_Comp(a,i_V3S) = trans_Comp(a,i_V3S) + age_y*trans_Comp(a-1,i_V3S) - age_y*trans_Comp(a,i_V3S);	
				out_Comp(a,i_V3I) = trans_Comp(a,i_V3I) + age_y*trans_Comp(a-1,i_V3I) - age_y*trans_Comp(a,i_V3I);	
				out_Comp(a,i_V3R) = trans_Comp(a,i_V3R) + age_y*trans_Comp(a-1,i_V3R) - age_y*trans_Comp(a,i_V3R);	
				out_Comp(a,i_V1F) = trans_Comp(a,i_V1F) + age_y*trans_Comp(a-1,i_V1F) - age_y*trans_Comp(a,i_V1F);	
			
				//Rcout << a+1 << " ";
			}
			//Rcout << "\nYearly age groups done\n";
			
			// ageing process: 2 week - 2 years old 
			for (int a = 155; a > 0; --a){
															
				out_Comp(a,i_M)   = trans_Comp(a,i_M)   + age_w*trans_Comp(a-1,i_M)   - age_w*trans_Comp(a,i_M)   - wane*trans_Comp(a,i_M);       
				out_Comp(a,i_S)   = trans_Comp(a,i_S)   + age_w*trans_Comp(a-1,i_S)   - age_w*trans_Comp(a,i_S)   + wane*trans_Comp(a,i_M);       
				out_Comp(a,i_I)   = trans_Comp(a,i_I)   + age_w*trans_Comp(a-1,i_I)   - age_w*trans_Comp(a,i_I)  ; 
				out_Comp(a,i_R)   = trans_Comp(a,i_R)   + age_w*trans_Comp(a-1,i_R)   - age_w*trans_Comp(a,i_R)  ; 
				out_Comp(a,i_V1S) = trans_Comp(a,i_V1S) + age_w*trans_Comp(a-1,i_V1S) - age_w*trans_Comp(a,i_V1S); 
				out_Comp(a,i_V1I) = trans_Comp(a,i_V1I) + age_w*trans_Comp(a-1,i_V1I) - age_w*trans_Comp(a,i_V1I); 
				out_Comp(a,i_V1R) = trans_Comp(a,i_V1R) + age_w*trans_Comp(a-1,i_V1R) - age_w*trans_Comp(a,i_V1R); 
				out_Comp(a,i_V2S) = trans_Comp(a,i_V2S) + age_w*trans_Comp(a-1,i_V2S) - age_w*trans_Comp(a,i_V2S);	
				out_Comp(a,i_V2I) = trans_Comp(a,i_V2I) + age_w*trans_Comp(a-1,i_V2I) - age_w*trans_Comp(a,i_V2I); 
				out_Comp(a,i_V2R) = trans_Comp(a,i_V2R) + age_w*trans_Comp(a-1,i_V2R) - age_w*trans_Comp(a,i_V2R);	
				out_Comp(a,i_V3S) = trans_Comp(a,i_V3S) + age_w*trans_Comp(a-1,i_V3S) - age_w*trans_Comp(a,i_V3S);	
				out_Comp(a,i_V3I) = trans_Comp(a,i_V3I) + age_w*trans_Comp(a-1,i_V3I) - age_w*trans_Comp(a,i_V3I);	
				out_Comp(a,i_V3R) = trans_Comp(a,i_V3R) + age_w*trans_Comp(a-1,i_V3R) - age_w*trans_Comp(a,i_V3R);	
				out_Comp(a,i_V1F) = trans_Comp(a,i_V1F) + age_w*trans_Comp(a-1,i_V1F) - age_w*trans_Comp(a,i_V1F);	
			
				//Rcout << a+1 << " ";
			}  
			//Rcout << "\nWeekly age groups done\n";
			
			// ageing process: 1 week at birth
			pop_fert_SR = 0.0;
			pop_fert_R = 0.0;
			for (int a = 171; a < 189; ++a){  // fertility age: 18-35 years old
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
			
			in_Comp = out_Comp;
		}	
		
		outp_incd(_, yr) = newinfect;
		outp_susc(_, yr) = in_Comp(_,i_S) + in_Comp(_,i_V1S) + in_Comp(_,i_V2S) + in_Comp(_,i_V3S);
		outp_reco(_, yr) = in_Comp(_,i_M) + in_Comp(_,i_R) + in_Comp(_,i_V1R) + in_Comp(_,i_V2R) + in_Comp(_,i_V3R);
    }
	
    outp["cases"] = outp_incd;  
    outp["susc"]  = outp_susc;
	outp["reco"]  = outp_reco;
	outp["out_Comp"] = in_Comp;
	
    return outp;
}