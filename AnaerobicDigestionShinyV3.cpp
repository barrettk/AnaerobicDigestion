$PROB Hydrolysis --> Acidogenesis --> Acetogenesis --> Methanogenesis
  
$PARAM 
  
  // Reactor Operating Parameters (must be in hours)
  scale_Day = 1 // 1 for Day scaling, 0 for hourly.
  Days = 4.17 // Number of dates
  flow = 810 // L/h (810 * 55 wells = 44550)
  mass_flow_per_well = 1726 // lb/h --> round(1866-22.11*348/55) (subtracted sand)
  BATCH = 1 // !1 for CSTR (i.e. 0)
  Temp = 30 // celsius
  vol_per_well = 786 // L --> round(810-3.78541*348/55) (subtracted sand)
  wells = 55
  // H2_perc_head_space = 0.15 // 0.015 // gas mole fraction
  // CO2_perc_head_space = 0.85 // 0.12 // gas mole fraction
  Head_Space_VolRatio = 2 // Fraction of liquid volume
  Pressure_Rest_IN = 101325 // Pressure (atmospheric or vacuum)

  // Component Inflow Rates/Concentrations
  WT_Perc_Water_IN = 0.875
  WT_Perc_Guar_IN = 0.0083
  WT_Perc_PEG_IN = 0.004
  WT_Perc_MeOL_IN = 0.004
  WT_Perc_ISO_IN = 0.004
  Conc_Hydro_Enzyme = 0.0002
  
  MW_PEG_In = 400
  molFracPEG9 = 0.725705329
  molFracPEG8 = 0.224809673
  molFracPEG7 = 0.046574116
  molFracPEG6 = 0.002686968
  molFracPEG5 = 0.000223914

  // PEG_In = 6974 // mol/h
  // HydroLight_IN = 6974 // mol/h

  // Bacetria (vol_per_well/volFromPaper in Liters --> 810/(20/1000) = 40500)
  Bact_ScaleFact_Acido = 0.3 
  Bact_ScaleFact_Aceto = 3
  Bact_ScaleFact_Meth = 0.5
  Bact_ScaleFact_Bact = 0.3

// Kinetic Data Import From Bacteria_kinetics.R
  // Acetogens
  // Vmax
    Intercept_Aceto_vmax = 201.7769
    a_Aceto_vmax = 10.35778
    b_Aceto_vmax = 1.858871
    c_Aceto_vmax = -0.05282864
  // Km
    Intercept_Aceto_km = 274.5769
    a_Aceto_km = -34.11182
    b_Aceto_km = 3.224341
    c_Aceto_km = -0.06061024
  // Methanogens
  // Vmax
    Intercept_Meth_vmax = -60.84599
    a_Meth_vmax = 44.76423
    b_Meth_vmax = -0.5325213
    c_Meth_vmax = -0.008983561
  // Km
    Intercept_Meth_km = 31.89081
    a_Meth_km = 3.502035
    b_Meth_km = 0.1449059
    c_Meth_km = -0.004059418

$FIXED
  
  // Henry's Constants mol/L-atm at 25C
  // https://www.atmos-chem-phys.net/15/4399/2015/acp-15-4399-2015.pdf
  Henry_0_H2 = 7.8e-4
  Henry_0_CO2 = 3.4e-2
  Henry_0_O2 = 1.3e-3
  Henry_0_N2 = 6.1e-4
  // Henry_0_CH4 reported in units of 1.4e-4 mol/m^3-Pa
  Henry_0_CH4 = 1.4e-2
  // Enthalpy (Delta) of Dissolution Divided by Gas Constant (-Hsol/R) - For calculating Henry's constant at different temperatures
  // (Units: K (temp))
  Hsol_H2 = 500 
  Hsol_CO2 = 2400
  Hsol_O2 = 1700 
  Hsol_N2 = 1300 
  Hsol_CH4 = 1600

  // Molecular Weights (g/mol)
  MW_guar = 1250000 // 800,000 - 8,000,000
  MW_glucose = 180
  MW_ethanol = 46
  MW_propA = 74
  MW_acetate = 59
  MW_H2 = 2
  MW_CO2 = 44
  MW_CH4 = 16
  MW_H2O = 18
  MW_methanol = 32
  MW_isopropanol = 60
  //MW_PEG400 = 400 // 380-420

  // Cleavable Bonds per mole guar gum: Guar Gum MW/Glucose MW (1,200,000/180 = 6666.667 ~ 6500)
  n_bonds = 6600 // 1200000/266-1 in paper
  // Density (g/mL)
  Density_guar = 0.9 
  Density_Water = 1 
  Density_Peg400 = 1.128 
  Density_Methanol = 0.791
  Density_Isopropanol = 0.785


$FIXED @annotated
  
  pi             : 3.1415926 : Pi (π)
  R_coeff        : 8314      : Gas Constant, Units of (L-Pa)/(K-mol)  
  cell_density   : 1         : Unis of g/mL
  N_cell_Meth_0  : 5.5e+07   : Methanogen Cells per mL (Typically 10^7 to 10^8 cells per ml)
  N_cell_Aceto_0 : 5.5e+07   : Acetogen Cells per mL (Typically 10^7 to 10^8 cells per ml)
  N_cell_Acido_0 : 5.5e+07   : Acetogen Cells per mL (Typically 10^7 to 10^8 cells per ml)
  N_cell_Bact_0  : 5.5e+07   : Bacetoides per mL (Assumed the same as other bacteria strains)
  
$PLUGIN Rcpp mrgx
  

$GLOBAL
  
  namespace {
  Rcpp::Environment env;
  }

  #include <Bacteria_param_functions.h>
  #define LOGIT(a) (log(a/(1-a)))
  #define ALOGIT(a) (exp(a)/(1+exp(a)))

$CMT @annotated
  
  GUAR                : Guar Gum (Substrate - moles)
  GLUCOSE             : Glucose (Substrate - moles)
  ETHANOL             : Ethanol (Substrate - moles)
  CO2_TOT             : CO2 (Gas and Liquid)
  H2_TOT              : H2 (Gas and Liquid)
  CH4_TOT             : CH4 (Gas and Liquid)
  PropAcid            : Propanoic Acid (Substrate - moles)
  ACETATE             : Acetic Acid (Substrate - moles)
  H2O                 : Water (moles)
  
  PEG_9               : PEG-9 (~PEG 414 - moles)
  PEG_8               : PEG-8 (~PEG 370 - moles)
  PEG_7               : PEG-7 (~PEG 326 - moles)
  PEG_6               : PEG-6 (~PEG 282 - moles)
  PEG_5               : PEG-5 (~PEG 238 - moles)
  PEG_4               : PEG-4 (~PEG 194 - moles)
  PEG_3               : PEG-3 (~PEG 150 - moles)
  DEG                 : PEG-2 (~PEG 106 - moles)
  EG                  : Ethylene Glycol (moles)
  Acetaldehyde        : Acetaldehyde (moles)
  
  METHANOL            : Methanol (moles)
  ISOPROPANOL         : Isopropanol (moles)
  
  ACIDOGEN            : Acidogen Growth (grams)
  ACETOGEN            : Acetogen Growth (grams)
  METHANOGEN          : Methanogen Growth (grams)
  BACTEROID           : Bacteroid Growth (grams)
  Temp2               : Temperature (C)
 
$OMEGA @name Bact_kinetics @annotated
  
  Eta_Vmax_Guar   : 0.06 : ETA on Vmax of Guar Gum
  Eta_Km_Guar     : 0.06 : ETA on Km of Guar Gum
  
  Eta_Umax_Gluc   : 0.06 : ETA on Umax of Glucose
  Eta_Ks_Gluc     : 0.06 : ETA on Ks of Glucose
  
  Eta_Umax_Eth    : 0.06 : ETA on Umax of Ethanol
  Eta_Ks_Eth      : 0.06 : ETA on Ks of Ethanol
  
  Eta_Umax_propA  : 0.06 : ETA on Umax of Propanoic Acid
  Eta_Ks_propA    : 0.06 : ETA on Ks of Propanoic Acid
  
  Eta_Umax_Acet   : 0.06 : ETA on Umax of Acetate
  Eta_Ks_Acet     : 0.06 : ETA on Ks of Acetate
  
  Eta_Kcat_PEG    : 0.05 : ETA on Kcat of PEG structures
  Eta_Km_PEG      : 0.05 : ETA on Km of PEG structures
  
  Eta_Kcat_MeOL   : 0.05 : ETA on Kcat of Methanol
  Eta_Kcat_ISO    : 0.05 : ETA on Kcat of Isopropanol
 
$OMEGA @name Bact_yields @annotated
  
  Eta_Y_Gluc   : 0.01 : ETA on Yield of Glucose
  
  Eta_Y_Eth    : 0.01 : ETA on Yield of Ethanol
  
  Eta_Y_propA  : 0.01 : ETA on Yield of Propanoic Acid
  
  Eta_Y_Acet   : 0.01 : ETA on Yield of Acetate
  
  Eta_Y_PEG    : 0.01 : ETA on Yield of PEG structures
  
  Eta_Y_MeOL   : 0.01 : ETA on Yield of Isopropanol
  
  Eta_Y_ISO    : 0.01 : ETA on Yield of Isopropanol
  
$MAIN
  
  // Sizing Parameters (i.e. head space volume & pressure) 
  if(scale_Day==0){
    double Day_Mult = 1;
  }else{
    Day_Mult = 24 * Days; //only running for 12 hours a day (not 24)
  }
  
  double vol = vol_per_well*wells*Day_Mult; // L
  double Head_Space_Vol = Head_Space_VolRatio * vol; //vol_per_well * Day_Mult; // L
  // PEG MW Calcs
  // C(2n)H(4n+2)O(n+1) + (n-1)H2O → nCH3COOH + nH2
  // int PegStoich_water = n_Peg-1;
  // int PegStoich_acetate = n_Peg;
  // int PegStoich_H2 = n_Peg; 
  
  // Influent moles
  double guar_IN = (WT_Perc_Guar_IN*mass_flow_per_well*453.6*wells*Day_Mult)/MW_guar; // mol/h or mole
  double water_IN = (WT_Perc_Water_IN*mass_flow_per_well*453.6*wells*Day_Mult)/MW_H2O; //mol/h or mole
  double PEG_IN = (WT_Perc_PEG_IN*mass_flow_per_well*453.6*wells*Day_Mult)/MW_PEG_In; //mol/h or mole
  double MeOL_IN = (WT_Perc_MeOL_IN*mass_flow_per_well*453.6*wells*Day_Mult)/MW_methanol; //mol/h or mole
  double ISO_IN = (WT_Perc_ISO_IN*mass_flow_per_well*453.6*wells*Day_Mult)/MW_isopropanol; //mol/h or mole
  // Volumes (L)
  double vol_guar_IN = (guar_IN*MW_guar/Density_guar)/1000;
  double vol_water_IN = (water_IN*MW_H2O/Density_Water)/1000;
  double vol_peg_IN = (PEG_IN*MW_PEG_In/Density_Peg400)/1000;
  double vol_methanol_IN = (MeOL_IN*MW_methanol/Density_Methanol)/1000;
  double vol_ISO_IN = (ISO_IN*MW_isopropanol/Density_Isopropanol)/1000;
  double vol_Rest = vol - vol_guar_IN - vol_water_IN - vol_peg_IN - vol_methanol_IN - vol_ISO_IN;
  
  
// Hydrolysis of Guar Gum
  
    //Can be sped up with temp --> rate doubled from 38 to 55C
    //See https://www.researchgate.net/publication/51811480_Relative_kinetics_of_anaerobic_digestion_under_thermophilic_and_mesophilic_conditions
  // Guar Gum Kinetic Parameters
  double Vmax_Guar_Hydro = (7.8e-10 * 1000 * 3600)*exp(Eta_Vmax_Guar); // mol/mL-s --> mol/L-h
  double Km_Guar_Hydro = (0.6 / 1000)*exp(Eta_Km_Guar); // mM --> M
  // Hydrolysis Enzyme
  double Conc_hydro_Enz_IN = (Conc_Hydro_Enzyme / 1e+06) * 1000 * 60 * vol_guar_IN; // units/mL (1 unit = 1 umol/min) --> units/L (1 unit = 1 mol/h)
        //reported as 1e-10 M
  
// Degradation of PEG by Bacteroides

  //MW's
  int MW_PEG_9 = MW_H2O + 44*9;
  int MW_PEG_8 = MW_H2O + 44*8;
  int MW_PEG_7 = MW_H2O + 44*7;
  int MW_PEG_6 = MW_H2O + 44*6;
  int MW_PEG_5 = MW_H2O + 44*5;
  int MW_PEG_4 = MW_H2O + 44*4;
  int MW_PEG_3 = MW_H2O + 44*3;
  int MW_DEG = MW_H2O + 44*2;
  int MW_EG = MW_H2O + 44*1;
  int MW_AcetHyde = 44;
  // Kcat's
  double kcat_PEG_9 = 0.264/140 * exp(Eta_Kcat_PEG);
  double kcat_DEG = 0.035/30;
  // Interpolate for others
  double Intercept = (kcat_PEG_9-((kcat_PEG_9-kcat_DEG)/(9-2))*9);
  double kcat_PEG_8 = ((kcat_PEG_9-kcat_DEG)/(9-2))*8 + Intercept;
  double kcat_PEG_7 = ((kcat_PEG_9-kcat_DEG)/(9-2))*7 + Intercept;
  double kcat_PEG_6 = ((kcat_PEG_9-kcat_DEG)/(9-2))*6 + Intercept;
  double kcat_PEG_5 = ((kcat_PEG_9-kcat_DEG)/(9-2))*5 + Intercept;
  double kcat_PEG_4 = ((kcat_PEG_9-kcat_DEG)/(9-2))*4 + Intercept;
  double kcat_PEG_3 = ((kcat_PEG_9-kcat_DEG)/(9-2))*3 + Intercept;
  
  double kcat_EG = abs(0.023 * exp(Eta_Kcat_PEG));
  double kcat_AcetHyde = abs(0.055 * exp(Eta_Kcat_PEG)); //come back later
  // Km's
  double km_PEG_22 = abs(0.0016 * exp(Eta_Km_PEG));
  double km_PEG_9 = abs(0.0007 * exp(Eta_Km_PEG));
  double km_DEG = abs(0.0065);
  double km_EG = abs(0.005 * exp(Eta_Km_PEG));
  double km_AcetHyde = abs(0.002 * exp(Eta_Km_PEG)); //come back later
  // Yields
  double Y_PEG_22 = abs(0.00365 * exp(Eta_Y_PEG));
  double Y_PEG_9 = abs(0.00556 * exp(Eta_Y_PEG));
  double Y_DEG = abs(0.03962 * exp(Eta_Y_PEG));
  double Y_EG = abs(0.08065 * exp(Eta_Y_PEG));
  double Y_AcetHyde = abs(0.04773 * exp(Eta_Y_PEG)); //verified
  
// Acidogenesis 
  
  // Glucose Kinetic Parameters
  //Used: https://ac-els-cdn-com.ezproxy2.library.drexel.edu/S0961953411004016/1-s2.0-S0961953411004016-main.pdf?_tid=1bde782c-a1a9-475e-9d65-95b54d158e66&acdnat=1550172665_09381fa91b974dd6780732ebec85643f
  double Umax_Glucose_Acido = abs(0.163*exp(Eta_Umax_Gluc)); //0.19; //0.0564; (1/h)
  double Ks_Glucose_Acido =  abs(((.042/1000)*MW_glucose)*exp(Eta_Ks_Gluc)); //0.0801; //(0.0801)*MW_glucose; // 0.005 Divided by MW Glucose
  double Y_Acido_Glucose = abs(((0.027*1000)/MW_glucose)*exp(Eta_Y_Gluc)); //16/MW_glucose; //0.86; // g Biomass/g Substrate
  double Y_Ace_Gluc = 0.789; //0.03
  double Y_Eth_Gluc = 0.249; //0.38
  double Y_PA_Gluc = 0.545; 
  double Y_H2_Gluc = 0.759; 
  double Y_CO2_Gluc = 0.816;
  double num_AcidogenReactions = 3; // 3 reactions utilize glucose as substrate
  
// Acetogenesis
  
    //Strains
    //https://onlinelibrary-wiley-com.ezproxy2.library.drexel.edu/doi/pdf/10.1002/0471468967
    
  // Method 1: Monod Equation
  // For Reactions:	
  // 1a. CH3CH2OH + 2H2O → CH3COO- + 3H2 + H+
  // 2a. 	CH3CH2OO-  + 3H2O → CH3COO- + H+ + 3H2 + HCO3- ; 2b. 	H+ + HCO3- → H2CO3 →  CO2 + H2O
  // 1. ---> Model 1a. ---> CH3CH2OH + 2H2O → CH3COOH + 3H2
  // 2. ---> Model 2a. + 2b. ---> CH3CH2OOH  + 3H2O → CH3COOH + 3H2 + CO2 + H2O 
  // 3. ---> 4CH3OH + CO2 → 3CH3COOH + 2H2O
  
  // Ethanol Kinetic Parameters
  // Used: https://www-jstor-org.ezproxy2.library.drexel.edu/stable/25164672?pq-origsite=summon&seq=6#metadata_info_tab_contents
  // Alternative: https://pdfs.semanticscholar.org/eb8a/50666189f5278c0cddec611f16dcb780baef.pdf
  double Umax_eth_Aceto = abs(0.8910*exp(Eta_Umax_Eth)); // (1/h) also 0.643 or 0.3052 or 7.2/24;
  double Ks_eth_Aceto = abs((0.2030/1000)*exp(Eta_Ks_Eth));  // also (7.8/1000)/MW_ethanol or 0.001; 
  double Y_Aceto_eth = abs((0.08911)*exp(Eta_Y_Eth)); // also 0.4/0.7 or 0.333; // g Biomass/g Substrate
  // Propanoic Acid Kinetic Parameters
  // Used: https://search-proquest-com.ezproxy2.library.drexel.edu/docview/1943342519?pq-origsite=summon
    // Bacteria Name: Desulfobulbus propionicus, grow in/use sulfate
  // Other: For Ks and y:https://onlinelibrary-wiley-com.ezproxy2.library.drexel.edu/doi/epdf/10.1021/bp00010a012
  double Umax_propA_Aceto = abs(0.1122*exp(Eta_Umax_propA));  //(1/h) also 1.66/24 or 7.2/24 or 0.4/24 or  0.274/24;
  double Ks_propA_Aceto = abs(((1.004e-05)/1000)*exp(Eta_Ks_propA)); // also 0.25/MW_propA or 0.001; 
  double Y_Aceto_propA = abs((0.03034)*exp(Eta_Y_propA)); // also 0.333; // g Biomass/g Substrate
  // Stoichiometric parameters of products in gas phase
  double Y_H2_Acido = 0.03;
  double Y_CO2_Acido = 0.03;
  double Y_H2_Ace = 0.03;
  double Y_CO2_Ace = 0.03;
  
  // Method 2: Michaelis-Menten Kinetics
  // Uses Kinetic Data Import From Bacteria_kinetics.R
  // For Reaction: 2CO2 + 4H2 → CH3COOH + 2H2O ---> Recast ---> (2/4)CO2 + H2 → (1/4)CH3COOH + (2/4)H2O
  
  // H2 Threshold
  double H2_Thresh_Aceto = 55; //Hydrogen Threshold (Pa)
  double H2_Thresh_Aceto_LIQ = H2_Thresh_Aceto * Henry_H2;
  // Decay
  // https://www-jstor-org.ezproxy2.library.drexel.edu/stable/pdf/25164672.pdf?refreqid=excelsior%3Ae05924d663e99ccb166f9aa85942883c
  // kd =0.024 1/d --> 0.001 1/h --> see other paper too
  // subtract -kd*X from bacteria to implement
  
  
// Methanogenesis

    //U sed: https://www-jstor-org.ezproxy2.library.drexel.edu/stable/25164672?pq-origsite=summon&seq=6#metadata_info_tab_contents
    // Also: methyltrophic bacteria for methanol: https://www.springer.com/cda/content/document/cda_downloaddocument/9789811041297-c2.pdf?SGWID=0-0-45-1602425-p180697263
    // Methanosarcina bakerii for Acetate, carbon dioxide, hydrogen, methanol
    // Methanol paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC242456/
  // Method 1: Monod Equation
  // For Reaction: CH3COOH → CH4 + CO2
  // Acetate Kinetic Parameters
  double Umax_acet_Meth = abs(0.245*exp(Eta_Umax_Acet)); //(1/h) also 0.44/24 or 0.138;
  double Ks_acet_Meth = abs((0.3/MW_acetate)*exp(Eta_Ks_Acet));  //also (58/1000)/MW_acetate or 0.3/MW_acetate or (20.51e-3) or 0.0208 or (0.13*10^-3)*MW_acetate 
  double Y_Meth_Acet = abs((0.06)*exp(Eta_Y_Acet)); // g Biomass/g Substrate also 3.295/MW_acetate 
  
  // Methanol Kinetic Parameters
  double Kcat_methanol = abs(0.002*exp(Eta_Kcat_MeOL)); //(1/h)
  double Ks_methanol = 20/1000; 
  double Y_methanol_X = abs(((16/MW_methanol)/2)*exp(Eta_Y_MeOL)); // g Biomass/g Substrate 
  double Y_H2_methanol_X = (16/MW_methanol)/2; // g Biomass/g Substrate 
  double Y_CH4_methanol = 0.85/2; //mol methane/mol methanol
  double Y_H2O_methanol = 0.85/2; //mol methane/mol methanol
  // Isopropanol Kinetic Parameters
  double Kcat_isopropanol = abs(0.0009*exp(Eta_Kcat_ISO)); //(1/h)
  double Ks_isopropanol = 40/1000; 
  double Y_isopropanol_X = abs((1.6/MW_isopropanol)*exp(Eta_Y_ISO)); // g Biomass/g Substrate 
  double Y_H2_isopropanol_X = (1.6/MW_isopropanol)/5; // g Biomass/g Substrate 
  double Y_CH4_isopropanol = 0.85/2; //mol methane/mol methanol
  double Y_H2O_isopropanol = 0.85/2; //mol methane/mol methanol
  // Method 2: Michaelis-Menten
  // Uses Kinetic Data Import From Bacteria_kinetics.R
  // For Reaction: CO2 + 4H2 → CH4 + 2H2O ---> Recast ---> (1/4)CO2 + H2 → (1/4)CH4 + (2/4)H2O

  // H2 Threshold
  double H2_Thresh_MSB = 1; //Hydrogen Threshold (Pa)Table
  double H2_Thresh_MSB_LIQ = H2_Thresh_MSB * Henry_H2;
  //Decay
  //https://www-jstor-org.ezproxy2.library.drexel.edu/stable/pdf/25164672.pdf?refreqid=excelsior%3Ae05924d663e99ccb166f9aa85942883c
  //kd =0.024 1/d --> 0.001 1/h
  //subtract -kd*X from bacteria to implement
  
// Biomass Calculations for all Strains -----------
  
  // Bacteria Count in paper 
  double vtot_bact_ini = 20; //Total volume in mL's - important for scaling vmax, which was reported per 20 mL's
  double V_Bact_ini = vtot_bact_ini / 1e+03; //Bacterial Suspension (20 mL = 0.02 liters)
  // Methanogenic Strain MSB Size Parameters
  double cell_width = 0.7 * 1e-06; //um
  double cell_length = 28 * 1e-06; //um
  // Cell Volume (mL) - Spirillum Shaped (curved rods) for Methanogenic Strain MSB --> V = (pi*r^2)*h
  double V_cell = (pi*pow((cell_width/2),2))*cell_length; // meters^3/cell
  
  double x_Meth = N_cell_Meth_0*1e+06*V_cell*(cell_density*1e+06); // Biomass concentration (gram/meters^3)
  double x_adj_Meth = (x_Meth)/(1000); // Biomass concentration (g/L)
  
  double x_Aceto = N_cell_Aceto_0*1e+06*V_cell*(cell_density*1e+06); // Biomass concentration (gram/meters^3)
  double x_adj_Aceto = (x_Aceto)/(1000); // Biomass concentration (g/L)
  
  double x_Acido = N_cell_Acido_0*1e+06*V_cell*(cell_density*1e+06); // Biomass concentration (gram/meters^3)
  double x_adj_Acido = (x_Aceto)/(1000); // Biomass concentration (g/L)
  x_adj_Acido = 0.3;
  x_adj_Aceto = 0.3;
  x_adj_Meth = 0.3;
  
  // Reactor Control
  if(BATCH==0){ // Not implemented
    double Flow = flow;
    // Concentrations (inflow in moles/h)
    double Conc_guar_IN = (guar_IN*n_bonds/Flow) ; // mol/L
    double Conc_water_IN = (water_IN/Flow);
  }else if(BATCH==1){
    Flow = 0; 
    // Concentrations (inflow in moles)
    H2O_0 = water_IN;
    GUAR_0 = guar_IN*n_bonds; // mol
    PEG_9_0 = PEG_IN*molFracPEG9; // mol
    PEG_8_0 = PEG_IN*molFracPEG8; // mol
    PEG_7_0 = PEG_IN*molFracPEG7; // mol
    PEG_6_0 = PEG_IN*molFracPEG6; // mol
    PEG_5_0 = PEG_IN*molFracPEG5; // mol
    METHANOL_0 = MeOL_IN; // mol
    ISOPROPANOL_0 = ISO_IN; // mol
  }
  
  // Temperature Control
  Temp2_0 = Temp; //works
  
  // Bacteria Control (g/L)
  ACIDOGEN_0 = Bact_ScaleFact_Acido * vol; 
  ACETOGEN_0 = Bact_ScaleFact_Aceto * vol;
  METHANOGEN_0 = Bact_ScaleFact_Meth * vol;
  BACTEROID_0 = Bact_ScaleFact_Bact * vol;
  
$ODE
    
  // Henry's Law for Gases
  // Calculate Henry's Constant at Temperature (in K) according to van 't Hoff equation (mol/L-Pa) (Good for small changes in temp)
  double TempK = Temp2 + 273.15;
  double Temp0K = 25 + 273.15;
  double Henry_H2 = (Henry_0_H2 * exp(Hsol_H2*((1/TempK)-1/Temp0K)))/101325;
  double Henry_CO2 = (Henry_0_CO2 * exp(Hsol_CO2*((1/TempK)-1/Temp0K)))/101325;
  double Henry_O2 = (Henry_0_O2 * exp(Hsol_O2*((1/TempK)-1/Temp0K)))/101325;
  double Henry_N2 = (Henry_0_N2 * exp(Hsol_N2*((1/TempK)-1/Temp0K)))/101325;
  double Henry_CH4 = (Henry_0_CH4 * exp(Hsol_CH4*((1/TempK)-1/Temp0K)))/101325;
  
  // Temperature Rate Control
  dxdt_Temp2 = 0;
  
  // Concentrations (mol/L)
    // Guar Gum Process
  double Conc_GUAR =  GUAR/vol;
  double Conc_GLUCOSE = GLUCOSE/vol;
  double Conc_ETHANOL = ETHANOL/vol;
  double Conc_CO2 = CO2_LIQ/vol;
  double Conc_H2 = H2_LIQ/vol;
  double Conc_PropAcid = PropAcid/vol;
  double Conc_ACETATE = ACETATE/vol;
  double Conc_CH4 = CH4_LIQ/vol;
  double Conc_H2O = H2O/vol;
  // PEG Process
  double Conc_PEG9 = PEG_9/vol;
  double Conc_PEG8 = PEG_8/vol;
  double Conc_PEG7 = PEG_7/vol;
  double Conc_PEG6 = PEG_6/vol;
  double Conc_PEG5 = PEG_5/vol;
  double Conc_PEG4 = PEG_4/vol;
  double Conc_PEG3 = PEG_3/vol;
  double Conc_DEG = DEG/vol;
  double Conc_EG = EG/vol;
  double Conc_AcetHyde = Acetaldehyde/vol;
  // Alcohols
  double Conc_METHANOL = METHANOL/vol;
  double Conc_ISOPROPANOL = ISOPROPANOL/vol;
  // Bacteria
  double Conc_ACIDOGEN = ACIDOGEN/vol;
  double Conc_ACETOGEN = ACETOGEN/vol;
  double Conc_METHANOGEN = METHANOGEN/vol;
  double Conc_BACTEROID = BACTEROID/vol;
  double C_TOT = Conc_GUAR/n_bonds + Conc_GLUCOSE + Conc_ETHANOL + Conc_PropAcid + Conc_CO2 + Conc_H2 +
    Conc_ACETATE + Conc_CH4 + Conc_H2O + Conc_PEG3 + Conc_PEG4 + Conc_PEG5 + Conc_PEG6 + Conc_PEG7 +
    Conc_PEG8 + Conc_PEG9 + Conc_DEG + Conc_EG + Conc_AcetHyde + Conc_METHANOL + Conc_ISOPROPANOL;
  
  // Rates
  
  // Hydrolysis of Guar Gum (Depolymerization)
  double kcat_G = Vmax_Guar_Hydro / Conc_hydro_Enz_IN;
  double vol_guar = ((Conc_GUAR/n_bonds)*vol*MW_guar/Density_guar)/1000; // L
  double H2O_Vol_L = (H2O*MW_H2O/1)/1000; // L
  double V_TOT = vol_Rest + vol_guar + H2O_Vol_L; //L
  double Head_Space_VTOT = Head_Space_VolRatio*V_TOT; //(V_TOT/wells);
  // Hydrolysis Enzyme
  double Conc_hydro_Enz = (Conc_Hydro_Enzyme / 1e+06) * 1000 * 60 * vol_guar; // units/mL (1 unit = 1 umol/min) --> units/L (1 unit = 1 mol/h) --> *volGuar
  // Rate Law
  double Hydro_Rate = (kcat_G * Conc_hydro_Enz) / (Km_Guar_Hydro + Conc_GUAR); // molar concentration of cleavable bonds in the system
  
  // Degradation of PEG Structures (Kcat and Km change at PEG 1000 and PEG 400)
  double Rate_PEG9 = (kcat_PEG_9 * Conc_PEG9)/(km_PEG_9 + Conc_PEG9);
  double Rate_PEG8 = (kcat_PEG_8 * Conc_PEG8)/(km_PEG_9 + Conc_PEG8);
  double Rate_PEG7 = (kcat_PEG_7 * Conc_PEG7)/(km_PEG_9 + Conc_PEG7);
  double Rate_PEG6 = (kcat_PEG_6 * Conc_PEG6)/(km_PEG_9 + Conc_PEG6);
  double Rate_PEG5 = (kcat_PEG_5 * Conc_PEG5)/(km_PEG_9 + Conc_PEG5);
  double Rate_PEG4 = (kcat_PEG_4 * Conc_PEG4)/(km_PEG_9 + Conc_PEG4);
  double Rate_PEG3 = (kcat_PEG_3 * Conc_PEG3)/(km_PEG_9 + Conc_PEG3);
  double Rate_DEG = (kcat_DEG * Conc_DEG)/(km_DEG + Conc_DEG);
  double Rate_EG = (kcat_EG * Conc_EG)/(km_EG + Conc_EG);
  double Rate_AcetHyde = (kcat_AcetHyde * Conc_AcetHyde)/(km_AcetHyde + Conc_AcetHyde);
  
  // Acidogenesis (Growth rates per substrate)
  double U_Acido_Glucose = (Umax_Glucose_Acido * Conc_GLUCOSE) / (Ks_Glucose_Acido + Conc_GLUCOSE);
  
  // Acetogenesis (Growth rates per substrate)
  // Method 1
  double U_Aceto_Eth = (Umax_eth_Aceto * Conc_ETHANOL) / (Ks_eth_Aceto + Conc_ETHANOL);
  double U_Aceto_PropA = (Umax_propA_Aceto * Conc_PropAcid) / (Ks_propA_Aceto + Conc_PropAcid);
  // Method 2
  double vmax_Aceto = param_fn::Aceto_vmax(Temp2,Intercept_Aceto_vmax,a_Aceto_vmax,b_Aceto_vmax,c_Aceto_vmax); // nmol/h
  double km_Aceto = param_fn::Aceto_km(Temp2,Intercept_Aceto_km,a_Aceto_km,b_Aceto_km,c_Aceto_km); // Pa
  double km_Aceto_LIQ = km_Aceto * Henry_H2; // mol/L
  // Convert to correct units
  double kcat_Aceto = (vmax_Aceto / (1e+09)) / x_adj_Aceto;// * 1000 * 4; // mol/h
  if(H2_Pressure >= H2_Thresh_Aceto && CO2_GAS > 0){
    // Normalize rate based on biomass concentration from paper (x_adj_Aceto)
    double Aceto_H2_gas = ((kcat_Aceto*(H2_Pressure - H2_Thresh_Aceto))/(km_Aceto + H2_Pressure - H2_Thresh_Aceto));// * (V_TOT/0.02); // (mol-gas/h)/(g/L)
    double Aceto_H2 = ((kcat_Aceto*(Conc_H2 - H2_Thresh_Aceto_LIQ))/(km_Aceto_LIQ + Conc_H2 - H2_Thresh_Aceto_LIQ));// * (V_TOT/0.02); // (mol-liquid/h)/(g/L)
  }else{
    Aceto_H2_gas = 0;
    //Rate_Aceto_H2_Gas = 0;
    Aceto_H2 = 0;
  }
  
  
  // Methanogenesis (Growth rates per substrate)
  // Method 1
  // Acetate
  double U_Meth_Acet = (Umax_acet_Meth * Conc_ACETATE) / (Ks_acet_Meth + Conc_ACETATE);
  // Alcholols
  double Rate_METHANOL = Kcat_methanol*(Conc_METHANOL/(METHANOL_0/vol)) * (Conc_METHANOGEN/0.5) * (((H2_TOT/4)/(H2_TOT_0+1))/1250000)*25; //(Kcat_methanol * (Conc_METHANOL/(METHANOL_0/vol))) / ((Ks_methanol/(METHANOL_0/vol)) + (Conc_METHANOL/(METHANOL_0/vol)));
  double Rate_ISOPROPANOL = Kcat_isopropanol*(Conc_ISOPROPANOL/(ISOPROPANOL_0/vol)) * (Conc_METHANOGEN/0.5) * (((H2_TOT/4)/(H2_TOT_0+1))/1250000)*25;
  // Method 2
  double vmax_Meth = param_fn::Meth_vmax(Temp2,Intercept_Meth_vmax,a_Meth_vmax,b_Meth_vmax,c_Meth_vmax); // nmol/h
  double km_Meth = param_fn::Meth_km(Temp2,Intercept_Meth_km,a_Meth_km,b_Meth_km,c_Meth_km); // Pa
  double km_Meth_LIQ = km_Meth * Henry_H2; // mol/L
  // Convert to correct units (multiplied by 1000*4 to compare to other papers)
  //0.5/16/24 = 0.001302, vmax_Meth2 = 0.0019 initially
  double kcat_Meth = (vmax_Meth / (1e+09))/ x_adj_Meth;// * 1000 * 4; // mol/h
  if(H2_Pressure >= H2_Thresh_MSB && CO2_GAS > 0){
    // Normalize rate based on biomass concentration from paper (x_adj_Meth)
    double Meth_H2_gas = (((kcat_Meth*(H2_Pressure-H2_Thresh_MSB))/(km_Meth + H2_Pressure - H2_Thresh_MSB)));// * (V_TOT/0.02); // (mol-gas/h)/(g/L)
    double Meth_H2 = (((kcat_Meth*(Conc_H2-H2_Thresh_MSB_LIQ))/(km_Meth_LIQ + Conc_H2 - H2_Thresh_MSB_LIQ)));// * (V_TOT/0.02); // (mol-liquid/h)/(g/L)
  }else{
    Meth_H2_gas = 0;
    //Rate_Meth_H2_Gas = 0;
    Meth_H2 = 0;
  }
  
  
  // Notes on Hydrolysis
  // Hydrolysis decreases guar gum's viscosity and increases solubility as it degrades (Mark–Houwink equation)
  // Three kinds of bonds in guar are susceptible to enzymatic hydrolysis: the endo- and exo-ß-1,4 linkages
  // on the D-mannose backbone and the å-1,6 linkage between the mannose unit and the galactose side-chain. 
  // The enzymes that cleave these bonds are respectively endo- and exo-ß-mannanase and œ-galactosidase.
  // å-galactosidase for natural guar - ORGANISMS produce it (ligand) - probably wont use that though
  // ß-mannanase has a molecular weight of 45,000 --> By assuming the enzyme molecule is a sphere in solution and
  // its density is 0.7 g/mL, which is a standard value for density of globular protein in solution, we can estimate
  // the size of the enzyme molecules --> membrane sizing to keep it in reactor (same with bacteria)
  
  // ODE's
  
  // Hydrolysis of Guar Gum (mol cleavable bonds/h) 
  dxdt_GUAR = Flow*(Conc_guar_IN*n_bonds - Conc_GUAR) - Hydro_Rate*Conc_GUAR*vol;
  
  dxdt_GLUCOSE = -Flow*Conc_GLUCOSE + (Hydro_Rate*Conc_GUAR - (U_Acido_Glucose*Conc_ACIDOGEN/(MW_glucose*Y_Acido_Glucose)) *
    num_AcidogenReactions) * vol;
  
  // Degradation of PEG Structures (mol/h) 
  dxdt_PEG_9 = (-Rate_PEG9*Conc_BACTEROID/(MW_PEG_9*Y_PEG_9)) * vol;
  dxdt_PEG_8 = (Rate_PEG9*Conc_BACTEROID/(MW_PEG_9*Y_PEG_9) - Rate_PEG8*Conc_BACTEROID/(MW_PEG_8*Y_PEG_9)) * vol;
  dxdt_PEG_7 = (Rate_PEG8*Conc_BACTEROID/(MW_PEG_8*Y_PEG_9) - Rate_PEG7*Conc_BACTEROID/(MW_PEG_7*Y_PEG_9)) * vol;
  dxdt_PEG_6 = (Rate_PEG7*Conc_BACTEROID/(MW_PEG_7*Y_PEG_9) - Rate_PEG6*Conc_BACTEROID/(MW_PEG_6*Y_PEG_9)) * vol;
  dxdt_PEG_5 = (Rate_PEG6*Conc_BACTEROID/(MW_PEG_6*Y_PEG_9) - Rate_PEG5*Conc_BACTEROID/(MW_PEG_5*Y_PEG_9)) * vol;
  dxdt_PEG_4 = (Rate_PEG5*Conc_BACTEROID/(MW_PEG_5*Y_PEG_9) - Rate_PEG4*Conc_BACTEROID/(MW_PEG_4*Y_PEG_9)) * vol;
  dxdt_PEG_3 = (Rate_PEG4*Conc_BACTEROID/(MW_PEG_4*Y_PEG_9) - Rate_PEG3*Conc_BACTEROID/(MW_PEG_3*Y_PEG_9)) * vol;
  // Acetogenesis of DEG and EG (mol/h) 
  dxdt_DEG = (Rate_PEG3*Conc_BACTEROID/(MW_PEG_3*Y_PEG_9) - Rate_DEG*Conc_ACETOGEN/(MW_DEG*Y_DEG)) * vol;
  dxdt_EG = (Rate_AcetHyde*Conc_BACTEROID/(MW_AcetHyde*Y_AcetHyde) - Rate_EG*Conc_ACETOGEN/(MW_EG*Y_EG)) * vol;
  // Conversion of Acetaldehyde (mol/h) 
  dxdt_Acetaldehyde = (Rate_PEG9*Conc_BACTEROID/(MW_PEG_9*Y_PEG_9) + Rate_PEG8*Conc_BACTEROID/(MW_PEG_8*Y_PEG_9) +
    Rate_PEG7*Conc_BACTEROID/(MW_PEG_7*Y_PEG_9) + Rate_PEG6*Conc_BACTEROID/(MW_PEG_6*Y_PEG_9) +
    Rate_PEG5*Conc_BACTEROID/(MW_PEG_5*Y_PEG_9) + Rate_PEG4*Conc_BACTEROID/(MW_PEG_4*Y_PEG_9) +
    Rate_PEG3*Conc_BACTEROID/(MW_PEG_3*Y_PEG_9) - Rate_AcetHyde*Conc_BACTEROID/(MW_AcetHyde*Y_AcetHyde)) * vol;
  // Alcohols (mol/h)
  dxdt_METHANOL = (-Rate_METHANOL)*vol;
  dxdt_ISOPROPANOL = (-Rate_ISOPROPANOL)*vol;
  // Bacteroides Degradation (g/h)
  dxdt_BACTEROID = (Rate_PEG9*Conc_BACTEROID + Rate_PEG8*Conc_BACTEROID + Rate_PEG7*Conc_BACTEROID + Rate_PEG6*Conc_BACTEROID +
    Rate_PEG5*Conc_BACTEROID + Rate_PEG4*Conc_BACTEROID + Rate_PEG3*Conc_BACTEROID + Rate_AcetHyde*Conc_BACTEROID) * vol;
  
  // Acidogenesis and Acetogenesis (g/h)
  dxdt_ACIDOGEN = ((U_Acido_Glucose*Conc_ACIDOGEN) * num_AcidogenReactions) * vol;
  
  dxdt_ACETOGEN = (U_Aceto_Eth*Conc_ACETOGEN + U_Aceto_PropA*Conc_ACETOGEN + Rate_DEG*Conc_ACETOGEN + Rate_EG*Conc_ACETOGEN) * vol + Aceto_H2*Conc_ACETOGEN*MW_H2*Y_H2_Acido +
    Aceto_H2*Conc_ACETOGEN*MW_CO2*Y_CO2_Acido ;
  
  // Components: Intermediate Products (Liquid) (mol/h)
  dxdt_ETHANOL = -Flow*Conc_ETHANOL + ((U_Acido_Glucose*Conc_ACIDOGEN/(MW_glucose*Y_Acido_Glucose))*2 - 
    (U_Aceto_Eth*Conc_ACETOGEN/(MW_ethanol*Y_Aceto_eth))) * vol;
  
  dxdt_PropAcid = -Flow*Conc_PropAcid + ((U_Acido_Glucose*Conc_ACIDOGEN/(MW_glucose*Y_Acido_Glucose))*2 - 
    (U_Aceto_PropA*Conc_ACETOGEN/(MW_propA*Y_Aceto_propA))) * vol;
  
  dxdt_ACETATE = -Flow*Conc_ACETATE + ((U_Acido_Glucose*Conc_ACIDOGEN/(MW_glucose*Y_Acido_Glucose))*3 + 
    (U_Aceto_Eth*Conc_ACETOGEN/(MW_ethanol*Y_Aceto_eth)) + (U_Aceto_PropA*Conc_ACETOGEN/(MW_propA*Y_Aceto_propA)) +
    Rate_EG*Conc_ACETOGEN/(MW_EG*Y_EG) + (Rate_DEG*Conc_ACETOGEN/(MW_DEG*Y_DEG))*2 - 
    (U_Meth_Acet*Conc_METHANOGEN/(MW_acetate*Y_Meth_Acet))) * vol + (Aceto_H2*Conc_ACETOGEN)/4;  
  
  dxdt_H2O = Flow * (Conc_water_IN - Conc_H2O) - ((U_Aceto_Eth*Conc_ACETOGEN/(MW_ethanol*Y_Aceto_eth))*2 -
    (U_Aceto_PropA*Conc_ACETOGEN/(MW_propA*Y_Aceto_propA))*3 - Rate_DEG*Conc_ACETOGEN/(MW_DEG*Y_DEG) -
    Rate_AcetHyde*Conc_BACTEROID/(MW_AcetHyde*Y_AcetHyde) + (Rate_METHANOL*Y_H2O_methanol) + 
    (Rate_ISOPROPANOL*Y_H2O_isopropanol)) * vol + (Aceto_H2*Conc_ACETOGEN)/2 + (Meth_H2*Conc_METHANOGEN)/2;
  
  // Methanogenesis (g/h)
  
  dxdt_METHANOGEN = Meth_H2*Conc_METHANOGEN*MW_H2*Y_H2_Ace + Meth_H2*Conc_METHANOGEN*MW_CO2*Y_CO2_Ace + 
    (U_Meth_Acet*Conc_METHANOGEN + (Rate_METHANOL*Y_methanol_X/MW_methanol) + 
    (Rate_ISOPROPANOL*Y_isopropanol_X/MW_isopropanol)) * vol;
  
  // Components: Final Products (Biogas) (mol/h)
  dxdt_CO2_TOT = ((U_Acido_Glucose*Conc_ACIDOGEN/(MW_glucose*Y_Acido_Glucose))*2 + (U_Acido_Glucose*Conc_ACIDOGEN/(MW_glucose*Y_Acido_Glucose))*2 +
    (U_Aceto_PropA*Conc_ACETOGEN/(MW_acetate*Y_Aceto_propA)) + (U_Meth_Acet*Conc_METHANOGEN/(MW_acetate*Y_Meth_Acet))) * vol -
    (Meth_H2_gas*Conc_METHANOGEN)/4 - (Aceto_H2_gas*Conc_ACETOGEN)/2; 
  
  dxdt_H2_TOT = (-(U_Acido_Glucose*Conc_ACIDOGEN/(MW_glucose*Y_Acido_Glucose))*2 + (U_Aceto_Eth*Conc_ACETOGEN/(Y_Aceto_eth*MW_ethanol))*3 + 
    (U_Aceto_PropA*Conc_ACETOGEN/(MW_acetate*Y_Aceto_propA)*3) + (Rate_EG*Conc_ACETOGEN/(MW_EG*Y_EG)) + 
    (Rate_DEG*Conc_ACETOGEN/(MW_DEG*Y_DEG))*2 - (Rate_METHANOL*Y_H2_methanol_X/Y_methanol_X) -
    (Rate_ISOPROPANOL*Y_H2_isopropanol_X/Y_isopropanol_X)*3) * vol - 
    Meth_H2_gas*Conc_METHANOGEN - Aceto_H2_gas*Conc_ACETOGEN; 
  
  dxdt_CH4_TOT =((U_Meth_Acet*Conc_METHANOGEN/(MW_acetate*Y_Meth_Acet)) + (Rate_METHANOL*Y_CH4_methanol) + 
    (Rate_ISOPROPANOL*Y_CH4_isopropanol)*3) * vol + (Meth_H2_gas*Conc_METHANOGEN)/4;
  
  // Components: Final Products (Biogas) (mol)
  double CO2_GAS = CO2_TOT/(1+(R_coeff*TempK*Henry_CO2*vol/Head_Space_VTOT));
  double H2_GAS = H2_TOT/(1+(R_coeff*TempK*Henry_H2*vol/Head_Space_VTOT));
  double CH4_GAS = CH4_TOT/(1+(R_coeff*TempK*Henry_CH4*vol/Head_Space_VTOT));
  // Components: Final Products (Biogas) (Pa)
  double CO2_Pressure = (CO2_GAS/Head_Space_VTOT) * R_coeff * TempK;
  double H2_Pressure = (H2_GAS/Head_Space_VTOT) * R_coeff * TempK;
  double CH4_Pressure = (CH4_GAS/Head_Space_VTOT) * R_coeff * TempK;
  
  // Components: Final Products (Liquid) (mol)
  double CO2_LIQ = CO2_Pressure * Henry_CO2 * vol; //CO2_GAS - CO2_Pressure * Henry_CO2 * vol;
  double H2_LIQ = H2_Pressure * Henry_H2 * vol; //H2_GAS - H2_Pressure * Henry_H2 * vol;
  double CH4_LIQ =  CH4_Pressure * Henry_CH4 * vol; //CH4_GAS - CH4_Pressure * Henry_CH4 * vol;
  
  // Guar (Actual Concentration in g/L)
  double GUAR_Conc = Conc_GUAR*MW_guar/n_bonds; // g/L
  
  //Pressure
  double Pressure_atm = (H2_Pressure + CO2_Pressure + CH4_Pressure + Pressure_Rest_IN)/101325; // Pressure in atm
  
  //Average MW of PEG Structures
  double molFrac_Peg9 = Conc_PEG9/C_TOT;
  double molFrac_Peg8 = Conc_PEG8/C_TOT;
  double molFrac_Peg7 = Conc_PEG7/C_TOT;
  double molFrac_Peg6 = Conc_PEG6/C_TOT;
  double molFrac_Peg5 = Conc_PEG5/C_TOT;
  double molFrac_Peg4 = Conc_PEG4/C_TOT;
  double molFrac_Peg3 = Conc_PEG3/C_TOT;
  double molFrac_DEG = Conc_DEG/C_TOT;
  
  double PEG_mol_frac = molFrac_DEG + molFrac_Peg3 + molFrac_Peg4 + molFrac_Peg5 + molFrac_Peg6 + molFrac_Peg7 + molFrac_Peg8 + molFrac_Peg9;
  double WT_AVG_PEG_MW = (molFrac_Peg9*MW_PEG_9 + molFrac_Peg8*MW_PEG_8 + molFrac_Peg7*MW_PEG_7 + molFrac_Peg6*MW_PEG_6 + 
    molFrac_Peg5*MW_PEG_5 + molFrac_Peg4*MW_PEG_4 + molFrac_Peg3*MW_PEG_3 + molFrac_DEG*MW_DEG)/PEG_mol_frac;
  
  double N_TOT = PEG_9 + PEG_8 + PEG_7 + PEG_6 + PEG_5 + PEG_4 + PEG_3 + DEG;
  double AVG_PEG_MW = (PEG_9*MW_PEG_9 + PEG_8*MW_PEG_8 + PEG_7*MW_PEG_7 + PEG_6*MW_PEG_6 + PEG_5*MW_PEG_5 + PEG_4*MW_PEG_4 + PEG_3*MW_PEG_3 + DEG*MW_DEG)/N_TOT;
  
$CAPTURE PEG_IN V_TOT Pressure_atm GUAR_Conc Conc_GLUCOSE Conc_ACIDOGEN Conc_ETHANOL Conc_PropAcid 
  Conc_ACETOGEN Conc_ACETATE Conc_METHANOGEN H2_LIQ H2_GAS CO2_LIQ CO2_GAS CH4_GAS CH4_LIQ 
  Conc_BACTEROID Conc_PEG9 Conc_PEG8 Conc_PEG7 Conc_PEG6 Conc_PEG5 Conc_PEG4 Conc_PEG3 Conc_DEG Conc_EG
  Conc_AcetHyde AVG_PEG_MW Conc_METHANOL Conc_ISOPROPANOL Rate_METHANOL Rate_ISOPROPANOL
    
$SET //req = s_(GUAR, H2_TOT, H2O)

