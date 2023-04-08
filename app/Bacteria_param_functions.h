// Date: Fri Jan  25 09:04:25 2019
#ifndef Bacteria_param_functions_h
#define Bacteria_param_functions_h
/*** R
source("/Users/kyleb/Desktop/Classwork/Senior Design/Modeling/R_Files/Model/Scripts/Bacteria_kinetics.r")
  coef_Aceto_vmax <<- coef(fit_Aceto_vmax)
  Intercept_Aceto_vmax <<- as.numeric(coef_Aceto_vmax[1])
  a_Aceto_vmax <<- as.numeric(coef_Aceto_vmax[2])
  b_Aceto_vmax <<- as.numeric(coef_Aceto_vmax[3])
  c_Aceto_vmax <<- as.numeric(coef_Aceto_vmax[4])

  coef_Aceto_km <<- coef(fit_Aceto_km)
  Intercept_Aceto_km <<- as.numeric(coef_Aceto_km[1])
  a_Aceto_km <<- as.numeric(coef_Aceto_km[2])
  b_Aceto_km <<- as.numeric(coef_Aceto_km[3])
  c_Aceto_km <<- as.numeric(coef_Aceto_km[4])
  
  coef_Meth_vmax <<- coef(fit_Meth_vmax)
  Intercept_Meth_vmax <<- as.numeric(coef_Meth_vmax[1])
  a_Meth_vmax <<- as.numeric(coef_Meth_vmax[2])
  b_Meth_vmax <<- as.numeric(coef_Meth_vmax[3])
  c_Meth_vmax <<- as.numeric(coef_Meth_vmax[4])
  
  coef_Meth_km <<- coef(fit_Meth_km)
  Intercept_Meth_km <<- as.numeric(coef_Meth_km[1])
  a_Meth_km <<- as.numeric(coef_Meth_km[2])
  b_Meth_km <<- as.numeric(coef_Meth_km[3])
  c_Meth_km <<- as.numeric(coef_Meth_km[4])
  */

namespace param_fn {
double polyn(double Temp2, double intercept, double a, double b, double c) {
  return(intercept + a*Temp2 + b*pow(Temp2,2) + c*pow(Temp2,3));
}

double Aceto_vmax(double Temp2, const double Intercept_Aceto_vmax,const double a_Aceto_vmax, const double b_Aceto_vmax, const double c_Aceto_vmax) {
  return(polyn(Temp2,Intercept_Aceto_vmax,a_Aceto_vmax,b_Aceto_vmax,c_Aceto_vmax));
}

double Aceto_km(double Temp2, const double Intercept_Aceto_km,const double a_Aceto_km, const double b_Aceto_km, const double c_Aceto_km) {
  return(polyn(Temp2,Intercept_Aceto_km,a_Aceto_km,b_Aceto_km,c_Aceto_km));
}

double Meth_vmax(double Temp2, const double Intercept_Meth_vmax,const double a_Meth_vmax, const double b_Meth_vmax, const double c_Meth_vmax) {
  return(polyn(Temp2,Intercept_Meth_vmax,a_Meth_vmax,b_Meth_vmax,c_Meth_vmax));
}

double Meth_km(double Temp2, const double Intercept_Meth_km,const double a_Meth_km, const double b_Meth_km, const double c_Meth_km) {
  return(polyn(Temp2,Intercept_Meth_km,a_Meth_km,b_Meth_km,c_Meth_km));
}

}

#endif