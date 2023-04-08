---
title: "Mathematical Model for Anaerobic Digestion in Batch Reactor"
date: "March 4th, 2019"
author: "Kyle Barrett"
header-includes: \usepackage{amsmath}
output:
  html_document:
   mathjax: "http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
   number_sections: true
params:
  set_title: "Mathematical Model for Anaerobic Digestion in Batch Reactor"
runtime: shiny
theme: cerulean
---
<style>
body {
text-align: justify}
</style>

<style type="text/css">

body{ /* Normal  */
      font-size: 15px;
  }
td {  /* Table  */
  font-size: 14px;
}
h1.title {
  font-size: 38px;
  color: DarkRed;
}
h1 { /* Header 1 */
  font-size: 26px;
  color: DarkBlue;
}
h2 { /* Header 2 */
    font-size: 19px;
  color: DarkRed;
}
h3 { /* Header 3 */
  font-size: 22px;
}
code.r{ /* Code block */
    font-size: 13px;
}
pre { /* Code block - determines code spacing between lines */
    font-size: 14px;
}
body .main-container {
max-width: 2000px;
}
</style>

# **Model Assumptions and Limitations**

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:center;"> Num </th>
   <th style="text-align:center;"> Comment </th>
   <th style="text-align:center;"> Effect of Addressing Comment </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> Ideally, each bacteria would be kept track of separately, and not as a category (i.e. acetogen, etc.) → Some strains participate in several included reactions, but have varying affinity per substrate. Others may only degrade a specific substrate. </td>
   <td style="text-align:center;"> Actual degredation rates could be much higher or lower, as some bacteria strains perform synergistically when together, while others may interfere with each other </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> Initial reactor temp is based on optimal conditions for only some of the bacteria (methanogens/some acetogen). However, all strains of bacteria can survive within the restricted temperature range. </td>
   <td style="text-align:center;"> Bacteria behavior would be more accurate. </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> Most kinetic parameters do not change with temperature due to lack of availability. </td>
   <td style="text-align:center;"> Model fit would be more reflective of temperature changes. </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> Kinetic parameters not fit to observational data. Based solely on kinetic paramaters from the literature that have been extrapolated to this environment. </td>
   <td style="text-align:center;"> Model would be MUCH more reflective of the process. </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> Some Parameters were adjusted by hand to more closely match the observed degradation times taken from the literature. </td>
   <td style="text-align:center;"> This is an approximation due to the lack of observational data and in some cases, kinietic paramaters. </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> Petroleum distillates breakdown into many components, including methanol and isopropanol. However to simplify the model, it was assumed that they only broke down to methanol and isopropanol. Thus the WT percents of each alcohol were increased to account for the presence of additional petroleum distillates. </td>
   <td style="text-align:center;"> Model would better represent the chosen contaminants present in hydraulic fracking fluid. </td>
  </tr>
</tbody>
</table>

# **Illustrations**
## **Proposed Pathway of Guar Gum**
<p align="center">
  <img src="www/Guar_Deg.png" alt="drawing" width="1050" height="760"/>
</p>

## **Proposed Pathway of PEG-400 (n ~ 8.7)**
<br>
<p align="center">
  <img src="www/PEG_Deg.png" alt="drawing" width="1050" height="760"/>
</p>

## **Proposed Pathway of Petroleum Distillates/Alcohols**
<br>
<p align="center">
  <img src="www/Alcohol_Deg.png" alt="drawing" width="910" height="560"/>
</p>

## **Overall Reaction Mechanism**
<br>
<p align="center">
  <img src="www/compartmentalModel.png" alt="drawing" width="1050" height="760"/>
</p>

# **Differential Equations (Guar Gum Only)**

  * Due to the large number of PDE's, only the equations relating to guar gum are provided mathematically below. See the model file for the rest (code).

## **Hydrolysis of Guar Gum**

<br>

$$
\mathbf{\frac{dGuar,bonds}{dt}}=\frac{dC_{Guar,bonds}}{dt} * vol=- \left(\frac{k_{cat}*C_e*C_{Guar,bonds}}{K_{m,G}+C_{Guar,bonds}}\right)*vol
$$

  * $C_{Guar, bonds}$ is molar concentration of *breakable bonds* of guar gum.
      
  * $C_{Guar,bonds}=C_{Guar}*N_{Bonds}$,    $\space\space\space N_{Bonds}=\frac{MW_{Guar}}{MW_{Glucose}} ≈ 6600\space bonds$,    $\space\space\space 	k_{cat}=\frac{v_{max,G}}{C_e}$

<br> 
     
$$
\mathbf{\frac{dGlucose}{dt}}= \left[\left(\frac{k_{cat}*C_e*C_{Guar,bonds}}{K_m+C_{Guar,bonds}}\right) -
Stoich \left(\frac{µ_{max,G}*C_{Glucose}}{Ks_G+C_{Glucose}}*\frac{X_{Acidogen}}{MW_{Glucose}*Y_{Acidogen, G}}\right)*3\space Reactions\right]*vol
$$
     
  * The acidogenic portion is **multiplied by 3** to account for the three acidogenic reactions that take place. The stoichiometry and Monod parameters happen to be the same in each case.
      
  *	Note: The parameters *$Stoich_{Acido, 1}$* and *$Stoich_{Acido, 2}$* in the **acidogenesis** section are technically equal, but are differentiated to draw attention to the fact that two reactions are utilizing/producing the specified component.
     
<br>

## **Acidogenesis and Acetogenesis** 

<br>

### **Liquids**
$$
\mathbf{\frac{dEthanol}{dt}}= \left[Stoich \left(\frac{µ_{max,G}*C_{Glucose}}{Ks_G+C_{Glucose}}*\frac{X_{Acidogen}}{MW_{Glucose}*Y_{Acidogen,   G}} \right) -Stoich \left(\frac{µ_{max,E}*C_{Ethanol}}{Ks_E+C_{Ethanol}}*\frac{X_{Acetogen}}{MW_{Ethanol}*Y_{Acetogen, E}}\right) \right]*vol
$$

<br>

$$
\mathbf{\frac{dPropA}{dt}}= \left[Stoich \left(\frac{µ_{max,G}*C_{Glucose}}{Ks_G+C_{Glucose}}*\frac{X_{Acidogen}}{MW_{Glucose}*Y_{Acidogen,   G}} \right) -Stoich \left(\frac{µ_{max,PA}*C_{PropA}}{Ks_{PA}+C_{PropA}}*\frac{X_{Acetogen}}{MW_{PropA}*Y_{Acetogen, PA}}\right) \right]*vol
$$

<br>

 

$$
\begin{align}
\mathbf{\frac{dAcetate}{dt}}= 
\begin{bmatrix}
Stoich \left(\frac{µ_{max,G}*C_{Glucose}}{Ks_G+C_{Glucose}}*\frac{X_{Acidogen}}{MW_{Glucose}*Y_{Acidogen,   G}} \right) \\\\+ 
Stoich \left(\frac{µ_{max,E}*C_{Ethanol}}{Ks_{E}+C_{Ethanol}}*\frac{X_{Acetogen}}{MW_{Ethanol}*Y_{Acetogen, E}}\right) \\\\+
Stoich \left(\frac{µ_{max,PA}*C_{PropA}}{Ks_{PA}+C_{PropA}}*\frac{X_{Acetogen}}{MW_{PropA}*Y_{Acetogen, PA}}\right) \\\\- 
Stoich \left(\frac{µ_{max,A}*C_{Acetate}}{Ks_{A}+C_{Acetate}}*\frac{X_{Acetogen}}{MW_{Acetate}*Y_{Methanogen, PA}}\right)\end{bmatrix}*vol + 
Stoich \left(\frac{v_{max,A}*\left(\left[H_2 \right]-\left[H_2 \right]^* \right)}{Km_A + \left[H_2 \right]-\left[H_2 \right]^* }*\frac{X_{Acetogen}}{X_{i,Acetogen}} \right)
\end{align}
$$

<br>

$$
\begin{align}
\mathbf{\frac{dWater}{dt}}= 
\begin{bmatrix}
-Stoich \left(\frac{µ_{max,E}*C_{Ethanol}}{Ks_{E}+C_{Ethanol}}*\frac{X_{Acetogen}}{MW_{Ethanol}*Y_{Acetogen, E}}\right) \\\\-
Stoich \left(\frac{µ_{max,PA}*C_{PropA}}{Ks_{PA}+C_{PropA}}*\frac{X_{Acetogen}}{MW_{PropA}*Y_{Acetogen, PA}}\right)  
\end{bmatrix}*vol\space + 
\begin{bmatrix}
Stoich \left(\frac{v_{max,A}*\left(\left[H_2 \right]-\left[H_2 \right]^* \right)}{Km_A + \left[H_2 \right]-\left[H_2 \right]^* }*\frac{X_{Acetogen}}{X_{i,Acetogen}} \right) +\\\\
Stoich \left(\frac{v_{max,M}*\left(\left[H_2 \right]-\left[H_2 \right]^* \right)}{Km_M + \left[H_2 \right]-\left[H_2 \right]^* }*\frac{X_{Methanogen}}{X_{i,Methanogen}} \right)\end{bmatrix}
\end{align}
$$
<br>

### **Gases**

  *	Note: The following PDEs keep track of total moles, i.e. *Liquid* + *Gas*. The partitioning of moles into each phase are determined by Henry's constant, and is described following the differential equations.

<br>

$$
\begin{align}
\mathbf{\frac{dCO_{2} (Tot)}{dt}}= 
\begin{bmatrix}
Stoich_{Acido,1} \left(\frac{µ_{max,G}*C_{Glucose}}{Ks_G+C_{Glucose}}*\frac{X_{Acidogen}}{MW_{Glucose}*Y_{Acidogen,   G}} \right) \\\\+ 
Stoich_{Acido,2} \left(\frac{µ_{max,G}*C_{Glucose}}{Ks_G+C_{Glucose}}*\frac{X_{Acidogen}}{MW_{Glucose}*Y_{Acidogen,   G}} \right) \\\\+ 
Stoich \left(\frac{µ_{max,PA}*C_{PropA}}{Ks_{PA}+C_{PropA}}*\frac{X_{Acetogen}}{MW_{PropA}*Y_{Acetogen, PA}}\right) \\\\+ 
Stoich \left(\frac{µ_{max,A}*C_{Acetate}}{Ks_{A}+C_{Acetate}}*\frac{X_{Acetogen}}{MW_{Acetate}*Y_{Methanogen, PA}}\right)\end{bmatrix}*vol\space - 
\begin{bmatrix}
Stoich \left(\frac{v_{max,A}*\left(\left[H_2 \right]-\left[H_2 \right]^* \right)}{Km_A + \left[H_2 \right]-\left[H_2 \right]^* }*\frac{X_{Acetogen}}{X_{i,Acetogen}} \right)+ \\\\
Stoich \left(\frac{v_{max,M}*\left(\left[H_2 \right]-\left[H_2 \right]^* \right)}{Km_M + \left[H_2 \right]-\left[H_2 \right]^* }*\frac{X_{Methanogen}}{X_{i,Methanogen}} \right)
\end{bmatrix}
\end{align}
$$
  *	See previous comment regarding the parameters *$Stoich_{Acido, 1}$* and *$Stoich_{Acido, 2}$*
  
<br>
  
$$
\begin{align}
\mathbf{\frac{dH_{2}(Tot)}{dt}}= 
\begin{bmatrix}
Stoich_{Acido,1} \left(\frac{µ_{max,G}*C_{Glucose}}{Ks_G+C_{Glucose}}*\frac{X_{Acidogen}}{MW_{Glucose}*Y_{Acidogen,   G}} \right) \\\\+ 
Stoich \left(\frac{µ_{max,E}*C_{Ethanol}}{Ks_{E}+C_{Ethanol}}*\frac{X_{Acetogen}}{MW_{Ethanol}*Y_{Acetogen, E}}\right) \\\\+
Stoich \left(\frac{µ_{max,PA}*C_{PropA}}{Ks_{PA}+C_{PropA}}*\frac{X_{Acetogen}}{MW_{PropA}*Y_{Acetogen, PA}}\right) \\\\
\end{bmatrix}*vol\space - 
\begin{bmatrix}
Stoich \left(\frac{v_{max,A}*\left(\left[H_2 \right]-\left[H_2 \right]^* \right)}{Km_A + \left[H_2 \right]-\left[H_2 \right]^* }*\frac{X_{Acetogen}}{X_{i,Acetogen}} \right)+ \\\\
Stoich \left(\frac{v_{max,M}*\left(\left[H_2 \right]-\left[H_2 \right]^* \right)}{Km_M + \left[H_2 \right]-\left[H_2 \right]^* }*\frac{X_{Methanogen}}{X_{i,Methanogen}} \right)
\end{bmatrix}
\end{align}
$$
  
<br>



### **Gas-Liquid Equilibrium**

$$
\begin{align}
\mathbf{CO_2 (Gas)}=\frac{CO_2 (Tot)}{\left[1+\frac{R*T*Henry_{CO_{2}}*vol}{vol_{Head Space} }\right]} ,&&
\mathbf{P_{CO_2}}=\left(\frac{CO_2 (Gas)}{vol_{Head Space}}\right)*R*T
\end{align}
$$

<br>
$$
\begin{align}
\mathbf{H_2 (Gas)}=\frac{H_2 (Tot)}{\left[1+\frac{R*T*Henry_{H_{2}}*vol}{vol_{Head Space} }\right]} ,&&
\mathbf{P_{H_2}}=\left(\frac{H_2 (Gas)}{vol_{Head Space}}\right)*R*T
\end{align}
$$

<br>
$$
\begin{align}
\mathbf{CO_2 (Liq)}=P_{CO_2}*Henry_{CO_{2}}*vol ,&&
\mathbf{H_2 (Liq)}=P_{H_2}*Henry_{H_{2}}*vol
\end{align}
$$

<br>

## **Methanogenesis**

<br>

$$
\mathbf{\frac{dCH_{4}(Tot)}{dt}}= Stoich \left(\frac{µ_{max,A}*C_{Acetate}}{Ks_{A}+C_{Acetate}}*\frac{X_{Acetogen}}{MW_{Acetate}*Y_{Methanogen, PA}}\right)*vol + 
Stoich \left(\frac{v_{max,M}*\left(\left[H_2 \right]-\left[H_2 \right]^* \right)}{Km_M + \left[H_2 \right]-\left[H_2 \right]^* }*\frac{X_{Methanogen}}{X_{i,Methanogen}} \right)
$$
<br>

### **Gas-Liquid Equilibrium** 

<br>
$$
\begin{align}
\mathbf{CH_4 (Gas)}=\frac{CH_4 (Tot)}{\left[1+\frac{R*T*Henry_{CH_{4}}*vol}{vol_{Head Space} }\right]} ,&&
\mathbf{P_{CH_4}}=\left(\frac{CH_4 (Gas)}{vol_{Head Space}}\right)*R*T
\end{align}
$$

<br>

$$
\mathbf{CH_4 (Liq)}=P_{CH_4}*Henry_{CH_4}*vol
$$
<br>

## **Bacteria Growth** 

<br>
$$
\begin{align}
\mathbf{\frac{dAcidogen}{dt}}= 3 *
\left[\frac{µ_{max,G}*C_{Glucose}}{Ks_G+C_{Glucose}}*X_{Acidogen} \right]*vol
\end{align}
$$
  *	The **'3'** accounts for each of the acidogenic reactions occuring on the subtrate (Glucose)

<br>
$$
\begin{align}
\mathbf{\frac{dAcetogen}{dt}}=
&\left[\frac{µ_{max,E}*C_{Ethanol}}{Ks_E+C_{Ethanol}}*X_{Acetogen}+
\frac{µ_{max,PA}*C_{PropA}}{Ks_{PA}+C_{PropA} }*X_{Acetogen} \right]*vol \\\\+
&\left[\frac{v_{max,A}*\left(\left[H_2 \right]-\left[H_2 \right]^* \right)}{Km_A + \left[H_2 \right]-\left[H_2 \right]^* }*\frac{X_{Acetogen}}{X_{i,Acetogen}} \right]*MW_{H_2}*Y_{Acetogen, H_2}
\end{align}
$$

<br>
$$
\begin{align}
\mathbf{\frac{dMethanogen}{dt}}=
\left[\frac{µ_{max,A}*C_{Acetate}}{Ks_A+C_{Acetate}}*X_{Methanogen} \right]*vol +
\left[\frac{v_{max,M}*\left(\left[H_2 \right]-\left[H_2 \right]^* \right)}{Km_M + \left[H_2 \right]-\left[H_2 \right]^* }*\frac{X_{Methanogen}}{X_{i,Methanogen}} \right]*MW_{H_2}*Y_{Methanogen, H_2}
\end{align}
$$

<br>
