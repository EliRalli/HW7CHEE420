%Nonisothermal PFR for an non-reversible gas phase elementary reaction

clear all; close all;

%Conversion range for the integration
volStart=0;%L
volEnd=3000;%L
volSpan=[volStart,volEnd];

%Initial conditions
initX=0;

%ODE Solver initialization:
[volOut,X]=ode45(@myFunction,volSpan,initX);

%Plotting Data
plot(volOut,X)

%Final Conversion of A
Xend=X(end,:)

%Oulet temperature at the End
DeltaH=-70000;%J/mol
To=500;%K
va=1;
vb=2;
vc=1;
vi=0;
Cpa=150;%J/(mol*K)
Cpb=150;%J/(mol*K)
Cpi=150;%J/(mol*K)
Cpc=450;%J/(mol*K)
yao=0.2;
ybo=0.6;
yco=0;
yio=0.2;
theta_a=1;
theta_b=ybo/yao;
theta_c=yco/yao;
theta_i=yio/yao;
Sum_Cpi_theta=Cpa+Cpb*theta_b+Cpc*theta_c+Cpi*theta_i;%Sum of the Cp data
Delta_Cp=vc*Cpc-(vb*Cpb+va*Cpa);%Cp(products)-Cp(reactants)
Tend=(Xend*(-DeltaH)+Sum_Cpi_theta*To+Xend*Delta_Cp)/(Sum_Cpi_theta+Xend*Delta_Cp)

%Inputting the differential equation(s) and parameters
function f = myFunction(vol,X);

%Reaction rate constant k and other parameters
E=500;%J/mol
R=8.314472;%J/(mol*K)
k_300K=0.0055;%L^2/(mol^2*sec)


%Initial volumetric flowrate
vo=50;%L/sec

%Cp Data
Cpa=150;%J/(mol*K)
Cpb=150;%J/(mol*K)
Cpi=150;%J/(mol*K)
Cpc=450;%J/(mol*K)

%Initial molar flowrate
Fto=5;%mol/sec

%Initial mole fraction of the gas phase
yao=0.2;
ybo=0.6;
yco=0;
yio=0.2;

%Initial molar flowrates
Fao=Fto*yao;
Fbo=Fto*ybo;
Fio=Fto*yio;
Fco=Fto*yco;
Fto=Fao+Fbo+Fco+Fio;

%Initial molar concentration
Cao=Fao/vo;
Cbo=Fbo/vo;
Cco=Fco/vo;
Cio=Fio/vo;

%The parameter
theta_a=1;
theta_b=ybo/yao;
theta_c=yco/yao;
theta_i=yio/yao;

%Change in the total number of moles per mole of reactant A
va=1;
vb=2;
vc=1;
vi=0;

delta=vc/va-(vb/va+va);

%Enthalply Heat of Reaction At Standard Conditions
DeltaH=-70000;%J/mol
%Outlet molar flowrate

Fa=Fao*(1-X);
Fb=Fao*(theta_b-vb/va*X);
Fc=Fao*(theta_c-vc/va*X);
Fi=Fao*(theta_i-vi/va*X);
Ft=Fa+Fb+Fc+Fi;

%Outlet volumentric flowrate
Z=1;
Zo=1;
Po=1;
P=1;
To=500;%K
%Oulet temperature as a function of conversion
Sum_Cpi_theta=Cpa+Cpb*theta_b+Cpc*theta_c+Cpi*theta_i;%Sum of the Cp data
Delta_Cp=vc*Cpc-(vb*Cpb+va*Cpa);%Cp(products)-Cp(reactants)
T=(X*(-DeltaH)+Sum_Cpi_theta*To+X*Delta_Cp)/(Sum_Cpi_theta+X*Delta_Cp);
epsilon=yao*delta;
v=vo*(Po/P)*(T/To)*(Z/Zo)*(1+epsilon*X);

%Outlet concentration
Ca=Fa/v;
Cb=Fb/v;
Cc=Fc/v;
Ci=Fi/v;

%Reaction rate constant as a function of temperature
k=k_300K*exp(E/R*(1/To-1/T));

%Rate law
ra=-k*Ca*Cb^2;
dXdV=-ra/Fao;
f=[dXdV];
end