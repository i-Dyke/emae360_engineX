function [problem,time]=CEA(varargin)

% August 16, 2024
% Version 1.3.0.0
%
% Program modified by Peter Nelson, Summer Intern, ESSCA
% Program modified by Christian Reynolds, Summer Intern, ESSCA
% Program modified by Chip Kopicz, MSFC ER13, ESSCA
%
% Program corrected to handle single fluid call for fluid property
% calculations. Works with gasses and liquids. Solids will not work.
% Graphical User Interface (GUI) provided with program.
% Added an explanation on how to use executable version and capability for
% the cell array input.
%
% Known Issues
%   If Trace and ions are specified in the problem call, the code will 
%   generate an error.
%
% July 22, 2024
% Version 1.2.0.0
%
% Program modified by Peter Nelson, Summer Intern, ESSCA
% Program modified by Chip Kopicz, MSFC ER13, ESSCA
%
% Program corrected to handle thermodynamic library with no condensed
% species.
% Program corrected to record weight percent of a reactant when the weight
% percent is entered by the program.
%
% June 1, 2024
% Version 1.1.0.0
%
% Program modified by Peter Nelson, Summer Intern, ESSCA
% Program modified by Austin Jones, Summer Intern, ESSCA
% Program modified by Chip Kopicz, MSFC ER13, ESSCA
%
% Program modified to create sublibrary for a propellant combination and is
% no longer limited to the pre-programmed libraries.
% Issues with processing condensed species have been addressed for solid
% propellant and low mixture ratio conditions.
% The species in the thermodynamic library have been rearranged with all
% gasses are grouped together and IFAZ values in order for MatLab to
% process. THE NEW THERMODYNAMIC LIBRARY PROVIDED WITH THIS RELEASE IS
% REQUIRED
% Code tracks each species input. Large number of inputs in very small
% values in previous version created issues with calculation roundoff.
% Multiple mixtures handled consistently.
% Number of Asub, Asup, and Pi/P increased beyond previous limit of 10
% total.
% Issues with calculating and printing transport properties corrected.
% Added unit distinction for kg-mol and gm-mol.
% Warnings are written to structured array. Errors still terminate program.
% SHOCK calculations were developed for shock tubes. The authors have
% observed good agreement with FORTRAN CEA up to Mach 22, but results were
% dependent on reactant combinations. Users should plan to verify
% calculations for velocities and Mach numbers in excess of Mach 10.
%
% February 7, 2017
% Version 1.0.0.0
%
% Program debugged, verified, and validated by
% Chip Kopicz
% ERC
%
% Program Translated from FORTRAN and Constructed in MATLAB by
% Sam Kobliska
% ERC., Jacobs ESSSA
%
%%%% Enthalpy-Pressure Examples
% HP1=CEA('problem','hp','equilibrium','o/f',6.0338,'case','CEAM-HP1','p,psia',3126.24,'reactants','fuel','H2(L)','H',2,'wt%',100,'h,cal/mol',-2154,'t(k)',20.17, 'oxid','O2(L)','O',2,'wt%',100,'h,cal/mol',-3102,'t(k)',90.18,'output','transport','end','screen');
% HP2=CEA('reac','oxid','Air','wtfrac',1,'t(k)',700,'fuel','C7H8(L)','wtfrac',0.4,'t(k)',298.15,'fuel','C8H18(L),n-octa','wtfrac',0.6,'t(k)',298.15,'prob','case','Example-3','hp','o/f',17,'p,bar',100,10,1,'omit','CH2CO,ketene','Air','HO(CO)2OH','CH2','CH3','CH2OH','CH3O','CH4','CH3OH','CCN','CNC','C2N2','C2O','C3H4,allene', 'C3H4,propyne','C3H4,cyclo-','C3','C3H5,allyl','C3H6,propylene','C3H6,cyclo-','C3H3,propargyl','C3H6O','C3H7,n-propyl','C3H7,i-propyl','Jet-A(g)','C3O2','C4','C4H2','C3H8O,2propanol','C4H4,1,3-cyclo-','C4H6,butadiene','C4H6,2butyne','C3H8O,1propanol', 'C4H8,tr2-butene','C4H8,isobutene','C4H8,cyclo-','C4H6,cyclo','(CH3COOH)2','C4H9,n-butyl','C4H9,i-butyl','C4H8,l-butene', 'C4H9,s-butyl','C4H9,t-butyl','C4H10,isobutane','C4H8,cis2-buten','C4H10,n-butane','C4N2','C5','C3H8','C5H6,l,3cyclo-', 'C5H8,cyclo-','C5H10,l-pentene','C10H21,n-decyl','C5H10,cyclo-','C5H11,pentyl','C5H11,t-pentyl','C12H10,biphenyl','C5H12,n-pentane','C5H12,i-pentane','CH3C(CH3)2CH3','C12H9,o-bipheny','C6H6','C6H5OH,phenol','C6H10,cyclo-','C6H2','C6H12,l-hexane','C6H12,cyclo-','C6H13,n-hexyl','C6H5,phenyl','C7H7,benzyl','C7H8','C7H8,cresol-mx','C6H5O,phenoxy','C7H14,l-heptane','C7H15,n-heptyl','C7H16,n-heptane','C10H8,azulene','C8H8,atyrene','C8H10,ethylbenz','C8H16,l-octene', 'C10H8,napthlene','C8H17,n-octyl','C8H18,isooctane','C8H18,n-octane','C9H19,n-octyl','C7H8(L)','C8H18(L),n-octa','Jet-A(L)','C6H6(L)','H2O(s)','H2O(L)','output','mks','trace',1e-15,'end','screen');
% HP3=CEA('reac','name','NH4CLO4(I)','wt%',72.06,'t(k)',298.15,'name','CHOS-Binder','C',1,'H',1.86955,'O',0.031256,'S',0.008415,'wt%',18.58, 'h,cal/mol',-2999.082,'t(k)',298.15,'name','AL(cr)','wt%',9,'t(k)',298.15,'name','MgO(s)','Mg',1,'O',1,'h,cal/mol',-143703,'wt%',0.2,'t(k)',298.15,'name','H2O(L)','wt%',0.16,'t(k)',298.15,'prob','case','Example-5','hp','p,psi',500,250,125,50,5, 'omit','CCN','CNC','C2N2','C2O', 'C3H4,allene','C3H4,propyne','C3H4,cyclo-','C3','C3H5,allyl','C3H6,propylene','C3H6,cyclo-','C3H3,propargyl','C3H6O','C3H7,n-propyl', 'C3H7,i-propyl','Jet-A(g)','C3O2','C4','C4H2','C3H8O,2propanol','C4H4,1,3-cyclo-','C4H6,butadiene','C4H6,2butyne','C3H8O,1propanol', 'C4H8,tr2-butene','C4H8,isobutene','C4H8,cyclo-','C4H6,cyclo','(CH3COOH)2','C4H9,n-butyl','C4H9,i-butyl','C4H8,l-butene','C4H9,s-butyl','C4H9,t-butyl','C4H10,isobutane','C4H8,cis2-buten','C4H10,n-butane','C4N2','C5','C3H8','C5H6,l,3cyclo-','C5H8,cyclo-','C5H10,l-pentene','C10H21,n-decyl','C5H10,cyclo-','C5H11,pentyl','C5H11,t-pentyl','C12H10,biphenyl','C5H12,n-pentane','C5H12,i-pentane', 'CH3C(CH3)2CH3','C12H9,o-bipheny','C6H6','C6H5OH,phenol','C6H10,cyclo-','C6H2','C6H12,l-hexane','C6H12,cyclo-','C6H13,n-hexyl','C6H5,phenyl','C7H7,benzyl','C7H8','C7H8,cresol-mx','C6H5O,phenoxy','C7H14,l-heptane','C7H15,n-heptyl','C7H16,n-heptane', 'C10H8,azulene','C8H8,atyrene','C8H10,ethylbenz','C8H16,l-octene','C10H8,napthlene','C8H17,n-octyl','C8H18,isooctane','C8H18,n-octane','C9H19,n-octyl','C7H8(L)','C8H18(L),n-octa','Jet-A(L)','C6H6(L)','H2O(s)','H2O(L)','output','calories','end','screen');
% HP4=CEA('reac','name','SR_AL','AL',1.00,'h,cal/mol',0.,'wt%',18.0,'t(k)',298.15,'name','AM_PERCL','N',1,'H',4,'O',4,'CL',1,'wt%',68.00,'h,cal/mol',-70700.,'t(k)',298.15,'name','Dio_Adi','C',22.,'H',42.,'O',4.,'wt%',2.0,'h,cal/mol',-296000.,'t(k)',298.15,'name','HTPB','C',7.332,'H',10.962,'O',0.058,'h,cal/mol',-250.,'wt%',11.00,'t(k)',298.15,'name','HX-752','C',14.,'H',16.,'N',2.0,'O',2.,'wt%',0.20,'h,cal/mol',-61000.,'t(k)',298.15,'name','stabilizer','C',38.,'H',66.,'N',2., 'wt%',0.10,'h,cal/mol',-206300.,'t(k)',298.15,'prob','hp','p,psia',1000,'outp','massf','transport','mks','end','screen');
% HP5=CEA('reac','name','AM_PERCL','N',1,'H',4,'O',4,'CL',1,'wt%',68.00,'h,cal/mol',-70700.,'t(k)',298.15,'name','SR_AL','AL',1.00,'h,cal/mol',0.,'wt%',18.0,'t(k)',298.15,'name','Dio_Adi','C',22.,'H',42.,'O',4.,'wt%',2.0,'h,cal/mol',-296000.,'t(k)',298.15,'name','HTPB','C',7.332,'H',10.962,'O',0.058,'h,cal/mol',-250.,'wt%',11.00,'t(k)',298.15,'name','HX-752','C',14.,'H',16.,'N',2.0,'O',2.,'wt%',0.20,'h,cal/mol',-61000.,'t(k)',298.15,'name','stabilizer','C',38.,'H',66.,'N',2.,'wt%',0.10,'h,cal/mol',-206300.,'t(k)',298.15, 'prob','hp','p,psia',1000,'outp','massf','mks','end','screen');
% HP6=CEA('problem','hp','o/f',136.9,'t(k)',3800','case','Example','p,psia',1828,'reactants','fuel','RPCust','wt%',100,'h,kj/mol',-20.57,'C',1,'H',1.95,'oxid','O2Cust','wt%',100,'h,kj/mol',-11.3927,'O',2,'end','screen');
%
%%%% Temperature-Pressure Examples
% TP1=CEA('problem','tp','equilibrium','o/f',6.0:.2:7.0,'case','CEAM-TP1','p,psi',1000.,2000,'t,R',5500:500:6500,'outp','transport','reactants','fuel','H2','wt%',100.,'t(k)',298.15,'oxid','O2','wt%',100,'t(k)',298.15,'end','screen');
% TP2=CEA('problem','case','Example-1','tp','p(atm)',1,0.1,0.01,'t(k)',3000,2000,'r.eq.ratio',1,1.5,'reac','fuel','H2','moles',1,'oxid','Air','moles',1,'only','Ar','C','CO','CO2','H','H2','H2O','HO2','HNO2','HNO3','N','NH','NO','N2','N2O3','O','O2','OH','O3','output','calories','end','screen');
% TP3=CEA('reactants','name','H2(L)','H',2,'moles',100,'h,J/mol',-9012,'name','O2(L)','O',2,'moles',60,'h,J/mol',-12979,'prob','tp','case','Example-14','p,atm',0.05,'t(k)',1000,500,350,305,304.3,304.2,304,300,'output','debug',5,'end','screen');
% TP4=CEA('problem','tp','equilibrium','t(k)', 527.07, 'p(atm)', 1,'phi', .625,'omit','Air', 'reactants','fuel','CH4','moles',.86,'t(k)',297.65,'fuel','N2','moles',.14,'ox','N2','moles',.769916,'ox','O2','moles',.206459,'ox','Ar','moles',.00920135,'ox','CO2','moles',.0003944,'ox','Ne','moles',1.79452E-5,'ox','He','moles',5.1272E-6,'ox','CH4','moles',1.479E-6,'ox','Kr','moles',1.0846E-6,'ox','H2','moles',4.93E-7,'ox','N2O','moles',2.958E-7,'ox','CO','moles',1.972E-7,'ox','Xe','moles',9.86E-8,'ox','H2O','moles',.014,'output','transport','end','screen');
% TP5=CEA('problem','tp','o/f',136.9,'t(k)',3800','case','Example','p,psia',1828,'reactants','fuel','RPCust','wt%',100,'h,kj/mol',-20.57,'C',1,'H',1.95,'oxid','O2Cust','wt%',100,'h,kj/mol',-11.3927,'O',2,'end','screen');
%
%%%% Entropy-Pressure Example
% SP1=CEA('reactants','fuel','H2(L)','H',2,'wt%',100.0,'oxid','O2(L)','O',2,'wt%',100.0,'only','H','H2','H2O','O','OH','O2','problem','sp','case','CEAM-SP1','o/f',6.00,'p,psi',628.960,'s/r',2.17959,'output','short','massf','end','screen');
% SP2=CEA('reactants','fuel','H2(L)','H',2,'wt%',100.0,'oxid','O2(L)','O',2,'wt%',100.0,'only','H','H2','H2O','O','OH','O2','problem','sp','o/f',6.00,'p,bar',43.36527,'s,kJ/kg-K',18.12328,'output','short','massf','end','screen');
%
%%%% Entropy-Specific Volume Example
% SV1=CEA('problem','sv','equilibrium','o/f',6.0338,'case','CEAM-SV1', 's,kj/kg-k',17.0993,'dens,kg/m3',9.823,'reactants', 'fuel','H2(L)','H',2,'wt%',100,'h,cal/mol',-2154,'t(k)',20.17, 'oxid','O2(L)','O',2,'wt%',100,'h,cal/mol',-3102,'t(k)',90.18, 'output','end','screen');
% SV2=CEA('problem','sv','equilibrium','o/f',6.0338,'case','CEAM-SV2', 's,kj/kg-k',17.0993,'v,m3/kg',1/9.823,'reactants', 'fuel','H2(L)','H',2,'wt%',100,'h,cal/mol',-2154,'t(k)',20.17, 'oxid','O2(L)','O',2,'wt%',100,'h,cal/mol',-3102,'t(k)',90.18, 'output','end','screen');
% SV3=CEA('problem','sv','equilibrium','o/f',6.0338,'case','CEAM-SV3', 's/r',2.056574,'dens,kg/m3',9.823,'reactants', 'fuel','H2(L)','H',2,'wt%',100,'h,cal/mol',-2154,'t(k)',20.17, 'oxid','O2(L)','O',2,'wt%',100,'h,cal/mol',-3102,'t(k)',90.18, 'output','end','screen');
%
%%%% Temperature-Specific Volume Example
% TV1=CEA('reac','fuel','H2','wt%',100,'oxid','Air','wt%',100,'problem','case','Example-2','phi,eq.ratio',1,'tv','t(k)',3000,'rho,g/cc',9.186e-5,8.0877e-6,6.6054e-7,'only', 'Ar','C','CO','CO2','H','H2','H2O','HO2','HNO2','HNO3','N','NH','NO','N2','N2O3','O','O2','OH','O3', 'output','transport','calories','end','screen');
% TV2=CEA('reac','fuel','H2','wt%',100,'oxid','Air','wt%',100,'problem','case','CEAM-TV2','phi,eq.ratio',1,'tv','t(k)',3000,'v,cc/g',1/9.186e-5,1/8.0877e-6,1/6.6054e-7,'only', 'Ar','C','CO','CO2','H','H2','H2O','HO2','HNO2','HNO3','N','NH','NO','N2','N2O3','O','O2','OH','O3', 'output','transport','calories','end','screen');
% TV3=CEA('reac','fuel','H2','wt%',100,'oxid','Air','wt%',100,'t(F)',800.33,'problem','case','CEAM-TV3','phi,eq.ratio',1,'tv','t,c',2726.85,'v,in3/kg',1/1.505315699e-6,'only', 'Ar','C','CO','CO2','H','H2','H2O','HO2','HNO2','HNO3','N','NH','NO','N2','N2O3','O','O2','OH','O3', 'output','transport','mks','end','screen');
%
% % Internal Energy-Specific Volume Example
% UV1=CEA('reac','oxid','Air','wtfrac',1,'t(k)',700,'fuel','C7H8(L)','wtfrac',0.4,'t(k)',298.15,'fuel','C8H18(L),n-octa','wtfrac',0.6,'t(k)', 298.15,'prob','case','Example-4','uv','o/f',17,'u/r',-45.1343,'rho,kg/m^3',14.428,'omit','CH2CO,ketene','Air','HO(CO)2OH','CH2','CH3','CH2OH','CH3O','CH4','CH3OH','CCN','CNC','C2N2','C2O','C3H4,allene', 'C3H4,propyne','C3H4,cyclo-','C3','C3H5,allyl','C3H6,propylene','C3H6,cyclo-','C3H3,propargyl','C3H6O','C3H7,n-propyl','C3H7,i-propyl','Jet-A(g)','C3O2','C4','C4H2','C3H8O,2propanol','C4H4,1,3-cyclo-','C4H6,butadiene','C4H6,2butyne','C3H8O,1propanol', 'C4H8,tr2-butene','C4H8,isobutene','C4H8,cyclo-','C4H6,cyclo','(CH3COOH)2','C4H9,n-butyl','C4H9,i-butyl','C4H8,l-butene', 'C4H9,s-butyl','C4H9,t-butyl','C4H10,isobutane','C4H8,cis2-buten','C4H10,n-butane','C4N2','C5','C3H8','C5H6,l,3cyclo-', 'C5H8,cyclo-','C5H10,l-pentene','C10H21,n-decyl','C5H10,cyclo-','C5H11,pentyl','C5H11,t-pentyl','C12H10,biphenyl','C5H12,n-pentane','C5H12,i-pentane','CH3C(CH3)2CH3','C12H9,o-bipheny','C6H6','C6H5OH,phenol','C6H10,cyclo-','C6H2','C6H12,l-hexane','C6H12,cyclo-','C6H13,n-hexyl','C6H5,phenyl','C7H7,benzyl','C7H8','C7H8,cresol-mx','C6H5O,phenoxy','C7H14,l-heptane','C7H15,n-heptyl','C7H16,n-heptane','C10H8,azulene','C8H8,atyrene','C8H10,ethylbenz','C8H16,l-octene', 'C10H8,napthlene','C8H17,n-octyl','C8H18,isooctane','C8H18,n-octane','C9H19,n-octyl','C7H8(L)','C8H18(L),n-octa','Jet-A(L)','C6H6(L)','H2O(s)','H2O(L)','output','mks','trace',1e-15,'end','screen');
% UV2=CEA('problem','uv','equilibrium','U,BTU/lbm',-2466.882733575299,'rho,lbm/ft3',0.236079691999656,'outp','transport','eng','reactants','fuel','H2','wt%',100.,'t(k)',298.15,'oxid','O2','wt%',100,'t(k)',298.15,'end','screen');
% UV3=CEA('problem','uv','equilibrium','U,BTU/lbm',-2466.882733575299,'v,ft3/lbm',1/0.236079691999656,'outp','transport','eng','reactants','fuel','H2','wt%',100.,'t(k)',298.15,'oxid','O2','wt%',100,'t(k)',298.15,'end','screen');
%
%%%% Rocket Examples
% RKT1=CEA('problem','rocket','equilibrium','fac','acat',3,'o/f',2.34,2.5,'case','CEAM-RKT1','p(psi)',633,800,'subsonic(ae/at)',2.59,1.01,'supsonic(ae/at)',1.01,15.0,30.0,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'h,cal/mol',-5430,'t(k)',300.0,'oxid','O2(L)','O',2,'wt%',100,'h,cal/mol',-3032,'t(k)',94.44,'output','transport','mks','end','screen');
% RKT2=CEA('reactants','name','HEHN','wt',44.500,'C',2,'H',9,'N',3,'O',4,'h,kcal/mol',-98.000,'rho,g/cc',1.428,'name', 'HAN','wt',41.830,'H',4,'N',2,'O',4,'h,kcal/mol',-95.300,'rho,g/cc',1.090, 'name','AN','wt',02.225,'H',4,'N',2,'O',3,'h,kcal/mol',-87.380, 'rho,g/cc',1.725, 'name','DIPY','wt',0.445,'C',10,'H',8,'N',2,'h,kcal/mol',44.478,'rho,g/cc',1.326, 'name','Water','wt',11.000,'H',2,'O',1,'h,kcal/mol',-68.313,'rho,g/cc',1.000, 'problem','rocket','frozen','p,psia',300,'sup,ae/at',10.0000,'output','massf','transport','end','screen');
% RKT3=CEA('problem','rocket','equilibrium','o/f',5.55157,'case','Example-8','p,bar',53.3172,'subar',1.58,'pi/p',10,100,1000,'supar',25,50,75,'reactants','fuel','H2(L)','wt%',100,'t(k)',20.27,'oxid','O2(L)','wt%',100,'t(k)',90.17,'output','mks','end','screen');
% RKT4=CEA('problem','rocket','equilibrium','o/f',5.55157,'case','Example-9', 'p,bar',53.3172,'fac','acat',1.58,'pi/p',10,100,1000,'supar',25:25:75,'reactants','fuel','H2(L)','wt%',100,'t(k)',20.27,'oxid','O2(L)','wt%',100,'t(k)',90.17,'output','mks','end','screen');
% RKT5=CEA('problem','rocket','equilibrium','o/f',5.55157,'case','Example-10','p,bar',53.3172,'fac','ma,kg/m^2',1333.9,'pi/p',10,100,1000,'supar',25:25:75,'reactants','fuel','H2(L)','t(k)',20.27,'oxid','O2(L)','t(k)',90.17,'output','mks','short','end','screen');
% RKT6=CEA('reactants','fuel','CH6N2(L)','rho,g/cc',0.874,'oxid','N2O4(L)','rho,g/cc',1.431,'problem','rocket','equilibirum','case','CEAM-RKT6','p,psi',1000,'pi/p',68.0457,'o/f',2.5,'eql','froz','nfz',2,'supar',5,10,25,50,75,100,150,200, 'only','C','CO','CO2','H','HNO','HNO2','HO2','H2','H2O','H2O2','N','NO','NO2','N2','N2O','O','OH','O2','HCO','NH','CH4','NH2','NH3','H2O(L)','C(gr)','output','mks','end','screen');
% RKT7=CEA('reac','name','SR_AL','AL',1.00,'h,cal/mol',0.,'wt%',18.0,'t(k)',298.15,'name','AM_PERCL','N',1,'H',4,'O',4,'CL',1,'wt%',68.00,'h,cal/mol',-70700.,'t(k)',298.15,'name','Dio_Adi','C',22.,'H',42.,'O',4.,'wt%',2.0,'h,cal/mol',-296000., 't(k)',298.15,'name','HTPB','C',7.332,'H',10.962,'O',0.058,'h,cal/mol',-250.,'wt%',11.00,'t(k)',298.15,'name','HX-752', 'C',14.,'H',16.,'N',2.0,'O',2.,'wt%',0.20,'h,cal/mol',-61000.,'t(k)',298.15,'name','stabilizer','C',38.,'H',66.,'N',2., 'wt%',0.10,'h,cal/mol',-206300.,'t(k)',298.15,'prob','rkt','ar','eql','p,psia',1000,'outp','massf','transport','mks','end','screen');
% RKT8=CEA('reactants','name','HEHN','wt',44.500,'C',2,'H',9,'N',3,'O',4,'h,kcal/mol',-98.000,'rho,g/cc',1.428,'name','HAN', 'wt',41.830,'H',4,'N',2,'O',4,'h,kcal/mol',-95.300,'rho,g/cc',1.090,'name','AN','wt',02.225,'H',4,'N',2,'O',3,'h,kcal/mol',-87.380,'rho,g/cc',1.725,'name','DIPY','wt',00.445,'C',10,'H',8,'N',2,'h,kcal/mol',44.478,'rho,g/cc',1.326,'name','Water','wt',11.000,'H',2,'O',1,'h,kcal/mol',-68.313,'rho,g/cc',1.000,'problem','rocket','frozen', 'p,psia',300,'sup,ae/at',10.0000,'output','massf','transport','end','screen');
% RKT9=CEA('reac','name','AM_PERCL','N',1,'H',4,'O',4,'CL',1,'wt%',68.00,'h,cal/mol',-70700.,'t(k)',298.15,'name','SR_AL', 'AL',1.00,'h,cal/mol',0.,'wt%',18.0,'t(k)',298.15, 'name','Dio_Adi','C',22.,'H',42.,'O',4.,'wt%',2.0,'h,cal/mol',-296000.,'t(k)',298.15, 'name','HTPB','C',7.332,'H',10.962,'O',0.058, 'h,cal/mol',-250.,'wt%',11.00,'t(k)', 298.15,'name','HX-752', 'C',14.,'H',16.,'N',2.0,'O',2.,'wt%',0.20,'h,cal/mol',-61000.,'t(k)',298.15,'name','stabilizer','C',38.,'H',66.,'N',2.,'wt%',0.10,'h,cal/mol',-206300.,'t(k)',298.15, 'prob','rocket','equilibrium','p,psia',1000,'outp','transport','massf','mks','end','screen');
% RKT10=CEA('problem','rocket','equilibrium','case','CEAM-RKT10','p,psi',1000,'ions','pi/p',68.0457,'sub',10,'supar',10,20,100,'reactants','fuel','Li(Cr)','moles',1,'t(k)',298.15,'oxid','F2(L)','moles',0.5556,'t(k)',85.02,'output','mks','end','screen');
% RKT11=CEA('reactants','fuel','N2H4(L)','wt%',80,'t(k)',298.15,'fuel','Be(a)','wt%',20,'t(k)',298.15,'oxid','H2O2(L)','wt%',100,'t(k)',298.15,'prob','rocket','case','CEAM-RKT11','p,psia',3000,'pi/p',3,10,30,300,'equilibrium','%fuel',67,'insert','BeO(L)','output','trace',1e-10,'calories','end','screen');
% RKT12=CEA('reactants','fuel','ADN','H',4','N',4,'O',4,'wt%',65,'h,cal/mol',-148,'t(k)',298.15,'fuel','Water','H',2,'O',1,'wt%',9,'h,cal/mol',-241.83,'t(k)',298.15,'oxid','Ammonia','H',3,'N',1,'wt%',6,'h,cal/mol',-45.94,'t(k)',298.15,'fuel','Methanol','C',1,'H',4,'O',1,'wt%',20,'h,cal/mol',-205.4,'t(k)',298.15,'prob','rocket','case','CEAM-RKT12','p,psia',3000,'pi/p',3,10,30,300,'equilibrium','nfz',2,'only','H2O','CO2','CO','N2','H2','end','screen');
% RKT13=CEA('reactants','fuel','ADN','H',4','N',4,'O',4,'wt%',65,'h,cal/mol',-148,'t(k)',298.15,'fuel','Water','H',2,'O',1,'wt%',9,'h,cal/mol',-285.80,'t(k)',298.15,'oxid','Ammonia','H',3,'N',1,'wt%',6,'h,cal/mol',-71.54,'t(k)',298.15,'fuel','Methanol','C',1,'H',4,'O',1,'wt%',20,'h,cal/mol',-238.6,'t(k)',298.15,'prob','rocket','case','CEAM-RKT13','p,psia',135.,'supar',50, 'equilibrium','nfz',2,'only','H2O','CO2','CO','N2','H2','end','screen');
% RKT14=CEA('prob','rkt','eql','p,psi',1375,'fac','ma,lbm/ft2',323.2,'o/f',6,'reac','fuel','H2(L)','oxid','O2(L)','output','eng','end','screen');
% RKT15=CEA('problem','rocket','equilibrium','o/f',2.72,'p(bar)',256.62,'supsonic(ae/at)',36.4,'reactants','fuel','RP-1','C',1,'H',1.95,'h,J/mol',-24717.7,'wt%',100,'oxid','O2(L)','wt%',100,'h,J/mol',-12979,'output','transport','massf','mks','short','end','screen');
% RKT16=CEA('problem','rocket','equilibrium','o/f',2.72,'p(bar)',256.62,'supsonic(ae/at)',36.4,'reactants','fuel','RP-1','C',1,'H',1.95,'h,J/mol',-24717.7,'wt%',100,'oxid','O2(L)','wt%',100,'output','transport','massf','mks','short','end','screen');
% RKT17=CEA('problem','rocket','equilibrium','o/f',2.72,'p(bar)',256.62,'supsonic(ae/at)',36.4,'reactants','fuel','RP-1','wt%',100,'oxid','O2(L)','wt%',100,'output','transport','mks','short','end','screen');
% RKT18=CEA('problem','rocket','equilibrium','o/f',2.72,'p(bar)',256.62,'supsonic(ae/at)',36.4,'reactants','fuel','RP-1','C',1,'H',1.95,'h,J/mol',-78700,'wt%',100,'oxid','O2(L)','wt%',100,'output','transport','massf','mks','short','end','screen');
% RKT19=CEA('problem','rocket','equilibrium','o/f',2.72,'p(bar)',256.62,'supsonic(ae/at)',36.4,'reactants','fuel','RP-1','C',1,'H',1.95,'h,J/mol',-29310,'wt%',100,'oxid','O2(L)','wt%',100,'output','transport','massf','mks','short','end','screen');
% RKT20=CEA('problem','rocket','equilibrium','o/f',5.55157,'case','CEAM-RKT20','p,bar',53.3172,'fac','mdot,gm/cm^2',13.339,'pi/p',10,100,1000,'supar',25:25:75,'reactants','fuel','H2(L)','t(k)',20.27,'oxid','O2(L)','t(k)',90.17,'output','transport','mks','short','end','screen');
% RKT21=CEA('problem','rocket','equilibrium','o/f',5.55157,'case','CEAM-RKT21','p,bar',53.3172,'fac','ma,kilog/m^2',133.39,'pi/p',10,100,1000,'supar',25:25:75,'reactants','fuel','H2(L)','t(k)',20.27,'oxid','O2(L)','t(k)',90.17,'output','transport','cgs','short','end','screen');
% RKT22=CEA('problem','rocket','equilibrium','p(psi)',200,'subar',3,'supar',5,10,25,50,100,'o/f',3.5,'reactants','fuel','CH4(L)','wt%',100,'t(k)',102,'oxid','O2(L)','wt%',100,'t(k)',100,'output','transport','SI','end','screen');
% RKT23=CEA('problem','rocket','equilibrium','p(atm)',1,'O/F',27.58164,'omit','Air','reactants','fuel','CH4','ox','Air','output','transport','end','screen');
% RKT24=CEA('problem','rocket','frozen','o/f',5,'p(psi)',1500,'supsonic(ae/at)',100,'reactants','fuel','H2(L)','wt%', 100,'oxid','O2(L)','wt%',100,'output','transport','end','screen');
% RKT25=CEA('problem','rocket','frozen','p(psi)',800,'supar',10,50,70,100,'o/f',.2,.5,3,'reactants','fuel','CH4(L)','wt%',100,'oxid','O2(L)','wt%',100,'output','transport','SI','end','screen');
% RKT26=CEA('problem','rocket','equilibrium','o/f',1.5:.25:2.75,'p,psia',400:50:1000,'subar',1.58,'supar',8,12,16,'reactants','fuel','AL(cr)','wt%',40,'t,k',298.15,'fuel','AL2O3','wt%',0,'t,k',298.15,'oxid','O2','wt%',100,'t,k',298.15,'oxid','N2','wt%',0,'t,k',298.15,'fuel','HTPB','wt%',52.32,'t,k',298.15,'C',7.332,'H',11,'O',0.058,'h,cal/mol',-250.0,'fuel','DioctylAdipate','wt%',7.68,'C',22,'H',42,'O',4,'t,k',298.15,'output','calories','massf','end','screen');
% RKT27=CEA('problem','rocket','equilibrium','o/f',2.1:.2:2.9,'p,psia',400:50:950,'supar',8,12,'reactants','fuel','AL(cr)','wt%',40,'t,k',298.15,'fuel','AL2O3','wt%',0,'t,k',298.15,'oxid','O2','wt%',66.66667,'t,k',298.15,'oxid','N2','wt%',33.33333,'t,k',298.15,'fuel','HTPB','wt%',52.32,'t,k',298.15,'C',7.332,'H',10.962,'O',0.058,'h,cal/mol',-250.0,'fuel','DioctylAdipate','wt%',7.68,'C',22,'H',42,'O',4,'t,k',298.15,'h,cal/mol',-296000,'output','calories','massf','transport','end','screen');
% RKT28=CEA('problem','rocket','equilibrium','o/f',1.3,'p,psia',400,'subar',1.58,'supar',8,12,16,'reactants','fuel','AL(cr)','wt%',40,'t,k',298.15,'fuel','AL2O3','wt%',0,'t,k',298.15,'oxid','O2','wt%',66.66667,'t,k',298.15,'oxid','N2','wt%',33.33333,'t,k',298.15,'fuel','HTPB','wt%',52.32,'t,k',298.15,'C',7.332,'H',10.962,'O',0.058,'h,cal/mol',-250.0,'fuel','DioctylAdipate','wt%',7.68,'C',22,'H',42,'O',4,'t,k',298.15,'h,cal/mol',-296000,'output','calories','massf','debug','end','screen');
% RKT29=CEA('problem','rocket','equilibrium','case','active','p(psi)',513.6973,'sub',10.4425,7.2208,5.3736,4.3403,3.4294,2.6256,1.9775,1.5578,1.2885,1.1211,1.0293,'sup',1.0972,1.4354,2.0794,2.7858,3.5417,4.3519,5.2454,6.2223,7.2825,8.4261,8.6382,'reactants','name','AL(cr)','wt',4140,'t(k)',368.15,'name','NH4CLO4(I)','wt',15451,'t(k)',368.15,'name','R-45HTLO','wt',1821.6,'t(k)',368.15,'h,kj/mol',-0.1243,'C',7.337,'H',10.982,'O',0.058,'name','DOA','wt',837.2,'t(k)',368.15,'h,kj/mol',-1303,'C',22,'H',42,'O',4,'output','massf','transport','end','screen');
% RKT30=CEA('problem','rocket','equilibrium','o/f',5.55157,'case','CEAM-RKT20','p,bar',53.3172,'pi/p',10,50,100,500,1000,'subar',3:2:11,'supar',25:25:125,'reactants','fuel','H2(L)','t(k)',20.27,'oxid','O2(L)','t(k)',90.17,'output','transport','mks','short','end','screen');
%
%%%% Shock Examples
% SHK1=CEA('reac','name','H2','moles',2.0,'name','O2','moles',1.0,'prob','sh','inc','fr','eql','u1',2837.1,'t',300.,'p,atm',1.0,'output','transport','end','screen');
% SHK2=CEA('reac','name','H2','moles',0.050,'t(k)',300.0,'name','O2','moles',0.050,'t(k)',300.0,'name','Ar','moles',0.900,'t(k)',300.0,'problem','case','Example-7','p,mmhg',10,20,'shock','u1,m/s',1000,1100,1200,1250,1300,1350,'incd','froz','eql','end','screen');
% SHK3=CEA('reac','name','H2','moles',2.0,'name','O2','moles',1.0,'prob','sh','inc','fr','eql','Mach1',5.2588453,'t',300.,'p,atm',1.0,'output','transport','end','screen');
% SHK4=CEA('reac','name','H2','moles',0.050,'t(k)',300.0,'name','O2','moles',0.050,'t(k)',300.0,'name','Ar','moles',0.900,'t(k)',300.0,'problem','case','CEAM-SHK4','p,mmhg',10,20,'shock','u1,m/s',1100,1200,1250,1300,1350,'incd','refl','froz','eql','end','screen');
% SHK5=CEA('reac','name','H2','moles',0.050,'t(k)',300.0,'name','O2','moles',0.050,'t(k)',300.0,'name','Ar','moles',0.900,'t(k)',300.0,'problem','case','Example-7','p,mmhg',10,20,'shock','u1,m/s',1000,1100,1200,1250,1300,1350,'incd','eql','output','transport','end','screen');
% SHK6=CEA('reac','name','H2','moles',0.050,'t(k)',300.0,'name','O2','moles',0.050,'t(k)',300.0,'name','Ar','moles',0.900,'t(k)',300.0,'problem','case','Example-7','p,mmhg',10,20,'shock','u1,m/s',1000,1100,1200,1250,1300,1350,'incd','frz','output','transport','end','screen');
% SHK7=CEA('problem','shock','inc','ions','eq','ref','u1',2000,'t,k',295,'p,psia',14.7,'react','name','Air','m',1,'t,k',298.15,'omit','Air','output','end','screen');
%
%%%% Detonation Examples
% DET1=CEA('reac','name','H2','moles',2.0,'name','O2','moles',1.0,'prob','det','t,k',300.,'p,atm',1.0,'output','transport','end','screen');
% DET2=CEA('reac','name','C2H4','moles',1.0000,'name','O2','moles',3.0000,'prob','det','t,k',298.,'p,atm',1.0,'output','transport','end','screen');
% DET3=CEA('reac','oxid','O2','wt%',100,'t(k)',298.15,'fuel','H2','wt%',100,'t(k)',298.15,'problem','det','case','Example-6','t(k)',298.15,500,'r,e',1,'p,bar',1,20,'output','trace',1e-5,'calories','transport','end','screen');
