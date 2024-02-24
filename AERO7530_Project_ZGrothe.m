% Zane Grothe
% AERO 7530
% Computer Project
% 12/10/21

clear all
close all
clc

% To run this script two .m function files and seven .txt text files are
% needed:
% CompEQs1.m
% CompEQs2.m
% JANAF_H.txt
% JANAF_H2.txt
% JANAF_H2O.txt
% JANAF_N2.txt
% JANAF_O.txt
% JANAF_O2.txt
% JANAF_OH.txt

% The seven .txt files are from the JANAF website. They are manually
% imported using the MATLAB Import Wizard (uiimport), because some of the
% import settings need to be adjusted.

% After running the script the Import Wizard window will appear for each
% .txt file. In the top right hand corner of the Import Wizard window is
% the "Number of text header lines" option. If this number is 3, just hit
% enter twice to import.
% The files for H, H2O, O, and OH will have 3. They will be imported first.
% If the number is 2, bump it up to 3 and hit enter twice.
% The files for H2, O2, and N2 will have 2.

J_H=uiimport('JANAF_H.txt');
JH=J_H.data;
J_H2O=uiimport('JANAF_H2O.txt');
JH2O=J_H2O.data;
J_O=uiimport('JANAF_O.txt');
JO=J_O.data;
J_OH=uiimport('JANAF_OH.txt');
JOH=J_OH.data;
J_H2=uiimport('JANAF_H2.txt');
JH2=J_H2.data;
J_O2=uiimport('JANAF_O2.txt');
JO2=J_O2.data;
J_N2=uiimport('JANAF_N2.txt');
JN2=J_N2.data;


% Choose values for Nitrogen content in oxidizer (xR_N2)----------
xR_N2=0;
while xR_N2<1
    
    n_ox=1/2/(1-xR_N2); % moles of oxidizer
    n_N2=n_ox*xR_N2; % moles of N2 in products
    
    % Guess the product tempurature (TempGuess)----------
    TempGuess=2500; % Kelvin
    R=8.314; % ideal gas constant (J/mol/K)
    
    % Initiate enthalpy contribution values for loop
    deltafH=0;
    deltaHsGuess=1;
    
    
    % Compute Kp values using Gibbs free energy of formation----------
    
    while deltafH<deltaHsGuess
        clc
        
        % Find Gibbs free energy of formation at guessed temperature
        deltafG_H=interp1(JH(:,1),JH(:,7),TempGuess,'spline');
        deltafG_H2O=interp1(JH2O(:,1),JH2O(:,7),TempGuess,'spline');
        deltafG_OH=interp1(JOH(:,1),JOH(:,7),TempGuess,'spline');
        deltafG_O=interp1(JO(:,1),JO(:,7),TempGuess,'spline');
        
        % Find Kp
        KpO=exp(-deltafG_O*1000/R/TempGuess);
        KpH=exp(-deltafG_H*1000/R/TempGuess);
        KpH2O=exp(-deltafG_H2O*1000/R/TempGuess);
        KpOH=exp(-deltafG_OH*1000/R/TempGuess);
        
        
        % Solve for mole fractions----------
        if xR_N2==0
            
            fun2=@(x)CompEQs2(x,n_ox,xR_N2,KpO,KpH,KpH2O,KpOH);
            x0=[.1,.1,.1,.1,.1,.1];
            x=fsolve(fun2,x0);
            
            X_H=x(1,1);
            X_H2=x(1,2);
            X_H2O=x(1,3);
            X_OH=x(1,4);
            X_O=x(1,5);
            X_O2=x(1,6);
            
            
            % Compute ntotal----------
            
            ntotal=2/(2*X_H2O+2*X_H2+X_H+X_OH);
            
            % Compute number of product moles
            n_H2O=X_H2O*ntotal;
            n_H2=X_H2*ntotal;
            n_H=X_H*ntotal;
            n_OH=X_OH*ntotal;
            n_O2=X_O2*ntotal;
            n_O=X_O*ntotal;
            
        else
            fun1=@(x)CompEQs1(x,n_ox,xR_N2,KpO,KpH,KpH2O,KpOH);
            x0=[.1,.1,.1,.1,.1,.1,.1];
            x=fsolve(fun1,x0);
            
            X_H=x(1,1);
            X_H2=x(1,2);
            X_H2O=x(1,3);
            X_OH=x(1,4);
            X_O=x(1,5);
            X_O2=x(1,6);
            X_N2=x(1,7);
            
            
            % Compute ntotal----------
            
            ntotal=n_N2/X_N2;
            
            % Compute number of product moles
            n_H2O=X_H2O*ntotal;
            n_H2=X_H2*ntotal;
            n_H=X_H*ntotal;
            n_OH=X_OH*ntotal;
            n_O2=X_O2*ntotal;
            n_O=X_O*ntotal;
            
        end
        
        % Compute product temperature----------
        
        % Find enthalpy of formations at 298.15K
        TempRef=298.15;
        deltafH_H=interp1(JH(:,1),JH(:,6),TempRef);
        deltafH_H2O=interp1(JH2O(:,1),JH2O(:,6),TempRef);
        deltafH_OH=interp1(JOH(:,1),JOH(:,6),TempRef);
        deltafH_O=interp1(JO(:,1),JO(:,6),TempRef);
        
        % Find enthalpy contribution from heats of formation
        deltafH=(n_H2O*deltafH_H2O+n_H*deltafH_H+n_OH*deltafH_OH ...
               +n_O*deltafH_O);
    
        % Find sensible enthalpy values at the guessed temperature
        deltaHs_H=interp1(JH(:,1),JH(:,5),TempGuess,'spline');
        deltaHs_H2O=interp1(JH2O(:,1),JH2O(:,5),TempGuess,'spline');
        deltaHs_OH=interp1(JOH(:,1),JOH(:,5),TempGuess,'spline');
        deltaHs_O=interp1(JO(:,1),JO(:,5),TempGuess,'spline');
        deltaHs_H2=interp1(JH2(:,1),JH2(:,5),TempGuess,'spline');
        deltaHs_O2=interp1(JO2(:,1),JO2(:,5),TempGuess,'spline');
        deltaHs_N2=interp1(JN2(:,1),JN2(:,5),TempGuess,'spline');
        
        % Find sensible enthalpy contribution from guessed temperature
        if xR_N2==0
            deltaHsGuess=n_H*deltaHs_H+n_H2O*deltaHs_H2O+n_OH*deltaHs_OH...
                        +n_O*deltaHs_O+n_H2*deltaHs_H2+n_O2*deltaHs_O2;
        else
            deltaHsGuess=n_H*deltaHs_H+n_H2O*deltaHs_H2O+n_OH*deltaHs_OH...
                        +n_O*deltaHs_O+n_H2*deltaHs_H2+n_O2*deltaHs_O2...
                        +n_N2*deltaHs_N2;
        end
        % repeat loop if sensible enthalpy contribution and enthalpy of
        % formation contribution do not compare
        TempGuess=TempGuess+10; % re-guess in steps of 10 degrees
    end
    clc
    % Collect data during loops
    if xR_N2==0
        AFT=TempGuess;
        NC=xR_N2;
        X_Hm=X_H;
        X_H2m=X_H2;
        X_H2Om=X_H2O;
        X_OHm=X_OH;
        X_Om=X_O;
        X_O2m=X_O2;
        X_N2m=0;
    else
        AFT=[AFT,TempGuess];
        NC=[NC,xR_N2];
        X_Hm=[X_Hm,X_H];
        X_H2m=[X_H2m,X_H2];
        X_H2Om=[X_H2Om,X_H2O];
        X_OHm=[X_OHm,X_OH];
        X_Om=[X_Om,X_O];
        X_O2m=[X_O2m,X_O2];
        X_N2m=[X_N2m,X_N2];
    end
    xR_N2=xR_N2+.2; % next concentration of N2 in reactants
end
AFT % show results

% Plot results----------

% Adiabatic Flame Temperature
NCx=linspace(0,1);
AFTy=interp1(NC,AFT,NCx,'spline');
figure(1)
plot(NC,AFT,'ro')
hold on
plot(NCx,AFTy,'b')
hold off
xlabel('N2 Concentration in Oxidizer')
ylabel('Adiabatic Flame Temp (K)')
title('Adiabatic Flame Temperature vs. Nitrogen Concentration in Oxidizer')

% Monatomic Hydrogen H
X_Hmy=interp1(NC,X_Hm,NCx,'spline');
figure(2)
plot(NC,X_Hm,'ro')
hold on
plot(NCx,X_Hmy,'b')
hold off
xlabel('N2 Concentration in Oxidizer')
ylabel('H Concentration in Products')
title('Monatomic Hydrogen Concentration in Products vs. Nitrogen Concentration in Oxidizer')

% Hydrogen H2
X_H2my=interp1(NC,X_H2m,NCx,'spline');
figure(3)
plot(NC,X_H2m,'ro')
hold on
plot(NCx,X_H2my,'b')
hold off
xlabel('N2 Concentration in Oxidizer')
ylabel('H2 Concentration in Products')
title('Hydrogen Concentration in Products vs. Nitrogen Concentration in Oxidizer')

% Water H2O
X_H2Omy=interp1(NC,X_H2Om,NCx,'spline');
figure(5)
plot(NC,X_H2Om,'ro')
hold on
plot(NCx,X_H2Omy,'b')
hold off
xlabel('N2 Concentration in Oxidizer')
ylabel('H2O Concentration in Products')
title('Water Concentration in Products vs. Nitrogen Concentration in Oxidizer')

% Hydroxide OH
X_OHmy=interp1(NC,X_OHm,NCx,'spline');
figure(6)
plot(NC,X_OHm,'ro')
hold on
plot(NCx,X_OHmy,'b')
hold off
xlabel('N2 Concentration in Oxidizer')
ylabel('OH Concentration in Products')
title('Hydroxide Concentration in Products vs. Nitrogen Concentration in Oxidizer')

% Monatomic Oxygen O
X_Omy=interp1(NC,X_Om,NCx,'spline');
figure(7)
plot(NC,X_Om,'ro')
hold on
plot(NCx,X_Omy,'b')
hold off
xlabel('N2 Concentration in Oxidizer')
ylabel('O Concentration in Products')
title('Monatomic Oxygen Concentration in Products vs. Nitrogen Concentration in Oxidizer')

% Oxygen O2
X_O2my=interp1(NC,X_O2m,NCx,'spline');
figure(8)
plot(NC,X_O2m,'ro')
hold on
plot(NCx,X_O2my,'b')
hold off
xlabel('N2 Concentration in Oxidizer')
ylabel('O2 Concentration in Products')
title('Oxygen Concentration in Products vs. Nitrogen Concentration in Oxidizer')

% Nitrogen N2
X_N2my=interp1(NC,X_N2m,NCx,'spline');
figure(9)
plot(NC,X_N2m,'ro')
hold on
plot(NCx,X_N2my,'b')
hold off
xlabel('N2 Concentration in Oxidizer')
ylabel('N2 Concentration in Products')
title('Nitrogen Concentration in Products vs. Nitrogen Concentration in Oxidizer')


