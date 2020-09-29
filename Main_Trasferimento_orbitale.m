%% Main_Trasferimento_orbitale.m
clear all
clear
clc

R_i = input('Raggio dell orbita iniziale [km]:');
R_f = input('Raggio dell orbita finale [km]:');
i_i = input('Inclinazione dell orbita iniziale [°]:');
i_f = input('Inclinazione dell orbita finale [°]:');

M = input('Massa del satellite [kg]:');
Isp = input('Impulso specifico [s]:');
g0=0.009822;

default =  398600.4415;
U = input(' Costante gravitazionale planetaria [km^3/s^2] (Terra=default): ');

% Chiamata alla funzione trasferimento_orbitale
[dV,T,a_h,e_h,theta1,theta2,V_i,V_f] = Trasferimento_orbitale(R_i,R_f,i_i,i_f,U);
dVt = sum(dV);

% Chiamata  alla funzione fuelusage
dM = fuelusage(M,Isp,g0,dVt);

% Stampa dei risultati
fprintf(' Raggio dell orbita iniziale: %.4f km \n',R_i)
fprintf(' Velocità dell orbita iniziale: %.4f km/s \n',V_i)
if isnumeric(i_i)
    fprintf(' Inclinazione dell orbita iniziale: %.4f ° \n',i_i)
end
fprintf('\n')
fprintf(' Raggio dell orbita finale: %.4f km \n',R_f)
fprintf(' Velocità dell orbita finale: %.4f km/s\n',V_f)
if isnumeric(i_f)
    fprintf(' Inclinazione dell orbita finale: %.4f ° \n',i_f)
end
fprintf('\n')
fprintf(' Semiasse maggiore dell orbita di trasferimento: %.4f km \n',a_h)
fprintf(' Eccentricità dell orbita di trasferimento: %.4f \n', e_h)
fprintf(' Cambio di velocità necessario al trasferimento: %.4f km/s \n',dVt)
fprintf('\n')
if isnumeric(theta1) && isnumeric(theta2) && theta1~=0
    fprintf(' Angolo di cambio piano al primo impulso: %.4f ° \n',theta1)
    fprintf(' Cambio di velocità necessario nel primo impulso: %.4f km/s \n',dV(1))
    fprintf(' Angolo di cambio piano al secondo impulso: %.4f ° \n',theta2)
    fprintf(' Cambio di velocità necessario nel secondo impulso: %.4f km/s \n',dV(2))
end
fprintf('\n')
fprintf(' Tempo di trasferimento: %.4f s \n',T)
fprintf(' Consumo di carburante: %.4f kg \n',dM)

disp('-----------------------------------------------------------------------------');
PH=input(' Eseguire una manovra di phasing) [s/n]: ','s');
if PH=='s'
    disp('Inserire l angolo tra le due posizioni dell orbita (negativo se la posizione finale precede quella iniziale)');
    Dtheta=input('Differenza di posizione [°]: '); 
    [dV_ph,a_ph,e_ph,T_ph] = Phasing(R_f,Dtheta,i_f,U);
    fprintf(' Semiasse maggiore dell orbita di phasing: %.4f km \n',a_ph)
    fprintf(' Eccentricità dell orbita di phasing: %.4f \n',e_ph)
    fprintf(' Periodo dell orbita di phasing: %.4f \n',T_ph)
    fprintf(' Cambio di velocità manovra di phasing: %.4f km/s \n',dV_ph)
end
