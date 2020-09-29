%% Trasferimento_orbitale.m
% Funzione per il calcolo dei parametri di un trasferimento alla Hohmann
%
%
% Inputs:  o R_i    - Raggio dell'orbita iniziale in km
%          o R_f    - Raggio dell'orbita finale in km
%          o U      - Costante grafitazionale del pianeta (km^3/s^2).
%                     Per la terra (default) è 398600.4415 km^3/s^2.
%          o i_i    - Inclinazione dell'orbita iniziale
%          o i_f    - Inclinazione dell'orbita finale
%
% Outputs: o dV   - Una matrice 1x2 dei vettori di cambio di velocità necessari allo "sparo" iniziale e finale in km/s
%          o T    - Tempo del trasferimento
%          o a_h  - Semiasse maggiore dell'orbita di trasferimento 
%

function [dV,T,a_h] = Trasferimento_orbitale(R_i,R_f,i_i,i_f,U)

if nargin < 2
    error('Numero insufficiente di imput')
elseif nargin > 5
    error('Numero eccessivo di input')
elseif nargin == 4 
    U = 398600.4415; % km^3/s^2 
    di = i_f - i_i; % Differenza di inclinazione tra le orbite iniziale e finale
elseif nargin == 2
    U = 398600.4415; % km^3/s^2
    di = 0; % °
else
    di = i_f - i_i; % Differenza di inclinazione tra le orbite iniziale e finale
end

if R_i<=6378
    error('Orbita iniziale troppo bassa')
end



    
% Velocità dell'orbita iniziale
V_i = (U/R_i)^.5; %km/s

% Velocità dell'orbita finale
V_f = (U/R_f)^.5; %km/s

a_h = (R_i + R_f)/2;  % Semiasse maggiore (km)

% Velocità dell'orbita di trasferimento al perigeo (a) e all'apogeo (b)
V_h_a = (2*U/R_i - U/a_h)^.5; %km/s
V_h_b = (2*U/R_f - U/a_h)^.5; %km/s

if di==0
% Caso di complanarità delle due orbite    
% Delta velocità necessari 
dVa = V_h_a - V_i; %km/s
dVb = V_f - V_h_b;  %km/s

dV = [dVa, dVb]; %km/s

elseif di~=0
disp('----------------------------------------------------------------------------------------------------------------------------------------------')
disp(' Tipologie di cambio piano')
disp(' 1 - Cambio di piano al perigeo')
disp(' 2 - Cambio di piano all apogeo')
disp(' 3 - Cambio di piano ripartito')
disp('----------------------------------------------------------------------------------------------------------------------------------------------')
cambio = input('Tipologia di cambio piano:');
    switch (cambio)
        case 1
        dVa = sqrt(V_i^2+V_h_a^2-(2*V_i*V_h_a*cos(di))); % km/s
        dVb = abs(V_f-V_h_b); %km/s
        dV = [dVa, dVb]; %km/s
        
        case 2
        dVa =  abs(V_h_a-V_i); % km/s
        dVb = sqrt(V_f^2+V_h_b^2-(2*V_f*V_h_b*cos(di))); %km/s
        dV = [dVa, dVb]; %km/s
        
        case 3
        [dV, dI] = dVdI(R_i,R_f,di,U,Tol);  
        
    end        
end

% Tempo di trasferimento
T = pi*(a_h^3/U)^.5; %sec

          
          

% Plot del trasferimento 
R_E = 6378.14;
theta = [0:0.001:2*pi]; %campo di valori dell'anomalia vera (da 0 a 360°) 

% Plottaggio del Pianeta Terra e applicazione delle textures sulla superficie della sfera 
TERRA = imread('planisphere.jpg','jpg'); 
props.FaceColor='texture'; 
props.EdgeColor='none'; 
props.FaceLighting='phong'; 
props.Cdata = TERRA; 
Center = [0; 0; 0]; 
[XX, YY, ZZ] = ellipsoid(Center(1),Center(2),Center(3),R_E,R_E,R_E,30); 
surface(-XX, -YY, -ZZ,props); 
hold on 

% Plottaggio dell'orbita del satellite nell'orbita iniziale
X_ORBIT_i = cos(theta) .* R_i; 
Y_ORBIT_i = sin(theta) .* R_i; 
Z_ORBIT_i = zeros(1,length(theta)) .* R_i; 
plot3(X_ORBIT_i,Y_ORBIT_i,Z_ORBIT_i,'LineWidth',2.0)
hold on 

% Parametri orbita di traferimento
theta_h = [0:0.001:pi];
e_h = (R_f - R_i)/(R_f + R_i);
R_h = (a_h*(1- e_h^2))./(1+e_h*cos(theta_h));

% Plottaggio dell'orbita di trasferimento 
X_ORBIT_h = cos(theta_h) .* R_h; 
Y_ORBIT_h = sin(theta_h) .* R_h; 
Z_ORBIT_h = theta_h .* 0.0 .* R_h; 
plot3(X_ORBIT_h,Y_ORBIT_h,Z_ORBIT_h,'g-','LineWidth',2.0) 

% Plottaggio dell'orbita finale
X_ORBIT_f = cos(theta) .* R_f; 
Y_ORBIT_f = sin(theta) .* R_f; 
Z_ORBIT_f = zeros(1,length(theta)) .* R_f;
plot3(X_ORBIT_f,Y_ORBIT_f,Z_ORBIT_f,'LineWidth',2.0)
hold on

grid on % Definizione degli assi su cui fare riferimento nell'immagine dell'orbita 
axis equal 
xlabel('X [km]') 
ylabel('Y [km]') 
zlabel('Z [km]') 
title('3D orbit') 
view( [1 1 1] )