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

function [dV,T,a_h,e_h,theta1,theta2,V_i,V_f] = Trasferimento_orbitale(R_i,R_f,i_i,i_f,U)

if nargin < 2 || nargin == 3
    error('Numero insufficiente di imput')
elseif nargin > 5
    error('Numero eccessivo di input')
elseif nargin>3 && nargin<=5
    U = 398600.4415; % km^3/s^2 
    di = abs(i_f - i_i); % Differenza di inclinazione tra le orbite iniziale e finale
elseif nargin == 2
    U = 398600.4415; % km^3/s^2
    di = 0; % °
end
RE = 6378.1363;
if R_i<=RE
    error('Orbita iniziale troppo bassa')
end
% Trasformazione da deg a rad
d2r = pi/180;
dI = di*d2r; % rad
r2d= 180/pi;
% Velocità dell'orbita iniziale
V_i = sqrt(U/R_i); %km/s

% Velocità dell'orbita finale
V_f = sqrt(U/R_f); %km/s

% Parametri orbita di trasferimento
a_h = (R_i + R_f)/2;  % Semiasse maggiore 
R_h_a = max(R_i,R_f); % raggio di apogeo 
R_h_p = min(R_i,R_f); % raggio di perigeo 
e_h = (R_h_a-R_h_p)/(R_i+R_f); % eccentricità 
i_h=0; %°

% Velocità dell'orbita di trasferimento al perigeo e all'apogeo
V_h_p = sqrt(2*U/R_h_p - U/a_h); %km/s
V_h_a = sqrt(2*U/R_h_p - U/a_h); %km/s

% Calcolo dei dV
if di==0 % Trasferimento tra orbite complanari 
    if R_i>R_f
        dVp = abs(V_f - V_h_p);
        dVa = abs(V_h_a - V_i);
        dV = [dVp, dVa]
    else
        dVp = abs(V_h_p - V_i); %km/s
        dVa = abs(V_f - V_h_a);  %km/s
        dV = [dVp, dVa]; %km/s
    end
elseif di~=0 % Trasferimento tra orbite non complanari
disp('----------------------------------------------------------------------------------------------------------------------------------------------')
disp(' Tipologie di cambio piano al nodo:')
disp(' 1 - Cambio di piano al perigeo')
disp(' 2 - Cambio di piano all apogeo')
disp(' 3 - Cambio di piano ripartito')
disp('----------------------------------------------------------------------------------------------------------------------------------------------')
    cambio = input('Tipologia di cambio piano:');
        switch (cambio)
            case 1
                dVp = sqrt(V_i^2+V_h_p^2-(2*V_i*V_h_p*cos(dI))); % km/s
                dVa = abs(V_f-V_h_a); %km/s
                dV = [dVp, dVa]; %km/s
                theta1=di; %°
                theta2=0; %°
                i_h=i_f; %°
            case 2
                dVp =  abs(V_h_p-V_i); % km/s
                dVa = sqrt(V_f^2+V_h_a^2-(2*V_f*V_h_a*cos(dI))); %km/s
                dV = [dVp, dVa]; %km/s
                theta1=0; %°
                theta2=di; %°
                i_h=i_i;
            case 3
                R1 = sqrt(2*R_f/(R_i+R_f));
                R2 = sqrt(2*R_i/R_f);
                R3 = sqrt(2*R_i/(R_i+R_f));
                R=[R1;R2;R3];
                    f1 = fopen('Rnorm.dat','w');
                    for i=1:3
                        fprintf(f1,'%.4f',R(i));
                        fprintf(f1,'\n');
                    end
                    fprintf(f1,'%.2f',dI);
                    fclose(f1)       
                %rtol=input('Tolleranza per l algoritmo di brent: ');
                rtol=10^(-8);
                [xroot, froot] = brent ('optinc',0,dI,rtol);
                THETA1 = abs(xroot); %rad
                THETA2 = dI - THETA1;       
                dV1 = V_i*sqrt(1+R1^2-2*R1*cos(THETA1));
                dV2 = V_f*sqrt(R2^2+(R1^2)*(R3^2)-2*(R2^2)*R3*cos(THETA2));
                dV = [dV1,dV2];
                theta1 = THETA1*r2d;
                theta2 = THETA2*r2d;
                if i_f>i_i
                    i_h = i_i + theta1;
                else
                    i_h = i_i - theta1;
                end
        end        
end

% Tempo di trasferimento
T = pi*(a_h^3/U)^.5; %sec


% Creazione vettore degli elementi orbitali
oevi(1) = R_i;
oevi(2) = 0.0;
oevi(3) = i_i*d2r;
oevi(4) = 0.0;
oevi(5) = 0.0;

% Determinazione corretta anomalia vera
if (R_f>R_i)
    oevi(6) = 0.0;
else  
    oevi(6) = 180.0 *d2r;
end

[ri, vi] = orb2eci(U, oevi);

oevti(1) = a_h;
oevti(2) = e_h;
oevti(3) = i_h*d2r;
oevti(4) = 0.0;
oevti(5) = 0.0;

if (R_f>R_i)
    oevti(6)=0.0;
else
    oevti(6)=180.0*d2r;
end

[rti,vti] = orb2eci(U,oevti);

oevtf(1) = a_h;
oevtf(2) = e_h;
oevtf(3) = i_h*d2r;
oevtf(4) = 0.0;
oevtf(5) = 0.0;

if (R_f>R_i)   
    oevtf(6) = 180.0*d2r;   
else
    oevtf(6) = 0.0;  
end

[rtf, vtf] = orb2eci(U, oevtf);

oevf(1) = R_f;
oevf(2) = 0.0;
oevf(3) = i_f*d2r;
oevf(4) = 0.0;
oevf(5) = 0.0;

if (R_f>R_i)
    oevf(6) = 180.0*d2r;   
else  
    oevf(6) = 0.0;   
end

[rf, vf] = orb2eci(U, oevf);

% Periodi orbitali
period1 = 2.0 * pi * oevi(1) * sqrt(oevi(1) / U);

period2 = 2.0 * pi * oevti(1) * sqrt(oevti(1) / U);

period3 = 2.0 * pi * oevf(1) * sqrt(oevf(1) / U);

deltat1 = period1 / 300;

simtime1 = -deltat1;

deltat2 = 0.5 * period2 / 300;

simtime2 = -deltat2;

deltat3 = period3 / 300;

simtime3 = -deltat3;

for i = 1:1:301
    
    simtime1 = simtime1 + deltat1;
    
    simtime2 = simtime2 + deltat2;
    
    simtime3 = simtime3 + deltat3;
    
    % Posizioni orbitali iniziali normalizzate
    
    [rwrk, vwrk] = twobody2 (U, simtime1, ri, vi);
    
    rp1_x(i) = rwrk(1) / RE;
    
    rp1_y(i) = rwrk(2) / RE;
    
    rp1_z(i) = rwrk(3) / RE;
    
    
    [rwrk, vwrk] = twobody2 (U, simtime2, rti, vti);
    
    rp2_x(i) = rwrk(1) / RE;
    
    rp2_y(i) = rwrk(2) / RE;
    
    rp2_z(i) = rwrk(3) / RE;
    
    [rwrk, vwrk] = twobody2 (U, simtime3, rf, vf);
    
    rp3_x(i) = rwrk(1) / RE;
    
    rp3_y(i) = rwrk(2) / RE;
    
    rp3_z(i) = rwrk(3) / RE;
    
end

figure(1);

% Creazione vettori degli assi

xaxisx = [1 1.5];
xaxisy = [0 0];
xaxisz = [0 0];

yaxisx = [0 0];
yaxisy = [1 1.5];
yaxisz = [0 0];

zaxisx = [0 0];
zaxisy = [0 0];
zaxisz = [1 1.5];

figure (1);
hold on;

grid on;

% plot earth
TERRA=imread('planisphere.jpg','jpg');
props.FaceColor='texture';
props.EdgeColor='none';
props.FaceLighting='phong';
props.Cdata=TERRA;
[x,y,z]=ellipsoid(0,0,0,1,1,1,24);
surface(-x,-y,-z,props);

% plot sistema di coordinate

plot3(xaxisx, xaxisy, xaxisz, '-g', 'LineWidth', 1);

plot3(yaxisx, yaxisy, yaxisz, '-r', 'LineWidth', 1);

plot3(zaxisx, zaxisy, zaxisz, '-b', 'LineWidth', 1);

% plot orbita iniziale

plot3(rp1_x, rp1_y, rp1_z, '-r', 'LineWidth', 1.5);

plot3(rp1_x(1), rp1_y(1), rp1_z(1), 'ob');

% plot orbita di trasferimento

plot3(rp2_x, rp2_y, rp2_z, '-b', 'LineWidth', 1.5);

plot3(rp2_x(end), rp2_y(end), rp2_z(end), 'ob');

% plot orbita finale

plot3(rp3_x, rp3_y, rp3_z, '-g', 'LineWidth', 1.5);

xlabel('X', 'FontSize', 12);

ylabel('Y', 'FontSize', 12);

zlabel('Z', 'FontSize', 12);

title('Trasferimento di Hohmann', 'FontSize', 16);

axis equal;

view(50, 20);

rotate3d on;

print -depsc -tiff -r300 hohmann1.eps











