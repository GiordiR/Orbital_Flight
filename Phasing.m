function [dV_ph,a_ph,e_ph,T_ph] = Phasing(R,Dtheta,i,U)
% Funzione che implementa una manovra di phasing per lo spostamento lungo
% l'orbita
%
% Inputs: R      - Raggio dell'orbita circolare [km]
%         Dtheta - Angolo tra le due posizioni dell'orbita [°]
%                 (negativo se antiorario, positivo se orario)
%         i      - Inclinazione orbita iniziale (°)
%         U      - Costante gravitazionale planetaria
%
% Outputs: dV_ph - Cambio di velocità necessario [km/s]
%          a_ph  - Semiasse maggiore orbita di phasing [km]
%          e_ph  - Eccentricità orbita di phasing 
%          T_ph  - Periodo dell'orbita phasing [s]
%

if nargin<2
    error('Numero insufficiente di input')
elseif nargin>4
    error('Numero eccessivo di input')
elseif nargin==2 || nargin==3
    U=398600.4415; 
end
    
d2r=pi/180;
RE=6378.1363;
% Periodo dell'orbita iniziale
T=(2*pi*R^(3/2))/sqrt(U);

% Differenza di tempo tra i due punti dell'orbita 
t=abs((T*Dtheta)/360);

% Periodo dell'rotbita di phasing
if Dtheta>0
    T_ph=T-t;
else
    T_ph=T+t;
end

% Semiasse maggiore dell'orbita di phasing
a_ph=((T_ph*sqrt(U))/(2*pi))^(2/3);

% Eccentricità dell'orbita di phasing
if Dtheta>0
    e_ph=(R/a_ph)-1;
    rp=a_ph*(1-e_ph);
    ra=R;
else
    e_ph=1-(R/a_ph);
    rp=R;
    ra=a_ph*(1+e_ph);
end

% Velocità orbita iniziale
V=sqrt(U/R);

% Velocità orbita di perigeo e apogeo dell'orbita di phasing
if Dtheta>0
    va=sqrt(U*(2/ra-1/a_ph));
else 
    vp=sqrt(U*(2/rp-1/a_ph));
end

% Delta velocità manovra
if Dtheta>0
    dV_ph=2*abs(va-V);
else
    dV_ph=2*abs(vp-V);
end


% Creazione vettore degli elementi orbitali
oevi(1)=R;
oevi(2)=0.0;
oevi(3)=i*d2r;
oevi(4)=0.0;
oevi(5)=0.0;
oevi(6)=0.0;

[ri,vi]=orb2eci(U,oevi);

oevph(1)=a_ph;
oevph(2)=e_ph;
oevph(3)=i*d2r;
oevph(4)=0.0;
oevph(5)=0.0;
oevph(6)=0.0;

[rph,vph]=orb2eci(U,oevph);

% Periodi orbitali
period1 = 2.0 * pi * oevi(1) * sqrt(oevi(1) / U);

period2 = 2.0 * pi * oevph(1) * sqrt(oevph(1) / U);

deltat1 = period1 / 300;

simtime1 = -deltat1;

deltat2 = period2 / 300;

simtime2 = -deltat2;

for i = 1:1:301
    
    simtime1 = simtime1 + deltat1;
    
    simtime2 = simtime2 + deltat2;
        
    % Posizioni orbitali iniziali normalizzate
    
    [rwrk, vwrk] = twobody2 (U, simtime1, ri, vi);
    
    rp1_x(i) = rwrk(1) / RE;
    
    rp1_y(i) = rwrk(2) / RE;
    
    rp1_z(i) = rwrk(3) / RE;
    
    
    [rwrk, vwrk] = twobody2 (U, simtime2, rph, vph);
    
    rp2_x(i) = rwrk(1) / RE;
    
    rp2_y(i) = rwrk(2) / RE;
    
    rp2_z(i) = rwrk(3) / RE;
end

figure(2);

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

figure (2);
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

% plot orbita di phasing

plot3(rp2_x, rp2_y, rp2_z, '-b', 'LineWidth', 1.5);

plot3(rp2_x(end), rp2_y(end), rp2_z(end), 'ob');


xlabel('X', 'FontSize', 12);

ylabel('Y', 'FontSize', 12);

zlabel('Z', 'FontSize', 12);

title('Manovra di phasing', 'FontSize', 16);

axis equal;

view(50, 20);

rotate3d on;



    


end

