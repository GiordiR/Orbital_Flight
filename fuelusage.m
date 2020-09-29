function dM = fuelusage(M,Isp,g0,dV)

if nargin < 2
    error('Numero insufficiente di imput')
elseif nargin > 4
    error('Numero eccessivo di input')
elseif nargin == 3
    if ischar(g0)
        Isp = 350; % s
    elseif ischar(Isp)
        g0=0.009822; % km/s^2
    end
elseif nargin == 3
    Isp = 350; % s
    g0 = 0.009822; % km/s^2
end

dM = M - M*exp(-dV/(Isp*g0));

end

