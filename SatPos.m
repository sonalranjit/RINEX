% Computing Satellite Position (ECEF)
function satcoord = SatPos(prn_sat,t,Ndata)

for i = 1:length(Ndata)
    IODE(i)     = Ndata(i).IODE;
    Crs(i)      = Ndata(i).Crs;
    Deltan(i)   = Ndata(i).Deltan;
    M0(i)       = Ndata(i).M0;
    Cuc(i)      = Ndata(i).Cuc;
    Ecc(i)      = Ndata(i).Ecc;
    Cus(i)      = Ndata(i).Cus;
    SqrtA(i)    = Ndata(i).SqrtA;
    ToE(i)      = Ndata(i).ToE;
    Cic(i)      = Ndata(i).Cic;
    Cis(i)      = Ndata(i).Cis;
    i0(i)       = Ndata(i).i0;
    Crc(i)      = Ndata(i).Crc;
    omega(i)    = Ndata(i).omega;
    OMEGA(i)    = Ndata(i).OMEGA;  
    OMEGADOT(i) = Ndata(i).OMEGADOT;
    IDOT(i)     = Ndata(i).IDOT;
end

IODE     = IODE([prn_sat])';
Crs      = Crs([prn_sat])';
Deltan   = Deltan([prn_sat])';
M0       = M0([prn_sat])';
Cuc      = Cuc([prn_sat])';
Ecc      = Ecc([prn_sat])';
Cus      = Cus([prn_sat])';
SqrtA    = SqrtA([prn_sat])';
ToE      = ToE([prn_sat])';
Cic      = Cic([prn_sat])';
Cis      = Cis([prn_sat])';
i0       = i0([prn_sat])';
Crc      = Crc([prn_sat])';
omega    = omega([prn_sat])';
OMEGA    = OMEGA([prn_sat])';
OMEGADOT = OMEGADOT([prn_sat])';
IDOT     = IDOT([prn_sat])';

mu = 3.986005e14; %value of earth's universal gravitational parameter
OMEGAe = 7.2921151467e-5; %values of earth's rotation rate
A = SqrtA.*SqrtA; %Semimajor axis
n0 = sqrt(mu./A.^3); %Computed mean motion
n = n0 + Deltan; %Corrected mean motion



tk = check_t(t-ToE);
Mk = M0 + n.*tk; %Mean Anomoly
Mk = rem(Mk+2*pi,2*pi);
Ek = Mk;
for i = 1:10; %Kepler's Equation for Eccentric Anomaly
    Ekold = Ek;
    Ek = Mk+Ecc.*sin(Ek);
    dE = rem(Ek-Ekold,2*pi);
    if abs(dE) < 1e-12;
        break
    end
end
Ek = rem(Ek+2*pi,2*pi);
Vk = atan2(sqrt(1-Ecc.^2).*sin(Ek),(cos(Ek)-Ecc)); %True Anomaly

Phik = Vk + omega; %Argument of Latitude
Phik = rem(Phik,2*pi);
dukPhik = Cus.*sin(2.*Phik)+Cuc.*cos(2.*Phik); %Argument of Latitude Correction
drk = Crs.*sin(2.*Phik)+Crc.*cos(2.*Phik); %Radius Correction
dik = Cis.*sin(2*Phik)+Cic.*cos(2*Phik); %Inclination Correction

uk = Phik + dukPhik; %Corrected Argument of Latitude
rk = A.*(1-Ecc.*cos(Ek))+drk; %Corrected Radius
ik = i0+dik+IDOT.*tk; %Corrected Inclination

xo = rk.*cos(uk); %x position in orbital plane
yo = rk.*sin(uk); %y position in orbital plane
%Corrected longitude of ascending node
OMEGAk = OMEGA+(OMEGADOT - OMEGAe).*tk-OMEGAe.*ToE; 
OMEGAk = rem(OMEGAk+2*pi,2*pi);
x = xo.*cos(OMEGAk)-yo.*cos(ik).*sin(OMEGAk); %x coordinate
y = xo.*sin(OMEGAk)+yo.*cos(ik).*cos(OMEGAk); %y coordinate
z = yo.*sin(ik); %z coordinate

satcoord = [x y z];
end







