%%  Determine Moho temperature from qs measurements, then compare to Schutt et al., 2018 Moho values
%   Updated on 8/21/2019
%   Created by Michael Berry, with significant contributions by Tony Lowry, geotherm adapted from XXXXX
clear all
close all

data        = load('data_file.dat');
len         = length(data(:,1));         %length of the vector
temp_obs    = data(:,4)+273.15;          %Schutt, 2018 Pn estimated Moho values
temp_pred   = zeros(len,1);              %empty vector for inputting modelled values.
res_map     = zeros(len,3);              %Array for map values
res_map(:,1)= data(:,1);                 %Long first
res_map(:,2)= data(:,2);                 %lat 2nd

%%Input values for each iteration
%Radiogenic heating values.
A0          = 1.11e-06;  %surface heat production (using a best fit)
lrad        = 17e3;      %radiogentic length scale [m]

%Conductive parameters
a           = 0.28;      %conductive parameter A
b           = 3.16e-4;   %conductive parameter B

%% Loop through the data for all lat long points
for j=1:len 
  %A0 = data(j,9) * 10^-6;
  H   = (data(j,5)*1000)+2e3;       %H     = Thickness of the crust [m] (+2 km to account for the seismic sampling zone)  
  Tr  = data(j,8)+273.15;           %Tr    = mantle potential temperature [K]
  Ts  = data(j,6)+273.15;           %Ts    = surface temperature [K]
  qs  = data(j,3)*1e-3;             %qs    = surface heat flow [w/m3]
%%
  rdl  = lrad;  % reset the radiogenic length 
  tA0  = A0;    % reset the radiogenic heat production
  t3   = 1e10;  % Dummy value for first calculation:
  lcon = -1;    % Dummy value for first calculation:
        while (lcon <= 0) || (lcon >= 4e5)
          qr = tA0*rdl*(1d0-exp(-H/rdl));
%% lcon from Qs prediction:
          lcon = (2d0/sqrt(pi))*(((Tr-Ts)/(a+b*Ts))-rdl*qr)/(qs-qr);
          t1   = (Ts+a*tA0*rdl*rdl*(1d0-exp(-H/rdl)-erf(H/lcon))+(Tr-Ts)*erf(H/lcon))...
                 /(1d0-b*tA0*rdl*rdl*(1d0-exp(-H/rdl)-erf(H/lcon)));
          
          t2   = (Ts+a*tA0*rdl*rdl*(1d0-exp(-(H-1d1)/rdl)-erf((H-1d1)/lcon))+(Tr-Ts)*erf((H-1d1)/lcon))...
                 /(1d0-b*tA0*rdl*rdl*(1d0-exp(-(H-1d1)/rdl)-erf((H-1d1)/lcon)));
          qm   = (t1-t2)/(1d1*(a+b*t1));
%% Reduce the crustal heat production; increase the heat flow if lcon is less than 0 or lcon > 400 km:
          if (lcon <= 0) || (lcon >= 4e5)
            if qs < 0.04
              qs = qs+0.001;
            else
              tA0 = .99*tA0;
              rdl = .95*rdl;
            end
          end
%% Reduce the crustal heat production; increase the heat flow if qm is less than 18 mW/m^2:          
          if (lcon >= 0) && (lcon <= 4e5) && (qm < .018) && (qs-qr < .025)
            t1 = .018-qm;
            if (t1 < t3)
              lcon = -1;
              rdl  = .95*rdl;
            elseif (t1 > t3)
              lkk  = lcon; % Temporary storage until we're out  of the loop
              lcon = 0.9999;
            end
            t3=t1;
          end
        end
        if (lcon == 0.9999)
          lcon = lkk;
        end
%% Defined MOHO Temp using EQN for crustal temp     
  temp_pred(j) = (Ts+a*tA0*rdl^2*(1-exp(-H/rdl)-erf(H/lcon))+(Tr-Ts)*erf(H/lcon))...
                 /(1-b*tA0*rdl^2*(1-exp(-H/rdl)-erf(H/lcon)));
  
  residual    = (temp_obs(j))-(temp_pred(j)); %observed minus predicited
  res_map(j,3)= residual;  
end 
%Write to ASCII file format: Long Lat Residual Temp [K] Pn Temp [K] Qs predicted temp [K]
dlmwrite('outfile.dat',[res_map temp_obs temp_pred],' ');


