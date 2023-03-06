

close all; close all; clc;

global mu

% 1 au = 149597870700 m
% 1 day = 86164.09054 s
au_km = 1.49597870700e8;% km per au
d_s = 86400; % seconds per day

mu = 1.327*10^11; %sun mu km3/s2

% all in km, km/s
r1I = au_km * [3.515868886595499e-2,-3.162046390773074, 4.493983111703389];
v1I = au_km / d_s * [-2.317577766980901e-3,9.843360903693031e-3,1.541856855538041e-2];

r2I = au_km * [7.249472033259724, 14.61063037906177, 14.24274452216359];
v2I = au_km / d_s * [-8.241709369476881e-3,-1.156219024581502e-2,1.317135977481448e-2];

rE = au_km * [-1.796136509111975e-1,9.667949206859814e-1,-3.668681017942158e-5];
vE = au_km / d_s * [-1.720038360888334e-2,-3.211186197806460e-3,7.927736735960840e-7];

[a,e,i,Omega,omega,f] = orbital_elements_from_R_V(r1I,v1I);
fprintf('1I orbital elements [a (au),e,i,Omega,omega,f (all deg)]:\n')
disp([a/au_km,e,i*180/pi,Omega*180/pi,omega*180/pi,f*180/pi])

[a,e,i,Omega,omega,f] = orbital_elements_from_R_V(r2I,v2I);
fprintf('2I orbital elements [a (au),e,i,Omega,omega,f (all deg)]:\n')
disp([a/au_km,e,i*180/pi,Omega*180/pi,omega*180/pi,f*180/pi])

[a,e,i,Omega,omega,f] = orbital_elements_from_R_V(rE,vE);
fprintf('Earth orbital elements [a (au),e,i,Omega,omega,f (all deg]:\n')
disp([a/au_km,e,i*180/pi,Omega*180/pi,omega*180/pi,f*180/pi])

% 3
n = 100;
launch = linspace(0,364,n); % launch window (d),2017/1/1-2017/12/31
arrive = linspace(212,760,n); % arrive window (d), 2017/8/1-2019/1/31
[launch, arrive] = meshgrid(launch,arrive);
Delta_V_tot = zeros(n,n);
for i=1:n
    for j = 1:n
        [R1, V1] = State_R0_V0_t(rE,vE,d_s.*launch(i,j));
        [R2, V2] = State_R0_V0_t(r1I,v1I,d_s.*arrive(i,j));
        [V1m, V2m] = LambertU_V1_V2_from_R1_R2_t(R1,R2,(d_s.*arrive(i,j)-d_s.*launch(i,j)));
        Delta_V_tot(i,j) = norm(V1m-V1) + norm(V2m-V2);
        %capping delta v
        if Delta_V_tot(i,j)>60
            Delta_V_tot(i,j)=NaN;
        end
    end
end
Delta_V_tot;
contour(launch,arrive,Delta_V_tot)

% 4
launch = linspace(0,1307*d_s); % launch window (s),2017/1/1-2020/7/31
arrive = linspace(881*d_s,1856*d_s); % arrive window (s), 2019/6/1-2022/1/31
[launch, arrive] = meshgrid(launch,arrive);
Delta_V_tot = zeros(n,n);
for i=1:n
    for j = 1:n
        [R1, V1] = State_R0_V0_t(rE,vE,d_s.*launch(i,j));
        [R2, V2] = State_R0_V0_t(r1I,v1I,d_s.*arrive(i,j));
        [V1m, V2m] = LambertU_V1_V2_from_R1_R2_t(R1,R2,(d_s.*arrive(i,j)-d_s.*launch(i,j)));
        Delta_V_tot(i,j) = norm(V1m-V1) + norm(V2m-V2);
        %capping delta v
        if Delta_V_tot(i,j)>60
            Delta_V_tot(i,j)=NaN;
        end
    end
end
Delta_V_tot;
contour(launch,arrive,Delta_V_tot)



