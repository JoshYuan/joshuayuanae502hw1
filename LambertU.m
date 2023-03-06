

% Joshua Yuan
% AE 502 HW1 Lambert solver
close all; close all; clc;

% Lambert's problem
% Determine elliptic trajectory
% that joins two points
% in a specified flight time

% specify time of flight, transfer angle theta
prompt = {'Specify time of flight (days):',...
    'Specify transfer angle \theta (in degrees):',...
    'Specify r2/r1:'};
dlgtitle = 'Inputs';
definput = {'365 115 500 115','75 75 285 285','1.524 1.524 1.524 1.524'}; 
dims = [1 40];
opts.Interpreter = 'tex';
inputs = inputdlg(prompt,dlgtitle,dims,definput,opts);

t_f_in = str2num(inputs{1});
theta_in = str2num(inputs{2});
r2in = str2num(inputs{3});
n = length(t_f_in);

r1=1;

% t_f_days = 115;
% theta_degrees = 75;

% normalize units
r = 1.495978*10^8;  % r input in km
mu = 1.327*10^11;   % mu of problem in km3/s2
ctu_s = sqrt(r^3/mu); % canonical time unit in seconds
ctu = ctu_s/(60*60*24); % canonical time unit wrt input time unit

for i = 1:n
    t_f_days = t_f_in(i);
    theta_degrees = theta_in(i);
    r2 = r2in(i);
    
    fprintf('For time of flight: tf = %3i days and transfer angle: theta = %3i degrees\n',...
        t_f_days,theta_degrees);
    
    t_f = t_f_days/ctu;    %tf in time units
    theta = deg2rad(theta_degrees);
    
%     r1 = 1;
    c = sqrt(r1^2+r2^2-2*r1*r2*cos(theta));
    s = (r1 + r2 + c)/2;
    a_m = s/2;
    
    alpha_m = pi; %2*asin(sqrt(s/(2*a_m)));
    beta_m = 2*asin(sqrt(1-c/s));%2*asin(sqrt((s-c)/(2*a_m)));
    t_m = sqrt(s^3/8)*(pi-beta_m+sin(beta_m));
    t_p = sqrt(2)/3*(s^(3/2)-sign(sin(theta))*(s-c)^(3/2));
    
    %syms a_s
    alpha0 = @(a) 2*asin(sqrt(s/(2*a)));
    beta0 = @(a) 2*asin(sqrt((s-c)/(2*a)));
    
    % choose beta
    if 0<=theta && theta<pi
        beta_fcn = beta0;
    elseif pi<=theta && theta<2*pi
        beta_fcn = @(a) -1*beta0(a);
    end
    
    % choose alpha
    if t_f <= t_m
        alpha_fcn = alpha0;
    elseif t_f > t_m
        alpha_fcn = @(a) 2*pi-alpha0(a);
    end
    
    f = @(a) t_f - ((a^(3/2))*(alpha_fcn(a)-beta_fcn(a)-sin(alpha_fcn(a))+sin(beta_fcn(a))));
    upper = 4*s/abs(sin(theta));
    % a only gets large if theta is very small, which is accounted for here
    a = fzero(f,[a_m upper]);
    alpha = alpha_fcn(a);
    beta = beta_fcn(a);
    
    p = 4*a*(s-r1)*(s-r2)/c^2*sin((alpha+beta)/2)^2;
    e = sqrt(1-p/a);
    
    r1v = [r1 0];
    r2v = r2*[cos(theta) sin(theta)];
    u1 = r1v/r1;
    u2 = r2v/r2;
    uc = (r2v-r1v)/c;
    A = sqrt(1/4/a)*cot(alpha/2);
    B = sqrt(1/4/a)*cot(beta/2);
    v1 = (B+A)*uc+(B-A)*u1;
    v2 = (B+A)*uc-(B-A)*u2;
    
    fprintf('In canonical units (mu=1, r=1):\ta = %4.3f\te = %1.3f\tv1 = %1.3fi + %1.3fj\n\n',a,e,v1);
end