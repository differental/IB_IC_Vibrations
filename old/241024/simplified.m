% A matlab program to plot the theoretical response of the building in the
% 1A vibration lab.
% Based on code written by Penny Cox, now maintained by Aidan Reillym, tidied up by Jim Woodhouse.
% Partially rewritten by Brian Chen, October 2024, with calculations done
% by Mihnea Butiurca
%
close all
clear all

%  USER MUST SET TWO VALUES TO DETERMINE PLOT OPTIONS


m = 1.83; % mass of one floor
L = 0.2; % length
N = 3; % number of degrees of freedom
b = 0.08; % width
E = 210E9; % Young's Modulus
d = 0.001; % thickness
I = b*d*d*d/12; % second moment of area
k = (24*E*I)/(L*L*L); % static stiffness for each floor




% Natural frequencies of the original structure without absorbers:
% f (Hz): 3.3933    9.5078   13.7391
% We want to tune the natural frequencies of each absorber
% sqrt(k/m)/2pi
% to one of the natural frequencies above (absorber mass negligible)
% k = (2 * pi * f) ** 2 * m

natural_frequencies = [3.3933, 9.5078, 13.7391];

n = 30; % number of absorbers
absorber_total_mass = 0.3;
absorber_m = absorber_total_mass / n * ones(n, 1); % array of absorber mass, length n
%absorber_location = [1, 1, 1, 1, 2, 3, 1, 2, 3, 1];
absorber_location = 3 * ones(n, 1);


absorber_k = zeros(n, 1); % initialise array of absorber stiffness, length n

for i = 1:n
    selected_frequency = natural_frequencies(mod(i, 3) + 1);
    absorber_k(i) = (2 * pi * selected_frequency)^2 * m; % tune absorber frequency to match one of the natural frequencies
    r = 0.95 + (1.05 - 0.95) * rand; % random number between 95% and 105%
    absorber_k(i) = absorber_k(i) * r;
end



M = m*eye(N+n);  % Create an identity matrix scaled by m

% Assign values from absorber_m to diagonal positions starting from M(4,4)
for i = 1:n
    M(i+N, i+N) = absorber_m(i);
end

%K = k*[2 -1 0;-1 2 -1;0 -1 1]; % create the stiffness matrix




K = zeros(N+n, N+n);  % Pre-allocate the stiffness matrix with zeros

K(1,1) = 2*k;
K(1,2) = -k;
K(2,1) = -k;
K(2,2) = 2*k;
K(3,3) = k;
K(2,3) = -k;
K(3,2) = -k;

for i = 1:n
    loc = absorber_location(i);
    K(loc, loc) = K(loc, loc) + absorber_k(i);
    K(3+i, 3+i) = absorber_k(i);
    K(loc, 3+i) = -absorber_k(i);
    K(3+i, loc) = -absorber_k(i);
end

disp(K);





% To include vibration absorbers, you will need to modify
%   the mass and stiffness matrices (above)

[V,D] = eig(K,M);
syms w;

for imode=1:N
  freqs(imode) = sqrt(D(imode,imode));
end

%  Print natural frequencies and mode vectors in command window
hertz = freqs/(2*pi);
modeshapes = V;


F = zeros(N+n, 1);
F(1) = 1;




B = K - ((w*w)*M); 
% harmonic solution for unit force at floor 1
disp = (inv(B))*F;

%start of ezplot section
  hold on
  
  ifloor=1;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line'),'Color','k')
  ifloor=2;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line','Color','b'),'Color','g')
  ifloor=3;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line','Color','b'),'Color','r')
  
  set(findobj('Type','line','Color','k'),'Color','b')
  
  set(findobj('Type','line'),'LineStyle','-')
  hold off

% Calculate frequency response functions
all_disp = [];
for w = 1:130
  B = K - ((w*w)*M); 
  % harmonic solution for unit force at floor 1
  disp = (inv(B))*F;
  all_disp = [all_disp disp];
end
