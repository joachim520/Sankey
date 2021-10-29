clear
close all
clc
###################################################################
% Simple Sankey Diagram Generator, See sankey.pdf for info
###################################################################
% User inputs:
W_0 = 300; % Input power [W]
L = (10/3) * W_0; % Horizontal span of diagram
efficiency = [0.98,0.92,0.7,0.85]; % Efficiencies of the different processes
###################################################################

% Calculated initial variables
n = length(efficiency); % Amount of (lossful) processes in the system
L_0 = 0.18 * L; % Fixed length before first processes
H_0 = W_0 / 2; % Initial arrow height declaration;

% Calculating first system variables (index = 1)
W_i = efficiency(1) * W_0;
W = [W_i];

Ls = [(1 * (L - L_0) )/n];

t_i = W_0 - W(1);
t = [t_i];

H = [H_0 + t_i];

DataValues = [W;Ls;t;H];

% Calculating first (unique) points (index = 1)
a = [ 0 ; 0 ];
b = [ 0 ; -W_0 ];
c = [ L_0 ; -W_0 ];
d = [ L_0 ; -W_0 - H_0 ];

e = [L_0 + t_i ; -W_0 - H_0 ];
f = [L_0 + ( t_i / 2 ) ; -W_0 - H_0 - ( t_i / 2 )];
g = [L_0 + t_i ; -W_i ];

q = [ L ; 0];

% Declaring sequence points (index > 1)
h = [;];
j = [;];
k = [;];
l = [;];
m = [;];
p = [;];

% Create the matrix to be plotted later
DiagramPoints = [a,b,c,d,f,e,g];


if n > 1
 for i = 2:n

   % Calculating system variables 
   W_i = efficiency(i) * W(i-1);
   L_i = (i * (L - L_0) )/n;
   t_i = W(i-1) - W_i;
   H_i = H(i-1) + t_i;
   % Storing
   W = [W,W_i];
   Ls = [Ls,L_i];
   t = [t,t_i];
   H = [H,H_i];

   DataValue = [W_i;L_i;t_i;H_i];
   DataValues = [DataValues,DataValue]; % Handy matrix to return
      
  end
  % Calcutating the point coordinates
  
  for i = 2:n
   i = i-1;
   % Calculating sequence points (index > 1)
   h_i = [ L_0 + Ls(i) ; -W(i) ];
   j_i = [ L_0 + Ls(i) ; -W(i) - H(i) ];
   k_i = [ L_0 + Ls(i) + t(i+1) ; -W(i) - H(i) ];
   l_i = [ L_0 + Ls(i) + ( t(i+1) / 2 ) ; -W(i) - H(i) - ( t(i+1) / 2 ) ];
   m_i = [ L_0 + Ls(i) + t(i+1) ; -W(i+1) ];
   p_i = [ L_0 + Ls(i+1) ; -W(i+1) ];
   % Storing
   h = [h,h_i];
   j = [j,j_i];
   k = [k,k_i];
   l = [l,l_i];
   m = [m,m_i];
   p = [p,p_i];
   DiagramPoints = [DiagramPoints,h_i,j_i,l_i,k_i,m_i,p_i]; % Adding sequence points to the point matrix
   
   if i == n-1
   % Calculating final point (index = n)
   r = [ L + ( W(i+1) / 2 ) ; ( -W(i+1) / 2 ) ];
   DiagramPoints = [DiagramPoints,r,q,a]; % Adding the final point to the point matrix
   end
  end
end

% Calculating total system efficiency
SystemEfficiency = 1;
for i = 1:n
  SystemEfficiency = SystemEfficiency * efficiency(i);
 end
SystemEfficiency % Return system efficiency

% Plotting the diagram 
figure
hold on
plot(DiagramPoints(1,1:length(DiagramPoints(1,:))),DiagramPoints(2,1:length(DiagramPoints(1,:))),'b', "linewidth", 3)
axis equal
hold off
