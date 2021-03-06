%MODEL NODAL ANALYSIS
% JASON PRASAD 100196970; 

% define all intial paramters Resistance to Conductance 
r_0= 1*10^3
r1=1; 
r2=2; 
r3=10;
r4=0.1; 

g_0=1/r_0; 
g1=1/r1; 
g2=1/r2; 
g3=1/r3; 
g4=1/r4; 

%initialize remaining paramteres for the circuit
cap=0.25
alpha=100;
l=0.2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A) create C & G matrix stuck with my orginal equations since there were
% matrices errors 

C = [0, 0, 0, 0, 0, 0;
    cap, -cap, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, l];

G = [1, 0, 0, 0, 0, 0;
    g1, -(g1 + g2), 0, 0, 0, -1;
   0, 0, -g3, 0, 0, 1;
   0, 0, -alpha, 1, 0, 0;
   0, 0, 0, g4, -(g4  + g_0), 0;
   0, -1, 1, 0, 0, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% B) DC SWEEP GRAPH DOESNT SEEM TO BE RIGHT????? 
F= zeros(6, 1); 
V_in=zeros(100,1);
V_out= zeros(100,1); 
V_3=zeros(100,1);
count=0; 
for V_in = -10:0.1:10
    count = count +1; 
    F(1)= V_in; 
    V_val= G\F; 
    
    V_in(count)=V_in; 
    V_3(count)=V_val(5);
    V_out(count)=V_val(6); 
    
end
figure (1)
plot(V_in,V_out);
plot(V_in,V_3);
hold on; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% C) AC SWEEP 

V_in=10;
w_min=0; 
w_max=100; 
m= 1000; 
omega_o=zeros(m,1);
V_out=zeros(m,1);
Unity_g=zeros(m,1);
count_2=0;
for w=linspace(w_min,w_max,m);
    count_2 =count_2 +1; 
    Temp_1= G +1*j*w*C; 
    V= Temp_1 \ F; 
    omega_o(count_2)=w; 
    V_out(count_2)=norm(V(5));
    Unity_g(count_2)=20*log10(norm(V(5))/V_in); 
    
end 
figure (2) 
plot(omega_o, V_out); 
hold on
plot(omega_o, Unity_g);
   
%random perturbations
n = 1000;
w = pi();
Crand = zeros(n, 1);
Grand = zeros(n, 1);

for index = 1:n
    ci = c + 0.5*randn();
    
    C = [0, 0, 0, 0, 0, 0;
    ci, -ci, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, L];

    Temp_2 = G + 1j*w*C;
    V = Temp_2\F;
    
    Crand(index) = ci;
    Grand(index) = 20*log10(norm(V(5))/V_in);
end

figure(4)
histogram(Crand)
figure(5)
histogram(Grand)

