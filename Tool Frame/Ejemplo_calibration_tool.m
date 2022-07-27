%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fecha: 30/01/2022
% Autor: Martin Oviedo
% Gmail: eliasoviedo1718@gmail.com
%
% Resumen :
%    Ejemplo de algoritmo para calibracion del tool frame con 4 puntos con
%    el robot Puma 560.
%
%       Seccion A) ==> Cargo datos del robot puma
%
%       Seccion B) ==> Simulo las muestras a tomar con diferentes orientaciones
%
%       Seccion C) ==> Aplico el metodo para estimar la traslacion de T6_H
%
% -----------
% Referencias
% ----------- 
%
%       [Katayama] Tohru Katayama, "Subspace Methods for System Identification", 
%                  Springer, 1nd Edition. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear all, close all

%% Constantes

mm = 10^(-3);

%% A) Cargo datos del robot puma

syms q1 q2 q3 q4 q5 q6 real

%Carga de datos del PUMA 560 - Toolbox Corke
mdl_puma560; 

%Tabla de Denavit-Hartenberg con Toolbox Simbolico Robot Puma 560     
p1.d=0       ; p1.a=0       ; p1.alfa=0       ;p1.q=q1;
p2.d=0       ; p2.a=0       ; p2.alfa=-pi/2   ;p2.q=q2;
p3.d=0.15    ; p3.a= 0.4318 ; p3.alfa=0       ;p3.q=q3;
p4.d=0.4318  ; p4.a=0.0203  ; p4.alfa=-pi/2   ;p4.q=q4;
p5.d=0       ; p5.a=0       ; p5.alfa=pi/2    ;p5.q=q5;
p6.d=0       ; p6.a=0       ; p6.alfa=-pi/2   ;p6.q=q6;

% La función transf_DHcraig calcula T {i} ==> {i-1} para cada fila de la tabla
T1r0=transf_DHcraig(p1);
T2r1=transf_DHcraig(p2);
T3r2=transf_DHcraig(p3);
T4r3=transf_DHcraig(p4);
T5r4=transf_DHcraig(p5);
T6r5=transf_DHcraig(p6);

% Utilizó el cambio de base para obtener la Matriz T {H} ==> {B}
T6r0=simplify(T1r0*T2r1*T3r2*T4r3*T5r4*T6r5);


%% B) Simulo las muestras a tomar con diferentes orientaciones

%Parametros que se pasan a metodo_de_piper para la cinematica inversa
DH=[ 0    , 0      , 0 ;  
    -pi/2 , 0      , 0;
      0   , 0.4318 , 0.15;
    -pi/2 , 0.0203 , 0.4318 ;
    pi/2  , 0      , 0;
    -pi/2 , 0      , 0;
];


%Tool a determinar
Tool = 300*mm;
T6_H = Matriz_TH_B( [0;0;Tool], [0,0,0] );


%Construccion de la matriz TH_B con diferentes orientaciones
ang_a_rad = pi/180;
tita = 2*mm;
p = 0.5*ones(3,1,4) + rand(3,1,4)*tita;

%Angulos de muestra
angulos_ZYX_a = [90, 70, 70]*ang_a_rad;
angulos_ZYX_b = [30, 0, 40]*ang_a_rad;
angulos_ZYX_c = [40, 40, -30]*ang_a_rad;
angulos_ZYX_d = [-20, 20, 30]*ang_a_rad;

TH_B_a = Matriz_TH_B( p(:,1,1), angulos_ZYX_a );
TH_B_b = Matriz_TH_B( p(:,1,2), angulos_ZYX_b );
TH_B_c = Matriz_TH_B( p(:,1,3), angulos_ZYX_c );
TH_B_d = Matriz_TH_B( p(:,1,4), angulos_ZYX_d );

TB_6_a = TH_B_a/T6_H;
TB_6_b = TH_B_b/T6_H;
TB_6_c = TH_B_c/T6_H;
TB_6_d = TH_B_d/T6_H;

%Conjunto de posibles soluciones
[conj_a, N_sol_a] = metodo_de_pieper(DH,TB_6_a); 
[conj_b, N_sol_b] = metodo_de_pieper(DH,TB_6_b);
[conj_c, N_sol_c] = metodo_de_pieper(DH,TB_6_c);
[conj_d, N_sol_d] = metodo_de_pieper(DH,TB_6_d);

%% C) Aplico el metodo para estimar la traslacion de T6_H

sol = zeros(4,4,4);
sol(:,:,1) = TB_6_a;
sol(:,:,2) = TB_6_b;
sol(:,:,3) = TB_6_c;
sol(:,:,4) = TB_6_d;

R = zeros(3,3,4);
P = zeros(3,1,4);

for i = 1:4
    R(:,:,i) = sol(1:3,1:3,i);
    P(:,:,i) = sol(1:3,4,i);
end

C = 6;
F = zeros(3*C, 3);
B = zeros(3*C, 1);

F(1:3,:) = R(:,:,1) - R(:,:,2);
F(4:6,:) = R(:,:,1) - R(:,:,3);
F(7:9,:) = R(:,:,1) - R(:,:,4);
F(10:12,:) = R(:,:,2) - R(:,:,3);
F(13:15,:) = R(:,:,2) - R(:,:,4);
F(16:18,:) = R(:,:,3) - R(:,:,4);


B(1:3,:) = -P(:,:,1) + P(:,:,2);
B(4:6,:) = -P(:,:,1) + P(:,:,3);
B(7:9,:) = -P(:,:,1) + P(:,:,4);
B(10:12,:) = -P(:,:,2) + P(:,:,3);
B(13:15,:) = -P(:,:,2) + P(:,:,4);
B(16:18,:) = -P(:,:,3) + P(:,:,4);


%Descomposicion en valores singulares
rango = rank(F);
[U,S,V] = svd(F);
U_aux = U(:,1:rango);
S_aux = S(1:rango,1:rango);
V_aux = V(:,1:rango);

Tool_est = V_aux/S_aux;
Tool_est = Tool_est*transpose(U_aux)*B;
Error = abs(Tool_est(3) - Tool);

disp('------------ Resultados ---------------')
disp(sprintf('Tool_est = %.4d [mm]', Tool_est(3)/mm))
disp(sprintf('Error = %.4d [mm]', Error/mm))
disp('---------------------------------------')

