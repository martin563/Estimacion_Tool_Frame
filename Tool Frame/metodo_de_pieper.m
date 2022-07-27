function [conj, N_sol] = metodo_de_pieper(DH, TH_B)
% ====================== Metodo de Pieper =================================
%
% Funcion que resuelve la cinematica inversa de un manipulador serial de 
% 6 articulaciones con las caracteristicas de que las 3 ultimas 
% articulaciones coincidan con su origen.
%
% La tabla correspondiente para que funcione el algoritmo debe ser :
%
%   i | alpha{i-1} | a{i-1} | d{i} | tita{i}
%   1 |     0      |    0   |   0  |   q1
%   2 |    alpha1  |   a1   |   d2 |   q2
%   3 |    alpha2  |   a2   |   d3 |   q3
%   4 |    alpha3  |   a3   |   d4 |   q4
%   5 |    alpha4  |   0    |   0  |   q5
%   6 |    alpha5  |   0    |   0  |   q6
%
% NOTAS: 
%   - Los alpha debe estar en radianes
%
%   - Solo esta el caso en el cual a1 = 0 y sin(alpha1) ~= 0.
%
% ---------
% Entradas
% ---------
%
%     ==> DH : matriz con los parametros de las primeras 3
%              columnas de la tabla de Denavit-Hartemberg ( alpha{i-1} | a{i-1} | d{i} )
%
%     ==> TH_B : matriz de posicion y orientacion deseadas para el manipulador                  
%            
% ---------
% Salidas
% ---------
%
%    ==> conj : Conjunto de posibles soluciones (en radianes)   
%
% -----------
% Referencias
% -----------
%
%       [Craig] John J. Craig, "Introduccion a la robotica", 
%               Pearson, 3nd Edicion, 2006.
%
%       [Piper] Donald Lee Pieper, "The Kinematics of manipulators under 
%               computer control", Stanford University, 1968. 
%
% =========================================================================

%Traspaso de datos:
alpha=DH(2:6,1);
a=DH(2:4,2);
d=DH(1:4,3);


%% Solucion para los titas
if a(1)==0 && sin(alpha(1))~=0
    
    [tita1,tita2,tita3]=caso1(alpha,a,d,TH_B);
    
    N3 = length(tita3);
    N2 = length(tita2);
    N1 = length(tita1);
    N_sol = 0;
    
    if N3 == 1 && N2 == 2 && N1 == 2
        N_sol = 4;
        conj = zeros(N_sol,6);
        [tita4_1,tita5_1,tita6_1]=calculo_euler(DH,tita1(1,1),tita2(1,1),tita3(1),TH_B);
        [tita4_2,tita5_2,tita6_2]=calculo_euler(DH,tita1(2,1),tita2(2,1),tita3(1),TH_B);
        
        tita4(1,:)=tita4_1(:); tita4(2,:)=tita4_2(:);
        tita5(1,:)=tita5_1(:); tita5(2,:)=tita5_2(:);
        tita6(1,:)=tita6_1(:); tita6(2,:)=tita6_2(:);
        conj(1,:) = [tita1(1,1),tita2(1,1),tita3(1),tita4(1,1),tita5(1,1),tita6(1,1)];
        conj(2,:) = [tita1(1,1),tita2(1,1),tita3(1),tita4(1,2),tita5(1,2),tita6(1,2)];
        conj(3,:) = [tita1(2,1),tita2(2,1),tita3(1),tita4(2,1),tita5(2,1),tita6(2,1)];
        conj(4,:) = [tita1(2,1),tita2(2,1),tita3(1),tita4(2,2),tita5(2,2),tita6(2,2)];
        
    elseif N3 == 2 && N2 == 1 && N3 == 2
        N_sol = 4;
        conj = zeros(N_sol,6);
        [tita4_1,tita5_1,tita6_1]=calculo_euler(DH,tita1(1,1),tita2(1),tita3(1),TH_B);
        [tita4_2,tita5_2,tita6_2]=calculo_euler(DH,tita1(2,1),tita2(1),tita3(2),TH_B);
        
        tita4(1,:)=tita4_1(:); tita4(2,:)=tita4_2(:);
        tita5(1,:)=tita5_1(:); tita5(2,:)=tita5_2(:);
        tita6(1,:)=tita6_1(:); tita6(2,:)=tita6_2(:);
        conj(1,:) = [tita1(1,1),tita2(1),tita3(1),tita4(1,1),tita5(1,1),tita6(1,1)];
        conj(2,:) = [tita1(1,1),tita2(1),tita3(1),tita4(1,2),tita5(1,2),tita6(1,2)];
        conj(3,:) = [tita1(2,1),tita2(1),tita3(2),tita4(2,1),tita5(2,1),tita6(2,1)];
        conj(4,:) = [tita1(2,1),tita2(1),tita3(2),tita4(2,2),tita5(2,2),tita6(2,2)];
        
    elseif N3 == 2 && N2 == 2 && N1 == 1
        N_sol = 4;
        conj = zeros(N_sol,6);
        [tita4_1,tita5_1,tita6_1]=calculo_euler(DH,tita1(1),tita2(1,1),tita3(1),TH_B);
        [tita4_2,tita5_2,tita6_2]=calculo_euler(DH,tita1(1),tita2(2,1),tita3(2),TH_B);
        
        tita4(1,:)=tita4_1(:); tita4(2,:)=tita4_2(:);
        tita5(1,:)=tita5_1(:); tita5(2,:)=tita5_2(:);
        tita6(1,:)=tita6_1(:); tita6(2,:)=tita6_2(:);
        conj(1,:) = [tita1(1),tita2(1,1),tita3(1),tita4(1,1),tita5(1,1),tita6(1,1)];
        conj(2,:) = [tita1(1),tita2(1,1),tita3(1),tita4(1,2),tita5(1,2),tita6(1,2)];
        conj(3,:) = [tita1(1),tita2(2,1),tita3(2),tita4(2,1),tita5(2,1),tita6(2,1)];
        conj(4,:) = [tita1(1),tita2(2,1),tita3(2),tita4(2,2),tita5(2,2),tita6(2,2)];
               
    elseif N3 == 2 && N2 == 2 && N1 == 2
        N_sol = 8;
        conj = zeros(N_sol,6);    
        [tita4_1,tita5_1,tita6_1]=calculo_euler(DH,tita1(1,1),tita2(1,1),tita3(1),TH_B);
        [tita4_2,tita5_2,tita6_2]=calculo_euler(DH,tita1(2,1),tita2(2,1),tita3(1),TH_B);
        [tita4_3,tita5_3,tita6_3]=calculo_euler(DH,tita1(1,2),tita2(1,2),tita3(2),TH_B);
        [tita4_4,tita5_4,tita6_4]=calculo_euler(DH,tita1(2,2),tita2(2,2),tita3(2),TH_B);
        
        tita4(1,:)=tita4_1(:); tita4(2,:)=tita4_2(:); tita4(3,:)=tita4_3(:); tita4(4,:)=tita4_4(:);
        tita5(1,:)=tita5_1(:); tita5(2,:)=tita5_2(:); tita5(3,:)=tita5_3(:); tita5(4,:)=tita5_4(:);
        tita6(1,:)=tita6_1(:); tita6(2,:)=tita6_2(:); tita6(3,:)=tita6_3(:); tita6(4,:)=tita6_4(:);
        conj(1,:) = [tita1(1,1),tita2(1,1),tita3(1),tita4(1,1),tita5(1,1),tita6(1,1)];
        conj(2,:) = [tita1(1,1),tita2(1,1),tita3(1),tita4(1,2),tita5(1,2),tita6(1,2)];
        conj(3,:) = [tita1(2,1),tita2(2,1),tita3(1),tita4(2,1),tita5(2,1),tita6(2,1)];
        conj(4,:) = [tita1(2,1),tita2(2,1),tita3(1),tita4(2,2),tita5(2,2),tita6(2,2)];
        conj(5,:) = [tita1(1,2),tita2(1,2),tita3(2),tita4(3,1),tita5(3,1),tita6(3,1)];
        conj(6,:) = [tita1(1,2),tita2(1,2),tita3(2),tita4(3,2),tita5(3,2),tita6(3,2)];
        conj(7,:) = [tita1(2,2),tita2(2,2),tita3(2),tita4(4,1),tita5(4,1),tita6(4,1)];
        conj(8,:) = [tita1(2,2),tita2(2,2),tita3(2),tita4(4,2),tita5(4,2),tita6(4,2)];
    end
    
    disp(sprintf('El numero de posibles soluciones es = %d ', N_sol));    
         
end


end

