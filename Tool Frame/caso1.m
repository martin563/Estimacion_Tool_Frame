function [tita1,tita2,tita3] = caso1(alpha,a,d,TH_B)
% =========================== Caso 1 ======================================
%
% Resuelve tita1, tita2 y tita3 para el caso a1=0 y sin(alpha1)~=0
%
% ---------
% Entradas
% ---------
%
%     ==> alpha : angulos alpha{i-1} de la tabla DH (radianes)
%
%     ==> a : a{i-1} de la tabla DH (metros)
%
%     ==> d : d{i} de la tabla DH (metros)
%
%     ==> TH_B : matriz de posicion y orientacion deseadas para el manipulador  
%            
% ---------
% Salidas
% ---------
%
%    ==> tita1 : Posibles soluciones para la 1° articulacion 
%
%    ==> tita2 : Posibles soluciones para la 2° articulacion            
%
%    ==> tita3 : Posibles soluciones para la 3° articulacion
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


%% Resolucion tita3

    cte=a(1)^(2)+a(2)^(2)+a(3)^(2)+d(2)^(2)+d(3)^(2)+d(4)^(2)+2*d(3)*d(4)*cos(alpha(3))+2*d(2)*(d(4)*cos(alpha(2))*cos(alpha(3))+d(3)*cos(alpha(2)) ); 
    r=TH_B(1,4)^(2)+TH_B(2,4)^(2)+TH_B(3,4)^(2) ;
    a3=2*a(2)*a(3)-2*d(2)*d(4)*sin(alpha(3))*sin(alpha(2));
    b3=2*a(2)*d(4)*sin(alpha(3))+2*d(2)*a(3)*sin(alpha(2));
    c3=r-cte;
        
    tita3 = solucion_algebraica_reduccion_polinomio( a3, b3, c3 );
    N_sol3 = length(tita3);
    
%% Resolucion tita2
%Parametros a definir:    
for i=1:1:N_sol3;    
    f1=a(3)*cos(tita3(i))+d(4)*sin(alpha(3))*sin(tita3(i))+a(2);
    f2=a(3)*cos(alpha(2))*sin(tita3(i))-d(4)*sin(alpha(3))*cos(alpha(2))*cos(tita3(i))-d(4)*sin(alpha(2))*cos(alpha(3))-d(3)*sin(alpha(2));
    f3=a(3)*sin(alpha(2))*sin(tita3(i))-d(4)*sin(alpha(3))*sin(alpha(2))*cos(tita3(i))+d(4)*cos(alpha(2))*cos(alpha(3))+d(3)*cos(alpha(2));

    k1=f1;
    k2=-f2;
    k3=f1^(2)+f2^(2)+f3^(2)+a(1)^(2)+d(2)^(2)+2*d(2)*f3;
    k4=f3*cos(alpha(1))+d(2)*cos(alpha(1));
    z=TH_B(3,4);

    a3=-k2;
    b3=k1;
    c3=(z-k4)/sin(alpha(1));
       
    tita2(:,i) = solucion_algebraica_reduccion_polinomio( a3, b3, c3 );
    N_sol2 = length(tita2(:,i));

%% Calculo tita1
for j=1:1:N_sol2
    g1=cos(tita2(j,i))*f1-sin(tita2(j,i))*f2+a(1);
    g2=sin(tita2(j,i))*cos(alpha(1))*f1+cos(tita2(j,i))*cos(alpha(1))*f2-sin(alpha(1))*f3-d(2)*sin(alpha(1));
    g3=sin(tita2(j,i))*sin(alpha(1))*f1+cos(tita2(j,i))*sin(alpha(1))*f2+cos(alpha(1))*f3+d(2)*cos(alpha(1));

    matrix=[g1,g2;-g2,g1];
    v=(1/(g1^(2)+g2^(2)))*matrix*TH_B(1:2,4);

    tita1(j,i)=atan2(v(2),v(1));

end

end

