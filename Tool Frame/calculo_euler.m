function [tita4,tita5,tita6] = calculo_euler(DH,tita1,tita2,tita3,TH_B)
% ======================== Calculo Euler ==================================
%
% Funcion que resuelve los posibles valores de las articulaciones 4, 5 y 6
% del robot una vez calculado los primeros tres.
%
% ---------
% Entradas
% ---------
%
%     ==> DH : matriz con los parametros de las primeras 3
%              columnas de la tabla de Denavit-Hartemberg ( alpha{i-1} | a{i-1} | d{i} )
%
%     ==> tita1 : Posibles soluciones para la 1° articulacion
%
%     ==> tita2 : Posibles soluciones para la 2° articulacion
%
%     ==> tita3 : Posibles soluciones para la 3° articulacion
%
%     ==> TH_B : matriz de posicion y orientacion deseadas para el manipulador 
%
% ---------
% Salidas
% ---------
%
%    ==> tita4 : Posibles soluciones para la 4° articulacion 
%
%    ==> tita5 : Posibles soluciones para la 5° articulacion 
%
%    ==> tita6 : Posibles soluciones para la 6° articulacion 
%
% -----------
% Referencias
% -----------
%
%       [Craig] John J. Craig, "Introduccion a la robotica", 
%               Pearson, 3nd Edicion, 2006.
%
% =========================================================================


%% Calculo de tita4, tita5 y tita6

p1.d=DH(1,3)     ;  p1.a=DH(1,2)   ; p1.alfa=DH(1,1)   ;p1.q=tita1;
p2.d=DH(2,3)     ;  p2.a=DH(2,2)   ; p2.alfa=DH(2,1)   ;p2.q=tita2;
p3.d=DH(3,3)     ;  p3.a=DH(3,2)   ; p3.alfa=DH(3,1)   ;p3.q=tita3;


%La función transf_DHcraig calcula T {i} ==> {i-1} para cada fila de la tabla
T1r0=transf_DHcraig(p1);
T2r1=transf_DHcraig(p2);
T3r2=transf_DHcraig(p3);

%Utilizó el cambio de base para obtener la Matriz T {3} ==> {0}
T3r0=T1r0*T2r1*T3r2;
R=T3r0\TH_B;

[tita4,tita5,tita6] =  solucion_angulos(R(1:3,1:3));

end

