function [tita4,tita5,tita6 ] = solucion_angulos(R)
% ======================== Solucion angulos ===============================
%
% Funcion que devuelve el conjunto de angulos que resulve la rotacion
% inversa (R6r3) que tenga la siguiente forma:
%
% R6r3 = 
%
%   [   c4*c5*c6 - s4*s6, - c6*s4 - c4*c5*s6, -c4*s5]
%   [              c6*s5,             -s5*s6,     c5]
%   [ - c4*s6 - c5*c6*s4,    c5*s4*s6 -c4*c6,  s4*s5]
%
%
% Siendo c4 = cos(q4), s4 = sin(s4), c5 = cos(q5), s5 = sin(s5), 
% c6 = cos(q6) y s6 = sin(s6).
% 
% ---------
% Entradas
% ---------
%
%     ==> R : Matriz de rotacion R6r3.
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

tita5(1)=atan2(sqrt(1-R(2,3)^(2)),R(2,3));
tita5(2)=atan2(-sqrt(1-R(2,3)^(2)),R(2,3));

tita4(1)=atan2(R(3,3),-R(1,3));
tita4(2)=atan2(-R(3,3),R(1,3));

tita6(1)=atan2(-R(2,2),R(2,1));
tita6(2)=atan2(R(2,2),-R(2,1));

end

