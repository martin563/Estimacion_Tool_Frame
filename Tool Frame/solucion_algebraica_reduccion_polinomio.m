function [ tita ] = solucion_algebraica_reduccion_polinomio( a, b, c )
% ====================== Metodo de Pieper =================================
%
% Funcion que resuelve la ecuacion a*cos(tita) + b*sin(tita) = c
% como se explica en [Craig, pag.114, ejemplo 4.3].
%
% ---------
% Entradas
% ---------
%
%     ==> a : constante a de la ecuacion a resolver.
%
%     ==> b : constante b de la ecuacion a resolver.
%
%     ==> c : constante c de la ecuacion a resolver.
%            
% ---------
% Salidas
% ---------
%
%    ==> tita : Posibles soluciones de la ecuacion a resolver. 
%
% -----------
% Referencias
% -----------
%
%       [Craig] John J. Craig, "Introduccion a la robotica", 
%               Pearson, 3nd Edicion, 2006.
%
% =========================================================================

if abs(a+c) > 1e-8
    
    %Solucion para la ecuacion de 2° grado
    raiz = b^(2)+a^(2)-c^(2);
    if raiz >= 0
        aux1=2*atan2( (b+sqrt(raiz)), (a+c) );
        aux2=2*atan2( (b-sqrt(raiz)), (a+c) );
        
        aux1 = atan2(sin(aux1),cos(aux1));
        aux2 = atan2(sin(aux2),cos(aux2));
        tita = [aux1;aux2];
    else
        aux1 = 0;
        aux2 = 0;
        tita = [aux1;aux2];
        disp('No hay solucion real para la configuracion dada - (SOLUCIÓN ALGEBRAICA POR REDUCCIÓN A POLINOMIO)');
    end
    
else
    
    %Solucion para la ecuacion de 1° grado 
%     aux1 = 2*atan2( (c - a), (2*b) );
%     aux2 = aux1;
    tita = 2*atan2( (c - a), (2*b) );
    
end

% aux1 = atan2(sin(aux1),cos(aux1));
% aux2 = atan2(sin(aux2),cos(aux2));
% tita = [aux1;aux2];


end

