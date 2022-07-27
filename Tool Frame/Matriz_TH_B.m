function [ TH_B ] = Matriz_TH_B( p, angulos_ZYX )
%Funcion para armar la matriz de la herramienta respecto a la base
%
% Entradas:
%       p     : Vector de traslacion de la matriz de roto-traslacion
%
%  angulos_ZYX: Angulos Z, Y, X para la matriz de orientacion segun convenio de Euler (radianes).     
%
% Salida:
%       TH_B: Matriz de rototraslacion de la herramienta respecto a la base.
%


Rot = rotz(angulos_ZYX(1))*roty(angulos_ZYX(2))*rotx(angulos_ZYX(3));

TH_B = zeros(4,4);
TH_B(1:3,4) = p;
TH_B(1:3,1:3) = Rot;
TH_B(4,4) = 1;


end

