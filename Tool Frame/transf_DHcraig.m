function [T] = transf_DHcraig(param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%clear
%syms q1 q2 q3 ;

alfa=param.alfa;
a=param.a;
d=param.d;
q=param.q;

cq=cos(q);
sq=sin(q);
calfa=cos(alfa);
salfa=sin(alfa);

if abs(calfa)<1e-10,calfa=0;end
if abs(salfa)<1e-10,salfa=0;end


T=[cq        -sq         0          a         ;...
   sq*calfa  cq*calfa  -salfa       -salfa*d  ;...
   sq*salfa  cq*salfa   calfa       calfa*d   ;...
   0         0                0      1];


end

