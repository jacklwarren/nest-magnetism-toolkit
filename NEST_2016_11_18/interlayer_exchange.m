function Cexg=interlayer_exchange(d,exg_terms,exg_exp_0,exg_exp_1,exg_exp_2,exg_cos);

% function calculates Cexg in Hexg=M*Cexg*M
% For surface to surface separation between layers = d then:
% Cexg=[C0+C1*exp(-d/D1)+C2*exp(-2*d/D2)]*cos((2*pi*d/lambda)+phi)
% C0, C1 and C2 are all per unit area.

% Input parameters:
% d = (m) [d] in separation dependent terms
% exg_terms=[logical,logical,logical,logical]; % include the following terms:
% exg_exp_0= (/m^2) [C0] in Hexg=C0*M
% exg_exp_1= [/m^2,m]: [C1,D1] in Hexg=C1*exp(-d/D1)*M
% exg_exp_2= [/m^2,m]: [C2,D2] in Hexg=C2*exp(-2*d/D2)*M
% exg_cos= [m,radians]: [lambda,phi] in Cexg*CC*cos((2*pi*d/lambda)+phi)*M

% Note that to use the cosine term at least one of C0, C1 or C2 must be used.
% The direction of the cosine can be inverted by the sign of lambda and 
% the cosine term can be made into a sin by use of the phase angle phi.

Cexg=0.0;

if exg_terms(1)
    Cexg=Cexg+exg_exp_0;
end
if exg_terms(2)
    Cexg=Cexg+exg_exp_1(1)*exp(0.0-d/exg_exp_1(2));
end
if exg_terms(3)
    Cexg=Cexg+exg_exp_2(1)*exp(0.0-2.0*d/exg_exp_2(2));
end
% Only multiply by the cosine term if there is a multiplying constant
if exg_terms(4) && (exg_terms(1) || exg_terms(2) ||exg_terms(3))
        Cexg=Cexg*cos((2*pi*dz/exg_cos(1))+exg_cos(2));
end

end


