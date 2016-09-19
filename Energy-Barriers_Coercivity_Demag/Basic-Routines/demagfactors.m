%% Demagnetising factors calculation 
% This function calculates demagnetising factors for truncated elliptic
% cones, elliptic cylinders and prisms. All units are assumed to be in SI
% units, and dimensions are specified in nm. The magnetisation in all the
% geometries here are assumed to be uniform. The method of calculating the
% demagnetising factors are given in Kalezhi et al. [1].
% The integral equations used to calculate demagnetising factors were
% determined using Fourier transform approach proposed by Beleggia et al.
% [2]. The thorough derivation of the magnetometric demagnetising factors
% is given by Kalezhi et al. [1]. For a truncated elliptic cone the
% demagnetisation factors are given by:
% 
% $$N_{xx}= \frac{ab}{2\pi V} \int\limits_{x=0}^{\pi}\cos x dx  \int\limits_{\phi_k = 0}^{2\pi} \frac{\cos^2\phi_k} {\sqrt{\cos^2\phi_k + \beta^{-2}\sin^2\phi_k}} d\phi_k \times\int\limits^{t}_{z=0} \int\limits_{z'=0}^{t} \frac{(1-\frac{z}{z_0})(1-\frac{z'}{z_0} )}{\sqrt{F(z,z',\phi_k,x)}}dzdz'$$
%
% where 
%
% $$ F(z,z',\phi_k,x) = (z' - z)^2(\cos^2\phi_k + \beta^{-2}\sin^2\phi_k) + a^2((1 - z/z_0)^2+(1 - z'/z_0)^2)
%                       - 2a^2(1 - z/z_0)(1 - z'/z_0)\cos x $$
%
% $$N_{yy}= \frac{ab}{2\pi V \beta^2} \int\limits_{x=0}^{\pi}\cos x dx  \int\limits_{\phi_k = 0}^{2\pi} \frac{\sin^2\phi_k} {\sqrt{\cos^2\phi_k + \beta^{-2}\sin^2\phi_k}} d\phi_k \times\int\limits^{t}_{z=0} \int\limits_{z'=0}^{t} \frac{(1-\frac{z}{z_0})(1-\frac{z'}{z_0} )}{\sqrt{F(z,z',\phi_k,x)}}dzdz'$$
%
% $$N_{zz} = 1 - N_{xx} - N_{yy}$$
%
% Where the surface of the truncated elliptic cone is described by
%
% $$ \frac{x^2}{a^2} + \frac{y^2}{b^2} = (1 - z/z_0)^2 $$
%
% Where $$z_o = \frac{t}{1 - a_t/a}$, $$0 <= z <= t, t<= z_0$,
% $$-a<=x<=a$, $$-b<=y<=b$ and $$\beta = b/a$. $$t$ is the height of the truncated cone, $$a_t$ is the
% semi-major axis at the top of the cone, $$a$ is the semi-major axis at the
% bottom, $$b$ is the semi-minor axis, and $$x$, $$y$, $$z$ are the axes.
%

%% Choose which shape to use
% There are three possible shapes that can have their demagnetisation
% factors determined by this function. Truncated elliptic cone, elliptic
% cylinder and prism. The dimensions of each must be specified in nm, to
% use this function the geometry must be specified as 1 = truncated
% elliptic cone, 2 = elliptic cylinder, or 3 = prism.
function [ nxx nyy nzz ] = demagfactors(a, b, t, alpha, islandgeo, demag_tol)
switch islandgeo
    case 1
%         disp('truncated elliptic cone')
        [ nxx nyy nzz ] = conedemagIntegral(a, b, t, alpha, demag_tol);
    case 2
%         disp('elliptic cylinder')
        [ nxx nyy nzz ] = cyldemagIntegral(a, b, t, demag_tol);
    case 3
%         disp('prism')
        [nxx nyy nzz] = prismdemagIntegral(a, b, t/2); % since c =t/2
    otherwise
        msg = 'Unknown geometry. Island geometry must be specified as 1 to 3';
        error(msg)
end
end

%% Truncated elliptic cone demagnetising factors 
% Uses a piecewise function to calculate the integral for the demagnetising
% factors $$N_{xx}$ , $$N_{yy}$ , $$N_{zz}$ for a truncated elliptic
% cone. 
function [ nxx nyy nzz ] = conedemagIntegral(a, b, t, alpha, tol)
% a is the semi-major axis 
% b is the semi-minor axis
% t is the height 
% alpha is the ratio of top to bottom axes, a_t and a.
% tol is the tolerance allowed for the integral
beta = b./a; % ratio of the semi-minor axis, b, and the semi-major axis, a.
V= pi*t.*a.*b.*((alpha + 1/2).^2 + 3/4)/3; % volume of the truncated cone
zo = t./(1-alpha); 
%remove constants from integral
cxx = a.*b.*zo./(2*pi*V); 
cyy = a.*b.*zo./(2*pi*V.*beta.^2);
czz = a.*b.*zo./(2*pi*V);

%max and min values for integration
phikmin = 0; phikmax = 2*pi;
xmax = pi; xmin = 0.0000005; % 0

%calculate double integrals over the entire shape
for j =1:length(a)
    Ixx(j) = dblquad(@integrndconexx,phikmin,phikmax,xmin,xmax,tol,[],a(j),b(j),t(j),alpha(j));
    Iyy(j) = dblquad(@integrndconeyy,phikmin,phikmax,xmin,xmax,tol,[],a(j),b(j),t(j),alpha(j));
    Izz(j) = dblquad(@integrndconezz,phikmin,phikmax,xmin,xmax,tol,[],a(j),b(j),t(j),alpha(j));

    nxxp(j) = cxx(j).*Ixx(j);
    nyyp(j) = cyy(j).*Iyy(j);
    nzzp(j) = czz(j).*Izz(j);
end
nzz = 1 -nzzp;
nyy = nyyp;
nxx = nxxp;
end

function out = integrndconexx(phik, x, a, b, t, alpha)
%calculates values for the integral over \phi_k for the N_{xx} component
beta = b./a;

for j = 1:length(phik)
zzpI(j) = zzpIntegral(phik(j), x, a, b, t, alpha);
end
 
out = cos(x).*cos(phik).^2./sqrt(cos(phik).^2 + sin(phik).^2./beta.^2 ).*zzpI;
end

function out = integrndconeyy(phik, x, a, b, t, alpha)
%calculates values for the integral over \phi_k for the N_{yy} component
beta = b./a;

for j = 1:length(phik)
zzpI(j) = zzpIntegral(phik(j), x, a, b, t, alpha);
end
 
out = cos(x).*sin(phik).^2./sqrt(cos(phik).^2 + sin(phik).^2./beta.^2 ).*zzpI;
end

function out = integrndconezz(phik, x, a, b, t, alpha)
%calculates values for the integral over \phi_k for the N_{zz} component
beta = b./a;

for j = 1:length(phik)
zzpI(j) = zzpIntegral(phik(j), x, a, b, t, alpha);
end
 
out = cos(x).*sqrt(cos(phik).^2 + sin(phik).^2./beta.^2 ).*zzpI;
end

function out = zzpIntegral(phik, x, al, bl, t, alpha)
%calculates a piecewise function of the integral for a truncated elliptic
%cone that can be integrated
beta = bl./al;
zo = t./(1-alpha);

Iphik = cos(phik).^2 + 1./beta.^2.*sin(phik).^2;
A=(al./zo).^2 + Iphik;
B=(al./zo).^2.*cos(x) + Iphik;
z = 1-t./zo;

I1t1 = 6.*log(A + sqrt(2*(A -B) ).*sqrt(A) -B) -6*log(A + sqrt(A).*sqrt(A.*z.^2 -2*B.*z + A) - B.*z).*z.^3;
I1t2 = -2 + 2*z.^3 - 3*(A.^2 -3.*B.^2)./A.^2.*log(-B + A + sqrt(A).*sqrt(2*(A - B))) + 3*(A.^2 - 3*B.^2)./A.^2.*log(-B + A.*z + sqrt(A).*sqrt(A.*z.^2 -2*B.*z + A) );
I1t3 = 3*(3*B + A).*sqrt(2*(A-B))./A.^(3/2) - 3*(3*B + A.*z).*sqrt(A.*z.^2 -2*B.*z + A)./A.^(3/2);
I1 = B./(18*A.^(3/2)).*( I1t1 +I1t2 + I1t3 );

I2t1 = -3*( A.^2 -3*B.^2)./A.^2.*log(-B.*z + A + sqrt(A).*sqrt(A.*(z.^2 +1) -2*B.*z)).*z.^3;
I2t2 = 3*(A.^2 - 3*B.^2)./A.^2.*log( -B.*z + A.*z + sqrt(A).*sqrt(2*A.*z.^2 - 2*B.*z.^2 )).*z.^3;
I2t3 = 3*( 3*B.*z + A ).*sqrt(A.*(z.^2 +1) -2*B.*z).*z./A.^(3/2) -2;
I2t4 = -3*(3*B.*z +A.*z ).*sqrt(2*(A - B).*z.^2 ).*z./A.^(3/2) + 2*z.^3;
I2t5 = 6*log( A.*z -B + sqrt(A).*sqrt(A.*(z.^2 +1) -2*B.*z));
I2t6 = -6*z.^3.*log( (A -B).*z + sqrt(A).*sqrt(2*(A-B).*z.^2) );
I2 = -B./(18*A.^(3/2)).*( I2t1 +I2t2 + I2t3 +I2t4 +I2t5 + I2t6 );

I3t1 = sqrt(A).*sqrt(2*(A-B)).*(4*A.^2 -B.*A -3*B.^2) + 3*B.*(A.^2 - B.^2).*log(-B + A + sqrt(A).*sqrt(2*(A-B)));
I3t2 = -sqrt(A).*sqrt(A.*z.^2 -2*B.*z + A).*( 2*(z.^2 +1).*A.^2 -B.*A.*z -3*B.^2 );
I3t3 = -3*B.*( A.^2 -B.^2).*log( -B + A.*z + sqrt(A).*sqrt(A.*z.^2 -2*B.*z + A) );
I3 =  1./(6*A.^(7/2)).*(I3t1 + I3t2 + I3t3);

I4t1 = 3*B.*(A.^2 -B.^2).*log(-B.*z + A + sqrt(A).*sqrt( A.*(z.^2 +1) - 2*B.*z )).*z.^3;
I4t2 = -3*B.*(A.^2 -B.^2).*log((A-B).*z + sqrt(A).*sqrt( 2*(A -B).*z.^2 )).*z.^3;
I4t3 = sqrt(A).*sqrt(A.*(z.^2 + 1) -2*B.*z ).*( 2*(z.^2 +1).*A.^2 -B.*A.*z -3*B.^2.*z.^2 );
I4t4 = -sqrt(A).*sqrt( 2*(A -B).*z.^2 ).*( 4*z.^2.*A.^2 -B.*A.*z.^2 -3*B.^2.*z.^2 );
I4 = -1./(6*A.^(7/2)).*(I4t1 + I4t2 + I4t3 + I4t4);

out = I1 + I2 + I3 + I4; 
end

%% Elliptic cylder demagnetising factor
% The calculation for the demagnetising factors for a elliptic cylinder are
% much the same as for a truncated elliptic cone, but inbuilt matlab
% functions can be used to calculate the integral over the shape. 
function [ nxx nyy nzz ] = cyldemagIntegral(a, b, t, tol)
% a is the semi-major axis 
% b is the semi-minor axis
% t is the height 
% tol is the tolerance allowed for the integral
beta = b./a; % ratio of the semi-minor axis, b, and the semi-major axis, a.
V= pi*t.*a.*b; % volume of the truncated cone

% separate constants from integral
cxx = a.*b./(2*pi*V);
cyy = a.*b./(2*pi*V.*beta.^2);
czz = a.*b./(2*pi*V);

% set limits of integration
phikmin = 0; phikmax = 2*pi;

%calculates double integral over shape
for j =1:length(a)
    Ixx(j) = quad(@integrndcylxx,phikmin,phikmax,tol,[],a(j),b(j),t(j));
    Iyy(j) = quad(@integrndcylyy,phikmin,phikmax,tol,[],a(j),b(j),t(j));
    Izz(j) = quad(@integrndcylzz,phikmin,phikmax,tol,[],a(j),b(j),t(j));

    nxxp(j) = cxx(j).*Ixx(j);
    nyyp(j) = cyy(j).*Iyy(j);
    nzzp(j) = czz(j).*Izz(j);
end
nzz = 1 -nzzp;
nyy = nyyp;
nxx = nxxp;
end
function out = integrndcylxx(phik, a, b,t)
%calculates values for the integral over \phi_k for the N_{xx} component
beta = b./a;
A = cos(phik).^2 + sin(phik).^2./beta.^2;

%return elliptic integral for given function
[elipk,elipe]=ellipke(4*a.^2./(A.*t.^2 +4*a.^2));
I = -8.*a./(3.*sqrt(A)) + t.*(4/3 -A.*t.^2./(3*a.^2)).*sqrt((A.*t.^2 + 4*a.^2)./(A.*t.^2)).*elipe + t.*(4/3 + A.*t.^2./(3*a.^2)).*sqrt(A.*t.^2./(A.*t.^2 + 4*a.^2)).*elipk;
out = cos(phik).^2./A.*I;
end
function out = integrndcylyy(phik, a, b,t)
%calculates values for the integral over \phi_k for the N_{yy} component
beta = b./a;
A = cos(phik).^2 + sin(phik).^2./beta.^2;

%return elliptic integral for given function
[elipk,elipe]=ellipke(4*a.^2./(A.*t.^2 +4*a.^2));
I = -8.*a./(3.*sqrt(A)) + t.*(4/3 -A.*t.^2./(3*a.^2)).*sqrt((A.*t.^2 + 4*a.^2)./(A.*t.^2)).*elipe + t.*(4/3 + A.*t.^2./(3*a.^2)).*sqrt(A.*t.^2./(A.*t.^2 + 4*a.^2)).*elipk;
out = sin(phik).^2./A.*I;
end
function out = integrndcylzz(phik, a, b,t)
%calculates values for the integral over \phi_k for the N_{zz} component
beta = b./a;
A = cos(phik).^2 + sin(phik).^2./beta.^2;

%return elliptic integral for given function
[elipk,elipe]=ellipke(4*a.^2./(A.*t.^2 +4*a.^2));
I = -8.*a./(3.*sqrt(A)) + t.*(4/3 -A.*t.^2./(3*a.^2)).*sqrt((A.*t.^2 + 4*a.^2)./(A.*t.^2)).*elipe + t.*(4/3 + A.*t.^2./(3*a.^2)).*sqrt(A.*t.^2./(A.*t.^2 + 4*a.^2)).*elipk;
out = I;
end


%% Prism demagnetising factors
% Calculates the demagnetising factors for a prism geometry. 
function [nxx nyy nzz] = prismdemagIntegral(a, b, c)
nzz = prismdemag(a, b, c);
nyy = prismdemag(c, a, b);
nxx = prismdemag(b, c, a);

end
function pdemag = prismdemag(a, b, c)
% Uses a piecewise function to determine the integral of the prism geometry
% within one function. 
t1 = log(((a.^2 + b.^2 + c.^2).^(1/2) - a)./((a.^2 + b.^2 + c.^2).^(1/2) + a) ).*(b.^2 - c.^2)./(2*b.*c);
t2 = log(((a.^2 + b.^2 + c.^2).^(1/2) - b)./((a.^2 + b.^2 + c.^2).^(1/2) + b) ).*(a.^2 - c.^2)./(2*a.*c);
t3 = log(((a.^2 + b.^2).^(1/2) + a)./((a.^2 + b.^2).^(1/2) - a) ).*b./(2*c);
t4 = log(((a.^2 + b.^2).^(1/2) + b)./((a.^2 + b.^2).^(1/2) - b) ).*a./(2*c);
t5 = log(((b.^2 + c.^2).^(1/2) - b)./((b.^2 + c.^2).^(1/2) + b) ).*c./(2*a);
t6 = log(((a.^2 + c.^2).^(1/2) - a)./((a.^2 + c.^2).^(1/2) + a) ).*c./(2*b);
t7 = 2*atan(a.*b./(c.*(a.^2 + b.^2 + c.^2).^(1/2)));
t8 = (a.^3 + b.^3 - 2*c.^3)./(3*a.*b.*c) + (a.^2 + b.^2 - 2*c.^2).*((a.^2 + b.^2 + c.^2).^(1/2))./(3*a.*b.*c);
t9 = ((a.^2 + c.^2).^(1/2) + (b.^2 + c.^2).^(1/2)).*c./(a.*b);
t10 = -((a.^2 + b.^2).^(3/2) + (b.^2 + c.^2).^(3/2) + (c.^2 + a.^2).^(3/2) )./(3*a.*b.*c);
pdemag = (t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10)/pi;
end
    
%% References
% [1] Kalezhi, Josephat, Jim J. Miles, and Branson D. Belle. "Dependence of
% switching fields on island shape in bit patterned media." Magnetics, IEEE
% Transactions on 45.10 (2009): 3531-3534.
%
% [2] M. Beleggia and M. De Graef, “On the computation of the
% demagnetization tensor field for an arbitrary particle shape using a
% Fourier space approach,”J. Magn. Magn. Mater., vol. 263, pp. L1–L9, Jul.
% 2003.
