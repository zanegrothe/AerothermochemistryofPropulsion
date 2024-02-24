% Zane Grothe
% AERO 7530
% Computer Project
% 12/10/21

% Composition equations *including* N2 in the products

function eq=CompEQs1(x,n_ox,xR_N2,KpO,KpH,KpH2O,KpOH)

eq(1)=n_ox*(1-xR_N2)*(2*x(3)+2*x(2)+x(1)+x(4))-x(3)-x(4)-2*x(6)-x(5);
eq(2)=xR_N2*(x(3)+x(4)+2*x(6)+x(5))-2*(1-xR_N2)*x(7);
eq(3)=x(3)+x(2)+x(1)+x(4)+x(6)+x(5)+x(7)-1;
eq(4)=KpO-x(5)^2/x(6);
eq(5)=KpH-x(1)^2/x(2);
eq(6)=KpH2O-x(2)*sqrt(x(6))/x(3);
eq(7)=KpOH-sqrt(x(2)*x(6))/x(4);

% H=x(1)
% H2=x(2)
% H2O=x(3)
% OH=x(4)
% O=x(5)
% O2=x(6)
% N2=x(7)
