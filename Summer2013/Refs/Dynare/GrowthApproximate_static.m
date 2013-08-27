function [residual, g1, g2] = GrowthApproximate_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 4, 1);

%
% Model equations
%

T15 = y(1)^params(2)*(1-y(3))^(1-params(2));
T19 = T15^(1-params(5))/y(1);
T31 = y(3)^(1-params(4));
T35 = 1+params(4)*exp(y(4))*y(2)^(params(4)-1)*T31-params(3);
T42 = exp(y(4))*(1-params(4))*params(2)/(1-params(2))*y(2)^params(4);
T45 = T42*y(3)^(-params(4));
T69 = (y(1)*(1-y(3))^(1-params(2))*getPowerDeriv(y(1),params(2),1)*getPowerDeriv(T15,1-params(5),1)-T15^(1-params(5)))/(y(1)*y(1));
T91 = getPowerDeriv(T15,1-params(5),1)*y(1)^params(2)*(-(getPowerDeriv(1-y(3),1-params(2),1)))/y(1);
lhs =T19;
rhs =T19*params(1)*T35;
residual(1)= lhs-rhs;
lhs =y(1);
rhs =(1-y(3))*T45;
residual(2)= lhs-rhs;
lhs =y(2);
rhs =T31*exp(y(4))*y(2)^params(4)-y(1)+y(2)*(1-params(3));
residual(3)= lhs-rhs;
lhs =y(4);
rhs =y(4)*params(6)+params(7)*x(1);
residual(4)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(4, 4);

  %
  % Jacobian matrix
  %

  g1(1,1)=T69-T35*params(1)*T69;
  g1(1,2)=(-(T19*params(1)*T31*params(4)*exp(y(4))*getPowerDeriv(y(2),params(4)-1,1)));
  g1(1,3)=T91-(T35*params(1)*T91+T19*params(1)*params(4)*exp(y(4))*y(2)^(params(4)-1)*getPowerDeriv(y(3),1-params(4),1));
  g1(1,4)=(-(T19*params(1)*params(4)*exp(y(4))*y(2)^(params(4)-1)*T31));
  g1(2,1)=1;
  g1(2,2)=(-((1-y(3))*y(3)^(-params(4))*exp(y(4))*(1-params(4))*params(2)/(1-params(2))*getPowerDeriv(y(2),params(4),1)));
  g1(2,3)=(-((-T45)+(1-y(3))*T42*getPowerDeriv(y(3),(-params(4)),1)));
  g1(2,4)=(-((1-y(3))*T45));
  g1(3,1)=1;
  g1(3,2)=1-(1-params(3)+T31*exp(y(4))*getPowerDeriv(y(2),params(4),1));
  g1(3,3)=(-(exp(y(4))*y(2)^params(4)*getPowerDeriv(y(3),1-params(4),1)));
  g1(3,4)=(-(T31*exp(y(4))*y(2)^params(4)));
  g1(4,4)=1-params(6);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],4,16);
end
end
