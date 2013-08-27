function [residual, g1, g2, g3] = GrowthApproximate_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(4, 1);
T14 = (1-y(5))^(1-params(2));
T15 = y(3)^params(2)*T14;
T25 = (1-y(8))^(1-params(2));
T26 = y(7)^params(2)*T25;
T29 = params(1)*T26^(1-params(5))/y(7);
T37 = params(4)*exp(y(9))*y(4)^(params(4)-1);
T39 = y(8)^(1-params(4));
T40 = T37*T39;
T43 = 1+T40-params(3);
T50 = (1-params(4))*params(2)/(1-params(2))*exp(y(6));
T53 = T50*y(1)^params(4);
T55 = y(5)^(-params(4));
T56 = T53*T55;
T59 = exp(y(6))*y(1)^params(4);
T60 = y(5)^(1-params(4));
T61 = T59*T60;
T76 = T14*getPowerDeriv(y(3),params(2),1);
T77 = getPowerDeriv(T15,1-params(5),1);
T78 = T76*T77;
T84 = T25*getPowerDeriv(y(7),params(2),1);
T85 = getPowerDeriv(T26,1-params(5),1);
T86 = T84*T85;
T91 = params(1)*(y(7)*T86-T26^(1-params(5)))/(y(7)*y(7));
T105 = T39*params(4)*exp(y(9))*getPowerDeriv(y(4),params(4)-1,1);
T110 = y(3)^params(2)*(-(getPowerDeriv(1-y(5),1-params(2),1)));
T124 = y(7)^params(2)*(-(getPowerDeriv(1-y(8),1-params(2),1)));
T127 = params(1)*T85*T124/y(7);
T129 = T37*getPowerDeriv(y(8),1-params(4),1);
lhs =T15^(1-params(5))/y(3);
rhs =T29*T43;
residual(1)= lhs-rhs;
lhs =y(3);
rhs =(1-y(5))*T56;
residual(2)= lhs-rhs;
lhs =y(4);
rhs =T61-y(3)+y(1)*(1-params(3));
residual(3)= lhs-rhs;
lhs =y(6);
rhs =params(6)*y(2)+params(7)*x(it_, 1);
residual(4)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(4, 10);

  %
  % Jacobian matrix
  %

  g1(1,3)=(y(3)*T78-T15^(1-params(5)))/(y(3)*y(3));
  g1(1,7)=(-(T43*T91));
  g1(1,4)=(-(T29*T105));
  g1(1,5)=T77*T110/y(3);
  g1(1,8)=(-(T43*T127+T29*T129));
  g1(1,9)=(-(T29*T40));
  g1(2,3)=1;
  g1(2,1)=(-((1-y(5))*T55*T50*getPowerDeriv(y(1),params(4),1)));
  g1(2,5)=(-((-T56)+(1-y(5))*T53*getPowerDeriv(y(5),(-params(4)),1)));
  g1(2,6)=(-((1-y(5))*T56));
  g1(3,3)=1;
  g1(3,1)=(-(1-params(3)+T60*exp(y(6))*getPowerDeriv(y(1),params(4),1)));
  g1(3,4)=1;
  g1(3,5)=(-(T59*getPowerDeriv(y(5),1-params(4),1)));
  g1(3,6)=(-T61);
  g1(4,2)=(-params(6));
  g1(4,6)=1;
  g1(4,10)=(-params(7));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(38,3);
  v2(1,1)=1;
  v2(1,2)=23;
  v2(1,3)=(y(3)*y(3)*(T78+y(3)*(T77*T14*getPowerDeriv(y(3),params(2),2)+T76*T76*getPowerDeriv(T15,1-params(5),2))-T78)-(y(3)*T78-T15^(1-params(5)))*(y(3)+y(3)))/(y(3)*y(3)*y(3)*y(3));
  v2(2,1)=1;
  v2(2,2)=67;
  v2(2,3)=(-(T43*params(1)*(y(7)*y(7)*(T86+y(7)*(T85*T25*getPowerDeriv(y(7),params(2),2)+T84*T84*getPowerDeriv(T26,1-params(5),2))-T86)-(y(7)*T86-T26^(1-params(5)))*(y(7)+y(7)))/(y(7)*y(7)*y(7)*y(7))));
  v2(3,1)=1;
  v2(3,2)=37;
  v2(3,3)=(-(T91*T105));
  v2(4,1)=1;
  v2(4,2)=64;
  v2(4,3)=  v2(3,3);
  v2(5,1)=1;
  v2(5,2)=34;
  v2(5,3)=(-(T29*T39*params(4)*exp(y(9))*getPowerDeriv(y(4),params(4)-1,2)));
  v2(6,1)=1;
  v2(6,2)=43;
  v2(6,3)=(y(3)*(T110*T76*getPowerDeriv(T15,1-params(5),2)+T77*getPowerDeriv(y(3),params(2),1)*(-(getPowerDeriv(1-y(5),1-params(2),1))))-T77*T110)/(y(3)*y(3));
  v2(7,1)=1;
  v2(7,2)=25;
  v2(7,3)=  v2(6,3);
  v2(8,1)=1;
  v2(8,2)=45;
  v2(8,3)=(T110*T110*getPowerDeriv(T15,1-params(5),2)+T77*y(3)^params(2)*getPowerDeriv(1-y(5),1-params(2),2))/y(3);
  v2(9,1)=1;
  v2(9,2)=77;
  v2(9,3)=(-(T43*params(1)*(y(7)*(T124*T84*getPowerDeriv(T26,1-params(5),2)+T85*getPowerDeriv(y(7),params(2),1)*(-(getPowerDeriv(1-y(8),1-params(2),1))))-T85*T124)/(y(7)*y(7))+T91*T129));
  v2(10,1)=1;
  v2(10,2)=68;
  v2(10,3)=  v2(9,3);
  v2(11,1)=1;
  v2(11,2)=74;
  v2(11,3)=(-(T105*T127+T29*params(4)*exp(y(9))*getPowerDeriv(y(4),params(4)-1,1)*getPowerDeriv(y(8),1-params(4),1)));
  v2(12,1)=1;
  v2(12,2)=38;
  v2(12,3)=  v2(11,3);
  v2(13,1)=1;
  v2(13,2)=78;
  v2(13,3)=(-(T127*T129+T43*params(1)*(T124*T124*getPowerDeriv(T26,1-params(5),2)+T85*y(7)^params(2)*getPowerDeriv(1-y(8),1-params(2),2))/y(7)+T127*T129+T29*T37*getPowerDeriv(y(8),1-params(4),2)));
  v2(14,1)=1;
  v2(14,2)=87;
  v2(14,3)=(-(T40*T91));
  v2(15,1)=1;
  v2(15,2)=69;
  v2(15,3)=  v2(14,3);
  v2(16,1)=1;
  v2(16,2)=84;
  v2(16,3)=(-(T29*T105));
  v2(17,1)=1;
  v2(17,2)=39;
  v2(17,3)=  v2(16,3);
  v2(18,1)=1;
  v2(18,2)=88;
  v2(18,3)=(-(T29*T129+T40*T127));
  v2(19,1)=1;
  v2(19,2)=79;
  v2(19,3)=  v2(18,3);
  v2(20,1)=1;
  v2(20,2)=89;
  v2(20,3)=(-(T29*T40));
  v2(21,1)=2;
  v2(21,2)=1;
  v2(21,3)=(-((1-y(5))*T55*T50*getPowerDeriv(y(1),params(4),2)));
  v2(22,1)=2;
  v2(22,2)=41;
  v2(22,3)=(-((-(T55*T50*getPowerDeriv(y(1),params(4),1)))+(1-y(5))*T50*getPowerDeriv(y(1),params(4),1)*getPowerDeriv(y(5),(-params(4)),1)));
  v2(23,1)=2;
  v2(23,2)=5;
  v2(23,3)=  v2(22,3);
  v2(24,1)=2;
  v2(24,2)=45;
  v2(24,3)=(-((-(T53*getPowerDeriv(y(5),(-params(4)),1)))+(-(T53*getPowerDeriv(y(5),(-params(4)),1)))+(1-y(5))*T53*getPowerDeriv(y(5),(-params(4)),2)));
  v2(25,1)=2;
  v2(25,2)=51;
  v2(25,3)=(-((1-y(5))*T55*T50*getPowerDeriv(y(1),params(4),1)));
  v2(26,1)=2;
  v2(26,2)=6;
  v2(26,3)=  v2(25,3);
  v2(27,1)=2;
  v2(27,2)=55;
  v2(27,3)=(-((-T56)+(1-y(5))*T53*getPowerDeriv(y(5),(-params(4)),1)));
  v2(28,1)=2;
  v2(28,2)=46;
  v2(28,3)=  v2(27,3);
  v2(29,1)=2;
  v2(29,2)=56;
  v2(29,3)=(-((1-y(5))*T56));
  v2(30,1)=3;
  v2(30,2)=1;
  v2(30,3)=(-(T60*exp(y(6))*getPowerDeriv(y(1),params(4),2)));
  v2(31,1)=3;
  v2(31,2)=41;
  v2(31,3)=(-(exp(y(6))*getPowerDeriv(y(1),params(4),1)*getPowerDeriv(y(5),1-params(4),1)));
  v2(32,1)=3;
  v2(32,2)=5;
  v2(32,3)=  v2(31,3);
  v2(33,1)=3;
  v2(33,2)=45;
  v2(33,3)=(-(T59*getPowerDeriv(y(5),1-params(4),2)));
  v2(34,1)=3;
  v2(34,2)=51;
  v2(34,3)=(-(T60*exp(y(6))*getPowerDeriv(y(1),params(4),1)));
  v2(35,1)=3;
  v2(35,2)=6;
  v2(35,3)=  v2(34,3);
  v2(36,1)=3;
  v2(36,2)=55;
  v2(36,3)=(-(T59*getPowerDeriv(y(5),1-params(4),1)));
  v2(37,1)=3;
  v2(37,2)=46;
  v2(37,3)=  v2(36,3);
  v2(38,1)=3;
  v2(38,2)=56;
  v2(38,3)=(-T61);
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),4,100);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],4,1000);
end
end
