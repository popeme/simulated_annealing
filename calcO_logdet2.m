function O = calcO_logdet2(COV)

% ------------------------------------------------------------------------------
% COV:      covariance matrix of X
% O:        O-information
% I:        integration (total correlation)
% D:        dual total correlation
% uses logdet formalism for numeric stability
% 2021
% ------------------------------------------------------------------------------

N = size(COV,1);
pie1 = 2*pi*exp(1);

% calculate H at level N
%H_n = 0.5*log(det(COV)*pie1^N);
H_n = logdet(COV,'chol')/2 + N*log(pie1)/2;

% calculate H as sum over individual elements
H_xi = 0.5*sum(log(pie1*diag(COV)));

% calculate H at level N-1
H_n1 = zeros(1,N);
for i=1:N
   vv = ones(1,N);
   vv(i) = 0;
   [~,b,~] = find(vv==1);
   %H_n1(i) = 0.5*log(det(COV(b,b))*pie1^(N-1));
   H_n1(i) = logdet(COV(b,b),'chol')/2 + (N-1)*log(pie1)/2;
end
% sum of H at level N-1
H_n1 = sum(H_n1);

% calculate o-information, integration (total correlation) and dual total correlation
O = (N-2)*H_n + H_xi - H_n1;
I = H_xi - H_n;
D = I - O;
