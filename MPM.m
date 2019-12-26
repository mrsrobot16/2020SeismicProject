function [wk,xk]=MP_fb(y,L,K)
% Matrix-Pencil method, based on the paper:
% Y. Hua and T.K. Sarkar, “Matrix pencil method for estimating 
% parameters of exponentially damped/undamped sinusoids in noise”, 
% Acoustics, Speech and Signal Processing, IEEE Transactions on, vol. 38, 
% no. 5, pp. 814–824, May 1990.


% y - samples
% L - predictor order, must be greather than 2K
% K - number of exponentials

% wk - estimated frequencies
% xk - estimated amplitudes

N=length(y);

if L<2*K
    error('Predictor order L must be greater or equal to 2K');
end

if L>=N
    error('Predictor order must be smaller than number of samples');
end    


% build matrices
Y0=zeros(N-L,L);
Y1=zeros(N-L,L);

for i=0:N-L-1
    Y0(i+1,:)=y((1:L)+i);
    Y1(i+1,:)=y((2:L+1)+i);
end

Y0FB=[Y0; conj(Y1(:,end:-1:1))];
Y1FB=[Y1; conj(Y0(:,end:-1:1))];

Y0=Y0FB;
Y1=Y1FB;

% reduce rank approximation
[U,S,V]=svd(Y1);
Sp=S(1:K,1:K); Up=U(:,1:K); Vp=V(:,1:K);

[U,S,V]=svd(Y0);
S=S(1:K,1:K); U=U(:,1:K); V=V(:,1:K);
Y0=U*S*V';

% G.E.V problem solution
Zl=Vp*inv(Sp)*Up'*Y0;
z=eig(Zl);

% retrieve exponents frequencies
[B,I]=sort(abs(z));
wk=sort(mod((-angle(z(I(end-K+1:end)))),2*pi)/(2*pi));

% recover amplitudes
[W,NN]=meshgrid(wk,0:N-1);
Mat=exp(j*2*pi*W.*NN) ;

xk = pinv(Mat)*y;

end