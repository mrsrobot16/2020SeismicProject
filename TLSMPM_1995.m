function [freq,amp,varargout] = matrixpencil(x,cut_off)

% function [freq,amp,flag,relres] = matrixpencil(x,cut_off) outputs the -->
% convergence "flag" and residual norm ("relres") as found from the -->
% least squares soltuion for amplitudes using the LSQR function.


% -- Total Least Squares Matrix Pencil Method --
% ------ Implementation by: Sara Makboul -------

% Based on formulations in Paper:
% Sarkar, T. K., & Pereira, O. (1995).
% Using the matrix pencil method to estimate the parameters
% of a sum of complex exponentials. IEEE Antennas and 
% Propagation Magazine, 37(1), 48–55. doi:10.1109/74.370583
% link: https://sci-hub.tw/10.1109/74.370583


% HOW TO:
% 0) Call on the current executing function in the command prompt:
% 0.5) "[freq,amp,varargout] = matrixpencil(x,cut_off)" and include:
% 1) input signal as vector "x"
% 2) input p (significant digits) for noise removal as the "cut_off" 


% NOTE:
% + The damping factor is ignored and not included in signal poles


%%

% checking # of output parameters specified in executing function call, -->
% "if _ OR _ " statement:

if (nargout > 4) | (nargout < 2)
    error('improper number of output parameters');
    return;
end


% setting default value (p=5) for "cut_off" if not specificed in the input:

if nargin == 1
    cut_off = 5;
end

%%

% setting "N": number of equally spaced points, as the length of the x -->
% signal vector:

N = length(x);

% setting "L": pencil parameter (anywhere between N/3 to N/2)
% ceil() rounds up to next larget integer:

L = ceil(N/3);

% OR (for larger L, use):

% L = (N/2) - 1; 

%%

% Hankel data matrix Y:
% ( N - L ) X ( L + 1)
% Hankel matrix: symmetric and constant across the anti-diagonals
% Hankel(c,r) makes first column: c and last row: r
% doesn't add the "__ -1"

Y = x(hankel(1:N-L, N-L:N));

%%

% Hankel Singular Value Decomposition (HSVD)

% Performs an SVD on the Hankel matrix Y

% Then, returns the (N-L) largest singular values = #.
% Note: S represents the sigma E diagonal matrix with the singular values.

[U,S,V] = svds(Y,N-L);

%%

% For M determination, and simultaneous denoising

% "diag()" returns a column vector of the main diagonal elements of S (E)
% these are the singular values:

singval = diag(S);

% The Max singular value is at the top (descending order):

maxsingval = singval(1);

% Assessing ratio of each singular value (in the column vector) to the
% largest singular value ("maxsingval")
% Then, uses the negative of the # of significant decimal digits in the 
% data (p), set as the "cut_off" during input, as the 'tolerance factor'.
% And if the if statement is true, the singular value is NOISE.

for indx = 2:L
    M = indx;
    if log10(abs(singval(indx)/maxsingval)) < -cut_off
        break;
    end
end

%%

% Finding [E'] (sigma), [V'], and [U']:
% [E'] has only M colums corresponding to M dominant singular values


Sp(1:M,1:M) = S(1:M,1:M);
clear S;

Vp = V(:,1:M);
clear V;

Up = U(:,1:M);
clear U;

%%

% The size of Vp is computed, and presented in separate variables: c and r:
% c = columns; r = rows

[c,r] = size(Vp)


% Removing last row from [V'] to create [V1']:

Vp1T = Vp(1:r-1,:);

% Removing first row from [V'] to create [V2']:

Vp2T = Vp(2:r,:);


% [Y1] and [Y2]
% created after SVD, from SVD's noise-filtered, and reconstructed matrices.
% Complex conjugate transpose (^H) to [V1'] and [V2'] denoted with (') 
% is applied to (Vp1T) and (Vp2T):

Y1 = Up*Sp*(Vp1T)';

Y2 = Up*Sp*(Vp2T)';

clear Up;
clear Sp;
clear Vp;
