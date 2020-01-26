function [freq,amp,varargout] = matrixpencil(x,cut_off)

% function [freq,amp,flag,relres] = matrixpencil(x,cut_off) -->
% --> Outputs the convergence "flag" and residual norm ("relres") -->
% --> as calculated from the least squares soltuion for amplitudes -->
% --> using the LSQR function.



% -- Total Least Squares Matrix Pencil Method --
% ------- Implementation by Sara Makboul -------

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
% + The modal order (M) is represented as "N".
% + There will be "N" values in the Frequencies ("freq") output vector
% + There will be "N" values in the Amplitudes ("amp") output vector
