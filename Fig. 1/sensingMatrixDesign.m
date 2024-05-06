function A = sensingMatrixDesign(N,J,L,type,mode)
% Sequences
if mode==1
    Ne=N;
end
if mode==2
    Ne=N*2^J;
end
if strcmp(type, 'Gaussian')
    A = (randn(L,Ne) + 1i*randn(L,Ne))*sqrt(0.5);


elseif strcmp(type, 'QAM')
    A = randi(4,L,Ne);
    A(A == 1) = 1+1j;
    A(A == 2) = 1-1j;
    A(A == 3) = -1+1j;
    A(A == 4) = -1-1j;
    A = A*sqrt(0.5);
    
else
    error('oooooops matrix');
end