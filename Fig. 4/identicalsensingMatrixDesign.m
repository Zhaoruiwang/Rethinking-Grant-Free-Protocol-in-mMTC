function A = identicalsensingMatrixDesign(N,J,L,type,mode)
% Sequences
if mode==1
    Ne=N;
end
if mode==2
    Ne=N*2^J;
end
if strcmp(type, 'Gaussian')
    A = (randn(L,1) + 1i*randn(L,1))*sqrt(0.5);
    A = repmat(A,1,Ne);

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