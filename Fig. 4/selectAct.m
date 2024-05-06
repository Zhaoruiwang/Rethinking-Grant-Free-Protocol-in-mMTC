function act_set = selectAct(gamma, g, k_iter, epsilon)

    N=length(gamma);
    % active set
    grad0 = g .* (gamma == 0);
    grad1 = g .* (gamma == 1);
    
    act_set0 = find( gamma == 0 & g < -0.8^(k_iter+2)*abs(max(grad0)) );
    act_set1 = find( gamma == 1 & g > 0.8^(k_iter+2)*abs(max(grad1)) );
    act_setm = find( gamma ~= 0 & gamma ~= 1 & abs(g) > max(0.1^k_iter, epsilon/sqrt(0.3*N)) );

    act_set = [act_set0; act_set1; act_setm];
end