function V_gamma=V1(N,set,gamma,g_act)

V_gamma=zeros(1,N);

for i=1:length(set)
    if gamma(i)==0
        V_gamma(set(i))=abs(min(max(gamma(i) - g_act(i), 0),inf) - gamma(i));
    elseif gamma(i)==1
        V_gamma(set(i))=abs(min(max(gamma(i) - g_act(i), -inf),1) - gamma(i));
    else
        V_gamma(set(i))=abs(min(max(gamma(i) - g_act(i), -inf),inf) - gamma(i));
    end
end


end