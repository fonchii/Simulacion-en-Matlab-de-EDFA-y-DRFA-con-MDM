function h_pcc = pcc_calc(beta,modos)
% beta: constante de propagación [1/m]
% modos: string de modos

Pmod=length(modos);
Kuv=0.01; % [1/m] field coupling coeficient, default value in VPI 0.01 m^-1
dcor=5e-2; % [m] coupling correlation length, default value VPI 5e-2
h_pcc=zeros(length(modos));

for i=1:length(beta)
    for u=1:Pmod
        betau=beta.(modos(u)){i};
        for v=1:Pmod
            betav=beta.(modos(v)){i};
            if (u~=v)
                h_pcc(u,v)=Kuv^2*(dcor/(1+(dcor*abs(betau-betav))^2)); %[1/m]
            end
        end
    end
end

% aux=find(DeltaBetauv>0);
% DeltaBref=[min(min(DeltaBetauv(aux))) max(max(DeltaBetauv(aux)))];
% 
% dcorv=linspace(2/DeltaBref(1),2/DeltaBref(2),10);
% 
% for j=1:length(dcorv)
%     h_pcc(j)=Kuv^2*dcorv(j)./(1+(dcorv(j)*DeltaBetauv(1,2)).^2);
% end


end