function [x]=fitkymo(t,data,guess)

%Fokker-Planck diffusion, d a the proportionality constant
fun=@(d) cost(data(:,2:end), spatialFRAP(t,data(:,2),d./data(:,1)/length(data(:,1))));
options=optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'InitialMeshSize',0.1);
[x,~,~,~] = patternsearch(fun,guess,[],[],[],[],0,[],[],options);

%---------------------------------------------------------------

function out=cost(A,B)
         out=1e5*immse(A,B);
        %out=10*immsre(A,B);
end

end