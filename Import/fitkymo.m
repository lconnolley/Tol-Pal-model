function [x,fval,x2,fval2]=fitkymo(t,data,guess)

%Fickian diffusion, d is the diffusion constant
fun=@(d) cost(data(:,2:end), spatialFRAP(t,data(:,2),d));
options=optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'InitialMeshSize',0.1);
[x,fval,exitflag,output] = patternsearch(fun,guess,[],[],[],[],0,[],[],options);


%Fokker-Planck diffusion, d a the proportionality constant
fun=@(d) cost(data(:,2:end), spatialFRAP(t,data(:,2),d./data(:,1)/length(data(:,1))));
%fun=@(d) JS(data(:,2:end), spatialFRAP(data(:,2),d./data(:,1)/length(data(:,1))));
options=optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'InitialMeshSize',0.1);
[x2,fval2,exitflag,output] = patternsearch(fun,guess,[],[],[],[],0,[],[],options);

%---------------------------------------------------------------

function out=JS(A,B)
         out=0;
         for i=1:size(A,2)
            out=out+sqrt(JSdivergence(A(:,i),B(:,i)));%data is normalised so as a pdf i.e. sum=1
         end
end

function out=cost(A,B)
    
         %out=1e5*immse(A,B);
         out=10*immsre(A,B);
         %out=1e2*JS(A,B);
end


end