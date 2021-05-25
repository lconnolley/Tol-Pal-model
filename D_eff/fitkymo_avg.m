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
    
         out=1e5*immse(A,B);
         %out=10*immsre(A,B);
         %out=1e2*JS(A,B);
end

function err = immsre(x, y)
%IMMSRE Mean-Squared RELATIVE Error
%based on matlabs immse
 

validateattributes(x,{'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', ...
    'single','double'},{'nonsparse'},mfilename,'A',1);
validateattributes(y,{'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', ...
    'single','double'},{'nonsparse'},mfilename,'B',1);

if ~isa(x,class(y))
    error(message('images:validate:differentClassMatrices','A','B'));
end
    
if ~isequal(size(x),size(y))
    error(message('images:validate:unequalSizeMatrices','A','B'));
end

if isempty(x) % If x is empty, y must also be empty
    err = [];
    return;
end

if isinteger(x)     
    x = double(x);
    y = double(y);
end

err = (norm(1-x(:)./y(:),2).^2)/numel(x);
end

end
