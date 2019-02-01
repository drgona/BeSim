function [Xdelayed] = delayData(X,numDelays,flipFlag)

% Xdelayed = delayData(X,numDelays)
% input data X in the format:
%           rows = samples
%           columns = features
% Outputs in the format:
%           rows = samples
%           columns = shifted features X(t,t-1,t-2,...,t-numDelays)              if flipFlag = false
%           columns = shifted features X(t-numDelays, t-numDelays+1,...,t-1,t)   if flipFlag = true
% uf numDelays = 1 than nothing happen 


% % initialize zero matrix
% Xdelayed = zeros(size(X,1)-numDelays, (numDelays+1)*size(X,2));
% % fill up the matrix with delays
% for i = 1:numDelays+1
%     Xdelayed(:, (i-1)*size(X,2)+1 : i*size(X,2) ) = X(i:end-(numDelays+1)+i,:);
% %       Xdelayed(:, (i-1)*size(X,2)+1 : i*size(X,2) ) = flip(X(i:end-(numDelays+1)+i,:),2);
% 
% end

% initialize zero matrix
Xdelayed = zeros(size(X,1)-(numDelays-1), numDelays*size(X,2));
% fill up the matrix with delays
for i = 1:numDelays
    Xdelayed(:, (i-1)*size(X,2)+1 : i*size(X,2) ) = X(i:end-numDelays+i,:);
%       Xdelayed(:, (i-1)*size(X,2)+1 : i*size(X,2) ) = flip(X(i:end-(numDelays+1)+i,:),2);

end

if flipFlag
    Xdelayed = flip(Xdelayed,1);
end

% figure
% plot(Xdelayed(1:100,[30,60,90]))

end