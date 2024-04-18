function [h2_array, lag_array] = ins_findmaxh2(H, tocompare)
%function to create directed weighted connection matrix
if nargin <2
    tocompare=0;
end
nbrchan = size(H.aw_h2, 1);
nbrwin = size(H.aw_h2, 3);

h2_array=zeros(nbrchan,nbrchan,nbrwin);
lag_array=zeros(nbrchan,nbrchan,nbrwin);
for i=1:nbrchan
    for j=1:nbrchan
        if i>j
            if tocompare == 1
                h2_array(i,j,:) = max(squeeze(H.aw_h2(i,j,:)),squeeze(H.aw_h2(j,i,:)));
            else
                for t=1:nbrwin
                    if H.aw_h2(i,j,t) > H.aw_h2(j,i,t)
                        if H.aw_lag(i,j,t) < 0
                            h2_array(i,j,t) = H.aw_h2(i,j,t);
                        else
                            h2_array(j,i,t) = H.aw_h2(i,j,t);
                        end
                    elseif H.aw_h2(j,i,t) > H.aw_h2(i,j,t)
                        if H.aw_lag(j,i,t) <0
                            h2_array(j,i, t) = H.aw_h2(j,i,t);
                        else
                            h2_array(i,j,t) = H.aw_h2(j,i,t);
                        end
                    elseif H.aw_h2(j,i,t) == H.aw_h2(i,j,t)
                        h2_array(i,j,t) = H.aw_h2(i,j,t);
                        h2_array(j,i,t) = H.aw_h2(i,j,t);
                    end
                end
            end
        end
        
    end
end

%%Method Aude with the delay
% h2val = max(H.aw_h2(i,j,t),H.aw_h2(j,i,t));
% delta = sign(H.aw_lag(i,j,t) - H.aw_lag(j,i,t));
% if delta == 0
%     h2_array(j,i, t) = h2val;
%     lag_array(j,i, t) = max(abs(H.aw_lag(i,j,t)), abs(H.aw_lag(j,i,t)));
% elseif delta > 0
%     h2_array(j,i, t) = h2val;
%     lag_array(j,i, t) = max(abs(H.aw_lag(i,j,t)), abs(H.aw_lag(j,i,t)));
% else
%     h2_array(i,j,t) = h2val;
%     lag_array(i,j,t) = max(abs(H.aw_lag(i,j,t)), abs(H.aw_lag(j,i,t)));
% end
% end

%si o mettre des deux côté
%regarder le signe du lag correspondant au H2max
%si R2 pas de normalisation tril
%si h2 basé sur direction ou non