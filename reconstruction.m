function x = reconstruction(s_seg, N_seg, alpha, beta)

N = sum(N_seg);
x = zeros(floor(N),1);
k = length(N_seg);

x(1:N_seg(1)) = s_seg{1};

for i = 2:k
    
    lam1 = N_seg(i-1);
    lam2 = N_seg(i);
    theta1 = lam1/10;
    theta2 = lam2/10;
    temp1 = s_seg{i-1};
    temp2 = s_seg{i};
    
    if lam1 <= alpha 
        
        if lam2 <= alpha
            
            a = mean(temp1);
            b = mean(temp2);
            
        elseif lam2 >= beta
            
            a = mean(temp1);
            b = mean(temp2(1:theta2));
            
        else
            
            a = mean(temp1);
            b = mean(temp2(1:alpha));
            
        end
                
    elseif lam1 >= beta
        
        if lam2 <= alpha
            
            a = mean(temp1((lam1-alpha):lam1));
            b = mean(temp2);
            
        elseif lam2 >= beta
            
            a = mean(temp1((lam1-alpha):lam1));
            b = mean(temp2(1:theta2));
            
        else
            
            a = mean(temp1((lam1-alpha):lam1));
            b = mean(temp2(1:alpha));
            
        end
         
    else
        
        if lam2 <= alpha
            
            a = mean(temp1((lam1-theta1):lam1));
            b = mean(temp2);
            
        elseif lam2 >= beta
            
            a = mean(temp1((lam1-theta1):lam1));
            b = mean(temp2(1:theta2));
            
        else
            
            a = mean(temp1(floor(lam1-theta1):lam1));
            b = mean(temp2(1:alpha));
            
        end
        
    end
    
    v = a - b;
    s_seg{i} = s_seg{i} + v;
    x((1:N_seg(i)) + sum(N_seg(1:i-1))) = s_seg{i};
    

end



end