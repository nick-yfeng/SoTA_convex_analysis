function [y_seg, N_seg, IsMA] = segmentation(y,s1)
N = length(s1);

y_seg = {};
N_seg = [];


s1(s1 ~= 0) = 1;
k = 1;

if s1(1) == 1
    s1 = ~s1;
    IsMA = 1;

    for i = 1:N
        

        if sum(s1((1+sum(N_seg)):i)) == 1
            y_seg{k} = y((1+sum(N_seg)):i-1);
            N_seg(k) = i-1-sum(N_seg);
            k = k + 1;
           
            s1 = ~s1;
        end
        
    end
    y_seg{k} = y((1+sum(N_seg)):N);
    N_seg(k) = N-sum(N_seg);
    
else
    IsMA = 0;
    
    for i = 1:N
        

        if sum(s1((1+sum(N_seg)):i)) == 1
            y_seg{k} = y((1+sum(N_seg)):i-1);
            N_seg(k) = i-1-sum(N_seg);
            k = k + 1;
           
            s1 = ~s1;
        end

    end
    y_seg{k} = y((1+sum(N_seg)):N);
    N_seg(k) = N-sum(N_seg);

end



end