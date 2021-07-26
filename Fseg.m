% FUNCTION Fseg(w)

function Fw = Fseg(w)

for i = 1 : length(w)

    %if w(i) > 1
    %    error( 'w(i) > 1 !!!' )
    %end
    
    if w(i) < 0.005
        Fw(i) = 600*w(i);
    else
        Fw(i) = 3;
    end
    
end