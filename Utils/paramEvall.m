function H = paramEval(A,th)
    if isa(th,"single")||isa(th,"double")
        H=permute(pagemtimes(th,permute(A,[3 1 2])),[2 3 1]);
    elseif isa(th,'sdpvar')
        H=0;
        for i=1:length(th)
            H=H+th(i)*A(:,:,i);
        end
    end
end

