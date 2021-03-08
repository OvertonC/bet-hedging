function compile
N=10000;
X = 0;
I = 0;
for i=1:N
    
    name = ['output',num2str(i),'.mat'];
    try 
        load(name);
        X = X + outputresults;

    catch
        V = 0;
        I = I + 1;
    end
    
end

X=X./(N-I);

save('output_compiled.mat','X');
