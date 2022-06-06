function [mse] = MSE(input1,input2)
    %assert (size(input1) == size(input2));
    mse=0;
    for i = 1: length(input1)
        mse = mse + (abs(input1(i)-input2(i)))^2/length(input1);
    end

    mse=sqrt(mse);
end