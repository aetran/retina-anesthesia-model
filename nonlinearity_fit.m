function [fun_nonlinear, pi_nonlinear] = nonlinearity_fit(method_nonlinear, xdata, ydata)

    switch method_nonlinear
        case 'sigmoid'
            fun_nonlinear = @(x, xdata) x(1)./(1+exp(-x(2)*(xdata-x(3))));
            x0 = [200, 1.5, 2];
        case 'relu'
            fun_nonlinear = @(x, xdata) ( max(x(1)*xdata-x(2), 0) );
            x0 = [1, 1];
    end

    if nargin == 3
        pi_nonlinear = lsqcurvefit(fun_nonlinear, x0, xdata, ydata);
    else
        pi_nonlinear = [];
    end

end