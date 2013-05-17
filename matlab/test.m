% A = [1 0 1; 0 1 0; 0 1 0];

% nodes = [1 2 3];

% [x, fit, exitflag, output] = MR1([2 4 8], eye(3), A);

% x = dlmread('x.dat', '\t');
% x = num2cell(x, 1)

% output = staSTATS(x);

% [x fit exitflag] = staMR(x);

% output = {output};

% [xStar Ffits EStar exitflag] = CMRv4(output, {[3 2 1] [4 3]})

% [x f] = staCMR(data{1})

load '../data/nakabayashi.mat'

CMRfits(2, data)