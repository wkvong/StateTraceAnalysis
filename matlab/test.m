A = [1 1 0; 0 0 0; 0 1 1];

nodes = [1 2 3];

[x, fit, exitflag, output] = MR1([2 4 8], eye(3), A);

% x = dlmread('../data/x.dat', '\t');

% output = staSTATS(x);
% [x fit exitflag] = staMR(x);

% output = {output};

% [xStar Ffits EStar exitflag] = CMRv4(output, {[3 2 1] [4 3]})

% [x f] = staCMR(data{1})

% load '../data/nakabayashi.mat'


% [p datafit fits] = CMRfits(100, data);
