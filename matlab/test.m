y = load('../data/delay.dat');
[x f] = staCMR(y);
staPLOT(y,x,{1:4 5:8},{'No Delay' 'Delay'},{'RB' 'II'},{[.3 .8] [.2 ...
                    .7]});
