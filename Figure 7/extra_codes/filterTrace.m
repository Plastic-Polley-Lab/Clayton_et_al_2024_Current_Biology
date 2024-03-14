function data = filterTrace(trace)
trace = double(trace);

Fs = 24414.06;
LowCutoff = 300;
[B,A] = butter(3,LowCutoff/(Fs/2),'high');
HighCutoff = 3000;
[D,C] = butter(2,HighCutoff/(Fs/2),'low');

data = filter(B,A,trace);
data = filter(D,C,data);
data = double(data);
end