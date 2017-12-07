function Ind = PeriodsAboveThr(signal,thr,period)
% PeriodsAboveThr finds all indices of consecutive segments of a signal which are
% greater or equal than thr and are equal or longer than period in sampling
% points

if size(signal,1) > size(signal,2)
    signal = signal';
end
if size(signal) ~= 1
    disp('input signal contains more than two channels')
    Ind = {};
    return
end

if ~isinteger(period)
    period = round(period);
end

Signal = signal(1,:) >= thr;
Signal = double(Signal);
Filter = ones(1,period);
Result = conv(Signal,Filter);
Indices = find(Result >= period);

if isempty(Indices)
    disp('No periods matching the criteria')
    Ind = {};
    return
end

if numel(Indices) == 1
    Ind{1} = Indices;
    D = zeros(1,period-1);
    for l=1:length(Ind)
        Ind{l} = horzcat(D,Ind{l});
        for k = 1:length(D)
            Ind{l}(1,k) = Ind{l}(1,length(D)+1)-length(D)+k-1;
        end
    end
    return
end

A = diff(Indices);
B = find(A > 1);

if isempty(B) == 1
    Ind{1} = Indices;
    D = zeros(1,period-1);
    for l=1:length(Ind)
        Ind{l} = horzcat(D,Ind{l});
        for k = 1:length(D)
            Ind{l}(1,k) = Ind{l}(1,length(D)+1)-length(D)+k-1;
        end
    end
    return
end

Dummy = [1 B+1];

Ind = cell(length(Dummy),1);
for l = 1:length(Dummy)-1
    Ind{l} = Indices(Dummy(l):Dummy(l+1)-1);
end
for l = length(Dummy)
    Ind{l} = Indices(Dummy(l):length(Indices));
end

D = zeros(1,period-1);
for l=1:length(Ind)
    Ind{l} = horzcat(D,Ind{l});
    for k = 1:length(D)
        Ind{l}(1,k) = Ind{l}(1,length(D)+1)-length(D)+k-1;
    end
end
















