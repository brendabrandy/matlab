%% INTRODUCTION TO STOCHASTICITY

% Electrical engineers often use MATLAB to simulate random process, and the
% reason that it is so popular is that a lot of problems in EE are
% stochastic problems. What happens when there's random noise in a channel,
% what is the probability of corrupting a signal to the point that it
% cannot be recovered? 

% Today, we are going to explore the idea of randomness with case studies.
% We are first going to show simple examples of stochastic processes, and
% will conclude with a more elaborate exploration of stochasticity in
% communication channels

%% EXAMPLE 1 : MONTY HALL 
% The classic monty hall problem -- do you switch the door or not. The
% correct answer is you always switch. But why?? Let's find out

% Let's do a simulation to experimentally see the results, let's say there
% are two contestants, the first contestant never switch, the second
% contestant switches. 

% There are different ways to do it. The laziest way to do it is to first
% assume that the contestants always chooses door 1. That way, if the prize
% is behind door 1, the first contestant wins, but if the prize not behind
% door 1.
N = 50000;
prize_door = unidrnd(3,1,N);
result = zeros(2,N);
result(1, prize_door == 1) = 1;
result(2, prize_door ~= 1) = 1;
win_rate = cumsum(result , 2);
x = 1:1:N;
win_rate = win_rate ./ x;

%% Plotting the results of example 1
figure;
plot(x, win_rate(1,:), 'DisplayName','no switch');
hold on;
plot(x, win_rate(2,:), 'DisplayName', 'switch');
legend('show');
title('Monty Hall Problem Simulation Results')
hold off;
%% EXAMPLE 2 : CASINO

% A casino offers a game in which a fair coin is tossed repeatedly until
% the first heads comes up, for a total of k+1 tosses. The casino will pay
% out 2^k dollars

% Let's see how much money you are going to get if you play the game
N = 5000;
num_tosses = zeros(1,N);
for i = 1:N
   res = randi([0,1]);
   while (res ~= 0)
      num_tosses(i) = num_tosses(i) + 1;
      res = randi([0,1]);
   end
end
figure;
hist(num_tosses, unique(num_tosses))

% Okay, now what if instead, the casino cuts off the payout
Pmax = 1:0.1:44;
Pmax = 2.^(Pmax);
payout = ceil(log2(Pmax))./2 + Pmax ./(2.^(ceil(log2(Pmax))));
figure;
plot(Pmax, payout);
title('Expected payout at different payout cutoffs');

%% Radar/Signal Detection Example
%% Generating theoretical signal.
t1 = linspace(0,1,1000);
t = linspace(0,2,2000);
s1 = sin(50*t1 + 300*t1.^2);
s2 = fliplr(s1);
s = [s1 s2];
plot(t,s)
%% Generating the actual signal
N = floor(exprnd(10000,1,10));
sentSignal = [];
for n = N
    sentSignal = [sentSignal s zeros(1,n)];
end
receivedSignal = sentSignal + 3*randn(size(sentSignal)) + 100;
subplot(211)
plot(sentSignal)
subplot(212)
plot(receivedSignal)
%% Finding our signal in the noise
x = conv(receivedSignal,fliplr(s),'valid');
figure
plot(x)

%% QAM Example
%% Set up constellation
n = 4;
z = exp((0:n-1)/4*2*pi*1j);
plot(z,'o')
axis([-2 2 -2 2])
%% Generate random signal (with apriori)
apriori = [0.05 0.20 0.25 0.5]; %Pick aprioris
cmf = [0 cumsum(apriori)];
r = rand(10000,1);
p = zeros(size(r));
accum = 0;
n = 0;
for n = 2:length(cmf)
    p(cmf(n-1) <= r & r < cmf(n)) = n-1;
end
%% Confirm our generation method worked
var = 0.3;
realSignal = z(p) + sqrt(var/2)*randn(size(z(p))) + ...
                    1j*sqrt(var/2)*randn(size(z(p)));
histogram(p,'Normalization','probability')
figure,histogram(p)
figure,plot(realSignal,'.'), hold on
axis([-2 2 -2 2])
%% Decision making (ML)
plot([-2,2],[-2,2],[-2 2],[2 -2],'r')
d = (realSignal - z.');
[~,D] = min(abs(d));
c = p(:) - D(:);
fractionWrong = sum(c ~= 0) / length(p)
%% Decision making (MAP)
d = (realSignal - z.');
[~,D] = max(normpdf(abs(d),0,sqrt(var)).*apriori.');
c = p(:) - D(:);
fractionWrong = sum(c ~= 0) / length(p)



