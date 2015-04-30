clf
clc
clear all
hold on
%%  Establish Independent Variable Range
r1=0;
r2=4;

%%
n=500;
for r=linspace(r1,r2,n)
%%  Establish Starting Point
    x=.1;
%%  Converge On Fixed Points
    for iter=1:100
    x=r*x*(1-x);
    end
%%  Plot Fixed Points
    for iter=1:50
    x=r*x*(1-x);
    plot(r,x)
    end
end
%%  Format Plot
xlabel('r')
ylabel('${\bar{x}}$','interpreter','latex')
xlim([0 4])
title('Homework 1 Problem 1')

