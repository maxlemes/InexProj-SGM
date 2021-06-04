function [A,B,c] = ObjFunction(n,m)

% Generating the objective function data: A and B are m x n real matrices 
% with m >= n, and c is a real positive number

A = [];
B = [];

A = 2 * rand(m,n) - 1;
B = 2 * rand(m,n) - 1;

c = 10;