clear
clc

load ./Networks/SceMINT.mat;

A = SceMINTnewNet;
LocalNodes = 20;

[Clusters] = ExactLocal(A, LocalNodes);