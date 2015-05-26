%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Kyle Booth
%Project: APS1022 - Financial Engineering II
%Date Due: Monday May 25, 2015 @ Noon
%Filename: dataRead.m (function file), 10 Assets in Portfolio
%Purpose: Provides raw data from YahooFinance to the main program file.

function data = dataRead(data);
data(:,:, 1) = csvread('data\BCE.csv');
data(:,:, 2) = csvread('data\BMO.csv');
data(:,:, 3) = csvread('data\CCT.csv');
data(:,:, 4) = csvread('data\IMO.csv');
data(:,:, 5) = csvread('data\MFC.csv');
data(:,:, 6) = csvread('data\POW.csv');
data(:,:, 7) = csvread('data\SAP.csv');
data(:,:, 8) = csvread('data\SNC.csv');
data(:,:, 9) = csvread('data\TDT.csv');
data(:,:, 10) = csvread('data\VRX.csv');


