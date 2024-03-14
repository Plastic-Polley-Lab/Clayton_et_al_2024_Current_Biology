clear;close all;clc
load('C:\Users\kx776\Dropbox\File_Transfer\Dan - Ke\Figure2\RawData\summaryData_all.mat')
% repeats.ko = {[6, 9], [7, 10], [8,11]};
% repeats.wt = {[8, 11], [9, 12], [10, 13]};
cutOff = 1;
summaryData = spont_movement(summaryData, cutOff);

additional_fig = 0;
plots = summary_plot_orofacial(summaryData, additional_fig, cutOff, repeats);