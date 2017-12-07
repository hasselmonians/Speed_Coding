function [z] = get_zScore_speed(self,varargin)

z = (self.svel - mean(self.svel))/std(self.svel);

