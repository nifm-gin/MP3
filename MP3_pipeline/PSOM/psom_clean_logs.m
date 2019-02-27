function [] = psom_clean_logs(path_logs)
% Clean up all temporary files in a log folder
%
% Syntax: PSOM_CLEAN(PATH_LOGS)
% PATH_LOGS (string) the logs folder 
%
% Copyright (c) Pierre Bellec, Centre de recherche de l'institut de
% griatrie de Montral, Dpartement d'informatique et de recherche
% oprationnelle, Universit de Montral, Canada, 2010-2012
% Maintainer : pierre.bellec@criugm.qc.ca
% See licensing information in the code.

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Seting up default arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
    error('Please specify PATH_LOGS')
end

list_files = { 'PIPE.end' , 'PIPE.running' , 'PIPE.failed' , 'PIPE.finished' , 'PIPE.exit' , 'PIPE.kill' , ...
 'heartbeat.mat' , 'PIPE.oqsub' , 'PIPE.eqsub' , 'news_feed.csv' };

for num_f = 1:length(list_files)
    file_name = [path_logs list_files{num_f}];
    if psom_exist(file_name)
        psom_clean(file_name,struct('flag_verbose',false));
    end
end
