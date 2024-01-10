function [data, tp_names, gene_names, struct_names]  = load_data(~)
%
%
%

  persistent local_data %#ok
  persistent local_tp_names %#ok
  persistent local_struct_names %#ok
  

  if ~isempty(local_data)
    data = local_data;
    tp_names = local_tp_names;
    struct_names = local_struct_names;
    return
  end
  
  % Load data normalized by mean over all genes
  filename = 'dev_data_raw.mat';
  dirname = './';
  fullname  = fullfile(dirname, filename);
  addpath(dirname);
  fprintf('Load data from file = "%s"\n', fullname);
  
  x = load(fullname, 'dev_data_raw', 'gene_names');
  data = x.dev_data_raw;
  gene_names = x.gene_names;  
  [num_tps, ~, num_structs] = size(data);  

  tp_names = timePoint_ind2name(1:num_tps);
  struct_names = struct_ind2name(1:num_structs);
  
