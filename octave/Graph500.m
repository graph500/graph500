SCALE = 10;
edgefactor = 16;
NBFS = 64;

rand ("seed", 103);

ij = kronecker_generator (SCALE, edgefactor);

tic;
G = kernel_1 (ij);
kernel_1_time = toc;

N = size (G, 1);
coldeg = full (spstats (G));
search_key = randperm (N);
search_key(coldeg(search_key) == 0) = [];
if length (search_key) > NBFS,
  search_key = search_key(1:NBFS);
else
  NBFS = length (search_key);
end
search_key = search_key - 1;  

kernel_2_time = Inf * ones (NBFS, 1);
kernel_2_nedge = zeros (NBFS, 1);

for k = 1:NBFS,
  tic;
  parent = kernel_2 (G, search_key(k));
  kernel_2_time(k) = toc;
  err = validate (parent, ij, search_key (k));
  if err <= 0,
    error (sprintf ("BFS %d from search key %d failed to validate: %d",
		    k, search_key(k), err));
  end
  slice = find (parent >= 0);
  kernel_2_nedge(k) = nnz (G(slice, slice))/2;
end

output (SCALE, edgefactor, NBFS, kernel_1_time, kernel_2_time, kernel_2_nedge);
