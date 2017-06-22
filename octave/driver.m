SCALE = 10;
edgefactor = 16;
NBFS = 64;

rand ("seed", 103);

ijw = kronecker_generator (SCALE, edgefactor);

tic;
G = kernel_1 (ijw);
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
kernel_3_time = Inf * ones (NBFS, 1);
kernel_3_nedge = zeros (NBFS, 1);

indeg = histc (ijw(:), 1:N); % For computing the number of edges

for k = 1:NBFS,
  tic;
  parent = kernel_2 (G, search_key(k));
  kernel_2_time(k) = toc;
  err = validate (parent, ijw, search_key(k), 0, false);
  if err <= 0,
    error (sprintf ("BFS %d from search key %d failed to validate: %d",
		    k, search_key(k), err));
  end
  kernel_2_nedge(k) = sum (indeg(parent >= 0))/2; % Volume/2

  tic;
  [parent, d] = kernel_3 (G, search_key(k));
  kernel_3_time(k) = toc;
  err = validate (parent, ijw, search_key (k), d, true);
  if err <= 0,
    error (sprintf ("SSSP %d from search key %d failed to validate: %d",
		    k, search_key(k), err));
  end
  kernel_3_nedge(k) = sum (indeg(parent >= 0))/2; % Volume/2
end

output (SCALE, NBFS, NBFS, kernel_1_time, kernel_2_time, kernel_2_nedge, kernel_3_time, kernel_3_nedge);
