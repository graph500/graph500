SCALE = 10;

edgefactor = 16;
maxweight = 255;
NROOT = 8;

NE = edgefactor * 2**SCALE;

ijw = edge_list (0, NE, SCALE, NE, maxweight);

tic;
G = kernel_1 (ijw);
kernel_1_time = toc;
NV = size (G, 1);

root = sample_roots (NV, NROOT, NE);

kernel_2_time = Inf * ones (NROOT, 1);
kernel_2_dmax = -ones (NROOT, 1);
kernel_2_verify_time = Inf * ones (NROOT, 1);
kernel_3_time = Inf * ones (NROOT, 1);
kernel_3_dmax = -ones (NROOT, 1);
kernel_3_verify_time = Inf * ones (NROOT, 1);

for k = 1:NROOT,
  tic;
  [parent, d] = kernel_2 (G, root(k));
  kernel_2_time(k) = toc;
  kernel_2_dmax(k) = max (d);
  tic;
  err = verify (parent, ijw, root (k), d, 1);
  kernel_2_verify_time(k) = toc;
  if err <= 0,
    error (sprintf (["BFS %d from search key %d"
		     " failed to validate: %d"],
		    k, root(k), err));
  end
end

for k = 1:NROOT,
  tic;
  [parent, d] = kernel_3 (G, root(k));
  kernel_3_time(k) = toc;
  kernel_3_dmax(k) = max (d);
  tic;
  err = verify (parent, ijw, root (k), d, 0);
  kernel_3_verify_time(k) = toc;
  if err <= 0,
    error (sprintf (["SSSP %d from search key %d"
		     " failed to validate: %d"],
		    k, root(k), err));
  end
end

output ("Jason's old laptop", SCALE, edgefactor,
	root, kernel_1_time,
	kernel_2_time, kernel_2_dmax, kernel_2_verify_time,
	kernel_3_time, kernel_3_dmax, kernel_3_verify_time,
        "An aging, soon-to-be-replaced laptop.");
