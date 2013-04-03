function [parent, d] = kernel_2 (G, root)
  ## Compute a breadth-first search tree starting
  ## from vertex root on the graph represented by
  ## the sparse matrix G

  N = size (G, 1);
  ## Adjust from zero labels.
  root = root + 1;
  parent = zeros (N, 1);
  parent (root) = root;
  d = zeros (N, 1);

  vlist = zeros (N, 1);
  vlist(1) = root;
  lastk = 1;
  for k = 1:N,
    v = vlist(k);
    if v == 0, break; end
    [I,J,V] = find (G(:, v));
    nxt = I(parent(I) == 0);
    parent(nxt) = v;
    d(nxt) = d(v) + 1;
    vlist(lastk + (1:length (nxt))) = nxt;
    lastk = lastk + length (nxt);
  end

  ## Adjust to zero labels.
  parent = parent - 1;
endfunction
