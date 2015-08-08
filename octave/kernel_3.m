function [parent, d] = kernel_3 (G, root)
  %% Compute the shortest path lengths and parent
  %% tree starting from vertex root on the graph
  %% represented by the sparse matrix G. Every
  %% vertex in G can be reached from root.

  N = size (G, 1);
  %% Adjust from zero labels.
  root = root + 1;
  d = inf * ones (N, 1);
  parent = zeros (N, 1);
  d (root) = 0;
  parent (root) = root;

  Q = 1:N;
  while length (Q) > 0
    [dist, idx] = min (d(Q));
    v = Q(idx);
    Q = setdiff (Q, v);
    [I, J, V] = find (G (:, v));
    for idx = 1:length(I),
      u = I(idx);
      dist_tmp = d(v) + V(idx);
      if dist_tmp < d(u),
        d(u) = dist_tmp;
        parent(u) = v;
      end
    end
  end

  %% Adjust back to zero labels.
  parent = parent - 1;

