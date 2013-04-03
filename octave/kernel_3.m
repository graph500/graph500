function [parent, d] = kernel_3 (G, root)
  ## Compute the shortest path lengths and parent
  ## tree starting from vertex root on the graph
  ## represented by the sparse matrix G. Every
  ## vertex in G can be reached from root.

  N = size (G, 1);
  ## Adjust from zero labels.
  root = root + 1;
  d = inf * ones (N, 1);
  parent = zeros (N, 1);
  d (root) = 0;
  parent (root) = root;

  ## Very inefficient version of Dijkstra's algorithm.
  Q = 1:N;
  while length (Q) > 0,
    [du, qk] = min (d(Q));
    u = Q(qk);
    Q = setdiff (Q, u);
    [V, J, W] = find (G (:, u));
    for vk = 1:length (V),
      v = V(vk);
      dtmp = d(u) + W(vk);
      if dtmp < d(v),
        d(v) = dtmp;
        parent(v) = u;
      end
    end
  end

  ## Adjust back to zero labels.
  parent -= 1;
endfunction
