function ij = kronecker_generator (SCALE, edgefactor)
%% Generate an edgelist according to the Graph500
%% parameters.  In this sample, the edge list is
%% returned in an array with two rows, where StartVertex
%% is first row and EndVertex is the second.  The vertex
%% labels start at zero.
%%
%% Example, creating a sparse matrix for viewing:
%%   ij = kronecker_generator (10, 16);
%%   G = sparse (ij(1,:)+1, ij(2,:)+1, ones (1, size (ij, 2)));
%%   spy (G);
%% The spy plot should appear fairly dense. Any locality
%% is removed by the final permutations.

  %% Set number of vertices.
  N = 2^SCALE;

  %% Set number of edges.
  M = edgefactor * N;

  %% Set initiator probabilities.
  [A, B, C] = deal (0.57, 0.19, 0.19);

  %% Create index arrays.
  ij = ones (2, M);
  %% Loop over each order of bit.
  ab = A + B;
  c_norm = C/(1 - (A + B));
  a_norm = A/(A + B);

  for ib = 1:SCALE,
    %% Compare with probabilities and set bits of indices.
    ii_bit = rand (1, M) > ab;
    jj_bit = rand (1, M) > ( c_norm * ii_bit + a_norm * not (ii_bit) );
    ij = ij + 2^(ib-1) * [ii_bit; jj_bit];
  end

  %% Permute vertex labels
  p = randperm (N);
  ij = p(ij);

  %% Permute the edge list
  p = randperm (M);
  ij = ij(:, p);

  %% Adjust to zero-based labels.
  ij = ij - 1;
