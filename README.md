This is a depth-range-diophantine CSL parser/matcher/enumerator, an original algorithm. The overall time complexity is polynomial.

Summary:
Generate a parse tree from the expression where groups nest. Generate equation fragments bottom-up, where
  leaf => add group size constant
  repetitions => multiply by free variable
  concatention => cartesian product of pair sums
  alternation => merge equation list

The resulting equations are a sum of products, and are rewritten as linear diophantine equations. Then for some match input, solve each rewritten equation in the root level equation list set equal to the input size using a linear diophantine equation solver. The general solution is then constrained using integer programming, with the trivial starting bounds 0 <= y[i] <= L/c[i], where y[i] is a free variable in the rewrite, c[i] is its constant factor, and L is the input size. Then use a k-factors factorization of each y[i] to generate solutions in terms of the original free variables before the rewrite, which we call x[i].

These x[i] solutions are plugged back into the parse tree to perform matching recursively.

With this set up, equations can be solved and matched in parallel, and during matching, a group's children can be matched in parallel.

Time Complexity Discussion:
Linear diophantine equation solving, integer programming min/max'ing, and k-factors factorization are known polynomial time algorithms. During equation generation, each alternation => merge creates n + m equations, each concatenation => cartesian product of pair sums creates n^2 equations, and the composition of these operations is polynomial.

Implemented:
1. parse tree
2. equation gen
3. solver
4. matcher

TODO:
1. backrefs
2. enumerator
3. paper
4. parallelize
5. CUDA
6. O(1) lookup for factorization

How to:
g++ -g -std=c++23 main3.cpp
./a.out

Relevant links:
1. https://swtch.com/~rsc/regexp/regexp1.html
2. https://perl.plover.com/NPC/NPC-3SAT.html

![image](https://github.com/user-attachments/assets/428639bb-6c5f-4739-b4b1-6e5c285fb3d5)
