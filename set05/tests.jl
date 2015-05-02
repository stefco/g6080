import Laplacian

# Test cycle
println("Testing `cycle`")
a = [1 2 3; 4 5 6; 7 8 9];
@assert Laplacian.cycle(a, 1, 2) == [2 3 1; 5 6 4; 8 9 7] "r(n) -> r(n+1) fail"
@assert Laplacian.cycle(a,-2, 2) == [2 3 1; 5 6 4; 8 9 7] "r(n) -> r(n-2) fail"
@assert Laplacian.cycle(a,-1, 2) == [3 1 2; 6 4 5; 9 7 8] "r(n) -> r(n-1) fail"
@assert Laplacian.cycle(a, 2, 2) == [3 1 2; 6 4 5; 9 7 8] "r(n) -> r(n+2) fail"
@assert Laplacian.cycle(a, 1, 1) == [4 5 6; 7 8 9; 1 2 3] "c(n) -> c(n+1) fail"
@assert Laplacian.cycle(a,-2, 1) == [4 5 6; 7 8 9; 1 2 3] "c(n) -> c(n-2) fail"
@assert Laplacian.cycle(a,-1, 1) == [7 8 9; 1 2 3; 4 5 6] "c(n) -> c(n-1) fail"
@assert Laplacian.cycle(a, 2, 1) == [7 8 9; 1 2 3; 4 5 6] "c(n) -> c(n+2) fail"
println("`cycle` passed tests.")

# Test delta
println("Testing `δ`")
@assert Laplacian.δ(i->i,j->j,3) * a == a "identity matrix operational fail"
@assert Laplacian.δ(i->i,j->j,10) == eye(10,10) "identity matrix equality fail"
@assert Laplacian.δ(i->i,j->mod1(j+1,5),5) * [1:5] == [5,1:4] "permute fail"
println("`δ` passed tests.")

# Test opsplit
println("Testing `opsplit`")
@assert Laplacian.opsplit(a, [2]) == ([1 0 3; 4 0 6; 7 0 9], [0 2 0; 0 5 0; 0 8 0]) "split fail"
println("`opsplit` passed tests")
