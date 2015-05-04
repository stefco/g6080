import Laplacian
import ConjugateGradient

# Test delta
a = [1 2 3; 4 5 6; 7 8 9];
println("Testing `δ`")
@assert Laplacian.δ(i->i,j->j,3) * a == a "identity matrix operational fail"
@assert Laplacian.δ(i->i,j->j,10) == eye(10,10) "identity matrix equality fail"
@assert Laplacian.δ(i->i,j->mod1(j+1,5),5) * [1:5] == [5,1:4] "permute fail"
println("`δ` passed tests.")

# Test opsplit
println("Testing `opsplit`")
@assert Laplacian.opsplit(a, [0,1,0]) == ([1 0 3; 4 0 6; 7 0 9], [0 2 0; 0 5 0; 0 8 0]) "split fail"
println("`opsplit` passed tests")

# Test conjgrad
println("Testing `conjgrad`")
c = [1 2 3; 2 4 5; 3 5 6]
for i in 1:3
    x = rand(3)
    b = c*x
    xx = ConjugateGradient.conjgrad(c,b)
    @assert dot(xx-x, xx-x) < 1e-10 "inaccurate conjugate gradient"
end
println("`conjgrad` passed tests")
