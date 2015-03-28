## PART 3

# Declare autocorrelation function
include("autocorrelation.jl");

# Create x and y arrays for plotting
xcorr = [50:50:5000];
ycorr = zeros( length( xcorr ), 5 );
for a in [1:5]
    ycorr(:,a) = autocorrelation(M, v, a, xcorr);
end
