function readarray(fileName)
    f=open(fileName);
    arr=Float64[];
    a = readline(f);    # Get the first line
    while a != ""
        push!(arr,parsefloat(Float64,a));
        a=readline(f);  # Get subsequent lines
    end
    return arr;
end
