% Decimal truncation function, ripped from RDM's cosine_recursion_trunc.m

function [r] = trunc(y)


    if ( eP < 1 + eps )

        r = y;

    else

        %  Preserve sign of y

        s = sign(y);

        %  First figure out how to make the value a number of
        %  order 1.

        fly = floor(log10(abs(y)));

        %  Make y a positive number of order 1

        y = abs(y * 10^(-fly));

        %  Now truncate y to P digits of precision
        %  floor() truncates towards zero

        y =  (floor(y * eP))/eP;

        %  Now multiply by 10 to the correct power

        r = s * y * 10^fly;
    end

end
