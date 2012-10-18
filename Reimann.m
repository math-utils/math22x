function y = Reimann(x,ul,ur)
if ul < 0 || ur < 0
    error('This program is only defined for positive "ul" and "ur". You need to rewrite this program to allow for negative values');
end
y = ul;
if x ~= 0
    if ul > ur
        if x >= (ul + ur)/2
            y = ur;
        end
    elseif ul < ur
        if x > ur
            y = ur;
        elseif x > ul && x < ur
            y = x;
        end
    end
end
end