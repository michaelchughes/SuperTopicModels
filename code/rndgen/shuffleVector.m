function xs = shuffleVector( xs )

if length(xs)>1
    xs = xs( randperm(length(xs)) );
end
    