%this function returns a nx3 array of n unique colors
function colors = generate_unique_colors(n)
    colors = rand(n,3);
    colors = unique(colors,'rows','first');
    if size(colors,1) < n
        uniqueness = false;
    else
        uniqueness = true;
    end
    while(~uniqueness)
        nbreUniqueColors = size(colors,1);
        for i=(nbreUniqueColors+1):n
            colors(i,:) = rand(1,3);
        end
        if size(colors,1) < n
            uniqueness = false;
        else
            uniqueness = true;
        end
    end
    colors = colors(randperm(size(colors, 1)), :);%shuffle rows randomly
end