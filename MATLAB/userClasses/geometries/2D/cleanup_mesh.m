% Function to clean up unused nodes and elements
function [clean_coordinates, clean_elements] = cleanup_mesh(coords, elems)
        used_nodes = unique(elems(:));
        clean_coordinates = coords(used_nodes, :);

        % Map old indices to new indices
        node_map = containers.Map(used_nodes, 1:length(used_nodes));
        clean_elements = arrayfun(@(x) node_map(x), elems);
end
