function g = floodgraph(g, basenode, subn_idx)
%FLOODGRAPH(g,node,subn_idx) assigns to all nodes linked to 'basenode' by any
%   path the label identified in 'subn_idx', identifying the subnetowrk to
%   which 'basenode' belongs. The graph returned has an additional field 
%   'subN' that contains the label of the subnetwork to which each node
%   belongs.

    g.Nodes(basenode,:).subN = subn_idx;
    neigh = neighbors(g,basenode);
    if all(g.Nodes(neigh,:).subN~=0)
        return;
    end
    for i=1:numel(neigh)
        if g.Nodes(neigh,:).subN(i) == 0
            g = floodgraph(g,neigh(i),subn_idx);
        end
    end
end