function c = connectivity(x,y,z,bw)
%CONNECTIVITY3 calculates the number of voxels in the binary volume 'bw'
%adjacent to the voxel of coordinates (x,y,z). Diagonal adjacency is not
%considered. If the voxel in (x,y,z) is of value 0, the connectivity
%returned is NaN
c = NaN;
if bw(x,y,z) == 1
    neighbx = bw([x+1 x-1],y,z);
    neighby = bw(x,[y+1 y-1],z);
    neighbz = bw(x,y,[z+1 z-1]);
    c = nnz(neighbx)+nnz(neighby)+nnz(neighbz);
end
end