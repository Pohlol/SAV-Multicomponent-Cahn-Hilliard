function createMesh()
    h=[0.2 0.15 0.1 0.07 0.05 0.04 0.03 0.025 0.02];
    fd = @(p) drectangle(p,-1,1,-1,1);
    for i=8:9
    [nodes,elements] = distmesh2d( fd, @huniform, h(i), [-1,-1;1,1], [-1,-1;-1,1;1,-1;1,1] );
    nodes = 5*nodes;
    
    save(['Nodes/Square_nodes',num2str(i),'.txt'],'nodes','-ASCII');
    save(['Nodes/Square_elements',num2str(i),'.txt'],'elements','-ASCII');
    end
end

