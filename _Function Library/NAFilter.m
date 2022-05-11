function F = NAFilter(F, x, f, NA)

    [x_mesh, y_mesh] = meshgrid(x);
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
    F(r_mesh > f*NA) = 0;    

end