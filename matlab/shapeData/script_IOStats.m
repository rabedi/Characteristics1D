dd2s = [0.9, 0.5, 0.1];
shapes = [1, 1.5, 2, 3, 4];
ndd2 = length(dd2s);
nshape = length(shapes);
for dd2i = 1:ndd2
    dd2 = dd2s(dd2i);
    for shapei = 1:nshape
        shape = shapes(shapei);
        fprintf(1, 'start dd2 = %g, ', dd2);
        fprintf(1, ' shape = %g\n', shape);
        IO_shape_dd2(shape, dd2);
    end
end