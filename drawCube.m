function drawCube(cube, color)
    % Extract cube coordinates and size
    x = cube(1);
    y = cube(2);
    z = cube(3);
    size = cube(4);
    
    % Define cube vertices
    vertices = [x, y, z;
                x+size, y, z;
                x+size, y+size, z;
                x, y+size, z;
                x, y, z+size;
                x+size, y, z+size;
                x+size, y+size, z+size;
                x, y+size, z+size];
    
    % Define cube faces
    faces = [1, 2, 3, 4;
             5, 6, 7, 8;
             1, 2, 6, 5;
             2, 3, 7, 6;
             3, 4, 8, 7;
             4, 1, 5, 8];
    
    % Plot cube
    patch('Faces', faces, 'Vertices', vertices, 'FaceColor', color, 'EdgeColor', 'k');
end