import os



from utils.tensor import tensor

def df2vtk(args, df, file):
    # input
    # dataframe containing tensor components, file
    # output 
    # vtk file which can be used in paraview


    # extract the filename of input.dat
    base=os.path.basename(file)
    filename = os.path.splitext(base)[0]
    filepath = os.path.join(args.vtk_dir, filename + '.vtk') # output file has the same name as the inout file but is .vtk
    
    
    # declare buffer and define file settings
    buffer = ['# vtk DataFile Version 5.1\n', 'vtk output\n', 'ASCII\n', 'DATASET POLYDATA\n']
    buffer.append('POINTS %u float\n' % (df.shape[0]))

    # extract particle locations
    for _, row in df.iterrows():
        buffer.append('%6f %6f %6f\n' % (row['rx'], row['ry'], row['rz']))

    # calculate tensor components to control the shape
    buffer.append('POINT_DATA %u\n' % (df.shape[0]))
    buffer.append('TENSORS tensors float\n')
    for _, row in df.iterrows():
        if row["visibility"] == "n":
            buffer.append('%6f %6f %6f %6f %6f %6f %6f %6f %6f\n' % (0,0,0,0,0,0,0,0,0))
        else:
            XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ = tensor(row['nx'], row['ny'], row['nz'], row['scale'], row['ratio'])
            buffer.append('%6f %6f %6f %6f %6f %6f %6f %6f %6f\n' % (XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ))

    # add scalar information - this helps coloring and partially visualising particles 
    # theta
    buffer.append('SCALARS theta float 1\n')
    buffer.append('LOOKUP_TABLE default\n')
    for _, row in df.iterrows():
        buffer.append('%f\n' % (row['theta']))
    # phi
    buffer.append('SCALARS phi float 1\n')
    buffer.append('LOOKUP_TABLE default\n')
    for _, row in df.iterrows():
        buffer.append('%f\n' % (row['phi']))
    # scale (corresponds to radius)
    buffer.append('SCALARS scale float 1\n')
    buffer.append('LOOKUP_TABLE default\n')
    for _, row in df.iterrows():
        buffer.append('%f\n' % (row['scale']))
    # ratio
    buffer.append('SCALARS ratio float 1\n')
    buffer.append('LOOKUP_TABLE default\n')
    for _, row in df.iterrows():
        buffer.append('%f\n' % (row['ratio']))


    # join strings in buffer
    vtk_string = "".join(buffer)

    # write it into a file
    with open(filepath, 'w') as f:
        f.write(vtk_string)
