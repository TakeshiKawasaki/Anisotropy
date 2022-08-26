import numpy as np
import pandas as pd

def custom_df(args, df):
    # this function adds additional information into the dataframe
    # the dataframe is later used in writing .vtk file   

    # add column 'scale' and 'ratio' 
    try:
        _ = df['scale'] + df['ratio']
    except:
        print(f"There is no column 'scale' and 'ratio' in the input data, adding 'scale = {args.scale}' and 'ratio = {args.ratio}' from the input arguments.") 
        df['scale'] = args.scale
        df['ratio'] = args.ratio

    # add column 'theta' from nx and ny
    df['theta'] = np.arctan(df['ny']/df['nx']) 
    # add column 'phi' from nz    
    df['phi'] = np.arccos(np.absolute(df['nz'])) 

    # convert from radian to degrees for more intuitive visualisation
    d = 180/np.pi
    df['theta'] = df['theta']*d
    df['phi'] = df['phi']*d

    # define a plane to extract cross section 
    try:
        a = args.plane_grad[0]
        b = args.plane_grad[1]
        c = args.plane_grad[2]
        d = args.plane_grad[3]
        w = args.plane_width
        x0 = args.plane_origin[0]
        y0 = args.plane_origin[1]
        z0 = args.plane_origin[2]
        print(f"plane equation is specified and is ({a})x + ({b})y + ({c})z + ({d}) = 0")
        print(f"plane goes through a point ({x0}, {y0}, {z0})")
        cross_section_flag = True  
    except:
        cross_section_flag = False
        w = None
     
    # if width of the planes are not provided via ars, it uses the value form args.scale
    if w is None:
        w = df.at[0,'scale']
    

    # extract particles to visualise
    # do this by adding column 'visibility' 
    for i, row in df.iterrows():
        
        if cross_section_flag ==True: # run this block if plane is defined
            # ensure the particles' origin satisfie the condition that they are inside 2 planes
            upper_bound = (a*(row['rx']-x0)+ b*(row['ry']-y0) + c*(row['rz']-z0) + d - w/2 < 0)
            lower_bound = (a*(row['rx']-x0)+ b*(row['ry']-y0) + c*(row['rz']-z0) + d + w/2 > 0)
            condition = upper_bound*lower_bound

            if condition == 1:
                df.at[i,'visibility'] = 'y'
            else:
                df.at[i,'visibility'] = 'n'
        else:
            df.at[i,'visibility'] = 'y'



