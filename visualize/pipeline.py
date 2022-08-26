# libraries 
import argparse
import pandas as pd
import os
from string import digits
import json

from utils.dat2df import dat2df
from utils.custom_df import custom_df
from utils.df2vtk import df2vtk



# input arguments
# Create the parser
def parse_agrs():
    parser = argparse.ArgumentParser(description='program to which takes the result from MD, then visualise it using Paraview')
    
    parser.add_argument('--input_dir', help='Provide a path to a directory containing input files', required=True)
    parser.add_argument('--vtk_dir', default='.', help='Provide a directory path to save processed vtk file')
    parser.add_argument('--save_csv', default='off', help='if "on", the program processed csv file')
    parser.add_argument('--video', default='off', help='if "on", the program creates .vtk.series file which can be animated in Paraview')
    parser.add_argument('--scale', default='1' , help='Scale determines the radius of particles. Every particles will have the same scale', type=float)
    parser.add_argument('--ratio', default='2' , help='Ratio determines the shape of particles. The parameter is the ratio of the long and short axis(long/short), which are the length of the eigenvectors. Every particles will have the same ratio', type=float)
    parser.add_argument("--plane_grad", nargs=4, help='parameter to define the gradient of the plane' , type=float)
    parser.add_argument("--plane_width", help='width of 2 plane cutting the 3D particle models', type=float)
    parser.add_argument("--plane_origin", nargs=3, default=[0, 0, 0], help='origin of the plane', type=float)
    args = parser.parse_args()
    return args


# single data to vtk
def vtk_generator(args, file):
    ## load .dat file, convert .dat to dataframe
    filename, file_extension = os.path.splitext(file)
    if file_extension == ".csv":
        df = pd.read_csv(file)
    else:
        df = dat2df(file)

    ## add more information into dataframe. e.g. angles, color etc
    custom_df(args, df)
    if args.save_csv == 'on':
        # save csv file
        csv = os.path.join(args.vtk_dir, os.path.basename(filename))
        df.to_csv(csv, index=False)
        print(f'saved csv {os.path.basename(csv)}')


    ## create vtk from dataframe
    df2vtk(args, df, file)

    print(f'processed {filename}\n')


# function to create .vtk.series from multiple 
def vtkseries_gen(vtk_dir):
    # takes vtk files in a list
    files = [os.path.join(vtk_dir, file) for file in os.listdir(vtk_dir) if file.endswith('.vtk')]
    # sort the list
    files.sort()
    print('combining following files to create .vtk.series file')
    print(files)

    # take filename from one of the files and remove number
    filename = os.path.basename(files[0]) # takes base filename
    filename = os.path.splitext(filename)[0] # remove extension
    filename = ''.join((x for x in filename if not x.isdigit())) # remove digits
    filepath = os.path.join(vtk_dir, filename + '.vtk.series') 

    # iterate through vtk files in the directory
    data=[]
    for time, file in enumerate(files):
        item = {"name": os.path.basename(file), "time": time}
        data.append(item)

    # save into .vtk.series as a json format
    string = {
        "file-series-version" : "1.0",
        "files" : data
    }

    jsonString = json.dumps(string, indent = 2)
    with open(filepath, "w") as jsonFile:
        jsonFile.write(jsonString)


def pipeline():
    # takes input data directory path
    args = parse_agrs()
    files = [os.path.join(args.input_dir, file) for file in os.listdir(args.input_dir)]

    # run pipeline for each file
    for file in files:
        vtk_generator(args, file)

    if args.video == 'on':
        vtkseries_gen(args.vtk_dir)
        print('saved .vtk.series file')
    
    print('\nCompleted pipeline')


if __name__ == '__main__':
    pipeline()

