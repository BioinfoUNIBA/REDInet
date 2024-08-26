##################
###
##################
import numpy as np

def from_2_to_3_dimensions(pandas_df, n_features):
    '''
    A function to convert the 2 dimensional table with events (16 rows-features per sample ) into 
    a 3-dimensional array to feed the neural network with a shape of:
        - n rows => n samples
        - m columns => m positions-timepoint
        - z third positions => z features per position-timepoint (16 features per position)
    '''
    # convert the loaded csv from pandas dataframe to numpy 2D array
    np_csv = pandas_df.values
    # convert into 3D array
    array_list = []
    for i in range(0, np_csv.shape[0], n_features):
        array_list.append(np_csv[i : i+n_features, :].T)
    array_3d = np.stack(array_list)
    return array_3d
