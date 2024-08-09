##################
### Code taken adapted and modified from first model prototype at:
### REDInet/Notebooks&Scripts/CNN_wavenet_kindney_first_model_prototype/CNN_inference.py
##################

# import basic modules
import pysam, gzip,os
import numpy as np
import pandas as pd
from datetime import datetime
from tensorflow.python.platform.tf_logging import warn
import tensorflow as tf
import keras
from sklearn.preprocessing import OneHotEncoder
from tqdm import tqdm
import argparse

### define utils function ###
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

### define main function ###
def CNN_inference(reditable_filepath, model_filepath, output_table_prefix_filepath=None, cov_threshold=50, AGfreq_threshold=0.01, minAGsubs=3):

    if output_table_prefix_filepath == None:
        output_table_prefix_filepath = reditable_filepath

    interval=101 # it depends on the model

    #### START ####
    
    starttime = datetime.now()

    editing = []
    # open reditools table (compressed and indexed via tabix)
    with gzip.open(reditable_filepath) as redi:
        for c,l in enumerate(redi):
            line = l.decode("utf-8").rstrip().split("\t")
            if line[2] == "A":
                if line[4] != "-":
                    if int(line[4]) >= cov_threshold:
                        if "AG" in line[7]:
                            AG_rna = eval(line[6])[2]/sum(eval(line[6]))
                            if AG_rna >= AGfreq_threshold:
                                if eval(line[6])[2] >= minAGsubs:
                                    editing.append(line)
            if c % 50000000 == 0:
                print(f"\t[{datetime.now()}] Sites evaluated: {c}", flush=True)
    print(f"[{datetime.now()}] Total evaluated rows:", c, flush=True)
    editing = pd.DataFrame(editing)
    print(f"[{datetime.now()}] Total extracted Candidates Editing sites for current sample:", editing.shape[0], flush=True)
    stoptime = datetime.now()
    print(f"[{datetime.now()}] Extraction of Editing Candidates finished for current sample. Elapsed time: {stoptime-starttime}.", flush=True)


    columns = ["Region", "Position", "Ref", "Strand", "Cov", "Qual", "[A,C,G,T]", "AllSubs", "Freq", "gCov", "gQual", "g[A,C,G,T]", "gAllSubs", "gFreq"]
    editing.columns = columns
    print("", flush=True)
    print(editing, flush=True)
    print(f"[{datetime.now()}] Starting extraction of intervals.", flush=True)

    # create ohe instance
    ohe = OneHotEncoder()
    ohe.fit(np.array(["A", "C", "G", "T"]).reshape(-1, 1))

    intervals = []
    starttime_preds = datetime.now()
    total_extracted = 0

    features_extracted_filepath = output_table_prefix_filepath + ".feature_vectors.tsv"
    features_extracted = open(features_extracted_filepath, "w")
    # mantain only non chrM regions
    df = editing.query("Region != 'chrM'")
    print(f"[{datetime.now()}] Loading reditable with tabix and pysam:", reditable_filepath, flush=True)
    start_time = datetime.now()
    srr = pysam.TabixFile(reditable_filepath)
    # extract pos examples
    with tqdm(total=df.shape[0], position=0, leave=True) as pbar:
        for site in df.itertuples():
            start = int(site.Position) - ((interval-1)/2)
            stop = int(site.Position) + ((interval-1)/2)
            AGrna = eval(site._7)[2]/sum(eval(site._7))
            srr_interval = []
            for s in srr.fetch(site.Region, start-1, stop):
                srr_interval.append(s.split("\t"))
            srr_interval = pd.DataFrame(srr_interval, columns=columns)
            # assess wheter interval is of the required length and if the entire interval is on the same strand
            if srr_interval.shape[0] == interval and len(set(srr_interval["Strand"])) == 1:
                strand = site.Strand
                #print("DEBUG - Strand type:", type(strand), "Strand:", strand)
                if strand in ["0","1"]: # avoid unstranded sites (0 --> minus strand; 1 --> plus strand)
                    intervals.append([site.Region, site.Position, site.Strand, AGrna, start, stop, stop-start + 1, srr_interval.shape[0]])
                    total_extracted += 1
                    # encode features vector and write to disk
                    seq = srr_interval.Ref.values.reshape(-1,1)
                    seq_ohe = ohe.transform(seq).toarray().T
                    vects_freqs = []
                    vects = []
                    for vect in srr_interval["[A,C,G,T]"]:
                        vect = np.array(eval(vect))
                        cov = sum(vect)
                        vect_freqs = vect / cov
                        vects_freqs.append(vect_freqs)
                        vects.append(vect)
                    vects_freqs = np.array(vects_freqs).T
                    vects = np.array(vects).T
                    site = pd.concat([pd.DataFrame(seq_ohe), pd.DataFrame(vects_freqs)])
                    if strand == 0: # flip horizontally if it is on reverse strand
                        site = pd.DataFrame(np.flip(site.values, axis=1))
                    # save to disk (append mode)
                    site.to_csv(features_extracted, mode="a", sep="\t", header = None, index=None)
            pbar.update(1)
    # create dataframe for editing sites candidates with complete intervaal
    intervals = pd.DataFrame(intervals)
    print(f"[{datetime.now()}] Total extracted Editing sites: {total_extracted}.", flush=True)
    stop_time_global = datetime.now()
    print(f"[{datetime.now()}] Features Extraction Finished. Elapsed time {datetime.now()-starttime_preds}.", flush=True)
    features_extracted.close()

    # start prediction step on the extracted features
    # loading features and preprocess these for inference
    ###--- START LOADING OF DATA ---###
    print(f"[{datetime.now()}] Loading features extracted of editing candidates from: {features_extracted_filepath}", flush=True)
    X = pd.read_table(features_extracted_filepath, header=None)
    # trasform into 3d tensor from series of 2d table from tsv input file
    X_3d = from_2_to_3_dimensions(X, 8)
    # create a log2 version on pseudo count to zoom differences among edited and non edited sites
    # create a tuple containing min and max log2 expected values useful to normalize frequency features
    log_range = (-13.28771238, 0.000144262)
    # go to pseudofrequency
    X_3d_log2 = X_3d.copy()
    X_3d_log2[:,:,4:] = np.log2(X_3d_log2[:,:,4:]+0.0001)
    # normalize log2 preudo-frequencies using log2 boudaries
    X_3d_log2[:,:,4:] = (X_3d_log2[:,:,4:]-log_range[0]) / (log_range[1] - log_range[0])
    # make prediction
    # load the model
    print(f"[{datetime.now()}] Loading model from:", model_filepath, flush=True)
    model = tf.keras.models.load_model(model_filepath)
    y_hat_proba = model.predict(X_3d_log2, batch_size=100)
    #print("DEBUG - y_hat_proba:\n", y_hat_proba)
    # adapt y_hat_proba for models with single node output (signmoid activation function)
    if y_hat_proba.shape[1] == 1:
        y_hat_proba_snp = 1.0 - y_hat_proba
        y_hat_proba_ed = y_hat_proba
        y_hat_proba = np.hstack([y_hat_proba_snp, y_hat_proba_ed])
    # produce output vector (0 --> SNP; 1 --> Editing)
    y_hat = np.array([np.argmax(i) for i in y_hat_proba])
    # append to intervals
    intervals["snp_proba"] = y_hat_proba[:,0]
    intervals["ed_proba"] = y_hat_proba[:,1]
    intervals["y_hat"] = y_hat
    intervals.columns = ["region", "position", "Strand", "FreqAGrna", "start", "stop", "int_len", "TabixLen", "snp_proba", "ed_proba", "y_hat"]
    # save to disk
    intervals.to_csv(output_table_prefix_filepath+".predictions.tsv", sep="\t", index=None)
    print(f"[{datetime.now()}] Predictions concluded. File saved to: {output_table_prefix_filepath+'.predictions.tsv'}", flush=True)
    print(f"[{datetime.now()}] Computation Finished. Total Elapsed time: {datetime.now()-starttime}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=f"""REDInet_Inference_light_ver.py""")
    parser.add_argument("-r",
                        "--reditable",
                        required=True,
                        type=str,
                        help="--reditable: \t a <str> with the fullpath for the input reditable file.")
    parser.add_argument("-m",
                        "--model",
                        required=False,
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),"REDInet.h5"),
                        type=str,
                        help="--model: \t a <str> with the fullpath for the model file. [REDIportal model]")
    parser.add_argument("-o",
                        "--output_table_prefix",
                        required=False,
                        default=None,
                        type=str,
                        help="--output_table_prefix: \t a <str> Indicating the prefix for the output files. [None]")
    parser.add_argument("-c",
                        "--cov_threshold",
                        required=False,
                        default=50,
                        type=int,
                        help="--cov_threshold: \t a <int> Indicating the minimum coverage to make inference. [50]")
    parser.add_argument("-f",
                        "--AGfreq_threshold",
                        required=False,
                        default=0.01,
                        type=float,
                        help="--AGfreq_threshold: \t a <float> Indicating the minimum AG substitution frequency to make inference. [0.01]")
    parser.add_argument("-s",
                        "--minAGsubs",
                        required=False,
                        default=3,
                        type=int,
                        help="--minAGsubs: \t a <int> Indicating the minimum AG substitutions to make inference. [3]")

    args = parser.parse_args()
    reditable = args.reditable
    model = args.model
    output_table_prefix = args.output_table_prefix
    cov_threshold = args.cov_threshold
    AGfreq_threshold = args.AGfreq_threshold
    minAGsubs = args.minAGsubs

    # CNN_inference(reditable_filepath, model_filepath, output_table_prefix_filepath=None, cov_threshold=50, AGfreq_threshold=0.01)
    # print some starting info related to version, used program and to the input arguments
    print(f"[{datetime.now()}] REDInet_Inference_light_ver.py:", flush=True)
    print(f"[{datetime.now()}] Input arguments:", flush=True)
    for argument in args.__dict__.keys():
        print(f"\t- {argument} --> {args.__dict__[argument]}", flush=True)

    # launch main function
    CNN_inference(reditable_filepath=reditable, 
                  model_filepath=model, 
                  output_table_prefix_filepath=output_table_prefix, 
                  cov_threshold=cov_threshold, 
                  AGfreq_threshold=AGfreq_threshold,
                  minAGsubs=minAGsubs)
