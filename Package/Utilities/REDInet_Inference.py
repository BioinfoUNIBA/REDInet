import pandas as pd
import numpy as np
import tensorflow as tf
import os, sys, pysam, gzip, subprocess, shlex, time, argparse
from datetime import datetime
from tqdm import tqdm
from glob import glob
from multiprocessing import Pool
from sklearn.preprocessing import OneHotEncoder
from  tensorflow import keras

class inference():
    
    def __init__(self, cov_threshold, AGfreq_threshold, AG_min, Multiprocessing, 
                 Processes, Utilities_path, Data_path, Results_path, Files_names):
        self.cov_threshold = cov_threshold
        self.AGfreq_threshold = AGfreq_threshold
        self.AG_min =  AG_min
        self.Multiprocessing = Multiprocessing
        self.Processes = Processes
        self.Utilities_path = Utilities_path
        self.Data_path = Data_path
        self.Results_path = Results_path
        self.Files_names = Files_names
        
    def from_2_to_3_dimensions(self, pandas_df):

        time_step = pandas_df.shape[1]
        n_matrices = int(pandas_df.shape[0]/8)

        return np.transpose(np.reshape(pandas_df.to_numpy(), (n_matrices, 8, time_step)), (0, 2, 1))

    def log_preprocessing(self, np_array):

        log_range = (-13.28771238, 0.000144262)
        np_array[:, :, 4:] = np.log2(np_array[:, :, 4:]+0.0001)
        np_array[:, :, 4:] = (np_array[:, :, 4:] - log_range[0]) / (log_range[1] - log_range[0])

        return np_array
      
    def extraction(self, name):    
        
        intervals = []
        features_matrices = pd.DataFrame()
        
        starttime = datetime.now()
        data = []
        data_discarded = []
        prefix = os.path.join(self.Data_path, name)
        print(f"Candidates sites identification in {name} start\n", flush=True)
        with gzip.open(prefix+".gz") as redi:
            for c,l in enumerate(redi):
                line = l.decode("utf-8").rstrip().split("\t")
                if line[0].find("chr") != -1:
                    if line[4] != "-":
                        if line[2] == "A":
                            if line[7] == "AG": 
                                if int(line[4]) >= self.cov_threshold:
                                    AG_rna = eval(line[6])[2]/sum(eval(line[6]))
                                    if AG_rna >= self.AGfreq_threshold:
                                        if eval(line[6])[2] >= self.AG_min:
                                            data.append(line)
                                        else:
                                            data_discarded.append(line)
                                    else:
                                        data_discarded.append(line)
                                else:
                                    data_discarded.append(line)
                                    
                if c % 50000000 == 0:
                    print(f"\tSites evaluated in {name}: {c}", flush=True)
        print(f"Total evaluated rows in {name}: ", c, flush=True)
        data = pd.DataFrame(data)
        data = data.iloc[:, 0:9]
        data_discarded =  pd.DataFrame(data_discarded)
        data_discarded = data_discarded.iloc[:, 0:9]
        print(f"Total candidates sites in {name}:", data.shape[0], flush=True)
        print(f"Total discarded sites in {name}:", data_discarded.shape[0], flush=True)
        stoptime = datetime.now()
        print(f"[{datetime.now()}] Candidates sites identification in {name} end. Elapsed time: {stoptime-starttime}.", flush=True)
        columns = ["Region", "Position", "Ref", "Strand", "Cov", "Qual", "Bases", "AllSubs", "Freq"]
        data.columns = columns
        data_discarded.columns = columns
        
        print(f"[{datetime.now()}] Features extraction in {name} start", flush=True)
        ohe = OneHotEncoder()
        ohe.fit(np.array(["A", "C", "G", "T"]).reshape(-1, 1))

        incompleates = []
        starttime_preds = datetime.now()
        total_extracted = 0
        
        df = data.query("Region != 'chrM'")
        del data
        df_dis = data_discarded.query("Region != 'chrM'")
        del data_discarded
        df_dis.to_csv(os.path.join(self.Results_path, f"{name}_discarded.txt"), sep="\t", index=None)
        del df_dis
        
        print(f"[{datetime.now()}] Loading table with tabix and pysam:", prefix, flush=True)
        start_time = datetime.now()
        srr = pysam.TabixFile(prefix+".gz")
        with tqdm(total=df.shape[0], position=0, leave=True) as pbar:
            for site in df.itertuples():
                start = int(site.Position) - ((101-1)/2)
                stop = int(site.Position) + ((101-1)/2)
                AGrna = eval(site.Bases)[2]/sum(eval(site.Bases))
                srr_interval = []
                for s in srr.fetch(site.Region, start-1, stop):
                    srr_interval.append(s.split("\t"))
                srr_interval = pd.DataFrame(srr_interval, columns=[f"{i}" for i in range(14)])
                if srr_interval.shape[0] == 101 and len(set(srr_interval["3"])) == 1:
                    intervals.append([site.Region, site.Position, site.Ref, site.Strand, 
                                      AGrna, site.Bases, start, stop, srr_interval.shape[0]])
                    total_extracted += 1
                    strand = site.Strand
                    seq = srr_interval["2"].values.reshape(-1,1)
                    seq_ohe = ohe.transform(seq).toarray().T
                    vects_freqs = []
                    strands = []
                    vects = []
                    for vect in srr_interval["6"]:
                        vect = np.array(eval(vect))
                        cov = sum(vect)
                        vect_freqs = vect / cov
                        vects_freqs.append(vect_freqs)
                        vects.append(vect)
                    vects_freqs = np.array(vects_freqs).T
                    vects = np.array(vects).T
                    encoded_site = pd.concat([pd.DataFrame(seq_ohe), pd.DataFrame(vects_freqs)])
                    encoded_site.reset_index(drop=True, inplace=True)
                    if strand == 0: 
                        encoded_site = pd.DataFrame(np.flip(encoded_site.values, axis=1))
                    features_matrices = pd.concat([features_matrices, encoded_site], axis=0)
                else:
                    incompleates.append(([site.Region, site.Position, site.Ref, site.Strand, 
                                          AGrna, site.Bases, start, stop, srr_interval.shape[0]]))
                pbar.update(1)
        del df       
        intervals = pd.DataFrame(intervals)
        incompleates = pd.DataFrame(incompleates)
        
        columns = ["Region", "Position", "Reference_Base", "Strand", "AG_Sub_Frequency", 
                             "Bases_Counts", "Start", "Stop", "TabixLen"]
        
        intervals.columns = columns
        incompleates.columns = columns
        incompleates.to_csv(os.path.join(self.Results_path, f"{name}_incompleates.txt"), sep="\t", index=None)
        del incompleates
        
        print(f"[{datetime.now()}] Total extracted sites in {name} sample: {total_extracted}.", flush=True)
        stop_time_global = datetime.now()
        print(f"[{datetime.now()}] Features extraction in {name} end. Elapsed time {datetime.now()-starttime_preds}.", flush=True)
        
        return intervals, features_matrices, name
                      
    def make_predictions(self)->None:
        
        model = tf.keras.models.load_model(os.path.join(self.Utilities_path, "REDInet_log_preprocessing_pos_zeros_snps_29_02_2024_22_25_45.h5"))
        
        if self.Multiprocessing == "yes":

            Inputs = []
            if type(self.Files_names) =="str":
                for Name in os.listdir(self.Data_path):
                    if Name.find(".tbi") == -1:
                        Inputs.append(Name.replace(".gz", ""))
            else:
                for Name in self.Files_names:
                    if Name.find(".gz") !=-1:
                    	if os.path.isfile(os.path.join(self.Data_path, Name)):
                            Inputs.append(Name.replace(".gz", ""))
                    	else:
                            sys.exit(f"{Name} file doesen't exist")
                    else:
                        sys.exit(f"{Name} is an incorrect file name")
      
            with Pool(processes) as pool:
                for metadata, features, Name in pool.map(self.extraction, Inputs):
                    
                    starttime_preds = datetime.now()
                    print(f"Inference on {Name} sample start", flush=True)
                    
                    features = self.from_2_to_3_dimensions(features)
                    features = self.log_preprocessing(features)
                    y_hat_proba = model.predict(features, batch_size=512, verbose=1)
                    
                    metadata["Not_Editing_Probability"] = 1.0 - y_hat_proba
                    metadata["Editing_Probability"] = y_hat_proba
                    predictions = []
                    for x in  y_hat_proba:
                        if x > 0.5:
                            predictions.append("Editing")
                        else:
                            predictions.append("Not_Editing")
                    metadata["Predicted_Class"] = predictions 
                    
                    metadata.to_csv(os.path.join(self.Results_path, f"{Name}_predictions.txt"), sep="\t", index=None)
                    print(f"[{datetime.now()}] Inference on {Name} sample ended. Elapsed time {datetime.now()-starttime_preds}.", flush=True)
                    
        else:
            
            if type(self.Files_names) =="str":
                Names = os.listdir(self.Data_path)
            else:
                Names = self.Files_names
    
            for Name in Names:
                if Name.find(".gz") !=-1:
                    if os.path.isfile(os.path.join(self.Data_path, Name)):

                        metadata, features, Name = self.extraction(Name.replace(".gz", ""))

                        starttime_preds = datetime.now()
                        print(f"Inference on {Name} sample start", flush=True)

                        features = self.from_2_to_3_dimensions(features)
                        features = self.log_preprocessing(features)
                        y_hat_proba = model.predict(features, batch_size=512, verbose=1)

                        metadata["Not_Editing_Probability"] = 1.0 - y_hat_proba
                        metadata["Editing_Probability"] = y_hat_proba

                        predictions = []
                        for x in  y_hat_proba:
                            if x > 0.5:
                                predictions.append("Editing")
                            else:
                                predictions.append("Not_Editing")
                        metadata["Predicted_Class"] = predictions

                        metadata.to_csv(os.path.join(self.Results_path, f"{Name}_predictions.txt"), sep="\t", index=None)
                        print(f"[{datetime.now()}] Inference on {Name} sample ended. Elapsed time {datetime.now()-starttime_preds}.", flush=True)
                    else:
                        sys.exit(f"{Name} file doesen't exist")
                else:
                    sys.exit(f"{Name} is an incorrect file name")


s = ("A-to-I RNA editing identification in RNAseq data from REDItools tabix-indexed files.\n"
     "\tOptions.\n"
     "\t--I: Input files folder path.\n"
     "\t--O: Output files folder path.\n"
     "\t--C: Minimun bases coverage.\n"
     "\t--F: Minimum AG substitutions rate.\n"
     "\t--M: Minimum number of detected guanosines.\n"
     "\t--N: List of files to be analyzed.\n"
     "\t--P: Use of multiprocessing.\n")

parser = argparse.ArgumentParser(description=s)
parser.add_argument("--I", type=str ,default="default", help=("Absolute path to the folder containing the files to be analyzed."
                                                              "By default is the Data folder in this package."))
parser.add_argument("--O", type=str ,default="default", help=("Absolute path to the folder in which the results are going to be stored."
                                                              "By default is the Results folder in this package."))
parser.add_argument("--C", type=int ,default=50, help=("Minimum number of RNAseq reads mapping in the positions to be analyzed."
                                                       "Value must be a natural number greater than zero."
                                                       "By defalult the value is 50."))
parser.add_argument("--F", type=float ,default=0.01, help=("Mimun percentage of AG substitutions in the positions to be analyzed."
                                                           "Value must be a floating point number grater than zero."
                                                           "By default the value is 0.01."))
parser.add_argument("--M", type=int, default=3, help=("Minimum detected guanosine number in place of adenines in the positions to be analyzed."
                                                      "Value must be a natural number greater than zero."
                                                      "By default the value is 3."))
parser.add_argument("--N", type=str, default="default", help=("Compleate list of .gz files to be analyzed."
                                                              "List must be included in box brackets and names must be separated by commas without spaces."
                                                              "Es:[data1.gz,data2.gz,data3.gz]."
                                                              "By default the tool analyze all the .gz files in the input directory"))
parser.add_argument("--P", type=str, default="no", choices=["yes", "no"], help=("Choose whether to work with multiprocessing or not." 
                                                                                "Possible choices are: yes or no"
                                                                                "By default is no"))
args = parser.parse_args()

utilities_path = os.path.dirname(os.path.abspath(__file__))

if args.I == "default":
    data_path = utilities_path.replace("Utilities", "Data")
else:
    data_path = args.I
    if not os.path.isdir(data_path):
        sys.exit("Invalid --I value.")
        
if args.O == "default":
    results_path = utilities_path.replace("Utilities", "Results")
else:
    results_path = args.O
    try:
        if not os.path.isdir(results_path):
            os.mkdir(results_path)
    except:
        sys.exit("Invalid --O value.")

coverage_threshold = args.C
if coverage_threshold <= 0:
    sys.exit("Invalid --C value")
frequence_threshold = args.F
if frequence_threshold <= 0:
    sys.exit("Invalid --F value.")
adenines_min = args.M
if adenines_min <= 0:
    sys.exit("Invalid --M value.")

multiprocessing = args.P

if multiprocessing == "yes":
    processes = multiprocessing.cpu_count()
else:
    processes = 1

if args.N != "default":
    files_names = args.N.replace("[", "").replace("]", "").split(",")
    state = False
    for file_name in files_names:
        if not os.path.isfile(os.path.join(data_path, file_name)):
            state = True
            print(f"{file_name} is missing in the input folder.", flush=True)
        if not os.path.isfile(os.path.join(data_path, file_name+".tbi")):
            state = True
            print(f"{file_name+".tbi"} is missing in the input folder.", flush=True)
    if state:
        sys.exit("Change the files list or upload the missing files in the input directory or change dhe input directory.")
    del state
else:
    gzs = glob(os.path.join(data_path, "*.gz"))
    if not gzs:
        sys.exit("No .gz file in the input folder.")
    else:
        for gz in gzs:
            if not os.path.isfile(os.path.join(data_path, gz+".tbi")):
                sys.exit(f"Missing {gz}.tbi file in the input folder.")
    del gzs
    
inference(coverage_threshold, frequence_threshold, adenines_min, multiprocessing, 
          processes, utilities_path, data_path, results_path, files_names).make_predictions()

print("\nExecution end.", flush=True)
