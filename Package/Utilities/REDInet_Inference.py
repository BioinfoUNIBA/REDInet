import numpy as np
import pandas as pd
import tensorflow as tf
import argparse, glob, gzip, multiprocessing, os, pysam, shlex, subprocess, sys, time
from datetime import datetime
from multiprocessing import Pool
from sklearn.preprocessing import OneHotEncoder
from tensorflow import keras
from tqdm import tqdm

class inference():
    
    def __init__(self, cov_threshold, AGfreq_threshold, AG_min, Multiprocessing, Assembly, Max_missings,
                 Imputations, Processes, Utilities_path, Data_path, Results_path, Files_names):
        self.cov_threshold = cov_threshold
        self.AGfreq_threshold = AGfreq_threshold
        self.AG_min =  AG_min
        self.Multiprocessing = Multiprocessing
        self.Assembly = Assembly
        self.Max_missings = Max_missings
        self.Imputations = Imputations
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
        print(f"[{datetime.now()}] Candidates sites identification in {name} sample start\n", flush=True)
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
                    print(f"\t[{datetime.now()}] Sites evaluated in {name} sample: {c}", flush=True)
        print(f"[{datetime.now()}] Total evaluated rows in {name} sample: ", c, flush=True)
        
        columns = ["Region", "Position", "Ref", "Strand", "Cov", "Qual", "Bases", "AllSubs", "Freq"]

        data_discarded =  pd.DataFrame(data_discarded)
        if data_discarded.shape[0] > 0:
            
            data_discarded = data_discarded.iloc[:, 0:9]
            data_discarded.columns = columns
            df_dis = data_discarded.query("Region != 'chrM'")
            del data_discarded
            print(f"[{datetime.now()}] Total discarded sites in {name} sample:", df_dis.shape[0], flush=True)
            df_dis.to_csv(os.path.join(self.Results_path, f"{name}_discarded.txt"), sep="\t", index=None)
            del df_dis

        data = pd.DataFrame(data)
        print(f"[{datetime.now()}] Total candidates sites in {name} sample:", data.shape[0], flush=True)
        print(f"[{datetime.now()}] Candidates sites identification in {name} sample end. Elapsed time: {datetime.now()-starttime}.", flush=True)
        
        if data.shape[0] > 0:
            
            data = data.iloc[:, 0:9]
            data.columns = columns
            print(f"[{datetime.now()}] Features extraction in {name} sample start", flush=True)
            starttime_preds = datetime.now()
            ohe = OneHotEncoder()
            ohe.fit(np.array(["A", "C", "G", "T"]).reshape(-1, 1))
            incompleates = []
            total_extracted = 0
            df = data.query("Region != 'chrM'")
            del data
            print(f"[{datetime.now()}] Loading table with tabix and pysam:", prefix, flush=True)
            start_time = datetime.now()
            srr = pysam.TabixFile(prefix+".gz")
            print(f"[{datetime.now()}] Loading human reference genome with pysam:", prefix, flush=True)
            
            if self.Assembly == "GRCh37":
                genome = pysam.FastaFile(os.path.join(self.Utilities_path, "GRCh37.primary_assembly.genome.fa"))
            else:
                genome = pysam.FastaFile(os.path.join(self.Utilities_path, "GRCh38.primary_assembly.genome.fa"))
                
            with tqdm(total=df.shape[0], position=0, leave=True) as pbar:
                
                for site in df.itertuples():
                    start = int(site.Position) - int((101-1)/2)
                    stop = int(site.Position) + int((101-1)/2)
                    AGrna = eval(site.Bases)[2]/sum(eval(site.Bases))
                    srr_interval = []
                    for s in srr.fetch(site.Region, start-1, stop):
                        info = s.split("\t")
                        srr_interval.append([info[0], info[1], info[2], info[3], info[6]])
                    srr_interval = pd.DataFrame(srr_interval, columns=[f"{i}" for i in range(5)])
                    srr_interval.iloc[:, 1] = srr_interval.iloc[:, 1].astype("int32")
                    n_missings = 101-int(srr_interval.shape[0])
                    temp = []
                   
                    if self.Imputations == "yes" and n_missings > 0 and n_missings <=self.Max_missings :
              
                        freqs_impud = {"A":"[1,0,0,0]", "C":"[0,1,0,0]", "G":"[0,0,1,0]", "T":"[0,0,0,1]"}
                        complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
                        missings = [i  for i in range(start, stop+1, 1) if i not in srr_interval.iloc[:, 1].tolist()]

                        if all(int(element) == 1 for element in srr_interval["3"]):
                            
                            for position in missings:
                                base = genome.fetch(site.Region, int(position)-1, int(position))
                                value = freqs_impud.get(base)
                                if value:
                                    temp.append([site.Region, position, base, "1", value])
                                else:
                                    pass
                                
                        elif all(int(element) == 0 for element in srr_interval["3"]):

                            for position in missings:
                                base = genome.fetch(site.Region, int(position)-1, int(position))
                                complement_base = complement.get(base)
                                if complement_base:
                                    value = freqs_impud.get(complement_base)
                                    if value:
                                        temp.append([site.Region, position, complement_base, "0", value])           
                                    else:
                                        pass
                                else:
                                    pass            
                        else:
                            pass

                    if len(temp) > 0 :

                        temp = pd.DataFrame(temp, columns=[f"{i}" for i in range(5)])
                        srr_interval = pd.concat([srr_interval, temp], axis=0)
                        srr_interval.sort_values(by=["1"], inplace=True, ignore_index=True)
                        srr_interval.reset_index(drop=True, inplace=True)

                    if srr_interval.shape[0] == 101 and len(set(srr_interval["3"])) == 1:

                        intervals.append([site.Region, site.Position, site.Ref, site.Strand, 
                                          AGrna, site.Bases, start, stop, n_missings])
                        total_extracted += 1
                        strand = site.Strand
                        seq = srr_interval["2"].values.reshape(-1,1)
                        seq_ohe = ohe.transform(seq).toarray().T
                        vects_freqs = []
                        strands = []
                        vects = []
                        for vect in srr_interval["4"]:
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
                                              AGrna, site.Bases, start, stop, n_missings]))
                    pbar.update(1)
                    
            del df       
                
            columns = ["Region", "Position", "Reference_Base", "Strand", "AG_Sub_Frequency", 
                                 "Bases_Counts", "Start", "Stop", "N. Missing Nucleotides"]
                
            incompleates = pd.DataFrame(incompleates)    
            
            if incompleates.shape[0] > 0:    
                incompleates.columns = columns
                print(f"[{datetime.now()}] Total sequence with eccessive missing nucleotides in {name} sample: {incompleates.shape[0]}", flush=True)
                incompleates.to_csv(os.path.join(self.Results_path, f"{name}_incompleates.txt"), sep="\t", index=None)
                del incompleates
            else:
                print(f"[{datetime.now()}] No sequences with gap for candidates sites in {name} sample", flush=True)
                
            intervals = pd.DataFrame(intervals)
            if intervals.shape[0]>0:
                intervals.columns = columns

            print(f"[{datetime.now()}] Total extracted sites in {name} sample: {total_extracted}.", flush=True)
            print(f"[{datetime.now()}] Features extraction in {name} sample end. Elapsed time {datetime.now()-starttime_preds}.", flush=True)
        
        else:
            del data
            intervals = pd.DataFrame()
            features_matrices = pd.DataFrame()
            
        return intervals, features_matrices, name
                      
    def make_predictions(self)->None:
        
        model = tf.keras.models.load_model(os.path.join(self.Utilities_path, "REDInet.h5"))
        
        if self.Multiprocessing == "yes":

            Inputs = []
            
            for Name in self.Files_names:
                Inputs.append(Name.replace(".gz", ""))
                
            with Pool(self.Processes) as pool:
                for metadata, features, Name in pool.map(self.extraction, Inputs):
                    
                    starttime_preds = datetime.now()
                    print(f"[{datetime.now()}] Inference on {Name} sample start", flush=True)

                    if metadata.shape[0]>0:
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
                        print(f"[{datetime.now()}] No Inference on {Name} sample could be carried out. Elapsed time {datetime.now()-starttime_preds}.", flush=True)
        else:
            
            Names = self.Files_names
            
            for Name in Names:
                
                metadata, features, Name = self.extraction(Name.replace(".gz", ""))
                
                starttime_preds = datetime.now()
                print(f"Inference on {Name} sample start", flush=True)
                
                if metadata.shape[0]>0:
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
                    print(f"[{datetime.now()}] No Inference on {Name} sample could be carried out. Elapsed time {datetime.now()-starttime_preds}.", flush=True)

                                                                        


s = ("A-to-I RNA editing identification in RNAseq data from REDItools tabix-indexed files.\n"
     "\tOptions.\n"
     "\t--I: Input files folder path.\n"
     "\t--O: Output files folder path.\n"
     "\t--C: Minimun bases coverage.\n"
     "\t--F: Minimum AG substitutions rate.\n"
     "\t--M: Minimum number of detected guanosines.\n"
     "\t--N: List of files to be analyzed.\n"
     "\t--P: Use of multiprocessing.\n"
     "\t--G: Accept sequences'gaps.\n"
     "\t--A: Human genome assembly to adopt.\n"
     "\t--S: Maximum number of gaps tollerated.\n")

parser = argparse.ArgumentParser(description=s)
parser.add_argument("--I", type=str, default="default", help=("Absolute path to the folder containing the files to be analyzed."
                                                              "By default is the Data folder in this package."))
parser.add_argument("--O", type=str, default="default", help=("Absolute path to the folder in which the results are going to be stored."
                                                              "By default is the Results folder in this package."))
parser.add_argument("--C", type=int, default=50, help=("Minimum number of RNAseq reads mapping in the positions to be analyzed."
                                                       "Value must be a natural number greater than zero."
                                                       "By defalult the value is 50."))
parser.add_argument("--F", type=float, default=0.01, help=("Mimun percentage of AG substitutions in the positions to be analyzed."
                                                           "Value must be a floating point number grater than zero."
                                                           "By default the value is 0.01."))
parser.add_argument("--M", type=int, default=3, help=("Minimum detected guanosine number in place of adenines in the positions to be analyzed."
                                                      "Value must be a natural number greater than zero."
                                                      "By default the value is 3."))
parser.add_argument("--N", type=str, default="default", help=("Compleate list of .gz files to be analyzed."
                                                              "List must be included in box brackets and names must be separated by commas without spaces."
                                                              "Es:[data1.gz,data2.gz,data3.gz]."
                                                              "By default the tool analyze all the .gz files in the input directory."))
parser.add_argument("--P", type=str, default="no", choices=["yes", "no"], help=("Choose whether to work with multiprocessing or not." 
                                                                                "Possible choices are: yes or no."
                                                                                "By default is no."))
parser.add_argument("--G", type=str, default="no", choices=["yes", "no"], help=("Choose to accept sequences with maximum 30 missing nucleotides." 
                                                                                "Possible choices are: yes or no."
                                                                                "By default is no."))
parser.add_argument("--A", type=str, default="GRCh37", choices=["GRCh37", "GRCh38"], help=("Human genome assembly to use in handling missing nucleotides." 
                                                                                           "Possible choices are: GRCh37 or GRCh38."
                                                                                           "By default is GRCh37."))
parser.add_argument("--S", type=int, default=10, help=("Number of missing nucleotides to be imputed" 
                                                       "Value must be an integer number grater than zero and equal or smaller than 30."
                                                       "By default the value is 10."))

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
try:
    if coverage_threshold <= 0:
        sys.exit("Invalid --C value")
except:
    sys.exit("Invalid --C value")
    
frequence_threshold = args.F
try:
    if frequence_threshold <= 0 and frequence_threshold >= 1.0:
        sys.exit("Invalid --F value.")
except:
    sys.exit("Invalid --F value.")
    
adenines_min = args.M
try:
    if adenines_min <= 0:
        sys.exit("Invalid --M value.")
except:
    sys.exit("Invalid --M value.")

max_missings = args.S
try:
    if max_missings < 0 and max_missings > 30:
        sys.exit("Invalid --S value.")
except:
    sys.exit("Invalid --S value.")

MultiProcessing = args.P

if MultiProcessing == "yes":
    processes = multiprocessing.cpu_count()
elif MultiProcessing == "no":
    processes = 1
else:
    sys.exit("Invalid --P value.")
    
imputations = args.G
if imputations != "yes" and imputations != "no":
    sys.exit("Invalid --G value.")
    
assembly = args.A
if assembly != "GRCh37" and assembly != "GRCh38":
    sys.exit("Invalid --A value.")

if args.N != "default":
    files_names = args.N.replace("[", "").replace("]", "").split(",")
    for file_name in files_names:
        if file_name.find(".gz") ==-1:
            sys.exit(f"{file_name} is an incorrect file name.")
        if not os.path.isfile(os.path.join(data_path, file_name)):
            sys.exit(f"{file_name} is missing in the input folder.")
        if not os.path.isfile(os.path.join(data_path, file_name+".tbi")):
            sys.exit(f"{file_name+'.tbi'} is missing in the input folder.")
else:
    files_names = [i.replace(f"{data_path}/", "") for i in glob.glob(os.path.join(data_path, "*.gz"))]
    if not files_names:
        sys.exit("No .gz file in the input folder.")
    else:
        for file_name in files_names:
            if not os.path.isfile(os.path.join(data_path, file_name+".tbi")):
                sys.exit(f"Missing {file_name}.tbi file in the input folder.")

inference(coverage_threshold, frequence_threshold, adenines_min, MultiProcessing, assembly, max_missings, 
          imputations, processes, utilities_path, data_path, results_path, files_names).make_predictions()

print("\n[{datetime.now()}] Execution end.", flush=True)