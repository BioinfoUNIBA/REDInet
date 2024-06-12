# Copyright (c) 2023-2024 Pietro Luca Mazzacuva <pietroluca.mazzacuva@unicampus.it>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os, math
import numpy as np
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.python.platform.tf_logging import warn
from tensorflow import keras
from keras import Input, Model, Sequential
from keras.layers import Conv1D, Activation, Multiply, BatchNormalization, Add, MaxPooling1D, Flatten, Dense
from keras.initializers import Orthogonal
from keras.regularizers import L2, L1L2
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

class Utils():
    
    def __init__(self, mode=None):
        self.mode = mode

    def from_2_to_3_dimensions(self, pandas_df, n_features):

        time_step = pandas_df.shape[1]
        n_matrices = int(pandas_df.shape[0]/n_features)

        return np.transpose(np.reshape(pandas_df.to_numpy(), (n_matrices, n_features, time_step)), (0, 2, 1))

    def log_preprocessing(self, np_array):

        log_range = (-13.28771238, 0.000144262)
        np_array[:, :, 4:] = np.log2(np_array[:, :, 4:]+0.0001)
        np_array[:, :, 4:] = (np_array[:, :, 4:] - log_range[0]) / (log_range[1] - log_range[0])

        return np_array


    def make_confusion_matrix(self, y_true, y_pred,
                              group_names=None,
                              categories='auto',
                              count=True,
                              percent=True,
                              cbar=True,
                              xyticks=True,
                              xyplotlabels=True,
                              sum_stats=True,
                              figsize=None,
                              cmap='Blues',
                              title=None,
                              path=None):

        '''
        ###############################################################################################
        CITATION: taken and modified from: https://github.com/DTrimarchi10/confusion_matrix/blob/master/cf_matrix.py
        ###############################################################################################
        '''
                
        group_names = ["True Neg", "False Pos", "False Neg", "True Pos"]
        categories = ["No-Editing", "Editing"]

        predictions = []
        for x in y_pred:
            if x > 0.5:
                predictions.append(1)
            else:
                predictions.append(0)
        y_pred = np.array(predictions)
        
        cf = confusion_matrix(y_true,  y_pred)

        blanks = ['' for i in range(cf.size)]

        if group_names and len(group_names)==cf.size:
            group_labels = ["{}\n".format(value) for value in group_names]
        else:
            group_labels = blanks

        if count:
            group_counts = ["{0:0.0f}\n".format(value) for value in cf.flatten()]
        else:
            group_counts = blanks

        if percent:
            group_percentages = ["{0:.2%}".format(value) for value in cf.flatten()/np.sum(cf)]
        else:
            group_percentages = blanks

        box_labels = [f"{v1}{v2}{v3}".strip() for v1, v2, v3 in zip(group_labels,group_counts,group_percentages)]
        box_labels = np.asarray(box_labels).reshape(cf.shape[0],cf.shape[1])

        if sum_stats:
            accuracy  = np.trace(cf) / float(np.sum(cf))

            if len(cf)==2:

                specificity = cf[0, 0]/sum(cf[0,:])
                precision = cf[1,1]/sum(cf[:,1])
                recall = cf[1,1]/sum(cf[1,:])
                f1_score = 2*precision*recall/(precision+recall)
                balanced_accuracy = (recall+specificity)/2
                g_mean = math.sqrt(precision*specificity)
                stats_text = "\n\nBalanced Accuracy={:0.3f}\n Precision={:0.3f}\nRecall={:0.3f}\nF1-Score={:0.3f}\nSpecificity={:0.3f}\nG-Mean={:0.3f}".format(balanced_accuracy, precision,                                                                                                                                                                                    recall, f1_score, specificity, g_mean)
            
            else:
                stats_text = "\n\nAccuracy={:0.3f}".format(accuracy)
        else:
            stats_text = ""

        if figsize==None:
            figsize = plt.rcParams.get('figure.figsize')

        if xyticks==False:
            categories=False

        plt.figure(figsize=figsize)
        sn.heatmap(cf,annot=box_labels,fmt="",cmap=cmap,cbar=cbar,xticklabels=categories,yticklabels=categories)

        if xyplotlabels:
            plt.ylabel('True label')
            plt.xlabel('Predicted label' + stats_text)
        else:
            plt.xlabel(stats_text)

        if title:
            plt.title(title)

        if path:
            plt.tight_layout()
            plt.savefig(path)
        plt.show()
            
    def make_learning_curve(self, history, path, name):
        f, ax = plt.subplots(1,2, figsize=(15,5))
        ax[0].plot(history.history["loss"], label="Train-Set Loss")
        ax[0].plot(history.history["val_loss"], label="Validation-Set Loss")
        ax[0].set_xlabel("Epochs")
        ax[0].set_ylabel("Loss")
        ax[0].set_title("Loss in Training and Validation Datasets")
        ax[0].set_yscale("log")
        ax[0].legend()
        ax[1].plot(history.history["accuracy"], label="Train-Set Accuracy")
        ax[1].plot(history.history["val_accuracy"], label="Validation-Set accuracy")
        ax[1].set_xlabel("Epochs")
        ax[1].set_ylabel("Accuracy")
        ax[1].set_ylim(-0.1 ,1.1)
        ax[1].set_title("Accuracy in Training and Validation Datasets")
        ax[1].legend()
        plt.tight_layout()
        plt.savefig(os.path.join(path, f"{name}_Loss_Accuracy.tiff"))
        plt.show()

class TCN():

    def __init__(self, input_shape=(101, 8), out_channels=35, up_out_channels=35, out_channels_mul = 3, max_rate=32, 
                 kernel_size=2, gain=1.0, l2=0.01, momentum=0.2, epsilon=2e-5, strides=1, pool_size=2, n_units=10):
        
        self.input_shape = input_shape
        self.out_channels = out_channels
        self.up_out_channels = up_out_channels
        self.out_channels_mul = out_channels_mul
        self.max_rate = max_rate
        self.kernel_size = kernel_size
        self.gain = gain
        self.l2 = l2
        self.momentum = momentum
        self.epsilon = epsilon
        self.strides = strides
        self.pool_size = pool_size
        self.n_units = n_units

    def dilated_conv_block(self, res_number):
        dilated_block = Sequential(name = f"Residual_Unit_{res_number}_Dilated_Convolutions_Block")
        for rate in [int(pow(2, i)) for i in range(0, int(np.log2(self.max_rate))+1, 1)]:
            if rate != self.max_rate:
                dilated_block.add(Conv1D(filters=self.out_channels , kernel_size=self.kernel_size, padding = "causal", activation="relu", 
                dilation_rate = rate, kernel_initializer = Orthogonal(gain=self.gain), name = f"Residual_Unit_{res_number}_Conv1D_DilationRate_{rate}"))
            else:
                dilated_block.add(Conv1D(filters=self.out_channels , kernel_size=self.kernel_size, padding = "causal", dilation_rate = rate,
                kernel_initializer = Orthogonal(gain=self.gain), name = f"Residual_Unit_{res_number}_Conv1D_DilationRate_{rate}"))
        return dilated_block

    def residual_unit(self, res_number):
        Input = tf.keras.Input(shape = (self.input_shape[0], self.out_channels), name=f"Residual_Unit_{res_number}_InputLayer")
        dilation_block = self.dilated_conv_block(res_number)(Input)
        sigmoid = Activation("sigmoid", name=f"Residual_Unit_{res_number}_Sigmoid_ActivationLayer")(dilation_block)
        tanh = Activation("tanh", name=f"Residual_Unit_{res_number}_Tanh_ActivationLayer")(dilation_block)
        gated_activation = Multiply(name=f"Residual_Unit_{res_number}_Gated_Activation_Unit")([sigmoid, tanh])
        conv = Conv1D(filters=self.up_out_channels , kernel_size=1, padding = "causal", strides=1,
                      kernel_initializer = Orthogonal(gain=self.gain), name=f"Residual_Unit_{res_number}_Conv1D_1X1")(gated_activation)
        norm = BatchNormalization(momentum=self.momentum, epsilon=self.epsilon, center=True, scale=True,
                                  beta_initializer="zeros",gamma_initializer="ones", name=f"Residual_Unit_{res_number}_BatchNorm")(conv)
        relu = Activation("relu", name=f"Residual_Unit_{res_number}_ReLU_ActivationLayer")(norm)
        res_connection = Add(name=f"Residual_Unit_{res_number}_Residual_Connection")([Input, relu])
        res_unit = Model(inputs=Input,  outputs=[res_connection, norm], name=f"Residual_Unit_{res_number}")
        return res_unit

    def maxpool_block(self):
        block = Sequential(name="MaxPooling_Block")
        block.add(Conv1D(filters=self.up_out_channels*self.out_channels_mul, kernel_size=self.kernel_size, padding = "causal", strides=self.strides,
                         kernel_initializer = Orthogonal(gain=self.gain), name="MaxPooling_Block_First_Conv1D"))
        block.add(MaxPooling1D(pool_size=self.pool_size, name="MaxPooling_Block_First_MaxPool"))
        block.add(Conv1D(filters=self.up_out_channels*self.out_channels_mul, kernel_size=self.kernel_size, padding = "causal", strides=self.strides,
                         kernel_initializer = Orthogonal(gain=self.gain), name="MaxPooling_Block_Second_Conv1D"))
        block.add(tf.keras.layers.MaxPooling1D(pool_size=self.pool_size, name="MaxPooling_Block_Second_MaxPool"))
        block.add(BatchNormalization(momentum=self.momentum, epsilon=self.epsilon, center=True, scale=True,
                                   beta_initializer="zeros",gamma_initializer="ones", name="MaxPooling_Block_BatchNorm"))
        block.add(Flatten(name="MaxPooling_Block_FlattenLayer"))
        return block

    def residual_stack(self):
        Input = tf.keras.Input(shape = (self.input_shape[0], self.out_channels), name="Residual_Units_Stack_InputLayer")
        for i in range(1, self.n_units+1, 1):
            if i == 1:
                globals()[f"res{i}"], globals()[f"skip{i}"] = self.residual_unit(i)(Input)
            else:
                globals()[f"res{i}"], globals()[f"skip{i}"] = self.residual_unit(i)(globals()[f"res{i-1}"])
        add1 = Add(name="Skip_Connections_AddLayer")([globals()[f"skip{i}"] for i in range(1, self.n_units+1, 1)])
        relu = Activation("relu", name="Skip_Connections_ReLU_ActivationLayer")(add1)
        conv_1x1_1 = Conv1D(filters=self.up_out_channels, kernel_size=1, padding = "causal", strides=1, activation = "relu",
                            kernel_initializer = Orthogonal(gain=self.gain), name="Skip_Connections_First_ConvID_1X1")(relu)
        conv_1x1_2 = Conv1D(filters=self.up_out_channels, kernel_size=1, padding = "causal", strides=1, activation = "relu",
                            kernel_initializer = Orthogonal(gain=self.gain), name="Skip_Connections_Second_ConvID_1X1")(conv_1x1_1)
        add2 = Add(name="Residual_Units_Stack_AddLayer")([conv_1x1_2, globals()[f"res{self.n_units}"]])      
        stack = keras.Model(inputs=Input, outputs=add2, name="Residual_Units_Stack")
        return stack

    def MLP(self):
        mlp = Sequential(name = "Multilayer_Perceptron") 
        mlp.add(Dense(self.up_out_channels*self.out_channels_mul, activation='relu', kernel_initializer = Orthogonal(gain=self.gain), name="Multilayer_Perceptron_First_DenseLayer"))
        mlp.add(BatchNormalization(momentum=self.momentum, epsilon=self.epsilon, center=True, scale=True,
                                   beta_initializer="zeros",gamma_initializer="ones", name="Multilayer_Perceptron_First_BatchNorm"))
        mlp.add(Dense(int((self.up_out_channels*self.out_channels_mul)/2), activation='relu', kernel_initializer = Orthogonal(gain=self.gain),
                      kernel_regularizer= L2(l2=self.l2), name="Multilayer_Perceptron_Second_DenseLayer"))
        mlp.add(BatchNormalization(momentum=self.momentum, epsilon=self.epsilon, center=True, scale=True,
                                   beta_initializer="zeros",gamma_initializer="ones", name="Multilayer_Perceptron_Second_BatchNorm"))
        mlp.add(Dense(int((self.up_out_channels*self.out_channels_mul)/4), activation='relu', kernel_initializer = Orthogonal(gain=self.gain),
                      kernel_regularizer= L2(l2=self.l2), name="Multilayer_Perceptron_Third_DenseLayer"))
        mlp.add(BatchNormalization(momentum=self.momentum, epsilon=self.epsilon, center=True, scale=True,
                                   beta_initializer="zeros",gamma_initializer="ones", name="Multilayer_Perceptron_Third_BatchNorm"))
        mlp.add(Dense(1, activation='sigmoid', name="Multilayer_Perceptron_Sigmoid"))
        return mlp

    def get_model(self):
        input = Input(shape=self.input_shape, name="REDInet_InputLayer")
        conv = Conv1D(filters = self.out_channels, kernel_size=self.kernel_size, padding = "same", activation="relu",
                      kernel_initializer = Orthogonal(gain=self.gain), name="Conv1D")(input)
        stack = self.residual_stack()(conv)
        block = self.maxpool_block()(stack)
        mlp = self.MLP()(block)
        redinet = Model(inputs=input, outputs=mlp, name="REDInet")
        return redinet
