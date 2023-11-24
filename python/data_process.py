import pandas as pd
import numpy as np
import os
import pickle
import h5py
from tqdm import tqdm
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
import smote_variants as sv
from sklearn.metrics import classification_report

class MyDataset():
    def __init__(self,data_dir,mode,group,chrom_list):
        data_path = None
        label_path = None
        X = None
        y = None
        for file in os.listdir(data_dir):
            if mode in file.lower() and file.endswith('.hdf5'):
            # if mode in file.lower() and file.endswith('.h5'):
                data_path = os.path.join(data_dir, file)
                X=h5py.File(data_path, 'r')
                print(f"Data file to load: {data_path}")
            if 'label' in file.lower() and file.endswith('.pickle'):
                label_path = os.path.join(data_dir, file)
                y=pickle.load(open(label_path, "rb"))
                print(f"Label file to load: {label_path}")
        cell_keys = list(filter(lambda key: 'cell_' in key, X[group].keys()))
        cell_list_len=len(cell_keys)
        #datasets = [X[group][f'cell_{i}'][()] for i in range(cell_list_len)]
        chr_list=[chrom.decode('utf-8') for chrom in X['compartment_raw']['bin']['chrom']]
        # chr_id=[]
        # for chrom in chrom_list:
        # 	chr_id.extend([index for (index,value) in enumerate(chr_list) if value == chrom])
        chr_id = [index for index, value in enumerate(chr_list) if value in chrom_list]
        # datasets = [X[group][f'cell_{i}'][chr_id] for i in range(cell_list_len)]
        ##多个pc时 [chr_id,0]即取pc1的值
        datasets = [X[group][f'cell_{i}'][chr_id,0] for i in range(cell_list_len)]
        data = pd.DataFrame({f'cell_{i}': dataset for i, dataset in enumerate(datasets)})
        #data = pd.DataFrame({cell_key: X[group][cell_key][:][np.isin(X[group]['bin']['chrom'][:], chrom_list)] for cell_key in cell_keys})
        label=pd.DataFrame()
        label['cell_id']= data.columns
        label['cell_type']=y['cell_type']

        self.data = data.T
        self.label = label
    def smote_data(self):
        num_rows_to_append = 1000
        df=self.data
        label=self.label
        original_index=df.index
        columns = df.columns
        zero_rows = pd.DataFrame([[0] * len(columns)] * num_rows_to_append, columns=columns,
                                index=['cell_{}'.format(i) for i in range(len(original_index), len(original_index) + num_rows_to_append)])
        df = pd.concat([df, zero_rows])

        additional_rows = pd.DataFrame({'cell_id': [f'cell_{i}' for i in range(len(label), len(label) + num_rows_to_append)],'cell_type': ['zero'] * num_rows_to_append})
        df_label = pd.concat([label, additional_rows], ignore_index=True)
        #df_label['cell_type'].value_counts().sort_index()
        label_keys = df_label['cell_type'].value_counts().sort_index().index.tolist()
        num_values = range(len(label_keys))
        label_num_dict = dict(zip(label_keys, num_values))

        X=(df)
        X.index=df_label['cell_type']
        y=(X.index).map(label_num_dict)
        oversampler = sv.MulticlassOversampling(oversampler='distance_SMOTE', oversampler_params={'n_jobs':2, 'random_state':42})

        # X_samp and y_samp contain the oversampled dataset
        X_samp, y_samp= oversampler.sample(X, y)
        df_smote = pd.DataFrame(data=X_samp, columns=df.columns)
        df_smote = df_smote.rename(index=lambda s: 'cell_' + str(s))
        label_num_dict_inv = {v: k for k, v in label_num_dict.items()}
        y_samp_str = []
        for num in y_samp:
            y_samp_str.append(label_num_dict_inv[num])
        df_smote_label=pd.DataFrame()
        df_smote_label['cell_id']=df_smote.index
        df_smote_label['cell_type']=y_samp_str
        df_filtered = df_smote.loc[~(df_smote == 0).all(axis=1)]
        df_label_filtered = df_smote_label[df_smote_label['cell_type'] != 'zero']
        self.data_smote = df_filtered
        self.label_smote = df_label_filtered
        print(self.label_smote['cell_type'].value_counts().sort_index()) 
        return self.data_smote, self.label_smote

    def get_data_and_label(self):
        print("data_shape：",self.data.shape)
        print("label_shape：",self.label.shape)
        print(self.label['cell_type'].value_counts().sort_index())
        
        return self.data, self.label

def process_data(data_dir, mode, group,chrom_list,smote):
    if smote is False:
        dataset = MyDataset(data_dir, mode, group,chrom_list)
        return dataset.get_data_and_label()
    else:
        dataset = MyDataset(data_dir, mode, group,chrom_list)
        return dataset.smote_data()

def load_dataset(data,label):
    data_np=np.array(data)
    label_encoder = LabelEncoder()
    indexed_labels = label_encoder.fit_transform(label['cell_type'].values)
    onehot_encoder = OneHotEncoder(sparse=False)
    encoded_labels = onehot_encoder.fit_transform(indexed_labels.reshape(-1, 1))
    encoded_labels_df = pd.DataFrame(encoded_labels, columns=[f'_{i}' for i in range(encoded_labels.shape[1])])
    X_train, X_test, y_train, y_test = train_test_split(data_np, encoded_labels_df, test_size=0.2, random_state=123)
    X_train_3d = np.expand_dims(X_train, axis=2)
    X_test_3d = np.expand_dims(X_test, axis=2)
    print("X_train.shape:", X_train_3d.shape)
    print("y_train.shape:", y_train.shape)
    print()
    print("X_test.shape:", X_test_3d.shape)
    print("y_test.shape:", y_test.shape)
    dataset=(X_train_3d, y_train, X_test_3d ,y_test)
    return dataset
def load_data_to_shap(data_dir, mode, group,chrom_list,model):
    dataset = MyDataset(data_dir, mode, group,chrom_list)
    data,label=dataset.get_data_and_label()
    data_np=np.array(data)
    label_encoder = LabelEncoder()
    indexed_labels = label_encoder.fit_transform(label['cell_type'].values)
    onehot_encoder = OneHotEncoder(sparse=False)
    encoded_labels = onehot_encoder.fit_transform(indexed_labels.reshape(-1, 1))
    encoded_labels_df = pd.DataFrame(encoded_labels, columns=[f'_{i}' for i in range(encoded_labels.shape[1])])
    X = np.expand_dims(data_np, axis=2)
    y = np.expand_dims(encoded_labels_df,axis=2)
    print("------测试网络------")
    predictions = model.predict(X, batch_size=16)
    print(classification_report(y.argmax(axis=1),
    predictions.argmax(axis=1), target_names=label_encoder.classes_))
    return X, label







