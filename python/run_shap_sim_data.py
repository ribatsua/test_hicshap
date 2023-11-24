from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd
from test import CNN_1D, train
import shap
import h5py
from tqdm import tqdm
from sklearn.metrics import classification_report

def load_sim_data(data,label):
    data_np=np.array(data)
    label_encoder = LabelEncoder()
    indexed_labels = label_encoder.fit_transform(label)
    onehot_encoder = OneHotEncoder(sparse=False)
    encoded_labels = onehot_encoder.fit_transform(indexed_labels.reshape(-1, 1))
    encoded_labels_df = pd.DataFrame(encoded_labels, columns=[f'_{i}' for i in range(encoded_labels.shape[1])])
    X_train, X_test, y_train, y_test = train_test_split(data_np.transpose(), encoded_labels_df, test_size=0.2, random_state=123)
    X_train_3d = np.expand_dims(X_train, axis=2)
    X_test_3d = np.expand_dims(X_test, axis=2)
    print("X_train.shape:", X_train_3d.shape)
    print("y_train.shape:", y_train.shape)
    print()
    print("X_test.shape:", X_test_3d.shape)
    print("y_test.shape:", y_test.shape)
    dataset=(X_train_3d, y_train, X_test_3d ,y_test)
    return dataset

def run_shap(model, data ,label):
    group='compartment_raw'
    data_np=np.array(data)
    label_encoder = LabelEncoder()
    indexed_labels = label_encoder.fit_transform(label)
    onehot_encoder = OneHotEncoder(sparse=False)
    encoded_labels = onehot_encoder.fit_transform(indexed_labels.reshape(-1, 1))
    encoded_labels_df = pd.DataFrame(encoded_labels, columns=[f'_{i}' for i in range(encoded_labels.shape[1])])
    X = np.expand_dims(data_np.transpose(), axis=2)
    y = np.expand_dims(encoded_labels_df,axis=2)
    print("------测试网络------")
    predictions = model.predict(X, batch_size=16)
    print(classification_report(y.argmax(axis=1),
    predictions.argmax(axis=1), target_names=label_encoder.classes_))
    background=X[np.random.choice(X.shape[0], 100, replace=False)]
    # a=X[0:343,:,:][np.random.choice(343, 50, replace=False)]
    # b=X[343:816,:,:][np.random.choice(473, 50, replace=False)]
    # c=X[816:1075,:,:][np.random.choice(259, 50, replace=False)]
    # background=np.concatenate((a,b,c),axis=0)
    e = shap.DeepExplainer(model, background)
    print("e.expected_value:",e.expected_value)
    out_list = []
    num_samples = np.shape(X)[0]
    print("num_samples:",num_samples)
    for sample in tqdm(range(num_samples)):
        # shap
        shap_values = e.shap_values(X[sample : sample + 1])
        out_list.append(shap_values)
    shap_arr = np.squeeze(np.array(out_list))
    # save_path=data_dir+'/'+group+'_shap_values_150_remove_chrX_ncb.h5'
    
    # save_path = os.path.join(data_dir,group)
    # file_path= os.path.join(save_path,'shap_values.h5')
    # if not os.path.exists(save_path):
    #     os.makedirs(save_path)
    file_path='sim_SAWB_shap_values_x3.h5'
    with h5py.File(file_path, 'w') as hf:
        grp = hf.create_group(group)
        for i in range(shap_arr.shape[2]):
            # cell_data = shap_arr[i].T  
            bin_data=shap_arr[:,:,i] 
            # grp.create_dataset(f'cell_{i + 1}', data=cell_data)
            grp.create_dataset(f'bin_{i + 1}', data=bin_data)
        label_shap=grp.create_group('label')
        label_shap.create_dataset('cell_type', data=[l.encode('utf8') for l in label],
                            dtype=h5py.special_dtype(vlen=str))
    
    hf.close()
    print("SHAP_value save to:",file_path)


if __name__ == '__main__':
    # data1_path = '/home/sim_data/meanA-1.csv'
    # data2_path = '/home/sim_data/var1.5.csv'
    data1_path = '/home/sim_data/meansA0.5.csv'
    data2_path = '/home/sim_data/meanwB2.csv'
    data3_path = '/home/sim_data/var1.5.csv'

    data1 = pd.read_csv(data1_path,index_col=0)
    data2 = pd.read_csv(data2_path,index_col=0)
    data3 = pd.read_csv(data3_path,index_col=0)
    a=data1.iloc[:,:-5]
    b=data2.iloc[:,:-5]
    c=data3.iloc[:,:-5]
    data=pd.merge(a,b,left_index=True,right_index=True)
    data=pd.merge(data,c,left_index=True,right_index=True)
    label1=['A']*a.shape[1]
    label2=['B']*a.shape[1]
    label3=['C']*a.shape[1]
    label=[]
    label.extend(label1)
    label.extend(label2)
    label.extend(label3)
    dataset=load_sim_data(data,label)
    model= CNN_1D(data,label)
    model = train(model,dataset,5,16,0.001)
    run_shap(model, data ,label)
