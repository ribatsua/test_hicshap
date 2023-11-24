import json
import shap
import h5py
import os
import numpy as np
from data_process import process_data, load_dataset, load_data_to_shap
from cnn_utils import  Simple1DCNN, train
from tqdm import tqdm

np.random.seed(123)


def load_config(config_path):
    with open(config_path, 'r') as config_file:
        config = json.load(config_file)
    return config

def run_shap(model, X, y ,data_dir,group):
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
    
    save_path = os.path.join(data_dir,group)
    file_path= os.path.join(save_path,'ori_shap_values.h5')
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    with h5py.File(file_path, 'w') as hf:
        grp = hf.create_group(group)
        for i in range(shap_arr.shape[2]):
            # cell_data = shap_arr[i].T  
            bin_data=shap_arr[:,:,i] 
            # grp.create_dataset(f'cell_{i + 1}', data=cell_data)
            grp.create_dataset(f'bin_{i + 1}', data=bin_data)
        label_shap=grp.create_group('label')
        label_shap.create_dataset('cell_type', data=[l.encode('utf8') for l in y['cell_type']],
                            dtype=h5py.special_dtype(vlen=str))
        label_shap.create_dataset('cell_id', data=[l.encode('utf8') for l in y['cell_id']],
                            dtype=h5py.special_dtype(vlen=str))
    
    hf.close()
    print("SHAP_value save to:",file_path)





def main():
    config_path='/home/python/higashi/notebook/config.json'
    config = load_config(config_path)
    data_dir = config['data_dir']
    mode = config['mode']
    group = config['group']
    epochs = config['epochs']
    batch_size = config['batch_size']
    lr = config['lr']
    chrom_list= config['chrom_list']
    smote = config['smote']
    print("data_dir:",data_dir,"\n","mode:",mode,"\n","group:",group,"\n","epochs:",epochs,"batch_size:",batch_size,"lr:",lr,"\n","chrom_list",chrom_list)

    data,label=process_data(data_dir, mode, group,chrom_list,smote)
    dataset=load_dataset(data,label)
    model = Simple1DCNN(data,label)
    model=train(model,dataset,epochs,batch_size,lr)

    X, y=load_data_to_shap(data_dir, mode, group,chrom_list,model)
    run_shap(model,X, y ,data_dir,group)


if __name__ == '__main__':
    main()







    


