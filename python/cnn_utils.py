from keras.models import Sequential
from keras.layers import Conv1D, MaxPooling1D, Flatten, Dense, Dropout
from keras.optimizers import Adam
from keras.callbacks import EarlyStopping

def Simple1DCNN(data,label):
    bin_num=data.shape[1]
    class_num=len(list(set(label['cell_type'])))
    model = Sequential()
    model.add(Conv1D(8,3,activation='relu',input_shape=(bin_num,1)))
    model.add(MaxPooling1D(2))
    model.add(Conv1D(16,3,activation='relu'))
    model.add(MaxPooling1D(2))
    model.add(Flatten())
    # model.add(Dense(16,activation='relu'))
    # model.add(Dropout(0.5))
    model.add(Dense(class_num, activation='softmax'))
    model.summary()
    return model

def train(model,dataset,epochs,batch_size,lr):
    Optimizer=Adam(learning_rate=lr)
    X_train_3d, y_train, X_test_3d ,y_test=dataset
    model.compile(loss='categorical_crossentropy', optimizer=Optimizer, 
                    metrics=['accuracy'])
    early_stopping = EarlyStopping(monitor='val_loss',  
                                mode='min',          
                                patience=5,          
                                verbose=1)

    model.fit(X_train_3d, y_train, epochs=epochs, batch_size=batch_size, validation_data=(X_test_3d, y_test), callbacks=[early_stopping])
    _, accuracy = model.evaluate(X_test_3d, y_test)
    print('Accuracy: %.2f%%' % (accuracy * 100))

    return model
