import tensorflow as tf
from tensorflow.keras import layers
import pickle
import numpy as np
import matplotlib.pyplot as plt
import HackatonUtils as utils


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    model = tf.keras.models.load_model('./TrainedNanoNet')
    model.summary()
    # predict
    import os
    input_6dlb = utils.generate_input(ref_path)

    predict_dist, _, _, _ = model.predict(np.asarray([input_6dlb]))

