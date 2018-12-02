import tensorflow as tf

import numpy as np


class Network:
    def __init__(self, data):
        self.data = data

    def Network(self):
        X = tf.placeholder(dtype="tf.float32", shape=(1, 2))
        b = tf.placeholder(dtype="tf.float32", shape=(1))
        Y = tf.placeholder(dtype="tf.string", shape=(1))

        W1 = tf.placeholder(shape=(2, 1), type="tf.float32")

        layer1 = tf.matmul(X * W1)