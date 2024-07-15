import tensorflow as tf

print("TensorFlow version:", tf.__version__)

hello = tf.constant('Hello, TensorFlow!')
print(hello.numpy())

a = tf.constant(2)
b = tf.constant(3)
c = a + b

print("The result of the computation a + b is:", c.numpy())
