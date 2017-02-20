import tensorflow as tf
import numpy as np

def example_1():

	np.set_printoptions(suppress=True, precision=4)

	param_range = 50
	param_increment = 1
	data_points = 100

	x_data = np.arange(1,data_points+1).reshape(1,1,data_points).astype(np.float32)
	y_data = np.ones((param_range,param_range,data_points))

	for W_param in range(1,1+param_range):
		for c_param in range(6,6+param_range):
			y_data[W_param-1, c_param-6, :] = np.log(x_data) * W_param/param_increment + 10/(x_data**2) + c_param/param_increment*np.sin(x_data) 

	# for W_param in range(5):
	# 	for c_param in range(5):
	# 		y = np.exp(x_data) * W[W_param, c_param] + 10/(x_data**2) + c[W_param, c_param]*np.sin(x_data) 

	W = tf.Variable(tf.random_uniform((param_range,param_range,1), -1.0, 1.0))
	# b = tf.Variable(tf.random_uniform([5,5], -1.0, 1.0))
	c = tf.Variable(tf.random_uniform((param_range,param_range,1), -1.0, 1.0))
	y = tf.Variable(tf.ones([param_range,param_range,1], tf.float32))

	y = np.log(x_data[:,:,:]) * W + 10/(x_data[:,:,:]*2) + c*np.sin(x_data[:,:,:]) 

	# # [.1,.3,.4]

	loss = tf.reduce_mean(tf.square(y - y_data))
	optimizer = tf.train.FtrlOptimizer(1)
	train = optimizer.minimize(loss)

	init = tf.initialize_all_variables()
	
	sess = tf.Session()
	sess.run(init)
	print sess.run(W[0,0,0])
	print sess.run(c[0,0,0])
	print sess.run(loss)

	# for W_param in range(1,1+param_range):
	# 	for c_param in range(6,6+param_range):
	# 		step = 0
	# 		# y[W_param-1, c_param-6, :] = np.log(x_data[:,:,:]) * W[W_param-1, c_param-6, :] + 10/(x_data[:,:,:]*2) + c[W_param-1, c_param-6, :]*np.sin(x_data[:,:,:])
	# 		loss = tf.reduce_mean(tf.square(y[W_param-1, c_param-6, :] - y_data[W_param-1, c_param-6, :]))
	# 		train = optimizer.minimize(loss)
	# 		sess.run(train)
	# 		while sess.run(loss) > .01 and step < 2000:
	# 			sess.run(train)
	# 			step += 1
	# 			# if step % 20 == 0:
	# 		print(step, sess.run(W[W_param-1, c_param-6, :]), sess.run(c[W_param-1, c_param-6, :]), sess.run(loss))

	for step in range(25001):
		sess.run(train)
		if step % 20 == 0:
			print(step, sess.run(W), sess.run(c), sess.run(loss))

	print sess.run(W)
	print sess.run(c)

if __name__ == '__main__':
	example_1()
