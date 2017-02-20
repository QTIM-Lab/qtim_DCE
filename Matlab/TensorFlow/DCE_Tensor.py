from __future__ import division
import nifti_util
import numpy as np
import nibabel as nib
import tensorflow as tf
import scipy.optimize
import matplotlib.pyplot as plt
import math
import random
import os

global_AIF = []

def calc_DCE_properties_single(filepath, T1_tissue=1000, T1_blood=1440, relaxivity=.0045, TR=5, TE=2.1, scan_time_seconds=(11*60), hematocrit=0.45, injection_start_time_seconds=60, flip_angle_degrees=30, label_file=[], label_suffix=[], label_value=1, mask_value=0, mask_threshold=0, T1_map_file=[], T1_map_suffix='-T1Map', AIF_label_file=[],  AIF_values_file=[], AIF_mode='label_average', AIF_label_suffix=[], AIF_label_value=1, label_mode='separate', param_file=[], default_population_AIF=False, intial_fitting_function_parameters=[-2,-2], outputs=['ktrans','ve','auc'], outfile_prefix=''):

	flip_angle_radians = flip_angle_degrees*np.pi/180
	image = nifti_util.nifti_2_numpy(filepath)
	dimension = len(image.shape)
	if dimension > 4:
		print 'Error: Images greater than dimension 4 are currently not supported. Skipping this volume...'
		return []
	time_interval_seconds = float(scan_time_seconds / image.shape[dimension-1])

	AIF_label_image, label_image, T1_image = retreive_labels_from_files(filepath, label_file, label_mode, label_suffix, label_value, AIF_label_file, AIF_label_value, AIF_mode, AIF_label_suffix, T1_map_file, T1_map_suffix)

	if T1_image != [] and (image.shape[0:-1] != T1_image.shape):
		print 'T1 map and DCE image are not the same shape. T1 map processing will be skipped.'
		T1_image = []

	AIF = generate_AIF(scan_time_seconds, injection_start_time_seconds, time_interval_seconds, image, AIF_label_image, AIF_mode, dimension, AIF_label_value)

	if AIF == []:
		print 'Problem calculating AIF. Skipping this volume...'
		return []

	contrast_image = generate_contrast_agent_concentration(image, T1_tissue, TR, flip_angle_degrees, injection_start_time_seconds, relaxivity, time_interval_seconds, hematocrit)
	contrast_AIF = generate_contrast_agent_concentration(AIF, T1_tissue, TR, flip_angle_degrees, injection_start_time_seconds, relaxivity, time_interval_seconds, hematocrit, T1_blood=1440)

	# Note-to-self: implement 'output' parameter functionality. Any custom-specified outputs
	# will currently break the program.
	parameter_maps = np.zeros((contrast_image.shape[0:-1] + (len(outputs),)), dtype=float)
	parameter_maps = simplex_optimize(contrast_image, contrast_AIF, time_interval_seconds, mask_value, mask_threshold, intial_fitting_function_parameters)

	for param_idx, param in enumerate(outputs):
		nifti_util.save_numpy_2_nifti(parameter_maps[:,:,param_idx], filepath, outfile_prefix + param + '.nii.gz')

def retreive_labels_from_files(filepath, label_file, label_mode, label_suffix, label_value, AIF_label_file, AIF_label_value, AIF_mode, AIF_label_suffix, T1_map_file, T1_map_suffix):

	if label_file != []:
		label_image = nifti_util.nifti_2_numpy(label_file)
	elif label_suffix != []:
		split_path = str.split(filepath, '.nii')
		if os.path.isfile(split_path[0] + label_suffix + '.nii' + split_path[1]):
			label_image = nifti_util.nifti_2_numpy(split_path[0] + label_suffix + '.nii' + split_path[1])
		else:
			print "No labelmap found at provided label suffix. Continuing without..."
			label_image = []
	else:
		label_image = []

	# Note: specific label value functionality not yet added.

	if AIF_mode == 'label_average':
		if AIF_label_file != []:
			AIF_label_image = nifti_util.nifti_2_numpy(AIF_label_file)
		elif AIF_label_suffix != []:
			split_path = str.split(filepath, '.nii')
			AIF_label_image = nifti_util.nifti_2_numpy(split_path[0] + AIF_label_suffix + '.nii' + split_path[1])
		elif label_mode == 'separate':
			print 'No label found for this AIF. If AIF label is in the same file as ROI, change the label_mode parameter to \'combined\'. Skipping this volume...'
			AIF_label_image = []
		elif label_mode == 'combined':
			if label_file != []:
				AIF_label_image = np.copy(label_image)
				AIF_label_image[label_image != AIF_label_value] = 0
			else:
				print 'No label found for this AIF. If the AIF label is in a separate file from the ROI, change the label_mode parameter to \'separate\'. If not, be sure that the AIF_label_value parameter matches the AIF label value in your ROI. Skipping this volume...'
				AIF_label_image = []

	elif AIF_mode == 'population':
		AIF_label_image = []

	if T1_map_file != []:
		T1_image = nifti_util.nifti_2_numpy(T1_map_file)
	elif T1_map_suffix != []:
		split_path = str.split(filepath, '.nii')
		if os.path.isfile(split_path[0] + T1_map_suffix + '.nii' + split_path[1]):
			T1_image = nifti_util.nifti_2_numpy(split_path[0] + T1_map_suffix + '.nii' + split_path[1])
		else:
			T1_image = []
			print 'No T1 map found at provided T1 map file suffix. Continuing without...'		
	else:
		T1_image = []

	return AIF_label_image, label_image, T1_image

def generate_contrast_agent_concentration(data_numpy, T1_tissue, TR, flip_angle_degrees, injection_start_time_seconds, relaxivity, time_interval_seconds, hematocrit, T1_blood=0, T1_map = []):

	flip_angle_radians = flip_angle_degrees*np.pi/180

	if T1_map != []:
		R1_pre = 1 / T1_map
		R1_pre = np.reshape(R1_pre.shape + (1,))
	elif T1_blood == 0:
		R1_pre = 1 / T1_tissue
	else:
		R1_pre = 1 / T1_blood

	a = np.exp(-1 * TR * R1_pre)
	relative_term = (1-a) / (1-a*np.cos(flip_angle_radians))

	# Note-to-self: Find a way for this to work regardless of array dimension.

	if len(data_numpy.shape) == 1:
		baseline = np.mean(data_numpy[0:int(np.round(injection_start_time_seconds/time_interval_seconds))])
		baseline = np.tile(baseline, data_numpy.shape[-1])
	if len(data_numpy.shape) == 2:
		baseline = np.mean(data_numpy[:,0:int(np.round(injection_start_time_seconds/time_interval_seconds))], axis=1)
		baseline = np.tile(np.reshape(baseline, (baseline.shape[0], 1)), (1,data_numpy.shape[-1]))
	if len(data_numpy.shape) == 3:
		baseline = np.mean(data_numpy[:,:,0:int(np.round(injection_start_time_seconds/time_interval_seconds))], axis=2)
		baseline = np.tile(np.reshape(baseline, (baseline.shape[0],baseline.shape[1], 1)), (1,1,data_numpy.shape[-1]))
	if len(data_numpy.shape) == 4:
		baseline = np.mean(data_numpy[:,:,:,0:int(np.round(injection_start_time_seconds/time_interval_seconds))], axis=3)
		baseline = np.tile(np.reshape(baseline, (baseline.shape[0],baseline.shape[1],baseline.shape[2], 1)), (1,1,1,data_numpy.shape[-1]))
	
	output_numpy = np.copy(data_numpy)

	output_numpy = output_numpy / baseline
	output_numpy = output_numpy * relative_term

	output_numpy = (output_numpy - 1) / (a * (output_numpy * np.cos(flip_angle_radians) - 1))
	output_numpy = -1 * (1 / (relaxivity * TR)) * np.log(output_numpy)

	output_numpy = np.nan_to_num(output_numpy)

	if T1_blood == 0:
		return output_numpy
	else:
		output_numpy = output_numpy / (1-hematocrit)
		return output_numpy

def tensor_optimize(batched_image_part, contrast_image_numpy, contrast_AIF_numpy, time_interval_seconds, mask_value=0, mask_threshold=0, intial_fitting_function_parameters=[1,1]):

	time_series = np.arange(0, contrast_AIF_numpy.size) / (60 / time_interval_seconds)

	image_batch_size = batched_image_part.shape[0]
	identity_batch = np.ones((image_batch_size, image_batch_size), dtype=float)

	tensor_AIF = tf.Variable(tf.constant(contrast_AIF_numpy, dtype=tf.float32), trainable=False)
	tensor_AIF_current = tf.placeholder("float", None)
	tensor_observed_current = tf.placeholder("float",None)
	tensor_AIF_previous = tf.Variable(tf.ones(shape=(image_batch_size, image_batch_size)), dtype=tf.float32, trainable=False)
	tensor_current_estimate = tf.Variable(tf.zeros(shape=(image_batch_size, image_batch_size)), trainable=False)	
	tensor_previous_estimate = tf.Variable(tf.zeros(shape=(image_batch_size, image_batch_size)), trainable=False)
	tensor_estimated = tf.Variable(tf.ones_like(batched_image_part,dtype=tf.float32))

	print 'sequence feeds made'
	tensor_observed = tf.constant(batched_image_part, dtype=tf.float32)
	tensor_time_interval = tf.constant(time_series[1], dtype=tf.float32)
	tensor_zero = tf.constant([0], dtype=tf.float32)
	tensor_indice = tf.placeholder(tf.int32, None)

	print 'constants made'
	tensor_ktrans = tf.Variable(tf.ones(shape=(image_batch_size, image_batch_size)),name='ktrans')
	tensor_ve = tf.Variable(tf.ones(shape=(image_batch_size, image_batch_size)), name='ve')

	print 'variables made'
	tensor_kep = tf.truediv(tensor_ktrans, tensor_ve)
	tensor_log_e = tf.mul(tf.constant(-1, dtype=tf.float32), tf.mul(tensor_kep, tensor_time_interval))
	tensor_capital_e = tf.exp(tensor_log_e)
	tensor_term_A = tf.mul(tensor_AIF_current, tf.sub(tensor_capital_e, tf.sub(tensor_log_e, 1)))
	tensor_term_B = tf.mul(tensor_AIF_previous, tf.sub(tensor_capital_e, tf.sub(tf.mul(tensor_log_e, tensor_capital_e), 1)))
	tensor_integral = tf.truediv(tf.sub(tensor_term_B, tensor_term_A), tf.pow(tensor_log_e, 2))
	tensor_concentration = tf.add(tf.mul(tensor_previous_estimate, tensor_capital_e), tf.mul(tensor_ktrans, tf.mul(tensor_time_interval, tensor_integral)))

	print 'operations made'
	tensor_update_AIF = tf.assign(tensor_AIF_previous, tensor_AIF_current)
	tensor_update_current_estimate = tf.assign(tensor_current_estimate, tensor_concentration)
	tensor_update_previous_estimate = tf.assign(tensor_previous_estimate, tensor_current_estimate)
	tensor_update_total_estimate = tf.assign(tensor_estimated[:,:,tensor_indice], tensor_concentration)

	print 'updaters made'
	tensor_compute_cost = tf.reduce_mean(tf.pow(tf.sub(tensor_observed, tensor_current_estimate),2))

	tensor_loss = tf.add(tensor_zero, tensor_zero)
	for i in xrange(0, np.shape(contrast_AIF_numpy)[-1]):
		tensor_loss = tf.add(tensor_loss, tf.reduce_mean(tf.square(tf.sub(tensor_observed[:,:,i], ))))

	tensor_loss = tf.add(tensor_zero, tensor_zero)
	# tensor_loss = tf.Variable([0], dtype=tf.float32, trainable=False)
	tensor_update_loss = tf.assign(tensor_loss, tf.add(tensor_compute_cost, tensor_loss))
	# tensor_increment_loss = tf.add(tensor_loss, )
	tensor_reset_loss = tf.assign(tensor_loss, tensor_zero)
	optimizer = tf.train.FtrlOptimizer(.7)
	train = optimizer.minimize(tensor_loss) 

	print 'optimizer made'
	trainable = tf.trainable_variables()
	init = tf.initialize_all_variables()
	print 'variables initalizes'
	sess = tf.Session()
	print 'session created'
	sess.run(init)

	print tf.trainable_variables()

	for training_runs in xrange(200):
		for i in xrange(0, np.shape(contrast_AIF_numpy)[-1]):
			
			sess.run(tensor_update_current_estimate, feed_dict={tensor_AIF_current: identity_batch*contrast_AIF_numpy[i]})
			sess.run(tensor_update_loss, feed_dict={tensor_observed: batched_image_part[:,:,i]})

			loss = sess.run(tf.add(tf.square(tf.sub(tensor_observed_current, tensor_current_estimate, feed_dict={tensor_observed_current: batched_image_part[:,:,i]}))))

			sess.run(tensor_update_AIF, feed_dict={tensor_AIF_current: identity_batch*contrast_AIF_numpy[i]})
			sess.run(tensor_update_previous_estimate)

		sess.run(train)
		if training_runs % 2 == 0:
			print (training_runs, sess.run(tensor_ktrans), sess.run(tensor_ve), sess.run(tensor_loss))
		sess.run(tensor_reset_loss)


	print 'escape!'

	fd = fg

	# for i in xrange(1, np.size(estimated_concentration)):

	# y_data = np.ones((image_batch_size, image_batch_size, contrast_image_numpy.shape[-1]), dtype=float)

	sess.close()

	## Tensor Estimate Concentration


def simplex_optimize(contrast_image_numpy, contrast_AIF_numpy, time_interval_seconds, mask_value=0, mask_threshold=0, intial_fitting_function_parameters=[1,1]):
	
	image_batch_size = 10
	# print np.ceil(contrast_image_numpy.shape[0]/image_batch_size)
	for row_batch in xrange(int(np.ceil(contrast_image_numpy.shape[0]/image_batch_size))):
		for col_batch in xrange(int(np.ceil(contrast_image_numpy.shape[1]/image_batch_size))):
			batched_image_part = contrast_image_numpy[row_batch:row_batch+image_batch_size, col_batch:col_batch+image_batch_size, :]
			tensor_optimize(batched_image_part, contrast_image_numpy, contrast_AIF_numpy, time_interval_seconds, mask_value=0, mask_threshold=0, intial_fitting_function_parameters=[1,1])

	# inexplicable minute conversion; investigate
	time_series = np.arange(0, contrast_AIF_numpy.size) / (60 / time_interval_seconds)

	fd = gf

	def cost_function(params):

		estimated_concentration = estimate_concentration(params, contrast_AIF_numpy, time_series)

		difference_term = observed_concentration - estimated_concentration
		difference_term = np.power(difference_term, 2)

		return np.sum(difference_term)


	result_grid = np.zeros((6,5,2),dtype=float)
	np.set_printoptions(threshold=np.nan)
	output_image = np.zeros((contrast_image_numpy.shape[0:-1] + (3,)), dtype=float)

	if len(contrast_image_numpy.shape) == 3:
		# for x in xrange(contrast_image_numpy.shape[0]):
			# for y in xrange(contrast_image_numpy.shape[1]):
		for x_idx, x in enumerate(range(0,50,1)):
			for y_idx, y in enumerate(range(10,70,1)):
				# if contrast_image_numpy[x,y,0] == mask_value or contrast_image_numpy[x,y,0] < mask_threshold:
					# continue

				print[x,y]
				observed_concentration = contrast_image_numpy[x,y,:]

				result_params = scipy.optimize.fmin(cost_function, intial_fitting_function_parameters, disp=1, ftol=.0001, xtol=1e-8)

				ktrans = np.exp(result_params[0]) #ktrans
				ve = 1 / (1 + np.exp(-result_params[1])) #ve
				# auc = np.trapz(observed_concentration, dx=time_interval_seconds) / np.trapz(contrast_AIF_numpy, dx=time_interval_seconds)

				print [ktrans, ve]
				# result_grid[y_idx, x_idx, 0] = ktrans
				# result_grid[y_idx, x_idx, 1] = ve

				# time_series = np.arange(0, contrast_AIF_numpy.size) / (60 / time_interval_seconds)
				# estimated_concentration = estimate_concentration(result_params, contrast_AIF_numpy, time_interval_seconds)
				# plt.plot(time_series, estimated_concentration, 'r--', time_series, observed_concentration, 'b--')
				# plt.show()

				output_image[x,y,0] = ktrans
				output_image[x,y,1] = ve
				# output_image[x,y,2] = auc

		output_image[output_image[:,:,0] > .7] = -.01
		output_image[output_image[:,:,1] > .9] = -.01
		output_image[output_image[:,:,1] < 1e-4] = -.01

		validation_image = np.copy(output_image)
		validation_image = np.ma.masked_equal(validation_image, -.01)
		for y in xrange(0,50,10):
			for x in xrange(10,70,10):
				result_grid[(x)/10,(y-1)/10,0] = np.mean((validation_image[x:x+10,y:y+10,0]).flatten())
				result_grid[(x)/10,(y-1)/10,1] = np.mean((validation_image[x:x+10,y:y+10,1]).flatten())

		print result_grid[:,:,0]
		print result_grid[:,:,1]

	elif len(contrast_image_numpy.shape) == 4:
		for x in xrange(contrast_image_numpy.shape[0]):
			for y in xrange(contrast_image_numpy.shape[1]):
				for z in xrange(contrast_image_numpy.shape[2]):
					if contrast_image_numpy[x,y,z,0] == mask_value or contrast_image_numpy[x,y,z,0] < mask_threshold:
						continue
					print[x,y,z]
					print contrast_image_numpy[x,y,z,0]
					observed_concentration = contrast_image_numpy[x,y,z,:]

					result_params = scipy.optimize.fmin(cost_function, intial_fitting_function_parameters, disp=1, ftol=.0001, xtol=1e-8)

					ktrans = np.exp(result_params[0]) #ktrans
					ve = 1 / (1 + np.exp(-result_params[1])) #ve
					# auc = np.trapz(observed_concentration, dx=time_interval_seconds) / np.trapz(contrast_AIF_numpy, dx=time_interval_seconds)

					print [ktrans, ve]
					output_image[x,y,z,0] = ktrans
					output_image[x,y,z,1] = ve
					# output_image[x,y,z,2] = auc
		output_image[output_image[:,:,:,0] > .7] = 0
		output_image[output_image[:,:,:,1] > .9] = 0

	else:
		print "Fitting not yet implemented for images with less than 3 dimensions or greater than four dimensions."
		return []


	return output_image

def estimate_concentration(params, contrast_AIF_numpy, time_series):

	estimated_concentration = np.zeros_like(contrast_AIF_numpy, dtype=float)
	if params[0] > 10 or params[1] > 10:
		return estimated_concentration

	ktrans = np.exp(params[0])
	ve = 1 / (1 + np.exp(-params[1]))
	kep = ktrans / ve

	# Notation is very inexact here. Clean it up later.

	for i in xrange(1, np.size(estimated_concentration)):
		log_e = -1 * kep * time_series[1]
		capital_E = np.exp(log_e)
		term_A = contrast_AIF_numpy[i] * (capital_E - log_e - 1)
		term_B = contrast_AIF_numpy[i-1] * (capital_E - (capital_E * log_e) - 1)
		integral_term = (term_A - term_B) / np.power(log_e, 2)
		estimated_concentration[i] = estimated_concentration[i-1]*capital_E + ktrans * time_series[1] * integral_term

	# Quick, error prone convolution method
	# res = np.exp(-1*kep*time_series)
	# estimated_concentration = ktrans * np.convolve(contrast_AIF_numpy, res) * time_series[1]
	# estimated_concentration = estimated_concentration[0:np.size(res)]

	return estimated_concentration

def validate_math():
	return

def validate_params():
	return

def validate_label():
	return

def generate_AIF(scan_time_seconds, injection_start_time_seconds, time_interval_seconds, image_numpy=[], AIF_label_numpy=[], AIF_mode='label_average', dimension=4, AIF_label_value=1):

	# It's an open question how to create labels for 2-D DCE phantoms. For now, I assume that people draw their label at time-point zero.

	if AIF_mode == 'label_average':
		if image_numpy != []:
			if AIF_label_numpy != []:


				# Note-to-self: Find a way for this to work regardless of array dimension.
				# Also find a better way to mask arrays.
				# Also this function is extremely messy. It will fail if an image has zero
				# values. This will probably occur, so it must be fixed soon.

				if dimension == 3:
					AIF_subregion = np.copy(image_numpy)
					label_mask = (AIF_label_numpy[:,:,0] != AIF_label_value).reshape((AIF_label_numpy.shape[0:-1] + (1,)))
					AIF_subregion = np.ma.array(AIF_subregion, mask=np.tile(label_mask, (1,)*(dimension-1) + (AIF_subregion.shape[-1],)))
					AIF_subregion = np.reshape(AIF_subregion, (np.product(AIF_subregion.shape[0:-1]), AIF_subregion.shape[-1]))
					AIF = AIF_subregion.mean(axis=0, dtype=np.float64)
					return AIF
				elif dimension == 4:
					AIF_subregion = np.copy(image_numpy)
					label_mask = (AIF_label_numpy[:,:,:] != AIF_label_value).reshape((AIF_label_numpy.shape[0:-1] + (1,)))
					AIF_subregion = np.ma.array(AIF_subregion, mask=np.tile(label_mask, (1,)*(dimension-1) + (AIF_subregion.shape[-1],)))
					AIF_subregion = np.reshape(AIF_subregion, (np.product(AIF_subregion.shape[0:-1]), AIF_subregion.shape[-1]))
					AIF = AIF_subregion.mean(axis=0, dtype=np.float64)
					return AIF_subregion
				else:
					print 'Error: too many or too few dimensions to calculate AIF currently. Unable to calculate AIF.'
					return []
			else:
				'Error: no AIF label detected. Unable to calculate AIF.'
				return []
		else:
			print 'No image provided to AIF function. Set AIF_mode to \'population\' to use a population AIF. Unable to calculate AIF.'
			return []

	elif AIF_mode == 'population':
		print 'Population AIF mode not yet completed. Unable to calculate AIF.'
	return []



def calc_DCE_properties_batch(folder, labels=True, param_file='', outputs=['ktrans','auc','ve'], T1_tissue=1000, T1_blood=1660, relaxivity=.0039, TR=5, TE=2.1, slice_time_interval=.5, hematocrit=0.4, injection_start_time_seconds=60, flip_angle_degrees=30, AIF_mode='population', AIF_label_suffix='-label-AIF', AIF_label_value=1):

	return

def create_4d_from_3d(filepath):
	nifti_3d = nib.load(filepath)
	numpy_3d = nifti_3d.get_data()
	numpy_4d = np.zeros((numpy_3d.shape[0], numpy_3d.shape[1], 5, numpy_3d.shape[2]), dtype=float)
	numpy_3d = np.reshape(numpy_3d, (numpy_3d.shape[0], numpy_3d.shape[1], 1, numpy_3d.shape[2]))

	numpy_4d = np.tile(numpy_3d, (1,1,5,1))
	t1_map = np.zeros((numpy_3d.shape[0], numpy_3d.shape[1], 5), dtype=float)
	t1_map[0:50,0:70,:] = 1000
	print t1_map
	t1_map[t1_map == 0] = 1440
	print t1_map

	nifti_util.save_numpy_2_nifti(numpy_4d, filepath, 'tofts_4d.nii')
	nifti_util.save_numpy_2_nifti(t1_map, filepath, 'tofts_t1map.nii')

if __name__ == '__main__':
	filepath = 'tofts_v6.nii.gz'
	# filepath = 'tofts_v9_5SNR.nii'
	# filepath = 'dce1_mc_ss.nii.gz'

	np.set_printoptions(suppress=True, precision=4)

	calc_DCE_properties_single(filepath, label_file=[], param_file=[], AIF_label_file=[], AIF_values_file=[], outputs=['ktrans','ve','auc'], T1_tissue=1000, T1_blood=1440, relaxivity=.0045, TR=5, TE=2.1, scan_time_seconds=(11*60), hematocrit=0.45, injection_start_time_seconds=60, flip_angle_degrees=30, AIF_mode='label_average', AIF_label_suffix='-AIF-label', AIF_label_value=1, label_mode='separate', default_population_AIF=False, intial_fitting_function_parameters=[.1,.1])
	# calc_DCE_properties_single(filepath, label_file=[], param_file=[], AIF_label_file=[], AIF_values_file=[], outputs=['ktrans','ve','auc'], T1_tissue=1000, T1_blood=1440, relaxivity=.0045, TR=5, TE=2.1, scan_time_seconds=(6*60), hematocrit=0.45, injection_start_time_seconds=60, flip_angle_degrees=30, AIF_mode='label_average', AIF_label_suffix='-AIF-label', AIF_label_value=1, label_mode='separate', default_population_AIF=False, intial_fitting_function_parameters=[-2,.1], outfile_prefix='tofts_v9_5SNR')
	# calc_DCE_properties_single(filepath, label_file=[], param_file=[], AIF_label_file=[], AIF_values_file=[], outputs=['ktrans','ve','auc'], T1_tissue=1000, T1_blood=1440, relaxivity=.0045, TR=5, TE=2.1, scan_time_seconds=(11*60), hematocrit=0.45, injection_start_time_seconds=60, flip_angle_degrees=30, AIF_mode='label_average', AIF_label_suffix='-AIF-label', AIF_label_value=1, label_mode='separate', default_population_AIF=False, intial_fitting_function_parameters=[1,1])


	# create_4d_from_3d(filepath)