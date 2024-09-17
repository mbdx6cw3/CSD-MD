#!/usr/bin/env python

'''
This module is for running a NN with a training set of data.
'''
from __future__ import print_function  # for tf printing
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
from tensorflow.keras.layers import Input, Dense, Layer
from tensorflow.keras.models import Model

class NuclearChargePairs(Layer):
    def __init__(self, _NC2, n_atoms, **kwargs):
        super(NuclearChargePairs, self).__init__()
        self._NC2 = _NC2
        self.n_atoms = n_atoms

    def call(self, atom_nc):
        a = tf.expand_dims(atom_nc, 2)
        b = tf.expand_dims(atom_nc, 1)
        c = a * b
        tri1 = tf.linalg.band_part(c, -1, 0)  # lower
        tri2 = tf.linalg.band_part(tri1, 0, 0)  # lower
        tri = tri1 - tri2
        nonzero_indices = tf.where(tf.not_equal(tri, tf.zeros_like(tri)))
        nonzero_values = tf.gather_nd(tri, nonzero_indices)
        nc_flat = tf.reshape(nonzero_values,
                             shape=(tf.shape(atom_nc)[0],
                                    self._NC2))  # reshape to _NC2
        return nc_flat


class CoordsToNRF(Layer):
    def __init__(self, max_NRF, _NC2, n_atoms, **kwargs):
        super(CoordsToNRF, self).__init__()
        self.max_NRF = max_NRF
        self._NC2 = _NC2
        self.n_atoms = n_atoms
        self.au2kcalmola = 627.5095 * 0.529177

    def compute_output_shape(self, input_shape):
        batch_size = input_shape[0]
        n_atoms = input_shape[1]
        return (batch_size, n_atoms, self._NC2)

    def call(self, coords_nc):
        coords, atom_nc = coords_nc
        a = tf.expand_dims(coords, 2)
        b = tf.expand_dims(coords, 1)
        diff = a - b
        diff2 = tf.reduce_sum(diff ** 2, axis=-1)  # get sqrd diff
        # flatten diff2 so that _NC2 values are left
        tri = tf.linalg.band_part(diff2, -1, 0)  # lower
        nonzero_indices = tf.where(tf.not_equal(tri, tf.zeros_like(tri)))
        nonzero_values = tf.gather_nd(tri, nonzero_indices)
        diff_flat = tf.reshape(nonzero_values,
                               shape=(tf.shape(tri)[0], -1))  # reshape to _NC2
        r = diff_flat ** 0.5
        recip_r2 = 1 / r ** 2
        _NRF = (((atom_nc * self.au2kcalmola) * recip_r2) /
                self.max_NRF)  # scaled
        _NRF = tf.reshape(_NRF,
                          shape=(
                          tf.shape(coords)[0], self._NC2))  # reshape to _NC2
        return _NRF


class Energy(Layer):
    def __init__(self, prescale, **kwargs):
        super(Energy, self).__init__()
        self.prescale = prescale

    def compute_output_shape(self, input_shape):
        batch_size = input_shape[0]
        return (batch_size, 1)

    def call(self, E_scaled):
        E = ((E_scaled - self.prescale[2]) /
             (self.prescale[3] - self.prescale[2]) *
             (self.prescale[1] - self.prescale[0]) + self.prescale[0])
        return E


class Eij(Layer):
    def __init__(self, _NC2, max_Eij, **kwargs):
        super(Eij, self).__init__()
        self._NC2 = _NC2
        self.max_Eij = max_Eij

    def compute_output_shape(self, input_shape):
        batch_size = input_shape[0]
        return (batch_size, self._NC2)

    def call(self, decomp_scaled):
        decomp_scaled = tf.reshape(decomp_scaled,
                                   shape=(tf.shape(decomp_scaled)[0], -1))
        decomp = (decomp_scaled - 0.5) * (2 * self.max_Eij)
        decomp = tf.reshape(decomp,
                            shape=(tf.shape(decomp_scaled)[0],
                                   self._NC2))  # reshape to _NC2
        return decomp


class ERecomposition(Layer):
    def __init__(self, n_atoms, _NC2, **kwargs):
        super(ERecomposition, self).__init__()
        self.n_atoms = n_atoms
        self._NC2 = _NC2

    def compute_output_shape(self, input_shape):
        batch_size = input_shape[0]
        return (batch_size, 1)

    def call(self, coords_decompFE):
        coords, decompFE = coords_decompFE
        decompFE = tf.reshape(decompFE, shape=(tf.shape(decompFE)[0], -1))
        a = tf.expand_dims(coords, 2)
        b = tf.expand_dims(coords, 1)
        diff = a - b
        diff2 = tf.reduce_sum(diff ** 2, axis=-1)  # get sqrd diff
        # flatten diff2 so that _NC2 values are left
        tri = tf.linalg.band_part(diff2, -1, 0)  # lower
        nonzero_indices = tf.where(tf.not_equal(tri, tf.zeros_like(tri)))
        nonzero_values = tf.gather_nd(tri, nonzero_indices)
        diff_flat = tf.reshape(nonzero_values,
                               shape=(tf.shape(tri)[0], -1))  # reshape to _NC2
        r_flat = diff_flat ** 0.5
        recip_r_flat = 1 / r_flat
        norm_recip_r = tf.reduce_sum(recip_r_flat ** 2, axis=1,
                                     keepdims=True) ** 0.5
        eij_E = recip_r_flat / norm_recip_r
        recompE = tf.einsum('bi, bi -> b', eij_E, decompFE)
        recompE = tf.reshape(recompE, shape=(tf.shape(coords)[0], 1))
        return recompE


class Force(Layer):
    def __init__(self, n_atoms, _NC2, **kwargs):
        super(Force, self).__init__()
        self.n_atoms = n_atoms
        self._NC2 = _NC2

    def compute_output_shape(self, input_shape):
        batch_size = input_shape[0]
        return (batch_size, self.n_atoms, 3)

    def call(self, E_coords):
        E, coords = E_coords
        gradients = tf.gradients(E, coords, unconnected_gradients='zero')
        return gradients[0] * -1


class Network(object):
    '''
    '''

    def __init__(self):
        self.model = None

    def build(self, n_atoms, ann_params, prescale):
        '''Input coordinates and z_types into model to get NRFS which then
        are used to predict decompFE, which are then recomposed to give
        Cart Fs and molecular E, both of which could be used in the loss
        function, could weight the E or Fs as required.
        '''

        # set variables
        n_pairs = int(n_atoms * (n_atoms - 1) / 2)
        activations = ann_params["activations"]
        n_layers = ann_params["n_layers"]
        n_nodes = ann_params["n_nodes"]
        if ann_params["n_nodes"] == "auto":
            n_nodes = [n_atoms * 30] * n_layers

        # set prescaling factors
        max_NRF = tf.constant(prescale[4], dtype=tf.float32)
        max_matFE = tf.constant(prescale[5], dtype=tf.float32)
        prescale = tf.constant(prescale, dtype=tf.float32)

        coords_layer = Input(shape=(n_atoms, 3), name='coords_layer')
        nuclear_charge_layer = Input(shape=(n_atoms),
                                     name='nuclear_charge_layer')

        nc_pairs_layer = NuclearChargePairs(n_pairs, n_atoms)(
            nuclear_charge_layer)

        # calculate scaled NRFs from coordinates and nuclear charges
        NRF_layer = CoordsToNRF(max_NRF, n_pairs, n_atoms, name='NRF_layer') \
            ([coords_layer, nc_pairs_layer])

        # define input layer as the NRFs
        connected_layer = NRF_layer

        # loop over hidden layers
        for l in range(n_layers):
            net_layer = Dense(units=n_nodes[l], activation=activations,
                              name='net_layerA{}'.format(l))(connected_layer)
            connected_layer = net_layer

        # output layer
        connected_layer = Dense(units=n_pairs, activation="linear",
                                name='net_layerQ')(connected_layer)

        # calculated unscaled q's
        unscale_qFE_layer = Eij(n_pairs, max_matFE, name='unscale_qF_layer') \
            (connected_layer)

        # calculate the scaled energy from the coordinates and unscaled qFE
        E_layer = ERecomposition(n_atoms, n_pairs, name='qFE_layer') \
            ([coords_layer, unscale_qFE_layer])

        # calculate the unscaled energy
        unscaleE_layer = Energy(prescale, name='unscaleE_layer')(E_layer)

        # obtain the forces by taking the gradient of the energy
        dE_dx = Force(n_atoms, n_pairs, name='dE_dx') \
            ([unscaleE_layer, coords_layer])

        # define the input layers and output layers used in the loss function
        model = Model(
            inputs=[coords_layer, nuclear_charge_layer],
            outputs=[dE_dx, unscale_qFE_layer, unscaleE_layer],
        )

        return model

    def read_params(input_file):
        try:
            param_file = open(input_file, "r")
        except FileNotFoundError:
            print("***ERROR: no input file in the current working directory")
            exit()
        params = {}
        for line in param_file:
            if line.startswith("#"):
                continue
            line = line.strip()
            key_word = line.split(" = ")
            if len(key_word) == 1:
                continue
            key_word_list = key_word[1].split(", ")
            if len(key_word_list) == 1:
                params[key_word[0].strip()] = key_word_list[0]
            if len(key_word_list) > 1:
                params[key_word[0].strip()] = key_word_list
        param_file.close()

        # check that all input is valid and convert types
        params["activations"] = str(params["activations"])
        accepted_strings = ["silu", "linear"]
        if params["activations"] not in accepted_strings:
            print("***ERROR: activation function type not accepted")
            exit()
        try:
            params["epochs"] = int(params["epochs"])
        except ValueError:
            print("***ERROR: Invalid number of epochs")
            exit()
        try:
            params["n_layers"] = int(params["n_layers"])
        except ValueError:
            print("***ERROR: Invalid number of layers")
            exit()
        params["nodes"] = str(params["n_nodes"])
        accepted_strings = ["auto"]
        if params["n_nodes"] not in accepted_strings:
            try:
                if params["n_layers"] == 1:
                    params["n_nodes"] = [int(params["n_nodes"])]
                elif params["n_layers"] > 1:
                    params["n_nodes"] = [eval(i) for i in params["n_nodes"]]
            except ValueError:
                print("***ERROR: Invalid number of nodes")
                exit()
        try:
            params["n_data"] = [eval(i) for i in params["n_data"]]
        except ValueError:
            print("***ERROR: Invalid training set size")
        try:
            params["loss_weights"] = [eval(i) for i in params["loss_weights"]]
        except ValueError:
            print("***ERROR: Invalid weights")
            exit()
        try:
            params["init_lr"] = float(params["init_lr"])
        except ValueError:
            print("***ERROR: Invalid initial learning rate")
            exit()
        try:
            params["min_lr"] = float(params["min_lr"])
        except ValueError:
            print("***ERROR: Invalid minimum learning rate")
            exit()
        try:
            params["lr_patience"] = int(params["lr_patience"])
        except ValueError:
            print("***ERROR: Invalid learning rate patience")
            exit()
        try:
            params["lr_factor"] = float(params["lr_factor"])
        except ValueError:
            print("***ERROR: Invalid learning rate factor")
            exit()
        return params

