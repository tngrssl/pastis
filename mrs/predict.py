#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
An attempt to calculate the optimal shim currents when performing localized MR spectroscopy on the spinal cord at a C3 cervical level at 7T using the cervical RF coil.

Knowing the volunteer's height, weight and the VOI X-Y-Z position in the scanner's reference frame, the script will output the 1st and 2nd order shim currents that should give a minimum water resonance linewidth. The calculation is done using a polynomial model that was and still is "trained" on a set of healthy volunteers. The optimal shim currents are written in a text file compatible with the FASTESTMAP windows script installed on the console.

All this is implemented as a "predict_shim" class that stores the knowledge and the multi-variable polynomial coeficients.

2019
@author: Tangi Roussel
"""

import mrs.reco as reco
from termcolor import cprint
import numpy as np
import matplotlib.pylab as plt
import datetime
import pickle

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline


class predict_shim():
    """The predict_shim class is an attempt to build a polynomial model to predict shims currents."""

    def __init__(self):

        self.data_filepaths = ""
        self.coil_nChannels = 8
        self._data_filepaths_list = []

        self._data_shim_values = []
        self._data_patient_heights = []
        self._data_patient_weights = []
        self._data_patient_ages = []
        self._data_patient_sex = []
        self._data_vrefs = []
        self._data_pos = []

        self.data_included_for_shim_prediction = {'Patient height': True, 'Patient weight': True, 'Patient age': False, 'Patient sex': False, 'Reference voltage': True, 'VOI position': True}
        self.data_included_for_ref_voltage_prediction = {'Patient height': True, 'Patient weight': True, 'Patient age': False, 'Patient sex': True, 'VOI position': False}

        # stored models
        self._model_poly_order = 3
        self._shim_model_list = []
        self._vref_model = []

        # stored predictions
        self._predicted_shims = []
        self._predicted_vref = []

    def build_library(self):
        """Read all data and extract parameters (patient name, height, weight, etc.) and stores it."""
        # hello
        cprint(">> mrs.predict.predict_shim.build_library:",
               'green', attrs=['bold'])

        print(" > reading data and extracting metrics...")

        # init
        self._data_shim_values = []
        self._data_patient_heights = []
        self._data_patient_weights = []
        self._data_patient_ages = []
        self._data_patient_sex = []
        self._data_vrefs = []
        self._data_pos = []

        # parsing file paths
        self._data_filepaths_list = reco._parse_string_into_list(self.data_filepaths)

        # for each twix file, read it and extract metrics
        for f in self._data_filepaths_list:
            s = reco.MRSData2(f, self.coil_nChannels, None)
            self._data_shim_values.append(s.shims)
            self._data_patient_heights.append(s.patient_height)
            self._data_patient_weights.append(s.patient_weight)
            self._data_patient_ages.append(s.patient_birthyear)
            self._data_patient_sex.append(s.patient_sex)
            self._data_vrefs.append(s.vref)
            self._data_pos.append(s.to_scanner(0, 0, 0))

        # convert to numpy
        self._data_shim_values = np.array(self._data_shim_values)
        self._data_patient_heights = np.array(self._data_patient_heights)
        self._data_patient_weights = np.array(self._data_patient_weights)
        self._data_patient_ages = np.array(self._data_patient_ages)
        self._data_patient_sex = np.array(self._data_patient_sex)
        self._data_vrefs = np.array(self._data_vrefs)
        self._data_pos = np.array(self._data_pos)

        # convert birthyears in age
        # actually, this could be a serious approximation if this thing works and is still used in several years
        # what impacts is the age on the day of the scan, not now...
        now = datetime.datetime.now()
        self._data_patient_ages = now.year - self._data_patient_ages

    def print_library_stats(self):
        """Print statistics about our library."""
        # hello
        cprint(">> mrs.predict.predict_shim.print_library_stats:", 'green', attrs=['bold'])

        print(" > statistics about current stored library...")
        print(" > Patient height = %.2f +/- %.2f" % (self._data_patient_heights.mean(), self._data_patient_heights.std()))
        print(" > Patient weight = %.2f +/- %.2f" % (self._data_patient_weights.mean(), self._data_patient_weights.std()))
        print(" > Patient age = %.2f +/- %.2f" % (self._data_patient_ages.mean(), self._data_patient_ages.std()))
        print(" > Patient sex = M=%d F=%d" % (np.sum(self._data_patient_sex == 0), np.sum(self._data_patient_sex == 1)))
        print("")
        print(" > Reference voltage = %.2f +/- %.2f" % (self._data_vrefs.mean(), self._data_vrefs.std()))
        print(" > Reference X VOI poisition = %.2f +/- %.2f" % (self._data_pos[:, 0].mean(), self._data_pos[:, 0].std()))
        print(" > Reference Y VOI poisition = %.2f +/- %.2f" % (self._data_pos[:, 1].mean(), self._data_pos[:, 1].std()))
        print(" > Reference Z VOI poisition = %.2f +/- %.2f" % (self._data_pos[:, 2].mean(), self._data_pos[:, 2].std()))

    def learn_from_library(self):
        """Build multi-variable polynomial model for shim and vref prediction."""
        # hello
        cprint(">> mrs.predict.predict_shim.learn_from_library:",
               'green', attrs=['bold'])

        print(" > learning and adjusting polynomial coefficients for shim prediction...")
        print(" > including: ", end="", flush=True)
        if(self.data_included_for_shim_prediction["Patient height"]):
            print("patient height", end="", flush=True)
        if(self.data_included_for_shim_prediction["Patient weight"]):
            print(", patient weight", end="", flush=True)
        if(self.data_included_for_shim_prediction["Patient age"]):
            print(", patient age", end="", flush=True)
        if(self.data_included_for_shim_prediction["Patient sex"]):
            print(", patient sex", end="", flush=True)
        if(self.data_included_for_shim_prediction["Reference voltage"]):
            print(", reference voltage", end="", flush=True)
        if(self.data_included_for_shim_prediction["VOI position"]):
            print(", VOI position", end="", flush=True)
        print("")

        # for each shim
        self._shim_model_list = []
        for iShim in range(8):
            # polynomial fit
            this_model = Pipeline([('poly', PolynomialFeatures(
                degree=self._model_poly_order)), ('linear', LinearRegression(fit_intercept=True))])
            # decide what to teach
            x = []
            if(self.data_included_for_shim_prediction["Patient height"]):
                x.append(self._data_patient_heights)
            if(self.data_included_for_shim_prediction["Patient weight"]):
                x.append(self._data_patient_weights)
            if(self.data_included_for_shim_prediction["Patient age"]):
                x.append(self._data_patient_ages)
            if(self.data_included_for_shim_prediction["Patient sex"]):
                x.append(self._data_patient_sex)
            if(self.data_included_for_shim_prediction["Reference voltage"]):
                x.append(self._data_vrefs)
            if(self.data_included_for_shim_prediction["VOI position"]):
                x.append(self._data_pos[:, 0])
                x.append(self._data_pos[:, 1])
                x.append(self._data_pos[:, 2])

            x = np.array(x).transpose()
            y = self._data_shim_values[:, iShim]
            this_model = this_model.fit(x, y)
            # store the model
            self._shim_model_list.append(this_model)

        print(" > learning and adjusting polynomial coefficients for Vref prediction...")
        # now for vref, decide what to teach
        x = []
        print(" > including: ", end="", flush=True)
        if(self.data_included_for_ref_voltage_prediction["Patient height"]):
            x.append(self._data_patient_heights)
            print("patient height", end="", flush=True)
        if(self.data_included_for_ref_voltage_prediction["Patient weight"]):
            x.append(self._data_patient_weights)
            print(", patient weight", end="", flush=True)
        if(self.data_included_for_ref_voltage_prediction["Patient age"]):
            x.append(self._data_patient_ages)
            print(", patient age", end="", flush=True)
        if(self.data_included_for_ref_voltage_prediction["Patient sex"]):
            x.append(self._data_patient_sex)
            print(", patient sex", end="", flush=True)
        if(self.data_included_for_ref_voltage_prediction["VOI position"]):
            x.append(self._data_pos[:, 0])
            x.append(self._data_pos[:, 1])
            x.append(self._data_pos[:, 2])
            print(", VOI position", end="", flush=True)
        print("")

        x = np.array(x).transpose()
        y = self._data_vrefs
        this_model = Pipeline([('poly', PolynomialFeatures(degree=self._model_poly_order)), ('linear', LinearRegression(fit_intercept=True))])
        self._vref_model = this_model.fit(x, y)

    def display(self):
        """Display library, 2D scatter correlation plots and predicted data."""
        # hello
        cprint(">> mrs.predict.predict_shim.display:", 'green', attrs=['bold'])

        shim_labels = ['X', 'Y', 'Z', 'Z2', 'ZX', 'ZY', 'Z2-Y2', 'XY', 'Fref (Hz)']
        first_order_shim_ranges = [-4000, 4000]
        second_order_shim_ranges = [-13000, 13000]

        plt.close("all")
        fig = plt.figure(1)
        fig.canvas.set_window_title("Shims vs. Patient heights")
        fig = plt.figure(2)
        fig.canvas.set_window_title("Shims vs. Patient weights")
        fig = plt.figure(3)
        fig.canvas.set_window_title("Shim vs. Patient age")
        fig = plt.figure(4)
        fig.canvas.set_window_title("Shim vs. Patient sex")
        fig = plt.figure(5)
        fig.canvas.set_window_title("Shim vs. Reference voltage")
        fig = plt.figure(6)
        fig.canvas.set_window_title("Shim vs. X positions")
        fig = plt.figure(7)
        fig.canvas.set_window_title("Shim vs. Y positions")
        fig = plt.figure(8)
        fig.canvas.set_window_title("Shim vs. Z positions")
        fig = plt.figure(9)
        fig.canvas.set_window_title("Patient height/weight vs. position")
        fig = plt.figure(10)
        fig.canvas.set_window_title("Patient height/weight vs. Reference voltage")

        plt.figure(1)
        for i_shim in range(8):
            plt.subplot(2, 4, i_shim + 1)
            plt.scatter(self._data_shim_values[:, i_shim], self._data_patient_heights)
            plt.xlabel(shim_labels[i_shim])
            plt.ylabel('Patient height (m)')
            plt.grid('on')
            if(i_shim < 3):
                plt.xlim(first_order_shim_ranges)
            else:
                plt.xlim(second_order_shim_ranges)

        plt.figure(2)
        for i_shim in range(8):
            plt.subplot(2, 4, i_shim + 1)
            plt.scatter(
                self._data_shim_values[:, i_shim], self._data_patient_weights)
            plt.xlabel(shim_labels[i_shim])
            plt.ylabel('Patient weight (m)')
            plt.grid('on')
            if(i_shim < 3):
                plt.xlim(first_order_shim_ranges)
            else:
                plt.xlim(second_order_shim_ranges)

        plt.figure(3)
        for i_shim in range(8):
            plt.subplot(2, 4, i_shim + 1)
            plt.scatter(
                self._data_shim_values[:, i_shim], self._data_patient_ages)
            plt.xlabel(shim_labels[i_shim])
            plt.ylabel('Patient age (y)')
            plt.grid('on')
            if(i_shim < 3):
                plt.xlim(first_order_shim_ranges)
            else:
                plt.xlim(second_order_shim_ranges)

        plt.figure(4)
        for i_shim in range(8):
            plt.subplot(2, 4, i_shim + 1)
            plt.scatter(
                self._data_shim_values[:, i_shim], self._data_patient_sex)
            plt.xlabel(shim_labels[i_shim])
            plt.ylabel('Patient sex (F=1, M=2)')
            plt.grid('on')
            if(i_shim < 3):
                plt.xlim(first_order_shim_ranges)
            else:
                plt.xlim(second_order_shim_ranges)

        plt.figure(5)
        for i_shim in range(8):
            plt.subplot(2, 4, i_shim + 1)
            plt.scatter(self._data_shim_values[:, i_shim], self._data_vrefs)
            plt.xlabel(shim_labels[i_shim])
            plt.ylabel('Reference voltage (V)')
            plt.grid('on')
            if(i_shim < 3):
                plt.xlim(first_order_shim_ranges)
            else:
                plt.xlim(second_order_shim_ranges)

        plt.figure(6)
        for i_shim in range(8):
            plt.subplot(2, 4, i_shim + 1)
            plt.scatter(
                self._data_shim_values[:, i_shim], self._data_pos[:, 0])
            plt.xlabel(shim_labels[i_shim])
            plt.ylabel('X position (mm)')
            plt.grid('on')
            if(i_shim < 3):
                plt.xlim(first_order_shim_ranges)
            else:
                plt.xlim(second_order_shim_ranges)

        plt.figure(7)
        for i_shim in range(8):
            plt.subplot(2, 4, i_shim + 1)
            plt.scatter(
                self._data_shim_values[:, i_shim], self._data_pos[:, 1])
            plt.xlabel(shim_labels[i_shim])
            plt.ylabel('Y position (mm)')
            plt.grid('on')
            if(i_shim < 3):
                plt.xlim(first_order_shim_ranges)
            else:
                plt.xlim(second_order_shim_ranges)

        plt.figure(8)
        for i_shim in range(8):
            plt.subplot(2, 4, i_shim + 1)
            plt.scatter(
                self._data_shim_values[:, i_shim], self._data_pos[:, 2])
            plt.xlabel(shim_labels[i_shim])
            plt.ylabel('Z position (mm)')
            plt.grid('on')
            if(i_shim < 3):
                plt.xlim(first_order_shim_ranges)
            else:
                plt.xlim(second_order_shim_ranges)

        plt.figure(9)
        plt.subplot(2, 3, 1)
        plt.scatter(self._data_patient_weights, self._data_pos[:, 0])
        plt.xlabel('Patient weight (kg)')
        plt.ylabel('X position (mm)')
        plt.grid('on')
        plt.subplot(2, 3, 2)
        plt.scatter(self._data_patient_weights, self._data_pos[:, 1])
        plt.xlabel('Patient weight (kg)')
        plt.ylabel('Y position (mm)')
        plt.grid('on')
        plt.subplot(2, 3, 3)
        plt.scatter(self._data_patient_weights, self._data_pos[:, 2])
        plt.xlabel('Patient weight (kg)')
        plt.ylabel('Z position (mm)')
        plt.grid('on')
        plt.subplot(2, 3, 4)
        plt.scatter(self._data_patient_heights, self._data_pos[:, 0])
        plt.xlabel('Patient height (m)')
        plt.ylabel('X position (mm)')
        plt.grid('on')
        plt.subplot(2, 3, 5)
        plt.scatter(self._data_patient_heights, self._data_pos[:, 1])
        plt.xlabel('Patient height (m)')
        plt.ylabel('Y position (mm)')
        plt.grid('on')
        plt.subplot(2, 3, 6)
        plt.scatter(self._data_patient_heights, self._data_pos[:, 2])
        plt.xlabel('Patient height (m)')
        plt.ylabel('Z position (mm)')
        plt.grid('on')

        plt.figure(10)
        plt.subplot(2, 3, 1)
        plt.scatter(self._data_patient_heights, self._data_vrefs)
        plt.xlabel('Patient height (m)')
        plt.ylabel('Reference voltage (V)')
        plt.grid('on')
        plt.subplot(2, 3, 2)
        plt.scatter(self._data_patient_weights, self._data_vrefs)
        plt.xlabel('Patient weight (kg)')
        plt.ylabel('Reference voltage (V)')
        plt.grid('on')
        plt.subplot(2, 3, 3)
        plt.scatter(self._data_patient_ages, self._data_vrefs)
        plt.xlabel('Patient age (y)')
        plt.ylabel('Reference voltage (V)')
        plt.grid('on')
        plt.subplot(2, 3, 4)
        plt.scatter(self._data_patient_sex, self._data_vrefs)
        plt.xlabel('Patient sex (F=1, M=2)')
        plt.ylabel('Reference voltage (V)')
        plt.grid('on')
        plt.subplot(2, 3, 5)
        plt.scatter(self._data_pos[:, 1], self._data_vrefs)
        plt.xlabel('Y position (mm)')
        plt.ylabel('Reference voltage (V)')
        plt.grid('on')
        plt.subplot(2, 3, 6)
        plt.scatter(self._data_pos[:, 2], self._data_vrefs)
        plt.xlabel('Z position (mm)')
        plt.ylabel('Reference voltage (V)')
        plt.grid('on')

        for i in range(1, 11):
            plt.figure(i)
            figManager = plt.get_current_fig_manager()
            figManager.window.showMaximized()
            plt.pause(0.01)
            plt.tight_layout()
            plt.pause(0.01)
            plt.tight_layout()

    def predict(self, twix_file_fullpath=None):
        """
        Predict shim currents and B1 voltage based on input parameters extracted from a MRS TWIX file (see input argument) or from user terminal inputs (if no argument passed).

        Parameters
        ----------
        twix_file_fullpath : string
            Full absolute file path pointing to the stored signal (TWIX file)

        """
        # hello
        cprint(">> mrs.predict.predict_shim.predict:", 'green', attrs=['bold'])

        # need to build input array
        x_for_shim = []
        x_for_vref = []
        if(twix_file_fullpath is None):
            print(" > ok, tell me more about your patient...")

            # ok no input TWIX file, let's ask questions
            if(self.data_included_for_shim_prediction["Patient height"]) or self.data_included_for_ref_voltage_prediction["Patient height"]:
                print("\n")
                print('* What is the height of the patient? [m]')
                this_data_patient_heights = input()

            if(self.data_included_for_shim_prediction["Patient weight"] or self.data_included_for_ref_voltage_prediction["Patient weight"]):
                print("\n")
                print('* What is the weight of the patient? [kg]')
                this_data_patient_weights = input()

            if(self.data_included_for_shim_prediction["Patient age"] or self.data_included_for_ref_voltage_prediction["Patient age"]):
                print("\n")
                print('* What is the age of the patient? [years]')
                this_data_patient_ages = input()

            if(self.data_included_for_shim_prediction["Patient sex"] or self.data_included_for_ref_voltage_prediction["Patient sex"]):
                print("\n")
                print(
                    '* What is the sex of the patient? [Female is 1, Male is 2]')
                this_data_patient_sex = input()

            if(self.data_included_for_shim_prediction["Reference voltage"]):
                print("\n")
                print('* What is the reference voltage of the patient? [V]')
                this_data_vrefs = input()

            if(self.data_included_for_shim_prediction["VOI position"] or self.data_included_for_ref_voltage_prediction["VOI position"]):
                print("\n")
                print(
                    "* What is the position of your VOI in the X direction (scanner's reference frame)? [mm]")
                this_data_pos_x = input()
                print("\n")
                print(
                    "* What is the position of your VOI in the Y direction (scanner's reference frame)? [mm]")
                this_data_pos_y = input()
                print("\n")
                print(
                    "* What is the position of your VOI in the Z direction (scanner's reference frame)? [mm]")
                this_data_pos_z = input()

            print("\n\n")
        else:
            print(" > extracting input parameters from TWIX file header...")
            s = reco.MRSData2(twix_file_fullpath, self.coil_nChannels, None)
            this_data_patient_heights = s.patient_height
            this_data_patient_weights = s.patient_weight
            now = datetime.datetime.now()
            this_data_patient_ages = now.year - s.patient_birthyear
            this_data_patient_sex = s.patient_sex
            this_data_vrefs = s.vref
            this_data_pos = s.to_scanner(0, 0, 0)
            this_data_pos_x = this_data_pos[0]
            this_data_pos_y = this_data_pos[1]
            this_data_pos_z = this_data_pos[2]

        # prepare prediction input arrays
        if(self.data_included_for_shim_prediction["Patient height"]):
            x_for_shim.append(this_data_patient_heights)
        if(self.data_included_for_ref_voltage_prediction["Patient height"]):
            x_for_vref.append(this_data_patient_heights)

        if(self.data_included_for_shim_prediction["Patient weight"]):
            x_for_shim.append(this_data_patient_weights)
        if(self.data_included_for_ref_voltage_prediction["Patient weight"]):
            x_for_vref.append(this_data_patient_weights)

        if(self.data_included_for_shim_prediction["Patient age"]):
            x_for_shim.append(this_data_patient_ages)
        if(self.data_included_for_ref_voltage_prediction["Patient age"]):
            x_for_vref.append(this_data_patient_ages)

        if(self.data_included_for_shim_prediction["Patient sex"]):
            x_for_shim.append(this_data_patient_sex)
        if(self.data_included_for_ref_voltage_prediction["Patient sex"]):
            x_for_vref.append(this_data_patient_sex)

        if(self.data_included_for_shim_prediction["Reference voltage"]):
            x_for_shim.append(this_data_vrefs)

        if(self.data_included_for_shim_prediction["VOI position"]):
            x_for_shim.append(this_data_pos_x)
            x_for_shim.append(this_data_pos_y)
            x_for_shim.append(this_data_pos_z)
        if(self.data_included_for_ref_voltage_prediction["VOI position"]):
            x_for_vref.append(this_data_pos_x)
            x_for_vref.append(this_data_pos_y)
            x_for_vref.append(this_data_pos_z)

        print(" > predicting shim...")
        self._predicted_shims = []
        for m in self._shim_model_list:
            x_for_shim = np.array(x_for_shim).transpose()
            this_shim_value = m.predict(x_for_shim.reshape(1, -1))
            self._predicted_shims.append(this_shim_value)
        self._predicted_shims = np.array([s[0] for s in self._predicted_shims])

        print(" > predicting reference voltage...")
        x_for_vref = np.array(x_for_vref).transpose()
        self._predicted_vref = self._vref_model.predict(
            x_for_vref.reshape(1, -1))
        self._predicted_vref = self._predicted_vref[0]

        # round results
        self._predicted_vref = self._predicted_vref.round()
        self._predicted_shims = self._predicted_shims.round()

        # results
        print("* Predicted shim = " + str(self._predicted_shims))
        print("* Predicted reference voltage = " + str(self._predicted_vref))
        print("\n")

    def save_shims_to_disk(self, shim_file_fullpath):
        """
        Save shims to a FASTESTMAP-compatible txt file.

        Parameters
        ----------
        shim_file_fullpath : string
            Full absolute file path pointing to the shim file
        """
        # hello
        cprint(">> mrs.predict.predict_shim.save_shims_to_disk:",
               'green', attrs=['bold'])

        # now let's start prediction
        print(" > saving predicted shims to [" + shim_file_fullpath + "]!")
        # need to check units used in this file, looks weird
        with open(shim_file_fullpath, 'w') as f:
            for s in self._predicted_shims:
                f.write('%6d ' % s)
            # add some dummy 1H frequency
            f.write('%9d' % 297205442.0)
        print(" > you can use this file with the FASTESTMAP script to update the shim voltages!")
        print(
            " > you will be off-resonance for sure! Do not forget to adjust the frequency!")

    def save_model_to_disk(self, pickle_file_fullpath):
        """
        Build multi-variable polynomial model for shim and vref prediction and predict shim currents and B1 voltage based on input parameters extracted from a MRS TWIX file (see input argument) or from user terminal inputs (if no argument passed).

        Parameters
        ----------
        twix_file_fullpath : string
            Full absolute file path pointing to the stored signal (TWIX file)

        """
        # hello
        cprint(">> mrs.predict.predict_shim.save_model_to_disk:",
               'green', attrs=['bold'])

        # now let's start prediction
        print(" > saving polynomial model to " + pickle_file_fullpath + " for later use!")
        with open(pickle_file_fullpath, 'wb') as f:
            pickle.dump([self], f)
