import sys
import numpy as np
import warnings
import os.path
from num_input import NumericalInput
from single_nucleon_form_factors import *

class WC_3_flavor:
    def __init__(self, coeff_dict, input_dict):
        """ 
        Class for Wilson coefficients in 3 flavor QCD + QED.

        coeff_dict (dictionary): dictionary containing the initial conditions
        of the dimension-five to dimension-seven three-flavor-QCD Wilson coefficients 
        of the form {'C51' : value, 'C52' : value, ...}. By default all values are
        set to zero. The Wilson coefficients should be specified in the MS-bar scheme 
        at 2 GeV.

        The possible names are (following the notation in the draft):
        Dirac fermion:  'C51', 'C52', 
                        'C61u', 'C61d', 'C61s', 'C62u', 'C62d', 'C62s',
                        'C63u', 'C63d', 'C63s', 'C64u', 'C64d', 'C64s',
                        'C65u', 'C65d', 'C65s', 'C66u', 'C66d', 'C66s',
                        'C67u', 'C67d', 'C67s', 'C68u', 'C68d', 'C68s',
                        'C69u', 'C69d', 'C69s', 'C610u', 'C610d', 'C610s',
                        'C71', 'C72', 'C73', 'C74', 'C75', 'C76', 'C78',
                        'C79u', 'C79d', 'C79s', 'C710u', 'C710d', 'C710s',
                        'C711u', 'C711d', 'C711s', 'C712u', 'C712d', 'C712s',
                        'C713u', 'C713d', 'C713s', 'C714u', 'C714d', 'C714s',
                        'C715u', 'C715d', 'C715s', 'C716u', 'C716d', 'C716s'

        input_dict (dictionary): dictionary containing desired input parameters.
        """
        # The dictionary of input parameters
        self.input_dict = input_dict

        # Quark level interaction operators
        self.wc_name_list = ['C51', 'C52', 
                             'C61u', 'C61d', 'C61s', 'C62u', 'C62d', 'C62s',
                             'C63u', 'C63d', 'C63s', 'C64u', 'C64d', 'C64s',
                             'C65u', 'C65d', 'C65s', 'C66u', 'C66d', 'C66s',
                             'C67u', 'C67d', 'C67s', 'C68u', 'C68d', 'C68s',
                             'C69u', 'C69d', 'C69s', 'C610u', 'C610d', 'C610s',
                             'C71', 'C72', 'C73', 'C74', 'C75', 'C76', 'C77', 'C78',
                             'C79u', 'C79d', 'C79s', 'C710u', 'C710d', 'C710s',
                             'C711u', 'C711d', 'C711s', 'C712u', 'C712d', 'C712s',
                             'C713u', 'C713d', 'C713s', 'C714u', 'C714d', 'C714s',
                             'C715u', 'C715d', 'C715s', 'C716u', 'C716d', 'C716s']

        # Initialize coefficient dictionary
        self.coeff_dict = {}

        # Issue a user warning if a key is not defined:
        for wc_name in coeff_dict.keys():
            if wc_name in self.wc_name_list:
                pass
            else:
                warnings.warn('The key ' + wc_name + ' is not a valid key.')

        # Create the dictionary.
        for wc_name in self.wc_name_list:
            if wc_name in coeff_dict.keys():
                self.coeff_dict[wc_name] = coeff_dict[wc_name]
            else:
                self.coeff_dict[wc_name] = 0.

        # Set default mL
        self.mL = (self.ip['mproton'] + self.ip['mneutron']) / 2.

    def match_nucleon_basis(self, q, mL = None):
        """
        Match the 3-flavor parton basis to the nucleon basis from 2208.07945

        q (numpy array): Outgoing 4-momentum of the electron
        mL (float): see 2208.07945

        Returns a dictionary contianing the Wilson coefficients of the nucleon operators
        """
        # Define the input parameters
        mpi = self.ip['mpi0']
        mp = self.ip['mproton']
        mn = self.ip['mneutron']
        mN = (mp+mn)/2

        alpha = 1/self.ip['alowinv']

        # Quark masses at 2GeV
        mu = self.ip['mu_at_2GeV']
        md = self.ip['md_at_2GeV']
        ms = self.ip['ms_at_2GeV']
        mtilde = 1/(1/mu + 1/md + 1/ms)
        
        # Quark charged
        qu = self.ip['qu']
        qd = self.ip['qd']
        qs = self.ip['qs']

        # Lepton masses
        me = self.ip['me']
        mmu = self.ip['mmu']

        if mL == None:
            mL = self.mL

        # Initialize the form factors to leading order in chiral counting
        F1up = F1('u', 'p', self.ip).value_zero_mom()
        F1dp = F1('d', 'p', self.ip).value_zero_mom()
        F1sp = F1('s', 'p', self.ip).value_zero_mom()

        F1un = F1('u', 'n', self.ip).value_zero_mom()
        F1dn = F1('d', 'n', self.ip).value_zero_mom()
        F1sn = F1('s', 'n', self.ip).value_zero_mom()

        F1spslope = F1('s', 'p', self.ip).first_deriv_zero_mom()
        F1snslope = F1('s', 'n', self.ip).first_deriv_zero_mom()

        F2up = F2('u', 'p', self.ip).value_zero_mom()
        F2dp = F2('d', 'p', self.ip).value_zero_mom()
        F2sp = F2('s', 'p', self.ip).value_zero_mom()

        F2un = F2('u', 'n', self.ip).value_zero_mom()
        F2dn = F2('d', 'n', self.ip).value_zero_mom()
        F2sn = F2('s', 'n', self.ip).value_zero_mom()

        FAup = FA('u', 'p', self.ip).value_zero_mom()
        FAdp = FA('d', 'p', self.ip).value_zero_mom()
        FAsp = FA('s', 'p', self.ip).value_zero_mom()

        FAun = FA('u', 'n', self.ip).value_zero_mom()
        FAdn = FA('d', 'n', self.ip).value_zero_mom()
        FAsn = FA('s', 'n', self.ip).value_zero_mom()

        FSup = FS('u', 'p', self.ip).value_zero_mom()
        FSdp = FS('d', 'p', self.ip).value_zero_mom()
        FSsp = FS('s', 'p', self.ip).value_zero_mom()

        FSun = FS('u', 'n', self.ip).value_zero_mom()
        FSdn = FS('d', 'n', self.ip).value_zero_mom()
        FSsn = FS('s', 'n', self.ip).value_zero_mom()

        FPup_pion = FP('u', 'p', self.ip).value_pion_pole()
        FPdp_pion = FP('d', 'p', self.ip).value_pion_pole()
        FPsp_pion = FP('s', 'p', self.ip).value_pion_pole()

        FPun_pion = FP('u', 'n', self.ip).value_pion_pole()
        FPdn_pion = FP('d', 'n', self.ip).value_pion_pole()
        FPsn_pion = FP('s', 'n', self.ip).value_pion_pole()

        FPup_eta = FP('u', 'p', self.ip).value_eta_pole()
        FPdp_eta = FP('d', 'p', self.ip).value_eta_pole()
        FPsp_eta = FP('s', 'p', self.ip).value_eta_pole()

        FPun_eta = FP('u', 'n', self.ip).value_eta_pole()
        FPdn_eta = FP('d', 'n', self.ip).value_eta_pole()
        FPsn_eta = FP('s', 'n', self.ip).value_eta_pole()

        FPup = FPup_pion + FPup_eta
        FPdp = FPdp_pion + FPdp_eta
        FPsp = FPsp_pion + FPsp_eta

        FPun = FPun_pion + FPun_eta
        FPdn = FPdn_pion + FPdn_eta
        FPsn = FPsn_pion + FPsn_eta

        FPpup_pion = FPprimed('u', 'p', self.ip).value_pion_pole()
        FPpdp_pion = FPprimed('d', 'p', self.ip).value_pion_pole()
        FPpsp_pion = FPprimed('s', 'p', self.ip).value_pion_pole()

        FPpun_pion = FPprimed('u', 'n', self.ip).value_pion_pole()
        FPpdn_pion = FPprimed('d', 'n', self.ip).value_pion_pole()
        FPpsn_pion = FPprimed('s', 'n', self.ip).value_pion_pole()

        FPpup_eta = FPprimed('u', 'p', self.ip).value_eta_pole()
        FPpdp_eta = FPprimed('d', 'p', self.ip).value_eta_pole()
        FPpsp_eta = FPprimed('s', 'p', self.ip).value_eta_pole()

        FPpun_eta = FPprimed('u', 'n', self.ip).value_eta_pole()
        FPpdn_eta = FPprimed('d', 'n', self.ip).value_eta_pole()
        FPpsn_eta = FPprimed('s', 'n', self.ip).value_eta_pole()

        FPpup = FPpup_pion + FPpup_eta
        FPpdp = FPpdp_pion + FPpdp_eta
        FPpsp = FPpsp_pion + FPpsp_eta

        FPpun = FPpun_pion + FPpun_eta
        FPpdn = FPpdn_pion + FPpdn_eta
        FPpsn = FPpsn_pion + FPpsn_eta

        FGp = FG('p', self.ip).value_zero_mom()
        FGn = FG('n', self.ip).value_zero_mom()

        FGtildep = FGtilde('p', self.ip).value_zero_mom()
        FGtilden = FGtilde('n', self.ip).value_zero_mom()

        FT0up = FT0('u', 'p', self.ip).value_zero_mom()
        FT0dp = FT0('d', 'p', self.ip).value_zero_mom()
        FT0sp = FT0('s', 'p', self.ip).value_zero_mom()

        FT0un = FT0('u', 'n', self.ip).value_zero_mom()
        FT0dn = FT0('d', 'n', self.ip).value_zero_mom()
        FT0sn = FT0('s', 'n', self.ip).value_zero_mom()

        FT1up = FT1('u', 'p', self.ip).value_zero_mom()
        FT1dp = FT1('d', 'p', self.ip).value_zero_mom()
        FT1sp = FT1('s', 'p', self.ip).value_zero_mom()

        FT1un = FT1('u', 'n', self.ip).value_zero_mom()
        FT1dn = FT1('d', 'n', self.ip).value_zero_mom()
        FT1sn = FT1('s', 'n', self.ip).value_zero_mom()

        Fgammap = Fgamma('p', self.ip).value_zero_mom()
        Fgamman = Fgamma('n', self.ip).value_zero_mom()

        Fgammatildep = Fgammatilde('p', self.ip).value_zero_mom()
        Fgammatilden = Fgammatilde('n', self.ip).value_zero_mom()

        #############################################################
        # NEEDS TO BE IMPLEMENTED IN single_nucleon_form_factors.py #
        #############################################################
        FT2up = 0.
        FT2dp = 0.
        FT2sp = 0.

        FT2un = 0.
        FT2dn = 0.
        FT2sn = 0.

        FT3up = 0.
        FT3dp = 0.
        FT3sp = 0.

        FT3un = 0.
        FT3dn = 0.
        FT3sn = 0.

        # Matching
        nucleon_coeff_dict = {
            'd1p':  + (1. / mu) * self.coeff_dict['C65u'] * FSup + (1. / md) * self.coeff_dict['C65d'] * FSdp + (1. / ms) * self.coeff_dict['C65s'] * FSsp\
                    + self.coeff_dict['C75'] * Fgammap + self.coeff_dict['C71'] * FGp,

            'd1n':  + (1. / mu) * self.coeff_dict['C65u'] * FSun + (1. / md) * self.coeff_dict['C65d'] * FSdn + (1. / ms) * self.coeff_dict['C65s'] * FSsn\
                    + self.coeff_dict['C75'] * Fgamman + self.coeff_dict['C71'] * FGn,

            'd2p':  + (1. / mu) * self.coeff_dict['C67u'] * FPup + (1. / md) * self.coeff_dict['C67d'] * FPdp + (1. / ms) * self.coeff_dict['C67s'] * FPsp\
                    + self.coeff_dict['C73'] * FGtildep + self.coeff_dict['C77'] * Fgammatildep\
                    + 1j * ((me - mmu) / (2 * mp)) * (self.coeff_dict['C63u'] * FPpup + self.coeff_dict['C63d'] * FPpdp + self.coeff_dict['C63s'] * FPpsp)\
                    + 4 * (me - mmu) * ((1 / mu) * self.coeff_dict['C714u'] * FT3up + (1 / md) * self.coeff_dict['C714d'] * FT3dp + (1 / ms) * self.coeff_dict['C714s'] * FT3sp)\
                    - 1j * ((mmu**2 - me**2) / (2 * mp)) * (self.coeff_dict['C711u'] * FPpup + self.coeff_dict['C711d'] * FPpdp + self.coeff_dict['C711s'] * FPpsp),

            'd2n':  + (1. / mu) * self.coeff_dict['C67u'] * FPun + (1. / md) * self.coeff_dict['C67d'] * FPdn + (1. / ms) * self.coeff_dict['C67s'] * FPsn\
                    + self.coeff_dict['C73'] * FGtilden + self.coeff_dict['C77'] * Fgammatilden\
                    + 1j * ((me - mmu) / (2 * mn)) * (self.coeff_dict['C63u'] * FPpun + self.coeff_dict['C63d'] * FPpdn + self.coeff_dict['C63s'] * FPpsn)\
                    + 4 * (me - mmu) * ((1 / mu) * self.coeff_dict['C714u'] * FT3un + (1 / md) * self.coeff_dict['C714d'] * FT3dn + (1 / ms) * self.coeff_dict['C714s'] * FT3sn)\
                    - 1j * ((mmu**2 - me**2) / (2 * mn)) * (self.coeff_dict['C711u'] * FPpun + self.coeff_dict['C711d'] * FPpdn + self.coeff_dict['C711s'] * FPpsn),

            'd3p':  + (1. / mu) * self.coeff_dict['C66u'] * FSup + (1. / md) * self.coeff_dict['C66d'] * FSdp + (1. / ms) * self.coeff_dict['C66s'] * FSsp\
                    + self.coeff_dict['C72'] * FGp + self.coeff_dict['C76'] * Fgammap,

            'd3n':  + (1. / mu) * self.coeff_dict['C66u'] * FSun + (1. / md) * self.coeff_dict['C66d'] * FSdn + (1. / ms) * self.coeff_dict['C66s'] * FSsn\
                    + self.coeff_dict['C72'] * FGn + self.coeff_dict['C76'] * Fgamman,

            'd4p':  + (1. / mu) * self.coeff_dict['C68u'] * FPup + (1. / md) * self.coeff_dict['C68d'] * FPdp + (1. / ms) * self.coeff_dict['C68s'] * FPsp\
                    + self.coeff_dict['C74'] * FGtildep + self.coeff_dict['C78'] * Fgammatildep\
                    + ((mmu + me) / (2 * mp)) * (self.coeff_dict['C64u'] * FPpup + self.coeff_dict['C64d'] * FPpdp + self.coeff_dict['C64s'] * FPpsp)\
                    - 1j * ((mmu**2 - me**2) / (2 * mp)) * (self.coeff_dict['C712u'] * FPpup + self.coeff_dict['C712d'] * FPpdp + self.coeff_dict['C712s'] * FPpsp)\
                    - 4 * 1j * (me + mmu) * ((1 / mu) * self.coeff_dict['C716u'] * FT3up + (1 / md) * self.coeff_dict['C716d'] * FT3dp + (1 / ms) * self.coeff_dict['C716s'] * FT3sp),

            'd4n':  + (1. / mu) * self.coeff_dict['C68u'] * FPun + (1. / md) * self.coeff_dict['C68d'] * FPdn + (1. / ms) * self.coeff_dict['C68s'] * FPsn\
                    + self.coeff_dict['C74'] * FGtilden + self.coeff_dict['C78'] * Fgammatilden\
                    + ((mmu + me) / (2 * mn)) * (self.coeff_dict['C64u'] * FPpun + self.coeff_dict['C64d'] * FPpdn + self.coeff_dict['C64s'] * FPpsn)\
                    - 1j * ((mmu**2 - me**2) / (2 * mn)) * (self.coeff_dict['C712u'] * FPpun + self.coeff_dict['C712d'] * FPpdn + self.coeff_dict['C712s'] * FPpsn)\
                    - 4 * 1j * (me + mmu) * ((1 / mu) * self.coeff_dict['C716u'] * FT3un + (1 / md) * self.coeff_dict['C716d'] * FT3dn + (1 / ms) * self.coeff_dict['C716s'] * FT3sn),

            'd5p':  + self.coeff_dict['C61u'] * F1up + self.coeff_dict['C61d'] * F1dp + self.coeff_dict['C61s'] * F1sp\
                    + (mmu + me) * (self.coeff_dict['C79u'] * F1up + self.coeff_dict['C79d'] * F1dp + self.coeff_dict['C79s'] * F1sp)\
                    - ((np.dot(q,q)) / (2 * mp)) * ((1 / mu) * self.coeff_dict['C713u'] * (FT1up - 4 * FT2up) + (1 / md) * self.coeff_dict['C713d'] * (FT1dp - 4 * FT2dp) + (1 / ms) * self.coeff_dict['C713s'] * (FT1sp - 4 * FT2sp)),

            'd5n':  + self.coeff_dict['C61u'] * F1un + self.coeff_dict['C61d'] * F1dn + self.coeff_dict['C61s'] * F1sn\
                    + (mmu + me) * (self.coeff_dict['C79u'] * F1un + self.coeff_dict['C79d'] * F1dn + self.coeff_dict['C79s'] * F1sn)\
                    - ((np.dot(q,q)) / (2 * mn)) * ((1 / mu) * self.coeff_dict['C713u'] * (FT1un - 4 * FT2un) + (1 / md) * self.coeff_dict['C713d'] * (FT1dn - 4 * FT2dn) + (1 / ms) * self.coeff_dict['C713s'] * (FT1sn - 4 * FT2sn)),

            'd6p':  - (1./2.) * (self.coeff_dict['C61u'] * F2up + self.coeff_dict['C61d'] * F2dp + self.coeff_dict['C61s'] * F2sp)\
                    - (1./2.) * (mmu + me) * (self.coeff_dict['C79u'] * F2up + self.coeff_dict['C79d'] * F2dp + self.coeff_dict['C79s'] * F2sp),

            'd6n':  - (1./2.) * (self.coeff_dict['C61u'] * F2un + self.coeff_dict['C61d'] * F2dn + self.coeff_dict['C61s'] * F2sn)\
                    - (1./2.) * (mmu + me) * (self.coeff_dict['C79u'] * F2un + self.coeff_dict['C79d'] * F2dn + self.coeff_dict['C79s'] * F2sn),

            'd7p':  + self.coeff_dict['C63u'] * FAup + self.coeff_dict['C63d'] * FAdp + self.coeff_dict['C63s'] * FAsp\
                    + (mmu + me) * (self.coeff_dict['C711u'] * FAup + self.coeff_dict['C711d'] * FAdp + self.coeff_dict['C711s'] * FAsp)\
                    # Missing form factor, check C714
                    + 2 * 1j * ((np.dot(q,q)) / mp) * ((1 / mu) * self.coeff_dict['C714u'] + (1 / md) * self.coeff_dict['C714d'] + (1 / ms) * self.coeff_dict['C714s']),

            'd7n':  + self.coeff_dict['C63u'] * FAun + self.coeff_dict['C63d'] * FAdn + self.coeff_dict['C63s'] * FAsn
                    + (mmu + me) * (self.coeff_dict['C711u'] * FAun + self.coeff_dict['C711d'] * FAdn + self.coeff_dict['C711s'] * FAsn)\
                    # Missing form factor, check C714
                    + 2 * 1j * ((np.dot(q,q)) / mn) * ((1 / mu) * self.coeff_dict['C714u'] + (1 / md) * self.coeff_dict['C714d'] + (1 / ms) * self.coeff_dict['C714s']),

            'd8p':  0.,
            'd8n':  0.,

            'd9p':  - (alpha / np.pi) * self.coeff_dict['C51'] * (mL / np.dot(q,q)) * (qu * F1up + qd * F1dp + qs * F1sp)\
                    - mL * (self.coeff_dict['C79u'] * F1up + self.coeff_dict['C79d'] * F1dp + self.coeff_dict['C79s'] * F1sp),

            'd9n':  - (alpha / np.pi) * self.coeff_dict['C51'] * (mL / np.dot(q,q)) * (qu * F1un + qd * F1dn + qs * F1sn)\
                    - mL * (self.coeff_dict['C79u'] * F1un + self.coeff_dict['C79d'] * F1dn + self.coeff_dict['C79s'] * F1sn),

            'd10p': + (alpha / (2 * np.pi)) * self.coeff_dict['C51'] * (mL / np.dot(q,q)) * (qu * F2up + qd * F2dp + qs * F2sp)\
                    + (mL / 2.) * (self.coeff_dict['C79u'] * F2up + self.coeff_dict['C79d'] * F2dp + self.coeff_dict['C79s'] * F2sp),

            'd10n': + (alpha / (2 * np.pi)) * self.coeff_dict['C51'] * (mL / np.dot(q,q)) * (qu * F2un + qd * F2dn + qs * F2sn)\
                    + (mL / 2.) * (self.coeff_dict['C79u'] * F2un + self.coeff_dict['C79d'] * F2dn + self.coeff_dict['C79s'] * F2sn),

            'd11p': - mL * (self.coeff_dict['C711u'] * FAup + self.coeff_dict['C711d'] * FAdp + self.coeff_dict['C711s'] * FAsp),

            'd11n': - mL * (self.coeff_dict['C711u'] * FAun + self.coeff_dict['C711d'] * FAdn + self.coeff_dict['C711s'] * FAsn),

            'd12p': 0.,
            'd12n': 0.,

            'd13p': + self.coeff_dict['C62u'] * F1up + self.coeff_dict['C62d'] * F1dp + self.coeff_dict['C62s'] * F1sp\
                    - 1j * (mmu - me) * (self.coeff_dict['C710u'] * F1up + self.coeff_dict['C710d'] * F1dp + self.coeff_dict['C710s'] * F1sp)\
                    - ((np.dot(q,q)) / (2 * mp)) * ((1 / mu) * self.coeff_dict['C715u'] * (FT1up - 4 * FT2up) + (1 / md) * self.coeff_dict['C715d'] * (FT1dp - 4 * FT2dp) + (1 / ms) * self.coeff_dict['C715s'] * (FT1sp - 4 * FT2sp)),

            'd13n': + self.coeff_dict['C62u'] * F1un + self.coeff_dict['C62d'] * F1dn + self.coeff_dict['C62s'] * F1sn\
                    - 1j * (mmu - me) * (self.coeff_dict['C710u'] * F1un + self.coeff_dict['C710d'] * F1dn + self.coeff_dict['C710s'] * F1sn)\
                    - ((np.dot(q,q)) / (2 * mn)) * ((1 / mu) * self.coeff_dict['C715u'] * (FT1un - 4 * FT2un) + (1 / md) * self.coeff_dict['C715d'] * (FT1dn - 4 * FT2dn) + (1 / ms) * self.coeff_dict['C715s'] * (FT1sn - 4 * FT2sn)),

            'd14p': - (1./2.) * (self.coeff_dict['C62u'] * F2up + self.coeff_dict['C62d'] * F2dp + self.coeff_dict['C62s'] * F2sp)\
                    + (1j / 2.) * (mmu - me) * (self.coeff_dict['C710u'] * F2up + self.coeff_dict['C710d'] * F2dp + self.coeff_dict['C710s'] * F2sp),

            'd14n': - (1./2.) * (self.coeff_dict['C62u'] * F2un + self.coeff_dict['C62d'] * F2dn + self.coeff_dict['C62s'] * F2sn)\
                    + (1j / 2.) * (mmu - me) * (self.coeff_dict['C710u'] * F2un + self.coeff_dict['C710d'] * F2dn + self.coeff_dict['C710s'] * F2sn),

            'd15p': + self.coeff_dict['C64u'] * FAup + self.coeff_dict['C64d'] * FAdp + self.coeff_dict['C64s'] * FAsp\
                    - 1j * (mmu - me) * (self.coeff_dict['C712u'] * FAup + self.coeff_dict['C712d'] * FAdp + self.coeff_dict['C712s'] * FAsp)\
                    + 2j * ((np.dot(q,q)) / mp) * ((1 / mu) * self.coeff_dict['C716u'] * FT3up + (1 / md) * self.coeff_dict['C716d'] * FT3dp + (1 / ms) * self.coeff_dict['C716s'] * FT3sp),

            'd15n': + self.coeff_dict['C64u'] * FAun + self.coeff_dict['C64d'] * FAdn + self.coeff_dict['C64s'] * FAsn\
                    - 1j * (mmu - me) * (self.coeff_dict['C712u'] * FAun + self.coeff_dict['C712d'] * FAdn + self.coeff_dict['C712s'] * FAsn)\
                    + 2j * ((np.dot(q,q)) / mn) * ((1 / mu) * self.coeff_dict['C716u'] * FT3un + (1 / md) * self.coeff_dict['C716d'] * FT3dn + (1 / ms) * self.coeff_dict['C716s'] * FT3sn),

            'd16p': 0.,
            'd16n': 0.,

            'd17p': + (alpha / np.pi) * self.coeff_dict['C52'] * (mL / np.dot(q,q)) * (qu * F1up + qd * F1dp + qs * F1sp)\
                    + mL * (self.coeff_dict['C710u'] * F1up + self.coeff_dict['C710d'] * F1dp + self.coeff_dict['C710s'] * F1sp),

            'd17n': + (alpha / np.pi) * self.coeff_dict['C52'] * (mL / np.dot(q,q)) * (qu * F1un + qd * F1dn + qs * F1sn)\
                    + mL * (self.coeff_dict['C710u'] * F1un + self.coeff_dict['C710d'] * F1dn + self.coeff_dict['C710s'] * F1sn),

            'd18p': - (alpha / (2 * np.pi)) * self.coeff_dict['C52'] * (mL / np.dot(q,q)) * (qu * F2up + qd * F2dp + qs * F2sp)\
                    - (mL / 2.) * (self.coeff_dict['C710u'] * F2up + self.coeff_dict['C710d'] * F2dp + self.coeff_dict['C710s'] * F2sp),

            'd18n': - (alpha / (2 * np.pi)) * self.coeff_dict['C52'] * (mL / np.dot(q,q)) * (qu * F2un + qd * F2dn + qs * F2sn)\
                    - (mL / 2.) * (self.coeff_dict['C710u'] * F2un + self.coeff_dict['C710d'] * F2dn + self.coeff_dict['C710s'] * F2sn),

            'd19p': + mL * (self.coeff_dict['C712u'] * FAup + self.coeff_dict['C712d'] * FAdp + self.coeff_dict['C712s'] * FAsp),

            'd19n': + mL * (self.coeff_dict['C712u'] * FAun + self.coeff_dict['C712d'] * FAdn + self.coeff_dict['C712s'] * FAsn),

            'd20p': 0,
            'd20n': 0,

            'd21p': + (1 / mu) * self.coeff_dict['C69u'] * FT0up + (1 / md) * self.coeff_dict['C69d'] * FT0dp + (1 / ms) * self.coeff_dict['C69s'] * FT0sp,

            'd21n': + (1 / mu) * self.coeff_dict['C69u'] * FT0un + (1 / md) * self.coeff_dict['C69d'] * FT0dn + (1 / ms) * self.coeff_dict['C69s'] * FT0sn,

            'd22p': - ((1 / mu) * self.coeff_dict['C69u'] * FT1up + (1 / md) * self.coeff_dict['C69d'] * FT1dp + (1 / ms) * self.coeff_dict['C69s'] * FT1sp),

            'd22n': - ((1 / mu) * self.coeff_dict['C69u'] * FT1un + (1 / md) * self.coeff_dict['C69d'] * FT1dn + (1 / ms) * self.coeff_dict['C69s'] * FT1sn),

            'd23p': - ((1 / mu) * self.coeff_dict['C69u'] * FT2up + (1 / md) * self.coeff_dict['C69d'] * FT2dp + (1 / ms) * self.coeff_dict['C69s'] * FT2sp),

            'd23n': - ((1 / mu) * self.coeff_dict['C69u'] * FT2un + (1 / md) * self.coeff_dict['C69d'] * FT2dn + (1 / ms) * self.coeff_dict['C69s'] * FT2sn),

            'd24p': - ((1 / mu) * self.coeff_dict['C69u'] * FT3up + (1 / md) * self.coeff_dict['C69d'] * FT3dp + (1 / ms) * self.coeff_dict['C69s'] * FT3sp),

            'd24n': - ((1 / mu) * self.coeff_dict['C69u'] * FT3un + (1 / md) * self.coeff_dict['C69d'] * FT3dn + (1 / ms) * self.coeff_dict['C69s'] * FT3sn),

            'd25p': + (1 / mu) * self.coeff_dict['C610u'] * FT0up + (1 / md) * self.coeff_dict['C610d'] * FT0dp + (1 / ms) * self.coeff_dict['C610s'] * FT0sp,

            'd25n': + (1 / mu) * self.coeff_dict['C610u'] * FT0un + (1 / md) * self.coeff_dict['C610d'] * FT0dn + (1 / ms) * self.coeff_dict['C610s'] * FT0sn,

            'd26p': - ((1 / mu) * self.coeff_dict['C610u'] * FT1up + (1 / md) * self.coeff_dict['C610d'] * FT1dp + (1 / ms) * self.coeff_dict['C610s'] * FT1sp),

            'd26n': - ((1 / mu) * self.coeff_dict['C610u'] * FT1un + (1 / md) * self.coeff_dict['C610d'] * FT1dn + (1 / ms) * self.coeff_dict['C610s'] * FT1sn),

            'd27p': - ((1 / mu) * self.coeff_dict['C610u'] * FT2up + (1 / md) * self.coeff_dict['C610d'] * FT2dp + (1 / ms) * self.coeff_dict['C610s'] * FT2sp),

            'd27n': - ((1 / mu) * self.coeff_dict['C610u'] * FT2un + (1 / md) * self.coeff_dict['C610d'] * FT2dn + (1 / ms) * self.coeff_dict['C610s'] * FT2sn),

            'd28p': - ((1 / mu) * self.coeff_dict['C610u'] * FT3up + (1 / md) * self.coeff_dict['C610d'] * FT3dp + (1 / ms) * self.coeff_dict['C610s'] * FT3sp),

            'd28n': - ((1 / mu) * self.coeff_dict['C610u'] * FT3un + (1 / md) * self.coeff_dict['C610d'] * FT3dn + (1 / ms) * self.coeff_dict['C610s'] * FT3sn),

            'd29p': - mL * ((1 / mu) * self.coeff_dict['C713u'] * (FT0up - (np.dot(q,q) / mp**2) * FT2up) + (1 / md) * self.coeff_dict['C713d'] * (FT0dp - (np.dot(q,q) / mp**2) * FT2dp) + (1 / ms) * self.coeff_dict['C713s'] * (FT0sp - (np.dot(q,q) / mp**2) * FT2sp)),

            'd29n': - mL * ((1 / mu) * self.coeff_dict['C713u'] * (FT0un - (np.dot(q,q) / mn**2) * FT2un) + (1 / md) * self.coeff_dict['C713d'] * (FT0dn - (np.dot(q,q) / mn**2) * FT2dn) + (1 / ms) * self.coeff_dict['C713s'] * (FT0sn - (np.dot(q,q) / mn**2) * FT2sn)),

            'd30p': - mL * ((1 / mu) * self.coeff_dict['C715u'] * (FT0up - (np.dot(q,q) / mp**2) * FT2up) + (1 / md) * self.coeff_dict['C715d'] * (FT0dp - (np.dot(q,q) / mp**2) * FT2dp) + (1 / ms) * self.coeff_dict['C715s'] * (FT0sp - (np.dot(q,q) / mp**2) * FT2sp)),

            'd30n': - mL * ((1 / mu) * self.coeff_dict['C715u'] * (FT0un - (np.dot(q,q) / mn**2) * FT2un) + (1 / md) * self.coeff_dict['C715d'] * (FT0dn - (np.dot(q,q) / mn**2) * FT2dn) + (1 / ms) * self.coeff_dict['C715s'] * (FT0sn - (np.dot(q,q) / mn**2) * FT2sn)),

            'd31p': - (1j / 4.) * mL * ((1 / mu) * self.coeff_dict['C716u'] * FT0up + (1 / md) * self.coeff_dict['C716d'] * FT0dp + (1 / ms) * self.coeff_dict['C716s'] * FT0sp),

            'd31n': - (1j / 4.) * mL * ((1 / mu) * self.coeff_dict['C716u'] * FT0un + (1 / md) * self.coeff_dict['C716d'] * FT0dn + (1 / ms) * self.coeff_dict['C716s'] * FT0sn), 

            'd32p': - (1j / 4.) * mL * ((1 / mu) * self.coeff_dict['C714u'] * FT0up + (1 / md) * self.coeff_dict['C714d'] * FT0dp + (1 / ms) * self.coeff_dict['C714s'] * FT0sp),

            'd32n': - (1j / 4.) * mL * ((1 / mu) * self.coeff_dict['C714u'] * FT0un + (1 / md) * self.coeff_dict['C714d'] * FT0dn + (1 / ms) * self.coeff_dict['C714s'] * FT0sn),
            }

        return nucleon_coeff_dict

    def get_isoscalar_isovector(self, q, mL):
        """
        Convert the Wilson coefficients into isoscalar and isovector components

        Returns dictionary of isoscalar and isovector Wilson coefficients
        """
        # Specify the total number of operators
        n_op = 20
        #  Generate nucleon operator basis
        wcN = self.match_nucleon_basis(q, mL)
        # Initialize a new dictionary for isoscalar and isovector pieces
        nucleon_coeff_dict_is_iv = {}
        for i in range(1, n_op + 1):
            key_p = 'd' + str(i) + 'p'
            key_n = 'd' + str(i) + 'n'
            nucleon_coeff_dict_is_iv['d'+str(i)+'is'] = (1./2.) * (wcN[key_p] + wcN[key_n])
            nucleon_coeff_dict_is_iv['d'+str(i)+'iv'] = (1./2.) * (wcN[key_p] - wcN[key_n])
            
        return nucleon_coeff_dict_is_iv

    def write_yaml(self, q, element, isotope, interaction, oscb, isochar, plots, mL = None, compute_rate = False, path = './wcN.yaml'):
        """
        Output .yaml file for input into rate computation. See 2208.07945 for more details.

        element (string): target element -- can use UI element spelling Sulfur or it's symbol S (also Britsh sp)
        isotope (int):    target element isotope -- an input of 0 averages over all isotopes
        interaction (string): type of interaction model -- bw is short for "Brown-Wildenthal"
        oscb:(int): Oscillator length scale (in fm). Omit or set to 0 to calc from average A (Abar) over isotopes
                    using:  b=Sqrt[41.467/(45*Abar^(-1/3) -25*Abar^(-2/3))]. Automatically fixed for some interactions.
        isochar (string or list): [1, 0]/isoscalar: isoscalar (1)
                                  [0, 1]/isovector: isovector (tau_3)
                                  [0.5, 0.5]/proton proton
                                  [0.5,-0.5]/neutron neutron
        plots (list, string): vcrm: Vector Charge (or S.I.) Response (M)
                              alsr: Axial Longitudinal Spin Response (Sigma'') 
                              alsr: Axial Transverse Spin Response (Sigma')
                              ssd:  Standard Spin-Dependent (or S.D.) Response
                              vtmr: Vector Transverse Magnetic Response (Delta)
                              vlr:  Vector Longitudinal Response (Phi'')
                              vter: Vector Transverse Electric Response (PhiT')
                              all:  All of the above
                              none: None of the above
        mL (float): 0 is non-relativistic, >0 relativistic

        compute_rates (bool): run computation for decay rate

        path (string): path to where the file should be saved. Defaults to './cwN.yaml'
        """
        if mL == None:
            mL = self.mL
        
        n_op = 20
        # Generate isoscalar and isovector Wilson coefficients
        wcN = self.get_isoscalar_isovector(q, mL)
        # Convert the Wilson coefficient dictionary into a nested list
        ds = []
        for i in range(1, n_op + 1):
            key_is = 'd' + str(i) + 'is'
            key_iv = 'd' + str(i) + 'iv'
            if np.imag(wcN[key_is] != 0.):
                warnings.warn('Key ', key_is, ' has a non-zero imaginary component --- should be real.')
            if np.imag(wcN[key_iv] != 0.):
                warnings.warn('Key ', key_iv, ' has a non-zero imaginary component --- should be real.')
            ds.append([i, np.real(wcN[key_is]), np.real(wcN[key_iv])])
        # cs is irrelevant for the expressions above
        cs = []
        # bs are also irrelevant for the expressions above
        bs = []
        # Format the Wilson coefficients for input into the yaml file
        cNR_model_yaml = 'Element: ' + element + '\n'\
                         + 'Isotope: ' + str(isotope) + '\n'\
                         + 'Interaction: ' + interaction + '\n'\
                         + 'oscb: ' + str(oscb) + '\n'\
                         + 'isochar: ' + isochar + '\n'\
                         + 'plots: ' + plots + '\n'\
                         + 'mL: ' + str(mL) + '\n'\
                         + 'ds: ' + str(ds) + '\n'\
                         + 'cs: ' + str(cs) + '\n'\
                         + 'bs: ' + str(bs) 

        output_file = str(os.path.expanduser(path))

        with open(output_file,'w') as f:
            f.write(cNR_model_yaml)

        if compute_rate == True:
            os.system("python mu2e_relativistic.py wcN.yaml")
