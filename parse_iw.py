"""
This project aims to parse the IMPACT World+ files from the original Microsoft access database to the different
implementation in the available LCA software (SimaPro, openLCA & brightway2). A developer version of IW+ is also parsed
and is made available in an Excel format. This version is aimed for entities wishing to integrated IW+ in their
tools/databases on their own.

In the end, IMPACT World+ will span over 5 different files, in different formats:
- the original Microsoft access database (referred to as the source version)
- the Excel version of IW+ regrouping characterization factors from the source version as well as additional
extrapolated CFs (referred to as the dev version)
- the IW+ version directly implementable in the SimaPro LCA software, in a .csv format (referred to as the SimaPro version)
- the IW+ version directly implementable in the openLCA LCA software, in a .zip format (referred to as the openLCA version)
- the IW+ version directly implementable in the brightway2 LCA software, in a .bw2package format (referred to as the
brightway2 version)

All equations modeling the AGTP of GHGs come from the work of Thomas Gasser (gasser@iiasa.ac.at) and Yue He (heyue@iiasa.ac.at)

file name: parse_iw.py
author: Maxime Agez
e-mail: maxime.agez@polymtl.ca
date created: 02-04-22
python version= 3.9
"""

import pandas as pd
import numpy as np
import os
import pkg_resources
import json
import country_converter as coco
import scipy.sparse
import brightway2 as bw2
import bw2io
import datetime
from datetime import datetime
import csv
import warnings
import uuid
import shutil
import zipfile
import logging
import sqlite3
import math
import molmass
from scipy.stats import gmean


class Parse:
    def __init__(self, path_access_db, version, new_path_access_db, bw2_project=None):
        """
        :param path_access_db: path to the Microsoft access database (source version)
        :param version: the version of IW+ to parse
        :param version: (optional) the name of a brightway2 project in which the database "biosphere3" is available

        Object instance variables:
        -------------------------
            - master_db : the master dataframe where basic IW CFs are stored (what is used to produce the dev.xlsx file)
            - ei35_iw   : the dataframe where IW CFs linked to ecoinvent v3.5 elementary flows are stored
            - ei36_iw   : the dataframe where IW CFs linked to ecoinvent v3.6 elementary flows are stored
            - ei371_iw  : the dataframe where IW CFs linked to ecoinvent v3.7.1 elementary flows are stored
            - ei38_iw   : the dataframe where IW CFs linked to ecoinvent v3.8 elementary flows are stored
            - ei391_iw  : the dataframe where IW CFs linked to ecoinvent v3.9.1 elementary flows are stored
            - ei310_iw  : the dataframe where IW CFs linked to ecoinvent v3.10 elementary flows are stored
            - iw_sp : the dataframe where IW CFs linked to SimaPro elementary flows are stored
            - olca_iw   : the dataframe where IW CFs linked to openLCA elementary flows are stored
            - exio_iw   : the dataframe where IW CFs linked to EXIOBASE elementary flows are stored

        Object insteance methods:
        -------------------------
            - load_cfs()
            - load_basic_cfs()
            - load_acid_eutro_cfs()
            - load_land_use_cfs()
            - load_particulates_cfs()
            - load_water_scarcity_cfs()
            - load_water_availability_eq_cfs()
            - load_water_availability_hh_cfs()
            - load_water_availability_terr_cfs()
            - load_thermally_polluted_water_cfs()
            - apply_rules()
            - create_not_regio_flows()
            - create_regio_flows_for_not_regio_ic()
            - order_things_around()
            - separate_regio_cfs()
            - link_to_ecoinvent()
            - export_to_bw2()
            - link_to_sp()
            - export_to_sp()
            - produce_files()
            - produce_files_hybrid_ecoinvent()

        """

        # ignoring some warnings
        warnings.filterwarnings(action='ignore', category=FutureWarning)
        warnings.filterwarnings(action='ignore', category=np.VisibleDeprecationWarning)
        warnings.filterwarnings(action='ignore', category=pd.errors.PerformanceWarning)
        warnings.filterwarnings(action='ignore', category=UserWarning)

        # set up logging tool
        self.logger = logging.getLogger('IW_Reborn')
        self.logger.setLevel(logging.INFO)
        self.logger.handlers = []
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.propagate = False

        self.path_access_db = path_access_db
        self.new_path_access_db = new_path_access_db
        self.version = str(version)
        self.bw2_project = bw2_project

        # OUTPUTs
        self.master_db = pd.DataFrame()
        self.master_db_not_regio = pd.DataFrame()
        self.ei35_iw = pd.DataFrame()
        self.ei36_iw = pd.DataFrame()
        self.ei371_iw = pd.DataFrame()
        self.ei38_iw = pd.DataFrame()
        self.ei39_iw = pd.DataFrame()
        self.ei310_iw = pd.DataFrame()
        self.simplified_version_ei35 = pd.DataFrame()
        self.simplified_version_ei36 = pd.DataFrame()
        self.simplified_version_ei371 = pd.DataFrame()
        self.simplified_version_ei38 = pd.DataFrame()
        self.simplified_version_ei39 = pd.DataFrame()
        self.simplified_version_ei310 = pd.DataFrame()
        self.ei35_iw_as_matrix = pd.DataFrame()
        self.ei36_iw_as_matrix = pd.DataFrame()
        self.ei371_iw_as_matrix = pd.DataFrame()
        self.ei38_iw_as_matrix = pd.DataFrame()
        self.ei39_iw_as_matrix = pd.DataFrame()
        self.ei310_iw_as_matrix = pd.DataFrame()
        self.iw_sp = pd.DataFrame()
        self.simplified_version_sp = pd.DataFrame()
        self.simplified_version_olca = pd.DataFrame()
        self.simplified_version_bw = pd.DataFrame()
        self.sp_data = {}
        self.olca_iw = pd.DataFrame()
        self.olca_data = {}
        self.olca_data_custom = {}
        self.exio_iw = pd.DataFrame()

        self.conn = sqlite3.connect(self.new_path_access_db)

# -------------------------------------------- Main methods ------------------------------------------------------------

    def load_cfs(self):
        """
        Load the characterization factors and stored them in master_db.
        :return: updated master_db
        """

        self.logger.info("Loading basic characterization factors...")
        self.load_basic_cfs()
        self.logger.info("Loading climate change characterization factors...")
        self.load_climate_change_cfs()
        self.logger.info("Loading ozone layer depletion characterization factors...")
        self.load_ozone_layer_depletion_cfs()
        self.logger.info("Loading photochemical ozone formation characterization factors...")
        self.load_photochemical_ozone_formation()
        self.logger.info("Loading acidification characterization factors...")
        self.load_freshwater_acidification_cfs()
        self.load_terrestrial_acidification_cfs()
        self.logger.info("Loading eutrophication characterization factors...")
        self.load_marine_eutrophication_cfs()
        self.load_freshwater_eutrophication_cfs()
        self.logger.info("Loading land use characterization factors...")
        self.load_land_use_cfs()
        self.logger.info("Loading particulate matter characterization factors...")
        self.load_particulates_cfs()
        self.logger.info("Loading water scarcity characterization factors...")
        self.load_water_scarcity_cfs()
        self.logger.info("Loading water availability characterization factors...")
        self.load_water_availability_eq_cfs()
        self.load_water_availability_hh_cfs()
        self.load_water_availability_terr_cfs()
        self.logger.info("Loading thermally polluted water characterization factors...")
        self.load_thermally_polluted_water_cfs()
        self.logger.info("Loading physical effect on biota characterization factors...")
        self.load_plastic_cfs()
        self.logger.info("Loading fisheries impact characterization factors...")
        self.load_fisheries_cfs()

        self.logger.info("Applying rules...")
        self.apply_rules()

        self.logger.info("Treating regionalized factors...")
        self.create_not_regio_flows()
        self.create_regio_flows_for_not_regio_ic()
        self.order_things_around()
        self.separate_regio_cfs()

        self.logger.info("Linking to ecoinvent elementary flows...")
        self.link_to_ecoinvent()

        self.logger.info("Linking to SimaPro elementary flows...")
        self.link_to_sp()

        self.logger.info("Linking to openLCA elementary flows...")
        self.link_to_olca()

        self.logger.info("Linking to exiobase environmental extensions...")
        self.link_to_exiobase()

        self.logger.info("Prepare the footprint version...")
        self.get_simplified_versions()

        self.get_total_hh_and_eq()

    def export_to_bw2(self, ei_flows_version=None):
        """
        This method creates a brightway2 method with the IW+ characterization factors.
        :param ei_flows_version: [str] Provide a specific ei version (e.g., 3.6) to be used to determine the elementary flows
                                 to be linked to iw+. Default values = eiv3.8 (in 2022)

        :return:
        """

        self.logger.info("Exporting to brightway2...")

        bw2.projects.set_current(self.bw2_project)

        bio = bw2.Database('biosphere3')

        # extract the uuid codes for biosphere flows
        bw_flows_with_codes = (
            pd.DataFrame(
                [(i.as_dict()['name'], i.as_dict()['categories'][0], i.as_dict()['categories'][1], i.as_dict()['code'])
                 if len(i.as_dict()['categories']) == 2
                 else (i.as_dict()['name'], i.as_dict()['categories'][0], 'unspecified', i.as_dict()['code'])
                 for i in bio],
                columns=['Elem flow name', 'Compartment', 'Sub-compartment', 'code'])
        )

        # merge with CF df to have the uuid codes and the corresponding CFs
        if ei_flows_version == '3.5':
            ei_in_bw = self.ei35_iw.merge(bw_flows_with_codes)
        elif ei_flows_version == '3.6':
            ei_in_bw = self.ei36_iw.merge(bw_flows_with_codes)
        elif ei_flows_version == '3.7.1':
            ei_in_bw = self.ei371_iw.merge(bw_flows_with_codes)
        elif ei_flows_version == '3.8':
            ei_in_bw = self.ei38_iw.merge(bw_flows_with_codes)
        elif ei_flows_version == '3.9.1':
            ei_in_bw = self.ei38_iw.merge(bw_flows_with_codes)
        else:
            ei_in_bw = self.ei39_iw.merge(bw_flows_with_codes)
        ei_in_bw_simple = self.simplified_version_bw.merge(bw_flows_with_codes)

        # create total HH and EQ categories
        ei_in_bw.set_index(['Impact category', 'CF unit', 'code'], inplace=True)
        total_hh = ei_in_bw.loc(axis=0)[:, 'DALY'].copy('deep')
        total_hh = total_hh.groupby('code').agg({'Compartment': 'first',
                                                 'Sub-compartment': 'first',
                                                 'Elem flow name': 'first',
                                                 'CAS number': 'first',
                                                 'CF value': sum,
                                                 'Elem flow unit': 'first',
                                                 'MP or Damage': 'first',
                                                 'Native geographical resolution scale': 'first'})
        total_hh.index = pd.MultiIndex.from_product([['Total human health'], ['DALY'], total_hh.index])
        total_eq = ei_in_bw.loc(axis=0)[:, 'PDF.m2.yr'].copy('deep')
        total_eq = total_eq.groupby('code').agg({'Compartment': 'first',
                                                 'Sub-compartment': 'first',
                                                 'Elem flow name': 'first',
                                                 'CAS number': 'first',
                                                 'CF value': sum,
                                                 'Elem flow unit': 'first',
                                                 'MP or Damage': 'first',
                                                 'Native geographical resolution scale': 'first'})
        total_eq.index = pd.MultiIndex.from_product([['Total ecosystem quality'], ['PDF.m2.yr'], total_eq.index])
        ei_in_bw = pd.concat([ei_in_bw, total_hh, total_eq])
        ei_in_bw.index.names = ['Impact category', 'CF unit', 'code']
        ei_in_bw = ei_in_bw.reset_index()

        ei_in_bw.set_index(['Impact category', 'CF unit'], inplace=True)
        ei_in_bw_simple.set_index(['Impact category', 'CF unit'], inplace=True)
        impact_categories = ei_in_bw.index.drop_duplicates()
        impact_categories_simple = ei_in_bw_simple.index.drop_duplicates()

        # -------------- For complete version of IW+ ----------------
        for ic in impact_categories:

            if ei_in_bw.loc[[ic], 'MP or Damage'].iloc[0] == 'Midpoint':
                mid_end = 'Midpoint'
                # create the name of the method
                name = ('IMPACT World+ ' + mid_end + ' ' + self.version, 'Midpoint', ic[0])
            else:
                mid_end = 'Damage'
                # create the name of the method
                if ic[1] == 'DALY':
                    name = ('IMPACT World+ ' + mid_end + ' ' + self.version, 'Human health', ic[0])
                else:
                    name = ('IMPACT World+ ' + mid_end + ' ' + self.version, 'Ecosystem quality', ic[0])

            # initialize the "Method" method
            new_method = bw2.Method(name)
            # register the new method
            new_method.register()
            # set its unit
            new_method.metadata["unit"] = ic[1]

            df = ei_in_bw.loc[[ic], ['code', 'CF value']].copy()
            df.set_index('code', inplace=True)

            data = []
            for stressor in df.index:
                data.append((('biosphere3', stressor), df.loc[stressor, 'CF value']))
            new_method.write(data)

        # -------------- For simplified version of IW+ ----------------
        for ic in impact_categories_simple:

            name = ('IMPACT World+ Footprint ' + self.version, ic[0])

            # initialize the "Method" method
            new_method = bw2.Method(name)
            # register the new method
            new_method.register()
            # set its unit
            new_method.metadata["unit"] = ic[1]

            df = ei_in_bw_simple.loc[[ic], ['code', 'CF value']].copy()
            df.set_index('code', inplace=True)

            data = []
            for stressor in df.index:
                data.append((('biosphere3', stressor), df.loc[stressor, 'CF value']))
            new_method.write(data)

    def export_to_sp(self):
        """
        This method creates the necessary information for the csv creation in SimaPro.
        :return:
        """

        self.logger.info("Exporting to SimaPro...")

        # csv accepts strings only
        self.iw_sp.loc[:, 'CF value'] = self.iw_sp.loc[:, 'CF value'].astype(str)
        self.simplified_version_sp.loc[:, 'CF value'] = self.simplified_version_sp.loc[:, 'CF value'].astype(str)

        # Metadata
        l = ['SimaPro 9.3.0.3', 'methods', 'Date: ' + datetime.now().strftime("%D"),
             'Time: ' + datetime.now().strftime("%H:%M:%S"),
             'Project: Methods', 'CSV Format version: 8.0.5', 'CSV separator: Semicolon',
             'Decimal separator: .', 'Date separator: -', 'Short date format: yyyy-MM-dd', 'Selection: Selection (1)',
             'Related objects (system descriptions, substances, units, etc.): Yes',
             'Include sub product stages and processes: No', "Open library: 'Methods'"]
        metadata = []
        for i in l:
            s = '{' + i + '}'
            metadata.append([s, '', '', '', '', ''])

        # metadata on the midpoint method
        midpoint_method_metadata = [['Method', '', '', '', '', ''], ['', '', '', '', '', ''],
                                    ['Name', '', '', '', '', ''],
                                    ['IMPACTWorld+ Midpoint ' + self.version, '', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                    [self.version.split('.')[0],self.version.split('.')[1], '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                    ['For more information on IMPACT World+: https://www.impactworldplus.org/en/.\n'
                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes', '', '', '', '', ''], ['', '', '', '', '', ''],
                                    ['Category', '', '', '', '', ''], ['Others', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Damage Assessment', '', '', '', '', ''], ['No', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Normalization', '', '', '', '', ''], ['No', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Weighting', '', '', '', '', ''], ['No', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Addition', '', '', '', '', ''], ['No', '', '', '', '', '']]
        # metadata on the damage method
        damage_method_metadata = [['Method', '', '', '', '', ''], ['', '', '', '', '', ''],
                                  ['Name', '', '', '', '', ''],
                                  ['IMPACTWorld+ Expert ' + self.version, '', '', '', '', ''],
                                  ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                  [self.version.split('.')[0],self.version.split('.')[1], '', '', '', ''],
                                  ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                  ['For more information on IMPACT World+: https://www.impactworldplus.org/en/.\n'
                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes', '', '', '', '', ''], ['', '', '', '', '', ''],
                                  ['Category', '', '', '', '', ''], ['Others', '', '', '', '', ''],
                                  ['', '', '', '', '', ''],
                                  ['Use Damage Assessment', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                  ['', '', '', '', '', ''],
                                  ['Use Normalization', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                  ['', '', '', '', '', ''],
                                  ['Use Weighting', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                  ['', '', '', '', '', ''],
                                  ['Use Addition', '', '', '', '', ''], ['Yes', '', '', '', '', '']]
        # metadata on the combined method
        combined_method_metadata = [['Method', '', '', '', '', ''], ['', '', '', '', '', ''],
                                    ['Name', '', '', '', '', ''],
                                    ['IMPACTWorld+ ' + self.version, '', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                    [self.version.split('.')[0], self.version.split('.')[1], '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                    ['For more information on IMPACT World+: https://www.impactworldplus.org/en/.\n'
                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes', '', '', '', '', ''], ['', '', '', '', '', ''],
                                    ['Category', '', '', '', '', ''], ['Others', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Damage Assessment', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Normalization', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Weighting', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Addition', '', '', '', '', ''], ['Yes', '', '', '', '', '']]
        # metadata on the simplified method
        simplified_method_metadata = [['Method', '', '', '', '', ''], ['', '', '', '', '', ''],
                                    ['Name', '', '', '', '', ''],
                                    ['IMPACTWorld+ Footprint ' + self.version, '', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                    [self.version.split('.')[0],self.version.split('.')[1], '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                    ['For more information on IMPACT World+: https://www.impactworldplus.org/en/.\n'
                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes', '', '', '', '', ''], ['', '', '', '', '', ''],
                                    ['Category', '', '', '', '', ''], ['Others', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Damage Assessment', '', '', '', '', ''], ['No', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Normalization', '', '', '', '', ''], ['No', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Weighting', '', '', '', '', ''], ['No', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Addition', '', '', '', '', ''], ['No', '', '', '', '', '']]

        # data for weighting and normalizing
        df = pd.read_csv(pkg_resources.resource_filename(
            __name__, '/Data/weighting_normalizing/weighting_and_normalization.csv'),
            header=None, delimiter=';').fillna('')
        weighting_info_damage = [[df.loc[i].tolist()[0], df.loc[i].tolist()[1], '', '', '', ''] for i in df.index]
        weighting_info_damage.insert(36, ['Physical effects on biota', '1.00E+00', '', '', '', ''])

        weighting_info_combined = weighting_info_damage.copy()
        weighting_info_combined[11] = ['Ozone layer depletion (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[12] = ['Particulate matter formation (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[13] = ['Photochemical oxidant formation (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[22] = ['Freshwater acidification (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[25] = ['Freshwater eutrophication (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[27] = ['Land occupation, biodiversity (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[28] = ['Land transformation, biodiversity (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[31] = ['Marine eutrophication (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[32] = ['Terrestrial acidification (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[36] = ['Physical effects on biota (damage)', '1.00E+00', '', '', '', '']

        # extracting midpoint CFs
        d_ic_unit = self.iw_sp.loc[self.iw_sp['MP or Damage'] == 'Midpoint',
                              ['Impact category', 'CF unit']].drop_duplicates().set_index('Impact category').iloc[:,
                    0].to_dict()
        midpoint_values = []
        for j in d_ic_unit.keys():
            midpoint_values.append(['', '', '', '', '', ''])
            midpoint_values.append(['Impact category', '', '', '', '', ''])
            midpoint_values.append([j, d_ic_unit[j], '', '', '', ''])
            midpoint_values.append(['', '', '', '', '', ''])
            midpoint_values.append(['Substances', '', '', '', '', ''])
            df = self.iw_sp[self.iw_sp['Impact category'] == j]
            df = df[df['CF unit'] == d_ic_unit[j]]
            df = df[['Compartment', 'Sub-compartment', 'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit']]
            for i in df.index:
                if type(df.loc[i, 'CAS number']) == float:
                    df.loc[i, 'CAS number'] = ''
                midpoint_values.append(df.loc[i].tolist())

        # extracting damage CFs
        d_ic_unit = self.iw_sp.loc[self.iw_sp['MP or Damage'] == 'Damage',
                              ['Impact category', 'CF unit']].drop_duplicates().set_index('Impact category').iloc[:,
                    0].to_dict()
        damage_values = []
        for j in d_ic_unit.keys():
            damage_values.append(['', '', '', '', '', ''])
            damage_values.append(['Impact category', '', '', '', '', ''])
            damage_values.append([j, d_ic_unit[j], '', '', '', ''])
            damage_values.append(['', '', '', '', '', ''])
            damage_values.append(['Substances', '', '', '', '', ''])
            df = self.iw_sp[self.iw_sp['Impact category'] == j]
            df = df[df['CF unit'] == d_ic_unit[j]]
            df = df[['Compartment', 'Sub-compartment', 'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit']]
            for i in df.index:
                if type(df.loc[i, 'CAS number']) == float:
                    df.loc[i, 'CAS number'] = ''
                damage_values.append(df.loc[i].tolist())

        # extracting combined CFs
        ic_unit = self.iw_sp.loc[:, ['Impact category', 'CF unit']].drop_duplicates()
        same_names = ['Freshwater acidification','Freshwater eutrophication','Land occupation, biodiversity',
                      'Land transformation, biodiversity','Marine eutrophication','Ozone layer depletion',
                      'Particulate matter formation','Photochemical oxidant formation','Terrestrial acidification',
                      'Physical effects on biota']
        combined_values = []
        for j in ic_unit.index:
            combined_values.append(['', '', '', '', '', ''])
            combined_values.append(['Impact category', '', '', '', '', ''])
            if ic_unit.loc[j,'Impact category'] in same_names:
                if ic_unit.loc[j,'CF unit'] in ['DALY','PDF.m2.yr']:
                    combined_values.append([ic_unit.loc[j,'Impact category']+' (damage)',
                                            ic_unit.loc[j,'CF unit'], '', '', '', ''])
                else:
                    combined_values.append([ic_unit.loc[j,'Impact category']+' (midpoint)',
                                            ic_unit.loc[j,'CF unit'], '', '', '', ''])
            else:
                combined_values.append([ic_unit.loc[j, 'Impact category'],
                                        ic_unit.loc[j, 'CF unit'], '', '', '', ''])
            combined_values.append(['', '', '', '', '', ''])
            combined_values.append(['Substances', '', '', '', '', ''])
            df = self.iw_sp.loc[[i for i in self.iw_sp.index if (
                    self.iw_sp.loc[i,'Impact category'] == ic_unit.loc[j,'Impact category'] and
                    self.iw_sp.loc[i, 'CF unit'] == ic_unit.loc[j, 'CF unit'])]]
            df = df[['Compartment', 'Sub-compartment', 'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit']]
            for i in df.index:
                if type(df.loc[i, 'CAS number']) == float:
                    df.loc[i, 'CAS number'] = ''
                combined_values.append(df.loc[i].tolist())

        # extracting simplified values
        ic_unit = self.simplified_version_sp.loc[:, ['Impact category', 'CF unit']].drop_duplicates()
        simplified_values = []
        for j in ic_unit.index:
            simplified_values.append(['', '', '', '', '', ''])
            simplified_values.append(['Impact category', '', '', '', '', ''])
            simplified_values.append([ic_unit.loc[j, 'Impact category'],
                                      ic_unit.loc[j, 'CF unit'], '', '', '', ''])
            simplified_values.append(['', '', '', '', '', ''])
            simplified_values.append(['Substances', '', '', '', '', ''])
            df = self.simplified_version_sp.loc[[i for i in self.simplified_version_sp.index if (
                    self.simplified_version_sp.loc[i, 'Impact category'] == ic_unit.loc[j, 'Impact category'] and
                    self.simplified_version_sp.loc[i, 'CF unit'] == ic_unit.loc[j, 'CF unit'])]]
            df = df[['Compartment', 'Sub-compartment', 'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit']]
            for i in df.index:
                simplified_values.append(df.loc[i].tolist())

        # dump everything in an attribute
        self.sp_data = {'metadata': metadata, 'midpoint_method_metadata': midpoint_method_metadata,
                        'damage_method_metadata': damage_method_metadata,
                        'combined_method_metadata': combined_method_metadata,
                        'simplified_method_metadata': simplified_method_metadata,
                        'weighting_info_damage': weighting_info_damage,
                        'weighting_info_combined': weighting_info_combined,
                        'midpoint_values': midpoint_values, 'damage_values': damage_values,
                        'combined_values': combined_values, 'simplified_values': simplified_values}

    def export_to_olca(self):
        """
        This method creates the necessary information for the creation of json files in openLCA.
        :return:
        """

        self.logger.info("Exporting to openLCA...")

        # --------------------- GENERAL METADATA OF IW+ ------------------------
        id_category = str(uuid.uuid4())

        category_metadata = {"@context": "http://greendelta.github.io/olca-schema/context.jsonld",
                             "@type": "Category",
                             "@id": id_category,
                             "name": "IMPACT World+",
                             "version": "0.00.000",
                             "modelType": "IMPACT_METHOD"}

        # -----------------------IW+ VERSION METADATA --------------------------
        category_names = list(set([(i[0], i[1]) for i in self.olca_iw.index]))
        category_names = [(i[0], i[1], str(uuid.uuid4())) for i in category_names]
        category_names = {(i[0], i[1]): i[2] for i in category_names}
        category_names_damage = {k: v for k, v in category_names.items() if k[1] in ['DALY', 'PDF.m2.yr']}
        category_names_midpoint = {k: v for k, v in category_names.items() if k[1] not in ['DALY', 'PDF.m2.yr']}
        category_names_footprint = list(set([(i[0], i[1]) for i in self.simplified_version_olca.index]))
        category_names_footprint = {(category_names_footprint[i][0], category_names_footprint[i][1]): str(uuid.uuid4())
                                    for i in range(len(category_names_footprint))}
        # need to differentiate midpoint from endpoint with the names of the categories
        category_names_combined = {}
        for category in category_names:
            if category[1] in ['DALY', 'PDF.m2.yr']:
                category_names_combined[(category[0] + ' (damage)', category[1])] = category_names[category]
            else:
                category_names_combined[(category[0] + ' (midpoint)', category[1])] = category_names[category]

        id_iw_damage = str(uuid.uuid4())
        id_iw_midpoint = str(uuid.uuid4())
        id_iw_footprint = str(uuid.uuid4())
        id_iw_combined = str(uuid.uuid4())
        norm_weight_id = str(uuid.uuid4())

        metadata_iw_damage = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "ImpactMethod",
            "@id": id_iw_damage,
            "name": "IMPACT World+ Expert v" + self.version,
            "lastChange": "2022-09-15T17:25:43.725-05:00",
            "category": {
                "@type": "Category",
                "@id": id_category,
                "name": "IMPACT World+",
                "categoryType": "ImpactMethod"},
            'impactCategories': [],
            'nwSets': [{
                "@type": "NwSet",
                "@id": norm_weight_id,
                "name": "IMPACT World+ (Stepwise 2006 values)"
            }]
        }

        metadata_iw_midpoint = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "ImpactMethod",
            "@id": id_iw_midpoint,
            "name": "IMPACT World+ Midpoint v" + self.version,
            "lastChange": "2022-09-15T17:25:43.725-05:00",
            "category": {
                "@type": "Category",
                "@id": id_category,
                "name": "IMPACT World+",
                "categoryType": "ImpactMethod"},
            'impactCategories': []
        }

        metadata_iw_footprint = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "ImpactMethod",
            "@id": id_iw_footprint,
            "name": "IMPACT World+ Footprint v" + self.version,
            "lastChange": "2022-09-15T17:25:43.725-05:00",
            "category": {
                "@type": "Category",
                "@id": id_category,
                "name": "CIRAIG methods",
                "categoryType": "ImpactMethod"},
            'impactCategories': []
        }

        metadata_iw_combined = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "ImpactMethod",
            "@id": id_iw_combined,
            "name": "IMPACT World+ Combined v" + self.version,
            "lastChange": "2022-09-15T17:25:43.725-05:00",
            "category": {
                "@type": "Category",
                "@id": id_category,
                "name": "IMPACT World+",
                "categoryType": "ImpactMethod"},
            'impactCategories': [],
            'nwSets': [{
                "@type": "NwSet",
                "@id": norm_weight_id,
                "name": "IMPACT World+ (Stepwise 2006 values)"
            }]
        }

        # ------------------------- IMPACT CATEGORIES METADATA -----------------
        for cat in category_names_damage:
            metadata_iw_damage['impactCategories'].append(
                {"@type": "ImpactCategory",
                 "@id": category_names_damage[cat],
                 "name": cat[0],
                 "refUnit": cat[1]}
            )

        for cat in category_names_midpoint:
            metadata_iw_midpoint['impactCategories'].append(
                {"@type": "ImpactCategory",
                 "@id": category_names_midpoint[cat],
                 "name": cat[0],
                 "refUnit": cat[1]}
            )

        for cat in category_names_footprint:
            metadata_iw_footprint['impactCategories'].append(
                {"@type": "ImpactCategory",
                 "@id": category_names_footprint[cat],
                 "name": cat[0],
                 "refUnit": cat[1]}
            )

        for cat in category_names_combined:
            metadata_iw_combined['impactCategories'].append(
                {"@type": "ImpactCategory",
                 "@id": category_names_combined[cat],
                 "name": cat[0],
                 "refUnit": cat[1]}
            )

        # --------------------------------- UNIT METADATA -------------------------
        # hardcoded, obtained from going inside the json files of an exported oLCA method
        unit_groups = {
            'kg': '93a60a57-a4c8-11da-a746-0800200c9a66',
            'm2*a': '93a60a57-a3c8-20da-a746-0800200c9a66',
            'm2': '93a60a57-a3c8-18da-a746-0800200c9a66',
            'kBq': '93a60a57-a3c8-16da-a746-0800200c9a66',
            'm3': '93a60a57-a3c8-12da-a746-0800200c9a66',
            'MJ': '93a60a57-a3c8-11da-a746-0800200c9a66'
        }
        flow_properties = {
            'm3': '93a60a56-a3c8-22da-a746-0800200c9a66',
            'm2*a': '93a60a56-a3c8-21da-a746-0800200c9a66',
            'm2': '93a60a56-a3c8-19da-a746-0800200c9a66',
            'kBq': '93a60a56-a3c8-17da-a746-0800200c9a66',
            'kg': '93a60a56-a3c8-11da-a746-0800200b9a66',
            'MJ': 'f6811440-ee37-11de-8a39-0800200c9a66'
        }

        # ------------------------- INTEGRATE CHARACTERIZATION FACTORS -----------------
        cf_dict_damage = {}
        for cat in category_names_damage:
            cf_values = {
                "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
                "@type": "ImpactCategory",
                "@id": category_names_damage[cat],
                "name": cat[0],
                "category": metadata_iw_damage["name"],
                "version": "0" + ''.join(self.version.split('.')) + "0.000",
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.olca_iw.loc[cat].copy()

            for flow_id in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[flow_id, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": flow_id,
                        "name": dff.loc[flow_id, 'flow_name'],
                        "category": "Elementary flows/"+eval(dff.loc[flow_id, 'comp'])[0]+'/'+eval(dff.loc[flow_id, 'comp'])[1],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[flow_id, 'unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[flow_id, 'unit']],
                        "name": dff.loc[flow_id, 'unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[flow_id, 'unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[flow_id, 'unit']
                    }
                })

            cf_dict_damage[cat] = cf_values

        cf_dict_midpoint = {}
        for cat in category_names_midpoint:
            cf_values = {
                "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
                "@type": "ImpactCategory",
                "@id": category_names_midpoint[cat],
                "name": cat[0],
                "category": metadata_iw_midpoint["name"],
                "version": "0" + ''.join(self.version.split('.')) + "0.000",
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.olca_iw.loc[cat].copy()

            for flow_id in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[flow_id, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": flow_id,
                        "name": dff.loc[flow_id, 'flow_name'],
                        "category": "Elementary flows/"+eval(dff.loc[flow_id, 'comp'])[0]+'/'+eval(dff.loc[flow_id, 'comp'])[1],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[flow_id, 'unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[flow_id, 'unit']],
                        "name": dff.loc[flow_id, 'unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[flow_id, 'unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[flow_id, 'unit']
                    }
                })

            cf_dict_midpoint[cat] = cf_values

        cf_dict_footprint = {}
        for cat in category_names_footprint:
            cf_values = {
                "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
                "@type": "ImpactCategory",
                "@id": category_names_footprint[cat],
                "name": cat[0],
                "category": metadata_iw_footprint["name"],
                "version": "0" + ''.join(self.version.split('.')) + "0.000",
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.simplified_version_olca.loc[cat].copy()

            for flow_id in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[flow_id, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": flow_id,
                        "name": dff.loc[flow_id, 'flow_name'],
                        "category": "Elementary flows/"+eval(dff.loc[flow_id, 'comp'])[0]+'/'+eval(dff.loc[flow_id, 'comp'])[1],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[flow_id, 'unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[flow_id, 'unit']],
                        "name": dff.loc[flow_id, 'unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[flow_id, 'unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[flow_id, 'unit']
                    }
                })

            cf_dict_footprint[cat] = cf_values

        cf_dict_combined = {}
        for cat in category_names_combined:
            cf_values = {
                "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
                "@type": "ImpactCategory",
                "@id": category_names_combined[cat],
                "name": cat[0],
                "category": metadata_iw_combined["name"],
                "version": "0" + ''.join(self.version.split('.')) + "0.000",
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.olca_iw.loc[(cat[0].split(' (')[0], cat[1])].copy()

            for flow_id in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[flow_id, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": flow_id,
                        "name": dff.loc[flow_id, 'flow_name'],
                        "category": "Elementary flows/"+eval(dff.loc[flow_id, 'comp'])[0]+'/'+eval(dff.loc[flow_id, 'comp'])[1],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[flow_id, 'unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[flow_id, 'unit']],
                        "name": dff.loc[flow_id, 'unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[flow_id, 'unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[flow_id, 'unit']
                    }
                })

            cf_dict_combined[cat] = cf_values

        # ------------------------- NORMALIZATION AND WEIGHTING -------------------
        norm = {}
        for i in category_names:
            if i[1] == 'DALY':
                norm[(i[0], i[1], category_names[i])] = {"normalisationFactor": 13.7, "weightingFactor": 5401.459854}
            elif i[1] == 'PDF.m2.yr':
                norm[(i[0], i[1], category_names[i])] = {"normalisationFactor": 1.01E-4, "weightingFactor": 1386.138614}

        normalization = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "NwSet",
            "@id": norm_weight_id,
            "name": "IMPACT World+ (Stepwise 2006 values)",
            "version": "00.00.000",
            "weightedScoreUnit": "EUR2003",
            "factors": []}

        for x in norm:
            normalization['factors'].append(
                {
                    "@type": "NwFactor",
                    "impactCategory": {
                        "@type": "ImpactCategory",
                        "@id": x[2],
                        "name": x[0],
                        "refUnit": x[1]
                    },
                    "normalisationFactor": norm[x]['normalisationFactor'],
                    "weightingFactor": norm[x]['weightingFactor']
                }
            )

        # ------------------------- STORE ALL DATA IN DICTIONARY ----------------------
        self.olca_data = {'category_metadata': category_metadata,
                          'metadata_iw_damage': metadata_iw_damage,
                          'metadata_iw_midpoint': metadata_iw_midpoint,
                          'metadata_iw_footprint': metadata_iw_footprint,
                          'metadata_iw_combined': metadata_iw_combined,
                          'cf_dict_damage': cf_dict_damage,
                          'cf_dict_midpoint': cf_dict_midpoint,
                          'cf_dict_footprint': cf_dict_footprint,
                          'cf_dict_combined': cf_dict_combined,
                          'normalization': normalization}

    def produce_files(self):
        """
        Function producing the different IW+ files for the different versions.
        :return: the IW+ files
        """

        self.logger.info("Creating all the files...")

        path = pkg_resources.resource_filename(__name__, '/Databases/Impact_world_' + self.version)

        # if the folders are not there yet, create them
        if not os.path.exists(path + '/Dev/'):
            os.makedirs(path + '/Dev/')
        if not os.path.exists(path + '/ecoinvent/'):
            os.makedirs(path + '/ecoinvent/')
        if not os.path.exists(path + '/exiobase/'):
            os.makedirs(path + '/exiobase/')
        if not os.path.exists(path + '/bw2/'):
            os.makedirs(path + '/bw2/')
        if not os.path.exists(path + '/SimaPro/'):
            os.makedirs(path + '/SimaPro/')
        if not os.path.exists(path + '/openLCA/'):
            os.makedirs(path + '/openLCA/')

        # Dev version
        self.master_db.to_excel(path + '/Dev/impact_world_plus_' + self.version + '_dev.xlsx')

        # ecoinvent versions in Excel format
        self.ei35_iw.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v35.xlsx')
        self.ei36_iw.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v36.xlsx')
        self.ei371_iw.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v371.xlsx')
        self.ei38_iw.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v38.xlsx')
        self.ei39_iw.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v39.xlsx')
        self.ei310_iw.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v310.xlsx')

        # ecoinvent version in DataFrame format
        self.ei35_iw_as_matrix.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v35_as_df.xlsx')
        self.ei36_iw_as_matrix.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v36_as_df.xlsx')
        self.ei371_iw_as_matrix.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v371_as_df.xlsx')
        self.ei38_iw_as_matrix.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v38_as_df.xlsx')
        self.ei39_iw_as_matrix.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v39_as_df.xlsx')
        self.ei310_iw_as_matrix.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v310_as_df.xlsx')

        # ecoinvent version in DataFrame format
        self.simplified_version_ei35.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v35.xlsx')
        self.simplified_version_ei36.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v36.xlsx')
        self.simplified_version_ei371.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v371.xlsx')
        self.simplified_version_ei38.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v38.xlsx')
        self.simplified_version_ei39.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v39.xlsx')
        self.simplified_version_ei310.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v310.xlsx')

        # exiobase version in DataFrame format
        self.exio_iw.to_excel(path + '/exiobase/impact_world_plus_' + self.version + '_expert_version_exiobase.xlsx')

        # brightway2 versions in bw2package format
        IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                                                                self.version in ic[0])]
        bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                             '_brightway2_expert_version', folder=path+'/bw2/')
        # bw2 footprint version
        IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' in ic[0] and
                                                                self.version in ic[0])]
        bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                             '_brightway2_footprint_version', folder=path+'/bw2/')

        # SimaPro version in csv format
        with open(path+'/SimaPro/impact_world_plus_'+self.version+'_midpoint_version_simapro_.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['midpoint_method_metadata'] +
                self.sp_data['midpoint_values'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/impact_world_plus_'+self.version+'_expert_version_simapro_.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['damage_method_metadata'] +
                self.sp_data['damage_values'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_damage'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/impact_world_plus_'+self.version+'_simapro.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['combined_method_metadata'] +
                self.sp_data['combined_values'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_combined'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/impact_world_plus_'+self.version+'_footprint_version_simapro_.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['simplified_method_metadata'] +
                self.sp_data['simplified_values'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])

        # create the openLCA version expert version (zip file)
        if not os.path.exists(path + '/openLCA/expert_version/'):
            os.makedirs(path + '/openLCA/expert_version/')
        if os.path.exists(path + '/openLCA/expert_version/impact_world_plus'+self.version+'_openLCA.zip'):
            os.remove(path + '/openLCA/expert_version/impact_world_plus'+self.version+'_openLCA.zip')
        if os.path.exists(path + '/openLCA/expert_version/oLCA_folders'):
            shutil.rmtree(path + '/openLCA/expert_version/oLCA_folders')
        zipObj = zipfile.ZipFile(path + '/openLCA/expert_version/impact_world_plus_'+self.version+'_openLCA.zip', 'w')

        if not os.path.exists(path + '/openLCA/expert_version/oLCA_folders/categories/'):
            os.makedirs(path + '/openLCA/expert_version/oLCA_folders/categories/')
        with open(path + '/openLCA/expert_version/oLCA_folders/categories/' + self.olca_data['category_metadata']['@id']
                  + '.json', 'w') as f:
            json.dump(self.olca_data['category_metadata'], f)
        zipObj.write(path + '/openLCA/expert_version/oLCA_folders/categories/' + self.olca_data['category_metadata'][
            '@id'] + '.json')

        if not os.path.exists(path + '/openLCA/expert_version/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/expert_version/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/expert_version/oLCA_folders/lcia_methods/' + self.olca_data['metadata_iw_damage']['@id'] +
                  '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw_damage'], f)
        zipObj.write(path + '/openLCA/expert_version/oLCA_folders/lcia_methods/' + self.olca_data['metadata_iw_damage']['@id']
                     + '.json')

        if not os.path.exists(path + '/openLCA/expert_version/oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/expert_version/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict_damage'].keys():
            with open(path + '/openLCA/expert_version/oLCA_folders/lcia_categories/' +
                      self.olca_data['cf_dict_damage'][cat]['@id'] + '.json', 'w') as f:
                json.dump(self.olca_data['cf_dict_damage'][cat], f)
            zipObj.write(path + '/openLCA/expert_version/oLCA_folders/lcia_categories/' + self.olca_data['cf_dict_damage'][
                cat]['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/expert_version/oLCA_folders/nw_sets/'):
            os.makedirs(path + '/openLCA/expert_version/oLCA_folders/nw_sets/')
        with open(path + '/openLCA/expert_version/oLCA_folders/nw_sets/' + self.olca_data['normalization']['@id'] +
                  '.json', 'w') as f:
            json.dump(self.olca_data['normalization'], f)
        zipObj.write(path + '/openLCA/expert_version/oLCA_folders/nw_sets/' + self.olca_data['normalization']['@id'] +
                     '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/expert_version/impact_world_plus_'+self.version+'_expert_version_openLCA', 'zip', path +
                            '/openLCA/expert_version/oLCA_folders/')

        # create the openLCA version midpoint version (zip file)
        if not os.path.exists(path + '/openLCA/midpoint_version/'):
            os.makedirs(path + '/openLCA/midpoint_version/')
        if os.path.exists(path + '/openLCA/midpoint_version/impact_world_plus'+self.version+'_openLCA.zip'):
            os.remove(path + '/openLCA/midpoint_version/impact_world_plus'+self.version+'_openLCA.zip')
        if os.path.exists(path + '/openLCA/midpoint_version/oLCA_folders'):
            shutil.rmtree(path + '/openLCA/midpoint_version/oLCA_folders')
        zipObj = zipfile.ZipFile(path + '/openLCA/midpoint_version/impact_world_plus_'+self.version+'_openLCA.zip', 'w')

        if not os.path.exists(path + '/openLCA/midpoint_version/oLCA_folders/categories/'):
            os.makedirs(path + '/openLCA/midpoint_version/oLCA_folders/categories/')
        with open(path + '/openLCA/midpoint_version/oLCA_folders/categories/' + self.olca_data['category_metadata']['@id']
                  + '.json', 'w') as f:
            json.dump(self.olca_data['category_metadata'], f)
        zipObj.write(path + '/openLCA/midpoint_version/oLCA_folders/categories/' + self.olca_data['category_metadata'][
            '@id'] + '.json')

        if not os.path.exists(path + '/openLCA/midpoint_version/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/midpoint_version/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/midpoint_version/oLCA_folders/lcia_methods/' + self.olca_data['metadata_iw_midpoint']['@id'] +
                  '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw_midpoint'], f)
        zipObj.write(path + '/openLCA/midpoint_version/oLCA_folders/lcia_methods/' + self.olca_data['metadata_iw_midpoint']['@id']
                     + '.json')

        if not os.path.exists(path + '/openLCA/midpoint_version/oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/midpoint_version/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict_midpoint'].keys():
            with open(path + '/openLCA/midpoint_version/oLCA_folders/lcia_categories/' +
                      self.olca_data['cf_dict_midpoint'][cat]['@id'] + '.json', 'w') as f:
                json.dump(self.olca_data['cf_dict_midpoint'][cat], f)
            zipObj.write(path + '/openLCA/midpoint_version/oLCA_folders/lcia_categories/' + self.olca_data['cf_dict_midpoint'][
                cat]['@id'] + '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/midpoint_version/impact_world_plus_'+self.version+'_midpoint_version_openLCA', 'zip', path +
                            '/openLCA/midpoint_version/oLCA_folders/')

        # create the openLCA version footprint version (zip file)
        if not os.path.exists(path + '/openLCA/footprint_version/'):
            os.makedirs(path + '/openLCA/footprint_version/')
        if os.path.exists(path + '/openLCA/footprint_version/impact_world_plus'+self.version+'_openLCA.zip'):
            os.remove(path + '/openLCA/footprint_version/impact_world_plus'+self.version+'_openLCA.zip')
        if os.path.exists(path + '/openLCA/footprint_version/oLCA_folders'):
            shutil.rmtree(path + '/openLCA/footprint_version/oLCA_folders')
        zipObj = zipfile.ZipFile(path + '/openLCA/footprint_version/impact_world_plus_'+self.version+'_openLCA.zip', 'w')

        if not os.path.exists(path + '/openLCA/footprint_version/oLCA_folders/categories/'):
            os.makedirs(path + '/openLCA/footprint_version/oLCA_folders/categories/')
        with open(path + '/openLCA/footprint_version/oLCA_folders/categories/' + self.olca_data['category_metadata']['@id']
                  + '.json', 'w') as f:
            json.dump(self.olca_data['category_metadata'], f)
        zipObj.write(path + '/openLCA/footprint_version/oLCA_folders/categories/' + self.olca_data['category_metadata'][
            '@id'] + '.json')

        if not os.path.exists(path + '/openLCA/footprint_version/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/footprint_version/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/footprint_version/oLCA_folders/lcia_methods/' +
                  self.olca_data['metadata_iw_footprint']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw_footprint'], f)
        zipObj.write(path + '/openLCA/footprint_version/oLCA_folders/lcia_methods/' +
                     self.olca_data['metadata_iw_footprint']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/footprint_version/oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/footprint_version/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict_footprint'].keys():
            with open(path + '/openLCA/footprint_version/oLCA_folders/lcia_categories/' +
                      self.olca_data['cf_dict_footprint'][cat]['@id'] + '.json', 'w') as f:
                json.dump(self.olca_data['cf_dict_footprint'][cat], f)
            zipObj.write(path + '/openLCA/footprint_version/oLCA_folders/lcia_categories/' +
                         self.olca_data['cf_dict_footprint'][cat]['@id'] + '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/footprint_version/impact_world_plus_' + self.version +
                            '_footprint_version_openLCA', 'zip', path + '/openLCA/footprint_version/oLCA_folders/')

        # create the openLCA version combined version (zip file)
        if not os.path.exists(path + '/openLCA/combined_version/'):
            os.makedirs(path + '/openLCA/combined_version/')
        if os.path.exists(path + '/openLCA/combined_version/impact_world_plus'+self.version+'_openLCA.zip'):
            os.remove(path + '/openLCA/combined_version/impact_world_plus'+self.version+'_openLCA.zip')
        if os.path.exists(path + '/openLCA/combined_version/oLCA_folders'):
            shutil.rmtree(path + '/openLCA/combined_version/oLCA_folders')
        zipObj = zipfile.ZipFile(path + '/openLCA/combined_version/impact_world_plus_'+self.version+'_openLCA.zip', 'w')

        if not os.path.exists(path + '/openLCA/combined_version/oLCA_folders/categories/'):
            os.makedirs(path + '/openLCA/combined_version/oLCA_folders/categories/')
        with open(path + '/openLCA/combined_version/oLCA_folders/categories/' + self.olca_data['category_metadata']['@id']
                  + '.json', 'w') as f:
            json.dump(self.olca_data['category_metadata'], f)
        zipObj.write(path + '/openLCA/combined_version/oLCA_folders/categories/' + self.olca_data['category_metadata'][
            '@id'] + '.json')

        if not os.path.exists(path + '/openLCA/combined_version/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/combined_version/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/combined_version/oLCA_folders/lcia_methods/' + self.olca_data['metadata_iw_combined']['@id'] +
                  '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw_combined'], f)
        zipObj.write(path + '/openLCA/combined_version/oLCA_folders/lcia_methods/' + self.olca_data['metadata_iw_combined']['@id']
                     + '.json')

        if not os.path.exists(path + '/openLCA/combined_version/oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/combined_version/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict_combined'].keys():
            with open(path + '/openLCA/combined_version/oLCA_folders/lcia_categories/' +
                      self.olca_data['cf_dict_combined'][cat]['@id'] + '.json', 'w') as f:
                json.dump(self.olca_data['cf_dict_combined'][cat], f)
            zipObj.write(path + '/openLCA/combined_version/oLCA_folders/lcia_categories/' + self.olca_data['cf_dict_combined'][
                cat]['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/combined_version/oLCA_folders/nw_sets/'):
            os.makedirs(path + '/openLCA/combined_version/oLCA_folders/nw_sets/')
        with open(path + '/openLCA/combined_version/oLCA_folders/nw_sets/' + self.olca_data['normalization']['@id'] +
                  '.json', 'w') as f:
            json.dump(self.olca_data['normalization'], f)
        zipObj.write(path + '/openLCA/combined_version/oLCA_folders/nw_sets/' + self.olca_data['normalization']['@id'] +
                     '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/combined_version/impact_world_plus_' + self.version +
                            '_combined_version_openLCA', 'zip', path + '/openLCA/combined_version/oLCA_folders/')

    def produce_files_hybrid_ecoinvent(self):
        """Specific method to create the files matching with hybrid-ecoinvent (pylcaio)."""

        path = pkg_resources.resource_filename(__name__, '/Databases/Impact_world_' + self.version)

        if not os.path.exists(path + '/for_hybrid_ecoinvent/ei35/'):
            os.makedirs(path + '/for_hybrid_ecoinvent/ei35/')
        if not os.path.exists(path + '/for_hybrid_ecoinvent/ei36/'):
            os.makedirs(path + '/for_hybrid_ecoinvent/ei36/')
        if not os.path.exists(path + '/for_hybrid_ecoinvent/ei371/'):
            os.makedirs(path + '/for_hybrid_ecoinvent/ei371/')
        if not os.path.exists(path + '/for_hybrid_ecoinvent/ei38/'):
            os.makedirs(path + '/for_hybrid_ecoinvent/ei38/')
        if not os.path.exists(path + '/for_hybrid_ecoinvent/ei39/'):
            os.makedirs(path + '/for_hybrid_ecoinvent/ei39/')

        scipy.sparse.save_npz(path+'/for_hybrid_ecoinvent/ei35/Ecoinvent_not_regionalized.npz',
                              scipy.sparse.csr_matrix(self.ei35_iw_as_matrix))
        scipy.sparse.save_npz(path+'/for_hybrid_ecoinvent/ei36/Ecoinvent_not_regionalized.npz',
                              scipy.sparse.csr_matrix(self.ei36_iw_as_matrix))
        scipy.sparse.save_npz(path+'/for_hybrid_ecoinvent/ei371/Ecoinvent_not_regionalized.npz',
                              scipy.sparse.csr_matrix(self.ei371_iw_as_matrix))
        scipy.sparse.save_npz(path+'/for_hybrid_ecoinvent/ei38/Ecoinvent_not_regionalized.npz',
                              scipy.sparse.csr_matrix(self.ei38_iw_as_matrix))
        scipy.sparse.save_npz(path+'/for_hybrid_ecoinvent/ei39/Ecoinvent_not_regionalized.npz',
                              scipy.sparse.csr_matrix(self.ei39_iw_as_matrix))

# ----------------------------------------- Secondary methods ----------------------------------------------------------

    def load_basic_cfs(self):
        """
        Loading the basic CFs. By basic we mean that these CFs do not require further treatment.

        Concerned impact categories:
            - Fossil and nuclear energy use
            - Freshwater ecotoxicity
            - Freshwater ecotoxicity, long term
            - Freshwater ecotoxicity, short term
            - Human toxicity cancer
            - Human toxicity cancer, long term
            - Human toxicity cancer, short term
            - Human toxicity non cancer
            - Human toxicity non cancer, long term
            - Human toxicity non cancer, short term
            - Ionizing radiation, ecosystem quality
            - Ionizing radiation, human health
            - Ionizing radiations
            - Marine acidification, short term
            = Marine acidification, long term
            - Mineral resources use
            - Ozone layer depletion
            - Photochemical oxidant formation

        :return: updated master_db
        """

        self.master_db = pd.concat([pd.read_sql(sql='SELECT * FROM [CF - not regionalized - IonizingRadiations]',
                                                con=self.conn),
                                    pd.read_sql(sql='SELECT * FROM [CF - not regionalized - MarAcid]',
                                                con=self.conn),
                                    pd.read_sql(sql='SELECT * FROM [CF - not regionalized - FossilResources]',
                                                con=self.conn),
                                    pd.read_sql(sql='SELECT * FROM [CF - not regionalized - MineralResources]',
                                                con=self.conn),
                                    pd.read_sql(sql='SELECT * FROM [CF - not regionalized - HumanTox]',
                                                con=self.conn),
                                    pd.read_sql(sql='SELECT * FROM [CF - not regionalized - EcotoxFW]',
                                                con=self.conn)])

        self.master_db = clean_up_dataframe(self.master_db)

    def load_climate_change_cfs(self):
        """
        Loading the CFs for the climate change impact categories.

        Concerned impact categories:
            - Climate change, short term
            - Climate change, long term
            - Climate change, human health, short term
            - Climate change, human health, long term
            - Climate change, ecosystem quality, short term
            - Climate change, ecosystem quality, long term

        :return: updated master_db
        """

        data = pd.read_sql('SELECT * FROM "CF - not regionalized - ClimateChange"', con=self.conn)
        mapping = pd.read_sql('SELECT * FROM "SI - Climate change - mapping GHGs"', con=self.conn)
        data = data.merge(mapping, left_on='Name', right_on='IPCC_name', how='inner').drop(['Name', 'IPCC_name'],axis=1)

        # ---------------------------- Climate change midpoint indicators ---------------------------------------------
        # Climate change, short term
        GWP_midpoint = data.loc[:, ['IW+_name', 'GWP-100', 'CAS number']]
        GWP_midpoint.columns = ['Elem flow name', 'CF value', 'CAS number']
        GWP_midpoint.loc[:, 'Impact category'] = 'Climate change, short term'
        GWP_midpoint.loc[:, 'CF unit'] = 'kg CO2 eq (short)'
        GWP_midpoint.loc[:, 'Compartment'] = 'Air'
        GWP_midpoint.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_midpoint.loc[:, 'Elem flow unit'] = 'kg'
        GWP_midpoint.loc[:, 'MP or Damage'] = 'Midpoint'
        GWP_midpoint.loc[:, 'Native geographical resolution scale'] = 'Global'

        # Climate change, long term
        GTP_midpoint = data.loc[:, ['IW+_name', 'GTP-100', 'CAS number']]
        GTP_midpoint.columns = ['Elem flow name', 'CF value', 'CAS number']
        GTP_midpoint.loc[:, 'Impact category'] = 'Climate change, long term'
        GTP_midpoint.loc[:, 'CF unit'] = 'kg CO2 eq (long)'
        GTP_midpoint.loc[:, 'Compartment'] = 'Air'
        GTP_midpoint.loc[:, 'Sub-compartment'] = '(unspecified)'
        GTP_midpoint.loc[:, 'Elem flow unit'] = 'kg'
        GTP_midpoint.loc[:, 'MP or Damage'] = 'Midpoint'
        GTP_midpoint.loc[:, 'Native geographical resolution scale'] = 'Global'

        # ---------------------------- Climate change damage indicators -----------------------------------------------

        # atmosphere properties
        M_ATMOS = 5.1352E18  # kg
        M_AIR = 28.97E-3  # kg/mol

        # ----------- IRFCO2 -----------#
        # Joos 2013
        aC1 = 0.2173
        aC2 = 0.2763
        aC3 = 0.2824
        aC4 = 0.2240
        tauC1 = 4.304
        tauC2 = 36.54
        tauC3 = 394.4

        # ----------- IRFclim -----------#
        # IPCC AR6 github
        kPulseT = 0.7634
        aT1 = 0.5809
        aT2 = 0.4191
        tauT1 = 3.444
        tauT2 = 285.1

        # ----------- IRFfdbk -----------#
        # Gasser 2017
        gamma = 3.015
        aS1 = 0.6368
        aS2 = 0.3322
        aS3 = 0.0310
        tauS1 = 2.376
        tauS2 = 30.14
        tauS3 = 490.1

        # calculate molecular mass of each GHG
        data.loc[:, 'Molecular Mass (kg/mol)'] = [molmass.Formula(i).mass / 1000 for i in data.Formula]
        # convert radiative efficiencies from W/m2/ppb to W/m2/kg for all GHGs
        data.loc[:, 'Radiative Efficiency (W/m2/kg)'] = data.loc[:, 'Radiative Efficiency (W/m2/ppb)'] / (
                    1E-9 * (data.loc[:, 'Molecular Mass (kg/mol)'] / M_AIR) * M_ATMOS)

        AACO2 = data.loc[data.loc[:, 'IW+_name'] == 'Carbon dioxide, fossil', 'Radiative Efficiency (W/m2/kg)'].iloc[0]

        # AGTP cumulative 100 years
        data.loc[:, 'AGTP cumulative 100 years'] = 0

        for t in range(0, 100):
            data.loc[data[data.loc[:, 'IW+_name'] == 'Carbon dioxide, fossil'].index[
                         0], 'AGTP cumulative 100 years'] += AGTPCO2(t, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3,
                                                                     kPulseT, aT1, tauT1, aT2, tauT2, AACO2)

            data.loc[data[data.loc[:, 'IW+_name'] == 'Methane, fossil'].index[
                         0], 'AGTP cumulative 100 years'] += AGTPCH4Fossil_Final(
                t, data.loc[data.loc[:, 'IW+_name'] == 'Methane, fossil', 'Lifetime (yr)'].iloc[0], kPulseT, aT1, tauT1,
                aT2, tauT2,
                data.loc[data.loc[:, 'IW+_name'] == 'Methane, fossil', 'Radiative Efficiency (W/m2/kg)'].iloc[0], aC1,
                aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

            data.loc[data[data.loc[:, 'IW+_name'] == 'Methane, biogenic'].index[
                         0], 'AGTP cumulative 100 years'] += AGTPCH4NonFossil_Final(
                t, data.loc[data.loc[:, 'IW+_name'] == 'Methane, biogenic', 'Lifetime (yr)'].iloc[0], kPulseT, aT1,
                tauT1, aT2, tauT2,
                data.loc[data.loc[:, 'IW+_name'] == 'Methane, biogenic', 'Radiative Efficiency (W/m2/kg)'].iloc[0], aC1,
                aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

            for ghg in data.index:
                if data.loc[ghg, 'IW+_name'] not in ['Carbon dioxide, fossil', 'Methane, fossil', 'Methane, biogenic']:
                    if data.loc[ghg, 'Lifetime (yr)'] != 0:
                        data.loc[ghg, 'AGTP cumulative 100 years'] += AGTPNonCO2_Final(
                            t, data.loc[ghg, 'Lifetime (yr)'], kPulseT, aT1, tauT1, aT2, tauT2,
                            data.loc[ghg, 'Radiative Efficiency (W/m2/kg)'], aC1, aC2, aC3, aC4, tauC1, tauC2,
                            tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

        # AGTP cumulative 500 years
        data.loc[:, 'AGTP cumulative 500 years'] = 0

        for t in range(0, 500):
            data.loc[data[data.loc[:, 'IW+_name'] == 'Carbon dioxide, fossil'].index[
                         0], 'AGTP cumulative 500 years'] += AGTPCO2(t, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3,
                                                                     kPulseT, aT1, tauT1, aT2, tauT2, AACO2)

            data.loc[data[data.loc[:, 'IW+_name'] == 'Methane, fossil'].index[
                         0], 'AGTP cumulative 500 years'] += AGTPCH4Fossil_Final(
                t, data.loc[data.loc[:, 'IW+_name'] == 'Methane, fossil', 'Lifetime (yr)'].iloc[0], kPulseT, aT1, tauT1,
                aT2, tauT2,
                data.loc[data.loc[:, 'IW+_name'] == 'Methane, fossil', 'Radiative Efficiency (W/m2/kg)'].iloc[0], aC1,
                aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

            data.loc[data[data.loc[:, 'IW+_name'] == 'Methane, biogenic'].index[
                         0], 'AGTP cumulative 500 years'] += AGTPCH4NonFossil_Final(
                t, data.loc[data.loc[:, 'IW+_name'] == 'Methane, biogenic', 'Lifetime (yr)'].iloc[0], kPulseT, aT1,
                tauT1, aT2, tauT2,
                data.loc[data.loc[:, 'IW+_name'] == 'Methane, biogenic', 'Radiative Efficiency (W/m2/kg)'].iloc[0], aC1,
                aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

            for ghg in data.index:
                if data.loc[ghg, 'IW+_name'] not in ['Carbon dioxide, fossil', 'Methane, fossil', 'Methane, biogenic']:
                    if data.loc[ghg, 'Lifetime (yr)'] != 0:
                        data.loc[ghg, 'AGTP cumulative 500 years'] += AGTPNonCO2_Final(
                            t, data.loc[ghg, 'Lifetime (yr)'], kPulseT, aT1, tauT1, aT2, tauT2,
                            data.loc[ghg, 'Radiative Efficiency (W/m2/kg)'],aC1, aC2, aC3, aC4, tauC1, tauC2,tauC3,
                            AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

        # get effect factors
        effect_factors = pd.read_sql('SELECT * FROM "SI - Climate change - effect factors"', con=self.conn).set_index('index')
        HH_effect_factor = effect_factors.loc['Total', 'DALY/K/yr']
        EQ_effect_factor = effect_factors.loc['Total', 'PDF.m2.yr/K/yr']

        # Climate change, human health, short term
        GWP_damage_HH_short = data.loc[:, ['IW+_name', 'AGTP cumulative 100 years', 'CAS number']].copy()
        GWP_damage_HH_short.columns = ['Elem flow name', 'CF value', 'CAS number']
        GWP_damage_HH_short.loc[:, 'Impact category'] = 'Climate change, human health, short term'
        GWP_damage_HH_short.loc[:, 'CF unit'] = 'DALY'
        GWP_damage_HH_short.loc[:, 'Compartment'] = 'Air'
        GWP_damage_HH_short.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_damage_HH_short.loc[:, 'Elem flow unit'] = 'kg'
        GWP_damage_HH_short.loc[:, 'MP or Damage'] = 'Damage'
        GWP_damage_HH_short.loc[:, 'Native geographical resolution scale'] = 'Global'
        GWP_damage_HH_short.loc[:, 'CF value'] *= HH_effect_factor

        # Climate change, human health, long term
        GWP_damage_HH_long = data.loc[:, ['IW+_name', 'AGTP cumulative 100 years', 'AGTP cumulative 500 years',
                                          'CAS number']].copy()
        GWP_damage_HH_long = GWP_damage_HH_long.rename(columns={'IW+_name': 'Elem flow name'})
        GWP_damage_HH_long.loc[:, 'Impact category'] = 'Climate change, human health, long term'
        GWP_damage_HH_long.loc[:, 'CF unit'] = 'DALY'
        GWP_damage_HH_long.loc[:, 'Compartment'] = 'Air'
        GWP_damage_HH_long.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_damage_HH_long.loc[:, 'Elem flow unit'] = 'kg'
        GWP_damage_HH_long.loc[:, 'MP or Damage'] = 'Damage'
        GWP_damage_HH_long.loc[:, 'Native geographical resolution scale'] = 'Global'
        GWP_damage_HH_long.loc[:, 'CF value'] = (GWP_damage_HH_long.loc[:,'AGTP cumulative 500 years'] -
                                                 GWP_damage_HH_long.loc[:,'AGTP cumulative 100 years']) * HH_effect_factor
        GWP_damage_HH_long = GWP_damage_HH_long.drop(['AGTP cumulative 100 years', 'AGTP cumulative 500 years'], axis=1)

        # Climate change, ecosystem quality, short term
        GWP_damage_EQ_short = data.loc[:, ['IW+_name', 'AGTP cumulative 100 years', 'CAS number']]
        GWP_damage_EQ_short.columns = ['Elem flow name', 'CF value', 'CAS number']
        GWP_damage_EQ_short.loc[:, 'Impact category'] = 'Climate change, ecosystem quality, short term'
        GWP_damage_EQ_short.loc[:, 'CF unit'] = 'PDF.m2.yr'
        GWP_damage_EQ_short.loc[:, 'Compartment'] = 'Air'
        GWP_damage_EQ_short.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_damage_EQ_short.loc[:, 'Elem flow unit'] = 'kg'
        GWP_damage_EQ_short.loc[:, 'MP or Damage'] = 'Damage'
        GWP_damage_EQ_short.loc[:, 'Native geographical resolution scale'] = 'Global'
        GWP_damage_EQ_short.loc[:, 'CF value'] *= EQ_effect_factor

        # Climate change, ecosystem quality, long term
        GWP_damage_EQ_long = data.loc[:, ['IW+_name', 'AGTP cumulative 100 years', 'AGTP cumulative 500 years',
                                          'CAS number']].copy()
        GWP_damage_EQ_long = GWP_damage_EQ_long.rename(columns={'IW+_name': 'Elem flow name'})
        GWP_damage_EQ_long.loc[:, 'Impact category'] = 'Climate change, ecosystem quality, long term'
        GWP_damage_EQ_long.loc[:, 'CF unit'] = 'PDF.m2.yr'
        GWP_damage_EQ_long.loc[:, 'Compartment'] = 'Air'
        GWP_damage_EQ_long.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_damage_EQ_long.loc[:, 'Elem flow unit'] = 'kg'
        GWP_damage_EQ_long.loc[:, 'MP or Damage'] = 'Damage'
        GWP_damage_EQ_long.loc[:, 'Native geographical resolution scale'] = 'Global'
        GWP_damage_EQ_long.loc[:, 'CF value'] = (GWP_damage_EQ_long.loc[:,'AGTP cumulative 500 years'] -
                                                 GWP_damage_EQ_long.loc[:,'AGTP cumulative 100 years']) * EQ_effect_factor
        GWP_damage_EQ_long = GWP_damage_EQ_long.drop(['AGTP cumulative 100 years', 'AGTP cumulative 500 years'], axis=1)

        self.master_db = pd.concat([self.master_db, GWP_midpoint, GTP_midpoint, GWP_damage_EQ_short, GWP_damage_EQ_long,
                                    GWP_damage_HH_short, GWP_damage_HH_long])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_ozone_layer_depletion_cfs(self):
        """
        Loading the CFs for the ozone layer depletion impact categories.

        Concerned impact categories:
            - Ozone layer depletion

        :return: updated master_db
        """

        data = pd.read_sql('SELECT * FROM "CF - not regionalized - OzoneLayerDepletion"', self.conn).set_index('index')

        # reformat
        data = data.loc[:, ['CAS', 'ODP (infinite)']].reset_index()
        data = data.rename(columns={'index': 'Elem flow name', 'CAS': 'CAS number', 'ODP (infinite)': 'CF value'})
        data.loc[:, 'Impact category'] = 'Ozone layer depletion'
        data.loc[:, 'Compartment'] = 'Air'
        data.loc[:, 'Sub-compartment'] = '(unspecified)'
        data.loc[:, 'CF unit'] = 'kg CFC-11 eq'
        data.loc[:, 'Elem flow unit'] = 'kg'
        data.loc[:, 'MP or Damage'] = 'Midpoint'
        data.loc[:, 'Native geographical resolution scale'] = 'Not regionalized'
        # map names to IW+ standard
        map = pd.read_sql('SELECT * FROM "SI - Ozone Layer Depletion - mapping"', self.conn)
        data.loc[:, 'Elem flow name'] = [dict(zip(map.WMO_name, map.loc[:, "IW+_name"]))[i] for i in
                                         data.loc[:, 'Elem flow name']]

        # originally coming from Hayashi 2006, says ReCiPe2016 report adopting Egalitarian perspective
        midpoint_to_damage_factor_long_term = 1.34e-3
        damage_data = data.copy('deep')
        damage_data.loc[:, 'CF value'] *= midpoint_to_damage_factor_long_term
        damage_data.loc[:, 'CF unit'] = 'DALY'
        damage_data.loc[:, 'MP or Damage'] = 'Damage'

        data = pd.concat([data, damage_data])
        data.index = [i for i in range(0, len(data.index))]

        # concat with master_db
        self.master_db = pd.concat([self.master_db, data])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_photochemical_ozone_formation(self):
        """
        Loading the CFs for the photochemical ozone formation impact categories.

        Concerned impact categories:
            - Photochemical ozone formation

        :return: updated master_db
        """

        data = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - PhotochemOxid]', con=self.conn)
        conc = pd.read_sql(sql='SELECT * FROM [SI - Photochemical Ozone Formation - mapping]', con=self.conn)
        effect_factor = pd.read_sql(sql='SELECT * FROM [SI - Photochemical Ozone Formation - effect factors]',
                                    con=self.conn)
        data = data.merge(conc, left_on=['Substance name'], right_on=['ReCiPe2016 name'])

        photochem_damage_hh = data.copy('deep')
        photochem_damage_hh = photochem_damage_hh.loc[:, ['CF HH (kg NOx-eq/kg)', 'IW+ name', 'CAS number']]
        photochem_damage_hh = photochem_damage_hh.rename(columns={'CF HH (kg NOx-eq/kg)': 'CF value',
                                                                  'IW+ name': 'Elem flow name'})
        photochem_damage_hh.loc[:, 'CF value'] *= effect_factor.loc[0, 'HH (DALY/kg NOx-eq)']
        photochem_damage_hh.loc[:, 'Impact category'] = 'Photochemical Ozone Formation, human health'
        photochem_damage_hh.loc[:, 'CF unit'] = 'DALY'
        photochem_damage_hh.loc[:, 'Compartment'] = 'Air'
        photochem_damage_hh.loc[:, 'Sub-compartment'] = '(unspecified)'
        photochem_damage_hh.loc[:, 'Elem flow unit'] = 'kg'
        photochem_damage_hh.loc[:, 'MP or Damage'] = 'Damage'
        photochem_damage_hh.loc[:, 'Native geographical resolution scale'] = 'Not regionalized'

        photochem_damage_eq = data.copy('deep')
        photochem_damage_eq = photochem_damage_eq.loc[:, ['CF EQ (kg NOx-eq/kg)', 'IW+ name', 'CAS number']]
        photochem_damage_eq = photochem_damage_eq.rename(columns={'CF EQ (kg NOx-eq/kg)': 'CF value',
                                                                  'IW+ name': 'Elem flow name'})
        photochem_damage_eq.loc[:, 'CF value'] *= (effect_factor.loc[0, 'EQ (species.yr/kg NOx-eq)'] /
                                                   effect_factor.loc[0, 'species density in ReCiPe (species/m2)'])
        photochem_damage_eq.loc[:, 'Impact category'] = 'Photochemical Ozone Formation, ecosystem quality'
        photochem_damage_eq.loc[:, 'CF unit'] = 'PDF.m2.yr'
        photochem_damage_eq.loc[:, 'Compartment'] = 'Air'
        photochem_damage_eq.loc[:, 'Sub-compartment'] = '(unspecified)'
        photochem_damage_eq.loc[:, 'Elem flow unit'] = 'kg'
        photochem_damage_eq.loc[:, 'MP or Damage'] = 'Damage'
        photochem_damage_eq.loc[:, 'Native geographical resolution scale'] = 'Not regionalized'

        # ReCiPe separates HH and EQ at midpoint. We only base our midpoint on  HH.
        photochem_midpoint = data.copy('deep')
        photochem_midpoint.loc[:, 'CF value'] = photochem_midpoint.loc[:, 'CF HH (kg NOx-eq/kg)']
        photochem_midpoint = photochem_midpoint.rename(columns={'IW+ name': 'Elem flow name'})
        photochem_midpoint = photochem_midpoint.loc[:, ['CF value', 'Elem flow name', 'CAS number']]
        photochem_midpoint.loc[:, 'Impact category'] = 'Photochemical Ozone Formation'
        photochem_midpoint.loc[:, 'CF unit'] = 'kgNOxeq'
        photochem_midpoint.loc[:, 'Compartment'] = 'Air'
        photochem_midpoint.loc[:, 'Sub-compartment'] = '(unspecified)'
        photochem_midpoint.loc[:, 'Elem flow unit'] = 'kg'
        photochem_midpoint.loc[:, 'MP or Damage'] = 'Midpoint'
        photochem_midpoint.loc[:, 'Native geographical resolution scale'] = 'Not regionalized'

        # concat with master_db
        self.master_db = pd.concat([self.master_db, photochem_midpoint, photochem_damage_hh, photochem_damage_eq])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_freshwater_acidification_cfs(self):
        """
        Loading the CFs for the freshwater acidification impact category. This includes CFs coming from the
        original article of IW+, as well as other CFs extrapolated from stoechiometric ratios.

        Concerned impact categories:
            - Freshwater acidification

        :return: updated master_db
        """

        # ------------------------------ LOADING DATA -----------------------------------
        cell_emissions = {
            'NH3': pd.read_sql('SELECT * FROM [CF - regionalized - AcidFW - native (NH3 emissions to air)]',
                               self.conn).set_index('cell'),
            'NOx': pd.read_sql('SELECT * FROM [CF - regionalized - AcidFW - native (NOx emissions to air)]',
                               self.conn).set_index('cell'),
            'SO2': pd.read_sql('SELECT * FROM [CF - regionalized - AcidFW - native (SO2 emissions to air)]',
                               self.conn).set_index('cell')}

        inter_country = pd.read_sql('SELECT * FROM [SI - Acidification/Eutrophication - Countries cell resolution]',
                                    self.conn).set_index('cell').T
        inter_continent = pd.read_sql('SELECT * FROM [SI - Acidification/Eutrophication - Continents cell resolution]',
                                      self.conn).set_index('cell').T

        # ------------------------- CALCULATING MIDPOINT AND DAMAGE CFS --------------------------------
        country_emissions = {}
        continent_emissions = {}
        for substance in cell_emissions:
            cell_emissions[substance].index.name = None
            country_emissions[substance] = inter_country.dot(cell_emissions[substance].Emission)
            continent_emissions[substance] = inter_continent.dot(cell_emissions[substance].Emission)

        cfs_midpoint = pd.DataFrame(None, cell_emissions.keys(), country_emissions['SO2'].index)

        for substance in cell_emissions:
            # for countries
            for country in country_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'midpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     country_emissions[substance].loc[country] *
                                     inter_country.loc[country]).sum()
                    cfs_midpoint.loc[substance, country] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # for continents
            for continent in continent_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'midpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     continent_emissions[substance].loc[continent] *
                                     inter_continent.loc[continent]).sum()
                    cfs_midpoint.loc[substance, continent] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # global
            GLO_value = (cell_emissions[substance].loc[:, 'midpoint'] *
                         cell_emissions[substance].loc[:, 'Emission'] /
                         cell_emissions[substance].Emission.sum()).sum()
            cfs_midpoint.loc[substance, 'GLO'] = GLO_value

        # divide by reference flow -> SO2
        cfs_midpoint /= cfs_midpoint.loc['SO2', 'GLO']

        cfs_damage = pd.DataFrame(None, cell_emissions.keys(), country_emissions['SO2'].index)

        for substance in cell_emissions:
            # for countries
            for country in country_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'endpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     country_emissions[substance].loc[country] *
                                     inter_country.loc[country]).sum()
                    cfs_damage.loc[substance, country] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # for continents
            for continent in continent_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'endpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     continent_emissions[substance].loc[continent] *
                                     inter_continent.loc[continent]).sum()
                    cfs_damage.loc[substance, continent] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # global
            GLO_value = (cell_emissions[substance].loc[:, 'endpoint'] *
                         cell_emissions[substance].loc[:, 'Emission'] /
                         cell_emissions[substance].Emission.sum()).sum()
            cfs_damage.loc[substance, 'GLO'] = GLO_value

        # ------------------------------ DATA FORMATTING -----------------------------------
        # convert country names to ISO 2-letter codes
        # lines of logger to avoid having warnings showing up
        coco_logger = coco.logging.getLogger()
        coco_logger.setLevel(coco.logging.CRITICAL)
        cfs_midpoint.columns = coco.convert(cfs_midpoint.iloc[:, :-7], to='ISO2') + ['RNA', 'RLA', 'RER', 'RAS', 'RAF',
                                                                                     'OCE', 'GLO']
        cfs_damage.columns = cfs_midpoint.columns
        # countries that don't exist anymore or that are not part of the UN are dropped
        cfs_midpoint.drop('not found', axis=1, inplace=True)
        cfs_damage.drop('not found', axis=1, inplace=True)

        cfs_midpoint = cfs_midpoint.stack().reset_index()
        cfs_midpoint.columns = ['Elem flow', 'Region code', 'CF value']
        cfs_damage = cfs_damage.stack().reset_index()
        cfs_damage.columns = ['Elem flow', 'Region code', 'CF value']
        cfs_midpoint.loc[:, 'CF unit'] = 'kg SO2 eq'
        cfs_midpoint.loc[:, 'MP or Damage'] = 'Midpoint'
        cfs_damage.loc[:, 'CF unit'] = 'PDF.m2.yr'
        cfs_damage.loc[:, 'MP or Damage'] = 'Damage'
        cfs = pd.concat([cfs_midpoint, cfs_damage])
        cfs.loc[:, 'Impact category'] = 'Freshwater acidification'
        cfs.loc[:, 'Compartment'] = 'Air'
        cfs.loc[:, 'Sub-compartment'] = '(unspecified)'
        cfs.loc[:, 'Elem flow unit'] = 'kg'
        cfs.loc[:, 'Native geographical resolution scale'] = 'Country'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NH3', 'CAS number'] = '007664-41-7'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NOx', 'CAS number'] = '011104-93-1'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'SO2', 'CAS number'] = '007446-09-5'
        cfs = clean_up_dataframe(cfs)
        cfs.loc[[i for i in cfs.index if cfs.loc[i, 'Region code'] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF', 'OCE', 'GLO']],
                'Native geographical resolution scale'] = 'Continent'
        cfs.loc[cfs.loc[:, 'Region code'] == 'GLO', 'Native geographical resolution scale'] = 'Global'

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - stoechiometry]', self.conn)
        for i in stoc.index:
            proxy = stoc.loc[i, 'Corresponding elem flow in IW']
            df = cfs[cfs.loc[:, 'Elem flow'] == proxy].copy('deep')
            if not df.empty:
                df.loc[:, 'Elem flow'] = stoc.loc[i, 'Elem flow name in Simapro']
                df.loc[:, 'CAS number'] = stoc.loc[i, 'CAS number']
                df.loc[:, 'CF value'] *= stoc.loc[i, 'Stoechiometric ratio']

                cfs = pd.concat([cfs, df])

        # ------------------------------------ FINAL FORMATTING ---------------------------------
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NH3', 'Elem flow'] = 'Ammonia'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NOx', 'Elem flow'] = 'Nitrogen oxides'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'SO2', 'Elem flow'] = 'Sulfur dioxide'
        # concat elem flow and region code in the same name
        cfs.loc[:, 'Elem flow name'] = [', '.join(i) for i in list(zip(cfs.loc[:, 'Elem flow'], cfs.loc[:, 'Region code']))]
        cfs = cfs.drop(['Elem flow', 'Region code'], axis=1)
        cfs = clean_up_dataframe(cfs)

        # concat with master_db
        self.master_db = pd.concat([self.master_db, cfs])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_terrestrial_acidification_cfs(self):
        """
        Loading the CFs for the terrestrial acidification impact category. This includes CFs coming from the
        original article of IW+, as well as other CFs extrapolated from stoechiometric ratios.

        Concerned impact categories:
            - Terrestrial acidification

        :return: updated master_db
        """

        # ------------------------------ LOADING DATA -----------------------------------
        cell_emissions = {
            'NH3': pd.read_sql('SELECT * FROM [CF - regionalized - AcidTerr - native (NH3 emissions to air)]',
                               self.conn).set_index('cell'),
            'NOx': pd.read_sql('SELECT * FROM [CF - regionalized - AcidTerr - native (NOx emissions to air)]',
                               self.conn).set_index('cell'),
            'SO2': pd.read_sql('SELECT * FROM [CF - regionalized - AcidTerr - native (SO2 emissions to air)]',
                               self.conn).set_index('cell')}

        inter_country = pd.read_sql('SELECT * FROM [SI - Acidification/Eutrophication - Countries cell resolution]',
                                    self.conn).set_index('cell').T
        inter_continent = pd.read_sql('SELECT * FROM [SI - Acidification/Eutrophication - Continents cell resolution]',
                                      self.conn).set_index('cell').T

        # ------------------------- CALCULATING MIDPOINT AND DAMAGE CFS --------------------------------
        country_emissions = {}
        continent_emissions = {}
        for substance in cell_emissions:
            cell_emissions[substance].index.name = None
            country_emissions[substance] = inter_country.dot(cell_emissions[substance].Emission)
            continent_emissions[substance] = inter_continent.dot(cell_emissions[substance].Emission)

        cfs_midpoint = pd.DataFrame(None, cell_emissions.keys(), country_emissions['SO2'].index)

        for substance in cell_emissions:
            # for countries
            for country in country_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'midpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     country_emissions[substance].loc[country] *
                                     inter_country.loc[country]).sum()
                    cfs_midpoint.loc[substance, country] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # for continents
            for continent in continent_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'midpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     continent_emissions[substance].loc[continent] *
                                     inter_continent.loc[continent]).sum()
                    cfs_midpoint.loc[substance, continent] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # global
            GLO_value = (cell_emissions[substance].loc[:, 'midpoint'] *
                         cell_emissions[substance].loc[:, 'Emission'] /
                         cell_emissions[substance].Emission.sum()).sum()
            cfs_midpoint.loc[substance, 'GLO'] = GLO_value

        # divide by reference flow -> SO2
        cfs_midpoint /= cfs_midpoint.loc['SO2', 'GLO']

        cfs_damage = pd.DataFrame(None, cell_emissions.keys(), country_emissions['SO2'].index)

        for substance in cell_emissions:
            # for countries
            for country in country_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'endpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     country_emissions[substance].loc[country] *
                                     inter_country.loc[country]).sum()
                    cfs_damage.loc[substance, country] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # for continents
            for continent in continent_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'endpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     continent_emissions[substance].loc[continent] *
                                     inter_continent.loc[continent]).sum()
                    cfs_damage.loc[substance, continent] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # global
            GLO_value = (cell_emissions[substance].loc[:, 'endpoint'] *
                         cell_emissions[substance].loc[:, 'Emission'] /
                         cell_emissions[substance].Emission.sum()).sum()
            cfs_damage.loc[substance, 'GLO'] = GLO_value

        # ------------------------------ DATA FORMATTING -----------------------------------
        # convert country names to ISO 2 letter codes
        # lines of logger to avoid having warnings showing up
        coco_logger = coco.logging.getLogger()
        coco_logger.setLevel(coco.logging.CRITICAL)
        cfs_midpoint.columns = coco.convert(cfs_midpoint.iloc[:, :-7], to='ISO2') + ['RNA', 'RLA', 'RER', 'RAS', 'RAF',
                                                                                     'OCE', 'GLO']
        cfs_damage.columns = cfs_midpoint.columns
        # countries that don't exist anymore or that are not part of the UN are dropped
        cfs_midpoint.drop('not found', axis=1, inplace=True)
        cfs_damage.drop('not found', axis=1, inplace=True)

        cfs_midpoint = cfs_midpoint.stack().reset_index()
        cfs_midpoint.columns = ['Elem flow', 'Region code', 'CF value']
        cfs_damage = cfs_damage.stack().reset_index()
        cfs_damage.columns = ['Elem flow', 'Region code', 'CF value']
        cfs_midpoint.loc[:, 'CF unit'] = 'kg SO2 eq'
        cfs_midpoint.loc[:, 'MP or Damage'] = 'Midpoint'
        cfs_damage.loc[:, 'CF unit'] = 'PDF.m2.yr'
        cfs_damage.loc[:, 'MP or Damage'] = 'Damage'
        cfs = pd.concat([cfs_midpoint, cfs_damage])
        cfs.loc[:, 'Impact category'] = 'Terrestrial acidification'
        cfs.loc[:, 'Compartment'] = 'Air'
        cfs.loc[:, 'Sub-compartment'] = '(unspecified)'
        cfs.loc[:, 'Elem flow unit'] = 'kg'
        cfs.loc[:, 'Native geographical resolution scale'] = 'Country'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NH3', 'CAS number'] = '007664-41-7'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NOx', 'CAS number'] = '011104-93-1'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'SO2', 'CAS number'] = '007446-09-5'
        cfs = clean_up_dataframe(cfs)
        cfs.loc[[i for i in cfs.index if cfs.loc[i, 'Region code'] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF', 'OCE', 'GLO']],
                'Native geographical resolution scale'] = 'Continent'
        cfs.loc[cfs.loc[:, 'Region code'] == 'GLO', 'Native geographical resolution scale'] = 'Global'

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - stoechiometry]', self.conn)
        for i in stoc.index:
            proxy = stoc.loc[i, 'Corresponding elem flow in IW']
            df = cfs[cfs.loc[:, 'Elem flow'] == proxy].copy('deep')
            if not df.empty:
                df.loc[:, 'Elem flow'] = stoc.loc[i, 'Elem flow name in Simapro']
                df.loc[:, 'CAS number'] = stoc.loc[i, 'CAS number']
                df.loc[:, 'CF value'] *= stoc.loc[i, 'Stoechiometric ratio']

                cfs = pd.concat([cfs, df])

        # ------------------------------------ FINAL FORMATTING ---------------------------------
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NH3', 'Elem flow'] = 'Ammonia'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NOx', 'Elem flow'] = 'Nitrogen oxides'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'SO2', 'Elem flow'] = 'Sulfur dioxide'
        # concat elem flow and region code in the same name
        cfs.loc[:, 'Elem flow name'] = [', '.join(i) for i in list(zip(cfs.loc[:, 'Elem flow'], cfs.loc[:, 'Region code']))]
        cfs = cfs.drop(['Elem flow', 'Region code'], axis=1)
        cfs = clean_up_dataframe(cfs)

        # concat with master_db
        self.master_db = pd.concat([self.master_db, cfs])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_marine_eutrophication_cfs(self):
        """
        Loading the CFs for the marine eutrophication impact category. This includes CFs coming from the
        original article of IW+, as well as other CFs extrapolated from stoechiometric ratios.

        Concerned impact categories:
            - Marine eutrophication

        :return: updated master_db
        """

        # ------------------------------ LOADING DATA -----------------------------------
        cell_emissions = {
            'NH3': pd.read_sql('SELECT * FROM [CF - regionalized - MarEutro - native (NH3 emissions to air)]',
                               self.conn).set_index('cell'),
            'NOx': pd.read_sql('SELECT * FROM [CF - regionalized - MarEutro - native (NOx emissions to air)]',
                               self.conn).set_index('cell'),
            'HNO3': pd.read_sql('SELECT * FROM [CF - regionalized - MarEutro - native (HNO3 emissions to air)]',
                               self.conn).set_index('cell')}

        inter_country = pd.read_sql('SELECT * FROM [SI - Acidification/Eutrophication - Countries cell resolution]',
                                    self.conn).set_index('cell').T
        inter_continent = pd.read_sql('SELECT * FROM [SI - Acidification/Eutrophication - Continents cell resolution]',
                                      self.conn).set_index('cell').T

        # ------------------------- CALCULATING MIDPOINT AND DAMAGE CFS --------------------------------
        country_emissions = {}
        continent_emissions = {}
        for substance in cell_emissions:
            cell_emissions[substance].index.name = None
            country_emissions[substance] = inter_country.dot(cell_emissions[substance].Emission)
            continent_emissions[substance] = inter_continent.dot(cell_emissions[substance].Emission)

        cfs_midpoint = pd.DataFrame(None, cell_emissions.keys(), country_emissions['NOx'].index)

        for substance in cell_emissions:
            # for countries
            for country in country_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'midpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     country_emissions[substance].loc[country] *
                                     inter_country.loc[country]).sum()
                    cfs_midpoint.loc[substance, country] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # for continents
            for continent in continent_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'midpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     continent_emissions[substance].loc[continent] *
                                     inter_continent.loc[continent]).sum()
                    cfs_midpoint.loc[substance, continent] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # global
            GLO_value = (cell_emissions[substance].loc[:, 'midpoint'] *
                         cell_emissions[substance].loc[:, 'Emission'] /
                         cell_emissions[substance].Emission.sum()).sum()
            cfs_midpoint.loc[substance, 'GLO'] = GLO_value

        # factor 0.7 is the reference (it's nitrogen GLO in water)
        cfs_midpoint /= 0.7

        cfs_damage = pd.DataFrame(None, cell_emissions.keys(), country_emissions['NOx'].index)

        for substance in cell_emissions:
            # for countries
            for country in country_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'endpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     country_emissions[substance].loc[country] *
                                     inter_country.loc[country]).sum()
                    cfs_damage.loc[substance, country] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # for continents
            for continent in continent_emissions[substance].index:
                try:
                    aggregated_cf = (cell_emissions[substance].loc[:, 'endpoint'] *
                                     cell_emissions[substance].loc[:, 'Emission'] /
                                     continent_emissions[substance].loc[continent] *
                                     inter_continent.loc[continent]).sum()
                    cfs_damage.loc[substance, continent] = aggregated_cf
                except ZeroDivisionError:
                    pass

            # global
            GLO_value = (cell_emissions[substance].loc[:, 'endpoint'] *
                         cell_emissions[substance].loc[:, 'Emission'] /
                         cell_emissions[substance].Emission.sum()).sum()
            cfs_damage.loc[substance, 'GLO'] = GLO_value

        # ------------------------------ DATA FORMATTING -----------------------------------
        # convert country names to ISO 2 letter codes
        # lines of logger to avoid having warnings showing up
        coco_logger = coco.logging.getLogger()
        coco_logger.setLevel(coco.logging.CRITICAL)
        cfs_midpoint.columns = coco.convert(cfs_midpoint.iloc[:, :-7], to='ISO2') + ['RNA', 'RLA', 'RER', 'RAS', 'RAF',
                                                                                     'OCE', 'GLO']
        cfs_damage.columns = cfs_midpoint.columns
        # countries that don't exist anymore or that are not part of the UN are dropped
        cfs_midpoint.drop('not found', axis=1, inplace=True)
        cfs_damage.drop('not found', axis=1, inplace=True)

        cfs_midpoint = cfs_midpoint.stack().reset_index()
        cfs_midpoint.columns = ['Elem flow', 'Region code', 'CF value']
        cfs_damage = cfs_damage.stack().reset_index()
        cfs_damage.columns = ['Elem flow', 'Region code', 'CF value']
        cfs_midpoint.loc[:, 'CF unit'] = 'kg N N-lim eq'
        cfs_midpoint.loc[:, 'MP or Damage'] = 'Midpoint'
        cfs_damage.loc[:, 'CF unit'] = 'PDF.m2.yr'
        cfs_damage.loc[:, 'MP or Damage'] = 'Damage'
        cfs = pd.concat([cfs_midpoint, cfs_damage])
        cfs.loc[:, 'Impact category'] = 'Marine eutrophication'
        cfs.loc[:, 'Compartment'] = 'Air'
        cfs.loc[:, 'Sub-compartment'] = '(unspecified)'
        cfs.loc[:, 'Elem flow unit'] = 'kg'
        cfs.loc[:, 'Native geographical resolution scale'] = 'Country'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NH3', 'CAS number'] = '007664-41-7'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NOx', 'CAS number'] = '011104-93-1'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'HNO3', 'CAS number'] = '007697-37-2'
        cfs = clean_up_dataframe(cfs)
        cfs.loc[[i for i in cfs.index if cfs.loc[i, 'Region code'] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF', 'OCE', 'GLO']],
                'Native geographical resolution scale'] = 'Continent'
        cfs.loc[cfs.loc[:, 'Region code'] == 'GLO', 'Native geographical resolution scale'] = 'Global'

        # add non-regionalized flows (water emissions)
        cfs = pd.concat(
            [cfs, pd.read_sql('SELECT * FROM [CF - not regionalized - MarEutro]', self.conn)])

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - stoechiometry]', self.conn)
        for i in stoc.index:
            if stoc.loc[i,'Formula'] not in list(cfs.loc[:,'Elem flow']):
                proxy = stoc.loc[i, 'Corresponding elem flow in IW']
                df = cfs[cfs.loc[:, 'Elem flow'] == proxy].copy('deep')
                if not df.empty:
                    df.loc[:, 'Elem flow'] = stoc.loc[i, 'Elem flow name in Simapro']
                    df.loc[:, 'CAS number'] = stoc.loc[i, 'CAS number']
                    df.loc[:, 'CF value'] *= stoc.loc[i, 'Stoechiometric ratio']

                    cfs = pd.concat([cfs, df])

        # ------------------------------------ FINAL FORMATTING ---------------------------------
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NH3', 'Elem flow'] = 'Ammonia'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NOx', 'Elem flow'] = 'Nitrogen oxides'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'HNO3', 'Elem flow'] = 'Nitric acid'
        # concat elem flow and region code in the same name
        cfs.loc[:, 'Elem flow name'] = [', '.join(i) if type(i[1]) == str else i[0]
                                        for i in list(zip(cfs.loc[:, 'Elem flow'], cfs.loc[:, 'Region code']))]
        cfs = cfs.drop(['Elem flow', 'Region code'], axis=1)
        cfs = clean_up_dataframe(cfs)

        # concat with master_db
        self.master_db = pd.concat([self.master_db, cfs])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_freshwater_eutrophication_cfs(self):
        """
        Loading the CFs for the freshwater eutrophication impact category. This includes CFs coming from the
        original article of IW+, as well as other CFs extrapolated from stoechiometric ratios.

        Concerned impact categories:
            - Freshwater eutrophication

        :return: updated master_db
        """

        # ------------------------------ LOADING DATA -----------------------------------
        cell_emissions = pd.read_sql(sql='SELECT * FROM [CF - regionalized - FWEutro - native]', con=self.conn)
        inter_country = pd.read_sql(
            sql='SELECT * FROM [SI - Freshwater Eutrophication - Countries/continents cell resolution]',
            con=self.conn)

        df = cell_emissions.merge(inter_country, left_on='cell', right_on='nindex').loc[
             :, ['cell', 'FF_P_yr', 'Population_per_cell', 'ISO_2DIGIT', 'CNTRY_NAME', 'CONT_NAME']]

        # 'NA' is converted to None by pandas. So we force it back to 'NA' ('NA' is for Namibia)
        df.ISO_2DIGIT = df.ISO_2DIGIT.fillna('NA')

        # ------------------------- CALCULATING MIDPOINT AND DAMAGE CFS --------------------------------
        countries = list(set(df.ISO_2DIGIT))
        countries.sort()
        cfs_midpoint = pd.DataFrame(None, countries, ['CF'])
        for country in countries:
            ponderation = df[df.ISO_2DIGIT == country].Population_per_cell / df[
                df.ISO_2DIGIT == country].Population_per_cell.sum()
            cfs_midpoint.loc[country, 'CF'] = (ponderation * df.loc[ponderation.index, 'FF_P_yr']).sum()

        # global value
        cfs_midpoint.loc['GLO', 'CF'] = (df.Population_per_cell / df.Population_per_cell.sum()).dot(
            df.loc[:, 'FF_P_yr'])

        continents = set(df.CONT_NAME)
        for continent in continents:
            ponderation = df[df.CONT_NAME == continent].Population_per_cell / df[
                df.CONT_NAME == continent].Population_per_cell.sum()
            cfs_midpoint.loc[continent, 'CF'] = (ponderation * df.loc[ponderation.index, 'FF_P_yr']).sum()

        # Damage: just multiply the value by 11.4
        cfs_damage = cfs_midpoint.copy('deep')
        cfs_damage *= 11.4

        # Ponderate by reference flow (Phosphate, GLO)
        cfs_midpoint /= cfs_midpoint.loc['GLO', 'CF']

        # ------------------------------ DATA FORMATTING -----------------------------------
        cfs_midpoint.loc[:, 'CF unit'] = 'kg PO4 P-lim eq'
        cfs_midpoint.loc[:, 'MP or Damage'] = 'Midpoint'
        cfs_damage.loc[:, 'CF unit'] = 'PDF.m2.yr'
        cfs_damage.loc[:, 'MP or Damage'] = 'Damage'
        cfs = pd.concat([cfs_midpoint, cfs_damage])
        cfs.loc[:, 'Impact category'] = 'Freshwater eutrophication'
        cfs.loc[:, 'Compartment'] = 'Water'
        cfs.loc[:, 'Sub-compartment'] = '(unspecified)'
        cfs.loc[:, 'Elem flow unit'] = 'kg'
        cfs.loc[:, 'Native geographical resolution scale'] = 'Country'
        cfs.loc[['RNA', 'RLA', 'RER', 'RAS', 'RAF', 'OCE'], 'Native geographical resolution scale'] = 'Continent'
        cfs.loc['GLO', 'Native geographical resolution scale'] = 'Global'
        cfs.loc[:, 'Elem flow'] = 'Phosphate'
        cfs.loc[:, 'CAS number'] = '014265-44-2'
        cfs = cfs.reset_index()
        cfs.columns = ['Region code', 'CF value', 'CF unit', 'MP or Damage', 'Impact category', 'Compartment',
                       'Sub-compartment', 'Elem flow unit', 'Native geographical resolution scale',
                       'Elem flow', 'CAS number']

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - stoechiometry]', self.conn)
        for i in stoc.index:
            proxy = stoc.loc[i, 'Corresponding elem flow in IW']
            df = cfs[cfs.loc[:, 'Elem flow'] == proxy].copy('deep')
            if not df.empty:
                df.loc[:, 'Elem flow'] = stoc.loc[i, 'Elem flow name in Simapro']
                df.loc[:, 'CAS number'] = stoc.loc[i, 'CAS number']
                df.loc[:, 'CF value'] *= stoc.loc[i, 'Stoechiometric ratio']

                cfs = pd.concat([cfs, df])

        # ------------------------------------ FINAL FORMATTING ---------------------------------
        cfs.loc[:, 'Elem flow name'] = [', '.join(i) for i in
                                        list(zip(cfs.loc[:, 'Elem flow'], cfs.loc[:, 'Region code']))]
        cfs = cfs.drop(['Elem flow', 'Region code'], axis=1)
        cfs = clean_up_dataframe(cfs)

        # concat with master_db
        self.master_db = pd.concat([self.master_db, cfs])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_land_use_cfs(self):
        """
        Loading the CFs for the land use impact categories.

        Concerned impact categories:
            - Land occupation, biodiversity
            - Land transformation, biodiversity

        :return: updated master_db
        """

        # ------------------------------ LOADING DATA -----------------------------------
        original_cf_occup = pd.read_sql(sql='SELECT * FROM [CF - regionalized - land occupation per biome]',
                                        con=self.conn)

        intersect_country = pd.read_sql(sql='SELECT * FROM [SI - Land occupation - countries cell resolution]',
                                        con=self.conn)

        land_use_type = pd.read_sql(sql='SELECT * FROM [SI - Land occupation - land use type per country/continent]',
                                    con=self.conn)

        recovery_times = pd.read_sql(sql='SELECT * FROM [SI - land transformation - recovery times]', con=self.conn)

        original_cf_occup = original_cf_occup.merge(intersect_country, left_on='Biome', right_on='WWF_MHTNUM',
                                                    how='outer')

        original_cf_occup = original_cf_occup.merge(land_use_type, on='ISO', how='outer')

        original_cf_occup = original_cf_occup.merge(recovery_times, left_on='WWF_MHTNUM', right_on='Biome_type')

        original_cf_occup = original_cf_occup.set_index(['ECO_CODE', 'ISO', 'Distribution'])
        original_cf_occup.index.names = None, None, None

        original_cf_occup = original_cf_occup.loc[:,
                            ['Forest/Grassland, not used', 'Secondary Vegetation', 'Forest, used',
                             'Pasture/meadow', 'Annual crops', 'Permanent crops',
                             'Agriculture, mosaic (Agroforestry)', 'Artificial areas', 'Proxy_area_AGRI',
                             'Proxy_area_PAST', 'Proxy_area_URB', 'Proxy_area_MANAGED_FOR', 'Proxy_area_FOR',
                             'ID_CONT_1', 'Time']]
        original_cf_occup = original_cf_occup.dropna(
            subset=['Forest/Grassland, not used', 'Secondary Vegetation', 'Forest, used',
                    'Pasture/meadow', 'Annual crops', 'Permanent crops',
                    'Agriculture, mosaic (Agroforestry)', 'Artificial areas'], how='all')
        original_cf_occup = original_cf_occup.drop(original_cf_occup[original_cf_occup.ID_CONT_1.isna()].index)

        original_cf_transfo = original_cf_occup.copy()

        for land_type in ['Forest/Grassland, not used', 'Secondary Vegetation', 'Forest, used', 'Pasture/meadow',
                          'Annual crops',
                          'Permanent crops', 'Agriculture, mosaic (Agroforestry)', 'Artificial areas']:
            original_cf_transfo.loc[:, land_type] *= original_cf_transfo.loc[:, 'Time'] * 0.5

        # ------------------------- CALCULATING MIDPOINT AND DAMAGE CFS --------------------------------
        # the proxy used for ponderation
        proxy = {'Forest/Grassland, not used': ('Proxy_area_FOR', 'surface_COUNT_FOR'),
                 'Secondary Vegetation': ('Proxy_area_MANAGED_FOR', 'surface_COUNT_FOR'),
                 'Forest, used': ('Proxy_area_MANAGED_FOR', 'surface_COUNT_FOR'),
                 'Pasture/meadow': ('Proxy_area_PAST', 'surface_COUNT_PAST'),
                 'Annual crops': ('Proxy_area_AGRI', 'surface_COUNT_AGRI'),
                 'Permanent crops': ('Proxy_area_AGRI', 'surface_COUNT_AGRI'),
                 'Agriculture, mosaic (Agroforestry)': ('Proxy_area_AGRI', 'surface_COUNT_AGRI'),
                 'Artificial areas': ('Proxy_area_URB', 'surface_COUNT_URB')}

        #                             - FOR LAND OCCUPATION -

        # for countries
        index = list(set([i[1] for i in original_cf_occup.index]))
        index.sort()
        cfs = pd.DataFrame(0, index, proxy)

        for cty in cfs.index:
            for land_type in cfs.columns:
                native_cfs = original_cf_occup.loc(axis=0)[:, cty, 'Median'].loc[:, land_type].dropna()

                ponderation = (original_cf_occup.loc[native_cfs.index, proxy[land_type][0]] /
                               original_cf_occup.loc[native_cfs.index, proxy[land_type][0]].sum())
                df = (native_cfs * ponderation).sum()

                cfs.loc[cty, land_type] = df

        # for continents
        for cont in set(original_cf_occup.ID_CONT_1):
            for land_type in cfs.columns:
                native_cfs = original_cf_occup[original_cf_occup.ID_CONT_1 == cont].loc(axis=0)[:, :, 'Median'].loc[:,
                             land_type].dropna()

                ponderation = (original_cf_occup.loc[native_cfs.index, proxy[land_type][0]] /
                               original_cf_occup.loc[native_cfs.index, proxy[land_type][0]].sum())
                df = (native_cfs * ponderation).sum()

                cfs.loc[cont, land_type] = df

        # for global value
        for land_type in cfs.columns:
            native_cfs = original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:, land_type].dropna()

            ponderation = (original_cf_occup.loc[native_cfs.index, proxy[land_type][0]] /
                           original_cf_occup.loc[native_cfs.index, proxy[land_type][0]].sum())
            cf = (native_cfs * ponderation).sum()

            cfs.loc['GLO', land_type] = cf

        # CFs obtained are already for damage categories
        occupation_damage_cfs = cfs.copy('deep')

        # to get midpoint CFS, simply normalize with reference
        occupation_midpoint_cfs = cfs / cfs.loc['GLO', 'Annual crops']

        # to get an unspecified value
        for country in set([i[1] for i in original_cf_occup.index]):

            CF_sum = 0

            for land_use in proxy.keys():
                if proxy[land_use][0] == 'Proxy_area_AGRI':
                    # gotta split agriculture area evenly between the three agriculture land types
                    CF_sum += (original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:, land_use] *
                               original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:, proxy[land_use][0]]).sum() / 3
                if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                    # gotta split forest area evenly between the two forest land types
                    CF_sum += (original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:, land_use] *
                               original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:, proxy[land_use][0]]).sum() / 2
                else:
                    CF_sum += (original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:, land_use] *
                               original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:, proxy[land_use][0]]).sum()

            total_area_covered = 0

            for land_use in proxy.keys():
                covered_or_not = original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:, proxy.keys()].mask(
                    original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:, proxy.keys()] > -10, 1)
                if proxy[land_use][0] == 'Proxy_area_AGRI':
                    # gotta split agriculture area evenly between the three agriculture land types
                    total_area_covered += (covered_or_not.loc[:, land_use] *
                                           original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:,
                                           proxy[land_use][0]]).sum() / 3
                if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                    # gotta split forest area evenly between the two forest land types
                    total_area_covered += (covered_or_not.loc[:, land_use] *
                                           original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:,
                                           proxy[land_use][0]]).sum() / 2
                else:
                    total_area_covered += (covered_or_not.loc[:, land_use] *
                                           original_cf_occup.loc(axis=0)[:, country, 'Median'].loc[:,
                                           proxy[land_use][0]]).sum()

            occupation_damage_cfs.loc[country, 'Unspecified'] = CF_sum / total_area_covered

        conts = ['OCE', 'RAF', 'RAS', 'RER', 'RLA', 'RNA', 'Antarctica']

        for cont in conts:

            cont_data = original_cf_occup.loc(axis=0)[:, :, 'Median'][
                original_cf_occup.loc(axis=0)[:, :, 'Median'].ID_CONT_1 == cont]

            CF_sum = 0

            for land_use in proxy.keys():
                if proxy[land_use][0] == 'Proxy_area_AGRI':
                    # gotta split agriculture area evenly between the three agriculture land types
                    CF_sum += (cont_data.loc[:, land_use] * cont_data.loc[:, proxy[land_use][0]]).sum() / 3
                if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                    # gotta split forest area evenly between the two forest land types
                    CF_sum += (cont_data.loc[:, land_use] * cont_data.loc[:, proxy[land_use][0]]).sum() / 2
                else:
                    CF_sum += (cont_data.loc[:, land_use] * cont_data.loc[:, proxy[land_use][0]]).sum()

            total_area_covered = 0

            for land_use in proxy.keys():
                covered_or_not = cont_data.loc[:, proxy.keys()].mask(cont_data.loc[:, proxy.keys()] > -10, 1)
                if proxy[land_use][0] == 'Proxy_area_AGRI':
                    # gotta split agriculture area evenly between the three agriculture land types
                    total_area_covered += (covered_or_not.loc[:, land_use] * cont_data.loc[:,
                                                                             proxy[land_use][0]]).sum() / 3
                if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                    # gotta split forest area evenly between the two forest land types
                    total_area_covered += (covered_or_not.loc[:, land_use] * cont_data.loc[:,
                                                                             proxy[land_use][0]]).sum() / 2
                else:
                    total_area_covered += (covered_or_not.loc[:, land_use] * cont_data.loc[:, proxy[land_use][0]]).sum()

            occupation_damage_cfs.loc[cont, 'Unspecified'] = CF_sum / total_area_covered

        # global
        CF_sum = 0

        for land_use in proxy.keys():
            if proxy[land_use][0] == 'Proxy_area_AGRI':
                # gotta split agriculture area evenly between the three agriculture land types
                CF_sum += (original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:, land_use] *
                           original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:, proxy[land_use][0]]).sum() / 3
            if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                # gotta split forest area evenly between the two forest land types
                CF_sum += (original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:, land_use] *
                           original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:, proxy[land_use][0]]).sum() / 2
            else:
                CF_sum += (original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:, land_use] *
                           original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:, proxy[land_use][0]]).sum()

        total_area_covered = 0

        for land_use in proxy.keys():
            covered_or_not = original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:, proxy.keys()].mask(
                original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:, proxy.keys()] > -10, 1)
            if proxy[land_use][0] == 'Proxy_area_AGRI':
                # gotta split agriculture area evenly between the three agriculture land types
                total_area_covered += (covered_or_not.loc[:, land_use] *
                                       original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:,
                                       proxy[land_use][0]]).sum() / 3
            if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                # gotta split forest area evenly between the two forest land types
                total_area_covered += (covered_or_not.loc[:, land_use] *
                                       original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:,
                                       proxy[land_use][0]]).sum() / 2
            else:
                total_area_covered += (covered_or_not.loc[:, land_use] *
                                       original_cf_occup.loc(axis=0)[:, :, 'Median'].loc[:, proxy[land_use][0]]).sum()

        occupation_damage_cfs.loc['GLO', 'Unspecified'] = CF_sum / total_area_covered

        occupation_midpoint_cfs.loc[:, 'Unspecified'] = (occupation_damage_cfs.loc[:, 'Unspecified'] /
                                                         occupation_damage_cfs.loc['GLO', 'Annual crops'])

        #                              - FOR LAND TRANSFORMATION -

        # for countries
        index = list(set([i[1] for i in original_cf_transfo.index]))
        index.sort()
        cfs = pd.DataFrame(0, index, proxy)

        for cty in cfs.index:
            for land_type in cfs.columns:
                native_cfs = original_cf_transfo.loc(axis=0)[:, cty, 'Median'].loc[:, land_type].dropna()

                ponderation = (original_cf_transfo.loc[native_cfs.index, proxy[land_type][0]] /
                               original_cf_transfo.loc[native_cfs.index, proxy[land_type][0]].sum())
                cf = (native_cfs * ponderation).sum()

                cfs.loc[cty, land_type] = cf

        # for continents
        for cont in set(original_cf_transfo.ID_CONT_1):
            for land_type in cfs.columns:
                native_cfs = original_cf_transfo[original_cf_transfo.ID_CONT_1 == cont].loc(axis=0)[:, :, 'Median'].loc[
                             :, land_type].dropna()

                ponderation = (original_cf_transfo.loc[native_cfs.index, proxy[land_type][0]] /
                               original_cf_transfo.loc[native_cfs.index, proxy[land_type][0]].sum())
                cf = (native_cfs * ponderation).sum()

                cfs.loc[cont, land_type] = cf

        for land_type in cfs.columns:
            native_cfs = original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:, land_type].dropna()

            ponderation = (original_cf_transfo.loc[native_cfs.index, proxy[land_type][0]] /
                           original_cf_transfo.loc[native_cfs.index, proxy[land_type][0]].sum())
            cf = (native_cfs * ponderation).sum()

            cfs.loc['GLO', land_type] = cf

        # CFs obtained are already for damage categories
        transformation_damage_cfs = cfs.copy('deep')

        # to get midpoint CFS, simply normalize with reference
        transformation_midpoint_cfs = cfs / cfs.loc['GLO', 'Annual crops']

        # to egt an unspecified value
        for country in set([i[1] for i in original_cf_transfo.index]):

            CF_sum = 0

            for land_use in proxy.keys():
                if proxy[land_use][0] == 'Proxy_area_AGRI':
                    # gotta split agriculture area evenly between the three agriculture land types
                    CF_sum += (original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:, land_use] *
                               original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:,
                               proxy[land_use][0]]).sum() / 3
                if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                    # gotta split forest area evenly between the two forest land types
                    CF_sum += (original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:, land_use] *
                               original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:,
                               proxy[land_use][0]]).sum() / 2
                else:
                    CF_sum += (original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:, land_use] *
                               original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:, proxy[land_use][0]]).sum()

            total_area_covered = 0

            for land_use in proxy.keys():
                covered_or_not = original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:, proxy.keys()].mask(
                    original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:, proxy.keys()] > -10, 1)
                if proxy[land_use][0] == 'Proxy_area_AGRI':
                    # gotta split agriculture area evenly between the three agriculture land types
                    total_area_covered += (covered_or_not.loc[:, land_use] *
                                           original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:,
                                           proxy[land_use][0]]).sum() / 3
                if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                    # gotta split forest area evenly between the two forest land types
                    total_area_covered += (covered_or_not.loc[:, land_use] *
                                           original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:,
                                           proxy[land_use][0]]).sum() / 2
                else:
                    total_area_covered += (covered_or_not.loc[:, land_use] *
                                           original_cf_transfo.loc(axis=0)[:, country, 'Median'].loc[:,
                                           proxy[land_use][0]]).sum()

            transformation_damage_cfs.loc[country, 'Unspecified'] = CF_sum / total_area_covered

        conts = ['OCE', 'RAF', 'RAS', 'RER', 'RLA', 'RNA', 'Antarctica']

        for cont in conts:

            cont_data = original_cf_transfo.loc(axis=0)[:, :, 'Median'][
                original_cf_transfo.loc(axis=0)[:, :, 'Median'].ID_CONT_1 == cont]

            CF_sum = 0

            for land_use in proxy.keys():
                if proxy[land_use][0] == 'Proxy_area_AGRI':
                    # gotta split agriculture area evenly between the three agriculture land types
                    CF_sum += (cont_data.loc[:, land_use] * cont_data.loc[:, proxy[land_use][0]]).sum() / 3
                if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                    # gotta split forest area evenly between the two forest land types
                    CF_sum += (cont_data.loc[:, land_use] * cont_data.loc[:, proxy[land_use][0]]).sum() / 2
                else:
                    CF_sum += (cont_data.loc[:, land_use] * cont_data.loc[:, proxy[land_use][0]]).sum()

            total_area_covered = 0

            for land_use in proxy.keys():
                covered_or_not = cont_data.loc[:, proxy.keys()].mask(cont_data.loc[:, proxy.keys()] > -10, 1)
                if proxy[land_use][0] == 'Proxy_area_AGRI':
                    # gotta split agriculture area evenly between the three agriculture land types
                    total_area_covered += (covered_or_not.loc[:, land_use] * cont_data.loc[:,
                                                                             proxy[land_use][0]]).sum() / 3
                if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                    # gotta split forest area evenly between the two forest land types
                    total_area_covered += (covered_or_not.loc[:, land_use] * cont_data.loc[:,
                                                                             proxy[land_use][0]]).sum() / 2
                else:
                    total_area_covered += (covered_or_not.loc[:, land_use] * cont_data.loc[:, proxy[land_use][0]]).sum()

            transformation_damage_cfs.loc[cont, 'Unspecified'] = CF_sum / total_area_covered

        # global
        CF_sum = 0

        for land_use in proxy.keys():
            if proxy[land_use][0] == 'Proxy_area_AGRI':
                # gotta split agriculture area evenly between the three agriculture land types
                CF_sum += (original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:, land_use] *
                           original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:, proxy[land_use][0]]).sum() / 3
            if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                # gotta split forest area evenly between the two forest land types
                CF_sum += (original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:, land_use] *
                           original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:, proxy[land_use][0]]).sum() / 2
            else:
                CF_sum += (original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:, land_use] *
                           original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:, proxy[land_use][0]]).sum()

        total_area_covered = 0

        for land_use in proxy.keys():
            covered_or_not = original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:, proxy.keys()].mask(
                original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:, proxy.keys()] > -10, 1)
            if proxy[land_use][0] == 'Proxy_area_AGRI':
                # gotta split agriculture area evenly between the three agriculture land types
                total_area_covered += (covered_or_not.loc[:, land_use] *
                                       original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:,
                                       proxy[land_use][0]]).sum() / 3
            if proxy[land_use][0] == 'Proxy_area_MANAGED_FOR':
                # gotta split forest area evenly between the two forest land types
                total_area_covered += (covered_or_not.loc[:, land_use] *
                                       original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:,
                                       proxy[land_use][0]]).sum() / 2
            else:
                total_area_covered += (covered_or_not.loc[:, land_use] *
                                       original_cf_transfo.loc(axis=0)[:, :, 'Median'].loc[:, proxy[land_use][0]]).sum()

        transformation_damage_cfs.loc['GLO', 'Unspecified'] = CF_sum / total_area_covered

        transformation_midpoint_cfs.loc[:, 'Unspecified'] = (transformation_damage_cfs.loc[:, 'Unspecified'] /
                                                             transformation_damage_cfs.loc['GLO', 'Annual crops'])

        # ------------------------------ DATA FORMATTING -----------------------------------
        occupation_damage_cfs = occupation_damage_cfs.stack().reset_index()
        occupation_midpoint_cfs = occupation_midpoint_cfs.stack().reset_index()
        transformation_damage_cfs = transformation_damage_cfs.stack().reset_index()
        transformation_midpoint_cfs = transformation_midpoint_cfs.stack().reset_index()
        occupation_damage_cfs.loc[:, 'Impact category'] = 'Land occupation, biodiversity'
        occupation_midpoint_cfs.loc[:, 'Impact category'] = 'Land occupation, biodiversity'
        transformation_damage_cfs.loc[:, 'Impact category'] = 'Land transformation, biodiversity'
        transformation_midpoint_cfs.loc[:, 'Impact category'] = 'Land transformation, biodiversity'
        occupation_damage_cfs.loc[:, 'MP or Damage'] = 'Damage'
        occupation_midpoint_cfs.loc[:, 'MP or Damage'] = 'Midpoint'
        transformation_damage_cfs.loc[:, 'MP or Damage'] = 'Damage'
        transformation_midpoint_cfs.loc[:, 'MP or Damage'] = 'Midpoint'
        occupation_damage_cfs.loc[:, 'CF unit'] = 'PDF.m2.yr'
        occupation_midpoint_cfs.loc[:, 'CF unit'] = 'm2 arable land eq .yr'
        transformation_damage_cfs.loc[:, 'CF unit'] = 'PDF.m2.yr'
        transformation_midpoint_cfs.loc[:, 'CF unit'] = 'm2 arable land eq'
        occupation_damage_cfs.loc[:, 'Elem flow unit'] = 'm2.yr'
        occupation_midpoint_cfs.loc[:, 'Elem flow unit'] = 'm2.yr'
        transformation_damage_cfs.loc[:, 'Elem flow unit'] = 'm2'
        transformation_midpoint_cfs.loc[:, 'Elem flow unit'] = 'm2'
        land_cfs = pd.concat(
            [occupation_damage_cfs, occupation_midpoint_cfs, transformation_damage_cfs, transformation_midpoint_cfs])
        land_cfs = clean_up_dataframe(land_cfs)
        continents = [i for i in land_cfs.index if
                      land_cfs.loc[i, 'level_0'] in ['RNA', 'RER', 'RLA', 'RAF', 'OCE', 'RAS', 'GLO']]

        continental_data = land_cfs.loc[continents].copy('deep')
        land_cfs = land_cfs.drop(continents)
        land_cfs.loc[:, 'level_0'] = coco.convert(land_cfs.loc[:, 'level_0'], to='ISO2')
        land_cfs = land_cfs.drop([i for i in land_cfs.index if land_cfs.loc[i, 'level_0'] == 'not found'])
        land_cfs = pd.concat([land_cfs, continental_data])
        land_cfs.level_1 = [i.lower() for i in land_cfs.level_1]
        land_cfs.loc[:, 'Elem flow name'] = [', '.join(i) for i in
                                             list(zip(land_cfs.loc[:, 'level_1'], land_cfs.loc[:, 'level_0']))]
        land_cfs.loc[land_cfs.loc[:, 'Impact category'] == 'Land occupation, biodiversity', 'Elem flow name'] = [
            'Occupation, ' + i[0].lower() + i[1:] for i in land_cfs.loc[land_cfs.loc[:, 'Impact category'] ==
                                                                        'Land occupation, biodiversity', 'Elem flow name']]
        land_cfs.loc[land_cfs.loc[:, 'Impact category'] == 'Land transformation, biodiversity', 'Elem flow name'] = [
            'Transformation, to ' + i[0].lower() + i[1:] for i in land_cfs.loc[land_cfs.loc[:, 'Impact category'] ==
                                                                               'Land transformation, biodiversity', 'Elem flow name']]
        for_negative = land_cfs.loc[land_cfs.loc[:, 'Impact category'] == 'Land transformation, biodiversity'].copy(
            'deep')
        for_negative.loc[:, 0] *= -1
        for_negative.loc[:, 'Elem flow name'] = ['Transformation, from' + i.split('Transformation, to')[1] for i in
                                                 for_negative.loc[:, 'Elem flow name']]
        land_cfs = pd.concat([land_cfs, for_negative])
        land_cfs = clean_up_dataframe(land_cfs)
        land_cfs.loc[:, 'Native geographical resolution scale'] = 'Country'
        land_cfs.loc[
            [i for i in land_cfs.index if land_cfs.loc[i, 'level_0'] in ['RNA', 'RER', 'RLA', 'RAF', 'OCE', 'RAS']],
            'Native geographical resolution scale'] = 'Continent'
        land_cfs.loc[land_cfs.level_0 == 'GLO', 'Native geographical resolution scale'] = 'Global'
        land_cfs.loc[:, 'CAS number'] = ''
        land_cfs.loc[:, 'Compartment'] = 'Raw'
        land_cfs.loc[:, 'Sub-compartment'] = 'land'
        land_cfs = land_cfs.drop(['level_0', 'level_1'], axis=1)
        land_cfs = land_cfs.rename(columns={0: 'CF value'})

        # concat with master_db
        self.master_db = pd.concat([self.master_db, land_cfs])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_particulates_cfs(self):
        """
        Load CFs for particulate mater formation.

        Concerned impact categories:
            - Particulate matter formation

        :return: updated master_db
        """

        data = pd.read_sql(sql='SELECT * FROM [CF - regionalized - ParticulateMatter - native]', con=self.conn)
        conc = pd.read_sql(sql='SELECT * FROM [SI - mapping with regions of ecoinvent]', con=self.conn).set_index(
            'Ecoinvent_short_name').PM

        # first we need to aggregate sub-regions to country level for a few countries which otherwise undefined
        countries_of_subregions_to_aggregate = ['Gabon', 'Kenya', 'Somalia', 'Uganda', 'South Africa', 'Russia',
                                                'Norway', 'Spain', 'United Kingdom', 'Mexico', 'Saudi Arabia', 'Canada',
                                                'United States', 'Australia','Brazil', 'China', 'India', 'Indonesia']

        for country in countries_of_subregions_to_aggregate:
            CF_urban = ((data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                i, 'Native geographical resolution scale'] == 'city')], 'CF urban'] *
                         data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                             i, 'Native geographical resolution scale'] == 'city')], 'Population urban'] /
                         data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                             i, 'Native geographical resolution scale'] == 'city')], 'Population urban'].sum())).sum()
            CF_rural = ((data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                i, 'Native geographical resolution scale'] == 'sub-regions')], 'CF rural'] *
                         data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                             i, 'Native geographical resolution scale'] == 'sub-regions')], 'Population rural'] /
                         data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                             i, 'Native geographical resolution scale'] == 'sub-regions')], 'Population rural'].sum())).sum()
            CF_unspecified = ((data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                i, 'Native geographical resolution scale'] == 'sub-regions')], 'CF unspecified'] *
                               (data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                                   i, 'Native geographical resolution scale'] == 'sub-regions')], 'Population urban'] +
                                data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                                    i, 'Native geographical resolution scale'] == 'sub-regions')], 'Population rural']) /
                               (data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                                   i, 'Native geographical resolution scale'] == 'sub-regions')], 'Population urban'] +
                                data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                                    i, 'Native geographical resolution scale'] == 'sub-regions')], 'Population rural']).sum())).sum()
            urban_pop = data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                i, 'Native geographical resolution scale'] == 'city')], 'Population urban'].sum()
            rural_pop = data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                i, 'Native geographical resolution scale'] == 'sub-regions')], 'Population rural'].sum()

            data = pd.concat([data, pd.DataFrame(
                [country, country, urban_pop, rural_pop, CF_urban, CF_rural, CF_unspecified, 'sub-regions'],
                index=['Country-Region', 'Country', 'Population urban', 'Population rural', 'CF urban', 'CF rural',
                       'CF unspecified', 'Native geographical resolution scale']).T])
            data = clean_up_dataframe(data)

        # we then need to determine intake fractions for countries (these are required for secondary PM)
        for country in set(data.Country):
            data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                i, 'Native geographical resolution scale'] == 'sub-regions')], 'iF urban'] = (
                ((data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                    i, 'Native geographical resolution scale'] == 'city')], 'iF urban'] *
                  data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                      i, 'Native geographical resolution scale'] == 'city')], 'Population urban'] /
                  data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                      i, 'Native geographical resolution scale'] == 'city')], 'Population urban'].sum())).sum()
            )
        for country in countries_of_subregions_to_aggregate:
            data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                i, 'Native geographical resolution scale'] == 'sub-regions')], 'iF rural'] = (
                ((data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                    i, 'Native geographical resolution scale'] == 'city')], 'iF rural'] *
                  data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                      i, 'Native geographical resolution scale'] == 'city')], 'Population urban'] /
                  data.loc[[i for i in data.index if (data.loc[i, 'Country'] == country and data.loc[
                      i, 'Native geographical resolution scale'] == 'city')], 'Population urban'].sum())).sum()
            )
        for subcont in set(data.loc[:, 'Sub-Continent']):
            data.loc[[i for i in data.index if (data.loc[i, 'Sub-Continent'] == subcont and data.loc[
                i, 'Native geographical resolution scale'] == 'sub-continent')], 'iF urban'] = (
                ((data.loc[[i for i in data.index if (data.loc[i, 'Sub-Continent'] == subcont and data.loc[
                    i, 'Native geographical resolution scale'] == 'city')], 'iF urban'] *
                  data.loc[[i for i in data.index if (data.loc[i, 'Sub-Continent'] == subcont and data.loc[
                      i, 'Native geographical resolution scale'] == 'city')], 'Population urban'] /
                  data.loc[[i for i in data.index if (data.loc[i, 'Sub-Continent'] == subcont and data.loc[
                      i, 'Native geographical resolution scale'] == 'city')], 'Population urban'].sum())).sum()
            )
        for cont in set(data.loc[:, 'Continent']):
            if cont != 'Global':
                data.loc[[i for i in data.index if (data.loc[i, 'Continent'] == cont and data.loc[
                    i, 'Native geographical resolution scale'] == 'continent')], 'iF urban'] = (
                    ((data.loc[[i for i in data.index if (data.loc[i, 'Continent'] == cont and data.loc[
                        i, 'Native geographical resolution scale'] == 'city')], 'iF urban'] *
                      data.loc[[i for i in data.index if (data.loc[i, 'Continent'] == cont and data.loc[
                          i, 'Native geographical resolution scale'] == 'city')], 'Population urban'] /
                      data.loc[[i for i in data.index if (data.loc[i, 'Continent'] == cont and data.loc[
                          i, 'Native geographical resolution scale'] == 'city')], 'Population urban'].sum())).sum()
                )
            else:
                data.loc[[i for i in data.index if (data.loc[i, 'Continent'] == cont and data.loc[
                    i, 'Native geographical resolution scale'] == 'global')], 'iF urban'] = (
                    ((data.loc[[i for i in data.index if
                                data.loc[i, 'Native geographical resolution scale'] == 'city'], 'iF urban'] *
                      data.loc[[i for i in data.index if
                                data.loc[i, 'Native geographical resolution scale'] == 'city'], 'Population urban'] /
                      data.loc[[i for i in data.index if data.loc[
                          i, 'Native geographical resolution scale'] == 'city'], 'Population urban'].sum())).sum()
                )

        # remove CFs for cities because we do not provide them in the dev version
        data = data[data.loc[:, 'Native geographical resolution scale'].isin(['sub-regions', 'continent', 'global'])]

        # third we need to determine unspecified intake fractions for everyone based on population
        data.loc[:, 'iF unspecified'] = (data.loc[:, 'iF urban'] * data.loc[:, 'Population urban'] / (
                    data.loc[:, 'Population urban'] + data.loc[:, 'Population rural']) +
                                         data.loc[:, 'iF rural'] * data.loc[:, 'Population rural'] / (
                                                     data.loc[:, 'Population urban'] + data.loc[:, 'Population rural']))

        # ------------------------------------------ PM 2.5 -----------------------------------------------------
        particulate_damage = pd.Series(dtype=float)

        for region in conc.index:
            if conc.loc[region] in data.loc[:, 'Country-Region'].tolist():
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(data.loc[data.loc[:, 'Country-Region'] == conc.loc[
                                                    region], 'CF urban'].iloc[0],
                                                          index=[('Particulates, < 2.5 um, ' + region, 'high. pop.')])])
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(data.loc[data.loc[:, 'Country-Region'] == conc.loc[
                                                    region], 'CF rural'].iloc[0],
                                                          index=[('Particulates, < 2.5 um, ' + region, 'low. pop.')])])
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(data.loc[data.loc[:, 'Country-Region'] == conc.loc[
                                                    region], 'CF unspecified'].iloc[0], index=[
                                                    ('Particulates, < 2.5 um, ' + region, '(unspecified)')])])
            elif conc.loc[region] in data.loc[:, 'Sub-Continent'].tolist():
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(data.loc[data.loc[:, 'Sub-Continent'] == conc.loc[
                                                    region], 'CF urban'].iloc[0],
                                                          index=[('Particulates, < 2.5 um, ' + region, 'high. pop.')])])
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(data.loc[data.loc[:, 'Sub-Continent'] == conc.loc[
                                                    region], 'CF rural'].iloc[0],
                                                          index=[('Particulates, < 2.5 um, ' + region, 'low. pop.')])])
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(data.loc[data.loc[:, 'Sub-Continent'] == conc.loc[
                                                    region], 'CF unspecified'].iloc[0], index=[
                                                    ('Particulates, < 2.5 um, ' + region, '(unspecified)')])])
            elif conc.loc[region] in data.loc[:, 'Continent'].tolist():
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(data.loc[data.loc[:, 'Continent'] == conc.loc[
                                                    region], 'CF urban'].iloc[0],
                                                          index=[('Particulates, < 2.5 um, ' + region, 'high. pop.')])])
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(data.loc[data.loc[:, 'Continent'] == conc.loc[
                                                    region], 'CF rural'].iloc[0],
                                                          index=[('Particulates, < 2.5 um, ' + region, 'low. pop.')])])
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(data.loc[data.loc[:, 'Continent'] == conc.loc[
                                                    region], 'CF unspecified'].iloc[0], index=[
                                                    ('Particulates, < 2.5 um, ' + region, '(unspecified)')])])

        # add GLO and RoW values
        particulate_damage = pd.concat([particulate_damage,
                                        pd.Series(data.loc[data.loc[:, 'Continent'] == 'Global', 'CF urban'].iloc[0],
                                                  index=[('Particulates, < 2.5 um, GLO', 'high. pop.')])])
        particulate_damage = pd.concat([particulate_damage,
                                        pd.Series(data.loc[data.loc[:, 'Continent'] == 'Global', 'CF rural'].iloc[0],
                                                  index=[('Particulates, < 2.5 um, GLO', 'low. pop.')])])
        particulate_damage = pd.concat([particulate_damage,
                                        pd.Series(
                                            data.loc[data.loc[:, 'Continent'] == 'Global', 'CF unspecified'].iloc[0],
                                            index=[('Particulates, < 2.5 um, GLO', '(unspecified)')])])
        particulate_damage = pd.concat([particulate_damage,
                                        pd.Series(data.loc[data.loc[:, 'Continent'] == 'Global', 'CF urban'].iloc[0],
                                                  index=[('Particulates, < 2.5 um, RoW', 'high. pop.')])])
        particulate_damage = pd.concat([particulate_damage,
                                        pd.Series(data.loc[data.loc[:, 'Continent'] == 'Global', 'CF rural'].iloc[0],
                                                  index=[('Particulates, < 2.5 um, RoW', 'low. pop.')])])
        particulate_damage = pd.concat([particulate_damage,
                                        pd.Series(
                                            data.loc[data.loc[:, 'Continent'] == 'Global', 'CF unspecified'].iloc[0],
                                            index=[('Particulates, < 2.5 um, RoW', '(unspecified)')])])

        particulate_damage.index = pd.MultiIndex.from_tuples(particulate_damage.index)
        particulate_damage = particulate_damage.reset_index().rename(
            columns={'level_0': 'Elem flow name', 'level_1': 'Sub-compartment', 0: 'CF value'})
        particulate_damage.loc[:, 'Impact category'] = 'Particulate matter formation'
        particulate_damage.loc[:, 'CF unit'] = 'DALY'
        particulate_damage.loc[:, 'Elem flow unit'] = 'kg'
        particulate_damage.loc[:, 'CAS number'] = None
        particulate_damage.loc[:, 'Compartment'] = 'Air'
        particulate_damage.loc[:, 'MP or Damage'] = 'Damage'

        # ------------------------------------------ PM 10 -----------------------------------------------------
        # source = SI - TableS2 of https://doi.org/10.1021/es103563z and also in Figure 4 of
        # https://www.frontiersin.org/articles/10.3389/fenvs.2021.692440/full
        pm2_5_to_10 = 0.6
        pm10_particulate_damage = particulate_damage.copy()
        pm10_particulate_damage.loc[:, 'CF value'] *= pm2_5_to_10
        pm10_particulate_damage.loc[:, 'Elem flow name'] = [i.replace('2.5 um', '10 um') for i in
                                                            pm10_particulate_damage.loc[:, 'Elem flow name']]

        # ------------------------------------------ Secondary PM -----------------------------------------------------
        secondary_pm_if = pd.read_sql(sql='SELECT * FROM [SI - ParticulateMatter - secondary PM intake fractions]',
                                      con=self.conn).set_index("precursor")

        secondary_pm_particulate_damage = pd.DataFrame()

        so2 = particulate_damage.copy()
        so2.loc[:, 'Elem flow name'] = [i.replace('Particulates, < 2.5 um', 'Sulfur dioxide') for i in
                                        so2.loc[:, 'Elem flow name']]

        for flow in so2.index:
            region = so2.loc[flow, 'Elem flow name'].split('Sulfur dioxide, ')[1]
            if region in conc.keys():
                if conc[region] in list(data.Country):
                    if so2.loc[flow, 'Sub-compartment'] == 'high. pop.':
                        so2.loc[flow, 'CF value'] = (
                                    so2.loc[flow, 'CF value'] / data.loc[data.Country == conc[region], 'iF urban'].iloc[
                                0] * secondary_pm_if.loc['SO2', 'urban'])
                    elif so2.loc[flow, 'Sub-compartment'] == 'low. pop.':
                        so2.loc[flow, 'CF value'] = (
                                    so2.loc[flow, 'CF value'] / data.loc[data.Country == conc[region], 'iF rural'].iloc[
                                0] * secondary_pm_if.loc['SO2', 'rural'])
                    elif so2.loc[flow, 'Sub-compartment'] == '(unspecified)':
                        so2.loc[flow, 'CF value'] = (so2.loc[flow, 'CF value'] /
                                                     data.loc[data.Country == conc[region], 'iF unspecified'].iloc[0] *
                                                     secondary_pm_if.loc['SO2', 'unspecified'])
                elif conc[region] in list(data.loc[:, 'Continent']):
                    if so2.loc[flow, 'Sub-compartment'] == 'high. pop.':
                        so2.loc[flow, 'CF value'] = (so2.loc[flow, 'CF value'] /
                                                     data.loc[data.Continent == conc[region], 'iF urban'].iloc[0] *
                                                     secondary_pm_if.loc['SO2', 'urban'])
                    elif so2.loc[flow, 'Sub-compartment'] == 'low. pop.':
                        so2.loc[flow, 'CF value'] = (so2.loc[flow, 'CF value'] /
                                                     data.loc[data.Continent == conc[region], 'iF rural'].iloc[0] *
                                                     secondary_pm_if.loc['SO2', 'rural'])
                    elif so2.loc[flow, 'Sub-compartment'] == '(unspecified)':
                        so2.loc[flow, 'CF value'] = (so2.loc[flow, 'CF value'] /
                                                     data.loc[data.Continent == conc[region], 'iF unspecified'].iloc[
                                                         0] * secondary_pm_if.loc['SO2', 'unspecified'])
            else:
                # for RoW and GLO
                if so2.loc[flow, 'Sub-compartment'] == 'high. pop.':
                    so2.loc[flow, 'CF value'] = (
                                so2.loc[flow, 'CF value'] / data.loc[data.Continent == 'Global', 'iF urban'].iloc[0] *
                                secondary_pm_if.loc['SO2', 'urban'])
                elif so2.loc[flow, 'Sub-compartment'] == 'low. pop.':
                    so2.loc[flow, 'CF value'] = (
                                so2.loc[flow, 'CF value'] / data.loc[data.Continent == 'Global', 'iF rural'].iloc[0] *
                                secondary_pm_if.loc['SO2', 'rural'])
                elif so2.loc[flow, 'Sub-compartment'] == '(unspecified)':
                    so2.loc[flow, 'CF value'] = (
                                so2.loc[flow, 'CF value'] / data.loc[data.Continent == 'Global', 'iF unspecified'].iloc[
                            0] * secondary_pm_if.loc['SO2', 'unspecified'])

        nh3 = particulate_damage.copy()
        nh3.loc[:, 'Elem flow name'] = [i.replace('Particulates, < 2.5 um', 'Ammonia') for i in
                                        nh3.loc[:, 'Elem flow name']]

        for flow in nh3.index:
            region = nh3.loc[flow, 'Elem flow name'].split('Ammonia, ')[1]
            if region in conc.keys():
                if conc[region] in list(data.Country):
                    if nh3.loc[flow, 'Sub-compartment'] == 'high. pop.':
                        nh3.loc[flow, 'CF value'] = (
                                    nh3.loc[flow, 'CF value'] / data.loc[data.Country == conc[region], 'iF urban'].iloc[
                                0] * secondary_pm_if.loc['NH3', 'urban'])
                    elif nh3.loc[flow, 'Sub-compartment'] == 'low. pop.':
                        nh3.loc[flow, 'CF value'] = (
                                    nh3.loc[flow, 'CF value'] / data.loc[data.Country == conc[region], 'iF rural'].iloc[
                                0] * secondary_pm_if.loc['NH3', 'rural'])
                    elif nh3.loc[flow, 'Sub-compartment'] == '(unspecified)':
                        nh3.loc[flow, 'CF value'] = (nh3.loc[flow, 'CF value'] /
                                                     data.loc[data.Country == conc[region], 'iF unspecified'].iloc[0] *
                                                     secondary_pm_if.loc['NH3', 'unspecified'])
                elif conc[region] in list(data.loc[:, 'Continent']):
                    if nh3.loc[flow, 'Sub-compartment'] == 'high. pop.':
                        nh3.loc[flow, 'CF value'] = (nh3.loc[flow, 'CF value'] /
                                                     data.loc[data.Continent == conc[region], 'iF urban'].iloc[0] *
                                                     secondary_pm_if.loc['NH3', 'urban'])
                    elif nh3.loc[flow, 'Sub-compartment'] == 'low. pop.':
                        nh3.loc[flow, 'CF value'] = (nh3.loc[flow, 'CF value'] /
                                                     data.loc[data.Continent == conc[region], 'iF rural'].iloc[0] *
                                                     secondary_pm_if.loc['NH3', 'rural'])
                    elif nh3.loc[flow, 'Sub-compartment'] == '(unspecified)':
                        nh3.loc[flow, 'CF value'] = (nh3.loc[flow, 'CF value'] /
                                                     data.loc[data.Continent == conc[region], 'iF unspecified'].iloc[
                                                         0] * secondary_pm_if.loc['NH3', 'unspecified'])
            else:
                # for RoW and GLO
                if nh3.loc[flow, 'Sub-compartment'] == 'high. pop.':
                    nh3.loc[flow, 'CF value'] = (
                                nh3.loc[flow, 'CF value'] / data.loc[data.Continent == 'Global', 'iF urban'].iloc[0] *
                                secondary_pm_if.loc['NH3', 'urban'])
                elif nh3.loc[flow, 'Sub-compartment'] == 'low. pop.':
                    nh3.loc[flow, 'CF value'] = (
                                nh3.loc[flow, 'CF value'] / data.loc[data.Continent == 'Global', 'iF rural'].iloc[0] *
                                secondary_pm_if.loc['NH3', 'rural'])
                elif nh3.loc[flow, 'Sub-compartment'] == '(unspecified)':
                    nh3.loc[flow, 'CF value'] = (
                                nh3.loc[flow, 'CF value'] / data.loc[data.Continent == 'Global', 'iF unspecified'].iloc[
                            0] * secondary_pm_if.loc['NH3', 'unspecified'])

        nox = particulate_damage.copy()
        nox.loc[:, 'Elem flow name'] = [i.replace('Particulates, < 2.5 um', 'Nitrogen oxides') for i in
                                        nox.loc[:, 'Elem flow name']]

        for flow in nox.index:
            region = nox.loc[flow, 'Elem flow name'].split('Nitrogen oxides, ')[1]
            if region in conc.keys():
                if conc[region] in list(data.Country):
                    if nox.loc[flow, 'Sub-compartment'] == 'high. pop.':
                        nox.loc[flow, 'CF value'] = (
                                    nox.loc[flow, 'CF value'] / data.loc[data.Country == conc[region], 'iF urban'].iloc[
                                0] * secondary_pm_if.loc['NOx', 'urban'])
                    elif nox.loc[flow, 'Sub-compartment'] == 'low. pop.':
                        nox.loc[flow, 'CF value'] = (
                                    nox.loc[flow, 'CF value'] / data.loc[data.Country == conc[region], 'iF rural'].iloc[
                                0] * secondary_pm_if.loc['NOx', 'rural'])
                    elif nox.loc[flow, 'Sub-compartment'] == '(unspecified)':
                        nox.loc[flow, 'CF value'] = (nox.loc[flow, 'CF value'] /
                                                     data.loc[data.Country == conc[region], 'iF unspecified'].iloc[0] *
                                                     secondary_pm_if.loc['NOx', 'unspecified'])
                elif conc[region] in list(data.loc[:, 'Continent']):
                    if nox.loc[flow, 'Sub-compartment'] == 'high. pop.':
                        nox.loc[flow, 'CF value'] = (nox.loc[flow, 'CF value'] /
                                                     data.loc[data.Continent == conc[region], 'iF urban'].iloc[0] *
                                                     secondary_pm_if.loc['NOx', 'urban'])
                    elif nox.loc[flow, 'Sub-compartment'] == 'low. pop.':
                        nox.loc[flow, 'CF value'] = (nox.loc[flow, 'CF value'] /
                                                     data.loc[data.Continent == conc[region], 'iF rural'].iloc[0] *
                                                     secondary_pm_if.loc['NOx', 'rural'])
                    elif nox.loc[flow, 'Sub-compartment'] == '(unspecified)':
                        nox.loc[flow, 'CF value'] = (nox.loc[flow, 'CF value'] /
                                                     data.loc[data.Continent == conc[region], 'iF unspecified'].iloc[
                                                         0] * secondary_pm_if.loc['NOx', 'unspecified'])
            else:
                # for RoW and GLO
                if nox.loc[flow, 'Sub-compartment'] == 'high. pop.':
                    nox.loc[flow, 'CF value'] = (
                                nox.loc[flow, 'CF value'] / data.loc[data.Continent == 'Global', 'iF urban'].iloc[0] *
                                secondary_pm_if.loc['NOx', 'urban'])
                elif nox.loc[flow, 'Sub-compartment'] == 'low. pop.':
                    nox.loc[flow, 'CF value'] = (
                                nox.loc[flow, 'CF value'] / data.loc[data.Continent == 'Global', 'iF rural'].iloc[0] *
                                secondary_pm_if.loc['NOx', 'rural'])
                elif nox.loc[flow, 'Sub-compartment'] == '(unspecified)':
                    nox.loc[flow, 'CF value'] = (
                                nox.loc[flow, 'CF value'] / data.loc[data.Continent == 'Global', 'iF unspecified'].iloc[
                            0] * secondary_pm_if.loc['NOx', 'unspecified'])

        # concat everything
        particulate_damage = clean_up_dataframe(pd.concat([particulate_damage, pm10_particulate_damage, so2, nh3, nox]))

        # determine the midpoint
        reference_value = particulate_damage.loc[[i for i in particulate_damage.index if (
                    particulate_damage.loc[i, 'Elem flow name'] == 'Particulates, < 2.5 um, GLO' and
                    particulate_damage.loc[i, 'Sub-compartment'] == '(unspecified)')], 'CF value'].iloc[0]
        particulate_midpoint = particulate_damage.copy()
        particulate_midpoint.loc[:, 'CF value'] /= reference_value
        particulate_midpoint.loc[:, 'MP or Damage'] = 'Midpoint'
        particulate_midpoint.loc[:, 'CF unit'] = 'kg PM2.5 eq'

        # concat midpoint and damage together
        particulate_cfs = clean_up_dataframe(pd.concat([particulate_midpoint, particulate_damage]))
        particulate_cfs.loc[:, 'CF value'] = particulate_cfs.loc[:, 'CF value'].fillna(0)

        # re-establish the native geographical resolution scale
        continents = list({k for k, v in conc.items() if v in set(data.loc[:, 'Continent'].dropna())})
        particulate_cfs.loc[:, 'Native geographical resolution scale'] = [
            'Global' if ('RoW' in i or 'GLO' in i) else 'Country' for i in particulate_cfs.loc[:, 'Elem flow name']]
        particulate_cfs.loc[[i for i in particulate_cfs.index if particulate_cfs.loc[i, 'Elem flow name'].split(', ')[
            -1] in continents], 'Native geographical resolution scale'] = 'Continent'

        self.master_db = pd.concat([self.master_db, particulate_cfs])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_availability_hh_cfs(self):
        """
        Load CFs for water availability human health.

        Concerned impact categories:
            - Water availability, human health

        :return: update master_db
        """

        # ------------------------------ LOADING DATA -----------------------------------
        data = pd.read_sql('SELECT * FROM "CF - regionalized - WaterAvailability_HH - native"', self.conn)

        # ------------------------- CALCULATING DAMAGE CFS --------------------------------
        data.loc[:, 'CF value'] = (data.loc[:, 'FATE - Scarcity'] *
                                   data.loc[:, 'EXP (1 - Adaptation Capacity)'] *
                                   data.loc[:, 'EF country (DALY/m3 deprived)'])

        water_hh_cfs = pd.DataFrame()
        # countries
        for country in set(data.ISO2):
            df = data.loc[data.ISO2 == country].copy('deep')
            df.loc[:, 'Water_use_HWCtot'] /= df.loc[:, 'Water_use_HWCtot'].sum()
            water_hh_cfs.loc[country, 'CF value'] = (df.loc[:, 'Water_use_HWCtot'] * df.loc[:, 'CF value']).sum()
        # continents
        for cont in set(data.Continent):
            df = data.loc[data.Continent == cont].copy('deep')
            df.loc[:, 'Water_use_HWCtot'] /= df.loc[:, 'Water_use_HWCtot'].sum()
            water_hh_cfs.loc[cont, 'CF value'] = (df.loc[:, 'Water_use_HWCtot'] * df.loc[:, 'CF value']).sum()
        # global value
        df = data.copy('deep')
        df.loc[:, 'Water_use_HWCtot'] /= df.loc[:, 'Water_use_HWCtot'].sum()
        water_hh_cfs.loc['GLO', 'CF value'] = (df.loc[:, 'Water_use_HWCtot'] * df.loc[:, 'CF value']).sum()

        # ------------------------------ DATA FORMATTING -----------------------------------
        water_hh_cfs.loc[:, 'Impact category'] = 'Water availability, human health'
        water_hh_cfs.loc[:, 'CF unit'] = 'DALY'
        water_hh_cfs.loc[:, 'Compartment'] = 'Raw'
        water_hh_cfs.loc[:, 'Sub-compartment'] = '(unspecified)'
        water_hh_cfs.loc[:, 'CAS number'] = '7732-18-5'
        water_hh_cfs.loc[:, 'Elem flow unit'] = 'm3'
        water_hh_cfs.loc[:, 'MP or Damage'] = 'Damage'
        water_hh_cfs.loc[:, 'Native geographical resolution scale'] = 'Country'
        water_hh_cfs.loc[:, 'Elem flow name'] = ['Water, unspecified natural origin, ' + i for i in water_hh_cfs.index]
        water_hh_cfs = water_hh_cfs.reset_index().drop('index', axis=1)

        # create other water flows
        df_cooling = water_hh_cfs.copy('deep')
        df_cooling.loc[:, 'Elem flow name'] = ['Water, cooling, unspecified natural origin' + i.split(
            'Water, unspecified natural origin')[1] for i in df_cooling.loc[:, 'Elem flow name']]
        df_lake = water_hh_cfs.copy('deep')
        df_lake.loc[:, 'Elem flow name'] = ['Water, lake' + i.split('Water, unspecified natural origin')[1]
                                            for i in df_lake.loc[:, 'Elem flow name']]
        df_river = water_hh_cfs.copy('deep')
        df_river.loc[:, 'Elem flow name'] = ['Water, river' + i.split('Water, unspecified natural origin')[1]
                                             for i in df_river.loc[:, 'Elem flow name']]
        df_well = water_hh_cfs.copy('deep')
        df_well.loc[:, 'Elem flow name'] = ['Water, well, in ground' + i.split('Water, unspecified natural origin')[1]
                                            for i in df_well.loc[:, 'Elem flow name']]
        df_water = water_hh_cfs.copy('deep')
        df_water.loc[:, 'Elem flow name'] = ['Water' + i.split('Water, unspecified natural origin')[1]
                                             for i in df_water.loc[:, 'Elem flow name']]
        df_water.loc[:, 'CF value'] *= -1
        df_water.loc[:, 'Compartment'] = 'Water'

        water_hh_cfs = pd.concat([water_hh_cfs, df_cooling, df_lake, df_river, df_well, df_water])
        water_hh_cfs = clean_up_dataframe(water_hh_cfs)

        # change the resolution to continents and global
        water_hh_cfs.loc[[i for i in water_hh_cfs.index if ('RER' in water_hh_cfs.loc[i, 'Elem flow name'] or
                                                            'RAS' in water_hh_cfs.loc[i, 'Elem flow name'] or
                                                            'RAF' in water_hh_cfs.loc[i, 'Elem flow name'] or
                                                            'RLA' in water_hh_cfs.loc[i, 'Elem flow name'] or
                                                            'OCE' in water_hh_cfs.loc[i, 'Elem flow name'] or
                                                            'RNA' in water_hh_cfs.loc[i, 'Elem flow name'])],
                         'Native geographical resolution scale'] = 'Continent'
        water_hh_cfs.loc[[i for i in water_hh_cfs.index if 'GLO' in water_hh_cfs.loc[i, 'Elem flow name']],
                         'Native geographical resolution scale'] = 'Global'

        # concat with master_db
        self.master_db = pd.concat([self.master_db, water_hh_cfs])
        self.master_db = clean_up_dataframe(self.master_db)

        # add the RoW geography based on the Global value
        df = self.master_db.loc[
            [i for i in self.master_db.index if (self.master_db.loc[i, 'Impact category'] ==
                                                 'Water availability, human health' and
                                                 'GLO' in self.master_db.loc[i, 'Elem flow name'])]]
        df['Elem flow name'] = [i.split(', GLO')[0] + ', RoW' for i in df['Elem flow name']]
        df['Native geographical resolution scale'] = 'Other region'
        self.master_db = pd.concat([self.master_db, df])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_availability_eq_cfs(self):
        """
        Load CFs for water availability freshwater ecosystem.

        Concerned impact categories:
            - Water availability, freshwater ecosystem

        :return: update master_db
        """

        # ------------------------------ LOADING DATA -----------------------------------
        original_cfs = pd.read_sql(sql='SELECT * FROM [CF - regionalized - WaterAvailability_FW - native]',
                                   con=self.conn)

        geos = pd.read_sql(sql='SELECT * FROM [SI - mapping with regions of ecoinvent]', con=self.conn).loc[
               :, 'Ecoinvent_short_name'].tolist()
        geos += ['GLO', 'RoW']

        water_flows = ['Water, lake', 'Water, river', 'Water, unspecified natural origin',
                       'Water, well, in ground', 'Water, cooling, unspecified natural origin']

        # ------------------------------ DATA FORMATTING -----------------------------------
        cfs = pd.DataFrame(original_cfs.loc[:, 'CF value'].median(),
                           index=pd.MultiIndex.from_product([water_flows, geos]),
                           columns=['CF value']).reset_index()
        cfs.columns = ['Elem flow', 'Region', 'CF value']

        cfs.loc[:, 'Impact category'] = 'Water availability, freshwater ecosystem'
        cfs.loc[:, 'Native geographical resolution scale'] = 'Country'
        cfs.loc[:, 'MP or Damage'] = 'Damage'
        cfs.loc[:, 'CAS number'] = '7732-18-5'
        cfs.loc[:, 'CF unit'] = 'PDF.m2.yr'
        cfs.loc[:, 'Elem flow unit'] = 'm3'
        cfs.loc[:, 'Compartment'] = 'Raw'
        cfs.loc[:, 'Sub-compartment'] = '(unspecified)'
        # pandas interprets "NA" (for Namibia) in sqlite as NaN
        cfs.loc[:, 'Region'] = cfs.loc[:, 'Region'].fillna('NA')

        # change resolution scale for GLO and RoW
        cfs.loc[cfs.loc[:, 'Region'].isin(['GLO']), 'Native geographical resolution scale'] = 'Global'
        cfs.loc[cfs.loc[:, 'Region'].isin(['RoW']), 'Native geographical resolution scale'] = 'Other region'

        water_comp = cfs.loc[cfs.loc[:, 'Elem flow'] == 'Water, lake'].copy('deep')
        water_comp.loc[:, 'Elem flow'] = 'Water'
        water_comp.loc[:, 'Compartment'] = 'Water'
        water_comp.loc[:, 'CF value'] *= -1

        cfs = clean_up_dataframe(pd.concat([cfs, water_comp]))
        cfs.loc[:, 'Elem flow name'] = [', '.join(i) for i in list(zip(cfs.loc[:, 'Elem flow'], cfs.loc[:, 'Region']))]
        cfs = cfs.drop(['Elem flow', 'Region'], axis=1)
        cfs = clean_up_dataframe(cfs)

        # concat with master_db
        self.master_db = pd.concat([self.master_db, cfs])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_scarcity_cfs(self):
        """
        Load CFs for water scarcity

        Concerned impact categories:
            - Water scarcity

        :return: updated master_db
        """

        data = pd.read_sql('SELECT * FROM "CF - regionalized - WaterScarcity - aggregated"', self.conn)
        conc = pd.read_sql('SELECT * FROM "SI - mapping with regions of ecoinvent"', self.conn).set_index(
            'AWARE').Ecoinvent_short_name

        # create the regionalized names (e.g., Water, AF)
        water_data = pd.DataFrame()

        for i in data.index:
            if data.ecoinvent_shortname[i] in conc.index:
                if type(conc.loc[data.ecoinvent_shortname[i]]) == str:
                    if data.loc[i, 'Water type'] == 'unspecified':
                        water_data.loc[i, 'Elem flow name'] = 'Water, ' + conc.loc[data.ecoinvent_shortname[i]]
                    elif data.loc[i, 'Water type'] == 'agri':
                        water_data.loc[i, 'Elem flow name'] = 'Water, agri, ' + conc.loc[data.ecoinvent_shortname[i]]
                    elif data.loc[i, 'Water type'] == 'non-agri':
                        water_data.loc[i, 'Elem flow name'] = 'Water, non-agri, ' + conc.loc[
                            data.ecoinvent_shortname[i]]
                    water_data.loc[i, 'CF value'] = data.loc[i, 'Annual']
                else:
                    for j in range(0, len(conc.loc[data.ecoinvent_shortname[i]])):
                        if data.loc[i, 'Water type'] == 'unspecified':
                            water_data.loc[str(i + j), 'Elem flow name'] = 'Water, ' + conc.loc[data.ecoinvent_shortname[i]].iloc[j]
                        elif data.loc[i, 'Water type'] == 'agri':
                            water_data.loc[str(i + j), 'Elem flow name'] = 'Water, agri, ' + conc.loc[data.ecoinvent_shortname[i]].iloc[j]
                        elif data.loc[i, 'Water type'] == 'non-agri':
                            water_data.loc[str(i + j), 'Elem flow name'] = 'Water, non-agri, ' + conc.loc[data.ecoinvent_shortname[i]].iloc[j]
                        water_data.loc[str(i + j), 'CF value'] = data.loc[i, 'Annual']
            else:
                if data.loc[i, 'Water type'] == 'unspecified':
                    water_data.loc[i, 'Elem flow name'] = 'Water, ' + data.ecoinvent_shortname[i]
                elif data.loc[i, 'Water type'] == 'agri':
                    water_data.loc[i, 'Elem flow name'] = 'Water, agri, ' + data.ecoinvent_shortname[i]
                elif data.loc[i, 'Water type'] == 'non-agri':
                    water_data.loc[i, 'Elem flow name'] = 'Water, non-agri, ' + data.ecoinvent_shortname[i]
                water_data.loc[i, 'CF value'] = data.loc[i, 'Annual']

        # formatting the data to IW+ format
        water_data.loc[:, 'Impact category'] = 'Water scarcity'
        water_data.loc[:, 'CF unit'] = 'm3 world-eq'
        water_data.loc[:, 'Compartment'] = 'Raw'
        water_data.loc[:, 'Sub-compartment'] = '(unspecified)'
        water_data.loc[:, 'CAS number'] = '7732-18-5'
        water_data.loc[:, 'Elem flow unit'] = 'm3'
        water_data.loc[:, 'MP or Damage'] = 'Midpoint'
        water_data.loc[:, 'Native geographical resolution scale'] = 'Country'

        water_data.loc[:, 'CF value'] = water_data.loc[:, 'CF value'].astype(float)

        # create the negative flows for the Water compartment
        water_extraction_data = water_data.copy('deep')
        water_extraction_data.loc[:, 'Compartment'] = 'Water'
        water_extraction_data.loc[:, 'CF value'] *= -1

        all_data = pd.concat([water_data, water_extraction_data])
        all_data = clean_up_dataframe(all_data)

        all_data.loc[[i for i in all_data.index if ('RER' in all_data.loc[i, 'Elem flow name'] or
                                                    'RAS' in all_data.loc[i, 'Elem flow name'] or
                                                    'RAF' in all_data.loc[i, 'Elem flow name'] or
                                                    'RLA' in all_data.loc[i, 'Elem flow name'] or
                                                    'OCE' in all_data.loc[i, 'Elem flow name'] or
                                                    'RNA' in all_data.loc[i, 'Elem flow name'])],
                         'Native geographical resolution scale'] = 'Continent'
        all_data.loc[[i for i in all_data.index if 'GLO' in all_data.loc[i, 'Elem flow name']],
                         'Native geographical resolution scale'] = 'Global'
        all_data.loc[[i for i in all_data.index if 'RoW' in all_data.loc[i, 'Elem flow name']],
                         'Native geographical resolution scale'] = 'Other region'

        self.master_db = pd.concat([self.master_db, all_data])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_availability_terr_cfs(self):
        """
        Load CFs for water availability terrestrial ecosystem.

        Concerned impact categories:
            - Water availability, terrestrial ecosystem

        :return: update master_db
        """

        original_cfs = pd.read_sql('SELECT * FROM "CF - not regionalized - WaterAvailabilityTerrestrial"',
                                   self.conn)

        geos = pd.read_sql(sql='SELECT * FROM [SI - mapping with regions of ecoinvent]', con=self.conn).loc[
               :,'Ecoinvent_short_name'].tolist()
        geos += ['GLO', 'RoW']

        list_flows = pd.MultiIndex.from_product([original_cfs.loc[:, 'Elem flow name'], geos])
        cfs = pd.DataFrame([i[0] + ', ' + i[1] for i in list_flows], columns=['Elem flow name'])

        cfs.loc[:, 'Impact category'] = 'Water availability, terrestrial ecosystem'
        cfs.loc[:, 'Native geographical resolution scale'] = 'Country'
        cfs.loc[:, 'MP or Damage'] = 'Damage'
        cfs.loc[:, 'CAS number'] = '7732-18-5'
        cfs.loc[:, 'CF unit'] = 'PDF.m2.yr'
        cfs.loc[:, 'Elem flow unit'] = 'm3'
        cfs.loc[:, 'Sub-compartment'] = 'in water'
        cfs.loc[:, 'Compartment'] = 'Raw'
        cfs.loc[:, 'CF value'] = 0

        cfs.loc[[i for i in cfs.index if 'GLO' in cfs.loc[i, 'Elem flow name']],
                'Native geographical resolution scale'] = 'Global'
        cfs.loc[[i for i in cfs.index if 'RoW' in cfs.loc[i, 'Elem flow name']],
                'Native geographical resolution scale'] = 'Other region'

        cfs.loc[[i for i in cfs.index if 'shallow' in cfs.loc[i, 'Elem flow name']], 'CF value'] = (
            original_cfs.loc[[i for i in original_cfs.index if
                              'shallow' in original_cfs.loc[i, 'Elem flow name']], 'CF value'].iloc[0]
        )

        cfs.loc[[i for i in cfs.index if 'shallow' not in cfs.loc[i, 'Elem flow name']], 'CF value'] = (
            original_cfs.loc[[i for i in original_cfs.index if
                              'shallow' not in original_cfs.loc[i, 'Elem flow name']], 'CF value'].iloc[0]
        )

        self.master_db = pd.concat([self.master_db, cfs])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_thermally_polluted_water_cfs(self):
        """
        Load CFs for thermally polluted water.

        Concerned impact categories:
            - Thermally polluted water

        :return: update master_db
        """

        db = pd.read_sql('SELECT * FROM "CF - not regionalized - ThermallyPollutedWater"', self.conn)
        geos = pd.read_sql(sql='SELECT * FROM [SI - mapping with regions of ecoinvent]', con=self.conn).loc[
               :, 'Ecoinvent_short_name'].tolist()
        geos += ['GLO', 'RoW']

        df = db.loc[0].copy('deep')
        for i, geo in enumerate(geos):
            db.loc[i + 1] = df
            db.loc[i + 1, 'Elem flow name'] = db.loc[0, 'Elem flow name'] + ', ' + geo
            if geo == 'GLO':
                db.loc[i + 1, 'Native geographical resolution scale'] = 'Global'
            elif geo == 'RoW':
                db.loc[i + 1, 'Native geographical resolution scale'] = 'Other region'
            else:
                db.loc[i + 1, 'Native geographical resolution scale'] = 'Country'
        db = db.drop(0)

        self.master_db = pd.concat([self.master_db, db])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_plastic_cfs(self):
        """
        Load CFs for plastics impact.

        Concerned impact categories:
            - Plastics physical effect on biota

        :return: update master_db
        """

        original_cfs = pd.read_sql('SELECT * from "CF - not regionalized - PhysicalImpactonBiota"',
                                   con=self.conn)

        original_cfs.drop(['Geometric st.dev.', 'Lower limit 95% CI (calculated)',
                           'Upper limit 95% CI (calculated)'], axis=1, inplace=True)
        original_cfs.loc[:, 'Impact category'] = 'Plastics physical effects on biota'
        original_cfs.loc[:, 'Native geographical resolution scale'] = 'Global'
        original_cfs.loc[:, 'MP or Damage'] = ['Midpoint' if original_cfs.loc[i, 'CF unit'] == 'CTUe' else 'Damage' for
                                               i in original_cfs.index]
        name_changes = {'EPS': 'Expanded Polystyrene (EPS)',
                        'HDPE': 'High density Polyethylene (HDPE)',
                        'LDPE': 'Low Density Polyethylene (LDPE)',
                        'PA (Nylon)': 'Polyamide/Nylon (PA)',
                        'PET': 'Polyethylene Terephtalate (PET)',
                        'PHA': 'Polyhydroxyalkanoates (PHA)',
                        'PLA': 'Polyactic acid (PLA)',
                        'PP': 'Polypropylene (PP)',
                        'PS': 'Polystyrene (PS)',
                        'PVC': 'Polyvinyl Chloride (PVC)',
                        'TRWP': 'Tyre and Road Wear Particles (TRWP)',
                        'PU/Spandex': 'Polyurethane (PU)',
                        'Acrylic': 'Acrylic polymers'}
        CAS = {'EPS': '9003-53-6',
               'HDPE': '9002-88-4',
               'LDPE': '9002-88-4',
               'PA (Nylon)': '25038-54-4',
               'PET': '25038-59-9',
               'PHA': '117068-64-1',
               'PLA': '26100-51-6',
               'PP': '9003-07-0',
               'PS': '9003-53-6',
               'PVC': '9002-86-2',
               'TRWP': None,
               'PU/Spandex': '9009-54-5',
               'Acrylic': '9065-11-6'}
        original_cfs.loc[:, 'CAS number'] = [CAS[i] for i in original_cfs.loc[:, 'Polymer type']]
        original_cfs.loc[:, 'Polymer type'] = [name_changes[i] for i in original_cfs.loc[:, 'Polymer type']]
        original_cfs = original_cfs.rename(columns={'Recommended CF (geometric mean)': 'CF value'})
        original_cfs.loc[:, 'Elem flow name'] = [
            original_cfs.loc[i, 'Polymer type'] + ' - ' + original_cfs.loc[i, 'Shape'] + ' - ' +
            str(original_cfs.loc[i, 'Size']) + 'um' for i in original_cfs.index]
        original_cfs.drop(['Polymer type', 'Size', 'Shape'], axis=1, inplace=True)

        self.master_db = pd.concat([self.master_db, original_cfs])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_fisheries_cfs(self):
        """
        Load CFs for fisheries impact.

        Concerned impact categories:
            - Fisheries impact

        :return: update master_db
        """

        cfs = pd.read_sql("SELECT * FROM [CF - regionalized - Fisheries]", self.conn)

        # only keep data that are classifies as the most robust, i.e., class I
        cfs = cfs.loc[cfs.loc[:, 'class'] == 'I']
        # change FAO_zone from numbers to strings
        cfs.loc[:, 'FAO_num'] = [str(i) for i in cfs.loc[:, 'FAO_num']]
        # convert the species/yr cf to PDF.m2.yr
        cfs.loc[:, 'CF (PDF.m2.yr)'] = (
                    cfs.loc[:, 'CF (species/yr)'] * cfs.loc[:, 'Area (m2)'] / cfs.loc[:, 'Species_num (nb sp.)'])
        # original cfs per tonnes of fish. We Want per kg
        cfs.loc[:, 'CF (PDF.m2.yr)'] /= 1000

        # aggregate impacts of fish per type of fish, i.e., pelagic or demersal
        cf_regions = pd.DataFrame(
            index=pd.MultiIndex.from_product([list(set(cfs.loc[:, 'FAO_num'])), ['Demersal', 'Pelagic']]),
            columns=['CF (PDF.m2.yr)'])
        for FAO_zone in set(cfs.loc[:, 'FAO_num']):
            biomass_demersal_in_zone = 0
            biomass_pelagic_in_zone = 0
            try:
                cf_regions.loc[[(FAO_zone, 'Demersal')], 'CF (PDF.m2.yr)'] = gmean(
                    cfs.loc[cfs.loc[:, 'FAO_num'] == FAO_zone].loc[cfs.loc[:, 'Type'] == 'Demersal', 'CF (PDF.m2.yr)'],
                    weights=cfs.loc[cfs.loc[:, 'FAO_num'] == FAO_zone].loc[
                        cfs.loc[:, 'Type'] == 'Demersal', 'B (tonnes)'])
                biomass_demersal_in_zone = cfs.loc[cfs.loc[:, 'FAO_num'] == FAO_zone].loc[
                    cfs.loc[:, 'Type'] == 'Demersal', 'B (tonnes)'].sum()
            except ZeroDivisionError:
                pass
            try:
                cf_regions.loc[[(FAO_zone, 'Pelagic')], 'CF (PDF.m2.yr)'] = gmean(
                    cfs.loc[cfs.loc[:, 'FAO_num'] == FAO_zone].loc[cfs.loc[:, 'Type'] == 'Pelagic', 'CF (PDF.m2.yr)'],
                    weights=cfs.loc[cfs.loc[:, 'FAO_num'] == FAO_zone].loc[
                        cfs.loc[:, 'Type'] == 'Pelagic', 'B (tonnes)'])
                biomass_pelagic_in_zone = cfs.loc[cfs.loc[:, 'FAO_num'] == FAO_zone].loc[
                    cfs.loc[:, 'Type'] == 'Pelagic', 'B (tonnes)'].sum()
            except ZeroDivisionError:
                pass

            total_biomass_in_zone = biomass_demersal_in_zone + biomass_pelagic_in_zone
            biomass_demersal_in_zone /= total_biomass_in_zone
            biomass_pelagic_in_zone /= total_biomass_in_zone

            # determine the Discard CF per region per type of fish
            cf_regions.loc[FAO_zone, 'Discard CF (PDF.m2.yr/t discarded fish)'] = math.exp(
                math.log(cf_regions.loc[[(FAO_zone, 'Pelagic')], 'CF (PDF.m2.yr)'].iloc[0]) * biomass_pelagic_in_zone +
                math.log(cf_regions.loc[[(FAO_zone, 'Demersal')], 'CF (PDF.m2.yr)'].iloc[0]) * biomass_demersal_in_zone)

        # determine global values for CF and discard CF
        glo = pd.DataFrame([gmean(cfs.loc[cfs.loc[:, 'Type'] == 'Demersal', 'CF (PDF.m2.yr)'],
                                  weights=cfs.loc[cfs.loc[:, 'Type'] == 'Demersal', 'B (tonnes)']),
                            gmean(cfs.loc[cfs.loc[:, 'Type'] == 'Pelagic', 'CF (PDF.m2.yr)'],
                                  weights=cfs.loc[cfs.loc[:, 'Type'] == 'Pelagic', 'B (tonnes)'])],
                           index=pd.MultiIndex.from_product([['GLO'], ['Demersal', 'Pelagic']]),
                           columns=['CF (PDF.m2.yr)'])

        glo.loc[:, 'Discard CF (PDF.m2.yr/t discarded fish)'] = (
            math.exp(math.log(glo.loc[[('GLO', 'Pelagic')], 'CF (PDF.m2.yr)'].iloc[0]) *
                     cfs.loc[cfs.loc[:, 'Type'] == 'Pelagic', 'B (tonnes)'].sum() /
                     (cfs.loc[cfs.loc[:, 'Type'] == 'Demersal', 'B (tonnes)'].sum() +
                      cfs.loc[cfs.loc[:, 'Type'] == 'Pelagic', 'B (tonnes)'].sum()) +
                     math.log(glo.loc[[('GLO', 'Demersal')], 'CF (PDF.m2.yr)'].iloc[0]) *
                     cfs.loc[cfs.loc[:, 'Type'] == 'Demersal', 'B (tonnes)'].sum() /
                     (cfs.loc[cfs.loc[:, 'Type'] == 'Demersal', 'B (tonnes)'].sum() +
                      cfs.loc[cfs.loc[:, 'Type'] == 'Pelagic', 'B (tonnes)'].sum())))

        # for FAO zones without data, use the global average
        cf_regions.loc[:, 'Discard CF (PDF.m2.yr/t discarded fish)'] = cf_regions.loc[:,
                                                                       'Discard CF (PDF.m2.yr/t discarded fish)'].fillna(
            glo.loc[:, 'Discard CF (PDF.m2.yr/t discarded fish)'].iloc[0])
        cf_regions.loc[[i for i in cf_regions.index if i[1] == 'Pelagic'], 'CF (PDF.m2.yr)'] = cf_regions.loc[
            [i for i in cf_regions.index if i[1] == 'Pelagic'], 'CF (PDF.m2.yr)'].fillna(
            glo.loc[('GLO', 'Pelagic'), 'CF (PDF.m2.yr)'])
        cf_regions.loc[[i for i in cf_regions.index if i[1] == 'Demersal'], 'CF (PDF.m2.yr)'] = cf_regions.loc[
            [i for i in cf_regions.index if i[1] == 'Demersal'], 'CF (PDF.m2.yr)'].fillna(
            glo.loc[('GLO', 'Demersal'), 'CF (PDF.m2.yr)'])

        data = pd.DataFrame()

        data.loc[:, 'Elem flow name'] = cfs.loc[:, 'ASFIS_spp_common'] + ', FAO zone ' + cfs.loc[:, 'FAO_num']
        data.loc[:, 'CF value'] = cfs.loc[:, 'CF (PDF.m2.yr)']

        data = pd.concat([data, pd.concat(
            [pd.DataFrame([i[1] + ' fish , FAO zone ' + i[0] for i in cf_regions.index], columns=['Elem flow name']),
             pd.DataFrame(cf_regions.loc[:, 'CF (PDF.m2.yr)'].values, columns=['CF value'])], axis=1)])

        data = pd.concat(
            [data, pd.concat([pd.DataFrame([i[1] + ' fish, GLO' for i in glo.index], columns=['Elem flow name']),
                              pd.DataFrame(glo.loc[:, 'CF (PDF.m2.yr)'].values, columns=['CF value'])], axis=1)])

        data = pd.concat([data, pd.concat([pd.DataFrame(
            [i[1] + ' fish, discarded, FAO zone ' + i[0] for i in cf_regions.index], columns=['Elem flow name']),
                                           pd.DataFrame(
                                               cf_regions.loc[:, 'Discard CF (PDF.m2.yr/t discarded fish)'].values,
                                               columns=['CF value'])], axis=1)])

        data = pd.concat(
            [data, pd.concat([pd.DataFrame([i[1] + ' fish, discarded, GLO' for i in glo.index], columns=['Elem flow name']),
                              pd.DataFrame(glo.loc[:, 'Discard CF (PDF.m2.yr/t discarded fish)'].values,
                                           columns=['CF value'])], axis=1)])
        data.loc[:, 'Impact category'] = 'Fisheries impact'
        data.loc[:, 'CF unit'] = 'PDF.m2.yr'
        data.loc[[i for i in data.index if 'discarded' not in data.loc[i, 'Elem flow name']], 'Compartment'] = 'Raw'
        data.loc[[i for i in data.index if 'discarded' in data.loc[i, 'Elem flow name']], 'Compartment'] = 'Water'
        data.loc[
            [i for i in data.index if 'discarded' not in data.loc[i, 'Elem flow name']], 'Sub-compartment'] = 'biotic'
        data.loc[[i for i in data.index if 'discarded' in data.loc[i, 'Elem flow name']], 'Sub-compartment'] = 'ocean'
        data.loc[:, 'Elem flow unit'] = 'kg'
        data.loc[:, 'MP or Damage'] = 'Damage'
        data.loc[:, 'Native geographical resolution scale'] = 'Country'
        data = clean_up_dataframe(data)
        data.loc[[i for i in data.index if ', GLO' in data.loc[i, 'Elem flow name']],
                 'Native geographical resolution scale'] = 'Global'

        self.master_db = pd.concat([self.master_db, data])
        self.master_db = clean_up_dataframe(self.master_db)

    def apply_rules(self):
        """
        Applying rules creating values for new sub compartments, for each substance.  The value is either equal to the
        unspecified subcomp value or is fixed to zero. If a value already exists for a sub-compartment to be created,
        this value is kept. Most of these subcomps are created to match with the subcomps used by the ecoinvent database.

        Created comp/subcomps:
            - Air/high. pop.
            - Air/low. pop.
            - Air/low. pop., long term
            - Air/stratosphere + troposphere
            - Air/indoor
            - Water/lake
            - Water/river
            - Water/ocean
            - Water/groundwater
            - Water/groundwater, long term
            - Soil/industrial
            - Raw/in ground
            - Raw/in water
            - Raw/biotic

        :return: updated master_db
        """

        # ---------------- Equal to unspecified -------------

        subcomps = {'Air': ['high. pop.', 'low. pop.', 'stratosphere + troposphere', 'indoor'],
                    'Water': ['lake', 'river'],
                    'Soil': ['industrial', 'agricultural'],
                    'Raw': ['in ground', 'in water', 'biotic']}

        # create a dataframe in which all new subcomps are created based on the unspecified value
        proxy = pd.DataFrame()
        for comp in subcomps.keys():
            if comp != 'Raw':
                unspecified = self.master_db.loc[
                    [i for i in self.master_db.index if (self.master_db.loc[i, 'Compartment'] == comp and
                                                         self.master_db.loc[i, 'Sub-compartment'] == '(unspecified)')]]
                dff = unspecified.copy()
                for subcomp in subcomps[comp]:
                    df = unspecified.copy()
                    df.loc[:, 'Sub-compartment'] = subcomp
                    dff = pd.concat([dff, df])
                proxy = pd.concat([proxy, dff])
            # Special case for Raw. The subcomp depends on the impact categories
            else:
                # For the Fossil and nuclear energy use category
                unspecified = self.master_db.loc[
                    [i for i in self.master_db.index if (self.master_db.loc[i, 'Compartment'] == 'Raw' and
                                                         self.master_db.loc[i, 'Sub-compartment'] == '(unspecified)' and
                                                         self.master_db.loc[
                                                             i, 'Impact category'] == 'Fossil and nuclear energy use')]]
                biotic = unspecified.loc[
                    [i for i in unspecified.index if ('wood' in unspecified.loc[i, 'Elem flow name'].lower() or
                                                      'peat' in unspecified.loc[
                                                          i, 'Elem flow name'].lower())]].copy()
                fossil = unspecified.loc[
                    [i for i in unspecified.index if ('wood' not in unspecified.loc[i, 'Elem flow name'].lower() and
                                                      'peat' not in unspecified.loc[
                                                          i, 'Elem flow name'].lower())]].copy()
                biotic.loc[:, 'Sub-compartment'] = 'biotic'
                fossil.loc[:, 'Sub-compartment'] = 'in ground'
                proxy = pd.concat([proxy, unspecified, biotic, fossil])

                # For the water flows
                unspecified = self.master_db.loc[
                    [i for i in self.master_db.index if (
                            self.master_db.loc[i, 'Compartment'] == 'Raw' and
                            self.master_db.loc[i, 'Sub-compartment'] == '(unspecified)' and
                            self.master_db.loc[i, 'Impact category'] in ['Water scarcity',
                                                                         'Thermally polluted water',
                                                                         'Water availability, terrestrial ecosystem',
                                                                         'Water availability, freshwater ecosystem',
                                                                         'Water availability, human health']
                    )]]
                dff = unspecified.copy()
                dff.loc[:, 'Sub-compartment'] = 'in water'
                proxy = pd.concat([proxy, dff])
                dff = unspecified.copy()
                dff.loc[:, 'Sub-compartment'] = 'in ground'
                proxy = pd.concat([proxy, dff])

            proxy = clean_up_dataframe(proxy)

        # only add new subcomps values if they do not already exist in master_db
        proxy.set_index(['Impact category', 'CF unit', 'Compartment', 'Sub-compartment', 'Elem flow name'],inplace=True)
        proxy.update(self.master_db.set_index(['Impact category', 'CF unit', 'Compartment',
                                               'Sub-compartment', 'Elem flow name']))
        proxy = proxy.reset_index()
        self.master_db = pd.concat([self.master_db, proxy]).drop_duplicates()

        self.master_db = clean_up_dataframe(self.master_db)

        # ------------ groundwater and ocean subcomps ------------

        # for the groundwater and ocean subcomps, only in some impact categories are the values equal to unspecified
        water_comp = self.master_db.loc[[i for i in self.master_db.index if self.master_db.loc[i, 'Compartment'] == 'Water']]
        to_unspecified = {'groundwater': ['Water availability, freshwater ecosystem',
                                          'Water availability, human health',
                                          'Water scarcity'],
                          'groundwater, long-term': ['Water availability, freshwater ecosystem',
                                                     'Water availability, human health',
                                                     'Water scarcity'],
                          'ocean': ['Marine eutrophication']}
        to_zero = {'groundwater': ['Freshwater ecotoxicity',
                                   'Freshwater ecotoxicity, long term',
                                   'Freshwater ecotoxicity, short term',
                                   'Freshwater eutrophication',
                                   'Human toxicity cancer',
                                   'Human toxicity cancer, long term',
                                   'Human toxicity cancer, short term',
                                   'Human toxicity non-cancer',
                                   'Human toxicity non-cancer, long term',
                                   'Human toxicity non-cancer, short term',
                                   'Ionizing radiation, ecosystem quality',
                                   'Ionizing radiation, human health',
                                   'Ionizing radiations',
                                   'Marine eutrophication'],
                   'groundwater, long-term': ['Freshwater ecotoxicity',
                                              'Freshwater ecotoxicity, long term',
                                              'Freshwater ecotoxicity, short term',
                                              'Freshwater eutrophication',
                                              'Human toxicity cancer',
                                              'Human toxicity cancer, long term',
                                              'Human toxicity cancer, short term',
                                              'Human toxicity non-cancer',
                                              'Human toxicity non-cancer, long term',
                                              'Human toxicity non-cancer, short term',
                                              'Ionizing radiation, ecosystem quality',
                                              'Ionizing radiation, human health',
                                              'Ionizing radiations',
                                              'Marine eutrophication'],
                   'ocean': ['Freshwater ecotoxicity',
                             'Freshwater ecotoxicity, long term',
                             'Freshwater ecotoxicity, short term',
                             'Freshwater eutrophication',
                             'Ionizing radiation, ecosystem quality',
                             'Human toxicity cancer',
                             'Human toxicity cancer, long term',
                             'Human toxicity cancer, short term',
                             'Human toxicity non-cancer',
                             'Human toxicity non-cancer, long term',
                             'Human toxicity non-cancer, short term',
                             'Water availability, freshwater ecosystem',
                             'Water availability, human health',
                             'Water scarcity']}
        # if an impact category already has specified groundwater flows, use that value instead of the unspecified subcomp value
        already_has_groundwater_values = set(
            water_comp[water_comp['Sub-compartment'] == 'groundwater']['Impact category'])

        for subcomp in to_unspecified:
            for cat in to_unspecified[subcomp]:
                if cat not in already_has_groundwater_values:
                    data = water_comp.loc[
                        [i for i in water_comp.index if (water_comp.loc[i, 'Impact category'] == cat and
                                                         water_comp.loc[
                                                             i, 'Sub-compartment'] == '(unspecified)')]]
                else:
                    data = water_comp.loc[
                        [i for i in water_comp.index if (water_comp.loc[i, 'Impact category'] == cat and
                                                         water_comp.loc[
                                                             i, 'Sub-compartment'] == 'groundwater')]]
                proxy = data.copy()
                proxy.loc[:, 'Sub-compartment'] = subcomp
                proxy.set_index(['Impact category', 'CF unit', 'Compartment', 'Sub-compartment', 'Elem flow name'],
                                inplace=True)
                proxy.update(self.master_db.set_index(['Impact category', 'CF unit', 'Compartment',
                                                       'Sub-compartment', 'Elem flow name']))
                proxy = proxy.reset_index()
                self.master_db = pd.concat([self.master_db, proxy]).drop_duplicates()
                self.master_db = clean_up_dataframe(self.master_db)

        for subcomp in to_zero:
            for cat in to_zero[subcomp]:
                data = water_comp.loc[[i for i in water_comp.index if (water_comp.loc[i, 'Impact category'] == cat and
                                                                       water_comp.loc[
                                                                           i, 'Sub-compartment'] == '(unspecified)')]]
                proxy = data.copy()
                proxy.loc[:, 'Sub-compartment'] = subcomp
                proxy.loc[:, 'CF value'] = 0
                proxy.set_index(['Impact category', 'CF unit', 'Compartment', 'Sub-compartment', 'Elem flow name'],
                                inplace=True)
                proxy.update(self.master_db.set_index(['Impact category', 'CF unit', 'Compartment',
                                                  'Sub-compartment', 'Elem flow name']))
                proxy = proxy.reset_index()
                self.master_db = pd.concat([self.master_db, proxy]).drop_duplicates()
                self.master_db = clean_up_dataframe(self.master_db)

        # -------------- low. pop., long-term --------------

        long_term_subcomps = {'low. pop., long-term': 'low. pop.'}

        long_term_cats = ['Climate change, ecosystem quality', 'Climate change, human health',
                          'Freshwater ecotoxicity', 'Human toxicity cancer',
                          'Human toxicity non-cancer', 'Marine acidification']

        for subcomp in long_term_subcomps.keys():
            # slice dataframe to only keep the corresponding subcomp (low. pop., long-term)
            data = self.master_db.loc[[i for i in self.master_db.index if (self.master_db.loc[i, 'Sub-compartment']
                                                                           == long_term_subcomps[subcomp])]].copy()
            # remove already existing subcomp values, these values will be redefined properly in this function
            self.master_db.drop(
                [i for i in self.master_db.index if (self.master_db.loc[i, 'Sub-compartment'] == subcomp and
                                                     ','.join(self.master_db.loc[i, 'Impact category'].split(',')[
                                                              :-1]) in long_term_cats)], inplace=True)
            for cat in long_term_cats:
                # slice dataframe to only keep corresponding impact category
                df = data.loc[[i for i in data.index if (cat in data.loc[i, 'Impact category'] and
                                                         cat != data.loc[i, 'Impact category'])]]
                # remove the "short term" and "long term" from impact category name
                df.loc[:, 'Impact category'] = [','.join(i.split(',')[:-1]) for i in df.loc[:, 'Impact category']]
                # now that they have the same category name, we can add their CF value by merging dataframes
                df = df.groupby(by=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment', 'Elem flow name',
                                    'CAS number', 'Elem flow unit', 'MP or Damage',
                                    # dropna=False because there are some NaNs with CAS numbers, which you remove some elementary flows from the dataframe
                                    'Native geographical resolution scale'], dropna=False).sum().reset_index()
                # rename the subcomp
                df.loc[:, 'Sub-compartment'] = subcomp
                # reimplement the "long term" in the imapct category name
                df.loc[:, 'Impact category'] = [i + ', long term' for i in df.loc[:, 'Impact category']]

                # concatenate new CFs to master_db
                self.master_db = pd.concat([self.master_db, df])
                # now change "long term" in impact category name for "short term"
                df.loc[:, 'Impact category'] = [i.replace('long term', 'short term') for i in
                                                df.loc[:, 'Impact category']]
                # fix the value for these CFs to zero
                df.loc[:, 'CF value'] = 0
                # and simply concatenate them to master_db as well
                self.master_db = pd.concat([self.master_db, df])
                # clean up index
                self.master_db = clean_up_dataframe(self.master_db)

        # now for midpoint categories
        for subcomp in long_term_subcomps.keys():
            data = self.master_db.loc[[i for i in self.master_db.index if (self.master_db.loc[i, 'Sub-compartment']
                                                                           == long_term_subcomps[subcomp] and
                                                                           self.master_db.loc[
                                                                               i, 'MP or Damage'] == 'Midpoint')]].copy()
            df = data.copy()
            df.loc[:, 'Sub-compartment'] = subcomp
            self.master_db = pd.concat([self.master_db, df])
            # clean up index
            self.master_db = clean_up_dataframe(self.master_db)

        # finally for the endpoint categories that have no long/short term differentiation
        ics = ['Marine eutrophication', 'Ozone layer depletion', 'Photochemical oxidant formation',
               'Terrestrial acidification', 'Particulate matter formation', 'Ionizing radiation, ecosystem quality',
               'Ionizing radiation, human health', 'Freshwater acidification']
        for subcomp in long_term_subcomps.keys():
            data = self.master_db.loc[[i for i in self.master_db.index if (self.master_db.loc[i, 'Sub-compartment']
                                                                           == long_term_subcomps[subcomp] and
                                                                           self.master_db.loc[
                                                                               i, 'Impact category'] in ics and
                                                                           self.master_db.loc[
                                                                               i, 'MP or Damage'] == 'Damage')]].copy()
            df = data.copy()
            df.loc[:, 'Sub-compartment'] = subcomp
            self.master_db = pd.concat([self.master_db, df])
            # clean up index
            self.master_db = clean_up_dataframe(self.master_db)

        # -------------- Shallow groundwater ----------------

        # water availability terrestrial ecosystem introduces a new water flow type: "water, shallow well, in ground"
        # its impact for all other water impact categories should be a copy of the typical groundwater CF




        # ----------------- Special cases --------------------

        # special case from/to soil or biomass flows should only be defined for short term damage categories so
        # remove the low. pop., long-term flow for them
        self.master_db = self.master_db.drop([i for i in self.master_db.index if (
                "soil or biomass" in self.master_db.loc[i,'Elem flow name'] and
                self.master_db.loc[i,'Impact category'] in [
                    'Climate change, ecosystem quality, long term',
                    'Climate change, human health, long term',
                    'Marine acidification, long term'])])

        # add zero flows for saline water to make it explicit for the user
        df = self.master_db[self.master_db['Elem flow name'] == 'Water, lake, GLO'].copy()
        df['Elem flow name'] = 'Water, salt, ocean'
        df['CF value'] = 0
        self.master_db = pd.concat([self.master_db, df])
        self.master_db = clean_up_dataframe(self.master_db)
        # add zero flows for saline water to make it explicit for the user
        df = self.master_db[self.master_db['Elem flow name'] == 'Water, lake, GLO'].copy()
        df['Elem flow name'] = 'Water, salt, sole'
        df['CF value'] = 0
        self.master_db = pd.concat([self.master_db, df])
        self.master_db = clean_up_dataframe(self.master_db)

    def create_not_regio_flows(self):
        """
        Method creates not regionalized flows (e.g., "Ammonia") from global values (e.g., "Ammonia, GLO"). Those flows
        are necessary for linkage to databases accepting regionalized flows, such as SimaPro.
        :return:
        """

        global_values = [i for i in self.master_db.index if ', GLO' in
                         self.master_db.loc[i, 'Elem flow name']]
        df = self.master_db.loc[global_values]
        df['Elem flow name'] = [i.split(', GLO')[0] for i in df['Elem flow name']]
        self.master_db = pd.concat([self.master_db, df])
        self.master_db = clean_up_dataframe(self.master_db)

    def create_regio_flows_for_not_regio_ic(self):
        """
        Some regionalized emissions (e.g., Ammonia) also impact non-regionalized impact categories (e.g.,
        Particulate matter). Hence, this method create non-regionalized CFs for these emissions.
        :return:
        """

        # identify regionalized flows, regions and impact categories

        regio_flows = set([', '.join(i.split(', ')[:-1]) for i in
                           self.master_db[self.master_db['Native geographical resolution scale'] == 'Country'].loc[
                           :, 'Elem flow name']])

        regio_regions = set(
            [i.split(', ')[-1] for i in self.master_db[self.master_db['Native geographical resolution scale'] ==
                                                       'Country'].loc[:, 'Elem flow name']])

        regio_ic = set(self.master_db[self.master_db['Native geographical resolution scale'] ==
                                      'Country'].loc[:, 'Impact category'])

        # identify which regionalized flows need to be characterized for non-regionalized impact categories

        flows_to_create = {}

        for ic in set(self.master_db.loc[:, 'Impact category']):
            if ic not in regio_ic:
                df = self.master_db[self.master_db['Impact category'] == ic].copy()
                for flow in regio_flows:
                    if flow in df['Elem flow name'].tolist():
                        # if the CF is equal to zero we don't care
                        if df.loc[[i for i in df.index if (df.loc[i, 'Elem flow name'] == flow and
                                                           df.loc[i, 'Impact category'] == ic)], 'CF value'].sum() != 0:
                            if ic in flows_to_create.keys():
                                flows_to_create[ic].append(flow)
                            else:
                                flows_to_create[ic] = [flow]

        # create CF for the previously identified flows

        for ic in flows_to_create.keys():
            df = self.master_db[self.master_db['Impact category'] == ic].copy()
            for flow in flows_to_create[ic]:
                dff = df[df['Elem flow name'] == flow]
                list_flow_added = ([flow + ', ' + i for i in regio_regions] * len(dff))
                list_flow_added.sort()
                dff = pd.concat([dff] * (len(regio_regions)))
                dff['Elem flow name'] = list_flow_added
                dff['Native geographical resolution scale'] = 'Country'

                self.master_db = pd.concat([self.master_db, dff])
                self.master_db = clean_up_dataframe(self.master_db)

    def order_things_around(self):
        """
        Simple method ordering the rows how we want them to appear, i.e., grouping midpoints together, DALY together
        and PDF together.
        :return:
        """

        self.master_db = self.master_db.sort_values(by=['Impact category', 'Elem flow name'])
        self.master_db = self.master_db.reset_index().drop('index', axis=1)

        DALY = self.master_db.loc[[i for i in self.master_db.index if self.master_db.loc[i, 'CF unit'] == 'DALY']]
        PDF = self.master_db.loc[[i for i in self.master_db.index if self.master_db.loc[i, 'CF unit'] == 'PDF.m2.yr']]
        MP = self.master_db.loc[
            [i for i in self.master_db.index if self.master_db.loc[i, 'MP or Damage'] == 'Midpoint']]

        self.master_db = clean_up_dataframe(pd.concat([MP, DALY, PDF]))

    def separate_regio_cfs(self):
        """
        Method to obtain two different versions of master_db. One with regionalized factors that will be used for
        SimaPro and openLCA versions. One with only non regionalized factors that will be used for brightway and
        ecoinvent versions.
        """

        self.master_db_not_regio = self.master_db.loc[[i for i in self.master_db.index if
                                                self.master_db.loc[i, 'Native geographical resolution scale'] not in [
                                                    'Continent', 'Country', 'Other region']]].copy()

        # dropping flow names with ", GLO" in them
        self.master_db_not_regio.drop([i for i in self.master_db_not_regio.index if ', GLO' in
                                       self.master_db_not_regio.loc[i, 'Elem flow name']], inplace=True)

    def link_to_ecoinvent(self):
        """
        Function that links names of substance from IW+ to the names of ecoinvent.
        :return: self.ei35_iw, self.ei36_iw, self.ei371_iw, self.ei38_iw
        """

        latest_ei_version = '3.10'

        ei_iw_db = self.master_db_not_regio.copy()

        # -------------- Mapping substances --------------

        mapping = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/mappings/ei'+
                                                                 latest_ei_version.replace('.','')+
                                                                 '/ei_iw_mapping.xlsx'))
        ei_mapping = mapping.loc[:, ['ecoinvent name', 'iw name']].dropna()
        not_one_for_one = ei_mapping[ei_mapping.loc[:, 'iw name'].duplicated(False)]
        one_for_one = ei_mapping[~ei_mapping.loc[:, 'iw name'].duplicated(False)]

        # for one_for_one it's easy! We just replace one for one
        ei_iw_db.loc[:, 'Elem flow name'] = ei_iw_db.loc[:, 'Elem flow name'].replace(
            one_for_one.loc[:, 'iw name'].tolist(), one_for_one.loc[:, 'ecoinvent name'].tolist())

        # for not_one_for_one it's harder, e.g., the "Zinc" substance from iw+ must be linked to multiple elementary flows in ecoinvent
        unique_not_one_for_one = set(not_one_for_one.loc[:, 'iw name'])
        for subst in unique_not_one_for_one:
            ei_df = not_one_for_one.loc[[i for i in not_one_for_one.index if not_one_for_one.loc[i, 'iw name'] == subst]]
            iw_df = ei_iw_db.loc[[i for i in ei_iw_db.index if ei_iw_db.loc[i, 'Elem flow name'] == subst]]
            new_df = pd.concat([iw_df] * len(ei_df))
            new_df = new_df.reset_index().drop('index', axis=1)
            for i, new_name in enumerate(ei_df.loc[:, 'ecoinvent name']):
                new_df.loc[len(iw_df) * i:len(iw_df) * (i + 1), 'Elem flow name'] = new_name
            ei_iw_db = pd.concat([ei_iw_db, new_df])
            ei_iw_db = clean_up_dataframe(ei_iw_db)

        # remove CFs from IW for substances that are not in ecoinvent
        ei_iw_db = ei_iw_db.loc[[i for i in ei_iw_db.index if
                                 ei_iw_db.loc[i, 'Elem flow name'] in ei_mapping.loc[:, 'ecoinvent name'].tolist()]]
        # CFs for minerals should only be for Mineral resources use and Fossil and nuclear energy use impact categories
        minerals = [i for i in ei_iw_db.index if (', in ground' in ei_iw_db.loc[i, 'Elem flow name'] and
                                                  ei_iw_db.loc[i, 'Impact category'] not in [
                                                      'Fossil and nuclear energy use',
                                                      'Mineral resources use'] and
                                                  'Water' not in ei_iw_db.loc[i, 'Elem flow name'])]
        ei_iw_db.drop(minerals, axis=0, inplace=True)
        # ions are only available in Water compartments! So remove those ions in air that don't make any sense.
        ions = [i for i in ei_iw_db.index if (', ion' in ei_iw_db.loc[i, 'Elem flow name'] and
                                              ei_iw_db.loc[i, 'Compartment'] != 'Water')]
        ei_iw_db.drop(ions, axis=0, inplace=True)

        # clean-up
        ei_iw_db = clean_up_dataframe(ei_iw_db)

        # sort
        ei_iw_db = ei_iw_db.sort_values(by=['Impact category', 'Elem flow name'])
        ei_iw_db = ei_iw_db.reset_index().drop('index', axis=1)

        # --------- Comp & subcomp shenanigans ------------

        with open(pkg_resources.resource_filename(__name__, "Data/mappings/ei"+
                                                                 latest_ei_version.replace('.','')+
                                                                 "/comps.json"), "r") as f:
            comps = json.load(f)
        with open(pkg_resources.resource_filename(__name__, "Data/mappings/ei"+
                                                                 latest_ei_version.replace('.','')+
                                                                 "/subcomps.json"), "r") as f:
            subcomps = json.load(f)

        ei_iw_db.Compartment = [comps[i] for i in ei_iw_db.Compartment]
        ei_iw_db.loc[:, 'Sub-compartment'] = [subcomps[i] if i in subcomps else None for i in
                                                ei_iw_db.loc[:, 'Sub-compartment']]

        # special cases: forestry subcomp = unspecified subcomp
        df = ei_iw_db.loc[
            [i for i in ei_iw_db.index if (ei_iw_db.loc[i, 'Sub-compartment'] == 'unspecified' and
                                           ei_iw_db.loc[i, 'Compartment'] == 'soil')]].copy()
        df.loc[:, 'Sub-compartment'] = 'forestry'
        ei_iw_db = pd.concat([ei_iw_db, df])
        ei_iw_db = clean_up_dataframe(ei_iw_db)

        # special cases: fossil well subcomp in water comp = ground- subcomp
        df = ei_iw_db.loc[[i for i in ei_iw_db.index if (ei_iw_db.loc[i, 'Sub-compartment'] == 'ground-' and
                                                             ei_iw_db.loc[i, 'Compartment'] == 'water')]].copy()
        df.loc[:, 'Sub-compartment'] = 'fossil well'
        ei_iw_db = pd.concat([ei_iw_db, df])
        ei_iw_db = clean_up_dataframe(ei_iw_db)

        # special cases: fossil well subcomp in raw comp = in ground subcomp
        df = ei_iw_db.loc[[i for i in ei_iw_db.index if (ei_iw_db.loc[i, 'Sub-compartment'] == 'in ground' and
                                                             ei_iw_db.loc[
                                                                 i, 'Compartment'] == 'natural resource')]].copy()
        df.loc[:, 'Sub-compartment'] = 'fossil well'
        ei_iw_db = pd.concat([ei_iw_db, df])
        ei_iw_db = clean_up_dataframe(ei_iw_db)

        # ----------- Unit shenanigans -------------
        ei_iw_db.loc[
            [i for i in ei_iw_db.index if ei_iw_db.loc[i, 'Elem flow unit'] == 'Bq'], 'CF value'] *= 1000
        ei_iw_db.loc[
            [i for i in ei_iw_db.index if ei_iw_db.loc[i, 'Elem flow unit'] == 'Bq'], 'Elem flow unit'] = 'kBq'
        ei_iw_db.loc[[i for i in ei_iw_db.index if
                        ei_iw_db.loc[i, 'Elem flow unit'] == 'm2.yr'], 'Elem flow unit'] = 'm2*year'

        # ---------- Odd cases ------------

        # some weird subcomps for energy flows
        energies = ['Energy, potential (in hydropower reservoir), converted',
                    'Energy, solar, converted',
                    'Energy, kinetic (in wind), converted',
                    'Energy, gross calorific value, in biomass, primary forest',
                    'Energy, gross calorific value, in biomass',
                    'Energy, geothermal, converted']
        subcomps_to_add = ['in air', 'in water', 'biotic']
        # does not make sense but it does not matter because it's corrected by linking to ecoinvent afterwards
        for energy in energies:
            for subcomp in subcomps_to_add:
                add = pd.DataFrame(
                    ['Fossil and nuclear energy use', 'MJ deprived', 'natural resource', subcomp, energy, None, 0.0,
                     'MJ', 'Midpoint', 'Global'],
                    ['Impact category', 'CF unit', 'Compartment', 'Sub-compartment', 'Elem flow name', 'CAS number',
                     'CF value','Elem flow unit', 'MP or Damage', 'Native geographical resolution scale']).T
                ei_iw_db = pd.concat([ei_iw_db, add])
                ei_iw_db = clean_up_dataframe(ei_iw_db)

        # fix to 0 "from soil/biomass" flows from long term impact categories
        ei_iw_db.loc[[i for i in ei_iw_db.index if (
                "soil or biomass" in ei_iw_db.loc[i, 'Elem flow name'] and
                ei_iw_db.loc[i, 'Impact category'] in [
                    'Climate change, ecosystem quality, long term',
                    'Climate change, human health, long term',
                    'Marine acidification, long term'])], 'CF value'] = 0
        # same for "methane, from soil or biomass" from resources
        ei_iw_db.loc[[i for i in ei_iw_db.index if (
                "Methane, from soil or biomass stock" == ei_iw_db.loc[i, 'Elem flow name'] and
                ei_iw_db.loc[i, 'Impact category'] == "Fossil and nuclear energy use")], 'CF value'] = 0

        # start with latest available version of ecoinvent
        self.ei310_iw = ei_iw_db.copy('deep')

        only_in_310 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == '3.10'].dropna(
            subset=['iw name']).loc[:, 'ecoinvent name'])

        self.ei39_iw = self.ei310_iw.drop([i for i in self.ei310_iw.index if self.ei310_iw.loc[i, 'Elem flow name'] in
                                          only_in_310]).copy('deep')

        only_in_39 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == 3.9].dropna(
            subset=['iw name']).loc[:, 'ecoinvent name'])

        self.ei38_iw = self.ei39_iw.drop([i for i in self.ei39_iw.index if self.ei39_iw.loc[i, 'Elem flow name'] in
                                          only_in_39]).copy('deep')

        only_in_38 = list(
            mapping[mapping.loc[:, 'introduced in ei v.'] == 3.8].dropna(subset=['iw name']).loc[:, 'ecoinvent name'])

        self.ei371_iw = self.ei38_iw.drop([i for i in self.ei38_iw.index if self.ei38_iw.loc[i, 'Elem flow name'] in
                                           only_in_38]).copy('deep')

        only_in_371 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == '3.7.1'].dropna(subset=['iw name']).loc[:,
                           'ecoinvent name'])

        self.ei36_iw = self.ei371_iw.drop([i for i in self.ei371_iw.index if self.ei371_iw.loc[i, 'Elem flow name'] in
                                           only_in_371]).copy('deep')

        only_in_36 = list(
            mapping[mapping.loc[:, 'introduced in ei v.'] == 3.6].dropna(subset=['iw name']).loc[:, 'ecoinvent name'])

        self.ei35_iw = self.ei36_iw.drop([i for i in self.ei36_iw.index if self.ei36_iw.loc[i, 'Elem flow name'] in
                                          only_in_36]).copy('deep')

        # ---------------------------- MATRIX VERSIONS -------------------------------

        # introducing UUID for stressors of ecoinvent
        stressors_ei35 = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/metadata/ei35/stressors.xlsx'))
        stressors_ei35 = stressors_ei35.set_index(['name', 'unit', 'comp', 'subcomp']).drop('cas', axis=1)
        stressors_ei35.index.names = (None, None, None, None)
        stressors_ei36 = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/metadata/ei36/stressors.xlsx'))
        stressors_ei36 = stressors_ei36.set_index(['name', 'unit', 'comp', 'subcomp']).drop('cas', axis=1)
        stressors_ei36.index.names = (None, None, None, None)
        stressors_ei371 = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/metadata/ei371/stressors.xlsx'))
        stressors_ei371 = stressors_ei371.set_index(['name', 'unit', 'comp', 'subcomp']).drop('cas', axis=1)
        stressors_ei371.index.names = (None, None, None, None)
        stressors_ei38 = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/metadata/ei38/stressors.xlsx'))
        stressors_ei38 = stressors_ei38.set_index(['name', 'unit', 'comp', 'subcomp']).drop('cas', axis=1)
        stressors_ei38.index.names = (None, None, None, None)
        stressors_ei39 = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/metadata/ei39/stressors.xlsx'))
        stressors_ei39 = stressors_ei39.set_index(['name', 'unit', 'comp', 'subcomp']).drop('cas', axis=1)
        stressors_ei39.index.names = (None, None, None, None)
        stressors_ei310 = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/metadata/ei310/stressors.xlsx'))
        stressors_ei310 = stressors_ei310.set_index(['name', 'unit', 'comp', 'subcomp']).drop('cas', axis=1)
        stressors_ei310.index.names = (None, None, None, None)

        df_iw = ei_iw_db.set_index(['Elem flow name', 'Elem flow unit', 'Compartment', 'Sub-compartment'])
        df_iw.index.names = (None, None, None, None)

        # matching with uuids
        df_ei35 = stressors_ei35.join(df_iw).dropna(subset=['Impact category', 'CF value'])
        df_ei35 = df_ei35.set_index('id').drop(['CAS number', 'MP or Damage',
                                                'Native geographical resolution scale'], axis=1)
        df_ei36 = stressors_ei36.join(df_iw).dropna(subset=['Impact category', 'CF value'])
        df_ei36 = df_ei36.set_index('id').drop(['CAS number', 'MP or Damage',
                                                'Native geographical resolution scale'], axis=1)
        df_ei371 = stressors_ei371.join(df_iw).dropna(subset=['Impact category', 'CF value'])
        df_ei371 = df_ei371.set_index('id').drop(['CAS number', 'MP or Damage',
                                                  'Native geographical resolution scale'], axis=1)
        df_ei38 = stressors_ei38.join(df_iw).dropna(subset=['Impact category', 'CF value'])
        df_ei38 = df_ei38.set_index('id').drop(['CAS number', 'MP or Damage',
                                                'Native geographical resolution scale'], axis=1)
        df_ei39 = stressors_ei39.join(df_iw).dropna(subset=['Impact category', 'CF value'])
        df_ei39 = df_ei39.set_index('id').drop(['CAS number', 'MP or Damage',
                                                'Native geographical resolution scale'], axis=1)
        df_ei310 = stressors_ei310.join(df_iw).dropna(subset=['Impact category', 'CF value'])
        df_ei310 = df_ei310.set_index('id').drop(['CAS number', 'MP or Damage',
                                                  'Native geographical resolution scale'], axis=1)
        # changing into a dataframe format readily available for matrix calculations
        self.ei35_iw_as_matrix = df_ei35.pivot_table(values='CF value', index='id',
                                                     columns=['Impact category', 'CF unit']).fillna(0)
        self.ei35_iw_as_matrix = self.ei35_iw_as_matrix.reindex(stressors_ei35.id).fillna(0)
        self.ei36_iw_as_matrix = df_ei36.pivot_table(values='CF value', index='id',
                                                     columns=['Impact category', 'CF unit']).fillna(0)
        self.ei36_iw_as_matrix = self.ei36_iw_as_matrix.reindex(stressors_ei36.id).fillna(0)
        self.ei371_iw_as_matrix = df_ei371.pivot_table(values='CF value', index='id',
                                                       columns=['Impact category', 'CF unit']).fillna(0)
        self.ei371_iw_as_matrix = self.ei371_iw_as_matrix.reindex(stressors_ei371.id).fillna(0)
        self.ei38_iw_as_matrix = df_ei38.pivot_table(values='CF value', index='id',
                                                     columns=['Impact category', 'CF unit']).fillna(0)
        self.ei38_iw_as_matrix = self.ei38_iw_as_matrix.reindex(stressors_ei38.id).fillna(0)
        self.ei39_iw_as_matrix = df_ei39.pivot_table(values='CF value', index='id',
                                                     columns=['Impact category', 'CF unit']).fillna(0)
        self.ei39_iw_as_matrix = self.ei39_iw_as_matrix.reindex(stressors_ei39.id).fillna(0)
        self.ei310_iw_as_matrix = df_ei310.pivot_table(values='CF value', index='id',
                                                       columns=['Impact category', 'CF unit']).fillna(0)
        self.ei310_iw_as_matrix = self.ei310_iw_as_matrix.reindex(stressors_ei310.id).fillna(0)

    def link_to_sp(self):
        """
        This method creates a SimaPro method with the IW+ characterization factors.
        :return:
        """

        # copy to not make changes on the original
        self.iw_sp = self.master_db.copy()

        # need to change the unit of land occupation flows to match SP nomenclature
        self.iw_sp.loc[[i for i in self.iw_sp.index if self.iw_sp.loc[i, 'Elem flow unit'] == 'm2.yr'], 'Elem flow unit'] = 'm2a'

        # SP does not like the "remote" subcomp of IW+ for particulate matter
        self.iw_sp.drop([i for i in self.iw_sp.index if self.iw_sp.loc[i, 'Sub-compartment'] == 'remote'], inplace=True)

        # Some water names are reserved names in SimaPro, so we modify it
        self.iw_sp.loc[self.iw_sp['Elem flow name'] == 'Water', 'Elem flow name'] = 'Water/m3'
        self.iw_sp.loc[self.iw_sp['Elem flow name'] == 'Water, agri', 'Elem flow name'] = 'Water/m3, agri'
        self.iw_sp.loc[self.iw_sp['Elem flow name'] == 'Water, non-agri', 'Elem flow name'] = 'Water/m3, non-agri'

        # now apply the mapping with the different SP flow names
        sp = pd.read_excel(pkg_resources.resource_filename(__name__, '/Data/mappings/SP/sp_mapping.xlsx'), None)
        sp = pd.concat([sp['Non regionalized'], sp['Regions existing in SP'], sp['Regions only in IW+']]).dropna()
        sp = clean_up_dataframe(sp)
        differences = sp.loc[[i for i in sp.index if sp.loc[i, 'SimaPro flows'] != sp.loc[i, 'IW+ flows']]]
        for diff in differences.index:
            if '%' not in sp.loc[diff, 'SimaPro flows']:
                df = self.iw_sp.loc[self.iw_sp.loc[:, 'Elem flow name'] == sp.loc[diff, 'IW+ flows']].copy()
                df.loc[:, 'Elem flow name'] = sp.loc[diff, 'SimaPro flows']
                self.iw_sp = pd.concat([self.iw_sp, df])
            # special case for minerals in ground, need to only apply their CFs to the Mineral resources use category
            else:
                df = self.iw_sp.loc[self.iw_sp.loc[:, 'Elem flow name'] == sp.loc[diff, 'IW+ flows']].copy()
                df = df[df['Impact category'] == 'Mineral resources use']
                df.loc[:, 'Elem flow name'] = sp.loc[diff, 'SimaPro flows']
                self.iw_sp = pd.concat([self.iw_sp, df])
        self.iw_sp = clean_up_dataframe(self.iw_sp)

        # need an unspecified subcomp for mineral resource uses, for some databases in SP (e.g., Industry2.0)
        df = self.iw_sp.loc[
            [i for i in self.iw_sp.index if self.iw_sp.loc[i, 'Impact category'] == 'Mineral resources use']].copy()
        df['Sub-compartment'] = '(unspecified)'
        self.iw_sp = pd.concat([self.iw_sp, df])
        self.iw_sp = clean_up_dataframe(self.iw_sp)

        # add some flows from SP that require conversions because of units
        new_flows = sp.loc[[i for i in sp.index if 'MJ' in sp.loc[i, 'SimaPro flows']], 'SimaPro flows'].values.tolist()
        for new_flow in new_flows:
            if 'Wood' not in new_flow:
                to_add = pd.DataFrame(['Fossil and nuclear energy use',
                                       'MJ deprived',
                                       'Raw',
                                       'in ground',
                                       new_flow,
                                       '',
                                       new_flow.split(' MJ')[0].split(', ')[-1],
                                       'kg',
                                       'Midpoint',
                                       'Global'], index=self.iw_sp.columns)
                self.iw_sp = pd.concat([self.iw_sp, to_add.T], axis=0)
                # same with other subcomp
                to_add = pd.DataFrame(['Fossil and nuclear energy use',
                                       'MJ deprived',
                                       'Raw',
                                       '(unspecified)',
                                       new_flow,
                                       '',
                                       new_flow.split(' MJ')[0].split(', ')[-1],
                                       'kg',
                                       'Midpoint',
                                       'Global'], index=self.iw_sp.columns)
                self.iw_sp = pd.concat([self.iw_sp, to_add.T], axis=0)
        self.iw_sp = clean_up_dataframe(self.iw_sp)
        # some other flows from SP that require conversions because of units
        new_flows = {
            'Gas, natural/kg': 'Gas, natural/m3', 'Uranium/kg': 'Uranium',
            'Gas, mine, off-gas, process, coal mining/kg':'Gas, mine, off-gas, process, coal mining/m3',
            'Water, cooling, unspecified natural origin/kg': 'Water, cooling, unspecified natural origin',
            'Water, process, unspecified natural origin/kg': 'Water, cooling, unspecified natural origin',
            'Water, unspecified natural origin/kg': 'Water, unspecified natural origin',
            'Water/kg': 'Water/m3', 'Wood, unspecified, standing/kg': 'Wood, unspecified, standing/m3'
        }
        for flow in new_flows:
            if 'Gas' in flow:
                density = 0.8  # kg/m3 (https://www.engineeringtoolbox.com/gas-density-d_158.html)
            elif 'Water' in flow:
                density = 1000  # kg/m3
            elif 'Wood' in flow:
                density = 600  # kg/m3 (https://www.engineeringtoolbox.com/wood-density-d_40.html)
            elif 'Uranium' in flow:
                density = 1  # kg/kg

            df = self.iw_sp[self.iw_sp['Elem flow name'] == new_flows[flow]].copy()
            df.loc[:, 'Elem flow name'] = flow
            df.loc[:, 'Elem flow unit'] = 'kg'
            df.loc[:, 'CF value'] /= density
            self.iw_sp = pd.concat([self.iw_sp, df])
        self.iw_sp = clean_up_dataframe(self.iw_sp)

        # some final adjustments for annoying flows...
        problems = {'Water, process, drinking': 'kg',
                    'Water, cooling, surface': 'kg',
                    'Water, process, surface': 'kg',
                    'Water, process, well': 'kg',
                    'Water, groundwater consumption': 'kg',
                    'Water, surface water consumption': 'kg',
                    'Water, Saline water consumption': 'kg'}

        for problem_child in problems:
            self.iw_sp.loc[self.iw_sp['Elem flow name'] == problem_child, 'Elem flow unit'] = 'kg'

    def link_to_olca(self):
        """
        This method creates an openLCA method with the IW+ characterization factors.
        :return:
        """

        self.olca_iw = self.master_db.copy()

        olca_mapping = pd.read_excel(pkg_resources.resource_filename(__name__, '/Data/mappings/oLCA/v1.10/oLCA_mapping.xlsx'))
        olca_flows = pd.read_excel(pkg_resources.resource_filename(__name__, '/Data/mappings/oLCA/v1.10/all_stressors.xlsx'))

        # before removing flows that are not matched, keep those related to energy that are treated later
        energy_flows = olca_mapping.loc[[i for i in olca_mapping.index if 'MJ' in olca_mapping.loc[i, 'oLCA name']],
                                        'oLCA name'].values.tolist()
        energy_flows = olca_flows.loc[[i for i in olca_flows.index if olca_flows.loc[i, 'flow_name'] in energy_flows]]

        # before removing flows that are not matched, keep those related to different units used
        flows_to_convert = {
            'Gas, natural/kg': 'Gas, natural/m3',
            'Gas, mine, off-gas, process, coal mining/kg': 'Gas, mine, off-gas, process, coal mining/m3',
            'Water, cooling, unspecified natural origin/kg': 'Water, cooling, unspecified natural origin',
            'Water, process, unspecified natural origin/kg': 'Water, cooling, unspecified natural origin',
            'Water, unspecified natural origin/kg': 'Water, unspecified natural origin',
            'Wood, unspecified, standing/kg': 'Wood, unspecified, standing/m3'
        }
        converting = olca_flows.loc[
            [i for i in olca_flows.index if olca_flows.loc[i, 'flow_name'] in flows_to_convert.keys()]]
        converting['reference'] = [flows_to_convert[i] for i in converting['flow_name']]
        assert converting.loc[31546, 'flow_name'] == 'Gas, mine, off-gas, process, coal mining/kg'
        converting.loc[31546, 'reference'] = 'Gas, mine, off-gas, process, coal mining'
        assert converting.loc[31563, 'flow_name'] == 'Gas, natural/kg'
        converting.loc[31563, 'reference'] = 'Gas, natural, in ground'
        assert converting.loc[60453, 'flow_name'] == 'Water, unspecified natural origin/kg'
        converting.loc[60453, 'reference'] = 'Water, unspecified natural origin/m3'

        # map oLCA and IW flow names + drop oLCA flows not linked to IW
        olca_flows = olca_flows.merge(olca_mapping, left_on=['flow_name'],
                                      right_on=['oLCA name'], how='left').dropna(subset=['iw name'])
        # split comp and subcomp
        olca_flows['Sub-compartment'] = [eval(i)[1] for i in olca_flows['comp']]
        olca_flows['Compartment'] = [eval(i)[0] for i in olca_flows['comp']]

        # map compartments between oLCA and IW
        with open(pkg_resources.resource_filename(__name__, '/Data/mappings/oLCA/v1.10/comps.json'), 'r') as f:
            comps = json.load(f)
        olca_flows['Compartment'] = [comps[i] for i in olca_flows['Compartment']]

        # map sub-compartments between oLCA and IW
        with open(pkg_resources.resource_filename(__name__, '/Data/mappings/oLCA/v1.10/subcomps.json'), 'r') as f:
            subcomps = json.load(f)

        # remove weird LT subcomps that only exist in openLCA flows
        olca_flows = olca_flows.drop([i for i in olca_flows.index if eval(olca_flows.loc[i, 'comp'])[1] in [
            'river, long-term', 'fresh water, long-term']])
        olca_flows['Sub-compartment'] = [subcomps[i] for i in olca_flows['Sub-compartment']]

        # switch iw names for olca names
        merge = olca_flows.merge(self.olca_iw, how='left',
                                 left_on=['iw name', 'Sub-compartment', 'Compartment'],
                                 right_on=['Elem flow name', 'Sub-compartment', 'Compartment'])
        # store flows of oLCA not linked yet
        only_in_olca = merge[merge.drop(['cas', 'CAS number'], axis=1).isna().any(axis=1)].copy('deep')
        # delete from original (full of NaNs)
        merge = merge.drop(merge[merge.drop(['cas', 'CAS number'], axis=1).isna().any(axis=1)].index)
        only_in_olca = only_in_olca.drop(['Impact category', 'CF unit', 'Elem flow name', 'CAS number', 'CF value',
                                          'Elem flow unit', 'MP or Damage', 'Native geographical resolution scale'],
                                         axis=1)
        # merge without the subcomp (that's what bothering)
        only_in_olca = only_in_olca.merge(self.olca_iw, how='left',
                                          left_on=['iw name', 'Compartment'],
                                          right_on=['Elem flow name', 'Compartment'])
        only_in_olca = only_in_olca.dropna(subset=['Impact category', 'CF unit', 'Elem flow name', 'CF value',
                                                   'Elem flow unit', 'MP or Damage',
                                                   'Native geographical resolution scale'])

        # specific case for radiative flows (should have proper values for ocean subcomp)
        only_in_olca = only_in_olca.drop(only_in_olca.loc[only_in_olca.unit == 'kBq'].loc[
                                             only_in_olca.comp == "('Emission to water', 'ocean')"].index)
        only_in_olca.drop('Sub-compartment_y', axis=1, inplace=True)
        only_in_olca.rename(columns={'Sub-compartment_x': 'Sub-compartment'}, inplace=True)
        self.olca_iw = clean_up_dataframe(pd.concat([merge, only_in_olca]))

        # CFs for water flows in kilograms to be converted
        self.olca_iw.loc[[i for i in self.olca_iw.index if (self.olca_iw.loc[i, 'unit'] == 'kg' and
                                                            'Water' in self.olca_iw.loc[i, 'flow_name'] and
                                                            '/kg' not in self.olca_iw.loc[
                                                                i, 'flow_name'])], 'CF value'] /= 1000

        # remove information from IW+ / only keep information relevant for oLCA
        self.olca_iw = self.olca_iw.drop(
            ['oLCA name', 'iw name', 'Sub-compartment', 'Compartment', 'Elem flow name', 'CAS number', 'Elem flow unit',
             'Native geographical resolution scale'], axis=1)
        # define unique index
        self.olca_iw = self.olca_iw.set_index(['Impact category', 'CF unit', 'flow_id'])
        # sort the multi-index
        self.olca_iw = self.olca_iw.sort_index()

        # oLCA flows are in kBq and iw+ flows in Bq. We convert CF values accordingly.
        self.olca_iw.loc[['Ionizing radiations',
                          'Ionizing radiation, ecosystem quality',
                          'Ionizing radiation, human health'], 'CF value'] *= 1000

        # flows which CF values change depending on their names, e.g., "Coal, 18MJ/kg"
        for id_ in energy_flows.index:
            if 'Wood' not in energy_flows.loc[id_, 'flow_name']:
                to_add = pd.DataFrame(['Fossil and nuclear energy use',
                                       'MJ deprived',
                                       energy_flows.loc[id_, 'flow_id'],
                                       energy_flows.loc[id_, 'flow_name'],
                                       energy_flows.loc[id_, 'comp'],
                                       energy_flows.loc[id_, 'cas'],
                                       energy_flows.loc[id_, 'unit'],
                                       energy_flows.loc[id_, 'flow_name'].split(' MJ')[0].split(', ')[-1],
                                       'Midpoint'],
                                      index=['Impact category', 'CF unit', 'flow_id', 'flow_name', 'comp', 'cas',
                                             'unit', 'CF value', 'MP or Damage'])
                to_add = to_add.T
                to_add.index = pd.MultiIndex.from_product([['Fossil and nuclear energy use'], ['MJ deprived'],
                                                           [energy_flows.loc[id_, 'flow_id']]])
                to_add.drop(['Impact category', 'CF unit', 'flow_id'], axis=1, inplace=True)

                self.olca_iw = pd.concat([self.olca_iw, to_add], axis=0)

        # flows which require conversions because of their unit, e.g., "Gas, natural/kg"
        for i in converting.index:
            df = self.olca_iw.loc[self.olca_iw['flow_name'] == converting.loc[i, 'reference']].copy()
            df.index.names = ['Impact category', 'CF unit', 'flow_id']
            df = df.loc[df['comp'] == converting.loc[i, 'comp']]
            df = df.reset_index()
            df['flow_id'] = converting.loc[i, 'flow_id']
            df['flow_name'] = converting.loc[i, 'flow_name']
            df = df.set_index(['Impact category', 'CF unit', 'flow_id'])
            if 'Gas' in converting.loc[i, 'flow_name']:
                df.loc[:, 'CF value'] /= 0.8  # kg/m3 (https://www.engineeringtoolbox.com/gas-density-d_158.html)
            elif 'Water' in converting.loc[i, 'flow_name']:
                df.loc[:, 'CF value'] /= 1000  # kg/m3
            self.olca_iw.update(df)

        self.olca_iw.index.names = ['Impact category', 'CF unit', 'flow_id']

    def link_to_exiobase(self):
        """
        This method creates an openLCA method with the IW+ characterization factors.
        :return:
        """

        EXIO_IW_concordance = pd.read_excel(pkg_resources.resource_filename(
            __name__, 'Data/mappings/exiobase/EXIO_IW_concordance.xlsx'))
        EXIO_IW_concordance.set_index('EXIOBASE', inplace=True)

        self.exio_iw = pd.DataFrame(0, EXIO_IW_concordance.index, list(set(list(zip(self.master_db.loc[:, 'Impact category'],
                                                                         self.master_db.loc[:, 'CF unit'])))))
        self.exio_iw.columns = pd.MultiIndex.from_tuples(self.exio_iw.columns, names=['Impact category', 'CF unit'])
        self.exio_iw = self.exio_iw.T.sort_index().T

        for flow in EXIO_IW_concordance.index:
            if not EXIO_IW_concordance.loc[flow].isna().iloc[0]:
                # identify all entries (any impact category, compartment, etc.) for given flow
                CF_flow = self.master_db.loc[self.master_db['Elem flow name'] == EXIO_IW_concordance.loc[flow, 'IW']].loc[:,
                          ['Impact category', 'CF unit', 'CF value', 'Compartment', 'Sub-compartment']]
                # name of the comp in lower case to match exiobase easily
                CF_flow.Compartment = [i.lower() for i in CF_flow.Compartment]
                # only keeping right compartments, always selecting unspecified sub-compartment
                if flow.split('- ')[-1] == 'soil':
                    # P - soil is not characterized in IW+
                    pass
                elif flow.split('- ')[-1] in ['air', 'water']:
                    CFs = pd.pivot(CF_flow, values='CF value', index=['Compartment', 'Sub-compartment'],
                                   columns=['Impact category', 'CF unit']).loc[
                        (flow.split('- ')[-1], '(unspecified)')].fillna(0)
                elif 'Occupation' in EXIO_IW_concordance.loc[flow, 'IW']:
                    CFs = pd.pivot(CF_flow, values='CF value', index=['Compartment', 'Sub-compartment'],
                                   columns=['Impact category', 'CF unit']).loc[('raw', 'land')].fillna(0)
                elif 'Fishery' in flow:
                    CFs = pd.pivot(CF_flow, values='CF value', index=['Compartment', 'Sub-compartment'],
                                   columns=['Impact category', 'CF unit']).loc[('raw', 'biotic')].fillna(0)
                else:
                    try:
                        CFs = pd.pivot(CF_flow, values='CF value', index=['Compartment', 'Sub-compartment'],
                                       columns=['Impact category', 'CF unit']).loc[('raw', '(unspecified)')].fillna(0)
                    except KeyError:
                        CFs = pd.pivot(CF_flow, values='CF value', index=['Compartment', 'Sub-compartment'],
                                       columns=['Impact category', 'CF unit']).loc[('raw', 'in ground')].fillna(0)
                CFs.name = flow
                # dumping the CF values in the C matrix
                self.exio_iw.update(pd.DataFrame(CFs).T)

        # EXIOBASE land occupation in km2 while IW in m2, so we convert
        self.exio_iw.loc[:, [i for i in self.exio_iw.columns if 'Land' in i[0]]] *= 1000000
        # EXIOBASE energy flows in TJ while IW in MJ, so we convert
        self.exio_iw.loc[:, 'Fossil and nuclear energy use'] = self.exio_iw.loc[:, 'Fossil and nuclear energy use'].values * 1000000
        # EXIOBASE mineral flows in kt while IW in kg, so we convert
        self.exio_iw.loc[:, 'Mineral resources use'] = self.exio_iw.loc[:, 'Mineral resources use'].values * 1000000
        # EXIOBASE water flows in Mm3 while IW in m3, so we convert
        self.exio_iw.loc[:, [i for i in self.exio_iw.columns if 'Water' in i[0]]] *= 1000000
        # more common to have impact categories as index
        self.exio_iw = self.exio_iw.T
        # keep same index format as previously
        self.exio_iw.index = [(i[0] + ' (' + i[1] + ')') for i in self.exio_iw.index]
        # HFC and PFC linked to "carbon dioxide" but should not impact marine acidification
        self.exio_iw.loc[['Marine acidification, long term (PDF.m2.yr)', 'Marine acidification, short term (PDF.m2.yr)'], [
            'HFC - air', 'PFC - air']] = 0

        self.special_case_minerals_exiobase()

    def special_case_minerals_exiobase(self):
        """
        Minerals in exiobase are accounted for as kilograms of ore, while IW+ operates with kilograms of metal. This
        function is here to take extra information from exiobase, such as average ore content, to create corresponding
        CFs for ores of exiobase.
        :return:
        """

        # loading the file with metal content information (obtained from the EXIOBASE team)
        metal_concentration_exiobase = pd.read_csv(pkg_resources.resource_filename(
            __name__, 'Data/metadata/exiobase/All_factors_applied_to_Exiobase_metals_minerals.csv'), sep=';')

        # extracting the average amount of metal per ore from this file
        average_gold_per_ore = metal_concentration_exiobase.loc[
            [i for i in metal_concentration_exiobase.index if
             'Gold' in metal_concentration_exiobase.CommodityName[i] and
             'Global average\r' in metal_concentration_exiobase.UsedComment[i]]].Concentration.mean()
        average_lead_per_ore = metal_concentration_exiobase.loc[
            [i for i in metal_concentration_exiobase.index if
             'Lead' in metal_concentration_exiobase.CommodityName[i]]].Concentration.mean()
        average_copper_per_ore = metal_concentration_exiobase.loc[
            [i for i in metal_concentration_exiobase.index if
             'opper' in metal_concentration_exiobase.CommodityName[i]]].Concentration.mean()
        average_silver_per_ore = metal_concentration_exiobase.loc[
            [i for i in metal_concentration_exiobase.index if
             'Silver' in metal_concentration_exiobase.CommodityName[i] and
             'Global average\r' in metal_concentration_exiobase.UsedComment[i]]].Concentration.mean()
        average_iron_per_ore = metal_concentration_exiobase.loc[
            [i for i in metal_concentration_exiobase.index if
             'Iron' in metal_concentration_exiobase.CommodityName[i] and
             'Global average\r' in metal_concentration_exiobase.UsedComment[i]]].Concentration.mean()
        average_nickel_per_ore = metal_concentration_exiobase.loc[
            [i for i in metal_concentration_exiobase.index if
             'Nickel' in metal_concentration_exiobase.CommodityName[i] and
             'Global average\r' in metal_concentration_exiobase.UsedComment[i]]].Concentration.mean()
        average_zinc_per_ore = metal_concentration_exiobase.loc[
            [i for i in metal_concentration_exiobase.index if
             'Zinc' in metal_concentration_exiobase.CommodityName[i]]].Concentration.mean()
        average_pgm_per_ore = metal_concentration_exiobase.loc[
            [i for i in metal_concentration_exiobase.index if
             'Platinum-group (PGM)' in metal_concentration_exiobase.CommodityName[i] and
             'By-products of other ore' in metal_concentration_exiobase.UsedComment[i]]].UsedFactor.mean()

        df = self.master_db.copy()
        df = df.set_index(['Elem flow name', 'Compartment', 'Sub-compartment'])

        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Domestic Extraction Used - Metal Ores - Gold ores'] = (
                average_gold_per_ore *  df.loc[('Gold','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Unused Domestic Extraction - Metal Ores - Gold ores'] = (
                average_gold_per_ore *  df.loc[('Gold','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Domestic Extraction Used - Metal Ores - Lead ores'] = (
                average_lead_per_ore * df.loc[('Lead','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Unused Domestic Extraction - Metal Ores - Lead ores'] = (
                average_lead_per_ore * df.loc[('Lead','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Domestic Extraction Used - Metal Ores - Copper ores'] = (
                average_copper_per_ore * df.loc[('Copper','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc[ 'Mineral resources use (kg deprived)', 'Unused Domestic Extraction - Metal Ores - Copper ores'] = (
                average_copper_per_ore * df.loc[('Copper','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Domestic Extraction Used - Metal Ores - Silver ores'] = (
                average_silver_per_ore * df.loc[('Silver','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Unused Domestic Extraction - Metal Ores - Silver ores'] = (
                average_silver_per_ore * df.loc[('Silver','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Domestic Extraction Used - Metal Ores - Iron ores'] = (
                average_iron_per_ore * df.loc[('Iron','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Unused Domestic Extraction - Metal Ores - Iron ores'] = (
                average_iron_per_ore * df.loc[('Iron','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Domestic Extraction Used - Metal Ores - Nickel ores'] = (
                average_nickel_per_ore * df.loc[('Nickel','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Unused Domestic Extraction - Metal Ores - Nickel ores'] = (
                average_nickel_per_ore * df.loc[('Nickel','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Domestic Extraction Used - Metal Ores - Zinc ores'] = (
                average_zinc_per_ore * df.loc[('Zinc','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Unused Domestic Extraction - Metal Ores - Zinc ores'] = (
                average_zinc_per_ore * df.loc[('Zinc','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Domestic Extraction Used - Metal Ores - PGM ores'] = (
                average_pgm_per_ore * df.loc[('Platinum','Raw','in ground'), 'CF value'].iloc[0] * 1000000)
        self.exio_iw.loc['Mineral resources use (kg deprived)', 'Unused Domestic Extraction - Metal Ores - PGM ores'] = (
                average_pgm_per_ore * df.loc[('Platinum','Raw','in ground'), 'CF value'].iloc[0] * 1000000)

        # loading the file describing which metals EXIOBASE includes in their other non-ferrous metals flow
        other_categories_composition = pd.read_excel(pkg_resources.resource_filename(
            __name__, 'Data/mappings/exiobase/Mineral_extension_exio_detailed_2016.xlsx'))

        # identify non ferrous metals among the list of mineral resources
        other_non_ferrous_metals_index = [i for i in other_categories_composition.index if
                                          'Other non-ferrous metal ores' in
                                          other_categories_composition.PhysicalTypeName[i]]
        # delete duplicated other non ferrous metals identified
        other_non_ferrous_metals = other_categories_composition.loc[other_non_ferrous_metals_index].groupby(
            other_categories_composition.CommodityName).head(1)
        # make a dataframe of it
        other_non_ferrous_metals = pd.DataFrame(0, other_non_ferrous_metals.CommodityName,
                                                ['Metal content', 'Ore abundance'])

        # iterate through the non-ferrous metals
        for metal in other_non_ferrous_metals.index:
            # if we can find a global average value pick it
            try:
                average_metal_content = metal_concentration_exiobase.loc[
                    [i for i in metal_concentration_exiobase.index if
                     metal in metal_concentration_exiobase.CommodityName[i] and
                     'Global average' in metal_concentration_exiobase.UsedComment[i]]].Concentration.mean()
            except TypeError:
                pass
            # if we didn't have a match
            if np.isnan(average_metal_content):
                # try with By-product of other ores
                try:
                    average_metal_content = metal_concentration_exiobase.loc[
                        [i for i in metal_concentration_exiobase.index if
                         metal in metal_concentration_exiobase.CommodityName[i] and
                         'By-product of other ores' in metal_concentration_exiobase.UsedComment[i]]].UsedFactor.mean()
                except TypeError:
                    pass
            # if for whatever reason we have a completely ridiculous concentration drop it
            if average_metal_content > 0.5:
                average_metal_content = None
            other_non_ferrous_metals.loc[metal, 'Metal content'] = average_metal_content

        # use 0.001 as default value
        other_non_ferrous_metals = other_non_ferrous_metals.fillna(0.001)

        abundance = pd.read_excel(pkg_resources.resource_filename(
            __name__, 'Data/metadata/exiobase/USGS_extraction_volumes.xlsx'), 'metals')
        abundance.set_index('Unnamed: 0', inplace=True)
        abundance /= abundance.sum()
        assert (other_non_ferrous_metals.index == abundance.index).all()
        other_non_ferrous_metals.loc[:, 'Ore abundance'] = abundance.values

        other_metal_concordance = pd.read_excel(pkg_resources.resource_filename(
            __name__, 'Data/mappings/exiobase/other_metals_matching.xlsx')).drop('comments', axis=1)
        other_metal_concordance.set_index('Unnamed: 0', inplace=True)
        other_metal_concordance.dropna(inplace=True)

        for flow in other_metal_concordance.index:
            CF_flow = self.master_db.loc[self.master_db['Elem flow name'] == other_metal_concordance.loc[flow].iloc[0]].loc[:,
                      ['Impact category', 'CF unit', 'CF value', 'Compartment', 'Sub-compartment']]
            CFs = pd.pivot(CF_flow, values='CF value', index=['Compartment', 'Sub-compartment'],
                           columns=['Impact category', 'CF unit']).loc[('Raw', 'in ground')].fillna(0)
            other_non_ferrous_metals.loc[flow, 'CF'] = CFs.loc['Mineral resources use'].iloc[0]

        CF = (other_non_ferrous_metals.loc[:, 'Metal content'] * other_non_ferrous_metals.loc[:,
                                                                 'Ore abundance'] * other_non_ferrous_metals.loc[:,
                                                                                    'CF']).sum()

        self.exio_iw.loc['Mineral resources use (kg deprived)',
                         'Domestic Extraction Used - Metal Ores - Other non-ferrous metal ores'] = CF * 1000000
        self.exio_iw.loc['Mineral resources use (kg deprived)',
                         'Unused Domestic Extraction - Metal Ores - Other non-ferrous metal ores'] = CF * 1000000

        # first we identify which mineral are part of "other minerals"
        other_minerals = other_categories_composition.loc[
            [i for i in other_categories_composition.index if
             'Other minerals' in other_categories_composition.PhysicalTypeName[i]]].groupby(
            other_categories_composition.CommodityName).head(1).drop(
            ['PhysicalTypeName', 'ProductTypeCode', 'CountryCode', 'ISOAlpha2', 'AccountingYear'], axis=1)
        other_minerals.index = other_minerals.CommodityName

        other_minerals.loc['Sillica sand (quartzsand)', 'CF'] = \
        df.loc[('Sand, quartz', 'Raw', 'in ground'), 'CF value'].iloc[0]
        other_minerals.loc['Feldspar', 'CF'] = df.loc[('Feldspar', 'Raw', 'in ground'), 'CF value'].iloc[0]
        other_minerals.loc['Graphite, natural', 'CF'] = \
        df.loc[('Metamorphous rock, graphite containing', 'Raw', 'in ground'), 'CF value'].iloc[0]
        other_minerals.loc['Magnesite', 'CF'] = df.loc[('Magnesite', 'Raw', 'in ground'), 'CF value'].iloc[0]
        other_minerals.loc['Talc (steatite, soapstone, pyrophyllite)', 'CF'] = \
        df.loc[('Talc', 'Raw', 'in ground'), 'CF value'].iloc[0]
        other_minerals.loc['Diatomite', 'CF'] = df.loc[('Diatomite', 'Raw', 'in ground'), 'CF value'].iloc[0]
        other_minerals.loc['Vermiculite', 'CF'] = df.loc[('Vermiculite', 'Raw', 'in ground'), 'CF value'].iloc[0]
        other_minerals.loc['Perlite', 'CF'] = df.loc[('Perlite', 'Raw', 'in ground'), 'CF value'].iloc[0]
        # arithmetic average because it does not matter (both are at 0)
        other_minerals.loc['Igneous rock (basalt, basaltic lava, diabase, granite, porphyry, etc.)',
                           'CF'] = (df.loc[('Basalt','Raw','in ground'), 'CF value'].iloc[0] +
                                    df.loc[('Granite','Raw','in ground'), 'CF value'].iloc[0]) / 2
        other_minerals.loc['Peat for agricultural use', 'CF'] = df.loc[('Peat', 'Raw', 'in ground'), 'CF value'].iloc[0]
        other_minerals.loc['Strontium minerals', 'CF'] = df.loc[('Strontium', 'Raw', 'in ground'), 'CF value'].iloc[0]
        other_minerals.loc['Abrasives, natural (puzzolan, pumice, volcanic cinder etc.)', 'CF'] = \
        df.loc[('Pumice', 'Raw', 'in ground'), 'CF value'].iloc[0]
        other_minerals.loc['Calcite', 'CF'] = df.loc[('Calcite', 'Raw', 'in ground'), 'CF value'].iloc[0]

        abundance_minerals = pd.read_excel(pkg_resources.resource_filename(
            __name__, 'Data/metadata/exiobase/USGS_extraction_volumes.xlsx'), 'minerals')

        # Include those in the dataframe containing all intel on other minerals
        other_minerals = pd.concat([other_minerals, abundance_minerals.set_index('CommodityName')], axis=1)
        # weighted average
        other_minerals.loc[:, 'Mineral extracted (kg)'] /= other_minerals.loc[:, 'Mineral extracted (kg)'].sum()
        # the CF for the aggregated mineral flow is thus:
        new_CF = (other_minerals.loc[:, 'CF'] * other_minerals.loc[:, 'Mineral extracted (kg)']).sum()

        self.exio_iw.loc['Mineral resources use (kg deprived)',
                         'Domestic Extraction Used - Non-Metallic Minerals - Other minerals'] = new_CF * 1000000
        self.exio_iw.loc[ 'Mineral resources use (kg deprived)',
                          'Unused Domestic Extraction - Non-Metallic Minerals - Other minerals'] = new_CF * 1000000

    def get_simplified_versions(self, ei_flows_version=None):

        # SimaPro
        self.simplified_version_sp = clean_up_dataframe(produce_simplified_version(self.iw_sp).reindex(
            self.iw_sp.columns, axis=1))

        # openLCA
        self.simplified_version_olca = produce_simplified_version_olca(self.olca_iw)
        assert not self.simplified_version_olca.index.duplicated().any()

        # ecoinvent
        self.simplified_version_ei35 = clean_up_dataframe(produce_simplified_version(self.ei35_iw).reindex(
            self.ei35_iw.columns, axis=1))

        self.simplified_version_ei36 = clean_up_dataframe(produce_simplified_version(self.ei36_iw).reindex(
            self.ei36_iw.columns, axis=1))

        self.simplified_version_ei371 = clean_up_dataframe(produce_simplified_version(self.ei371_iw).reindex(
            self.ei371_iw.columns, axis=1))

        self.simplified_version_ei38 = clean_up_dataframe(produce_simplified_version(self.ei38_iw).reindex(
            self.ei38_iw.columns, axis=1))

        self.simplified_version_ei39 = clean_up_dataframe(produce_simplified_version(self.ei39_iw).reindex(
            self.ei39_iw.columns, axis=1))

        self.simplified_version_ei310 = clean_up_dataframe(produce_simplified_version(self.ei310_iw).reindex(
            self.ei310_iw.columns, axis=1))

        # brightway2
        if ei_flows_version == '3.5':
            self.simplified_version_bw = self.simplified_version_ei35.copy('deep')
        elif ei_flows_version == '3.6':
            self.simplified_version_bw = self.simplified_version_ei36.copy('deep')
        elif ei_flows_version == '3.7.1':
            self.simplified_version_bw = self.simplified_version_ei371.copy('deep')
        elif ei_flows_version == '3.8':
            self.simplified_version_bw = self.simplified_version_ei38.copy('deep')
        elif ei_flows_version == '3.9.1':
            self.simplified_version_bw = self.simplified_version_ei39.copy('deep')
        else:
            self.simplified_version_bw = self.simplified_version_ei39.copy('deep')

    def get_total_hh_and_eq(self):
        """
        OpenLCA doesn't allow for reliable contribution analyses for total damage categories (unlike
        SimaPro). So we create two additional impact categories "Total human health" and "Total ecosystem quality".
        :return:
        """

        total_hh = self.olca_iw.loc(axis=0)[:, 'DALY'].groupby(axis=0, level=2).agg({'flow_name': 'first',
                                                                                     'comp': 'first',
                                                                                     'cas': 'first',
                                                                                     'unit': 'first',
                                                                                     'CF value': sum,
                                                                                     'MP or Damage': 'first'})
        total_hh.index = pd.MultiIndex.from_tuples([('Total human health', 'DALY', i) for i in total_hh.index])

        total_eq = self.olca_iw.loc(axis=0)[:, 'PDF.m2.yr'].groupby(axis=0, level=2).agg({'flow_name': 'first',
                                                                                          'comp': 'first',
                                                                                          'cas': 'first',
                                                                                          'unit': 'first',
                                                                                          'CF value': sum,
                                                                                          'MP or Damage': 'first'})
        total_eq.index = pd.MultiIndex.from_tuples([('Total ecosystem quality', 'PDF.m2.yr', i) for i in total_eq.index])

        self.olca_iw = pd.concat([self.olca_iw, total_hh, total_eq])


# -------------- Support modules -------------------

MCO2 = 44.01
MCH4 = 16.043
MC = 12.011
Y = 0.75
tauOH = 9.7

def AGTPCO2(t, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3, kPulseT, aT1, tauT1, aT2, tauT2, AACO2):
    term1 = aC1 * aT1 * (1 - np.exp(-t / tauT1)) + aC1 * aT2 * (1 - np.exp(-t / tauT2))
    term2 = (aC2 * aT1 * (np.exp(-t / tauT1) - np.exp(-t / tauC1)) * tauC1) / (tauT1 - tauC1)
    term3 = (aC2 * aT2 * (np.exp(-t / tauT2) - np.exp(-t / tauC1)) * tauC1) / (tauT2 - tauC1)
    term4 = (aC3 * aT1 * (np.exp(-t / tauT1) - np.exp(-t / tauC2)) * tauC2) / (tauT1 - tauC2)
    term5 = (aC3 * aT2 * (np.exp(-t / tauT2) - np.exp(-t / tauC2)) * tauC2) / (tauT2 - tauC2)
    term6 = (aC4 * aT1 * (np.exp(-t / tauT1) - np.exp(-t / tauC3)) * tauC3) / (tauT1 - tauC3)
    term7 = (aC4 * aT2 * (np.exp(-t / tauT2) - np.exp(-t / tauC3)) * tauC3) / (tauT2 - tauC3)

    return AACO2 * kPulseT * (term1 + term2 + term3 + term4 + term5 + term6 + term7)


def AGTPNonCO2(t, tauNonCO2, kPulseT, aT1, tauT1, aT2, tauT2, AANonCO2):
    term1 = (aT1 * (np.exp(-t / tauT1) - np.exp(-t / tauNonCO2))) / (tauT1 - tauNonCO2)
    term2 = (aT2 * (np.exp(-t / tauT2) - np.exp(-t / tauNonCO2))) / (tauT2 - tauNonCO2)

    return AANonCO2 * kPulseT * (term1 + term2) * tauNonCO2


def DAGTPNonCO2(t, tauNonCO2, kPulseT, aT1, tauT1, aT2, tauT2, AANonCO2, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2,
                gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3):
    term1 = -aC1 * aS1 * aT1 ** 2 * tauT1 ** 2 * tauNonCO2 * (-np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                (tauT1 - tauNonCO2) ** 2 * (tauT1 - tauS1)) + aC1 * aS1 * aT1 ** 2 * tauT1 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS1)) + aC1 * aS1 * aT1 ** 2 * tauT1 * t * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS1)) + aC1 * aS1 * aT1 ** 2 * tauT1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauNonCO2) ** 2 * (tauNonCO2 - tauS1)) + aC1 * aS1 * aT1 ** 2 * tauT1 * (
                        -tauT1 * tauS1 + (tauT1 * (t + tauS1) - t * tauS1) * np.exp(
                    t * (-1 / tauS1 + 1 / tauT1))) * np.exp(-t / tauT1) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS1) ** 2) - aC1 * aS1 * aT1 ** 2 * tauT1 * (
                        tauS1 - tauS1 * np.exp(-t / tauS1)) / (
                        tauS1 * (tauT1 - tauNonCO2)) - aC1 * aS1 * aT1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauNonCO2) * (tauNonCO2 - tauS1)) + aC1 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauT2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauT2) * (tauT1 - tauS1) * (
                tauT2 - tauNonCO2)) + aC1 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauT2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS1)) + aC1 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauNonCO2) * (tauT1 - tauS1) * (
                -tauT2 + tauNonCO2)) + aC1 * aS1 * aT1 * tauT1 ** 2 * aT2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS1)) - aC1 * aS1 * aT1 * tauT1 * aT2 * tauT2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / ((tauT1 - tauT2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS1)) - aC1 * aS1 * aT1 * tauT1 * aT2 * tauT2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / ((tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                tauT2 - tauS1)) - aC1 * aS1 * aT1 * tauT1 * aT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                tauNonCO2 - tauS1)) - aC1 * aS1 * aT1 * tauT1 * aT2 * (1 - np.exp(-t / tauS1)) / (
                        tauT1 - tauNonCO2) + aC1 * aS1 * aT1 * aT2 * tauT2 ** 2 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / ((-tauT1 + tauNonCO2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS1)) + aC1 * aS1 * aT1 * aT2 * tauT2 ** 2 * (-np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS1)) - aC1 * aS1 * aT1 * aT2 * tauT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / ((-tauT1 + tauNonCO2) * (tauT2 - tauNonCO2) * (
                tauNonCO2 - tauS1)) - aC1 * aS1 * aT1 * aT2 * tauT2 * (1 - np.exp(-t / tauS1)) / (
                        tauT2 - tauNonCO2) - aC1 * aS1 * aT1 * aT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauNonCO2) * (tauNonCO2 - tauS1)) - aC1 * aS1 * aT1 * aT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauNonCO2) * (tauNonCO2 - tauS1)) + aC1 * aS1 * aT1 * aT2 * tauNonCO2 * (
                        1 - np.exp(-t / tauS1)) / (tauT2 - tauNonCO2) + aC1 * aS1 * aT1 * aT2 * tauNonCO2 * (
                        1 - np.exp(-t / tauS1)) / (
                        tauT1 - tauNonCO2) - aC1 * aS1 * aT2 ** 2 * tauT2 ** 2 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) ** 2 * (tauT2 - tauS1)) + aC1 * aS1 * aT2 ** 2 * tauT2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS1)) + aC1 * aS1 * aT2 ** 2 * tauT2 * t * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS1)) + aC1 * aS1 * aT2 ** 2 * tauT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauNonCO2) ** 2 * (tauNonCO2 - tauS1)) + aC1 * aS1 * aT2 ** 2 * tauT2 * (
                        -tauT2 * tauS1 + (tauT2 * (t + tauS1) - t * tauS1) * np.exp(
                    t * (-1 / tauS1 + 1 / tauT2))) * np.exp(-t / tauT2) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS1) ** 2) - aC1 * aS1 * aT2 ** 2 * tauT2 * (
                        tauS1 - tauS1 * np.exp(-t / tauS1)) / (
                        tauS1 * (tauT2 - tauNonCO2)) - aC1 * aS1 * aT2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauNonCO2) * (tauNonCO2 - tauS1)) - aC1 * aS2 * aT1 ** 2 * tauT1 ** 2 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) ** 2 * (tauT1 - tauS2)) + aC1 * aS2 * aT1 ** 2 * tauT1 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS2)) + aC1 * aS2 * aT1 ** 2 * tauT1 * t * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS2)) + aC1 * aS2 * aT1 ** 2 * tauT1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauNonCO2) ** 2 * (tauNonCO2 - tauS2)) + aC1 * aS2 * aT1 ** 2 * tauT1 * (
                        -tauT1 * tauS2 + (tauT1 * (t + tauS2) - t * tauS2) * np.exp(
                    t * (-1 / tauS2 + 1 / tauT1))) * np.exp(-t / tauT1) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS2) ** 2) - aC1 * aS2 * aT1 ** 2 * tauT1 * (
                        tauS2 - tauS2 * np.exp(-t / tauS2)) / (
                        tauS2 * (tauT1 - tauNonCO2)) - aC1 * aS2 * aT1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauNonCO2) * (tauNonCO2 - tauS2)) + aC1 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauT2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauT2) * (tauT1 - tauS2) * (
                tauT2 - tauNonCO2)) + aC1 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauT2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS2)) + aC1 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauNonCO2) * (tauT1 - tauS2) * (
                -tauT2 + tauNonCO2)) + aC1 * aS2 * aT1 * tauT1 ** 2 * aT2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS2)) - aC1 * aS2 * aT1 * tauT1 * aT2 * tauT2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT1 - tauT2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS2)) - aC1 * aS2 * aT1 * tauT1 * aT2 * tauT2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                tauT2 - tauS2)) - aC1 * aS2 * aT1 * tauT1 * aT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                tauNonCO2 - tauS2)) - aC1 * aS2 * aT1 * tauT1 * aT2 * (1 - np.exp(-t / tauS2)) / (tauT1 - tauNonCO2)

    term2 = aC1 * aS2 * aT1 * aT2 * tauT2 ** 2 * tauNonCO2 * (-np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                (-tauT1 + tauNonCO2) * (tauT2 - tauNonCO2) * (tauT2 - tauS2)) + aC1 * aS2 * aT1 * aT2 * tauT2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS2)) - aC1 * aS2 * aT1 * aT2 * tauT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / ((-tauT1 + tauNonCO2) * (tauT2 - tauNonCO2) * (
                tauNonCO2 - tauS2)) - aC1 * aS2 * aT1 * aT2 * tauT2 * (1 - np.exp(-t / tauS2)) / (
                        tauT2 - tauNonCO2) - aC1 * aS2 * aT1 * aT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauNonCO2) * (tauNonCO2 - tauS2)) - aC1 * aS2 * aT1 * aT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauNonCO2) * (tauNonCO2 - tauS2)) + aC1 * aS2 * aT1 * aT2 * tauNonCO2 * (
                        1 - np.exp(-t / tauS2)) / (tauT2 - tauNonCO2) + aC1 * aS2 * aT1 * aT2 * tauNonCO2 * (
                        1 - np.exp(-t / tauS2)) / (
                        tauT1 - tauNonCO2) - aC1 * aS2 * aT2 ** 2 * tauT2 ** 2 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) ** 2 * (tauT2 - tauS2)) + aC1 * aS2 * aT2 ** 2 * tauT2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS2)) + aC1 * aS2 * aT2 ** 2 * tauT2 * t * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS2)) + aC1 * aS2 * aT2 ** 2 * tauT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauNonCO2) ** 2 * (tauNonCO2 - tauS2)) + aC1 * aS2 * aT2 ** 2 * tauT2 * (
                        -tauT2 * tauS2 + (tauT2 * (t + tauS2) - t * tauS2) * np.exp(
                    t * (-1 / tauS2 + 1 / tauT2))) * np.exp(-t / tauT2) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS2) ** 2) - aC1 * aS2 * aT2 ** 2 * tauT2 * (
                        tauS2 - tauS2 * np.exp(-t / tauS2)) / (
                        tauS2 * (tauT2 - tauNonCO2)) - aC1 * aS2 * aT2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauNonCO2) * (tauNonCO2 - tauS2)) - aC1 * aS3 * aT1 ** 2 * tauT1 ** 2 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) ** 2 * (tauT1 - tauS3)) + aC1 * aS3 * aT1 ** 2 * tauT1 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS3)) + aC1 * aS3 * aT1 ** 2 * tauT1 * t * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS3)) + aC1 * aS3 * aT1 ** 2 * tauT1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauNonCO2) ** 2 * (tauNonCO2 - tauS3)) + aC1 * aS3 * aT1 ** 2 * tauT1 * (
                        -tauT1 * tauS3 + (tauT1 * (t + tauS3) - t * tauS3) * np.exp(
                    t * (-1 / tauS3 + 1 / tauT1))) * np.exp(-t / tauT1) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS3) ** 2) - aC1 * aS3 * aT1 ** 2 * tauT1 * (
                        tauS3 - tauS3 * np.exp(-t / tauS3)) / (
                        tauS3 * (tauT1 - tauNonCO2)) - aC1 * aS3 * aT1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauNonCO2) * (tauNonCO2 - tauS3)) + aC1 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauT2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauT2) * (tauT1 - tauS3) * (
                tauT2 - tauNonCO2)) + aC1 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauT2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS3)) + aC1 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauNonCO2) * (tauT1 - tauS3) * (
                -tauT2 + tauNonCO2)) + aC1 * aS3 * aT1 * tauT1 ** 2 * aT2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauNonCO2) * (tauT1 - tauS3)) - aC1 * aS3 * aT1 * tauT1 * aT2 * tauT2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / ((tauT1 - tauT2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS3)) - aC1 * aS3 * aT1 * tauT1 * aT2 * tauT2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / ((tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                tauT2 - tauS3)) - aC1 * aS3 * aT1 * tauT1 * aT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                tauNonCO2 - tauS3)) - aC1 * aS3 * aT1 * tauT1 * aT2 * (1 - np.exp(-t / tauS3)) / (
                        tauT1 - tauNonCO2) + aC1 * aS3 * aT1 * aT2 * tauT2 ** 2 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / ((-tauT1 + tauNonCO2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS3)) + aC1 * aS3 * aT1 * aT2 * tauT2 ** 2 * (-np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS3)) - aC1 * aS3 * aT1 * aT2 * tauT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / ((-tauT1 + tauNonCO2) * (tauT2 - tauNonCO2) * (
                tauNonCO2 - tauS3)) - aC1 * aS3 * aT1 * aT2 * tauT2 * (1 - np.exp(-t / tauS3)) / (
                        tauT2 - tauNonCO2) - aC1 * aS3 * aT1 * aT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauNonCO2) * (tauNonCO2 - tauS3)) - aC1 * aS3 * aT1 * aT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauNonCO2) * (tauNonCO2 - tauS3)) + aC1 * aS3 * aT1 * aT2 * tauNonCO2 * (
                        1 - np.exp(-t / tauS3)) / (tauT2 - tauNonCO2) + aC1 * aS3 * aT1 * aT2 * tauNonCO2 * (
                        1 - np.exp(-t / tauS3)) / (
                        tauT1 - tauNonCO2) - aC1 * aS3 * aT2 ** 2 * tauT2 ** 2 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) ** 2 * (tauT2 - tauS3)) + aC1 * aS3 * aT2 ** 2 * tauT2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS3)) + aC1 * aS3 * aT2 ** 2 * tauT2 * t * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS3)) + aC1 * aS3 * aT2 ** 2 * tauT2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauNonCO2) ** 2 * (tauNonCO2 - tauS3)) + aC1 * aS3 * aT2 ** 2 * tauT2 * (
                        -tauT2 * tauS3 + (tauT2 * (t + tauS3) - t * tauS3) * np.exp(
                    t * (-1 / tauS3 + 1 / tauT2))) * np.exp(-t / tauT2) / (
                        (tauT2 - tauNonCO2) * (tauT2 - tauS3) ** 2) - aC1 * aS3 * aT2 ** 2 * tauT2 * (
                        tauS3 - tauS3 * np.exp(-t / tauS3)) / (
                        tauS3 * (tauT2 - tauNonCO2)) - aC1 * aS3 * aT2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauNonCO2) * (tauNonCO2 - tauS3)) - aC1 * aT1 ** 2 * tauT1 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / (tauT1 - tauNonCO2) ** 2 + aC1 * aT1 ** 2 * tauT1 * tauNonCO2 * np.exp(-t / tauT1) / (
                        tauT1 - tauNonCO2) ** 2 + aC1 * aT1 ** 2 * tauT1 / (
                        tauT1 - tauNonCO2) - aC1 * aT1 ** 2 * tauT1 * np.exp(-t / tauT1) / (
                        tauT1 - tauNonCO2) - aC1 * aT1 ** 2 * t * np.exp(-t / tauT1) / (
                        tauT1 - tauNonCO2) - aC1 * aT1 ** 2 * tauNonCO2 / (
                        tauT1 - tauNonCO2) + aC1 * aT1 ** 2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        tauT1 - tauNonCO2) + aC1 * aT1 * tauT1 * aT2 * tauT2 * np.exp(-t / tauT2) / (
                        (tauT1 - tauT2) * (tauT2 - tauNonCO2)) - aC1 * aT1 * tauT1 * aT2 * tauT2 * np.exp(
        -t / tauT1) / ((tauT1 - tauT2) * (tauT2 - tauNonCO2)) + aC1 * aT1 * tauT1 * aT2 * tauT2 * np.exp(-t / tauT2) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2)) - aC1 * aT1 * tauT1 * aT2 * tauT2 * np.exp(
        -t / tauT1) / ((tauT1 - tauT2) * (tauT1 - tauNonCO2)) + aC1 * aT1 * tauT1 * aT2 * tauNonCO2 * np.exp(
        -t / tauT1) / ((tauT1 - tauNonCO2) * (tauT2 - tauNonCO2)) + aC1 * aT1 * tauT1 * aT2 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2)) + aC1 * aT1 * tauT1 * aT2 / (
                        tauT1 - tauNonCO2) - aC1 * aT1 * tauT1 * aT2 * np.exp(-t / tauT1) / (
                        tauT1 - tauNonCO2) + aC1 * aT1 * aT2 * tauT2 * tauNonCO2 * np.exp(-t / tauT2) / (
                        (tauT1 - tauNonCO2) * (tauT2 - tauNonCO2)) + aC1 * aT1 * aT2 * tauT2 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((-tauT1 + tauNonCO2) * (tauT2 - tauNonCO2)) + aC1 * aT1 * aT2 * tauT2 / (
                        tauT2 - tauNonCO2) - aC1 * aT1 * aT2 * tauT2 * np.exp(-t / tauT2) / (
                        tauT2 - tauNonCO2) - aC1 * aT1 * aT2 * tauNonCO2 / (
                        tauT2 - tauNonCO2) + aC1 * aT1 * aT2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        tauT2 - tauNonCO2) - aC1 * aT1 * aT2 * tauNonCO2 / (
                        tauT1 - tauNonCO2) + aC1 * aT1 * aT2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        tauT1 - tauNonCO2) - aC1 * aT2 ** 2 * tauT2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        tauT2 - tauNonCO2) ** 2 + aC1 * aT2 ** 2 * tauT2 * tauNonCO2 * np.exp(-t / tauT2) / (
                        tauT2 - tauNonCO2) ** 2 + aC1 * aT2 ** 2 * tauT2 / (
                        tauT2 - tauNonCO2) - aC1 * aT2 ** 2 * tauT2 * np.exp(-t / tauT2) / (
                        tauT2 - tauNonCO2) - aC1 * aT2 ** 2 * t * np.exp(-t / tauT2) / (
                        tauT2 - tauNonCO2) - aC1 * aT2 ** 2 * tauNonCO2 / (
                        tauT2 - tauNonCO2) + aC1 * aT2 ** 2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        tauT2 - tauNonCO2) + aC2 * aS1 * aT1 ** 2 * tauT1 ** 2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauC1) ** 2 * (tauT1 - tauNonCO2) * (
                tauT1 - tauS1)) + aC2 * aS1 * aT1 ** 2 * tauT1 ** 2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) ** 2 * (
                tauT1 - tauS1)) - aC2 * aS1 * aT1 ** 2 * tauT1 * t * tauC1 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS1)) - aC2 * aS1 * aT1 ** 2 * tauT1 * tauC1 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC1)) / ((tauT1 - tauC1) ** 2 * (tauT1 - tauNonCO2) * (
                tauC1 - tauS1)) - aC2 * aS1 * aT1 ** 2 * tauT1 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS1)) - aC2 * aS1 * aT1 ** 2 * tauT1 * tauC1 * (
                        -tauT1 * tauS1 + (tauT1 * (t + tauS1) - t * tauS1) * np.exp(
                    t * (-1 / tauS1 + 1 / tauT1))) * np.exp(-t / tauT1) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS1) ** 2) - aC2 * aS1 * aT1 ** 2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauC1 - tauS1)) + aC2 * aS1 * aT1 ** 2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) - aC2 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC1 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT1 - tauS1) * (
                            tauT2 - tauC1)) - aC2 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC1 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC1) * (tauT1 - tauS1) * (
                            tauT2 - tauNonCO2)) - aC2 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauT1 - tauS1) * (
                            -tauT2 + tauC1)) - aC2 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauT1 - tauS1) * (
                            -tauT2 + tauNonCO2)) + aC2 * aS1 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC1 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT2 - tauC1) * (
                            tauT2 - tauS1)) + aC2 * aS1 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC1 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS1)) + aC2 * aS1 * aT1 * tauT1 * aT2 * tauC1 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (-tauT2 + tauC1) * (
                            tauC1 - tauS1)) + aC2 * aS1 * aT1 * tauT1 * aT2 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                            tauNonCO2 - tauS1)) - aC2 * aS1 * aT1 * aT2 * tauT2 ** 2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauC1) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS1)) - aC2 * aS1 * aT1 * aT2 * tauT2 ** 2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS1)) + aC2 * aS1 * aT1 * aT2 * tauT2 * tauC1 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC1)) / (
                        (-tauT1 + tauC1) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauC1 - tauS1)) + aC2 * aS1 * aT1 * aT2 * tauT2 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) + aC2 * aS1 * aT1 * aT2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC1)) / (
                        (tauT1 - tauC1) * (-tauT2 + tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauC1 - tauS1)) + aC2 * aS1 * aT1 * aT2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC1)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauC1 - tauNonCO2) * (
                            tauC1 - tauS1)) - aC2 * aS1 * aT1 * aT2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC1) * (-tauT2 + tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) - aC2 * aS1 * aT1 * aT2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) + aC2 * aS1 * aT2 ** 2 * tauT2 ** 2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / ((tauT2 - tauC1) ** 2 * (tauT2 - tauNonCO2) * (
                tauT2 - tauS1)) + aC2 * aS1 * aT2 ** 2 * tauT2 ** 2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) ** 2 * (
                tauT2 - tauS1)) - aC2 * aS1 * aT2 ** 2 * tauT2 * t * tauC1 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS1)) - aC2 * aS1 * aT2 ** 2 * tauT2 * tauC1 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC1)) / ((tauT2 - tauC1) ** 2 * (tauT2 - tauNonCO2) * (
                tauC1 - tauS1)) - aC2 * aS1 * aT2 ** 2 * tauT2 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS1)) - aC2 * aS1 * aT2 ** 2 * tauT2 * tauC1 * (
                        -tauT2 * tauS1 + (tauT2 * (t + tauS1) - t * tauS1) * np.exp(
                    t * (-1 / tauS1 + 1 / tauT2))) * np.exp(-t / tauT2) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS1) ** 2) - aC2 * aS1 * aT2 ** 2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC1)) / (
                        (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauC1 - tauS1)) + aC2 * aS1 * aT2 ** 2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) + aC2 * aS2 * aT1 ** 2 * tauT1 ** 2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauC1) ** 2 * (tauT1 - tauNonCO2) * (
                tauT1 - tauS2)) + aC2 * aS2 * aT1 ** 2 * tauT1 ** 2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) ** 2 * (
                tauT1 - tauS2)) - aC2 * aS2 * aT1 ** 2 * tauT1 * t * tauC1 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS2)) - aC2 * aS2 * aT1 ** 2 * tauT1 * tauC1 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC1)) / ((tauT1 - tauC1) ** 2 * (tauT1 - tauNonCO2) * (
                tauC1 - tauS2)) - aC2 * aS2 * aT1 ** 2 * tauT1 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS2)) - aC2 * aS2 * aT1 ** 2 * tauT1 * tauC1 * (
                        -tauT1 * tauS2 + (tauT1 * (t + tauS2) - t * tauS2) * np.exp(
                    t * (-1 / tauS2 + 1 / tauT1))) * np.exp(-t / tauT1) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS2) ** 2) - aC2 * aS2 * aT1 ** 2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauC1 - tauS2)) + aC2 * aS2 * aT1 ** 2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) - aC2 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC1 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT1 - tauS2) * (
                            tauT2 - tauC1)) - aC2 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC1 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC1) * (tauT1 - tauS2) * (
                            tauT2 - tauNonCO2)) - aC2 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauT1 - tauS2) * (
                            -tauT2 + tauC1)) - aC2 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauT1 - tauS2) * (
                            -tauT2 + tauNonCO2)) + aC2 * aS2 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC1 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT2 - tauC1) * (
                            tauT2 - tauS2)) + aC2 * aS2 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC1 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS2)) + aC2 * aS2 * aT1 * tauT1 * aT2 * tauC1 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (-tauT2 + tauC1) * (
                            tauC1 - tauS2)) + aC2 * aS2 * aT1 * tauT1 * aT2 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                            tauNonCO2 - tauS2)) - aC2 * aS2 * aT1 * aT2 * tauT2 ** 2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauC1) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS2)) - aC2 * aS2 * aT1 * aT2 * tauT2 ** 2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS2)) + aC2 * aS2 * aT1 * aT2 * tauT2 * tauC1 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC1)) / (
                        (-tauT1 + tauC1) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauC1 - tauS2)) + aC2 * aS2 * aT1 * aT2 * tauT2 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) + aC2 * aS2 * aT1 * aT2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC1)) / (
                        (tauT1 - tauC1) * (-tauT2 + tauNonCO2) * (tauC1 - tauNonCO2) * (tauC1 - tauS2))

    term3 = aC2 * aS2 * aT1 * aT2 * tauC1 ** 3 * tauNonCO2 * (-np.exp(-t / tauS2) + np.exp(-t / tauC1)) / (
                (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauC1 - tauNonCO2) * (
                    tauC1 - tauS2)) - aC2 * aS2 * aT1 * aT2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC1) * (-tauT2 + tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) - aC2 * aS2 * aT1 * aT2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) + aC2 * aS2 * aT2 ** 2 * tauT2 ** 2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT2 - tauC1) ** 2 * (tauT2 - tauNonCO2) * (
                tauT2 - tauS2)) + aC2 * aS2 * aT2 ** 2 * tauT2 ** 2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) ** 2 * (
                tauT2 - tauS2)) - aC2 * aS2 * aT2 ** 2 * tauT2 * t * tauC1 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS2)) - aC2 * aS2 * aT2 ** 2 * tauT2 * tauC1 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC1)) / ((tauT2 - tauC1) ** 2 * (tauT2 - tauNonCO2) * (
                tauC1 - tauS2)) - aC2 * aS2 * aT2 ** 2 * tauT2 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS2)) - aC2 * aS2 * aT2 ** 2 * tauT2 * tauC1 * (
                        -tauT2 * tauS2 + (tauT2 * (t + tauS2) - t * tauS2) * np.exp(
                    t * (-1 / tauS2 + 1 / tauT2))) * np.exp(-t / tauT2) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS2) ** 2) - aC2 * aS2 * aT2 ** 2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC1)) / (
                        (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauC1 - tauS2)) + aC2 * aS2 * aT2 ** 2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) + aC2 * aS3 * aT1 ** 2 * tauT1 ** 2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauC1) ** 2 * (tauT1 - tauNonCO2) * (
                tauT1 - tauS3)) + aC2 * aS3 * aT1 ** 2 * tauT1 ** 2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) ** 2 * (
                tauT1 - tauS3)) - aC2 * aS3 * aT1 ** 2 * tauT1 * t * tauC1 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS3)) - aC2 * aS3 * aT1 ** 2 * tauT1 * tauC1 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC1)) / ((tauT1 - tauC1) ** 2 * (tauT1 - tauNonCO2) * (
                tauC1 - tauS3)) - aC2 * aS3 * aT1 ** 2 * tauT1 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS3)) - aC2 * aS3 * aT1 ** 2 * tauT1 * tauC1 * (
                        -tauT1 * tauS3 + (tauT1 * (t + tauS3) - t * tauS3) * np.exp(
                    t * (-1 / tauS3 + 1 / tauT1))) * np.exp(-t / tauT1) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS3) ** 2) - aC2 * aS3 * aT1 ** 2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauC1 - tauS3)) + aC2 * aS3 * aT1 ** 2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) - aC2 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC1 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT1 - tauS3) * (
                            tauT2 - tauC1)) - aC2 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC1 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC1) * (tauT1 - tauS3) * (
                            tauT2 - tauNonCO2)) - aC2 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauT1 - tauS3) * (
                            -tauT2 + tauC1)) - aC2 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (tauT1 - tauS3) * (
                            -tauT2 + tauNonCO2)) + aC2 * aS3 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC1 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT2 - tauC1) * (
                            tauT2 - tauS3)) + aC2 * aS3 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC1 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS3)) + aC2 * aS3 * aT1 * tauT1 * aT2 * tauC1 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC1)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (-tauT2 + tauC1) * (
                            tauC1 - tauS3)) + aC2 * aS3 * aT1 * tauT1 * aT2 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                            tauNonCO2 - tauS3)) - aC2 * aS3 * aT1 * aT2 * tauT2 ** 2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauC1) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS3)) - aC2 * aS3 * aT1 * aT2 * tauT2 ** 2 * tauC1 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS3)) + aC2 * aS3 * aT1 * aT2 * tauT2 * tauC1 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC1)) / (
                        (-tauT1 + tauC1) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauC1 - tauS3)) + aC2 * aS3 * aT1 * aT2 * tauT2 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) + aC2 * aS3 * aT1 * aT2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC1)) / (
                        (tauT1 - tauC1) * (-tauT2 + tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauC1 - tauS3)) + aC2 * aS3 * aT1 * aT2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC1)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauC1 - tauNonCO2) * (
                            tauC1 - tauS3)) - aC2 * aS3 * aT1 * aT2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC1) * (-tauT2 + tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) - aC2 * aS3 * aT1 * aT2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) + aC2 * aS3 * aT2 ** 2 * tauT2 ** 2 * tauC1 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauC1) ** 2 * (tauT2 - tauNonCO2) * (tauT2 - tauS3))

    term4 = aC2 * aS3 * aT2 ** 2 * tauT2 ** 2 * tauC1 * tauNonCO2 * (-np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                (tauT2 - tauC1) * (tauT2 - tauNonCO2) ** 2 * (
                    tauT2 - tauS3)) - aC2 * aS3 * aT2 ** 2 * tauT2 * t * tauC1 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS3)) - aC2 * aS3 * aT2 ** 2 * tauT2 * tauC1 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC1)) / ((tauT2 - tauC1) ** 2 * (tauT2 - tauNonCO2) * (
                tauC1 - tauS3)) - aC2 * aS3 * aT2 ** 2 * tauT2 * tauC1 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS3)) - aC2 * aS3 * aT2 ** 2 * tauT2 * tauC1 * (
                        -tauT2 * tauS3 + (tauT2 * (t + tauS3) - t * tauS3) * np.exp(
                    t * (-1 / tauS3 + 1 / tauT2))) * np.exp(-t / tauT2) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS3) ** 2) - aC2 * aS3 * aT2 ** 2 * tauC1 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC1)) / (
                        (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauC1 - tauS3)) + aC2 * aS3 * aT2 ** 2 * tauC1 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (tauC1 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) + aC2 * aT1 ** 2 * tauT1 * tauC1 ** 2 * np.exp(-t / tauC1) / (
                        (tauT1 - tauC1) ** 2 * (tauT1 - tauNonCO2)) - aC2 * aT1 ** 2 * tauT1 * tauC1 ** 2 * np.exp(
        -t / tauT1) / ((tauT1 - tauC1) ** 2 * (
                tauT1 - tauNonCO2)) + aC2 * aT1 ** 2 * tauT1 * tauC1 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (tauT1 - tauC1) * (
                            tauT1 - tauNonCO2) ** 2) - aC2 * aT1 ** 2 * tauT1 * tauC1 * tauNonCO2 * np.exp(
        -t / tauT1) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) ** 2) + aC2 * aT1 ** 2 * t * tauC1 * np.exp(-t / tauT1) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2)) - aC2 * aT1 ** 2 * tauC1 ** 2 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                tauC1 - tauNonCO2)) + aC2 * aT1 ** 2 * tauC1 ** 2 * tauNonCO2 * np.exp(-t / tauC1) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                            tauC1 - tauNonCO2)) - aC2 * aT1 * tauT1 * aT2 * tauT2 * tauC1 * np.exp(-t / tauT2) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauC1)) + aC2 * aT1 * tauT1 * aT2 * tauT2 * tauC1 * np.exp(-t / tauT1) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauC1)) - aC2 * aT1 * tauT1 * aT2 * tauT2 * tauC1 * np.exp(-t / tauT2) / (
                        (tauT1 - tauT2) * (tauT1 - tauC1) * (
                            tauT2 - tauNonCO2)) + aC2 * aT1 * tauT1 * aT2 * tauT2 * tauC1 * np.exp(-t / tauT1) / (
                        (tauT1 - tauT2) * (tauT1 - tauC1) * (
                            tauT2 - tauNonCO2)) + aC2 * aT1 * tauT1 * aT2 * tauC1 ** 2 * np.exp(-t / tauC1) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauC1)) + aC2 * aT1 * tauT1 * aT2 * tauC1 ** 2 * np.exp(-t / tauT1) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                            -tauT2 + tauC1)) + aC2 * aT1 * tauT1 * aT2 * tauC1 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauNonCO2)) + aC2 * aT1 * tauT1 * aT2 * tauC1 * tauNonCO2 * np.exp(-t / tauT1) / (
                        (tauT1 - tauC1) * (tauT1 - tauNonCO2) * (
                            -tauT2 + tauNonCO2)) + aC2 * aT1 * aT2 * tauT2 * tauC1 ** 2 * np.exp(-t / tauC1) / (
                        (tauT1 - tauC1) * (tauT2 - tauC1) * (
                            tauT2 - tauNonCO2)) + aC2 * aT1 * aT2 * tauT2 * tauC1 ** 2 * np.exp(-t / tauT2) / (
                        (-tauT1 + tauC1) * (tauT2 - tauC1) * (
                            tauT2 - tauNonCO2)) + aC2 * aT1 * aT2 * tauT2 * tauC1 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((tauT1 - tauNonCO2) * (tauT2 - tauC1) * (
                tauT2 - tauNonCO2)) + aC2 * aT1 * aT2 * tauT2 * tauC1 * tauNonCO2 * np.exp(-t / tauT2) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC1) * (
                            tauT2 - tauNonCO2)) + aC2 * aT1 * aT2 * tauC1 ** 2 * tauNonCO2 * np.exp(-t / tauC1) / (
                        (-tauT1 + tauNonCO2) * (-tauT2 + tauC1) * (
                            tauC1 - tauNonCO2)) + aC2 * aT1 * aT2 * tauC1 ** 2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (-tauT1 + tauNonCO2) * (-tauT2 + tauC1) * (
                            -tauC1 + tauNonCO2)) + aC2 * aT1 * aT2 * tauC1 ** 2 * tauNonCO2 * np.exp(-t / tauC1) / (
                        (-tauT1 + tauC1) * (-tauT2 + tauNonCO2) * (
                            tauC1 - tauNonCO2)) + aC2 * aT1 * aT2 * tauC1 ** 2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (-tauT1 + tauC1) * (-tauT2 + tauNonCO2) * (
                            -tauC1 + tauNonCO2)) + aC2 * aT2 ** 2 * tauT2 * tauC1 ** 2 * np.exp(-t / tauC1) / (
                        (tauT2 - tauC1) ** 2 * (tauT2 - tauNonCO2)) - aC2 * aT2 ** 2 * tauT2 * tauC1 ** 2 * np.exp(
        -t / tauT2) / ((tauT2 - tauC1) ** 2 * (
                tauT2 - tauNonCO2)) + aC2 * aT2 ** 2 * tauT2 * tauC1 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (tauT2 - tauC1) * (
                            tauT2 - tauNonCO2) ** 2) - aC2 * aT2 ** 2 * tauT2 * tauC1 * tauNonCO2 * np.exp(
        -t / tauT2) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) ** 2) + aC2 * aT2 ** 2 * t * tauC1 * np.exp(-t / tauT2) / (
                        (tauT2 - tauC1) * (tauT2 - tauNonCO2)) - aC2 * aT2 ** 2 * tauC1 ** 2 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                tauC1 - tauNonCO2)) + aC2 * aT2 ** 2 * tauC1 ** 2 * tauNonCO2 * np.exp(-t / tauC1) / (
                        (tauT2 - tauC1) * (tauT2 - tauNonCO2) * (
                            tauC1 - tauNonCO2)) + aC3 * aS1 * aT1 ** 2 * tauT1 ** 2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauC2) ** 2 * (tauT1 - tauNonCO2) * (
                tauT1 - tauS1)) + aC3 * aS1 * aT1 ** 2 * tauT1 ** 2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) ** 2 * (
                tauT1 - tauS1)) - aC3 * aS1 * aT1 ** 2 * tauT1 * t * tauC2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS1)) - aC3 * aS1 * aT1 ** 2 * tauT1 * tauC2 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC2)) / ((tauT1 - tauC2) ** 2 * (tauT1 - tauNonCO2) * (
                tauC2 - tauS1)) - aC3 * aS1 * aT1 ** 2 * tauT1 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS1)) - aC3 * aS1 * aT1 ** 2 * tauT1 * tauC2 * (
                        -tauT1 * tauS1 + (tauT1 * (t + tauS1) - t * tauS1) * np.exp(
                    t * (-1 / tauS1 + 1 / tauT1))) * np.exp(-t / tauT1) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS1) ** 2) - aC3 * aS1 * aT1 ** 2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS1)) + aC3 * aS1 * aT1 ** 2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) - aC3 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT1 - tauS1) * (
                            tauT2 - tauC2)) - aC3 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC2) * (tauT1 - tauS1) * (
                            tauT2 - tauNonCO2)) - aC3 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauT1 - tauS1) * (
                            -tauT2 + tauC2)) - aC3 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauT1 - tauS1) * (
                            -tauT2 + tauNonCO2)) + aC3 * aS1 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT2 - tauC2) * (
                            tauT2 - tauS1)) + aC3 * aS1 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS1)) + aC3 * aS1 * aT1 * tauT1 * aT2 * tauC2 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (-tauT2 + tauC2) * (
                            tauC2 - tauS1)) + aC3 * aS1 * aT1 * tauT1 * aT2 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                            tauNonCO2 - tauS1)) - aC3 * aS1 * aT1 * aT2 * tauT2 ** 2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauC2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS1)) - aC3 * aS1 * aT1 * aT2 * tauT2 ** 2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS1)) + aC3 * aS1 * aT1 * aT2 * tauT2 * tauC2 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC2)) / (
                        (-tauT1 + tauC2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauC2 - tauS1)) + aC3 * aS1 * aT1 * aT2 * tauT2 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) + aC3 * aS1 * aT1 * aT2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC2)) / (
                        (tauT1 - tauC2) * (-tauT2 + tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS1)) + aC3 * aS1 * aT1 * aT2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS1)) - aC3 * aS1 * aT1 * aT2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC2) * (-tauT2 + tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) - aC3 * aS1 * aT1 * aT2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) + aC3 * aS1 * aT2 ** 2 * tauT2 ** 2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT2 - tauC2) ** 2 * (tauT2 - tauNonCO2) * (tauT2 - tauS1))

    term5 = aC3 * aS1 * aT2 ** 2 * tauT2 ** 2 * tauC2 * tauNonCO2 * (-np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                (tauT2 - tauC2) * (tauT2 - tauNonCO2) ** 2 * (
                    tauT2 - tauS1)) - aC3 * aS1 * aT2 ** 2 * tauT2 * t * tauC2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS1)) - aC3 * aS1 * aT2 ** 2 * tauT2 * tauC2 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC2)) / ((tauT2 - tauC2) ** 2 * (tauT2 - tauNonCO2) * (
                tauC2 - tauS1)) - aC3 * aS1 * aT2 ** 2 * tauT2 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS1)) - aC3 * aS1 * aT2 ** 2 * tauT2 * tauC2 * (
                        -tauT2 * tauS1 + (tauT2 * (t + tauS1) - t * tauS1) * np.exp(
                    t * (-1 / tauS1 + 1 / tauT2))) * np.exp(-t / tauT2) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS1) ** 2) - aC3 * aS1 * aT2 ** 2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC2)) / (
                        (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS1)) + aC3 * aS1 * aT2 ** 2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) + aC3 * aS2 * aT1 ** 2 * tauT1 ** 2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauC2) ** 2 * (tauT1 - tauNonCO2) * (
                tauT1 - tauS2)) + aC3 * aS2 * aT1 ** 2 * tauT1 ** 2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) ** 2 * (
                tauT1 - tauS2)) - aC3 * aS2 * aT1 ** 2 * tauT1 * t * tauC2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS2)) - aC3 * aS2 * aT1 ** 2 * tauT1 * tauC2 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC2)) / ((tauT1 - tauC2) ** 2 * (tauT1 - tauNonCO2) * (
                tauC2 - tauS2)) - aC3 * aS2 * aT1 ** 2 * tauT1 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS2)) - aC3 * aS2 * aT1 ** 2 * tauT1 * tauC2 * (
                        -tauT1 * tauS2 + (tauT1 * (t + tauS2) - t * tauS2) * np.exp(
                    t * (-1 / tauS2 + 1 / tauT1))) * np.exp(-t / tauT1) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS2) ** 2) - aC3 * aS2 * aT1 ** 2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS2)) + aC3 * aS2 * aT1 ** 2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) - aC3 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT1 - tauS2) * (
                            tauT2 - tauC2)) - aC3 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC2) * (tauT1 - tauS2) * (
                            tauT2 - tauNonCO2)) - aC3 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauT1 - tauS2) * (
                            -tauT2 + tauC2)) - aC3 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauT1 - tauS2) * (
                            -tauT2 + tauNonCO2)) + aC3 * aS2 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT2 - tauC2) * (
                            tauT2 - tauS2)) + aC3 * aS2 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS2)) + aC3 * aS2 * aT1 * tauT1 * aT2 * tauC2 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (-tauT2 + tauC2) * (
                            tauC2 - tauS2)) + aC3 * aS2 * aT1 * tauT1 * aT2 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                            tauNonCO2 - tauS2)) - aC3 * aS2 * aT1 * aT2 * tauT2 ** 2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauC2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS2)) - aC3 * aS2 * aT1 * aT2 * tauT2 ** 2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS2)) + aC3 * aS2 * aT1 * aT2 * tauT2 * tauC2 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC2)) / (
                        (-tauT1 + tauC2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauC2 - tauS2)) + aC3 * aS2 * aT1 * aT2 * tauT2 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) + aC3 * aS2 * aT1 * aT2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC2)) / (
                        (tauT1 - tauC2) * (-tauT2 + tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS2)) + aC3 * aS2 * aT1 * aT2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS2)) - aC3 * aS2 * aT1 * aT2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC2) * (-tauT2 + tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) - aC3 * aS2 * aT1 * aT2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) + aC3 * aS2 * aT2 ** 2 * tauT2 ** 2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT2 - tauC2) ** 2 * (tauT2 - tauNonCO2) * (
                tauT2 - tauS2)) + aC3 * aS2 * aT2 ** 2 * tauT2 ** 2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) ** 2 * (
                tauT2 - tauS2)) - aC3 * aS2 * aT2 ** 2 * tauT2 * t * tauC2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS2)) - aC3 * aS2 * aT2 ** 2 * tauT2 * tauC2 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC2)) / ((tauT2 - tauC2) ** 2 * (tauT2 - tauNonCO2) * (
                tauC2 - tauS2)) - aC3 * aS2 * aT2 ** 2 * tauT2 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS2)) - aC3 * aS2 * aT2 ** 2 * tauT2 * tauC2 * (
                        -tauT2 * tauS2 + (tauT2 * (t + tauS2) - t * tauS2) * np.exp(
                    t * (-1 / tauS2 + 1 / tauT2))) * np.exp(-t / tauT2) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS2) ** 2) - aC3 * aS2 * aT2 ** 2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC2)) / (
                        (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS2)) + aC3 * aS2 * aT2 ** 2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) + aC3 * aS3 * aT1 ** 2 * tauT1 ** 2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauC2) ** 2 * (tauT1 - tauNonCO2) * (
                tauT1 - tauS3)) + aC3 * aS3 * aT1 ** 2 * tauT1 ** 2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) ** 2 * (
                tauT1 - tauS3)) - aC3 * aS3 * aT1 ** 2 * tauT1 * t * tauC2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS3)) - aC3 * aS3 * aT1 ** 2 * tauT1 * tauC2 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC2)) / ((tauT1 - tauC2) ** 2 * (tauT1 - tauNonCO2) * (
                tauC2 - tauS3)) - aC3 * aS3 * aT1 ** 2 * tauT1 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS3)) - aC3 * aS3 * aT1 ** 2 * tauT1 * tauC2 * (
                        -tauT1 * tauS3 + (tauT1 * (t + tauS3) - t * tauS3) * np.exp(
                    t * (-1 / tauS3 + 1 / tauT1))) * np.exp(-t / tauT1) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS3) ** 2) - aC3 * aS3 * aT1 ** 2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS3)) + aC3 * aS3 * aT1 ** 2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) - aC3 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT1 - tauS3) * (
                            tauT2 - tauC2)) - aC3 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC2) * (tauT1 - tauS3) * (
                            tauT2 - tauNonCO2)) - aC3 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauT1 - tauS3) * (
                            -tauT2 + tauC2)) - aC3 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (tauT1 - tauS3) * (
                            -tauT2 + tauNonCO2)) + aC3 * aS3 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT2 - tauC2) * (
                            tauT2 - tauS3)) + aC3 * aS3 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS3)) + aC3 * aS3 * aT1 * tauT1 * aT2 * tauC2 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (-tauT2 + tauC2) * (
                            tauC2 - tauS3)) + aC3 * aS3 * aT1 * tauT1 * aT2 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                            tauNonCO2 - tauS3)) - aC3 * aS3 * aT1 * aT2 * tauT2 ** 2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauC2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS3)) - aC3 * aS3 * aT1 * aT2 * tauT2 ** 2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS3)) + aC3 * aS3 * aT1 * aT2 * tauT2 * tauC2 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC2)) / (
                        (-tauT1 + tauC2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauC2 - tauS3)) + aC3 * aS3 * aT1 * aT2 * tauT2 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) + aC3 * aS3 * aT1 * aT2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC2)) / (
                        (tauT1 - tauC2) * (-tauT2 + tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS3)) + aC3 * aS3 * aT1 * aT2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS3)) - aC3 * aS3 * aT1 * aT2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC2) * (-tauT2 + tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) - aC3 * aS3 * aT1 * aT2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) + aC3 * aS3 * aT2 ** 2 * tauT2 ** 2 * tauC2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / ((tauT2 - tauC2) ** 2 * (tauT2 - tauNonCO2) * (
                tauT2 - tauS3)) + aC3 * aS3 * aT2 ** 2 * tauT2 ** 2 * tauC2 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) ** 2 * (
                tauT2 - tauS3)) - aC3 * aS3 * aT2 ** 2 * tauT2 * t * tauC2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS3)) - aC3 * aS3 * aT2 ** 2 * tauT2 * tauC2 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC2)) / ((tauT2 - tauC2) ** 2 * (tauT2 - tauNonCO2) * (
                tauC2 - tauS3)) - aC3 * aS3 * aT2 ** 2 * tauT2 * tauC2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS3)) - aC3 * aS3 * aT2 ** 2 * tauT2 * tauC2 * (
                        -tauT2 * tauS3 + (tauT2 * (t + tauS3) - t * tauS3) * np.exp(
                    t * (-1 / tauS3 + 1 / tauT2))) * np.exp(-t / tauT2) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS3) ** 2) - aC3 * aS3 * aT2 ** 2 * tauC2 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC2)) / (
                        (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauC2 - tauS3)) + aC3 * aS3 * aT2 ** 2 * tauC2 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (tauC2 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) + aC3 * aT1 ** 2 * tauT1 * tauC2 ** 2 * np.exp(-t / tauC2) / (
                        (tauT1 - tauC2) ** 2 * (tauT1 - tauNonCO2)) - aC3 * aT1 ** 2 * tauT1 * tauC2 ** 2 * np.exp(
        -t / tauT1) / ((tauT1 - tauC2) ** 2 * (
                tauT1 - tauNonCO2)) + aC3 * aT1 ** 2 * tauT1 * tauC2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (tauT1 - tauC2) * (
                            tauT1 - tauNonCO2) ** 2) - aC3 * aT1 ** 2 * tauT1 * tauC2 * tauNonCO2 * np.exp(
        -t / tauT1) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) ** 2) + aC3 * aT1 ** 2 * t * tauC2 * np.exp(-t / tauT1) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2)) - aC3 * aT1 ** 2 * tauC2 ** 2 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                tauC2 - tauNonCO2)) + aC3 * aT1 ** 2 * tauC2 ** 2 * tauNonCO2 * np.exp(-t / tauC2) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                            tauC2 - tauNonCO2)) - aC3 * aT1 * tauT1 * aT2 * tauT2 * tauC2 * np.exp(-t / tauT2) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauC2)) + aC3 * aT1 * tauT1 * aT2 * tauT2 * tauC2 * np.exp(-t / tauT1) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauC2)) - aC3 * aT1 * tauT1 * aT2 * tauT2 * tauC2 * np.exp(-t / tauT2) / (
                        (tauT1 - tauT2) * (tauT1 - tauC2) * (
                            tauT2 - tauNonCO2)) + aC3 * aT1 * tauT1 * aT2 * tauT2 * tauC2 * np.exp(-t / tauT1) / (
                        (tauT1 - tauT2) * (tauT1 - tauC2) * (
                            tauT2 - tauNonCO2)) + aC3 * aT1 * tauT1 * aT2 * tauC2 ** 2 * np.exp(-t / tauC2) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauC2)) + aC3 * aT1 * tauT1 * aT2 * tauC2 ** 2 * np.exp(-t / tauT1) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                            -tauT2 + tauC2)) + aC3 * aT1 * tauT1 * aT2 * tauC2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauNonCO2)) + aC3 * aT1 * tauT1 * aT2 * tauC2 * tauNonCO2 * np.exp(-t / tauT1) / (
                        (tauT1 - tauC2) * (tauT1 - tauNonCO2) * (
                            -tauT2 + tauNonCO2)) + aC3 * aT1 * aT2 * tauT2 * tauC2 ** 2 * np.exp(-t / tauC2) / (
                        (tauT1 - tauC2) * (tauT2 - tauC2) * (
                            tauT2 - tauNonCO2)) + aC3 * aT1 * aT2 * tauT2 * tauC2 ** 2 * np.exp(-t / tauT2) / (
                        (-tauT1 + tauC2) * (tauT2 - tauC2) * (
                            tauT2 - tauNonCO2)) + aC3 * aT1 * aT2 * tauT2 * tauC2 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((tauT1 - tauNonCO2) * (tauT2 - tauC2) * (
                tauT2 - tauNonCO2)) + aC3 * aT1 * aT2 * tauT2 * tauC2 * tauNonCO2 * np.exp(-t / tauT2) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC2) * (
                            tauT2 - tauNonCO2)) + aC3 * aT1 * aT2 * tauC2 ** 2 * tauNonCO2 * np.exp(-t / tauC2) / (
                        (-tauT1 + tauNonCO2) * (-tauT2 + tauC2) * (
                            tauC2 - tauNonCO2)) + aC3 * aT1 * aT2 * tauC2 ** 2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (-tauT1 + tauNonCO2) * (-tauT2 + tauC2) * (
                            -tauC2 + tauNonCO2)) + aC3 * aT1 * aT2 * tauC2 ** 2 * tauNonCO2 * np.exp(-t / tauC2) / (
                        (-tauT1 + tauC2) * (-tauT2 + tauNonCO2) * (
                            tauC2 - tauNonCO2)) + aC3 * aT1 * aT2 * tauC2 ** 2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (-tauT1 + tauC2) * (-tauT2 + tauNonCO2) * (
                            -tauC2 + tauNonCO2)) + aC3 * aT2 ** 2 * tauT2 * tauC2 ** 2 * np.exp(-t / tauC2) / (
                        (tauT2 - tauC2) ** 2 * (tauT2 - tauNonCO2)) - aC3 * aT2 ** 2 * tauT2 * tauC2 ** 2 * np.exp(
        -t / tauT2) / ((tauT2 - tauC2) ** 2 * (
                tauT2 - tauNonCO2)) + aC3 * aT2 ** 2 * tauT2 * tauC2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (tauT2 - tauC2) * (
                            tauT2 - tauNonCO2) ** 2) - aC3 * aT2 ** 2 * tauT2 * tauC2 * tauNonCO2 * np.exp(
        -t / tauT2) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) ** 2) + aC3 * aT2 ** 2 * t * tauC2 * np.exp(-t / tauT2) / (
                        (tauT2 - tauC2) * (tauT2 - tauNonCO2)) - aC3 * aT2 ** 2 * tauC2 ** 2 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                tauC2 - tauNonCO2)) + aC3 * aT2 ** 2 * tauC2 ** 2 * tauNonCO2 * np.exp(-t / tauC2) / (
                        (tauT2 - tauC2) * (tauT2 - tauNonCO2) * (
                            tauC2 - tauNonCO2)) + aC4 * aS1 * aT1 ** 2 * tauT1 ** 2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauC3) ** 2 * (tauT1 - tauNonCO2) * (
                tauT1 - tauS1)) + aC4 * aS1 * aT1 ** 2 * tauT1 ** 2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) ** 2 * (
                tauT1 - tauS1)) - aC4 * aS1 * aT1 ** 2 * tauT1 * t * tauC3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS1)) - aC4 * aS1 * aT1 ** 2 * tauT1 * tauC3 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC3)) / ((tauT1 - tauC3) ** 2 * (tauT1 - tauNonCO2) * (
                tauC3 - tauS1)) - aC4 * aS1 * aT1 ** 2 * tauT1 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS1)) - aC4 * aS1 * aT1 ** 2 * tauT1 * tauC3 * (
                        -tauT1 * tauS1 + (tauT1 * (t + tauS1) - t * tauS1) * np.exp(
                    t * (-1 / tauS1 + 1 / tauT1))) * np.exp(-t / tauT1) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS1) ** 2) - aC4 * aS1 * aT1 ** 2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC3)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS1)) + aC4 * aS1 * aT1 ** 2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) - aC4 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT1 - tauS1) * (
                            tauT2 - tauC3)) - aC4 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC3) * (tauT1 - tauS1) * (
                            tauT2 - tauNonCO2)) - aC4 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauT1 - tauS1) * (
                            -tauT2 + tauC3)) - aC4 * aS1 * aT1 * tauT1 ** 2 * aT2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauT1 - tauS1) * (
                            -tauT2 + tauNonCO2)) + aC4 * aS1 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT2 - tauC3) * (
                            tauT2 - tauS1)) + aC4 * aS1 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS1)) + aC4 * aS1 * aT1 * tauT1 * aT2 * tauC3 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC3)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (-tauT2 + tauC3) * (
                            tauC3 - tauS1)) + aC4 * aS1 * aT1 * tauT1 * aT2 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                            tauNonCO2 - tauS1)) - aC4 * aS1 * aT1 * aT2 * tauT2 ** 2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauC3) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS1)) - aC4 * aS1 * aT1 * aT2 * tauT2 ** 2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS1)) + aC4 * aS1 * aT1 * aT2 * tauT2 * tauC3 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC3)) / (
                        (-tauT1 + tauC3) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauC3 - tauS1)) + aC4 * aS1 * aT1 * aT2 * tauT2 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) + aC4 * aS1 * aT1 * aT2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC3)) / (
                        (tauT1 - tauC3) * (-tauT2 + tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS1)) + aC4 * aS1 * aT1 * aT2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC3)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS1)) - aC4 * aS1 * aT1 * aT2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC3) * (-tauT2 + tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) - aC4 * aS1 * aT1 * aT2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) + aC4 * aS1 * aT2 ** 2 * tauT2 ** 2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / ((tauT2 - tauC3) ** 2 * (tauT2 - tauNonCO2) * (
                tauT2 - tauS1)) + aC4 * aS1 * aT2 ** 2 * tauT2 ** 2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) ** 2 * (
                tauT2 - tauS1)) - aC4 * aS1 * aT2 ** 2 * tauT2 * t * tauC3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauT2)) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS1)) - aC4 * aS1 * aT2 ** 2 * tauT2 * tauC3 ** 3 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC3)) / ((tauT2 - tauC3) ** 2 * (tauT2 - tauNonCO2) * (
                tauC3 - tauS1)) - aC4 * aS1 * aT2 ** 2 * tauT2 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS1)) - aC4 * aS1 * aT2 ** 2 * tauT2 * tauC3 * (
                        -tauT2 * tauS1 + (tauT2 * (t + tauS1) - t * tauS1) * np.exp(
                    t * (-1 / tauS1 + 1 / tauT2))) * np.exp(-t / tauT2) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS1) ** 2) - aC4 * aS1 * aT2 ** 2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauC3)) / (
                        (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS1)) + aC4 * aS1 * aT2 ** 2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS1) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS1)) + aC4 * aS2 * aT1 ** 2 * tauT1 ** 2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauC3) ** 2 * (tauT1 - tauNonCO2) * (
                tauT1 - tauS2)) + aC4 * aS2 * aT1 ** 2 * tauT1 ** 2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) ** 2 * (
                tauT1 - tauS2)) - aC4 * aS2 * aT1 ** 2 * tauT1 * t * tauC3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS2)) - aC4 * aS2 * aT1 ** 2 * tauT1 * tauC3 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC3)) / ((tauT1 - tauC3) ** 2 * (tauT1 - tauNonCO2) * (
                tauC3 - tauS2)) - aC4 * aS2 * aT1 ** 2 * tauT1 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS2)) - aC4 * aS2 * aT1 ** 2 * tauT1 * tauC3 * (
                        -tauT1 * tauS2 + (tauT1 * (t + tauS2) - t * tauS2) * np.exp(
                    t * (-1 / tauS2 + 1 / tauT1))) * np.exp(-t / tauT1) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS2) ** 2) - aC4 * aS2 * aT1 ** 2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC3)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS2)) + aC4 * aS2 * aT1 ** 2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) - aC4 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT1 - tauS2) * (
                            tauT2 - tauC3)) - aC4 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC3) * (tauT1 - tauS2) * (
                            tauT2 - tauNonCO2)) - aC4 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauT1 - tauS2) * (
                            -tauT2 + tauC3)) - aC4 * aS2 * aT1 * tauT1 ** 2 * aT2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauT1 - tauS2) * (
                            -tauT2 + tauNonCO2)) + aC4 * aS2 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT2 - tauC3) * (
                            tauT2 - tauS2)) + aC4 * aS2 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS2)) + aC4 * aS2 * aT1 * tauT1 * aT2 * tauC3 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC3)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (-tauT2 + tauC3) * (
                            tauC3 - tauS2)) + aC4 * aS2 * aT1 * tauT1 * aT2 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                            tauNonCO2 - tauS2)) - aC4 * aS2 * aT1 * aT2 * tauT2 ** 2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauC3) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS2)) - aC4 * aS2 * aT1 * aT2 * tauT2 ** 2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS2)) + aC4 * aS2 * aT1 * aT2 * tauT2 * tauC3 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC3)) / (
                        (-tauT1 + tauC3) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauC3 - tauS2)) + aC4 * aS2 * aT1 * aT2 * tauT2 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) + aC4 * aS2 * aT1 * aT2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC3)) / (
                        (tauT1 - tauC3) * (-tauT2 + tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS2)) + aC4 * aS2 * aT1 * aT2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC3)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS2)) - aC4 * aS2 * aT1 * aT2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC3) * (-tauT2 + tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) - aC4 * aS2 * aT1 * aT2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) + aC4 * aS2 * aT2 ** 2 * tauT2 ** 2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT2 - tauC3) ** 2 * (tauT2 - tauNonCO2) * (
                tauT2 - tauS2)) + aC4 * aS2 * aT2 ** 2 * tauT2 ** 2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) ** 2 * (
                tauT2 - tauS2)) - aC4 * aS2 * aT2 ** 2 * tauT2 * t * tauC3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauT2)) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS2)) - aC4 * aS2 * aT2 ** 2 * tauT2 * tauC3 ** 3 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC3)) / ((tauT2 - tauC3) ** 2 * (tauT2 - tauNonCO2) * (
                tauC3 - tauS2)) - aC4 * aS2 * aT2 ** 2 * tauT2 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS2)) - aC4 * aS2 * aT2 ** 2 * tauT2 * tauC3 * (
                        -tauT2 * tauS2 + (tauT2 * (t + tauS2) - t * tauS2) * np.exp(
                    t * (-1 / tauS2 + 1 / tauT2))) * np.exp(-t / tauT2) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS2) ** 2) - aC4 * aS2 * aT2 ** 2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauC3)) / (
                        (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS2)) + aC4 * aS2 * aT2 ** 2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS2) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS2)) + aC4 * aS3 * aT1 ** 2 * tauT1 ** 2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauC3) ** 2 * (tauT1 - tauNonCO2) * (
                tauT1 - tauS3)) + aC4 * aS3 * aT1 ** 2 * tauT1 ** 2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) ** 2 * (
                tauT1 - tauS3)) - aC4 * aS3 * aT1 ** 2 * tauT1 * t * tauC3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS3)) - aC4 * aS3 * aT1 ** 2 * tauT1 * tauC3 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC3)) / ((tauT1 - tauC3) ** 2 * (tauT1 - tauNonCO2) * (
                tauC3 - tauS3)) - aC4 * aS3 * aT1 ** 2 * tauT1 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS3)) - aC4 * aS3 * aT1 ** 2 * tauT1 * tauC3 * (
                        -tauT1 * tauS3 + (tauT1 * (t + tauS3) - t * tauS3) * np.exp(
                    t * (-1 / tauS3 + 1 / tauT1))) * np.exp(-t / tauT1) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                tauT1 - tauS3) ** 2) - aC4 * aS3 * aT1 ** 2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC3)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS3)) + aC4 * aS3 * aT1 ** 2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) - aC4 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT1 - tauS3) * (
                            tauT2 - tauC3)) - aC4 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauT2 * tauC3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC3) * (tauT1 - tauS3) * (
                            tauT2 - tauNonCO2)) - aC4 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauT1 - tauS3) * (
                            -tauT2 + tauC3)) - aC4 * aS3 * aT1 * tauT1 ** 2 * aT2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT1)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (tauT1 - tauS3) * (
                            -tauT2 + tauNonCO2)) + aC4 * aS3 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (tauT2 - tauC3) * (
                            tauT2 - tauS3)) + aC4 * aS3 * aT1 * tauT1 * aT2 * tauT2 ** 2 * tauC3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (tauT1 - tauT2) * (tauT1 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS3)) + aC4 * aS3 * aT1 * tauT1 * aT2 * tauC3 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC3)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (-tauT2 + tauC3) * (
                            tauC3 - tauS3)) + aC4 * aS3 * aT1 * tauT1 * aT2 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (-tauT2 + tauNonCO2) * (
                            tauNonCO2 - tauS3)) - aC4 * aS3 * aT1 * aT2 * tauT2 ** 2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauC3) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS3)) - aC4 * aS3 * aT1 * aT2 * tauT2 ** 2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauT2 - tauS3)) + aC4 * aS3 * aT1 * aT2 * tauT2 * tauC3 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC3)) / (
                        (-tauT1 + tauC3) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauC3 - tauS3)) + aC4 * aS3 * aT1 * aT2 * tauT2 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) + aC4 * aS3 * aT1 * aT2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC3)) / (
                        (tauT1 - tauC3) * (-tauT2 + tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS3)) + aC4 * aS3 * aT1 * aT2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC3)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS3)) - aC4 * aS3 * aT1 * aT2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT1 - tauC3) * (-tauT2 + tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) - aC4 * aS3 * aT1 * aT2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) + aC4 * aS3 * aT2 ** 2 * tauT2 ** 2 * tauC3 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / ((tauT2 - tauC3) ** 2 * (tauT2 - tauNonCO2) * (
                tauT2 - tauS3)) + aC4 * aS3 * aT2 ** 2 * tauT2 ** 2 * tauC3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) ** 2 * (
                tauT2 - tauS3)) - aC4 * aS3 * aT2 ** 2 * tauT2 * t * tauC3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauT2)) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS3)) - aC4 * aS3 * aT2 ** 2 * tauT2 * tauC3 ** 3 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC3)) / ((tauT2 - tauC3) ** 2 * (tauT2 - tauNonCO2) * (
                tauC3 - tauS3)) - aC4 * aS3 * aT2 ** 2 * tauT2 * tauC3 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) ** 2 * (
                tauNonCO2 - tauS3)) - aC4 * aS3 * aT2 ** 2 * tauT2 * tauC3 * (
                        -tauT2 * tauS3 + (tauT2 * (t + tauS3) - t * tauS3) * np.exp(
                    t * (-1 / tauS3 + 1 / tauT2))) * np.exp(-t / tauT2) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                tauT2 - tauS3) ** 2) - aC4 * aS3 * aT2 ** 2 * tauC3 ** 3 * tauNonCO2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauC3)) / (
                        (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauC3 - tauS3)) + aC4 * aS3 * aT2 ** 2 * tauC3 ** 2 * tauNonCO2 ** 2 * (
                        -np.exp(-t / tauS3) + np.exp(-t / tauNonCO2)) / (
                        (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (tauC3 - tauNonCO2) * (
                            tauNonCO2 - tauS3)) + aC4 * aT1 ** 2 * tauT1 * tauC3 ** 2 * np.exp(-t / tauC3) / (
                        (tauT1 - tauC3) ** 2 * (tauT1 - tauNonCO2)) - aC4 * aT1 ** 2 * tauT1 * tauC3 ** 2 * np.exp(
        -t / tauT1) / ((tauT1 - tauC3) ** 2 * (
                tauT1 - tauNonCO2)) + aC4 * aT1 ** 2 * tauT1 * tauC3 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (tauT1 - tauC3) * (
                            tauT1 - tauNonCO2) ** 2) - aC4 * aT1 ** 2 * tauT1 * tauC3 * tauNonCO2 * np.exp(
        -t / tauT1) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) ** 2) + aC4 * aT1 ** 2 * t * tauC3 * np.exp(-t / tauT1) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2)) - aC4 * aT1 ** 2 * tauC3 ** 2 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                tauC3 - tauNonCO2)) + aC4 * aT1 ** 2 * tauC3 ** 2 * tauNonCO2 * np.exp(-t / tauC3) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                            tauC3 - tauNonCO2)) - aC4 * aT1 * tauT1 * aT2 * tauT2 * tauC3 * np.exp(-t / tauT2) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauC3)) + aC4 * aT1 * tauT1 * aT2 * tauT2 * tauC3 * np.exp(-t / tauT1) / (
                        (tauT1 - tauT2) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauC3)) - aC4 * aT1 * tauT1 * aT2 * tauT2 * tauC3 * np.exp(-t / tauT2) / (
                        (tauT1 - tauT2) * (tauT1 - tauC3) * (
                            tauT2 - tauNonCO2)) + aC4 * aT1 * tauT1 * aT2 * tauT2 * tauC3 * np.exp(-t / tauT1) / (
                        (tauT1 - tauT2) * (tauT1 - tauC3) * (
                            tauT2 - tauNonCO2)) + aC4 * aT1 * tauT1 * aT2 * tauC3 ** 2 * np.exp(-t / tauC3) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauC3)) + aC4 * aT1 * tauT1 * aT2 * tauC3 ** 2 * np.exp(-t / tauT1) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                            -tauT2 + tauC3)) + aC4 * aT1 * tauT1 * aT2 * tauC3 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                            tauT2 - tauNonCO2)) + aC4 * aT1 * tauT1 * aT2 * tauC3 * tauNonCO2 * np.exp(-t / tauT1) / (
                        (tauT1 - tauC3) * (tauT1 - tauNonCO2) * (
                            -tauT2 + tauNonCO2)) + aC4 * aT1 * aT2 * tauT2 * tauC3 ** 2 * np.exp(-t / tauC3) / (
                        (tauT1 - tauC3) * (tauT2 - tauC3) * (
                            tauT2 - tauNonCO2)) + aC4 * aT1 * aT2 * tauT2 * tauC3 ** 2 * np.exp(-t / tauT2) / (
                        (-tauT1 + tauC3) * (tauT2 - tauC3) * (
                            tauT2 - tauNonCO2)) + aC4 * aT1 * aT2 * tauT2 * tauC3 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((tauT1 - tauNonCO2) * (tauT2 - tauC3) * (
                tauT2 - tauNonCO2)) + aC4 * aT1 * aT2 * tauT2 * tauC3 * tauNonCO2 * np.exp(-t / tauT2) / (
                        (-tauT1 + tauNonCO2) * (tauT2 - tauC3) * (
                            tauT2 - tauNonCO2)) + aC4 * aT1 * aT2 * tauC3 ** 2 * tauNonCO2 * np.exp(-t / tauC3) / (
                        (-tauT1 + tauNonCO2) * (-tauT2 + tauC3) * (
                            tauC3 - tauNonCO2)) + aC4 * aT1 * aT2 * tauC3 ** 2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (-tauT1 + tauNonCO2) * (-tauT2 + tauC3) * (
                            -tauC3 + tauNonCO2)) + aC4 * aT1 * aT2 * tauC3 ** 2 * tauNonCO2 * np.exp(-t / tauC3) / (
                        (-tauT1 + tauC3) * (-tauT2 + tauNonCO2) * (
                            tauC3 - tauNonCO2)) + aC4 * aT1 * aT2 * tauC3 ** 2 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (-tauT1 + tauC3) * (-tauT2 + tauNonCO2) * (
                            -tauC3 + tauNonCO2)) + aC4 * aT2 ** 2 * tauT2 * tauC3 ** 2 * np.exp(-t / tauC3) / (
                        (tauT2 - tauC3) ** 2 * (tauT2 - tauNonCO2)) - aC4 * aT2 ** 2 * tauT2 * tauC3 ** 2 * np.exp(
        -t / tauT2) / ((tauT2 - tauC3) ** 2 * (
                tauT2 - tauNonCO2)) + aC4 * aT2 ** 2 * tauT2 * tauC3 * tauNonCO2 * np.exp(-t / tauNonCO2) / (
                        (tauT2 - tauC3) * (
                            tauT2 - tauNonCO2) ** 2) - aC4 * aT2 ** 2 * tauT2 * tauC3 * tauNonCO2 * np.exp(
        -t / tauT2) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) ** 2) + aC4 * aT2 ** 2 * t * tauC3 * np.exp(-t / tauT2) / (
                        (tauT2 - tauC3) * (tauT2 - tauNonCO2)) - aC4 * aT2 ** 2 * tauC3 ** 2 * tauNonCO2 * np.exp(
        -t / tauNonCO2) / ((tauT2 - tauC3) * (tauT2 - tauNonCO2) * (
                tauC3 - tauNonCO2)) + aC4 * aT2 ** 2 * tauC3 ** 2 * tauNonCO2 * np.exp(-t / tauC3) / (
                        (tauT2 - tauC3) * (tauT2 - tauNonCO2) * (tauC3 - tauNonCO2)) + (
                        aC1 * aS1 * aT2 ** 2 * tauNonCO2 - aC1 * aS1 * aT2 ** 2 * tauNonCO2 * np.exp(-t / tauS1)) / (
                        tauT2 - tauNonCO2) + (
                        aC1 * aS2 * aT2 ** 2 * tauNonCO2 - aC1 * aS2 * aT2 ** 2 * tauNonCO2 * np.exp(-t / tauS2)) / (
                        tauT2 - tauNonCO2) + (
                        aC1 * aS3 * aT2 ** 2 * tauNonCO2 - aC1 * aS3 * aT2 ** 2 * tauNonCO2 * np.exp(-t / tauS3)) / (
                        tauT2 - tauNonCO2) + (
                        aC1 * aS1 * aT1 ** 2 * tauNonCO2 - aC1 * aS1 * aT1 ** 2 * tauNonCO2 * np.exp(-t / tauS1)) / (
                        tauT1 - tauNonCO2) + (
                        aC1 * aS2 * aT1 ** 2 * tauNonCO2 - aC1 * aS2 * aT1 ** 2 * tauNonCO2 * np.exp(-t / tauS2)) / (
                        tauT1 - tauNonCO2) + (
                        aC1 * aS3 * aT1 ** 2 * tauNonCO2 - aC1 * aS3 * aT1 ** 2 * tauNonCO2 * np.exp(-t / tauS3)) / (
                        tauT1 - tauNonCO2)

    return (MCO2 / MC) * 1E12 * AANonCO2 * AACO2 * gamma * kPulseT ** 2 * tauNonCO2 * (
                term1 + term2 + term3 + term4 + term5)


def AGTPNonCO2_Final(t, tauNonCO2, kPulseT, aT1, tauT1, aT2, tauT2, AANonCO2, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3,
                     AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3):
    return AGTPNonCO2(t, tauNonCO2, kPulseT, aT1, tauT1, aT2, tauT2, AANonCO2) + DAGTPNonCO2(t, tauNonCO2, kPulseT, aT1,
                                                                                             tauT1, aT2, tauT2,
                                                                                             AANonCO2, aC1, aC2, aC3,
                                                                                             aC4, tauC1, tauC2, tauC3,
                                                                                             AACO2, gamma, aS1, aS2,
                                                                                             aS3, tauS1, tauS2, tauS3)


def DAGWPCH4toCO2(t, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2):
    return AACO2*MCO2*Y*(aC1*tauOH*(t + tauOH*(-1 + np.exp(-t/tauOH))) - aC2*tauC1**2*tauOH*(-np.exp(-t/tauOH) + np.exp(-t/tauC1))/(tauC1 - tauOH) + aC2*tauC1*(tauOH - tauOH*np.exp(-t/tauOH)) - aC3*tauC2**2*tauOH*(-np.exp(-t/tauOH) + np.exp(-t/tauC2))/(tauC2 - tauOH) + aC3*tauC2*(tauOH - tauOH*np.exp(-t/tauOH)) - aC4*tauC3**2*tauOH*(-np.exp(-t/tauOH) + np.exp(-t/tauC3))/(tauC3 - tauOH) + aC4*tauC3*(tauOH - tauOH*np.exp(-t/tauOH)))/(MCH4*tauOH)


def DAGTPCH4toCO2(t, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3, kPulseT, aT1, tauT1, aT2, tauT2, AACO2):
    return AACO2*MCO2*Y*kPulseT*(-aC1*aT1*tauOH*tauT1*(-np.exp(-t/tauT1) + np.exp(-t/tauOH))/(tauOH - tauT1) + aC1*aT1*(tauOH - tauOH*np.exp(-t/tauOH)) - aC1*aT2*tauOH*tauT2*(-np.exp(-t/tauT2) + np.exp(-t/tauOH))/(tauOH - tauT2) + aC1*aT2*(tauOH - tauOH*np.exp(-t/tauOH)) - aC2*aT1*tauC1**2*tauOH*(-np.exp(-t/tauOH) + np.exp(-t/tauC1))/((-tauC1 + tauT1)*(tauC1 - tauOH)) + aC2*aT1*tauC1*tauOH*tauT1*(-np.exp(-t/tauT1) + np.exp(-t/tauOH))/((-tauC1 + tauT1)*(tauOH - tauT1)) - aC2*aT2*tauC1**2*tauOH*(-np.exp(-t/tauOH) + np.exp(-t/tauC1))/((-tauC1 + tauT2)*(tauC1 - tauOH)) + aC2*aT2*tauC1*tauOH*tauT2*(-np.exp(-t/tauT2) + np.exp(-t/tauOH))/((-tauC1 + tauT2)*(tauOH - tauT2)) - aC3*aT1*tauC2**2*tauOH*(-np.exp(-t/tauOH) + np.exp(-t/tauC2))/((-tauC2 + tauT1)*(tauC2 - tauOH)) + aC3*aT1*tauC2*tauOH*tauT1*(-np.exp(-t/tauT1) + np.exp(-t/tauOH))/((-tauC2 + tauT1)*(tauOH - tauT1)) - aC3*aT2*tauC2**2*tauOH*(-np.exp(-t/tauOH) + np.exp(-t/tauC2))/((-tauC2 + tauT2)*(tauC2 - tauOH)) + aC3*aT2*tauC2*tauOH*tauT2*(-np.exp(-t/tauT2) + np.exp(-t/tauOH))/((-tauC2 + tauT2)*(tauOH - tauT2)) - aC4*aT1*tauC3**2*tauOH*(-np.exp(-t/tauOH) + np.exp(-t/tauC3))/((-tauC3 + tauT1)*(tauC3 - tauOH)) + aC4*aT1*tauC3*tauOH*tauT1*(-np.exp(-t/tauT1) + np.exp(-t/tauOH))/((-tauC3 + tauT1)*(tauOH - tauT1)) - aC4*aT2*tauC3**2*tauOH*(-np.exp(-t/tauOH) + np.exp(-t/tauC3))/((-tauC3 + tauT2)*(tauC3 - tauOH)) + aC4*aT2*tauC3*tauOH*tauT2*(-np.exp(-t/tauT2) + np.exp(-t/tauOH))/((-tauC3 + tauT2)*(tauOH - tauT2)))/(MCH4*tauOH)


def AGTPCH4Fossil_Final(t, tauNonCO2, kPulseT, aT1, tauT1, aT2, tauT2, AANonCO2, aC1, aC2, aC3, aC4, tauC1, tauC2,
                        tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3):

    return AGTPNonCO2_Final(t, tauNonCO2, kPulseT, aT1, tauT1, aT2, tauT2, AANonCO2, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3) + DAGTPCH4toCO2(t, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3, kPulseT, aT1, tauT1, aT2, tauT2, AACO2)


def AGTPCH4NonFossil_Final(t, tauNonCO2, kPulseT, aT1, tauT1, aT2, tauT2, AANonCO2, aC1, aC2, aC3, aC4, tauC1, tauC2,
                           tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3):
    return AGTPNonCO2_Final(t, tauNonCO2, kPulseT, aT1, tauT1, aT2, tauT2, AANonCO2, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3) + DAGTPCH4toCO2(t, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3, kPulseT, aT1, tauT1, aT2, tauT2, AACO2) - (MCO2/MCH4)*AGTPCO2(t, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3, kPulseT, aT1, tauT1, aT2, tauT2, AACO2)


def produce_simplified_version(complete_dataframe):
    """
    Method producing the simplified version of IW+ in which there are only 6 indicators:
    - Climate change, short term (=GWP100)
    - Water scarcity
    - Fossil resources
    - Remaining Human health damage
    - Remaining Ecosystem quality damage
    The damage impacts of climate change and water use on the Human health and Ecosystem quality area of protections
    are removed, hence only leaving "Remaining [...] damage". These are removed to avoid potential misuse by users
    where the impact of climate change and water use would be double-counted. Once as a midpoint category, and
    another time as a damage category.
    :return:
    """

    simplified_version = complete_dataframe.copy('deep')

    # midpoint categories excluded from simplified version
    midpoint_drop = ['Climate change, long term', 'Freshwater acidification', 'Freshwater ecotoxicity',
                     'Freshwater eutrophication', 'Mineral resources use', 'Human toxicity cancer',
                     'Human toxicity non-cancer', 'Ionizing radiations', 'Land occupation, biodiversity',
                     'Land transformation, biodiversity', 'Marine eutrophication', 'Ozone layer depletion',
                     'Particulate matter formation', 'Photochemical oxidant formation', 'Terrestrial acidification']
    # endpoint categories excluded from simplified version
    endpoint_drop = ['Climate change, ecosystem quality, long term',
                     'Climate change, ecosystem quality, short term',
                     'Climate change, human health, long term', 'Climate change, human health, short term',
                     'Water availability, freshwater ecosystem', 'Water availability, human health',
                     'Water availability, terrestrial ecosystem',
                     'Human toxicity cancer, long term', 'Human toxicity non-cancer, long term',
                     'Freshwater ecotoxicity, long term', 'Marine acidification, long term']
    # dropping midpoint_drop
    simplified_version.drop([i for i in simplified_version.index if (
            simplified_version.loc[i, 'Impact category'] in midpoint_drop and
            simplified_version.loc[i, 'MP or Damage'] == 'Midpoint')],
                            inplace=True)
    # dropping endpoint_drop
    simplified_version.drop([i for i in simplified_version.index if
                             simplified_version.loc[i, 'Impact category'] in endpoint_drop], inplace=True)
    # storing the cas number to put them back at the end
    cas = simplified_version[['Elem flow name', 'CAS number']].drop_duplicates().set_index('Elem flow name').to_dict()[
        'CAS number']

    # setting index to allow use of groupby later
    simplified_version = simplified_version.set_index(['CF unit', 'Compartment', 'Sub-compartment',
                                                       'Elem flow name', 'Elem flow unit', 'MP or Damage'])
    # isolate and group HH CFs
    hh_simplified = simplified_version.loc['DALY'].copy()
    hh_simplified = hh_simplified.drop(['Impact category', 'CAS number',
                                        'Native geographical resolution scale'], axis=1).groupby(
        hh_simplified.index).sum()
    hh_simplified.index = pd.MultiIndex.from_tuples(hh_simplified.index)
    # isolate and group EQ CFs
    eq_simplified = simplified_version.loc['PDF.m2.yr'].copy()
    eq_simplified = eq_simplified.drop(['Impact category', 'CAS number',
                                        'Native geographical resolution scale'], axis=1).groupby(
        eq_simplified.index).sum()
    eq_simplified.index = pd.MultiIndex.from_tuples(eq_simplified.index)
    # delete HH and EQ CFs from original df
    simplified_version.drop(['DALY', 'PDF.m2.yr'], inplace=True)
    simplified_version = simplified_version.reset_index()
    # make hh_simplified respect the format of self.simplified_version_sp for concatenation
    hh_simplified = hh_simplified.reset_index()
    hh_simplified.columns = ['Compartment', 'Sub-compartment', 'Elem flow name', 'Elem flow unit',
                             'MP or Damage', 'CF value']
    hh_simplified.loc[:, 'CF unit'] = 'DALY'
    hh_simplified.loc[:, 'Impact category'] = 'Remaining Human health damage'
    # make eq_simplified respect the format of self.simplified_version_sp for concatenation
    eq_simplified = eq_simplified.reset_index()
    eq_simplified.columns = ['Compartment', 'Sub-compartment', 'Elem flow name', 'Elem flow unit',
                             'MP or Damage', 'CF value']
    eq_simplified.loc[:, 'CF unit'] = 'PDF.m2.yr'
    eq_simplified.loc[:, 'Impact category'] = 'Remaining Ecosystem quality damage'
    # concat
    simplified_version = pd.concat([simplified_version, hh_simplified, eq_simplified])
    # put back the CAS numbers
    simplified_version['CAS number'] = [cas[i] for i in simplified_version['Elem flow name']]

    simplified_version = clean_up_dataframe(simplified_version)

    simplified_version.loc[[i for i in simplified_version.index if simplified_version.loc[i, 'Impact category'] ==
                            'Climate change, short term'], 'Impact category'] = 'Carbon footprint'
    simplified_version.loc[[i for i in simplified_version.index if simplified_version.loc[i, 'Impact category'] ==
                            'Water scarcity'], 'Impact category'] = 'Water scarcity footprint'

    return simplified_version


def produce_simplified_version_olca(complete_dataframe):

    simplified = complete_dataframe.copy('deep')

    # midpoint categories excluded from simplified version
    midpoint_drop = ['Climate change, long term', 'Freshwater acidification', 'Freshwater ecotoxicity',
                     'Freshwater eutrophication', 'Mineral resources use', 'Human toxicity cancer',
                     'Human toxicity non-cancer', 'Ionizing radiations', 'Land occupation, biodiversity',
                     'Land transformation, biodiversity', 'Marine eutrophication', 'Ozone layer depletion',
                     'Particulate matter formation', 'Photochemical oxidant formation', 'Terrestrial acidification']
    # endpoint categories excluded from simplified version
    endpoint_drop = ['Climate change, ecosystem quality, long term',
                     'Climate change, ecosystem quality, short term',
                     'Climate change, human health, long term', 'Climate change, human health, short term',
                     'Water availability, freshwater ecosystem', 'Water availability, human health',
                     'Water availability, terrestrial ecosystem',
                     'Human toxicity cancer, long term', 'Human toxicity non-cancer, long term',
                     'Freshwater ecotoxicity, long term', 'Marine acidification, long term']

    simplified = simplified.drop([ic for ic in simplified.index if (ic[0] in midpoint_drop and
                                                                    ic[1] not in ['DALY', 'PDF.m2.yr'])])
    simplified = simplified.drop([ic for ic in simplified.index if (ic[0] in endpoint_drop and
                                                                    ic[1] in ['DALY', 'PDF.m2.yr'])])
    hh_simplified = simplified.loc(axis=0)[:, 'DALY'].groupby(axis=0, level=2).agg({'flow_name': 'first',
                                                                                    'comp': 'first',
                                                                                    'cas': 'first',
                                                                                    'unit': 'first',
                                                                                    'CF value': sum,
                                                                                    'MP or Damage': 'first'})
    hh_simplified.index = pd.MultiIndex.from_product([['Remaining Human health damage'], ['DALY'], hh_simplified.index])

    eq_simplified = simplified.loc(axis=0)[:, 'PDF.m2.yr'].groupby(axis=0, level=2).agg({'flow_name': 'first',
                                                                                         'comp': 'first',
                                                                                         'cas': 'first',
                                                                                         'unit': 'first',
                                                                                         'CF value': sum,
                                                                                         'MP or Damage': 'first'})
    eq_simplified.index = pd.MultiIndex.from_product(
        [['Remaining Ecosystem quality damage'], ['PDF.m2.yr'], eq_simplified.index])

    simplified = simplified.drop('DALY', axis=0, level=1)
    simplified = simplified.drop('PDF.m2.yr', axis=0, level=1)

    simplified = pd.concat([simplified, hh_simplified, eq_simplified])

    simplified.index = [('Carbon footprint', i[1], i[2]) if i[0] == 'Climate change, short term' else i
                        for i in simplified.index]
    simplified.index = [('Water scarcity footprint', i[1], i[2]) if i[0] == 'Water scarcity' else i
                        for i in simplified.index]
    simplified.index = pd.MultiIndex.from_tuples(simplified.index)

    return simplified


def clean_up_dataframe(df):
    # remove duplicates
    df = df.drop_duplicates()
    # fix index
    df = df.reset_index().drop('index', axis=1)
    return df


def convert_country_codes(db):
    """
    Converts 3 letters ISO codes to 2 letters ISO codes, which are used by ecoinvent.
    :param db:
    :return:
    """
    # lines of logger to avoid having warnings showing up
    coco_logger = coco.logging.getLogger()
    coco_logger.setLevel(coco.logging.CRITICAL)
    # converting 3 letters ISO codes to 2 letters ISO codes, which are used by ecoinvent
    db.loc[[i for i in db.index if db.loc[i, 'Resolution'] == 'Country'], 'Region code'] = coco.convert(
        db.loc[[i for i in db.index if db.loc[i, 'Resolution'] == 'Country'], 'Region code'], to='ISO2', not_found=None)

    return db


def mapping_with_ei35():
    """
    Support function kept for transparency on how the mapping with ecoinvent was performed.
    :return: nothing
    """
    with gzip.open( '.../ecoinvent3.5.cutoffPandas_symmNorm.gz.pickle','rb') as f:
        eco_35 = pd.read_pickle(f)

    # check which substances have a direct match with their names
    matching_by_name = self.master_db.loc[:, ['Elem flow name']].merge(eco_35['STR'].loc[:, ['name']],
                                                                     left_on='Elem flow name',
                                                                     right_on='name', how='inner').drop_duplicates()
    # harmonize CAS numbers
    eco_35['STR'].cas = [i.lstrip('0') if i != None else None for i in eco_35['STR'].cas]
    corrected_cas = []
    for i in self.master_db.loc[:, 'CAS number']:
        if type(i) != float and i != None:
            corrected_cas.append(i.lstrip('0'))
        else:
            corrected_cas.append(None)
    self.master_db.loc[:, 'CAS number'] = corrected_cas
    # merge dataframes
    matching_by_cas = self.master_db.loc[:, ['CAS number', 'Elem flow name']].merge(eco_35['STR'].loc[:, ['name', 'cas']],
                                                                                  left_on='CAS number',
                                                                                  right_on='cas',
                                                                                  how='inner').drop_duplicates().dropna()
    # if there is a duplicate it means 2 different ecoinvent names are linked to a single IW name
    # can't have that, it will require matching by hand to sort it out
    issues_matching_cas = list(
        set(matching_by_cas[matching_by_cas.loc[:, 'Elem flow name'].duplicated()].loc[:, 'Elem flow name']))
    matching_by_cas = matching_by_cas.loc[
        [i for i in matching_by_cas.index if matching_by_cas.loc[i, 'Elem flow name'] not in issues_matching_cas]]

    auto_matching = pd.concat([matching_by_name, matching_by_cas])
    auto_matching.drop(['CAS number', 'cas'], axis=1, inplace=True)

    # plus add manual mapping later on


def mapping_with_sp():
    # determine all substance within SP (multiple databases)
    def add_flows_dbs(sp_emissions, path_db_file):
        db_flows = pd.read_excel(path_db_file)

        resources_indexes = db_flows.loc[
            [i for i in db_flows.index if db_flows.loc[i, 'SimaPro 9.3.0.3'] == 'Resources']].index
        materials_indexes = db_flows.loc[
            [i for i in db_flows.index if db_flows.loc[i, 'SimaPro 9.3.0.3'] == 'Materials/fuels']].index
        air_indexes = db_flows.loc[
            [i for i in db_flows.index if db_flows.loc[i, 'SimaPro 9.3.0.3'] == 'Emissions to air']].index
        water_indexes = db_flows.loc[
            [i for i in db_flows.index if db_flows.loc[i, 'SimaPro 9.3.0.3'] == 'Emissions to water']].index
        soil_indexes = db_flows.loc[
            [i for i in db_flows.index if db_flows.loc[i, 'SimaPro 9.3.0.3'] == 'Emissions to soil']].index
        waste_indexes = db_flows.loc[
            [i for i in db_flows.index if db_flows.loc[i, 'SimaPro 9.3.0.3'] == 'Final waste flows']].index

        for i in range(0, len(resources_indexes)):
            sp_emissions = pd.concat([sp_emissions, db_flows.loc[
                                                    air_indexes.values.tolist()[i] + 1:water_indexes.values.tolist()[
                                                                                           i] - 1].loc[:,
                                                    'SimaPro 9.3.0.3']])
            sp_emissions = pd.concat([sp_emissions, db_flows.loc[
                                                    water_indexes.values.tolist()[i] + 1:soil_indexes.values.tolist()[
                                                                                             i] - 1].loc[:,
                                                    'SimaPro 9.3.0.3']])
            sp_emissions = pd.concat([sp_emissions, db_flows.loc[
                                                    soil_indexes.values.tolist()[i] + 1:waste_indexes.values.tolist()[
                                                                                            i] - 1].loc[:,
                                                    'SimaPro 9.3.0.3']])
            sp_emissions = pd.concat([sp_emissions, db_flows.loc[
                                                    resources_indexes.values.tolist()[i] + 1:
                                                    materials_indexes.values.tolist()[i] - 1].loc[:, 'SimaPro 9.3.0.3']])

        return sp_emissions

    sp_emissions = pd.DataFrame()
    sp_emissions = add_flows_dbs(sp_emissions, 'C://Users/11max/PycharmProjects/IW_Reborn/Data/mappings/SP/ecoinvent.XLSX')
    sp_emissions = add_flows_dbs(sp_emissions, 'C://Users/11max/PycharmProjects/IW_Reborn/Data/mappings/SP/agribalyse.XLSX')
    sp_emissions = add_flows_dbs(sp_emissions,
                                 'C://Users/11max/PycharmProjects/IW_Reborn/Data/mappings/SP/agrifootprint.XLSX')
    sp_emissions = add_flows_dbs(sp_emissions, 'C://Users/11max/PycharmProjects/IW_Reborn/Data/mappings/SP/ELCD.XLSX')
    sp_emissions = add_flows_dbs(sp_emissions,
                                 'C://Users/11max/PycharmProjects/IW_Reborn/Data/mappings/SP/Industry2.0.XLSX')
    sp_emissions = add_flows_dbs(sp_emissions, 'C://Users/11max/PycharmProjects/IW_Reborn/Data/mappings/SP/US-ei2.2.XLSX')
    sp_emissions = add_flows_dbs(sp_emissions, 'C://Users/11max/PycharmProjects/IW_Reborn/Data/mappings/SP/USLCI.XLSX')
    sp_emissions = add_flows_dbs(sp_emissions, 'C://Users/11max/PycharmProjects/IW_Reborn/Data/mappings/SP/WFDB.XLSX')

    sp_emissions = sp_emissions.dropna().drop_duplicates()


def extract_olca_flowsv1():
    """
    Works for openLCA v1.x versions
    :return:
    """
    # export relevant information from json files into a single dict
    path = 'C://Users/11max/PycharmProjects/IW_Reborn/Data/metadata/OLCA/ei38/flows/'
    directory = os.fsencode(path)

    all_flows = []
    for file in os.listdir(directory):
        filename = os.fsdecode(file)

        with open(path + filename, 'r') as f:
            read = json.load(f)

        if 'cas' in read.keys():
            flow = {'flow_name': read['name'], 'flow_id': read['@id'], 'comp': read['category']['@id'],
                    'cas': read['cas'], 'unit': read['flowProperties'][0]['flowProperty']['@id']}
        else:
            flow = {'flow_name': read['name'], 'flow_id': read['@id'], 'comp': read['category']['@id'],
                    'cas': None, 'unit': read['flowProperties'][0]['flowProperty']['@id']}

        all_flows.append(flow)

    # now we go fetch information on categories and units from the other json files
    path = 'C://Users/11max/PycharmProjects/IW_Reborn/Data/metadata/OLCA/ei38/categories/'
    directory = os.fsencode(path)

    categories = {}

    for file in os.listdir(directory):
        filename = os.fsdecode(file)

        with open(path + filename, 'r') as f:
            read = json.load(f)

        if 'category' in read.keys():
            categories[read['@id']] = (read['category']['name'], read['name'])

    # hardcoded because it's too much of a mess inside openLCA files
    units = {'4e10f566-0358-489a-8e3a-d687b66c50e6': 'kg*a', 'f6811440-ee37-11de-8a39-0800200c9a66': 'MJ',
             '93a60a56-a3c8-22da-a746-0800200c9a66': 'm3', '93a60a56-a3c8-17da-a746-0800200c9a66': 'kBq',
             '93a60a56-a3c8-19da-a746-0800200c9a66': 'm2', '93a60a56-a3c8-21da-a746-0800200c9a66': 'm2*a',
             '441238a3-ba09-46ec-b35b-c30cfba746d1': 'm3*a', '93a60a56-a3c8-11da-a746-0800200b9a66': 'kg'}

    all_flows = pd.DataFrame(all_flows)
    all_flows['comp'] = [categories[i] for i in all_flows['comp']]
    all_flows['unit'] = [units[i] for i in all_flows['unit']]

    return all_flows


def extract_olca_flowsv2():
    """
    Works for openLCA v2.x
    :return:
    """
    # export relevant information from json files into a single dict
    path = 'C://Users/11max/PycharmProjects/IW_Reborn/Data/metadata/OLCA/v2.0/ei38/flows/'
    directory = os.fsencode(path)

    all_flows = []
    for file in os.listdir(directory):
        filename = os.fsdecode(file)

        with open(path + filename, 'r') as f:
            read = json.load(f)

        if 'cas' in read.keys():
            flow = {'flow_name': read['name'], 'flow_id': read['@id'], 'comp': read['category'],
                    'cas': read['cas'], 'unit': read['flowProperties'][0]['flowProperty']['refUnit']}
        else:
            flow = {'flow_name': read['name'], 'flow_id': read['@id'], 'comp': read['category'],
                    'cas': None, 'unit': read['flowProperties'][0]['flowProperty']['refUnit']}

        all_flows.append(flow)

    return pd.DataFrame(all_flows)
