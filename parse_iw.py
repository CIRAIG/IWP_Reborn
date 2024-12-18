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
from tqdm import tqdm


class Parse:
    def __init__(self, path_access_db, version, bw2_projects):
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
            - load_water_availability_fw_cfs()
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
        self.version = str(version)
        self.bw2_projects = bw2_projects

        # OUTPUTs
        self.master_db = pd.DataFrame()
        self.master_db_carbon_neutrality = pd.DataFrame()
        self.master_db_not_regio = pd.DataFrame()
        self.master_db_not_regio_carbon_neutrality = pd.DataFrame()
        self.ei38_iw = pd.DataFrame()
        self.ei38_iw_carbon_neutrality = pd.DataFrame()
        self.ei39_iw = pd.DataFrame()
        self.ei39_iw_carbon_neutrality = pd.DataFrame()
        self.ei310_iw = pd.DataFrame()
        self.ei310_iw_carbon_neutrality = pd.DataFrame()
        self.simplified_version_ei38 = pd.DataFrame()
        self.simplified_version_ei39 = pd.DataFrame()
        self.simplified_version_ei310 = pd.DataFrame()
        self.iw_sp = pd.DataFrame()
        self.iw_sp_carbon_neutrality = pd.DataFrame()
        self.simplified_version_sp = pd.DataFrame()
        self.simplified_version_olca = pd.DataFrame()
        self.simplified_version_bw = pd.DataFrame()
        self.sp_data = {}
        self.olca_iw = pd.DataFrame()
        self.olca_iw_carbon_neutrality = pd.DataFrame()
        self.olca_data = {}
        self.olca_data_custom = {}
        self.exio_iw = pd.DataFrame()

        self.conn = sqlite3.connect(self.path_access_db)

# -------------------------------------------- Main methods ------------------------------------------------------------

    def load_cfs(self):
        """
        Load the characterization factors and stored them in master_db.
        :return: updated master_db
        """

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
        self.load_water_availability_fw_cfs()
        self.load_water_availability_hh_cfs()
        self.load_water_availability_terr_cfs()
        self.logger.info("Loading thermally polluted water characterization factors...")
        self.load_thermally_polluted_water_cfs()
        self.logger.info("Loading plastic physical effects on biota characterization factors...")
        self.load_plastic_cfs()
        self.logger.info("Loading fisheries impact characterization factors...")
        self.load_fisheries_cfs()

        self.logger.info("Harmonizing regionalized substances across indicators...")
        self.harmonize_regionalized_substances()

        self.logger.info("Applying rules...")
        self.apply_rules()

        self.logger.info("Treating regionalized factors...")
        self.create_not_regio_flows()
        self.create_regio_flows_for_not_regio_ic()
        self.order_things_around()

        self.logger.info("Managing biogenic carbon shenanigans...")
        self.deal_with_biogenic_carbon()
        self.deal_with_temporary_storage_of_carbon()
        self.separate_ghg_indicators()

        self.logger.info("Create non-regionalized version for ecoinvent...")
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

    def export_to_bw2(self):
        """
        This method creates a brightway2 method with the IW+ characterization factors.
        :param ei_flows_version: [str] Provide a specific ei version (e.g., 3.6) to be used to determine the elementary flows
                                 to be linked to iw+. Default values = eiv3.8 (in 2022)

        :return:
        """

        self.logger.info("Exporting to brightway2...")

        for project in self.bw2_projects:
            bw2.projects.set_current(project)
            bio = bw2.Database('biosphere3')
            ei_version = project.split('ecoinvent')[1]
            bw_flows_with_codes = (
                pd.DataFrame(
                    [(i.as_dict()['name'], i.as_dict()['categories'][0], i.as_dict()['categories'][1],
                      i.as_dict()['code'])
                     if len(i.as_dict()['categories']) == 2
                     else (i.as_dict()['name'], i.as_dict()['categories'][0], 'unspecified', i.as_dict()['code'])
                     for i in bio],
                    columns=['Elem flow name', 'Compartment', 'Sub-compartment', 'code'])
            )
            if project == 'ecoinvent3.8':
                ei_in_bw_normal = self.ei38_iw.merge(bw_flows_with_codes)
                ei_in_bw_carbon_neutrality = self.ei38_iw_carbon_neutrality.merge(bw_flows_with_codes)
                ei_in_bw_simple = self.simplified_version_ei38.merge(bw_flows_with_codes)
            elif project == 'ecoinvent3.9':
                ei_in_bw_normal = self.ei39_iw.merge(bw_flows_with_codes)
                ei_in_bw_carbon_neutrality = self.ei39_iw_carbon_neutrality.merge(bw_flows_with_codes)
                ei_in_bw_simple = self.simplified_version_ei39.merge(bw_flows_with_codes)
            elif project == 'ecoinvent3.10':
                ei_in_bw_normal = self.ei310_iw.merge(bw_flows_with_codes)
                ei_in_bw_carbon_neutrality = self.ei310_iw_carbon_neutrality.merge(bw_flows_with_codes)
                ei_in_bw_simple = self.simplified_version_ei310.merge(bw_flows_with_codes)

            for ei_in_bw_format in ['normal', 'carbon neutrality']:
                if ei_in_bw_format == 'normal':
                    ei_in_bw = ei_in_bw_normal
                elif ei_in_bw_format == 'carbon neutrality':
                    ei_in_bw = ei_in_bw_carbon_neutrality
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
                impact_categories = ei_in_bw.index.drop_duplicates()

                # -------------- For complete version of IW+ ----------------
                for ic in impact_categories:
                    if ei_in_bw.loc[[ic], 'MP or Damage'].iloc[0] == 'Midpoint':
                        mid_end = 'Midpoint'
                        if ei_in_bw_format == 'normal':
                            name = ('IMPACT World+ ' + mid_end + ' ' + self.version + ' for ecoinvent v' +
                                    ei_version + ' (incl. CO2 uptake)', 'Midpoint', ic[0])
                        elif ei_in_bw_format == 'carbon neutrality':
                            name = ('IMPACT World+ ' + mid_end + ' ' + self.version + ' for ecoinvent v' +
                                    ei_version, 'Midpoint', ic[0])
                    else:
                        mid_end = 'Damage'
                        if ic[1] == 'DALY':
                            if ei_in_bw_format == 'normal':
                                name = ('IMPACT World+ ' + mid_end + ' ' + self.version + ' for ecoinvent v' +
                                        ei_version + ' (incl. CO2 uptake)', 'Human health', ic[0])
                            elif ei_in_bw_format == 'carbon neutrality':
                                name = ('IMPACT World+ ' + mid_end + ' ' + self.version + ' for ecoinvent v' +
                                        ei_version, 'Human health', ic[0])
                        else:
                            if ei_in_bw_format == 'normal':
                                name = ('IMPACT World+ ' + mid_end + ' ' + self.version + ' for ecoinvent v' +
                                        ei_version + ' (incl. CO2 uptake)', 'Ecosystem quality', ic[0])
                            elif ei_in_bw_format == 'carbon neutrality':
                                name = ('IMPACT World+ ' + mid_end + ' ' + self.version + ' for ecoinvent v' +
                                        ei_version, 'Ecosystem quality', ic[0])

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

            ei_in_bw_simple.set_index(['Impact category', 'CF unit'], inplace=True)
            impact_categories_simple = ei_in_bw_simple.index.drop_duplicates()

            for ic in impact_categories_simple:

                name = ('IMPACT World+ Footprint ' + self.version + ' for ecoinvent v' + ei_version, ic[0])

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
        self.iw_sp_carbon_neutrality.loc[:, 'CF value'] = self.iw_sp_carbon_neutrality.loc[:, 'CF value'].astype(str)
        self.simplified_version_sp.loc[:, 'CF value'] = self.simplified_version_sp.loc[:, 'CF value'].astype(str)

        # Metadata
        l = ['SimaPro 9.6', 'methods', 'Date: ' + datetime.now().strftime("%D"),
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
                                    ['IMPACT World+ Midpoint ' + self.version + ' (incl. CO2 uptake)', '', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                    ['2','1', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                    ['IMPACT World+ Midpoint ' + self.version + ' (incl. CO2 uptake)' + chr(int("007F", 16)) + chr(int("007F", 16)) +
                                     'New category:' + chr(int("007F", 16)) +
                                     '- Plastics physical effect on biota' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'Updated categories:' + chr(int("007F", 16)) +
                                     '- Climate change indicators with -1/+1 approach biogenic carbon' + chr(int("007F", 16)) +
                                     '- Fossil and nuclear energy use'+chr(int("007F", 16)) +
                                     '- Ozone layer depletion' + chr(int("007F", 16)) +
                                     '- Particulate matter formation' + chr(int("007F", 16)) +
                                     '- Photochemical ozone formation' + chr(int("007F", 16)) +
                                     '- Water scarcity' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(int("007F", 16)) +
                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes',
                                     '', '', '', '', ''], ['', '', '', '', '', ''],
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
                                  ['IMPACT World+ Expert ' + self.version + ' (incl. CO2 uptake)', '', '', '', '', ''],
                                  ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                  ['2','1', '', '', '', ''],
                                  ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                  ['IMPACT World+ Expert ' + self.version + ' (incl. CO2 uptake)' + chr(int("007F", 16)) + chr(int("007F", 16)) +
                                   'New categories:' + chr(int("007F", 16)) + '- Marine ecotoxicity' + chr(int("007F", 16)) +
                                   '- Terrestrial ecotoxicity' + chr(int("007F", 16)) + '- Fisheries' + chr(int("007F", 16)) +
                                   '- Plastics physical effect on biota' + chr(int("007F", 16)) +
                                   '- Phochemical ozone formation, ecosystem quality' + chr(int("007F", 16)) + chr(int("007F", 16)) +
                                   'Updated categories:' + chr(int("007F", 16)) +
                                   '- Climate change indicators with -1/+1 approach biogenic carbon' + chr(int("007F", 16)) +
                                   '- Climate change damage indicators' + chr(int("007F", 16)) +
                                   '- Ozone layer depletion' + chr(int("007F", 16)) +
                                   '- Particulate matter formation' + chr(int("007F", 16)) +
                                   '- Photochemical ozone formation' + chr(int("007F", 16)) +
                                   '- Water availability, human health' + chr(int("007F", 16)) +
                                   '- Water availability, terrestrial ecosystems' + chr(int("007F", 16))
                                   + chr(int("007F", 16)) +
                                   'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(int("007F", 16)) +
                                   'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes',
                                   '', '', '', '', ''], ['', '', '', '', '', ''],
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
                                    ['IMPACT World+ ' + self.version + ' (incl. CO2 uptake)', '', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                    ['2','1', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                    ['IMPACT World+ ' + self.version  + ' (incl. CO2 uptake)' + chr(int("007F", 16)) +
                                     'New categories:' + chr(int("007F", 16)) +
                                     '- Marine ecotoxicity' + chr(int("007F", 16)) +
                                     '- Terrestrial ecotoxicity' + chr(int("007F", 16)) +
                                     '- Fisheries' + chr(int("007F", 16)) +
                                     '- Plastics physical effect on biota' + chr(int("007F", 16)) +
                                     '- Phochemical ozone formation, ecosystem quality' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'Updated categories:' + chr(int("007F", 16)) +
                                     '- Climate change indicators with -1/+1 approach biogenic carbon' + chr(int("007F", 16)) +
                                     '- Climate change damage indicators' + chr(int("007F", 16)) +
                                     '- Fossil and nuclear energy use' + chr(int("007F", 16)) +
                                     '- Ozone layer depletion' + chr(int("007F", 16)) +
                                     '- Particulate matter formation' + chr(int("007F", 16)) +
                                     '- Photochemical ozone formation' + chr(int("007F", 16)) +
                                     '- Water availability, human health' + chr(int("007F", 16)) +
                                     '- Water availability, terrestrial ecosystems' + chr(int("007F", 16)) +
                                     '- Water scarcity' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(int("007F", 16)) +
                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes',
                                     '', '', '', '', ''], ['', '', '', '', '', ''],
                                    ['Category', '', '', '', '', ''], ['Others', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Damage Assessment', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Normalization', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Weighting', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                    ['', '', '', '', '', ''],
                                    ['Use Addition', '', '', '', '', ''], ['Yes', '', '', '', '', '']]
        # metadata on the midpoint method
        midpoint_method_metadata_carboneutrality = [['Method', '', '', '', '', ''], ['', '', '', '', '', ''],
                                    ['Name', '', '', '', '', ''],
                                    ['IMPACT World+ Midpoint ' + self.version, '', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                    ['2','1', '', '', '', ''],
                                    ['IMPACT World+ Midpoint ' + self.version + chr(int("007F", 16)) + chr(int("007F", 16)) +
                                     'New category:' + chr(int("007F", 16)) +
                                     '- Plastics physical effect on biota' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'Updated categories:' + chr(int("007F", 16)) +
                                     '- Fossil and nuclear energy use'+chr(int("007F", 16)) +
                                     '- Ozone layer depletion' + chr(int("007F", 16)) +
                                     '- Particulate matter formation' + chr(int("007F", 16)) +
                                     '- Photochemical ozone formation' + chr(int("007F", 16)) +
                                     '- Water scarcity' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(int("007F", 16)) +
                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes',
                                     '', '', '', '', ''], ['', '', '', '', '', ''],
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
        damage_method_metadata_carboneutrality = [['Method', '', '', '', '', ''], ['', '', '', '', '', ''],
                                  ['Name', '', '', '', '', ''],
                                  ['IMPACT World+ Expert ' + self.version, '', '', '', '', ''],
                                  ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                  ['2','1', '', '', '', ''],
                                  ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                  ['IMPACT World+ Expert ' + self.version + chr(int("007F", 16)) + chr(int("007F", 16)) +
                                   'New categories:' + chr(int("007F", 16)) + '- Marine ecotoxicity' + chr(int("007F", 16)) +
                                   '- Terrestrial ecotoxicity' + chr(int("007F", 16)) + '- Fisheries' + chr(int("007F", 16)) +
                                   '- Plastics physical effect on biota' + chr(int("007F", 16)) +
                                   '- Phochemical ozone formation, ecosystem quality' + chr(int("007F", 16)) + chr(int("007F", 16)) +
                                   'Updated categories:' + chr(int("007F", 16)) +
                                   '- Climate change damage indicators' + chr(int("007F", 16)) +
                                   '- Ozone layer depletion' + chr(int("007F", 16)) +
                                   '- Particulate matter formation' + chr(int("007F", 16)) +
                                   '- Photochemical ozone formation' + chr(int("007F", 16)) +
                                   '- Water availability, human health' + chr(int("007F", 16)) +
                                   '- Water availability, terrestrial ecosystems' + chr(int("007F", 16))
                                   + chr(int("007F", 16)) +
                                   'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(int("007F", 16)) +
                                   'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes',
                                   '', '', '', '', ''], ['', '', '', '', '', ''],
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
        combined_method_metadata_carboneutrality = [['Method', '', '', '', '', ''], ['', '', '', '', '', ''],
                                    ['Name', '', '', '', '', ''],
                                    ['IMPACT World+ ' + self.version, '', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                    ['2','1', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                    ['IMPACT World+ ' + self.version + chr(int("007F", 16)) +
                                     'New categories:' + chr(int("007F", 16)) +
                                     '- Marine ecotoxicity' + chr(int("007F", 16)) +
                                     '- Terrestrial ecotoxicity' + chr(int("007F", 16)) +
                                     '- Fisheries' + chr(int("007F", 16)) +
                                     '- Plastics physical effect on biota' + chr(int("007F", 16)) +
                                     '- Phochemical ozone formation, ecosystem quality' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'Updated categories:' + chr(int("007F", 16)) +
                                     '- Climate change damage indicators' + chr(int("007F", 16)) +
                                     '- Fossil and nuclear energy use' + chr(int("007F", 16)) +
                                     '- Ozone layer depletion' + chr(int("007F", 16)) +
                                     '- Particulate matter formation' + chr(int("007F", 16)) +
                                     '- Photochemical ozone formation' + chr(int("007F", 16)) +
                                     '- Water availability, human health' + chr(int("007F", 16)) +
                                     '- Water availability, terrestrial ecosystems' + chr(int("007F", 16)) +
                                     '- Water scarcity' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(int("007F", 16)) +
                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes',
                                     '', '', '', '', ''], ['', '', '', '', '', ''],
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
                                    ['IMPACT World+ Footprint ' + self.version, '', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                    ['2','1', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                    ['IMPACT World+ Footprint ' + self.version + chr(int("007F", 16)) +
                                     'Updated categories:' + chr(int("007F", 16)) +
                                     '- Fossil and nuclear energy use' + chr(int("007F", 16)) +
                                     '- Water footprint - Scarcity' + chr(int("007F", 16)) +
                                     '- Human health (residual)' + chr(int("007F", 16)) +
                                     '- Ecosystem quality (residual)' + chr(int("007F", 16)) +
                                     'For details on what the footprint version of IW+ entails, please consult this page: '
                                     'https://www.impactworldplus.org/version-2-0-1/' + chr(int("007F", 16)) +
                                     'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.'  + chr(int("007F", 16)) +
                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes',
                                     '', '', '', '', ''], ['', '', '', '', '', ''],
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
        weighting_info_damage_carboneutrality = [[df.loc[i].tolist()[0], df.loc[i].tolist()[1], '', '', '', ''] for i in df.index]
        weighting_info_damage_carboneutrality[10] = ['Ionizing radiations, human health', '1.00E+00', '', '', '', '']
        weighting_info_damage_carboneutrality[13] = ['Photochemical ozone formation, human health', '1.00E+00', '', '', '', '']
        weighting_info_damage_carboneutrality.insert(22, ['Fisheries impact', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality[27] = ['Ionizing radiations, ecosystem quality', '1.00E+00', '', '', '', '']
        weighting_info_damage_carboneutrality.insert(32, ['Marine ecotoxicity, long term', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(33, ['Marine ecotoxicity, short term', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(35, ['Photochemical ozone formation, ecosystem quality', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(36, ['Plastics physical effects on biota', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(38, ['Terrestrial ecotoxicity, long term', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(39, ['Terrestrial ecotoxicity, short term', '1.00E+00', '', '', '', ''])

        weighting_info_combined_carboneutrality = weighting_info_damage_carboneutrality.copy()
        weighting_info_combined_carboneutrality[11] = ['Ozone layer depletion (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[12] = ['Particulate matter formation (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[23] = ['Freshwater acidification (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[26] = ['Freshwater eutrophication (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[28] = ['Land occupation, biodiversity (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[29] = ['Land transformation, biodiversity (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[34] = ['Marine eutrophication (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[36] = ['Plastics physical effects on biota (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[37] = ['Terrestrial acidification (damage)', '1.00E+00', '', '', '', '']

        weighting_info_damage = [[df.loc[i].tolist()[0], df.loc[i].tolist()[1], '', '', '', ''] for i in df.index]
        weighting_info_damage[4] = ['Climate change, HH, LT, fossil', '1.00E+00', '', '', '', '']
        weighting_info_damage[5] = ['Climate change, HH, ST, fossil', '1.00E+00', '', '', '', '']
        weighting_info_damage[10] = ['Ionizing radiations, human health', '1.00E+00', '', '', '', '']
        weighting_info_damage[13] = ['Photochemical ozone formation, human health', '1.00E+00', '', '', '', '']
        weighting_info_damage[20] = ['Climate change, EQ, LT, fossil', '1.00E+00', '', '', '', '']
        weighting_info_damage[21] = ['Climate change, EQ, ST, fossil', '1.00E+00', '', '', '', '']
        weighting_info_damage.insert(22, ['Fisheries impact', '1.00E+00', '', '', '', ''])
        weighting_info_damage[27] = ['Ionizing radiations, ecosystem quality', '1.00E+00', '', '', '', '']
        weighting_info_damage.insert(32, ['Marine ecotoxicity, long term', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(33, ['Marine ecotoxicity, short term', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(35, ['Photochemical ozone formation, ecosystem quality', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(36, ['Plastics physical effects on biota', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(38, ['Terrestrial ecotoxicity, long term', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(39, ['Terrestrial ecotoxicity, short term', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(6, ['Climate change, HH, LT, biogenic', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(7, ['Climate change, HH, ST, biogenic', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(8, ['Climate change, HH, LT, CO2 uptake', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(9, ['Climate change, HH, ST, CO2 uptake', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(10, ['Climate change, HH, LT, land transformation', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(11, ['Climate change, HH, ST, land transformation', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(28, ['Climate change, EQ, LT, biogenic', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(29, ['Climate change, EQ, ST, biogenic', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(30, ['Climate change, EQ, LT, CO2 uptake', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(31, ['Climate change, EQ, ST, CO2 uptake', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(32, ['Climate change, EQ, LT, land transformation', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(33, ['Climate change, EQ, ST, land transformation', '1.00E+00', '', '', '', ''])

        weighting_info_combined = weighting_info_damage.copy()
        weighting_info_combined[17] = ['Ozone layer depletion (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[18] = ['Particulate matter formation (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[35] = ['Freshwater acidification (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[38] = ['Freshwater eutrophication (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[40] = ['Land occupation, biodiversity (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[41] = ['Land transformation, biodiversity (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[46] = ['Marine eutrophication (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[48] = ['Plastics physical effects on biota (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[49] = ['Terrestrial acidification (damage)', '1.00E+00', '', '', '', '']

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

        d_ic_unit = self.iw_sp_carbon_neutrality.loc[self.iw_sp_carbon_neutrality['MP or Damage'] == 'Midpoint',
                              ['Impact category', 'CF unit']].drop_duplicates().set_index('Impact category').iloc[:,
                    0].to_dict()
        midpoint_values_carboneutrality = []
        for j in d_ic_unit.keys():
            midpoint_values_carboneutrality.append(['', '', '', '', '', ''])
            midpoint_values_carboneutrality.append(['Impact category', '', '', '', '', ''])
            midpoint_values_carboneutrality.append([j, d_ic_unit[j], '', '', '', ''])
            midpoint_values_carboneutrality.append(['', '', '', '', '', ''])
            midpoint_values_carboneutrality.append(['Substances', '', '', '', '', ''])
            df = self.iw_sp_carbon_neutrality[self.iw_sp_carbon_neutrality['Impact category'] == j]
            df = df[df['CF unit'] == d_ic_unit[j]]
            df = df[['Compartment', 'Sub-compartment', 'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit']]
            for i in df.index:
                if type(df.loc[i, 'CAS number']) == float:
                    df.loc[i, 'CAS number'] = ''
                midpoint_values_carboneutrality.append(df.loc[i].tolist())

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

        d_ic_unit = self.iw_sp_carbon_neutrality.loc[self.iw_sp_carbon_neutrality['MP or Damage'] == 'Damage',
                              ['Impact category', 'CF unit']].drop_duplicates().set_index('Impact category').iloc[:,
                    0].to_dict()
        damage_values_carboneutrality = []
        for j in d_ic_unit.keys():
            damage_values_carboneutrality.append(['', '', '', '', '', ''])
            damage_values_carboneutrality.append(['Impact category', '', '', '', '', ''])
            damage_values_carboneutrality.append([j, d_ic_unit[j], '', '', '', ''])
            damage_values_carboneutrality.append(['', '', '', '', '', ''])
            damage_values_carboneutrality.append(['Substances', '', '', '', '', ''])
            df = self.iw_sp_carbon_neutrality[self.iw_sp_carbon_neutrality['Impact category'] == j]
            df = df[df['CF unit'] == d_ic_unit[j]]
            df = df[['Compartment', 'Sub-compartment', 'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit']]
            for i in df.index:
                if type(df.loc[i, 'CAS number']) == float:
                    df.loc[i, 'CAS number'] = ''
                damage_values_carboneutrality.append(df.loc[i].tolist())

        # extracting combined CFs
        ic_unit = self.iw_sp.loc[:, ['Impact category', 'CF unit']].drop_duplicates()
        same_names = ['Freshwater acidification','Freshwater eutrophication','Land occupation, biodiversity',
                      'Land transformation, biodiversity','Marine eutrophication','Ozone layer depletion',
                      'Particulate matter formation','Terrestrial acidification',
                      'Plastics physical effects on biota']
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

        ic_unit = self.iw_sp_carbon_neutrality.loc[:, ['Impact category', 'CF unit']].drop_duplicates()
        same_names = ['Freshwater acidification','Freshwater eutrophication','Land occupation, biodiversity',
                      'Land transformation, biodiversity','Marine eutrophication','Ozone layer depletion',
                      'Particulate matter formation','Terrestrial acidification',
                      'Plastics physical effects on biota']
        combined_values_carboneutrality = []
        for j in ic_unit.index:
            combined_values_carboneutrality.append(['', '', '', '', '', ''])
            combined_values_carboneutrality.append(['Impact category', '', '', '', '', ''])
            if ic_unit.loc[j,'Impact category'] in same_names:
                if ic_unit.loc[j,'CF unit'] in ['DALY','PDF.m2.yr']:
                    combined_values_carboneutrality.append([ic_unit.loc[j,'Impact category']+' (damage)',
                                            ic_unit.loc[j,'CF unit'], '', '', '', ''])
                else:
                    combined_values_carboneutrality.append([ic_unit.loc[j,'Impact category']+' (midpoint)',
                                            ic_unit.loc[j,'CF unit'], '', '', '', ''])
            else:
                combined_values_carboneutrality.append([ic_unit.loc[j, 'Impact category'],
                                        ic_unit.loc[j, 'CF unit'], '', '', '', ''])
            combined_values_carboneutrality.append(['', '', '', '', '', ''])
            combined_values_carboneutrality.append(['Substances', '', '', '', '', ''])
            df = self.iw_sp_carbon_neutrality.loc[[i for i in self.iw_sp_carbon_neutrality.index if (
                    self.iw_sp_carbon_neutrality.loc[i,'Impact category'] == ic_unit.loc[j,'Impact category'] and
                    self.iw_sp_carbon_neutrality.loc[i, 'CF unit'] == ic_unit.loc[j, 'CF unit'])]]
            df = df[['Compartment', 'Sub-compartment', 'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit']]
            for i in df.index:
                if type(df.loc[i, 'CAS number']) == float:
                    df.loc[i, 'CAS number'] = ''
                combined_values_carboneutrality.append(df.loc[i].tolist())

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
                        'weighting_info_damage_carboneutrality': weighting_info_damage_carboneutrality,
                        'weighting_info_combined_carboneutrality': weighting_info_combined_carboneutrality,
                        'midpoint_values': midpoint_values, 'damage_values': damage_values,
                        'combined_values': combined_values, 'simplified_values': simplified_values,
                        'midpoint_method_metadata_carboneutrality': midpoint_method_metadata_carboneutrality,
                        'damage_method_metadata_carboneutrality': damage_method_metadata_carboneutrality,
                        'combined_method_metadata_carboneutrality': combined_method_metadata_carboneutrality,
                        'midpoint_values_carboneutrality': midpoint_values_carboneutrality,
                        'damage_values_carboneutrality': damage_values_carboneutrality,
                        'combined_values_carboneutrality': combined_values_carboneutrality}

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
                             "version": "2.1",
                             "modelType": "IMPACT_METHOD"}

        # -----------------------IW+ VERSION METADATA --------------------------
        category_names = {(i[0], i[1]): str(uuid.uuid4()) for i in
                          set(list(zip(self.olca_iw.loc[:, 'Impact category'], self.olca_iw.loc[:, 'CF unit'])))}
        category_names_damage = {k: v for k, v in category_names.items() if k[1] in ['DALY', 'PDF.m2.yr']}
        category_names_midpoint = {k: v for k, v in category_names.items() if k[1] not in ['DALY', 'PDF.m2.yr']}
        category_names_footprint = {(i[0], i[1]): str(uuid.uuid4()) for i in set(list(
            zip(self.simplified_version_olca.loc[:, 'Impact category'],
                self.simplified_version_olca.loc[:, 'CF unit'])))}
        category_names_carboneutrality = {(i[0], i[1]): str(uuid.uuid4()) for i in
                                          set(list(zip(self.olca_iw_carbon_neutrality.loc[:, 'Impact category'],
                                                       self.olca_iw_carbon_neutrality.loc[:, 'CF unit'])))}
        category_names_damage_carboneutrality = {k: v for k, v in category_names_carboneutrality.items() if
                                                 k[1] in ['DALY', 'PDF.m2.yr']}
        category_names_midpoint_carboneutrality = {k: v for k, v in category_names_carboneutrality.items() if
                                                   k[1] not in ['DALY', 'PDF.m2.yr']}
        # need to differentiate midpoint from endpoint with the names of the categories
        category_names_combined = {}
        for category in category_names:
            if category[1] in ['DALY', 'PDF.m2.yr']:
                category_names_combined[(category[0] + ' (damage)', category[1])] = category_names[category]
            else:
                category_names_combined[(category[0] + ' (midpoint)', category[1])] = category_names[category]

        category_names_combined_carboneutrality = {}
        for category in category_names_carboneutrality:
            if category[1] in ['DALY', 'PDF.m2.yr']:
                category_names_combined_carboneutrality[(category[0] + ' (damage)', category[1])] = \
                category_names_carboneutrality[category]
            else:
                category_names_combined_carboneutrality[(category[0] + ' (midpoint)', category[1])] = \
                category_names_carboneutrality[category]

        id_iw_damage = str(uuid.uuid4())
        id_iw_midpoint = str(uuid.uuid4())
        id_iw_footprint = str(uuid.uuid4())
        id_iw_combined = str(uuid.uuid4())
        norm_weight_id = str(uuid.uuid4())
        norm_weight_id_carboneutrality = str(uuid.uuid4())
        id_iw_damage_carboneutrality = str(uuid.uuid4())
        id_iw_midpoint_carboneutrality = str(uuid.uuid4())
        id_iw_combined_carboneutrality = str(uuid.uuid4())

        metadata_iw_damage = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "ImpactMethod",
            "@id": id_iw_damage,
            "name": "IMPACT World+ Expert v" + self.version + ' (incl. CO2 uptake)',
            "lastChange": "2024-09-15T17:25:43.725-05:00",
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
            "name": "IMPACT World+ Midpoint v" + self.version + ' (incl. CO2 uptake)',
            "lastChange": "2024-09-15T17:25:43.725-05:00",
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
            "lastChange": "2024-09-15T17:25:43.725-05:00",
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
            "name": "IMPACT World+ Combined v" + self.version + ' (incl. CO2 uptake)',
            "lastChange": "2024-09-15T17:25:43.725-05:00",
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

        metadata_iw_damage_carboneutrality = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "ImpactMethod",
            "@id": id_iw_damage_carboneutrality,
            "name": "IMPACT World+ Expert v" + self.version,
            "lastChange": "2024-09-15T17:25:43.725-05:00",
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

        metadata_iw_midpoint_carboneutrality = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "ImpactMethod",
            "@id": id_iw_midpoint_carboneutrality,
            "name": "IMPACT World+ Midpoint v" + self.version,
            "lastChange": "2024-09-15T17:25:43.725-05:00",
            "category": {
                "@type": "Category",
                "@id": id_category,
                "name": "IMPACT World+",
                "categoryType": "ImpactMethod"},
            'impactCategories': []
        }

        metadata_iw_combined_carboneutrality = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "ImpactMethod",
            "@id": id_iw_combined_carboneutrality,
            "name": "IMPACT World+ Combined v" + self.version,
            "lastChange": "2024-09-15T17:25:43.725-05:00",
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

        for cat in category_names_damage_carboneutrality:
            metadata_iw_damage_carboneutrality['impactCategories'].append(
                {"@type": "ImpactCategory",
                 "@id": category_names_damage_carboneutrality[cat],
                 "name": cat[0],
                 "refUnit": cat[1]}
            )

        for cat in category_names_midpoint_carboneutrality:
            metadata_iw_midpoint_carboneutrality['impactCategories'].append(
                {"@type": "ImpactCategory",
                 "@id": category_names_midpoint_carboneutrality[cat],
                 "name": cat[0],
                 "refUnit": cat[1]}
            )

        for cat in category_names_combined_carboneutrality:
            metadata_iw_combined_carboneutrality['impactCategories'].append(
                {"@type": "ImpactCategory",
                 "@id": category_names_combined_carboneutrality[cat],
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
            'MJ': '93a60a57-a3c8-11da-a746-0800200c9a66',
            'kg*a': '59f191d6-5dd3-4553-af88-1a32accfe308'
        }
        flow_properties = {
            'm3': '93a60a56-a3c8-22da-a746-0800200c9a66',
            'm2*a': '93a60a56-a3c8-21da-a746-0800200c9a66',
            'm2': '93a60a56-a3c8-19da-a746-0800200c9a66',
            'kBq': '93a60a56-a3c8-17da-a746-0800200c9a66',
            'kg': '93a60a56-a3c8-11da-a746-0800200b9a66',
            'MJ': 'f6811440-ee37-11de-8a39-0800200c9a66',
            'kg*a': '4e10f566-0358-489a-8e3a-d687b66c50e6'
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
                "version": '2.1',
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.olca_iw.loc[self.olca_iw.loc[:, 'Impact category'] == cat[0]].loc[
                self.olca_iw.loc[:, 'CF unit'] == cat[1]].copy()

            for i in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[i, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": dff.loc[i, 'flow_id'],
                        "name": dff.loc[i, 'Elem flow name'],
                        "category": "Elementary flows/" + dff.loc[i, 'Compartment'] + '/' + dff.loc[
                            i, 'Sub-compartment'],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[i, 'Elem flow unit']],
                        "name": dff.loc[i, 'Elem flow unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[i, 'Elem flow unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[i, 'Elem flow unit']
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
                "version": '2.1',
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.olca_iw.loc[self.olca_iw.loc[:, 'Impact category'] == cat[0]].loc[
                self.olca_iw.loc[:, 'CF unit'] == cat[1]].copy()

            for i in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[i, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": dff.loc[i, 'flow_id'],
                        "name": dff.loc[i, 'Elem flow name'],
                        "category": "Elementary flows/" + dff.loc[i, 'Compartment'] + '/' + dff.loc[
                            i, 'Sub-compartment'],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[i, 'Elem flow unit']],
                        "name": dff.loc[i, 'Elem flow unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[i, 'Elem flow unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[i, 'Elem flow unit']
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
                "version": '2.1',
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.simplified_version_olca.loc[self.simplified_version_olca.loc[:, 'Impact category'] == cat[0]].loc[
                self.simplified_version_olca.loc[:, 'CF unit'] == cat[1]].copy()

            for i in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[i, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": dff.loc[i, 'flow_id'],
                        "name": dff.loc[i, 'Elem flow name'],
                        "category": "Elementary flows/" + dff.loc[i, 'Compartment'] + '/' + dff.loc[
                            i, 'Sub-compartment'],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[i, 'Elem flow unit']],
                        "name": dff.loc[i, 'Elem flow unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[i, 'Elem flow unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[i, 'Elem flow unit']
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
                "version": '2.1',
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.olca_iw.loc[self.olca_iw.loc[:, 'Impact category'] == cat[0].split(' (')[0]].loc[
                self.olca_iw.loc[:, 'CF unit'] == cat[1]].copy()

            for i in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[i, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": dff.loc[i, 'flow_id'],
                        "name": dff.loc[i, 'Elem flow name'],
                        "category": "Elementary flows/" + dff.loc[i, 'Compartment'] + '/' + dff.loc[
                            i, 'Sub-compartment'],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[i, 'Elem flow unit']],
                        "name": dff.loc[i, 'Elem flow unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[i, 'Elem flow unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    }
                })

            cf_dict_combined[cat] = cf_values

        cf_dict_damage_carboneutrality = {}
        for cat in category_names_damage_carboneutrality:
            cf_values = {
                "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
                "@type": "ImpactCategory",
                "@id": category_names_damage_carboneutrality[cat],
                "name": cat[0],
                "category": metadata_iw_damage_carboneutrality["name"],
                "version": '2.1',
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.olca_iw_carbon_neutrality.loc[self.olca_iw_carbon_neutrality.loc[:, 'Impact category'] == cat[0]].loc[
                self.olca_iw_carbon_neutrality.loc[:, 'CF unit'] == cat[1]].copy()

            for i in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[i, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": dff.loc[i, 'flow_id'],
                        "name": dff.loc[i, 'Elem flow name'],
                        "category": "Elementary flows/" + dff.loc[i, 'Compartment'] + '/' + dff.loc[
                            i, 'Sub-compartment'],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[i, 'Elem flow unit']],
                        "name": dff.loc[i, 'Elem flow unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[i, 'Elem flow unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    }
                })

            cf_dict_damage_carboneutrality[cat] = cf_values

        cf_dict_midpoint_carboneutrality = {}
        for cat in category_names_midpoint_carboneutrality:
            cf_values = {
                "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
                "@type": "ImpactCategory",
                "@id": category_names_midpoint_carboneutrality[cat],
                "name": cat[0],
                "category": metadata_iw_midpoint_carboneutrality["name"],
                "version": '2.1',
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.olca_iw_carbon_neutrality.loc[self.olca_iw_carbon_neutrality.loc[:, 'Impact category'] == cat[0]].loc[
                self.olca_iw_carbon_neutrality.loc[:, 'CF unit'] == cat[1]].copy()

            for i in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[i, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": dff.loc[i, 'flow_id'],
                        "name": dff.loc[i, 'Elem flow name'],
                        "category": "Elementary flows/" + dff.loc[i, 'Compartment'] + '/' + dff.loc[
                            i, 'Sub-compartment'],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[i, 'Elem flow unit']],
                        "name": dff.loc[i, 'Elem flow unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[i, 'Elem flow unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    }
                })

            cf_dict_midpoint_carboneutrality[cat] = cf_values

        cf_dict_combined_carboneutrality = {}
        for cat in category_names_combined_carboneutrality:
            cf_values = {
                "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
                "@type": "ImpactCategory",
                "@id": category_names_combined_carboneutrality[cat],
                "name": cat[0],
                "category": metadata_iw_combined_carboneutrality["name"],
                "version": '2.1',
                "refUnit": cat[1],
                "impactFactors": []
            }

            dff = self.olca_iw_carbon_neutrality.loc[self.olca_iw_carbon_neutrality.loc[:, 'Impact category'] == cat[0].split(' (')[0]].loc[
                self.olca_iw_carbon_neutrality.loc[:, 'CF unit'] == cat[1]].copy()

            for i in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[i, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id": dff.loc[i, 'flow_id'],
                        "name": dff.loc[i, 'Elem flow name'],
                        "category": "Elementary flows/" + dff.loc[i, 'Compartment'] + '/' + dff.loc[
                            i, 'Sub-compartment'],
                        "flowType": "ELEMENTARY_FLOW",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    },
                    "unit": {
                        "@type": "Unit",
                        "@id": unit_groups[dff.loc[i, 'Elem flow unit']],
                        "name": dff.loc[i, 'Elem flow unit']
                    },
                    "flowProperty": {
                        "@type": "FlowProperty",
                        "@id": flow_properties[dff.loc[i, 'Elem flow unit']],
                        "category": "Technical flow properties",
                        "refUnit": dff.loc[i, 'Elem flow unit']
                    }
                })

            cf_dict_combined_carboneutrality[cat] = cf_values

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
            "version": "0.00.000",
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

        norm_carboneutrality = {}
        for i in category_names_carboneutrality:
            if i[1] == 'DALY':
                norm_carboneutrality[(i[0], i[1], category_names_carboneutrality[i])] = {"normalisationFactor": 13.7,
                                                                                         "weightingFactor": 5401.459854}
            elif i[1] == 'PDF.m2.yr':
                norm_carboneutrality[(i[0], i[1], category_names_carboneutrality[i])] = {"normalisationFactor": 1.01E-4,
                                                                                         "weightingFactor": 1386.138614}

        normalization_carboneutrality = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "NwSet",
            "@id": norm_weight_id_carboneutrality,
            "name": "IMPACT World+ (Stepwise 2006 values)",
            "version": "0.00.000",
            "weightedScoreUnit": "EUR2003",
            "factors": []}

        for x in norm_carboneutrality:
            normalization_carboneutrality['factors'].append(
                {
                    "@type": "NwFactor",
                    "impactCategory": {
                        "@type": "ImpactCategory",
                        "@id": x[2],
                        "name": x[0],
                        "refUnit": x[1]
                    },
                    "normalisationFactor": norm_carboneutrality[x]['normalisationFactor'],
                    "weightingFactor": norm_carboneutrality[x]['weightingFactor']
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
                          'normalization': normalization,
                          'metadata_iw_damage_carboneutrality': metadata_iw_damage_carboneutrality,
                          'metadata_iw_midpoint_carboneutrality': metadata_iw_midpoint_carboneutrality,
                          'metadata_iw_combined_carboneutrality': metadata_iw_combined_carboneutrality,
                          'cf_dict_damage_carboneutrality': cf_dict_damage_carboneutrality,
                          'cf_dict_midpoint_carboneutrality': cf_dict_midpoint_carboneutrality,
                          'cf_dict_combined_carboneutrality': cf_dict_combined_carboneutrality,
                          'normalization_carboneutrality': normalization_carboneutrality}

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
        self.master_db.to_excel(path + '/Dev/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_dev.xlsx')
        self.master_db_carbon_neutrality.to_excel(path + '/Dev/impact_world_plus_' + self.version + '_dev.xlsx')

        # ecoinvent versions in Excel format
        self.ei38_iw.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_expert_version_ecoinvent_v38.xlsx')
        self.ei39_iw.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_expert_version_ecoinvent_v39.xlsx')
        self.ei310_iw.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_expert_version_ecoinvent_v310.xlsx')
        self.ei311_iw.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_expert_version_ecoinvent_v311.xlsx')
        self.ei38_iw_carbon_neutrality.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v38.xlsx')
        self.ei39_iw_carbon_neutrality.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v39.xlsx')
        self.ei310_iw_carbon_neutrality.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v310.xlsx')
        self.ei311_iw_carbon_neutrality.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v311.xlsx')

        # ecoinvent version in DataFrame format
        self.simplified_version_ei38.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v38.xlsx')
        self.simplified_version_ei39.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v39.xlsx')
        self.simplified_version_ei310.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v310.xlsx')
        self.simplified_version_ei311.to_excel(path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v311.xlsx')

        # exiobase version in DataFrame format
        self.exio_iw.to_excel(path + '/exiobase/impact_world_plus_' + self.version + '_expert_version_exiobase.xlsx')

        # brightway2 versions in bw2package format
        for project in self.bw2_projects:
            bw2.projects.set_current(project)
            if '3.8' in project:
                IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                                                                        self.version in ic[0] and "for ecoinvent" in ic[0] and
                                                                        ' (incl. CO2 uptake)' in ic[0])]
                bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                     ' (incl. CO2 uptake)_brightway2_expert_version_ei38',
                                                     folder=path+'/bw2/')
                IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                                                                        self.version in ic[0] and "for ecoinvent" in ic[0] and
                                                                        ' (incl. CO2 uptake)' not in ic[0])]
                bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                     '_brightway2_expert_version_ei38', folder=path+'/bw2/')
                IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' in ic[0] and
                                                                        self.version in ic[0] and "for ecoinvent" in ic[0])]
                bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                     '_brightway2_footprint_version_ei38', folder=path+'/bw2/')
            elif '3.9' in project:
                IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                                                                        self.version in ic[0] and "for ecoinvent" in ic[0] and
                                                                        ' (incl. CO2 uptake)' in ic[0])]
                bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                     ' (incl. CO2 uptake)_brightway2_expert_version_ei39',
                                                     folder=path+'/bw2/')
                IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                                                                        self.version in ic[0] and "for ecoinvent" in ic[0] and
                                                                        ' (incl. CO2 uptake)' not in ic[0])]
                bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                     '_brightway2_expert_version_ei39', folder=path+'/bw2/')
                IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' in ic[0] and
                                                                        self.version in ic[0] and "for ecoinvent" in ic[0])]
                bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                     '_brightway2_footprint_version_ei39', folder=path+'/bw2/')
            elif '3.10' in project:
                IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                                                                        self.version in ic[0] and "for ecoinvent" in ic[0] and
                                                                        ' (incl. CO2 uptake)' in ic[0])]
                bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                     ' (incl. CO2 uptake)_brightway2_expert_version_ei310',
                                                     folder=path+'/bw2/')
                IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                                                                        self.version in ic[0] and "for ecoinvent" in ic[0] and
                                                                        ' (incl. CO2 uptake)' not in ic[0])]
                bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                     '_brightway2_expert_version_ei310', folder=path+'/bw2/')
                IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Footprint' in ic[0] and
                                                                        self.version in ic[0] and "for ecoinvent" in ic[0])]
                bw2io.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                     '_brightway2_footprint_version_ei310', folder=path+'/bw2/')

        # SimaPro version in csv format
        with open(path+'/SimaPro/impact_world_plus_'+self.version+' (incl. CO2 uptake)_midpoint_version_simapro.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['midpoint_method_metadata'] +
                self.sp_data['midpoint_values'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/impact_world_plus_'+self.version+' (incl. CO2 uptake)_expert_version_simapro.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['damage_method_metadata'] +
                self.sp_data['damage_values'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_damage'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/impact_world_plus_'+self.version+' (incl. CO2 uptake)_simapro.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['combined_method_metadata'] +
                self.sp_data['combined_values'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_combined'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/impact_world_plus_'+self.version+'_footprint_version_simapro.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['simplified_method_metadata'] +
                self.sp_data['simplified_values'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/impact_world_plus_'+self.version+'_midpoint_version_simapro.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['midpoint_method_metadata_carboneutrality'] +
                self.sp_data['midpoint_values_carboneutrality'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/impact_world_plus_'+self.version+'_expert_version_simapro.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['damage_method_metadata_carboneutrality'] +
                self.sp_data['damage_values_carboneutrality'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_damage_carboneutrality'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/impact_world_plus_'+self.version+'_simapro.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['combined_method_metadata_carboneutrality'] +
                self.sp_data['combined_values_carboneutrality'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_combined_carboneutrality'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])

        # create the openLCA version expert version (zip file)
        if not os.path.exists(path + '/openLCA/expert_version (incl. CO2 uptake)/'):
            os.makedirs(path + '/openLCA/expert_version (incl. CO2 uptake)/')
        if os.path.exists(path + '/openLCA/expert_version (incl. CO2 uptake)/impact_world_plus'+self.version+
                          ' (incl. CO2 uptake)_openLCA.zip'):
            os.remove(path + '/openLCA/expert_version (incl. CO2 uptake)/impact_world_plus'+self.version+
                      ' (incl. CO2 uptake)_openLCA.zip')
        if os.path.exists(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders'):
            shutil.rmtree(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders')
        zipObj = zipfile.ZipFile(path + '/openLCA/expert_version (incl. CO2 uptake)/impact_world_plus_'+self.version+
                                 ' (incl. CO2 uptake)_openLCA.zip', 'w')

        if not os.path.exists(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/categories/'):
            os.makedirs(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/categories/')
        with open(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/categories/' +
                  self.olca_data['category_metadata']['@id']  + '.json', 'w') as f:
            json.dump(self.olca_data['category_metadata'], f)
        zipObj.write(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/categories/' +
                     self.olca_data['category_metadata']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/' +
                  self.olca_data['metadata_iw_damage']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw_damage'], f)
        zipObj.write(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/' +
                     self.olca_data['metadata_iw_damage']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict_damage'].keys():
            with open(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/' +
                      self.olca_data['cf_dict_damage'][cat]['@id'] + '.json', 'w') as f:
                json.dump(self.olca_data['cf_dict_damage'][cat], f)
            zipObj.write(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/' +
                         self.olca_data['cf_dict_damage'][cat]['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/nw_sets/'):
            os.makedirs(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/nw_sets/')
        with open(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/nw_sets/' +
                  self.olca_data['normalization']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['normalization'], f)
        zipObj.write(path + '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/nw_sets/' +
                     self.olca_data['normalization']['@id'] + '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/expert_version (incl. CO2 uptake)/impact_world_plus_'+self.version+
                            ' (incl. CO2 uptake)_expert_version_openLCA', 'zip', path +
                            '/openLCA/expert_version (incl. CO2 uptake)/oLCA_folders/')

        # create the openLCA version midpoint version (zip file)
        if not os.path.exists(path + '/openLCA/midpoint_version (incl. CO2 uptake)/'):
            os.makedirs(path + '/openLCA/midpoint_version (incl. CO2 uptake)/')
        if os.path.exists(path + '/openLCA/midpoint_version (incl. CO2 uptake)/impact_world_plus'+self.version+
                          ' (incl. CO2 uptake)_openLCA.zip'):
            os.remove(path + '/openLCA/midpoint_version (incl. CO2 uptake)/impact_world_plus'+self.version+
                      ' (incl. CO2 uptake)_openLCA.zip')
        if os.path.exists(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders'):
            shutil.rmtree(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders')
        zipObj = zipfile.ZipFile(path + '/openLCA/midpoint_version (incl. CO2 uptake)/impact_world_plus_'+self.version+
                                 ' (incl. CO2 uptake)_openLCA.zip', 'w')

        if not os.path.exists(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/categories/'):
            os.makedirs(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/categories/')
        with open(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/categories/' +
                  self.olca_data['category_metadata']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['category_metadata'], f)
        zipObj.write(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/categories/' +
                     self.olca_data['category_metadata']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/' +
                  self.olca_data['metadata_iw_midpoint']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw_midpoint'], f)
        zipObj.write(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/' +
                     self.olca_data['metadata_iw_midpoint']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict_midpoint'].keys():
            with open(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/' +
                      self.olca_data['cf_dict_midpoint'][cat]['@id'] + '.json', 'w') as f:
                json.dump(self.olca_data['cf_dict_midpoint'][cat], f)
            zipObj.write(path + '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/' +
                         self.olca_data['cf_dict_midpoint'][cat]['@id'] + '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/midpoint_version (incl. CO2 uptake)/impact_world_plus_'+self.version+
                            ' (incl. CO2 uptake)_midpoint_version_openLCA', 'zip', path +
                            '/openLCA/midpoint_version (incl. CO2 uptake)/oLCA_folders/')

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
        if not os.path.exists(path + '/openLCA/combined_version (incl. CO2 uptake)/'):
            os.makedirs(path + '/openLCA/combined_version (incl. CO2 uptake)/')
        if os.path.exists(path + '/openLCA/combined_version (incl. CO2 uptake)/impact_world_plus'+self.version+
                          ' (incl. CO2 uptake)_openLCA.zip'):
            os.remove(path + '/openLCA/combined_version (incl. CO2 uptake)/impact_world_plus'+self.version+
                      ' (incl. CO2 uptake)_openLCA.zip')
        if os.path.exists(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders'):
            shutil.rmtree(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders')
        zipObj = zipfile.ZipFile(path + '/openLCA/combined_version (incl. CO2 uptake)/impact_world_plus_'+self.version+
                                 ' (incl. CO2 uptake)_openLCA.zip', 'w')

        if not os.path.exists(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/categories/'):
            os.makedirs(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/categories/')
        with open(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/categories/' +
                  self.olca_data['category_metadata']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['category_metadata'], f)
        zipObj.write(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/categories/' +
                     self.olca_data['category_metadata']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/' +
                  self.olca_data['metadata_iw_combined']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw_combined'], f)
        zipObj.write(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/lcia_methods/' +
                     self.olca_data['metadata_iw_combined']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict_combined'].keys():
            with open(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/' +
                      self.olca_data['cf_dict_combined'][cat]['@id'] + '.json', 'w') as f:
                json.dump(self.olca_data['cf_dict_combined'][cat], f)
            zipObj.write(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/lcia_categories/' +
                         self.olca_data['cf_dict_combined'][cat]['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/nw_sets/'):
            os.makedirs(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/nw_sets/')
        with open(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/nw_sets/' +
                  self.olca_data['normalization']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['normalization'], f)
        zipObj.write(path + '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/nw_sets/' +
                     self.olca_data['normalization']['@id'] + '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/combined_version (incl. CO2 uptake)/impact_world_plus_' + self.version +
                            ' (incl. CO2 uptake)_combined_version_openLCA', 'zip', path +
                            '/openLCA/combined_version (incl. CO2 uptake)/oLCA_folders/')

        # create the openLCA version carboneutrality expert version (zip file)
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
        zipObj.write(path + '/openLCA/expert_version/oLCA_folders/categories/' +
                     self.olca_data['category_metadata']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/expert_version/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/expert_version/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/expert_version/oLCA_folders/lcia_methods/' +
                  self.olca_data['metadata_iw_damage_carboneutrality']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw_damage_carboneutrality'], f)
        zipObj.write(path + '/openLCA/expert_version/oLCA_folders/lcia_methods/' +
                     self.olca_data['metadata_iw_damage_carboneutrality']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/expert_version/oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/expert_version/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict_damage_carboneutrality'].keys():
            with open(path + '/openLCA/expert_version/oLCA_folders/lcia_categories/' +
                      self.olca_data['cf_dict_damage_carboneutrality'][cat]['@id'] + '.json', 'w') as f:
                json.dump(self.olca_data['cf_dict_damage_carboneutrality'][cat], f)
            zipObj.write(path + '/openLCA/expert_version/oLCA_folders/lcia_categories/' +
                         self.olca_data['cf_dict_damage_carboneutrality'][cat]['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/expert_version/oLCA_folders/nw_sets/'):
            os.makedirs(path + '/openLCA/expert_version/oLCA_folders/nw_sets/')
        with open(path + '/openLCA/expert_version/oLCA_folders/nw_sets/' +
                  self.olca_data['normalization_carboneutrality']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['normalization_carboneutrality'], f)
        zipObj.write(path + '/openLCA/expert_version/oLCA_folders/nw_sets/' +
                     self.olca_data['normalization_carboneutrality']['@id'] + '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/expert_version/impact_world_plus_'+self.version+
                            '_expert_version_openLCA', 'zip', path +
                            '/openLCA/expert_version/oLCA_folders/')

        # create the openLCA version carboneutrality midpoint version (zip file)
        if not os.path.exists(path + '/openLCA/midpoint_version/'):
            os.makedirs(path + '/openLCA/midpoint_version/')
        if os.path.exists(path + '/openLCA/midpoint_version/impact_world_plus'+self.version+'_openLCA.zip'):
            os.remove(path + '/openLCA/midpoint_version/impact_world_plus'+self.version+'_openLCA.zip')
        if os.path.exists(path + '/openLCA/midpoint_version/oLCA_folders'):
            shutil.rmtree(path + '/openLCA/midpoint_version/oLCA_folders')
        zipObj = zipfile.ZipFile(path + '/openLCA/midpoint_version/impact_world_plus_'+self.version+'_openLCA.zip', 'w')

        if not os.path.exists(path + '/openLCA/midpoint_version/oLCA_folders/categories/'):
            os.makedirs(path + '/openLCA/midpoint_version/oLCA_folders/categories/')
        with open(path + '/openLCA/midpoint_version/oLCA_folders/categories/' +
                  self.olca_data['category_metadata']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['category_metadata'], f)
        zipObj.write(path + '/openLCA/midpoint_version/oLCA_folders/categories/' +
                     self.olca_data['category_metadata']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/midpoint_version/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/midpoint_version/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/midpoint_version/oLCA_folders/lcia_methods/' +
                  self.olca_data['metadata_iw_midpoint_carboneutrality']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw_midpoint_carboneutrality'], f)
        zipObj.write(path + '/openLCA/midpoint_version/oLCA_folders/lcia_methods/' +
                     self.olca_data['metadata_iw_midpoint_carboneutrality']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/midpoint_version /oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/midpoint_version/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict_midpoint_carboneutrality'].keys():
            with open(path + '/openLCA/midpoint_version/oLCA_folders/lcia_categories/' +
                      self.olca_data['cf_dict_midpoint_carboneutrality'][cat]['@id'] + '.json', 'w') as f:
                json.dump(self.olca_data['cf_dict_midpoint_carboneutrality'][cat], f)
            zipObj.write(path + '/openLCA/midpoint_version/oLCA_folders/lcia_categories/' +
                         self.olca_data['cf_dict_midpoint_carboneutrality'][cat]['@id'] + '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/midpoint_version/impact_world_plus_'+self.version+
                            '_midpoint_version_openLCA', 'zip', path +
                            '/openLCA/midpoint_version/oLCA_folders/')

        # create the openLCA version carboneutrality combined version (zip file)
        if not os.path.exists(path + '/openLCA/combined_version/'):
            os.makedirs(path + '/openLCA/combined_version/')
        if os.path.exists(path + '/openLCA/combined_version/impact_world_plus'+self.version+'_openLCA.zip'):
            os.remove(path + '/openLCA/combined_version/impact_world_plus'+self.version+'_openLCA.zip')
        if os.path.exists(path + '/openLCA/combined_version/oLCA_folders'):
            shutil.rmtree(path + '/openLCA/combined_version/oLCA_folders')
        zipObj = zipfile.ZipFile(path + '/openLCA/combined_version/impact_world_plus_'+self.version+'_openLCA.zip', 'w')

        if not os.path.exists(path + '/openLCA/combined_version/oLCA_folders/categories/'):
            os.makedirs(path + '/openLCA/combined_version/oLCA_folders/categories/')
        with open(path + '/openLCA/combined_version/oLCA_folders/categories/' +
                  self.olca_data['category_metadata']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['category_metadata'], f)
        zipObj.write(path + '/openLCA/combined_version/oLCA_folders/categories/' +
                     self.olca_data['category_metadata']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/combined_version/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/combined_version/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/combined_version/oLCA_folders/lcia_methods/' +
                  self.olca_data['metadata_iw_combined_carboneutrality']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw_combined_carboneutrality'], f)
        zipObj.write(path + '/openLCA/combined_version/oLCA_folders/lcia_methods/' +
                     self.olca_data['metadata_iw_combined_carboneutrality']['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/combined_version/oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/combined_version/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict_combined_carboneutrality'].keys():
            with open(path + '/openLCA/combined_version/oLCA_folders/lcia_categories/' +
                      self.olca_data['cf_dict_combined_carboneutrality'][cat]['@id'] + '.json', 'w') as f:
                json.dump(self.olca_data['cf_dict_combined_carboneutrality'][cat], f)
            zipObj.write(path + '/openLCA/combined_version/oLCA_folders/lcia_categories/' +
                         self.olca_data['cf_dict_combined_carboneutrality'][cat]['@id'] + '.json')

        if not os.path.exists(path + '/openLCA/combined_version /oLCA_folders/nw_sets/'):
            os.makedirs(path + '/openLCA/combined_version/oLCA_folders/nw_sets/')
        with open(path + '/openLCA/combined_version/oLCA_folders/nw_sets/' +
                  self.olca_data['normalization_carboneutrality']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['normalization_carboneutrality'], f)
        zipObj.write(path + '/openLCA/combined_version/oLCA_folders/nw_sets/' +
                     self.olca_data['normalization_carboneutrality']['@id'] + '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/combined_version/impact_world_plus_' + self.version +
                            '_combined_version_openLCA', 'zip', path +
                            '/openLCA/combined_version/oLCA_folders/')

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
            - Marine ecotoxicity, short term
            - Marine ecotoxicity, long term
            - Terrestrial ecotoxicity, short term
            - Terrestrial ecotoxicity, long term
            - Human toxicity cancer
            - Human toxicity cancer, long term
            - Human toxicity cancer, short term
            - Human toxicity non cancer
            - Human toxicity non cancer, long term
            - Human toxicity non cancer, short term
            - Ionizing radiations, ecosystem quality
            - Ionizing radiations, human health
            - Ionizing radiations
            - Marine acidification, short term
            = Marine acidification, long term
            - Mineral resources use

        :return: updated master_db
        """

        self.logger.info("Loading ionizing radiations characterization factors...")
        ionizing = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - IonizingRadiations]', con=self.conn)

        self.logger.info("Loading marine acidification characterization factors...")
        mar_acid = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - MarineAcidification]', con=self.conn)

        self.logger.info("Loading fossil resources characterization factors...")
        fossils = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - FossilResources]', con=self.conn)

        self.logger.info("Loading mineral resources characterization factors...")
        minerals = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - MineralResources]', con=self.conn)

        elem_flow_list = pd.read_sql(sql='SELECT * FROM [SI - Mapping with elementary flows]', con=self.conn)

        self.logger.info("Loading toxicity characterization factors...")
        toxicity = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - HumanTox]', con=self.conn)
        toxicity = toxicity.merge(elem_flow_list.loc[:, ['Name IW+', 'CAS-Usetox2_FW']], left_on=['CAS number'],
                                  right_on=['CAS-Usetox2_FW'], how='left').drop_duplicates()
        toxicity = toxicity.drop(['Elem flow name', 'CAS number'], axis=1)
        toxicity = toxicity.rename(columns={'Name IW+': 'Elem flow name', 'CAS-Usetox2_FW': 'CAS number'})

        self.logger.info("Loading freshwater ecotoxicity characterization factors...")
        fw_ecotoxicity = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - FreshwaterEcotox]', con=self.conn)
        fw_ecotoxicity = fw_ecotoxicity.merge(elem_flow_list.loc[:, ['Name IW+', 'CAS-Usetox2_FW']], left_on=['CAS number'],
                                        right_on=['CAS-Usetox2_FW'], how='left').drop_duplicates()
        fw_ecotoxicity = fw_ecotoxicity.drop(['Elem flow name', 'CAS number'], axis=1)
        fw_ecotoxicity = fw_ecotoxicity.rename(columns={'Name IW+': 'Elem flow name', 'CAS-Usetox2_FW': 'CAS number'})

        self.logger.info("Loading marine ecotoxicity characterization factors...")
        mar_ecotoxicity = pd.read_sql('SELECT * FROM [CF - not regionalized - MarineEcotox]', self.conn)
        mar_ecotoxicity = mar_ecotoxicity.merge(elem_flow_list.loc[:, ['Name IW+', 'CAS-Usetox2_Mar_Terr']],
                                                left_on=['CAS number'],
                                                right_on=['CAS-Usetox2_Mar_Terr'], how='left').drop_duplicates()
        mar_ecotoxicity = mar_ecotoxicity.drop(['Elem flow name', 'CAS number'], axis=1)
        mar_ecotoxicity = mar_ecotoxicity.rename(
            columns={'Name IW+': 'Elem flow name', 'CAS-Usetox2_Mar_Terr': 'CAS number'})

        self.logger.info("Loading terrestrial ecotoxicity characterization factors...")
        terr_ecotoxicity = pd.read_sql('SELECT * FROM [CF - not regionalized - TerrestrialEcotox]', self.conn)
        terr_ecotoxicity = terr_ecotoxicity.merge(elem_flow_list.loc[:, ['Name IW+', 'CAS-Usetox2_Mar_Terr']],
                                                left_on=['CAS number'],
                                                right_on=['CAS-Usetox2_Mar_Terr'], how='left').drop_duplicates()
        terr_ecotoxicity = terr_ecotoxicity.drop(['Elem flow name', 'CAS number'], axis=1)
        terr_ecotoxicity = terr_ecotoxicity.rename(
            columns={'Name IW+': 'Elem flow name', 'CAS-Usetox2_Mar_Terr': 'CAS number'})

        self.master_db = pd.concat([ionizing, mar_acid, fossils, minerals, toxicity,
                                    fw_ecotoxicity, mar_ecotoxicity, terr_ecotoxicity])

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
        # add carbon monoxide, which is based on the (C) stoechiometric ratio between CO2 and CO (1.57)
        monoxide = data.loc[data.Name == 'Carbon dioxide'].copy()
        monoxide.Name = 'Carbon monoxide'
        monoxide.Formula = 'CO'
        monoxide.loc[:, 'Lifetime (yr)'] = 2 / 12  # 2 months
        for indicator in ['Radiative Efficiency (W/m2/ppb)', 'AGWP-20 (pW/m2/yr/kg)', 'GWP-20',
                          'AGWP-100 (pW/m2/yr/kg)', 'GWP-100', 'AGWP-500 (pW/m2/yr/kg)', 'GWP-500',
                          'AGTP-50 (pK/kg)', 'GTP-50', 'AGTP-100 (pK/kg)', 'GTP-100']:
            monoxide.loc[monoxide.index[0], indicator] = float(monoxide.loc[monoxide.index[0], indicator]) * 1.57
        data = clean_up_dataframe(pd.concat([data, monoxide]))

        mapping = pd.read_sql('SELECT * FROM "SI - Mapping with elementary flows"', con=self.conn)
        data = data.merge(mapping.loc[:, ['Name IW+', 'Name-ipcc', 'CAS IW+']].dropna(subset='Name-ipcc'),
                          left_on='Name', right_on='Name-ipcc', how='inner').drop(['Name', 'Name-ipcc'], axis=1)

        # ---------------------------- Climate change midpoint indicators ---------------------------------------------
        # Climate change, short term
        GWP_midpoint = data.loc[:, ['Name IW+', 'GWP-100', 'CAS IW+']]
        GWP_midpoint.columns = ['Elem flow name', 'CF value', 'CAS number']
        GWP_midpoint.loc[:, 'Impact category'] = 'Climate change, short term'
        GWP_midpoint.loc[:, 'CF unit'] = 'kg CO2 eq (short)'
        GWP_midpoint.loc[:, 'Compartment'] = 'Air'
        GWP_midpoint.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_midpoint.loc[:, 'Elem flow unit'] = 'kg'
        GWP_midpoint.loc[:, 'MP or Damage'] = 'Midpoint'
        GWP_midpoint.loc[:, 'Native geographical resolution scale'] = 'Global'

        # Climate change, long term
        GTP_midpoint = data.loc[:, ['Name IW+', 'GTP-100', 'CAS IW+']]
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

        AACO2 = data.loc[data.loc[:, 'Name IW+'] == 'Carbon dioxide, fossil', 'Radiative Efficiency (W/m2/kg)'].iloc[0]

        # AGTP cumulative 100 years
        data.loc[:, 'AGTP cumulative 100 years'] = 0

        for t in range(0, 100):
            data.loc[data[data.loc[:, 'Name IW+'] == 'Carbon dioxide, fossil'].index[
                         0], 'AGTP cumulative 100 years'] += AGTPCO2(t, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3,
                                                                     kPulseT, aT1, tauT1, aT2, tauT2, AACO2)

            data.loc[data[data.loc[:, 'Name IW+'] == 'Methane, fossil'].index[
                         0], 'AGTP cumulative 100 years'] += AGTPCH4Fossil_Final(
                t, data.loc[data.loc[:, 'Name IW+'] == 'Methane, fossil', 'Lifetime (yr)'].iloc[0], kPulseT, aT1, tauT1,
                aT2, tauT2,
                data.loc[data.loc[:, 'Name IW+'] == 'Methane, fossil', 'Radiative Efficiency (W/m2/kg)'].iloc[0], aC1,
                aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

            data.loc[data[data.loc[:, 'Name IW+'] == 'Methane, biogenic'].index[
                         0], 'AGTP cumulative 100 years'] += AGTPCH4NonFossil_Final(
                t, data.loc[data.loc[:, 'Name IW+'] == 'Methane, biogenic', 'Lifetime (yr)'].iloc[0], kPulseT, aT1,
                tauT1, aT2, tauT2,
                data.loc[data.loc[:, 'Name IW+'] == 'Methane, biogenic', 'Radiative Efficiency (W/m2/kg)'].iloc[0], aC1,
                aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

            for ghg in data.index:
                if data.loc[ghg, 'Name IW+'] not in ['Carbon dioxide, fossil', 'Methane, fossil', 'Methane, biogenic']:
                    if data.loc[ghg, 'Lifetime (yr)'] != 0:
                        data.loc[ghg, 'AGTP cumulative 100 years'] += AGTPNonCO2_Final(
                            t, data.loc[ghg, 'Lifetime (yr)'], kPulseT, aT1, tauT1, aT2, tauT2,
                            data.loc[ghg, 'Radiative Efficiency (W/m2/kg)'], aC1, aC2, aC3, aC4, tauC1, tauC2,
                            tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

        # AGTP cumulative 500 years
        data.loc[:, 'AGTP cumulative 500 years'] = 0

        for t in range(0, 500):
            data.loc[data[data.loc[:, 'Name IW+'] == 'Carbon dioxide, fossil'].index[
                         0], 'AGTP cumulative 500 years'] += AGTPCO2(t, aC1, aC2, aC3, aC4, tauC1, tauC2, tauC3,
                                                                     kPulseT, aT1, tauT1, aT2, tauT2, AACO2)

            data.loc[data[data.loc[:, 'Name IW+'] == 'Methane, fossil'].index[
                         0], 'AGTP cumulative 500 years'] += AGTPCH4Fossil_Final(
                t, data.loc[data.loc[:, 'Name IW+'] == 'Methane, fossil', 'Lifetime (yr)'].iloc[0], kPulseT, aT1, tauT1,
                aT2, tauT2,
                data.loc[data.loc[:, 'Name IW+'] == 'Methane, fossil', 'Radiative Efficiency (W/m2/kg)'].iloc[0], aC1,
                aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

            data.loc[data[data.loc[:, 'Name IW+'] == 'Methane, biogenic'].index[
                         0], 'AGTP cumulative 500 years'] += AGTPCH4NonFossil_Final(
                t, data.loc[data.loc[:, 'Name IW+'] == 'Methane, biogenic', 'Lifetime (yr)'].iloc[0], kPulseT, aT1,
                tauT1, aT2, tauT2,
                data.loc[data.loc[:, 'Name IW+'] == 'Methane, biogenic', 'Radiative Efficiency (W/m2/kg)'].iloc[0], aC1,
                aC2, aC3, aC4, tauC1, tauC2, tauC3, AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

            for ghg in data.index:
                if data.loc[ghg, 'Name IW+'] not in ['Carbon dioxide, fossil', 'Methane, fossil', 'Methane, biogenic']:
                    if data.loc[ghg, 'Lifetime (yr)'] != 0:
                        data.loc[ghg, 'AGTP cumulative 500 years'] += AGTPNonCO2_Final(
                            t, data.loc[ghg, 'Lifetime (yr)'], kPulseT, aT1, tauT1, aT2, tauT2,
                            data.loc[ghg, 'Radiative Efficiency (W/m2/kg)'],aC1, aC2, aC3, aC4, tauC1, tauC2,tauC3,
                            AACO2, gamma, aS1, aS2, aS3, tauS1, tauS2, tauS3)

        # get effect factors
        effect_factors = pd.read_sql('SELECT * FROM "SI - Climate change - effect factors"', con=self.conn).set_index('index')
        HH_effect_factor = effect_factors.loc['Total', 'DALY/°C/yr']
        EQ_effect_factor = effect_factors.loc['Total', 'PDF.m2.yr/K/yr']

        # Climate change, human health, short term
        GWP_damage_HH_short = data.loc[:, ['Name IW+', 'AGTP cumulative 100 years', 'CAS IW+']].copy()
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
        GWP_damage_HH_long = data.loc[:, ['Name IW+', 'AGTP cumulative 100 years', 'AGTP cumulative 500 years',
                                          'CAS IW+']].copy()
        GWP_damage_HH_long = GWP_damage_HH_long.rename(columns={'Name IW+': 'Elem flow name','CAS IW+':'CAS number'})
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
        GWP_damage_EQ_short = data.loc[:, ['Name IW+', 'AGTP cumulative 100 years', 'CAS IW+']]
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
        GWP_damage_EQ_long = data.loc[:, ['Name IW+', 'AGTP cumulative 100 years', 'AGTP cumulative 500 years',
                                          'CAS IW+']].copy()
        GWP_damage_EQ_long = GWP_damage_EQ_long.rename(columns={'Name IW+': 'Elem flow name','CAS IW+':'CAS number'})
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
        mapping = pd.read_sql('SELECT * FROM "SI - Mapping with elementary flows"', con=self.conn)
        data.loc[:, 'Elem flow name'] = [dict(zip(mapping.loc[:, 'Name-ODP'], mapping.loc[:, "Name IW+"]))[i] for i in
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
        mapping = pd.read_sql('SELECT * FROM "SI - Mapping with elementary flows"', con=self.conn)
        effect_factor = pd.read_sql(sql='SELECT * FROM [SI - Photochemical ozone formation - effect factors]',
                                    con=self.conn)
        data = data.merge(mapping, left_on=['Substance name'], right_on=['Name-photochem'])

        photochem_damage_hh = data.copy('deep')
        photochem_damage_hh = photochem_damage_hh.loc[:, ['CF HH (kg NOx-eq/kg)', 'Name IW+', 'CAS IW+']]
        photochem_damage_hh = photochem_damage_hh.rename(columns={'CF HH (kg NOx-eq/kg)': 'CF value',
                                                                  'Name IW+': 'Elem flow name',
                                                                  'CAS IW+': 'CAS number'})
        photochem_damage_hh.loc[:, 'CF value'] *= effect_factor.loc[0, 'HH (DALY/kg NOx-eq)']
        photochem_damage_hh.loc[:, 'Impact category'] = 'Photochemical ozone formation, human health'
        photochem_damage_hh.loc[:, 'CF unit'] = 'DALY'
        photochem_damage_hh.loc[:, 'Compartment'] = 'Air'
        photochem_damage_hh.loc[:, 'Sub-compartment'] = '(unspecified)'
        photochem_damage_hh.loc[:, 'Elem flow unit'] = 'kg'
        photochem_damage_hh.loc[:, 'MP or Damage'] = 'Damage'
        photochem_damage_hh.loc[:, 'Native geographical resolution scale'] = 'Not regionalized'

        photochem_damage_eq = data.copy('deep')
        photochem_damage_eq = photochem_damage_eq.loc[:, ['CF EQ (kg NOx-eq/kg)', 'Name IW+', 'CAS IW+']]
        photochem_damage_eq = photochem_damage_eq.rename(columns={'CF EQ (kg NOx-eq/kg)': 'CF value',
                                                                  'Name IW+': 'Elem flow name',
                                                                  'CAS IW+': 'CAS number'})
        photochem_damage_eq.loc[:, 'CF value'] *= (effect_factor.loc[0, 'EQ (species.yr/kg NOx-eq)'] /
                                                   effect_factor.loc[0, 'species density in ReCiPe (species/m2)'])
        photochem_damage_eq.loc[:, 'Impact category'] = 'Photochemical ozone formation, ecosystem quality'
        photochem_damage_eq.loc[:, 'CF unit'] = 'PDF.m2.yr'
        photochem_damage_eq.loc[:, 'Compartment'] = 'Air'
        photochem_damage_eq.loc[:, 'Sub-compartment'] = '(unspecified)'
        photochem_damage_eq.loc[:, 'Elem flow unit'] = 'kg'
        photochem_damage_eq.loc[:, 'MP or Damage'] = 'Damage'
        photochem_damage_eq.loc[:, 'Native geographical resolution scale'] = 'Not regionalized'

        # ReCiPe separates HH and EQ at midpoint. We only base our midpoint on  HH.
        photochem_midpoint = data.copy('deep')
        photochem_midpoint.loc[:, 'CF value'] = photochem_midpoint.loc[:, 'CF HH (kg NOx-eq/kg)']
        photochem_midpoint = photochem_midpoint.rename(columns={'Name IW+': 'Elem flow name',
                                                                'CAS IW+': 'CAS number'})
        photochem_midpoint = photochem_midpoint.loc[:, ['CF value', 'Elem flow name', 'CAS number']]
        photochem_midpoint.loc[:, 'Impact category'] = 'Photochemical ozone formation'
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
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NH3', 'CAS number'] = '7664-41-7'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NOx', 'CAS number'] = '11104-93-1'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'SO2', 'CAS number'] = '7446-09-5'
        cfs = clean_up_dataframe(cfs)
        cfs.loc[[i for i in cfs.index if cfs.loc[i, 'Region code'] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF', 'OCE', 'GLO']],
                'Native geographical resolution scale'] = 'Continent'
        cfs.loc[cfs.loc[:, 'Region code'] == 'GLO', 'Native geographical resolution scale'] = 'Global'

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - Stoechiometry]', self.conn)
        stoc = stoc.loc[stoc.loc[:, 'Impact category'] == 'Freshwater acidification']
        for i in stoc.index:
            proxy = stoc.loc[i, 'Proxy molecule']
            comp = stoc.loc[i, 'Compartment']
            df = cfs[cfs.loc[:, 'Elem flow'] == proxy].loc[cfs.loc[:, 'Compartment'] == comp].copy('deep')
            if not df.empty:
                df.loc[:, 'Elem flow'] = stoc.loc[i, 'Elem flow name']
                df.loc[:, 'CAS number'] = stoc.loc[i, 'CAS number']
                df.loc[:, 'CF value'] *= stoc.loc[i, 'Proxy ratio']

                cfs = clean_up_dataframe(pd.concat([cfs, df]))

        # ------------------------------------ FINAL FORMATTING ---------------------------------
        cfs = cfs.drop(cfs.loc[cfs.loc[:, 'Elem flow'].isin(['NH3', 'NOx', 'SO2'])].index)
        # concat elem flow and region code in the same name
        cfs.loc[:, 'Elem flow name'] = [', '.join(i) for i in list(zip(cfs.loc[:, 'Elem flow'], cfs.loc[:, 'Region code']))]
        cfs = cfs.drop(['Elem flow', 'Region code'], axis=1)
        cfs = clean_up_dataframe(cfs)

        # add value for RME (Middle East) based on RAS (Asia)
        df = cfs.loc[cfs.loc[:, 'Elem flow name'].str.contains(', RAS')].copy()
        df.loc[:, 'Elem flow name'] = [i.split(', RAS')[0] + ', RME' for i in df.loc[:, 'Elem flow name']]
        cfs = clean_up_dataframe(pd.concat([cfs, df]))
        # add RoW based on GLO
        df = cfs.loc[cfs.loc[:, 'Elem flow name'].str.contains(', GLO')].copy()
        df.loc[:, 'Elem flow name'] = [i.split(', GLO')[0] + ', RoW' for i in df.loc[:, 'Elem flow name']]
        df.loc[:, 'Native geographical resolution scale'] = 'Other region'
        cfs = clean_up_dataframe(pd.concat([cfs, df]))

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
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NH3', 'CAS number'] = '7664-41-7'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NOx', 'CAS number'] = '11104-93-1'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'SO2', 'CAS number'] = '7446-09-5'
        cfs = clean_up_dataframe(cfs)
        cfs.loc[[i for i in cfs.index if cfs.loc[i, 'Region code'] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF', 'OCE', 'GLO']],
                'Native geographical resolution scale'] = 'Continent'
        cfs.loc[cfs.loc[:, 'Region code'] == 'GLO', 'Native geographical resolution scale'] = 'Global'

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - Stoechiometry]', self.conn)
        stoc = stoc.loc[stoc.loc[:, 'Impact category'] == 'Terrestrial acidification']
        for i in stoc.index:
            proxy = stoc.loc[i, 'Proxy molecule']
            comp = stoc.loc[i, 'Compartment']
            df = cfs[cfs.loc[:, 'Elem flow'] == proxy].loc[cfs.loc[:, 'Compartment'] == comp].copy('deep')
            if not df.empty:
                df.loc[:, 'Elem flow'] = stoc.loc[i, 'Elem flow name']
                df.loc[:, 'CAS number'] = stoc.loc[i, 'CAS number']
                df.loc[:, 'CF value'] *= stoc.loc[i, 'Proxy ratio']

                cfs = clean_up_dataframe(pd.concat([cfs, df]))

        # ------------------------------------ FINAL FORMATTING ---------------------------------
        cfs = cfs.drop(cfs.loc[cfs.loc[:, 'Elem flow'].isin(['NH3', 'NOx', 'SO2'])].index)
        # concat elem flow and region code in the same name
        cfs.loc[:, 'Elem flow name'] = [', '.join(i) for i in list(zip(cfs.loc[:, 'Elem flow'], cfs.loc[:, 'Region code']))]
        cfs = cfs.drop(['Elem flow', 'Region code'], axis=1)
        cfs = clean_up_dataframe(cfs)

        # add value for RME (Middle East) based on RAS (Asia)
        df = cfs.loc[cfs.loc[:, 'Elem flow name'].str.contains(', RAS')].copy()
        df.loc[:, 'Elem flow name'] = [i.split(', RAS')[0] + ', RME' for i in df.loc[:, 'Elem flow name']]
        cfs = clean_up_dataframe(pd.concat([cfs, df]))
        # add RoW based on GLO
        df = cfs.loc[cfs.loc[:, 'Elem flow name'].str.contains(', GLO')].copy()
        df.loc[:, 'Elem flow name'] = [i.split(', GLO')[0] + ', RoW' for i in df.loc[:, 'Elem flow name']]
        df.loc[:, 'Native geographical resolution scale'] = 'Other region'
        cfs = clean_up_dataframe(pd.concat([cfs, df]))

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
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NH3', 'CAS number'] = '7664-41-7'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'NOx', 'CAS number'] = '11104-93-1'
        cfs.loc[cfs.loc[:, 'Elem flow'] == 'HNO3', 'CAS number'] = '7697-37-2'
        cfs = clean_up_dataframe(cfs)
        cfs.loc[[i for i in cfs.index if cfs.loc[i, 'Region code'] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF', 'OCE', 'GLO']],
                'Native geographical resolution scale'] = 'Continent'
        cfs.loc[cfs.loc[:, 'Region code'] == 'GLO', 'Native geographical resolution scale'] = 'Global'

        # add non-regionalized flows (water emissions)
        cfs = clean_up_dataframe(pd.concat(
            [cfs, pd.read_sql('SELECT * FROM [CF - not regionalized - MarEutro]', self.conn)]))

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - Stoechiometry]', self.conn)
        stoc = stoc.loc[stoc.loc[:, 'Impact category'] == 'Marine eutrophication']
        for i in stoc.index:
            proxy = stoc.loc[i, 'Proxy molecule']
            comp = stoc.loc[i, 'Compartment']
            df = cfs[cfs.loc[:, 'Elem flow'] == proxy].loc[cfs.loc[:, 'Compartment'] == comp].copy('deep')
            if not df.empty:
                df.loc[:, 'Elem flow'] = stoc.loc[i, 'Elem flow name']
                df.loc[:, 'CAS number'] = stoc.loc[i, 'CAS number']
                df.loc[:, 'CF value'] *= stoc.loc[i, 'Proxy ratio']

                cfs = clean_up_dataframe(pd.concat([cfs, df]))

        # ------------------------------------ FINAL FORMATTING ---------------------------------
        cfs = cfs.drop(cfs.loc[cfs.loc[:, 'Elem flow'].isin(['NH3', 'NOx', 'SO2'])].index)
        # concat elem flow and region code in the same name
        cfs.loc[:, 'Elem flow name'] = [', '.join(i) if type(i[1]) == str else i[0]
                                        for i in list(zip(cfs.loc[:, 'Elem flow'], cfs.loc[:, 'Region code']))]
        cfs = cfs.drop(['Elem flow', 'Region code'], axis=1)
        cfs = clean_up_dataframe(cfs)

        # add value for RME (Middle East) based on RAS (Asia)
        df = cfs.loc[cfs.loc[:, 'Elem flow name'].str.contains(', RAS')].copy()
        df.loc[:, 'Elem flow name'] = [i.split(', RAS')[0] + ', RME' for i in df.loc[:, 'Elem flow name']]
        cfs = clean_up_dataframe(pd.concat([cfs, df]))
        # add RoW based on GLO
        df = cfs.loc[cfs.loc[:, 'Elem flow name'].str.contains(', GLO')].copy()
        df.loc[:, 'Elem flow name'] = [i.split(', GLO')[0] + ', RoW' for i in df.loc[:, 'Elem flow name']]
        df.loc[:, 'Native geographical resolution scale'] = 'Other region'
        cfs = clean_up_dataframe(pd.concat([cfs, df]))

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
        cfs.loc[:, 'CAS number'] = '14265-44-2'
        cfs = cfs.reset_index()
        cfs.columns = ['Region code', 'CF value', 'CF unit', 'MP or Damage', 'Impact category', 'Compartment',
                       'Sub-compartment', 'Elem flow unit', 'Native geographical resolution scale',
                       'Elem flow', 'CAS number']

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - Stoechiometry]', self.conn)
        stoc = stoc.loc[stoc.loc[:, 'Impact category'] == 'Freshwater eutrophication']
        for i in stoc.index:
            proxy = stoc.loc[i, 'Proxy molecule']
            comp = stoc.loc[i, 'Compartment']
            df = cfs[cfs.loc[:, 'Elem flow'] == proxy].loc[cfs.loc[:, 'Compartment'] == comp].copy('deep')
            if not df.empty:
                df.loc[:, 'Elem flow'] = stoc.loc[i, 'Elem flow name']
                df.loc[:, 'CAS number'] = stoc.loc[i, 'CAS number']
                df.loc[:, 'CF value'] *= stoc.loc[i, 'Proxy ratio']

                cfs = clean_up_dataframe(pd.concat([cfs, df]))

        # ------------------------------------ FINAL FORMATTING ---------------------------------
        cfs.loc[:, 'Elem flow name'] = [', '.join(i) for i in
                                        list(zip(cfs.loc[:, 'Elem flow'], cfs.loc[:, 'Region code']))]
        cfs = cfs.drop(['Elem flow', 'Region code'], axis=1)
        cfs = clean_up_dataframe(cfs)

        # add value for RME (Middle East) based on RAS (Asia)
        df = cfs.loc[cfs.loc[:, 'Elem flow name'].str.contains(', RAS')].copy()
        df.loc[:, 'Elem flow name'] = [i.split(', RAS')[0] + ', RME' for i in df.loc[:, 'Elem flow name']]
        cfs = clean_up_dataframe(pd.concat([cfs, df]))
        # add RoW based on GLO
        df = cfs.loc[cfs.loc[:, 'Elem flow name'].str.contains(', GLO')].copy()
        df.loc[:, 'Elem flow name'] = [i.split(', GLO')[0] + ', RoW' for i in df.loc[:, 'Elem flow name']]
        df.loc[:, 'Native geographical resolution scale'] = 'Other region'
        cfs = clean_up_dataframe(pd.concat([cfs, df]))

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
        original_cf_occup = pd.read_sql(sql='SELECT * FROM [CF - regionalized - Land occupation per biome]',
                                        con=self.conn)

        intersect_country = pd.read_sql(sql='SELECT * FROM [SI - Land occupation - countries cell resolution]',
                                        con=self.conn)

        land_use_type = pd.read_sql(sql='SELECT * FROM [SI - Land occupation - land use type per country/continent]',
                                    con=self.conn)

        recovery_times = pd.read_sql(sql='SELECT * FROM [SI - Land transformation - recovery times]', con=self.conn)

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

        # assign global default value to zero values (except for "not used" land type)
        for indicator in ['Occupation, ', 'Transformation, from ', 'Transformation, to ']:
            for land_type in ['secondary vegetation', 'forest, used', 'pasture/meadow', 'annual crops',
                              'permanent crops',
                              'agriculture, mosaic (agroforestry)', 'artificial areas', 'unspecified']:
                for level in ['Midpoint', 'Damage']:
                    to_adapt_to_global_default = (land_cfs.loc[land_cfs.loc[:, 'CF value'] == 0].loc[
                                                      land_cfs.loc[:, 'MP or Damage'] == level].loc[
                                                      land_cfs.loc[:, 'Elem flow name'].str.contains(
                                                          indicator + land_type.split('(')[0])].index)
                    land_cfs.loc[to_adapt_to_global_default, 'CF value'] = (
                        land_cfs.loc[land_cfs.loc[:, 'MP or Damage'] == level].loc[
                        land_cfs.loc[:, 'Elem flow name'] == indicator + land_type + ', GLO', 'CF value'].iloc[0])

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
        conc = pd.read_sql(sql='SELECT * FROM [SI - Mapping with regions of ecoinvent]', con=self.conn).set_index(
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

        # add ISO 2-letter codes to countries
        coco_logger = coco.logging.getLogger()
        coco_logger.setLevel(coco.logging.CRITICAL)
        data.loc[:, 'Country_code'] = coco.convert(data.Country, to='ISO2')

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

        # add other regions (e.g., Andorra)
        for region in data.Country_code:
            if region not in conc.index and region != 'not found':
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(
                                                    data.loc[data.loc[:, 'Country_code'] == region, 'CF urban'].iloc[0],
                                                    index=[('Particulates, < 2.5 um, ' + region, 'high. pop.')])])
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(
                                                    data.loc[data.loc[:, 'Country_code'] == region, 'CF rural'].iloc[0],
                                                    index=[('Particulates, < 2.5 um, ' + region, 'low. pop.')])])
                particulate_damage = pd.concat([particulate_damage,
                                                pd.Series(data.loc[data.loc[:,
                                                                   'Country_code'] == region, 'CF unspecified'].iloc[0],
                                                          index=[
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

        so2 = particulate_damage.copy()
        so2.loc[:, 'Elem flow name'] = [i.replace('Particulates, < 2.5 um', 'Sulfur dioxide') for i in
                                        so2.loc[:, 'Elem flow name']]
        so2.loc[:, 'CAS number'] = '007446-09-5'

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
            elif region in data.Country_code:
                if so2.loc[flow, 'Sub-compartment'] == 'high. pop.':
                    so2.loc[flow, 'CF value'] = (
                                so2.loc[flow, 'CF value'] / data.loc[data.Country_code == region, 'iF urban'].iloc[0] *
                                secondary_pm_if.loc['SO2', 'urban'])
                elif so2.loc[flow, 'Sub-compartment'] == 'low. pop.':
                    so2.loc[flow, 'CF value'] = (
                                so2.loc[flow, 'CF value'] / data.loc[data.Country_code == region, 'iF rural'].iloc[0] *
                                secondary_pm_if.loc['SO2', 'rural'])
                elif so2.loc[flow, 'Sub-compartment'] == '(unspecified)':
                    so2.loc[flow, 'CF value'] = (so2.loc[flow, 'CF value'] /
                                                 data.loc[data.Country_code == region, 'iF unspecified'].iloc[0] *
                                                 secondary_pm_if.loc['SO2', 'unspecified'])
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
        nh3.loc[:, 'CAS number'] = '007664-41-7'

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
            elif region in data.Country_code:
                if nh3.loc[flow, 'Sub-compartment'] == 'high. pop.':
                    nh3.loc[flow, 'CF value'] = (nh3.loc[flow, 'CF value'] /
                                                 data.loc[data.Country_code == region, 'iF urban'].iloc[0] *
                                                 secondary_pm_if.loc['NH3', 'urban'])
                elif nh3.loc[flow, 'Sub-compartment'] == 'low. pop.':
                    nh3.loc[flow, 'CF value'] = (nh3.loc[flow, 'CF value'] /
                                                 data.loc[data.Country_code == region, 'iF rural'].iloc[0] *
                                                 secondary_pm_if.loc['NH3', 'rural'])
                elif nh3.loc[flow, 'Sub-compartment'] == '(unspecified)':
                    nh3.loc[flow, 'CF value'] = (nh3.loc[flow, 'CF value'] /
                                                 data.loc[data.Country_code == region, 'iF unspecified'].iloc[
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
        nox.loc[:, 'CAS number'] = '011104-93-1'

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
            elif region in data.Country_code:
                if nox.loc[flow, 'Sub-compartment'] == 'high. pop.':
                    nox.loc[flow, 'CF value'] = (
                                nox.loc[flow, 'CF value'] / data.loc[data.Country_code == region, 'iF urban'].iloc[0] *
                                secondary_pm_if.loc['NOx', 'urban'])
                elif nox.loc[flow, 'Sub-compartment'] == 'low. pop.':
                    nox.loc[flow, 'CF value'] = (
                                nox.loc[flow, 'CF value'] / data.loc[data.Country_code == region, 'iF rural'].iloc[0] *
                                secondary_pm_if.loc['NOx', 'rural'])
                elif nox.loc[flow, 'Sub-compartment'] == '(unspecified)':
                    nox.loc[flow, 'CF value'] = (nox.loc[flow, 'CF value'] /
                                                 data.loc[data.Country_code == region, 'iF unspecified'].iloc[0] *
                                                 secondary_pm_if.loc['NOx', 'unspecified'])
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
        particulate_cfs.loc[:, 'Native geographical resolution scale'] = 'Country'
        particulate_cfs.loc[[i for i in particulate_cfs.index if particulate_cfs.loc[i, 'Elem flow name'].split(', ')[
            -1] in continents], 'Native geographical resolution scale'] = 'Other region'
        particulate_cfs.loc[[i for i in particulate_cfs.index if
                             particulate_cfs.loc[i, 'Elem flow name'].split(', ')[-1] in [
                                 'RER', 'RLA', 'RNA', 'RAS', 'RAF', 'OCE', 'RME']],
                            'Native geographical resolution scale'] = 'Continent'
        particulate_cfs.loc[[i for i in particulate_cfs.index if
                             particulate_cfs.loc[i, 'Elem flow name'].split(', ')[-1] == 'GLO'],
                            'Native geographical resolution scale'] = 'Global'
        particulate_cfs.loc[[i for i in particulate_cfs.index if
                             particulate_cfs.loc[i, 'Elem flow name'].split(', ')[-1] == 'RoW'],
                            'Native geographical resolution scale'] = 'Other region'

        # add zero values for PMs above 2.5um
        big_pms = particulate_cfs.loc[[i for i in particulate_cfs.index if (
                particulate_cfs.loc[i, 'Elem flow name'] == 'Particulates, < 2.5 um, GLO' and
                particulate_cfs.loc[i, 'Sub-compartment'] == '(unspecified)')]].copy()
        big_pms = pd.concat([big_pms] * 2)
        big_pms.loc[:, 'CF value'] = 0
        big_pms.iloc[0, big_pms.columns.get_loc('Elem flow name')] = 'Particulates, > 10 um'
        big_pms.iloc[1, big_pms.columns.get_loc('Elem flow name')] = 'Particulates, > 10 um'
        big_pms.iloc[2, big_pms.columns.get_loc('Elem flow name')] = 'Particulates, > 2.5 um, and < 10um'
        big_pms.iloc[3, big_pms.columns.get_loc('Elem flow name')] = 'Particulates, > 2.5 um, and < 10um'

        self.master_db = pd.concat([self.master_db, particulate_cfs, big_pms])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_scarcity_cfs(self):
        """
        Load CFs for water scarcity

        Concerned impact categories:
            - Water scarcity

        :return: updated master_db
        """

        data = pd.read_sql('SELECT * FROM "CF - regionalized - WaterScarcity - aggregated"', self.conn)
        conc = pd.read_sql('SELECT * FROM "SI - Mapping with regions of ecoinvent"', self.conn).set_index(
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
        all_data.loc[[i for i in all_data.index if ('RoW' in all_data.loc[i, 'Elem flow name'] or
                                                    'Akrotiri' in all_data.loc[i, 'Elem flow name'] or
                                                    'Asia without China' in all_data.loc[i, 'Elem flow name'] or
                                                    'Australia, including overseas territories' in all_data.loc[i, 'Elem flow name'] or
                                                    'BALTSO' in all_data.loc[i, 'Elem flow name'] or
                                                    'CENTREL' in all_data.loc[i, 'Elem flow name'] or
                                                    'CUSMA/T-MEC/USMCA' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta and Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canary Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'Central Asia' in all_data.loc[i, 'Elem flow name'] or
                                                    'China w/o Inner Mongol' in all_data.loc[i, 'Elem flow name'] or
                                                    'Crimea' in all_data.loc[i, 'Elem flow name'] or
                                                    'Cyprus No Mans Area' in all_data.loc[i, 'Elem flow name'] or
                                                    'Dhekelia Base' in all_data.loc[i, 'Elem flow name'] or
                                                    'ENTSO-E' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Austria' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without NORDEL (NCPA)' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland and Austria' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland and France' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe, without Russia and Türkiye' in all_data.loc[i, 'Elem flow name'] or
                                                    'FSU' in all_data.loc[i, 'Elem flow name'] or
                                                    'France, including overseas territories' in all_data.loc[i, 'Elem flow name'] or
                                                    'Guantanamo Bay' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Africa' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Asia, without China and GCC' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, EU27 & EFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Gulf Cooperation Council' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, North America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, North America, without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Russia & RER w/o EU27 & EFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, South America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IN-Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'MRO' in all_data.loc[i, 'Elem flow name'] or
                                                    'NAFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'NORDEL' in all_data.loc[i, 'Elem flow name'] or
                                                    'NPCC' in all_data.loc[i, 'Elem flow name'] or
                                                    'North America without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Northern Cyprus' in all_data.loc[i, 'Elem flow name'] or
                                                    'Québec, HQ distribution network' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o AT+BE+CH+DE+FR+IT' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o CH+DE' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o DE+NL+NO' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o DE+NL+NO+RU' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o DE+NL+RU' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o RU' in all_data.loc[i, 'Elem flow name'] or
                                                    'Russia (Asia)' in all_data.loc[i, 'Elem flow name'] or
                                                    'Russia (Europe)' in all_data.loc[i, 'Elem flow name'] or
                                                    'SAS' in all_data.loc[i, 'Elem flow name'] or
                                                    'Siachen Glacier' in all_data.loc[i, 'Elem flow name'] or
                                                    'Somaliland' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE without France' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE without Germany' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE without Germany and France' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-AMERICAS' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-ASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-AUSTRALIANZ' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-CAMERICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-CARIBBEAN' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EEUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-MAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-MELANESIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-MICRONESIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-NAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-NEUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-OCEANIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-POLYNESIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SAMERICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SEASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SEUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-WAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-WASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'United States of America, including overseas territories' in all_data.loc[i, 'Elem flow name'] or
                                                    'WECC' in all_data.loc[i, 'Elem flow name'] or
                                                    'WEU' in all_data.loc[i, 'Elem flow name']
                                                    )],
                         'Native geographical resolution scale'] = 'Other region'

        # adding the different other water flows (lake, river, well, etc.)
        df_lake = all_data.loc[[i for i in all_data.index if
                                'agri' not in all_data.loc[i, 'Elem flow name'] and all_data.loc[
                                    i, 'Compartment'] == 'Raw']].copy()
        df_lake.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, lake') for i in df_lake.loc[:, 'Elem flow name']]
        df_river = all_data.loc[[i for i in all_data.index if
                                 'agri' not in all_data.loc[i, 'Elem flow name'] and all_data.loc[
                                     i, 'Compartment'] == 'Raw']].copy()
        df_river.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, river') for i in
                                             df_river.loc[:, 'Elem flow name']]
        df_unspe = all_data.loc[[i for i in all_data.index if
                                 'agri' not in all_data.loc[i, 'Elem flow name'] and all_data.loc[
                                     i, 'Compartment'] == 'Raw']].copy()
        df_unspe.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, unspecified natural origin') for i in
                                             df_unspe.loc[:, 'Elem flow name']]
        df_well = all_data.loc[[i for i in all_data.index if
                                'agri' not in all_data.loc[i, 'Elem flow name'] and all_data.loc[
                                    i, 'Compartment'] == 'Raw']].copy()
        df_well.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, well, in ground') for i in
                                            df_well.loc[:, 'Elem flow name']]
        df_cooling = all_data.loc[[i for i in all_data.index if
                                   'agri' not in all_data.loc[i, 'Elem flow name'] and all_data.loc[
                                       i, 'Compartment'] == 'Raw']].copy()
        df_cooling.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, cooling, unspecified natural origin') for i in
                                               df_cooling.loc[:, 'Elem flow name']]

        # drop "Water" flow in Raw comp, that flow is only for the water comp
        all_data = all_data.drop([i for i in all_data.index if
                       'agri' not in all_data.loc[i, 'Elem flow name'] and all_data.loc[i, 'Compartment'] == 'Raw'])

        all_data = clean_up_dataframe(pd.concat([all_data, df_lake, df_river, df_unspe, df_well, df_cooling]))

        self.master_db = pd.concat([self.master_db, all_data])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_availability_hh_cfs(self):
        """
        Load CFs for water availability human health.

        Concerned impact categories:
            - Water availability, human health

        :return: update master_db
        """

        data = pd.read_sql('SELECT * FROM "CF - regionalized - WaterAvailability_HH - aggregated"', self.conn).loc[
               :, ['ecoinvent_shortname', 'CF_tot']]
        conc = pd.read_sql('SELECT * FROM "SI - Mapping with regions of ecoinvent"', self.conn).set_index(
            'AWARE').Ecoinvent_short_name

        # create the regionalized names (e.g., Water, AF)
        water_data = pd.DataFrame()
        for i in data.index:
            if data.ecoinvent_shortname[i] in conc.index:
                if type(conc.loc[data.ecoinvent_shortname[i]]) == str:
                    water_data.loc[i, 'Elem flow name'] = 'Water, ' + conc.loc[data.ecoinvent_shortname[i]]
                    water_data.loc[i, 'CF value'] = data.loc[i, 'CF_tot']
                else:
                    for j in range(0, len(conc.loc[data.ecoinvent_shortname[i]])):
                        water_data.loc[str(i + j), 'Elem flow name'] = 'Water, ' + conc.loc[data.ecoinvent_shortname[i]].iloc[j]
                        water_data.loc[str(i + j), 'CF value'] = data.loc[i, 'CF_tot']
            else:
                water_data.loc[i, 'Elem flow name'] = 'Water, ' + data.ecoinvent_shortname[i]
                water_data.loc[i, 'CF value'] = data.loc[i, 'CF_tot']

        # formatting the data to IW+ format
        water_data.loc[:, 'Impact category'] = 'Water availability, human health'
        water_data.loc[:, 'CF unit'] = 'DALY'
        water_data.loc[:, 'Compartment'] = 'Raw'
        water_data.loc[:, 'Sub-compartment'] = '(unspecified)'
        water_data.loc[:, 'CAS number'] = '7732-18-5'
        water_data.loc[:, 'Elem flow unit'] = 'm3'
        water_data.loc[:, 'MP or Damage'] = 'Damage'
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
        all_data.loc[[i for i in all_data.index if ('RoW' in all_data.loc[i, 'Elem flow name'] or
                                                    'Akrotiri' in all_data.loc[i, 'Elem flow name'] or
                                                    'Asia without China' in all_data.loc[i, 'Elem flow name'] or
                                                    'Australia, including overseas territories' in all_data.loc[i, 'Elem flow name'] or
                                                    'BALTSO' in all_data.loc[i, 'Elem flow name'] or
                                                    'CENTREL' in all_data.loc[i, 'Elem flow name'] or
                                                    'CUSMA/T-MEC/USMCA' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta and Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canary Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'Central Asia' in all_data.loc[i, 'Elem flow name'] or
                                                    'China w/o Inner Mongol' in all_data.loc[i, 'Elem flow name'] or
                                                    'Crimea' in all_data.loc[i, 'Elem flow name'] or
                                                    'Cyprus No Mans Area' in all_data.loc[i, 'Elem flow name'] or
                                                    'Dhekelia Base' in all_data.loc[i, 'Elem flow name'] or
                                                    'ENTSO-E' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Austria' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without NORDEL (NCPA)' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland and Austria' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland and France' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe, without Russia and Türkiye' in all_data.loc[i, 'Elem flow name'] or
                                                    'FSU' in all_data.loc[i, 'Elem flow name'] or
                                                    'France, including overseas territories' in all_data.loc[i, 'Elem flow name'] or
                                                    'Guantanamo Bay' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Africa' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Asia, without China and GCC' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, EU27 & EFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Gulf Cooperation Council' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, North America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, North America, without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Russia & RER w/o EU27 & EFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, South America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IN-Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'MRO' in all_data.loc[i, 'Elem flow name'] or
                                                    'NAFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'NORDEL' in all_data.loc[i, 'Elem flow name'] or
                                                    'NPCC' in all_data.loc[i, 'Elem flow name'] or
                                                    'North America without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Northern Cyprus' in all_data.loc[i, 'Elem flow name'] or
                                                    'Québec, HQ distribution network' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o AT+BE+CH+DE+FR+IT' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o CH+DE' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o DE+NL+NO' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o DE+NL+NO+RU' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o DE+NL+RU' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o RU' in all_data.loc[i, 'Elem flow name'] or
                                                    'Russia (Asia)' in all_data.loc[i, 'Elem flow name'] or
                                                    'Russia (Europe)' in all_data.loc[i, 'Elem flow name'] or
                                                    'SAS' in all_data.loc[i, 'Elem flow name'] or
                                                    'Siachen Glacier' in all_data.loc[i, 'Elem flow name'] or
                                                    'Somaliland' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE without France' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE without Germany' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE without Germany and France' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-AMERICAS' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-ASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-AUSTRALIANZ' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-CAMERICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-CARIBBEAN' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EEUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-MAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-MELANESIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-MICRONESIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-NAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-NEUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-OCEANIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-POLYNESIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SAMERICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SEASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SEUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-WAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-WASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'United States of America, including overseas territories' in all_data.loc[i, 'Elem flow name'] or
                                                    'WECC' in all_data.loc[i, 'Elem flow name'] or
                                                    'WEU' in all_data.loc[i, 'Elem flow name']
                                                    )],
                         'Native geographical resolution scale'] = 'Other region'

        # adding the different other water flows (lake, river, well, etc.)
        df_lake = all_data[all_data.Compartment == 'Raw'].copy()
        df_lake.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, lake') for i in df_lake.loc[:, 'Elem flow name']]
        df_river = all_data[all_data.Compartment == 'Raw'].copy()
        df_river.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, river') for i in
                                             df_river.loc[:, 'Elem flow name']]
        df_unspe = all_data[all_data.Compartment == 'Raw'].copy()
        df_unspe.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, unspecified natural origin') for i in
                                             df_unspe.loc[:, 'Elem flow name']]
        df_well = all_data[all_data.Compartment == 'Raw'].copy()
        df_well.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, well, in ground') for i in
                                            df_well.loc[:, 'Elem flow name']]
        df_cooling = all_data[all_data.Compartment == 'Raw'].copy()
        df_cooling.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, cooling, unspecified natural origin') for i in
                                               df_cooling.loc[:, 'Elem flow name']]

        # drop "Water" flow in Raw comp, that flow is only for the water comp
        all_data = all_data[all_data.Compartment != 'Raw']

        all_data = clean_up_dataframe(pd.concat([all_data, df_lake, df_river, df_unspe, df_well, df_cooling]))

        # for missing CFs, forced value to zero
        all_data.loc[:, 'CF value'] = all_data.loc[:, 'CF value'].fillna(0)

        self.master_db = pd.concat([self.master_db, all_data])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_availability_fw_cfs(self):
        """
        Load CFs for water availability freshwater ecosystem.

        Concerned impact categories:
            - Water availability, freshwater ecosystem

        :return: update master_db
        """

        data = pd.read_sql(sql='SELECT * FROM [CF - regionalized - WaterAvailability_EQ_fw - native]', con=self.conn)
        conc = pd.read_sql('SELECT * FROM "SI - Mapping with regions of ecoinvent"', self.conn).set_index(
            'AWARE').Ecoinvent_short_name
        geos = pd.read_sql(sql='SELECT * FROM [CF - regionalized - WaterScarcity - aggregated]', con=self.conn).loc[
               :, 'ecoinvent_shortname'].tolist()

        CF_value = data.loc[:, 'CF value'].median()
        # create the regionalized names (e.g., Water, AF)
        water_data = pd.DataFrame()
        for i, geo in enumerate(geos):
            if geo in conc.index:
                if type(conc.loc[geo]) == str:
                    water_data.loc[i, 'Elem flow name'] = 'Water, ' + conc.loc[geo]
                    water_data.loc[i, 'CF value'] = CF_value
                else:
                    for j in range(0, len(conc.loc[geo])):
                        water_data.loc[str(i + j), 'Elem flow name'] = 'Water, ' + conc.loc[geo].iloc[j]
                        water_data.loc[str(i + j), 'CF value'] = CF_value
            else:
                water_data.loc[i, 'Elem flow name'] = 'Water, ' + geo
                water_data.loc[i, 'CF value'] = CF_value

        # formatting the data to IW+ format
        water_data.loc[:, 'Impact category'] = 'Water availability, freshwater ecosystem'
        water_data.loc[:, 'CF unit'] = 'PDF.m2.yr'
        water_data.loc[:, 'Compartment'] = 'Raw'
        water_data.loc[:, 'Sub-compartment'] = '(unspecified)'
        water_data.loc[:, 'CAS number'] = '7732-18-5'
        water_data.loc[:, 'Elem flow unit'] = 'm3'
        water_data.loc[:, 'MP or Damage'] = 'Damage'
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
        all_data.loc[[i for i in all_data.index if ('RoW' in all_data.loc[i, 'Elem flow name'] or
                                                    'Akrotiri' in all_data.loc[i, 'Elem flow name'] or
                                                    'Asia without China' in all_data.loc[i, 'Elem flow name'] or
                                                    'Australia, including overseas territories' in all_data.loc[i, 'Elem flow name'] or
                                                    'BALTSO' in all_data.loc[i, 'Elem flow name'] or
                                                    'CENTREL' in all_data.loc[i, 'Elem flow name'] or
                                                    'CUSMA/T-MEC/USMCA' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta and Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canary Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'Central Asia' in all_data.loc[i, 'Elem flow name'] or
                                                    'China w/o Inner Mongol' in all_data.loc[i, 'Elem flow name'] or
                                                    'Crimea' in all_data.loc[i, 'Elem flow name'] or
                                                    'Cyprus No Mans Area' in all_data.loc[i, 'Elem flow name'] or
                                                    'Dhekelia Base' in all_data.loc[i, 'Elem flow name'] or
                                                    'ENTSO-E' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Austria' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without NORDEL (NCPA)' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland and Austria' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland and France' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe, without Russia and Türkiye' in all_data.loc[i, 'Elem flow name'] or
                                                    'FSU' in all_data.loc[i, 'Elem flow name'] or
                                                    'France, including overseas territories' in all_data.loc[i, 'Elem flow name'] or
                                                    'Guantanamo Bay' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Africa' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Asia, without China and GCC' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, EU27 & EFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Gulf Cooperation Council' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, North America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, North America, without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Russia & RER w/o EU27 & EFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, South America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IN-Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'MRO' in all_data.loc[i, 'Elem flow name'] or
                                                    'NAFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'NORDEL' in all_data.loc[i, 'Elem flow name'] or
                                                    'NPCC' in all_data.loc[i, 'Elem flow name'] or
                                                    'North America without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Northern Cyprus' in all_data.loc[i, 'Elem flow name'] or
                                                    'Québec, HQ distribution network' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o AT+BE+CH+DE+FR+IT' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o CH+DE' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o DE+NL+NO' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o DE+NL+NO+RU' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o DE+NL+RU' in all_data.loc[i, 'Elem flow name'] or
                                                    'RER w/o RU' in all_data.loc[i, 'Elem flow name'] or
                                                    'Russia (Asia)' in all_data.loc[i, 'Elem flow name'] or
                                                    'Russia (Europe)' in all_data.loc[i, 'Elem flow name'] or
                                                    'SAS' in all_data.loc[i, 'Elem flow name'] or
                                                    'Siachen Glacier' in all_data.loc[i, 'Elem flow name'] or
                                                    'Somaliland' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE without France' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE without Germany' in all_data.loc[i, 'Elem flow name'] or
                                                    'UCTE without Germany and France' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-AMERICAS' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-ASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-AUSTRALIANZ' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-CAMERICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-CARIBBEAN' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EEUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-EUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-MAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-MELANESIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-MICRONESIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-NAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-NEUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-OCEANIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-POLYNESIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SAMERICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SEASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-SEUROPE' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-WAFRICA' in all_data.loc[i, 'Elem flow name'] or
                                                    'UN-WASIA' in all_data.loc[i, 'Elem flow name'] or
                                                    'United States of America, including overseas territories' in all_data.loc[i, 'Elem flow name'] or
                                                    'WECC' in all_data.loc[i, 'Elem flow name'] or
                                                    'WEU' in all_data.loc[i, 'Elem flow name']
                                                    )],
                         'Native geographical resolution scale'] = 'Other region'

        # adding the different other water flows (lake, river, well, etc.)
        df_lake = all_data[all_data.Compartment == 'Raw'].copy()
        df_lake.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, lake') for i in df_lake.loc[:, 'Elem flow name']]
        df_river = all_data[all_data.Compartment == 'Raw'].copy()
        df_river.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, river') for i in
                                             df_river.loc[:, 'Elem flow name']]
        df_unspe = all_data[all_data.Compartment == 'Raw'].copy()
        df_unspe.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, unspecified natural origin') for i in
                                             df_unspe.loc[:, 'Elem flow name']]
        df_well = all_data[all_data.Compartment == 'Raw'].copy()
        df_well.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, well, in ground') for i in
                                            df_well.loc[:, 'Elem flow name']]
        df_cooling = all_data[all_data.Compartment == 'Raw'].copy()
        df_cooling.loc[:, 'Elem flow name'] = [i.replace('Water', 'Water, cooling, unspecified natural origin') for i in
                                               df_cooling.loc[:, 'Elem flow name']]

        # drop "Water" flow in Raw comp, that flow is only for the water comp
        all_data = all_data[all_data.Compartment != 'Raw']

        all_data = clean_up_dataframe(pd.concat([all_data, df_lake, df_river, df_unspe, df_well, df_cooling]))

        # concat with master_db
        self.master_db = pd.concat([self.master_db, all_data])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_availability_terr_cfs(self):
        """
        Load CFs for water availability terrestrial ecosystem.

        Concerned impact categories:
            - Water availability, terrestrial ecosystem

        :return: update master_db
        """

        data = pd.read_sql('SELECT * FROM "CF - regionalized - WaterAvailability_EQ_terr - aggregated"', self.conn)
        conc = pd.read_sql('SELECT * FROM "SI - Mapping with regions of ecoinvent"', self.conn).set_index(
            'AWARE').Ecoinvent_short_name

        # create the regionalized names (e.g., Water, AF)
        water_data = pd.DataFrame()
        for i in data.index:
            if data.ecoinvent_shortname[i] in conc.index:
                if type(conc.loc[data.ecoinvent_shortname[i]]) == str:
                    water_data.loc[i, 'Elem flow name'] = 'Water, well, in ground, ' + conc.loc[data.ecoinvent_shortname[i]]
                    water_data.loc[i, 'CF value'] = data.loc[i, 'CF (PDF.m2.yr/m3)']
                else:
                    for j in range(0, len(conc.loc[data.ecoinvent_shortname[i]])):
                        water_data.loc[str(i + j), 'Elem flow name'] = 'Water, well, in ground, ' + conc.loc[data.ecoinvent_shortname[i]].iloc[j]
                        water_data.loc[str(i + j), 'CF value'] = data.loc[i, 'CF (PDF.m2.yr/m3)']
            else:
                water_data.loc[i, 'Elem flow name'] = 'Water, well, in ground, ' + data.ecoinvent_shortname[i]
                water_data.loc[i, 'CF value'] = data.loc[i, 'CF (PDF.m2.yr/m3)']

        water_data.loc[:, 'Impact category'] = 'Water availability, terrestrial ecosystem'
        water_data.loc[:, 'CF unit'] = 'PDF.m2.yr'
        water_data.loc[:, 'Compartment'] = 'Raw'
        water_data.loc[:, 'Sub-compartment'] = '(unspecified)'
        water_data.loc[:, 'CAS number'] = '7732-18-5'
        water_data.loc[:, 'Elem flow unit'] = 'm3'
        water_data.loc[:, 'MP or Damage'] = 'Damage'
        water_data.loc[:, 'Native geographical resolution scale'] = 'Country'
        water_data.loc[:, 'CF value'] = water_data.loc[:, 'CF value'].astype(float)

        water_data.loc[[i for i in water_data.index if ('RER' in water_data.loc[i, 'Elem flow name'] or
                                            'RAS' in water_data.loc[i, 'Elem flow name'] or
                                            'RAF' in water_data.loc[i, 'Elem flow name'] or
                                            'RLA' in water_data.loc[i, 'Elem flow name'] or
                                            'OCE' in water_data.loc[i, 'Elem flow name'] or
                                            'RNA' in water_data.loc[i, 'Elem flow name'])],
                         'Native geographical resolution scale'] = 'Continent'
        water_data.loc[[i for i in water_data.index if 'GLO' in water_data.loc[i, 'Elem flow name']],
                         'Native geographical resolution scale'] = 'Global'
        water_data.loc[[i for i in water_data.index if ('RoW' in water_data.loc[i, 'Elem flow name'] or
                                            'Akrotiri' in water_data.loc[i, 'Elem flow name'] or
                                            'Asia without China' in water_data.loc[i, 'Elem flow name'] or
                                            'Australia, including overseas territories' in water_data.loc[i, 'Elem flow name'] or
                                            'BALTSO' in water_data.loc[i, 'Elem flow name'] or
                                            'CENTREL' in water_data.loc[i, 'Elem flow name'] or
                                            'CUSMA/T-MEC/USMCA' in water_data.loc[i, 'Elem flow name'] or
                                            'Canada without Alberta' in water_data.loc[i, 'Elem flow name'] or
                                            'Canada without Alberta and Quebec' in water_data.loc[i, 'Elem flow name'] or
                                            'Canada without Quebec' in water_data.loc[i, 'Elem flow name'] or
                                            'Canary Islands' in water_data.loc[i, 'Elem flow name'] or
                                            'Central Asia' in water_data.loc[i, 'Elem flow name'] or
                                            'China w/o Inner Mongol' in water_data.loc[i, 'Elem flow name'] or
                                            'Crimea' in water_data.loc[i, 'Elem flow name'] or
                                            'Cyprus No Mans Area' in water_data.loc[i, 'Elem flow name'] or
                                            'Dhekelia Base' in water_data.loc[i, 'Elem flow name'] or
                                            'ENTSO-E' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe without Austria' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe without NORDEL (NCPA)' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe without Switzerland' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe without Switzerland and Austria' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe without Switzerland and France' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe, without Russia and Türkiye' in water_data.loc[i, 'Elem flow name'] or
                                            'FSU' in water_data.loc[i, 'Elem flow name'] or
                                            'France, including overseas territories' in water_data.loc[i, 'Elem flow name'] or
                                            'Guantanamo Bay' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, Africa' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, Asia, without China and GCC' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, EU27 & EFTA' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, Gulf Cooperation Council' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, North America' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, North America, without Quebec' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, Russia & RER w/o EU27 & EFTA' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, South America' in water_data.loc[i, 'Elem flow name'] or
                                            'IN-Islands' in water_data.loc[i, 'Elem flow name'] or
                                            'MRO' in water_data.loc[i, 'Elem flow name'] or
                                            'NAFTA' in water_data.loc[i, 'Elem flow name'] or
                                            'NORDEL' in water_data.loc[i, 'Elem flow name'] or
                                            'NPCC' in water_data.loc[i, 'Elem flow name'] or
                                            'North America without Quebec' in water_data.loc[i, 'Elem flow name'] or
                                            'Northern Cyprus' in water_data.loc[i, 'Elem flow name'] or
                                            'Québec, HQ distribution network' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o AT+BE+CH+DE+FR+IT' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o CH+DE' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o DE+NL+NO' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o DE+NL+NO+RU' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o DE+NL+RU' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o RU' in water_data.loc[i, 'Elem flow name'] or
                                            'Russia (Asia)' in water_data.loc[i, 'Elem flow name'] or
                                            'Russia (Europe)' in water_data.loc[i, 'Elem flow name'] or
                                            'SAS' in water_data.loc[i, 'Elem flow name'] or
                                            'Siachen Glacier' in water_data.loc[i, 'Elem flow name'] or
                                            'Somaliland' in water_data.loc[i, 'Elem flow name'] or
                                            'UCTE' in water_data.loc[i, 'Elem flow name'] or
                                            'UCTE without France' in water_data.loc[i, 'Elem flow name'] or
                                            'UCTE without Germany' in water_data.loc[i, 'Elem flow name'] or
                                            'UCTE without Germany and France' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-AMERICAS' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-ASIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-AUSTRALIANZ' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-CAMERICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-CARIBBEAN' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-EAFRICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-EASIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-EEUROPE' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-EUROPE' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-MAFRICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-MELANESIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-MICRONESIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-NAFRICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-NEUROPE' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-OCEANIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-POLYNESIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-SAMERICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-SASIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-SEASIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-SEUROPE' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-WAFRICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-WASIA' in water_data.loc[i, 'Elem flow name'] or
                                            'United States of America, including overseas territories' in water_data.loc[i, 'Elem flow name'] or
                                            'WECC' in water_data.loc[i, 'Elem flow name'] or
                                            'WEU' in water_data.loc[i, 'Elem flow name'])],
                         'Native geographical resolution scale'] = 'Other region'

        self.master_db = pd.concat([self.master_db, water_data])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_thermally_polluted_water_cfs(self):
        """
        Load CFs for thermally polluted water.

        Concerned impact categories:
            - Thermally polluted water

        :return: update master_db
        """

        data = pd.read_sql('SELECT * FROM "CF - not regionalized - ThermallyPollutedWater"', self.conn)
        conc = pd.read_sql('SELECT * FROM "SI - Mapping with regions of ecoinvent"', self.conn).set_index(
            'AWARE').Ecoinvent_short_name
        geos = pd.read_sql(sql='SELECT * FROM [CF - regionalized - WaterScarcity - aggregated]', con=self.conn).loc[
               :, 'ecoinvent_shortname'].tolist()

        CF_value = data.loc[:, 'CF value'].iloc[0]
        # create the regionalized names (e.g., Water, AF)
        water_data = pd.DataFrame()
        for i, geo in enumerate(geos):
            if geo in conc.index:
                if type(conc.loc[geo]) == str:
                    water_data.loc[i, 'Elem flow name'] = 'Water, cooling, unspecified natural origin, ' + conc.loc[geo]
                    water_data.loc[i, 'CF value'] = CF_value
                else:
                    for j in range(0, len(conc.loc[geo])):
                        water_data.loc[str(i + j), 'Elem flow name'] = 'Water, cooling, unspecified natural origin, ' + conc.loc[geo].iloc[j]
                        water_data.loc[str(i + j), 'CF value'] = CF_value
            else:
                water_data.loc[i, 'Elem flow name'] = 'Water, cooling, unspecified natural origin, ' + geo
                water_data.loc[i, 'CF value'] = CF_value

        # formatting the data to IW+ format
        water_data.loc[:, 'Impact category'] = 'Thermally polluted water'
        water_data.loc[:, 'CF unit'] = 'PDF.m2.yr'
        water_data.loc[:, 'Compartment'] = 'Raw'
        water_data.loc[:, 'Sub-compartment'] = '(unspecified)'
        water_data.loc[:, 'CAS number'] = '7732-18-5'
        water_data.loc[:, 'Elem flow unit'] = 'm3'
        water_data.loc[:, 'MP or Damage'] = 'Damage'
        water_data.loc[:, 'Native geographical resolution scale'] = 'Country'
        water_data.loc[:, 'CF value'] = water_data.loc[:, 'CF value'].astype(float)

        water_data.loc[[i for i in water_data.index if ('RER' in water_data.loc[i, 'Elem flow name'] or
                                            'RAS' in water_data.loc[i, 'Elem flow name'] or
                                            'RAF' in water_data.loc[i, 'Elem flow name'] or
                                            'RLA' in water_data.loc[i, 'Elem flow name'] or
                                            'OCE' in water_data.loc[i, 'Elem flow name'] or
                                            'RNA' in water_data.loc[i, 'Elem flow name'])],
                         'Native geographical resolution scale'] = 'Continent'
        water_data.loc[[i for i in water_data.index if 'GLO' in water_data.loc[i, 'Elem flow name']],
                         'Native geographical resolution scale'] = 'Global'
        water_data.loc[[i for i in water_data.index if ('RoW' in water_data.loc[i, 'Elem flow name'] or
                                            'Akrotiri' in water_data.loc[i, 'Elem flow name'] or
                                            'Asia without China' in water_data.loc[i, 'Elem flow name'] or
                                            'Australia, including overseas territories' in water_data.loc[i, 'Elem flow name'] or
                                            'BALTSO' in water_data.loc[i, 'Elem flow name'] or
                                            'CENTREL' in water_data.loc[i, 'Elem flow name'] or
                                            'CUSMA/T-MEC/USMCA' in water_data.loc[i, 'Elem flow name'] or
                                            'Canada without Alberta' in water_data.loc[i, 'Elem flow name'] or
                                            'Canada without Alberta and Quebec' in water_data.loc[i, 'Elem flow name'] or
                                            'Canada without Quebec' in water_data.loc[i, 'Elem flow name'] or
                                            'Canary Islands' in water_data.loc[i, 'Elem flow name'] or
                                            'Central Asia' in water_data.loc[i, 'Elem flow name'] or
                                            'China w/o Inner Mongol' in water_data.loc[i, 'Elem flow name'] or
                                            'Crimea' in water_data.loc[i, 'Elem flow name'] or
                                            'Cyprus No Mans Area' in water_data.loc[i, 'Elem flow name'] or
                                            'Dhekelia Base' in water_data.loc[i, 'Elem flow name'] or
                                            'ENTSO-E' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe without Austria' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe without NORDEL (NCPA)' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe without Switzerland' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe without Switzerland and Austria' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe without Switzerland and France' in water_data.loc[i, 'Elem flow name'] or
                                            'Europe, without Russia and Türkiye' in water_data.loc[i, 'Elem flow name'] or
                                            'FSU' in water_data.loc[i, 'Elem flow name'] or
                                            'France, including overseas territories' in water_data.loc[i, 'Elem flow name'] or
                                            'Guantanamo Bay' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, Africa' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, Asia, without China and GCC' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, EU27 & EFTA' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, Gulf Cooperation Council' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, North America' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, North America, without Quebec' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, Russia & RER w/o EU27 & EFTA' in water_data.loc[i, 'Elem flow name'] or
                                            'IAI Area, South America' in water_data.loc[i, 'Elem flow name'] or
                                            'IN-Islands' in water_data.loc[i, 'Elem flow name'] or
                                            'MRO' in water_data.loc[i, 'Elem flow name'] or
                                            'NAFTA' in water_data.loc[i, 'Elem flow name'] or
                                            'NORDEL' in water_data.loc[i, 'Elem flow name'] or
                                            'NPCC' in water_data.loc[i, 'Elem flow name'] or
                                            'North America without Quebec' in water_data.loc[i, 'Elem flow name'] or
                                            'Northern Cyprus' in water_data.loc[i, 'Elem flow name'] or
                                            'Québec, HQ distribution network' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o AT+BE+CH+DE+FR+IT' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o CH+DE' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o DE+NL+NO' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o DE+NL+NO+RU' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o DE+NL+RU' in water_data.loc[i, 'Elem flow name'] or
                                            'RER w/o RU' in water_data.loc[i, 'Elem flow name'] or
                                            'Russia (Asia)' in water_data.loc[i, 'Elem flow name'] or
                                            'Russia (Europe)' in water_data.loc[i, 'Elem flow name'] or
                                            'SAS' in water_data.loc[i, 'Elem flow name'] or
                                            'Siachen Glacier' in water_data.loc[i, 'Elem flow name'] or
                                            'Somaliland' in water_data.loc[i, 'Elem flow name'] or
                                            'UCTE' in water_data.loc[i, 'Elem flow name'] or
                                            'UCTE without France' in water_data.loc[i, 'Elem flow name'] or
                                            'UCTE without Germany' in water_data.loc[i, 'Elem flow name'] or
                                            'UCTE without Germany and France' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-AMERICAS' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-ASIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-AUSTRALIANZ' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-CAMERICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-CARIBBEAN' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-EAFRICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-EASIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-EEUROPE' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-EUROPE' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-MAFRICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-MELANESIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-MICRONESIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-NAFRICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-NEUROPE' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-OCEANIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-POLYNESIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-SAMERICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-SASIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-SEASIA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-SEUROPE' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-WAFRICA' in water_data.loc[i, 'Elem flow name'] or
                                            'UN-WASIA' in water_data.loc[i, 'Elem flow name'] or
                                            'United States of America, including overseas territories' in water_data.loc[i, 'Elem flow name'] or
                                            'WECC' in water_data.loc[i, 'Elem flow name'] or
                                            'WEU' in water_data.loc[i, 'Elem flow name'])],
                         'Native geographical resolution scale'] = 'Other region'

        # concat with master_db
        self.master_db = pd.concat([self.master_db, water_data])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_plastic_cfs(self):
        """
        Load CFs for plastics impact.

        Concerned impact categories:
            - Plastics physical effect on biota

        :return: update master_db
        """

        original_cfs = pd.read_sql('SELECT * from "CF - not regionalized - PlasticPhysicalImpactonBiota"',
                                   con=self.conn)

        original_cfs.drop(['Geometric st.dev.', 'Lower limit 95% CI (calculated)',
                           'Upper limit 95% CI (calculated)'], axis=1, inplace=True)
        original_cfs.loc[:, 'Impact category'] = 'Plastics physical effects on biota'
        original_cfs.loc[:, 'Native geographical resolution scale'] = 'Global'
        original_cfs.loc[:, 'MP or Damage'] = ['Midpoint' if original_cfs.loc[i, 'CF unit'] == 'CTUe' else 'Damage' for
                                               i in original_cfs.index]
        CAS = {'EPS': '9003-53-6',
               'HDPE': '9002-88-4',
               'LDPE': '9002-88-4',
               'PA/Nylon': '25038-54-4',
               'PET': '25038-59-9',
               'PHA': '117068-64-1',
               'PLA': '26100-51-6',
               'PP': '9003-07-0',
               'PS': '9003-53-6',
               'PVC': '9002-86-2',
               'TRWP': None,
               'PU/Spandex': '9009-54-5',
               'Acrylic': '9065-11-6'}
        shape_names = {'Beads/spheres': 'Microplastic beads',
                       'Film fragments': 'Microplastic film fragments',
                       'Microfibers/cylinders': 'Plastic microfibers'}
        original_cfs.loc[:, 'CAS number'] = [CAS[i] for i in original_cfs.loc[:, 'Polymer type']]
        original_cfs.loc[:, 'Shape'] = [shape_names[i] for i in original_cfs.loc[:, 'Shape']]
        original_cfs = original_cfs.rename(columns={'Recommended CF (geometric mean)': 'CF value'})
        original_cfs.loc[:, 'Elem flow name'] = [
            original_cfs.loc[i, 'Shape'] + ' - ' + original_cfs.loc[i, 'Polymer type'] + ' (' +
            str(original_cfs.loc[i, 'Size']) + ' µm diameter)' for i in original_cfs.index]

        # we create CFs for default sizes which depend on the shape of the microplastics
        default_sizes = {'Microplastic beads': 1000,
                         'Microplastic film fragments': 100,
                         'Plastic microfibers': 10}
        for shape in default_sizes.keys():
            df = original_cfs.loc[[i for i in original_cfs.index if
                                   shape == original_cfs.loc[i, 'Shape'] and default_sizes[shape] == original_cfs.loc[
                                       i, 'Size']]].copy()
            df.loc[:, 'Elem flow name'] = [i.split(str(default_sizes[shape]) + ' µm diameter')[0] + 'default)' for i in
                                           df.loc[:, 'Elem flow name']]
            original_cfs = clean_up_dataframe(pd.concat([original_cfs, df]))

        original_cfs.drop(['Polymer type', 'Size', 'Shape'], axis=1, inplace=True)

        self.master_db = pd.concat([self.master_db, original_cfs])
        self.master_db = clean_up_dataframe(self.master_db)

        # have an unspecified sub-compartment in the water compartment for plastics
        df = self.master_db.loc[self.master_db.loc[:, 'Impact category'] == 'Plastics physical effects on biota'].loc[
            self.master_db.loc[:, 'Sub-compartment'] == 'ocean'].copy()
        df.loc[:, 'Sub-compartment'] = '(unspecified)'
        self.master_db = clean_up_dataframe(pd.concat([self.master_db, df]))

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
            [pd.DataFrame([i[1] + ' fish, FAO zone ' + i[0] for i in cf_regions.index], columns=['Elem flow name']),
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

    def harmonize_regionalized_substances(self):
        """
        Between different regionalized impact categories, the spatialization precision is not the same. One category
        can define the region "AD" (Andorra) for Ammonia, while another does not provide a CF for AD. Not solving these
        disparities means that when using Ammonia, AD, it is not characterized on all the indicators. So we identify
        these problemes and assign the continental value to each (e.g., RER for AD)
        :return: updated self.master_db
        """

        map = pd.read_sql('SELECT * FROM [SI - Mapping countries to continents]', self.conn).set_index('country')

        # keep both as separate loops, first one has to run before the second one
        for indicator in ['Particulate matter formation', 'Freshwater acidification', 'Terrestrial acidification',
                          'Marine eutrophication']:
            for substance in ['Ammonia', 'Nitrogen oxides', 'Sulfur dioxide']:
                if indicator == 'Particulate matter formation':
                    df = self.master_db.loc[self.master_db.loc[:, 'Impact category'] == indicator].loc[
                        self.master_db.loc[:, 'Elem flow name'] == substance + ', UN-OCEANIA'].copy('deep')
                    df.loc[:, 'Elem flow name'] = substance + ', OCE'
                    self.master_db = clean_up_dataframe(pd.concat([self.master_db, df]))

                    df = self.master_db.loc[self.master_db.loc[:, 'Impact category'] == indicator].loc[
                        self.master_db.loc[:, 'Elem flow name'] == substance + ', US-PR'].copy('deep')
                    df.loc[:, 'Elem flow name'] = substance + ', PR'
                    self.master_db = clean_up_dataframe(pd.concat([self.master_db, df]))
                else:
                    df = self.master_db.loc[self.master_db.loc[:, 'Impact category'] == indicator].loc[
                        self.master_db.loc[:, 'Elem flow name'] == substance + ', OCE'].copy('deep')
                    df.loc[:, 'Elem flow name'] = substance + ', UN-OCEANIA'
                    self.master_db = clean_up_dataframe(pd.concat([self.master_db, df]))

                    df = self.master_db.loc[self.master_db.loc[:, 'Impact category'] == indicator].loc[
                        self.master_db.loc[:, 'Elem flow name'] == substance + ', PR'].copy('deep')
                    df.loc[:, 'Elem flow name'] = substance + ', US-PR'
                    self.master_db = clean_up_dataframe(pd.concat([self.master_db, df]))

        for substance in ['Ammonia', 'Nitrogen oxides', 'Sulfur dioxide']:
            all_existing_geos = set(self.master_db.loc[self.master_db.loc[:, 'Elem flow name'].str.contains(
                substance + ', '), 'Elem flow name'])
            all_existing_geos = [i.split(substance + ', ')[1] for i in all_existing_geos]

            for indicator in ['Particulate matter formation', 'Freshwater acidification', 'Terrestrial acidification',
                              'Marine eutrophication']:
                existing_geos = set(self.master_db.loc[self.master_db.loc[:, 'Impact category'] == indicator].loc[
                                        self.master_db.loc[:, 'Elem flow name'].str.contains(
                                            substance + ', '), 'Elem flow name'])
                existing_geos = [i.split(substance + ', ')[1] for i in existing_geos]

                to_create = set(all_existing_geos) - set(existing_geos)

                for new_geo in to_create:
                    if new_geo in map.index:
                        df = self.master_db.loc[self.master_db.loc[:, 'Impact category'] == indicator].loc[
                            self.master_db.loc[:, 'Elem flow name'] == substance + ', ' + map.loc[
                                new_geo, 'continent']].copy('deep')
                        df.loc[:, 'Elem flow name'] = substance + ', ' + new_geo
                        if ('-' in new_geo and new_geo not in ['ENTSO-E', 'UN-SEASIA']) or len(new_geo) == 2:
                            df.loc[:, 'Native geographical resolution scale'] = 'Country'
                        else:
                            df.loc[:, 'Native geographical resolution scale'] = 'Other region'
                        self.master_db = clean_up_dataframe(pd.concat([self.master_db, df]))

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
                                   'Marine ecotoxicity, long term',
                                   'Marine ecotoxicity, short term',
                                   'Terrestrial ecotoxicity, long term',
                                   'Terrestrial ecotoxicity, short term',
                                   'Freshwater eutrophication',
                                   'Human toxicity cancer',
                                   'Human toxicity cancer, long term',
                                   'Human toxicity cancer, short term',
                                   'Human toxicity non-cancer',
                                   'Human toxicity non-cancer, long term',
                                   'Human toxicity non-cancer, short term',
                                   'Ionizing radiations, ecosystem quality',
                                   'Ionizing radiations, human health',
                                   'Ionizing radiations',
                                   'Marine eutrophication'],
                   'groundwater, long-term': ['Freshwater ecotoxicity',
                                              'Freshwater ecotoxicity, long term',
                                              'Freshwater ecotoxicity, short term',
                                              'Marine ecotoxicity, long term',
                                              'Marine ecotoxicity, short term',
                                              'Terrestrial ecotoxicity, long term',
                                              'Terrestrial ecotoxicity, short term',
                                              'Freshwater eutrophication',
                                              'Human toxicity cancer',
                                              'Human toxicity cancer, long term',
                                              'Human toxicity cancer, short term',
                                              'Human toxicity non-cancer',
                                              'Human toxicity non-cancer, long term',
                                              'Human toxicity non-cancer, short term',
                                              'Ionizing radiations, ecosystem quality',
                                              'Ionizing radiations, human health',
                                              'Ionizing radiations',
                                              'Marine eutrophication'],
                   'ocean': ['Freshwater ecotoxicity',
                             'Freshwater ecotoxicity, long term',
                             'Freshwater ecotoxicity, short term',
                             'Marine ecotoxicity, long term',
                             'Marine ecotoxicity, short term',
                             'Terrestrial ecotoxicity, long term',
                             'Terrestrial ecotoxicity, short term',
                             'Freshwater eutrophication',
                             'Ionizing radiations, ecosystem quality',
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

        long_term_cats = ['Climate change, ecosystem quality', 'Climate change, human health', 'Freshwater ecotoxicity',
                          'Marine ecotoxicity', 'Terrestrial ecotoxicity', 'Human toxicity cancer',
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
        ics = ['Marine eutrophication', 'Ozone layer depletion', 'Photochemical ozone formation, human health',
               'Photochemical ozone formation, ecosystem quality',
               'Terrestrial acidification', 'Particulate matter formation', 'Ionizing radiations, ecosystem quality',
               'Ionizing radiations, human health', 'Freshwater acidification']
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

        # ----------------- Saline water --------------------

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

        regio_flows_per_ic = {}

        for ic in set(self.master_db.loc[:, 'Impact category']):
            df = self.master_db[self.master_db['Impact category'] == ic].copy()
            regio_flows_per_ic[ic] = set([', '.join(i.split(', ')[:-1]) for i in df['Elem flow name'] if ', FR' in i])

        flows_to_create = {}

        for ic in set(self.master_db.loc[:, 'Impact category']):
            df = self.master_db[self.master_db['Impact category'] == ic].copy()
            for flow in set.union(*[set(i) for i in list(regio_flows_per_ic.values())]):
                if flow in df['Elem flow name'].tolist():
                    # if the CF is equal to zero we don't care
                    if df.loc[[i for i in df.index if (df.loc[i, 'Elem flow name'] == flow and
                                                       df.loc[i, 'Impact category'] == ic)], 'CF value'].sum() != 0:
                        if ic in flows_to_create.keys():
                            flows_to_create[ic].append(flow)
                        else:
                            flows_to_create[ic] = [flow]

        flows_to_create = {k: set(v) for k, v in flows_to_create.items()}

        for ic in flows_to_create.keys():
            flows_to_create[ic] = flows_to_create[ic] - regio_flows_per_ic[ic]

        flows_to_create = {k: v for k, v in flows_to_create.items() if v != set()}

        for substance in set.union(*[set(i) for i in list(flows_to_create.values())]):
            df = self.master_db.loc[
                self.master_db[self.master_db['Elem flow name'].str.contains(substance, na=False)].index.tolist()]
            df = df[df.loc[:, 'Native geographical resolution scale'].isin(['Country', 'Continent', 'Other region'])]
            regions = set([i.split(substance + ', ')[-1] for i in df.loc[:, 'Elem flow name'] if substance + ', ' in i])
            ic_flows_to_create = {k for k, v in flows_to_create.items() if substance in v}
            for ic in ic_flows_to_create:
                dff = self.master_db.loc[self.master_db.loc[:, 'Elem flow name'] == substance].loc[
                    self.master_db.loc[:, 'Impact category'] == ic].copy()
                list_flow_added = ([substance + ', ' + i for i in regions] * len(dff))
                list_flow_added.sort()
                dff = pd.concat([dff] * (len(regions)))
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

    def deal_with_biogenic_carbon(self):
        """
        Biogenic carbon can be followed either with the carbon neutrality approach (where e.g., CO2 bio = 0 kgCO2eq and
        Methane, bio = 27 kgCO2eq) or with the +/-1 approach. Here we deal with both. The default option (i.e.,
        self.master_db) follows the +/-1 approach while we create a copy self.master_db_carbon_neutrality where the
        carbon neutrality approach is followed.
        :return:
        """

        self.master_db_carbon_neutrality = self.master_db.copy()

        co2_bio_release = self.master_db.loc[
            self.master_db.loc[:, 'Elem flow name'].str.contains('Carbon dioxide')].copy()
        co2_bio_release.loc[:, 'Elem flow name'] = 'Carbon dioxide, biogenic, release'

        co2_bio_uptake = self.master_db.loc[
            self.master_db.loc[:, 'Elem flow name'].str.contains('Carbon dioxide')].copy()
        co2_bio_uptake.loc[:, 'Elem flow name'] = 'Carbon dioxide, biogenic, uptake'
        co2_bio_uptake.loc[:, 'CF value'] = -co2_bio_uptake.loc[:, 'CF value']

        co_bio_release = self.master_db.loc[
            self.master_db.loc[:, 'Elem flow name'].str.contains('Carbon monoxide')].copy()
        co_bio_release.loc[:, 'Elem flow name'] = 'Carbon monoxide, biogenic'

        self.master_db = self.master_db.drop(
            self.master_db.loc[self.master_db.loc[:, 'Elem flow name'].str.contains('Methane, biogenic')].index)
        ch4_bio_release = self.master_db.loc[
            self.master_db.loc[:, 'Elem flow name'].str.contains('Methane, fossil')].copy()
        ch4_bio_release.loc[:, 'Elem flow name'] = 'Methane, biogenic'

        co2_to_soil = self.master_db.loc[
            self.master_db.loc[:, 'Elem flow name'].str.contains('Carbon dioxide')].loc[
            self.master_db.loc[:, 'Sub-compartment'] == '(unspecified)'].copy()
        co2_to_soil.loc[:, 'Elem flow name'] = 'Carbon dioxide, to soil or biomass stock'
        co2_to_soil.loc[:, 'Compartment'] = 'Soil'
        co2_to_soil.loc[:, 'CF value'] = -co2_to_soil.loc[:, 'CF value']
        df = co2_to_soil.copy()
        df.loc[:, 'Sub-compartment'] = 'industrial'
        co2_to_soil = clean_up_dataframe(pd.concat([co2_to_soil, df]))
        df = co2_to_soil.copy()
        df.loc[:, 'Sub-compartment'] = 'agricultural'
        co2_to_soil = clean_up_dataframe(pd.concat([co2_to_soil, df]))
        df = co2_to_soil.copy()
        df.loc[:, 'Sub-compartment'] = 'forestry'
        co2_to_soil = clean_up_dataframe(pd.concat([co2_to_soil, df]))

        self.master_db = clean_up_dataframe(pd.concat([self.master_db, co2_bio_release, co2_bio_uptake,
                                                       co_bio_release, co2_to_soil, ch4_bio_release]))

        co2_bio_release = self.master_db_carbon_neutrality.loc[
            self.master_db_carbon_neutrality.loc[:, 'Elem flow name'].str.contains('Carbon dioxide')].copy()
        co2_bio_release.loc[:, 'Elem flow name'] = 'Carbon dioxide, biogenic, release'
        co2_bio_release.loc[:, 'CF value'] = 0

        co2_bio_uptake = self.master_db_carbon_neutrality.loc[
            self.master_db_carbon_neutrality.loc[:, 'Elem flow name'].str.contains('Carbon dioxide')].copy()
        co2_bio_uptake.loc[:, 'Elem flow name'] = 'Carbon dioxide, biogenic, uptake'
        co2_bio_uptake.loc[:, 'CF value'] = 0

        co_bio_release = self.master_db_carbon_neutrality.loc[
            self.master_db_carbon_neutrality.loc[:, 'Elem flow name'].str.contains('Carbon monoxide')].copy()
        co_bio_release.loc[:, 'Elem flow name'] = 'Carbon monoxide, biogenic'
        co_bio_release.loc[:, 'CF value'] = 0

        co2_to_soil = self.master_db_carbon_neutrality.loc[
            self.master_db_carbon_neutrality.loc[:, 'Elem flow name'].str.contains('Carbon dioxide')].loc[
            self.master_db_carbon_neutrality.loc[:, 'Sub-compartment'] == '(unspecified)'].copy()
        co2_to_soil.loc[:, 'Elem flow name'] = 'Carbon dioxide, to soil or biomass stock'
        co2_to_soil.loc[:, 'Compartment'] = 'Soil'
        co2_to_soil.loc[:, 'CF value'] = -co2_to_soil.loc[:, 'CF value']
        df = co2_to_soil.copy()
        df.loc[:, 'Sub-compartment'] = 'industrial'
        co2_to_soil = clean_up_dataframe(pd.concat([co2_to_soil, df]))
        df = co2_to_soil.copy()
        df.loc[:, 'Sub-compartment'] = 'agricultural'
        co2_to_soil = clean_up_dataframe(pd.concat([co2_to_soil, df]))
        df = co2_to_soil.copy()
        df.loc[:, 'Sub-compartment'] = 'forestry'
        co2_to_soil = clean_up_dataframe(pd.concat([co2_to_soil, df]))

        self.master_db_carbon_neutrality = clean_up_dataframe(
            pd.concat([self.master_db_carbon_neutrality, co2_bio_release, co2_bio_uptake,
                       co_bio_release, co2_to_soil]))

        self.master_db_carbon_neutrality.loc[[i for i in self.master_db_carbon_neutrality.index if (
                'Marine acidification' in self.master_db_carbon_neutrality.loc[i, 'Impact category'] and
                self.master_db_carbon_neutrality.loc[i, 'Elem flow name'] == 'Methane, biogenic')], 'CF value'] = 0

    def deal_with_temporary_storage_of_carbon(self):
        """
        Some LCI databases cover flows of temporary storage of carbon (in kgy). The associated CF is simply 1/100 of the
        normal CF.
        :return:
        """

        def temporary_storage_ghg_cf(origin_flow, name_storage_flow, master_db_format):
            """Function takes in the origin GHG flow name (e.g., Carbon dioxide, fossil) and a new name for the
            storage flow (e.g., Correction flow for delayed emission of biogenic carbon dioxide). It then creates CFs
            for the temporary storage flows."""
            df = master_db_format.loc[master_db_format.loc[:, 'Elem flow name'] == origin_flow].loc[~
            master_db_format.loc[:, 'Impact category'].isin(['Fossil and nuclear energy use',
                                                             'Climate change, short term',
                                                             'Climate change, long term'])].copy()
            df.loc[:, 'Elem flow name'] = name_storage_flow
            df.loc[:, 'Elem flow unit'] = 'kgy'
            df.loc[:, 'CF value'] /= -100
            df.loc[df.loc[:, 'Impact category'] == 'Climate change, long term', 'CF value'] = 0
            df.loc[df.loc[:, 'Impact category'] == 'Climate change, ecosystem quality, long term', 'CF value'] = -(
                df.loc[df.loc[:, 'Impact category'] == 'Climate change, ecosystem quality, short term', 'CF value']
            ).iloc[0]
            df.loc[df.loc[:, 'Impact category'] == 'Climate change, human health, long term', 'CF value'] = -(
                df.loc[df.loc[:, 'Impact category'] == 'Climate change, human health, short term', 'CF value']
            ).iloc[0]
            if origin_flow in ['Carbon dioxide, fossil', 'Methane, fossil', 'Carbon monoxide, fossil']:
                df.loc[df.loc[:, 'Impact category'] == 'Marine acidification, long term', 'CF value'] = -(
                    df.loc[df.loc[:, 'Impact category'] == 'Marine acidification, short term', 'CF value']
                ).iloc[0]
            return df

        self.master_db = clean_up_dataframe(pd.concat([
            self.master_db,
            temporary_storage_ghg_cf('Carbon dioxide, biogenic, release',
                                     'Correction flow for delayed emission of biogenic carbon dioxide',
                                     self.master_db),
            temporary_storage_ghg_cf('Carbon dioxide, fossil',
                                     'Correction flow for delayed emission of fossil carbon dioxide',
                                     self.master_db),
            temporary_storage_ghg_cf('Methane, biogenic',
                                     'Correction flow for delayed emission of biogenic methane',
                                     self.master_db),
            temporary_storage_ghg_cf('Methane, fossil',
                                     'Correction flow for delayed emission of fossil methane',
                                     self.master_db),
            temporary_storage_ghg_cf('Dinitrogen monoxide',
                                     'Correction flow for delayed emission of nitrous oxide',
                                     self.master_db),
            temporary_storage_ghg_cf('Sulfur hexafluoride',
                                     'Correction flow for delayed emission of sulphur hexafluoride',
                                     self.master_db)]))

        self.master_db_carbon_neutrality = clean_up_dataframe(pd.concat([
            self.master_db_carbon_neutrality,
            temporary_storage_ghg_cf('Carbon dioxide, biogenic, release',
                                     'Correction flow for delayed emission of biogenic carbon dioxide',
                                     self.master_db_carbon_neutrality),
            temporary_storage_ghg_cf('Carbon dioxide, fossil',
                                     'Correction flow for delayed emission of fossil carbon dioxide',
                                     self.master_db_carbon_neutrality),
            temporary_storage_ghg_cf('Methane, biogenic',
                                     'Correction flow for delayed emission of biogenic methane',
                                     self.master_db_carbon_neutrality),
            temporary_storage_ghg_cf('Methane, fossil',
                                     'Correction flow for delayed emission of fossil methane',
                                     self.master_db_carbon_neutrality),
            temporary_storage_ghg_cf('Dinitrogen monoxide',
                                     'Correction flow for delayed emission of nitrous oxide',
                                     self.master_db_carbon_neutrality),
            temporary_storage_ghg_cf('Sulfur hexafluoride',
                                     'Correction flow for delayed emission of sulphur hexafluoride',
                                     self.master_db_carbon_neutrality)]))

    def separate_ghg_indicators(self):
        """
        It is recommended by standards to separate fossil, biogenic, carbon uptake and land transformation flows for
        the GWP100 indicator. Since IW+ has additional GHG-related indicators, we apply the same logic the those as
        well. This only applies to the -/+1 approach, so only to self.master_db and not to self.master_db_carbon_neutrality
        :return:
        """

        biogenic = ['Carbon dioxide, biogenic, release',
                    'Correction flow for delayed emission of biogenic carbon dioxide',
                    'Methane, biogenic', 'Correction flow for delayed emission of biogenic methane',
                    'Carbon monoxide, biogenic']
        land_use = ['Carbon dioxide, to soil or biomass stock']
        CO2_uptake = ['Carbon dioxide, biogenic, uptake']

        self.master_db.loc[
            self.master_db['Elem flow name'].isin(biogenic) &
            self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category'] = [
            i+', biogenic' for i in self.master_db.loc[
                self.master_db['Elem flow name'].isin(biogenic) &
                self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category']]

        self.master_db.loc[
            self.master_db['Elem flow name'].isin(land_use) &
            self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category'] = [
            i+', land transformation' for i in self.master_db.loc[
                self.master_db['Elem flow name'].isin(land_use) &
                self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category']]

        self.master_db.loc[
            self.master_db['Elem flow name'].isin(CO2_uptake) &
            self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category'] = [
            i+', CO2 uptake' for i in self.master_db.loc[
                self.master_db['Elem flow name'].isin(CO2_uptake) &
                self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category']]

        self.master_db.loc[
            ~self.master_db['Elem flow name'].isin(CO2_uptake + biogenic + land_use) &
            self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category'] = [
            i+', fossil' for i in self.master_db.loc[
                ~self.master_db['Elem flow name'].isin(CO2_uptake + biogenic + land_use) &
                self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category']]

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

        self.master_db_not_regio_carbon_neutrality = self.master_db_carbon_neutrality.loc[
            [i for i in self.master_db_carbon_neutrality.index if
             self.master_db_carbon_neutrality.loc[i, 'Native geographical resolution scale'] not in [
                 'Continent', 'Country', 'Other region']]].copy()

        # dropping flow names with ", GLO" in them
        self.master_db_not_regio_carbon_neutrality.drop(
            [i for i in self.master_db_not_regio_carbon_neutrality.index if ', GLO' in
             self.master_db_not_regio_carbon_neutrality.loc[i, 'Elem flow name']], inplace=True)

    def link_to_ecoinvent(self):
        """
        Function that links names of substance from IW+ to the names of ecoinvent.
        :return: self.ei35_iw, self.ei36_iw, self.ei371_iw, self.ei38_iw
        """

        latest_ei_version = '3.11'

        for db_format in ['normal', 'carbon neutrality']:
            if db_format == 'normal':
                ei_iw_db = self.master_db_not_regio.copy()
            elif db_format == 'carbon neutrality':
                ei_iw_db = self.master_db_not_regio_carbon_neutrality.copy()

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
            for subst in tqdm(unique_not_one_for_one, leave=True):
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

            # fix comp/subcomp of pesky biogenic carbon elementary flows
            ei_iw_db = ei_iw_db.drop(
                ei_iw_db.loc[ei_iw_db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                                     'Carbon dioxide, in air'])].loc[
                    ei_iw_db.loc[:, 'Sub-compartment'] != 'unspecified'].index)
            ei_iw_db.loc[ei_iw_db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                                 'Carbon dioxide, in air']), 'Compartment'] = 'natural resource'
            ei_iw_db.loc[ei_iw_db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                                 'Carbon dioxide, in air']), 'Sub-compartment'] = 'in air'

            if db_format == 'normal':
                # recategorize these specific flows
                ei_iw_db.loc[
                    ei_iw_db['Elem flow name'].isin(['Carbon dioxide, from soil or biomass stock',
                                                     'Carbon monoxide, from soil or biomass stock',
                                                     'Methane, from soil or biomass stock']) &
                    ei_iw_db['Impact category'].str.contains(', fossil', na=False), 'Impact category'] = [
                    i.replace(', fossil', ', land transformation') for i in ei_iw_db.loc[
                        ei_iw_db['Elem flow name'].isin(['Carbon dioxide, from soil or biomass stock',
                                                     'Carbon monoxide, from soil or biomass stock',
                                                     'Methane, from soil or biomass stock']) &
                        ei_iw_db['Impact category'].str.contains(', fossil', na=False), 'Impact category']]

                # start with latest available version of ecoinvent
                self.ei311_iw = ei_iw_db.copy('deep')
                only_in_311 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == '3.11'].dropna(
                    subset=['iw name']).loc[:, 'ecoinvent name'])

                self.ei310_iw = self.ei311_iw.drop([i for i in self.ei311_iw.index if self.ei311_iw.loc[i, 'Elem flow name'] in
                                                  only_in_311]).copy('deep')

                only_in_310 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == '3.10'].dropna(
                    subset=['iw name']).loc[:, 'ecoinvent name'])

                self.ei39_iw = self.ei310_iw.drop([i for i in self.ei310_iw.index if self.ei310_iw.loc[i, 'Elem flow name'] in
                                                  only_in_310]).copy('deep')

                only_in_39 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == 3.9].dropna(
                    subset=['iw name']).loc[:, 'ecoinvent name'])

                self.ei38_iw = self.ei39_iw.drop([i for i in self.ei39_iw.index if self.ei39_iw.loc[i, 'Elem flow name'] in
                                                  only_in_39]).copy('deep')

                metals_in_ei38 = {'Aluminium': 'Aluminium III',
                                  'Antimony': 'Antimony ion',
                                  'Arsenic': 'Arsenic ion',
                                  'Barium': 'Barium II',
                                  'Beryllium': 'Beryllium II',
                                  'Cadmium': 'Cadmium II',
                                  'Chromium': 'Chromium III',
                                  'Cobalt': 'Cobalt II',
                                  'Copper': 'Copper ion',
                                  'Iron': 'Iron ion',
                                  'Lead': 'Lead II',
                                  'Manganese': 'Manganese II',
                                  'Mercury': 'Mercury II',
                                  'Molybdenum': 'Molybdenum VI',
                                  'Nickel': 'Nickel II',
                                  'Selenium': 'Selenium IV',
                                  'Silver': 'Silver I',
                                  'Strontium': 'Strontium I',
                                  'Thallium': 'Thallium I',
                                  'Tin': 'Tin ion',
                                  'Vanadium': 'Vanadium V',
                                  'Zinc': 'Zinc II'}

                # special case for ei3.8, since metal names in ei38 are shared between ions and their metallic form
                for metal in metals_in_ei38:
                    df = self.ei311_iw.loc[self.ei311_iw.loc[:, 'Elem flow name'] == metals_in_ei38[metal]].copy()
                    df.loc[:, 'Elem flow name'] = metal
                    self.ei38_iw = clean_up_dataframe(pd.concat([self.ei38_iw, df]))

            elif db_format == 'carbon neutrality':

                self.ei311_iw_carbon_neutrality = ei_iw_db.copy('deep')

                only_in_311 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == '3.11'].dropna(
                    subset=['iw name']).loc[:, 'ecoinvent name'])

                self.ei310_iw_carbon_neutrality = self.ei311_iw_carbon_neutrality.drop(
                    [i for i in self.ei311_iw_carbon_neutrality.index if self.ei311_iw_carbon_neutrality.loc[i, 'Elem flow name'] in
                     only_in_311]).copy('deep')

                only_in_310 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == '3.10'].dropna(
                    subset=['iw name']).loc[:, 'ecoinvent name'])

                self.ei39_iw_carbon_neutrality = self.ei310_iw_carbon_neutrality.drop(
                    [i for i in self.ei310_iw_carbon_neutrality.index if self.ei310_iw_carbon_neutrality.loc[i, 'Elem flow name'] in
                     only_in_310]).copy('deep')

                only_in_39 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == 3.9].dropna(
                    subset=['iw name']).loc[:, 'ecoinvent name'])

                self.ei38_iw_carbon_neutrality = self.ei39_iw_carbon_neutrality.drop(
                    [i for i in self.ei39_iw_carbon_neutrality.index if self.ei39_iw_carbon_neutrality.loc[i, 'Elem flow name'] in
                     only_in_39]).copy('deep')

                metals_in_ei38 = {'Aluminium': 'Aluminium III',
                                  'Antimony': 'Antimony ion',
                                  'Arsenic': 'Arsenic ion',
                                  'Barium': 'Barium II',
                                  'Beryllium': 'Beryllium II',
                                  'Cadmium': 'Cadmium II',
                                  'Chromium': 'Chromium III',
                                  'Cobalt': 'Cobalt II',
                                  'Copper': 'Copper ion',
                                  'Iron': 'Iron ion',
                                  'Lead': 'Lead II',
                                  'Manganese': 'Manganese II',
                                  'Mercury': 'Mercury II',
                                  'Molybdenum': 'Molybdenum VI',
                                  'Nickel': 'Nickel II',
                                  'Selenium': 'Selenium IV',
                                  'Silver': 'Silver I',
                                  'Strontium': 'Strontium I',
                                  'Thallium': 'Thallium I',
                                  'Tin': 'Tin ion',
                                  'Vanadium': 'Vanadium V',
                                  'Zinc': 'Zinc II'}

                # special case for ei3.8, since metal names in ei38 are shared between ions and their metallic form
                for metal in metals_in_ei38:
                    df = self.ei311_iw_carbon_neutrality.loc[self.ei311_iw_carbon_neutrality.loc[:, 'Elem flow name'] == metals_in_ei38[metal]].copy()
                    df.loc[:, 'Elem flow name'] = metal
                    self.ei38_iw_carbon_neutrality = clean_up_dataframe(pd.concat([self.ei38_iw_carbon_neutrality, df]))

    def link_to_sp(self):
        """
        This method creates a SimaPro method with the IW+ characterization factors.
        :return:
        """

        def linking(db, carboneutral):

            carboneutral = carboneutral

            # -------------------------------- MAPPING -------------------------------------

            # apply the mapping with the different SP flow names
            sp = pd.read_excel(pkg_resources.resource_filename(__name__, '/Data/mappings/SP/sp_mapping.xlsx'), None)
            sp = clean_up_dataframe(pd.concat([sp['Non regionalized'], sp['Regionalized']]))
            sp = sp.drop('Unnamed: 0', axis=1).loc[:, ['Name', 'Name IW+']].dropna()
            differences = sp.loc[sp.Name != sp.loc[:, 'Name IW+']].set_index('Name IW+')
            double_iw_flow = sp.loc[sp.loc[:, 'Name IW+'].duplicated(), 'Name IW+'].tolist()
            sp = sp.set_index('Name IW+')

            # go to dictionaries because it's waaaay faster to process than dataframes
            iw_sp_dict = dict(zip(list(zip(db.loc[:, 'Elem flow name'],
                                           db.loc[:, 'Impact category'],
                                           db.loc[:, 'CF unit'],
                                           db.loc[:, 'Compartment'],
                                           db.loc[:, 'Sub-compartment'],
                                           db.loc[:, 'CAS number'],
                                           db.loc[:, 'Elem flow unit'],
                                           db.loc[:, 'MP or Damage'],
                                           db.loc[:, 'Native geographical resolution scale'],
                                           )),
                                  db.loc[:, 'CF value']))

            for diff in tqdm(differences.index, leave=True):
                if diff not in double_iw_flow:
                    # simply rename the key
                    for k in list(iw_sp_dict.keys()):
                        if k[0] == diff:
                            new_key = (differences.loc[diff, 'Name'],) + k[1:]
                            iw_sp_dict[new_key] = iw_sp_dict.pop(k)
                else:
                    for i in range(len(sp.loc[diff, 'Name'])):
                        # here we loop through the multiple SP names connected to one IW+ name
                        if i == len(sp.loc[diff, 'Name']):
                            # if it's the final SP name, we pop the original key
                            for k in list(iw_sp_dict.keys()):
                                if k[0] == diff:
                                    new_key = (sp.loc[diff, 'Name'].iloc[i],) + k[1:]
                                    iw_sp_dict[new_key] = iw_sp_dict.pop(k)
                        else:
                            for k in list(iw_sp_dict.keys()):
                                # otherwise we just add a new key to the dict
                                if k[0] == diff:
                                    new_key = (sp.loc[diff, 'Name'].iloc[i],) + k[1:]
                                    iw_sp_dict[new_key] = iw_sp_dict[k]

            # create the dataframe from the dictionary
            db = pd.DataFrame.from_dict(iw_sp_dict, orient='index')
            db.index = pd.MultiIndex.from_tuples(db.index)
            db = db.reset_index()
            db.columns = ['Elem flow name', 'Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                          'CAS number', 'Elem flow unit', 'MP or Damage', 'Native geographical resolution scale',
                          'CF value']

            # ------------------------------- SUBCOMPS ----------------------------------

            # need to change the unit of land occupation flows to match SP nomenclature
            db.loc[db.loc[:, 'Elem flow unit'] == 'm2.yr', 'Elem flow unit'] = 'm2a'

            # Some water names are reserved names in SimaPro, so we modify it
            db.loc[db['Elem flow name'] == 'Water', 'Elem flow name'] = 'Water/m3'
            db.loc[db['Elem flow name'] == 'Water, agri', 'Elem flow name'] = 'Water/m3, agri'
            db.loc[db['Elem flow name'] == 'Water, non-agri', 'Elem flow name'] = 'Water/m3, non-agri'

            # need an unspecified subcomp for mineral resource uses, for some databases in SP (e.g., Industry2.0)
            df = db.loc[
                [i for i in db.index if db.loc[i, 'Impact category'] == 'Mineral resources use']].copy()
            df['Sub-compartment'] = '(unspecified)'
            db = clean_up_dataframe(pd.concat([db, df]))

            # fix comp/subcomp of pesky biogenic carbon elementary flows
            db = db.drop(
                db.loc[db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                         'Carbon dioxide, in air'])].loc[
                    db.loc[:, 'Sub-compartment'] != '(unspecified)'].index)
            db.loc[db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                     'Carbon dioxide, in air']), 'Compartment'] = 'Raw'
            db.loc[db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                     'Carbon dioxide, in air']), 'Sub-compartment'] = 'in air'
            db = db.drop(db.loc[db.loc[:, 'Elem flow name'] == 'Carbon dioxide, biogenic, uptake'].index)

            # --------------------------- UNIT CONVERSIONS ----------------------------------

            mj_flows = sp.loc[sp.loc[:, 'Name'].str.contains('MJ'), 'Name']
            for flow in mj_flows:
                db.loc[db.loc[:, 'Elem flow name'] == flow, 'CF value'] = float(
                    flow.split(' MJ')[0].split(', ')[-1])

            gj_flows = sp.loc[sp.loc[:, 'Name'].str.contains('GJ'), 'Name']
            for flow in gj_flows:
                if 'IN-GJ' not in flow:
                    db.loc[[i for i in db.index if db.loc[i, 'Elem flow name'] == flow and
                            db.loc[i, 'Impact category'] == 'Fossil and nuclear energy use'], 'CF value'] = float(
                        flow.split(' GJ')[0].split(', ')[-1]) * 1000

            # ensure units are coherent for energy flows
            db.loc[db.loc[:, 'Elem flow name'].isin([i for i in mj_flows if 'per kg' in i]), 'Elem flow unit'] = 'kg'
            db.loc[db.loc[:, 'Elem flow name'].isin([i for i in mj_flows if 'per m3' in i]), 'Elem flow unit'] = 'm3'

            # some other flows from SP that require conversions because of units
            new_flows = {
                'Gas, natural/kg': 'Gas, natural/m3',
                'Gas, mine, off-gas, process, coal mining/kg': 'Gas, mine, off-gas, process, coal mining/m3',
                'Wood, unspecified, standing/kg': 'Wood, unspecified, standing/m3',
                'Wood (16.9 MJ/kg)': 'Wood, unspecified, standing/m3'
            }
            for flow in new_flows:
                if 'Gas' in flow:
                    density = 0.8  # kg/m3 (https://www.engineeringtoolbox.com/gas-density-d_158.html)
                elif 'Wood' in flow:
                    density = 600  # kg/m3 (https://www.engineeringtoolbox.com/wood-density-d_40.html)

                df = db[db['Elem flow name'] == new_flows[flow]].copy()
                df.loc[:, 'Elem flow name'] = flow
                df.loc[:, 'Elem flow unit'] = 'kg'
                df.loc[:, 'CF value'] /= density
                db = pd.concat([db, df])
            db = clean_up_dataframe(db)

            # some water flows are in kilograms instead of cubic meters...
            problems = ['Water',
                        'Water, barrage',
                        'Water, cooling, drinking',
                        'Water, cooling, salt, ocean',
                        'Water, cooling, unspecified natural origin/kg',
                        'Water, process, unspecified natural origin/kg',
                        'Water, unspecified natural origin/kg',
                        'Water, process, drinking',
                        'Water, cooling, surface',
                        'Water, cooling, well',
                        'Water, process, surface',
                        'Water, process, salt, ocean',
                        'Water, process, well',
                        'Water, groundwater consumption',
                        'Water, surface water consumption',
                        'Water, Saline water consumption',
                        'Water, thermoelectric groundwater consumption',
                        'Water, thermoelectric saline water consumption',
                        'Water, thermoelectric surface water consumption',
                        'Thermally polluted water',
                        'Turbined water/kg']

            for problem_child in problems:
                db.loc[db['Elem flow name'] == problem_child, 'Elem flow unit'] = 'kg'
                db.loc[db['Elem flow name'] == problem_child, 'CF value'] /= 1000

            # finally, SimaPro limits to 12 characters the size of the CF unit. We rename some to avoid issues
            db.loc[db.loc[:, 'CF unit'] == 'kg CFC-11 eq', 'CF unit'] = 'kg CFC11 eq'
            db.loc[db.loc[:, 'CF unit'] == 'm2 arable land eq', 'CF unit'] = 'm2 ar ld eq'
            db.loc[db.loc[:, 'CF unit'] == 'm2 arable land eq .yr', 'CF unit'] = 'm2 ar ld.yr eq'

            # ------------------------------ SUB-CATEGORIES SHENANIGANS ----------------------------

            if not carboneutral:
                db.loc[db['Elem flow name'].isin(['Carbon dioxide, land transformation',
                                                  'Carbon monoxide, land transformation',
                                                  'Methane, land transformation']) &
                    db['Impact category'].str.contains(', fossil', na=False), 'Impact category'] = [
                    i.replace(', fossil', ', land transformation') for i in db.loc[
                        db['Elem flow name'].isin(['Carbon dioxide, land transformation',
                                                   'Carbon monoxide, land transformation',
                                                   'Methane, land transformation']) &
                        db['Impact category'].str.contains(', fossil', na=False), 'Impact category']]

                # the names of sub-categories exceed 40 characters (i.e., the limit of SimaPro), so we rename them
                db.loc[db.loc[:, 'Impact category'].str.contains(
                    'Climate change, ecosystem quality, short term'), 'Impact category'] = [
                    'Climate change, EQ, ST' + i.split('Climate change, ecosystem quality, short term')[1] for i in
                    db.loc[db.loc[:, 'Impact category'].str.contains(
                            'Climate change, ecosystem quality, short term'), 'Impact category']]

                db.loc[db.loc[:, 'Impact category'].str.contains(
                    'Climate change, ecosystem quality, long term'), 'Impact category'] = [
                    'Climate change, EQ, LT' + i.split('Climate change, ecosystem quality, long term')[1] for i in
                    db.loc[db.loc[:, 'Impact category'].str.contains(
                            'Climate change, ecosystem quality, long term'), 'Impact category']]

                db.loc[db.loc[:, 'Impact category'].str.contains(
                    'Climate change, human health, short term'), 'Impact category'] = [
                    'Climate change, HH, ST' + i.split('Climate change, human health, short term')[1] for i in
                    db.loc[db.loc[:, 'Impact category'].str.contains(
                            'Climate change, human health, short term'), 'Impact category']]

                db.loc[db.loc[:, 'Impact category'].str.contains(
                    'Climate change, human health, long term'), 'Impact category'] = [
                    'Climate change, HH, LT' + i.split('Climate change, human health, long term')[1] for i in
                    db.loc[db.loc[:, 'Impact category'].str.contains(
                            'Climate change, human health, long term'), 'Impact category']]

            # ------------------------------ DOING SIMAPRO'S JOB ----------------------------

            # Préconsultants can't seem to be able to harmonize their own substance list therefore some pollutants do not
            # have the same spelling depending on the compartment. We "fix" that here by hardcoding the necessary changes.
            db.loc[[i for i in db.index if
                            db.loc[i, 'Elem flow name'] == '2-(Chloromethyl)-3-Chloro-1-Propene' and
                            db.Compartment[i] in ['Water', 'Soil']], 'Elem flow name'] = '2-(Chloromethyl)-3-chloro-1-propene'
            db.loc[[i for i in db.index if
                            db.loc[i, 'Elem flow name'] == '2-Methyl-2,4-Pentanediol' and
                            db.Compartment[i] in ['Water', 'Soil']], 'Elem flow name'] = '2-Methyl-2,4-pentanediol'
            db.loc[[i for i in db.index if
                            db.loc[i, 'Elem flow name'] == 'Alpha-Naphthylamine' and
                            db.Compartment[i] == 'Air'], 'Elem flow name'] = 'alpha-Naphthylamine'
            db.loc[[i for i in db.index if
                            db.loc[i, 'Elem flow name'] == 'Alpha-pinene' and
                            db.Compartment[i] == 'Air'], 'Elem flow name'] = 'alpha-Pinene'
            db.loc[[i for i in db.index if
                            db.loc[i, 'Elem flow name'] == 'Phosphorodithioic acid, O,O-diethyl ester' and
                            db.Compartment[i] == 'Air'], 'Elem flow name'] = 'Phosphorodithioic acid, o,o-diethyl ester'
            db.loc[[i for i in db.index if
                            db.loc[i, 'Elem flow name'] == 'Phosphorodithioic acid, O,O-dimethyl ester' and
                            db.Compartment[i] == 'Air'], 'Elem flow name'] = 'Phosphorodithioic acid, o,o-dimethyl ester'
            db.loc[[i for i in db.index if
                            db.loc[i, 'Elem flow name'] == 't-Butyl acetate' and
                            db.Compartment[i] in ['Water', 'Soil']], 'Elem flow name'] = 'T-butyl acetate'

            return db

        self.iw_sp = linking(self.master_db, carboneutral=False)
        self.iw_sp_carbon_neutrality = linking(self.master_db_carbon_neutrality, carboneutral=True)

    def link_to_olca(self):
        """
        This method creates an openLCA method with the IW+ characterization factors.
        :return:
        """

        def linking(db):

            # -------------------------------- MAPPING -------------------------------------

            olca = pd.read_excel(pkg_resources.resource_filename(
                __name__, '/Data/mappings/oLCA/v2.1.1/oLCA_mapping.xlsx'), None)
            olca = clean_up_dataframe(pd.concat([olca['Non regionalized'], olca['Regionalized']]))
            olca = olca.drop('Unnamed: 0', axis=1).loc[:, ['Name', 'Name IW+']].dropna()
            differences = olca.loc[olca.Name != olca.loc[:, 'Name IW+']].set_index('Name IW+')
            double_iw_flow = olca.loc[olca.loc[:, 'Name IW+'].duplicated(), 'Name IW+'].tolist()
            olca = olca.set_index('Name IW+')

            # go to dictionaries because it's waaaay faster to process than dataframes
            iw_olca_dict = dict(zip(list(zip(db.loc[:, 'Elem flow name'],
                                             db.loc[:, 'Impact category'],
                                             db.loc[:, 'CF unit'],
                                             db.loc[:, 'Compartment'],
                                             db.loc[:, 'Sub-compartment'],
                                             db.loc[:, 'CAS number'],
                                             db.loc[:, 'Elem flow unit'],
                                             db.loc[:, 'MP or Damage'],
                                             db.loc[:, 'Native geographical resolution scale'],
                                             )),
                                    db.loc[:, 'CF value']))

            for diff in tqdm(differences.index, leave=True):
                if diff not in double_iw_flow:
                    # simply rename the key
                    for k in list(iw_olca_dict.keys()):
                        if k[0] == diff:
                            new_key = (differences.loc[diff, 'Name'],) + k[1:]
                            iw_olca_dict[new_key] = iw_olca_dict.pop(k)
                else:
                    for i in range(len(olca.loc[diff, 'Name'])):
                        # here we loop through the multiple oLCA names connected to one IW+ name
                        if i == len(olca.loc[diff, 'Name']):
                            # if it's the final oLCA name, we pop the original key
                            for k in list(iw_olca_dict.keys()):
                                if k[0] == diff:
                                    new_key = (olca.loc[diff, 'Name'].iloc[i],) + k[1:]
                                    iw_olca_dict[new_key] = iw_olca_dict.pop(k)
                        else:
                            for k in list(iw_olca_dict.keys()):
                                # otherwise we just add a new key to the dict
                                if k[0] == diff:
                                    new_key = (olca.loc[diff, 'Name'].iloc[i],) + k[1:]
                                    iw_olca_dict[new_key] = iw_olca_dict[k]

            # create the dataframe from the dictionary
            db = pd.DataFrame.from_dict(iw_olca_dict, orient='index')
            db.index = pd.MultiIndex.from_tuples(db.index)
            db = db.reset_index()
            db.columns = ['Elem flow name', 'Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                    'CAS number', 'Elem flow unit', 'MP or Damage', 'Native geographical resolution scale',
                                    'CF value']

            # ------------------------------------ UNITS --------------------------------

            # Unit name changes and others
            db.loc[db.loc[:, 'Elem flow unit'] == 'Bq', 'CF value'] *= 1000
            db.loc[db.loc[:, 'Elem flow unit'] == 'Bq', 'Elem flow unit'] = 'kBq'
            db.loc[db.loc[:, 'Elem flow unit'] == 'm2.yr', 'Elem flow unit'] = 'm2*a'
            db.loc[db.loc[:, 'Elem flow unit'] == 'kgy', 'Elem flow unit'] = 'kg*a'

            mj_flows = olca.loc[olca.loc[:, 'Name'].str.contains('MJ'), 'Name']
            for flow in mj_flows:
                db.loc[db.loc[:, 'Elem flow name'] == flow, 'CF value'] = float(
                    flow.split(' MJ')[0].split(', ')[-1])

            gj_flows = olca.loc[olca.loc[:, 'Name'].str.contains('GJ'), 'Name']
            for flow in gj_flows:
                if 'IN-GJ' not in flow:
                    db.loc[[i for i in db.index if db.loc[i, 'Elem flow name'] == flow and
                            db.loc[i, 'Impact category'] == 'Fossil and nuclear energy use'], 'CF value'] = float(
                        flow.split(' GJ')[0].split(', ')[-1]) * 1000

            # some other flows from oLCA that require conversions because of units
            new_flows = {
                'Gas, natural/kg': 'Gas, natural/m3',
                'Gas, mine, off-gas, process, coal mining/kg': 'Gas, mine, off-gas, process, coal mining/m3'
            }
            for flow in new_flows:
                if 'Gas' in flow:
                    density = 0.8  # kg/m3 (https://www.engineeringtoolbox.com/gas-density-d_158.html)

                df = db[db['Elem flow name'] == new_flows[flow]].copy()
                df.loc[:, 'Elem flow name'] = flow
                df.loc[:, 'Elem flow unit'] = 'kg'
                df.loc[:, 'CF value'] /= density
                db = pd.concat([db, df])
            db = clean_up_dataframe(db)

            # some water flows are in kilograms instead of cubic meters...
            problems = ['Water/kg',
                        'Water (fresh water)',
                        'Water (river water from technosphere turbined)',
                        'Water (river water from technosphere, turbined)',
                        'Water (with river silt)',
                        'Water, with river silt',
                        'Water, barrage',
                        'Water, cooling, drinking',
                        'Water, cooling, salt, ocean',
                        'Water, cooling, surface',
                        'Water, cooling, well',
                        'Water, cooling, well, in ground',
                        'Water, ground',
                        'Water, groundwater consumption',
                        'Water, process, drinking',
                        'Water, process, salt, ocean',
                        'Water, process, surface',
                        'Water, process, well',
                        'Water, Saline water consumption',
                        'Water, surface',
                        'Water, surface water consumption',
                        'Water, Surface water consumption',
                        'Water, thermoelectric groundwater consumption',
                        'Water, thermoelectric saline water consumption',
                        'Water, thermoelectric surface water consumption',
                        'Water, cooling, unspecified natural origin/kg',
                        'Water, process, unspecified natural origin/kg',
                        'Water, unspecified natural origin/kg',
                        'Thermally polluted water'
                        ]

            for problem_child in problems:
                db.loc[db['Elem flow name'] == problem_child, 'Elem flow unit'] = 'kg'
                db.loc[db['Elem flow name'] == problem_child, 'CF value'] /= 1000

            # specific problem for "Water" which is both defined in kg and in m3
            df = db.loc[db.loc[:, 'Elem flow name'] == 'Water'].copy()
            df.loc[:, 'Elem flow unit'] = 'kg'
            df.loc[:, 'CF value'] /= 1000
            db = clean_up_dataframe(pd.concat([db, df]))

            # --------------------------- COMPS AND SUBCOMPS --------------------------------
            with open(pkg_resources.resource_filename(__name__, '/Data/mappings/oLCA/v2.1.1/comps.json'), 'r') as f:
                comps = json.load(f)
            db.Compartment = [{v: k for k, v in comps.items()}[i] for i in db.Compartment]

            db.loc[db.loc[:, 'Sub-compartment'] == '(unspecified)', 'Sub-compartment'] = 'unspecified'
            db.loc[db.loc[:, 'Sub-compartment'] == 'groundwater', 'Sub-compartment'] = 'ground water'
            db.loc[db.loc[:, 'Sub-compartment'] == 'groundwater, long-term', 'Sub-compartment'] = 'ground water, long-term'
            db.loc[db.loc[:, 'Sub-compartment'] == 'high. pop.', 'Sub-compartment'] = 'high population density'
            db.loc[db.loc[:, 'Sub-compartment'] == 'low. pop.', 'Sub-compartment'] = 'low population density'
            db.loc[db.loc[:, 'Sub-compartment'] == 'low. pop., long-term',
                             'Sub-compartment'] = 'low population density, long-term'
            db.loc[db.loc[:, 'Sub-compartment'] == 'stratosphere + troposphere',
                             'Sub-compartment'] = 'lower stratosphere + upper troposphere'

            df = db.loc[db.loc[:, 'Sub-compartment'] == 'lake'].copy()
            df.loc[:, 'Sub-compartment'] = 'surface water'
            db = clean_up_dataframe(pd.concat([db, df]))

            # fix comp/subcomp of pesky biogenic carbon elementary flows
            db = db.drop(
                db.loc[db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                         'Carbon dioxide, in air'])].loc[
                    db.loc[:, 'Sub-compartment'] != 'unspecified'].index)
            db.loc[db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                     'Carbon dioxide, in air']), 'Compartment'] = 'Resource'
            db.loc[db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                     'Carbon dioxide, in air']), 'Sub-compartment'] = 'in air'
            db = db.drop(db.loc[db.loc[:, 'Elem flow name'] == 'Carbon dioxide, biogenic, uptake'].index)

            # ------------------------------ SUB-CATEGORIES SHENANIGANS ----------------------------

            db.loc[db['Elem flow name'].isin(['Carbon dioxide, land transformation',
                                              'Carbon monoxide, land transformation',
                                              'Methane, land transformation',
                                              'Carbon dioxide, from soil or biomass stock',
                                              'Carbon monoxide, from soil or biomass stock',
                                              'Methane, from soil or biomass stock']) &
                db['Impact category'].str.contains(', fossil', na=False), 'Impact category'] = [
                i.replace(', fossil', ', land transformation') for i in db.loc[
                    db['Elem flow name'].isin(['Carbon dioxide, land transformation',
                                               'Carbon monoxide, land transformation',
                                               'Methane, land transformation',
                                               'Carbon dioxide, from soil or biomass stock',
                                               'Carbon monoxide, from soil or biomass stock',
                                               'Methane, from soil or biomass stock']) &
                    db['Impact category'].str.contains(', fossil', na=False), 'Impact category']]

            # --------------------------- ADD OLCA UUIDS ------------------------------------
            olca_flows = pd.read_excel(pkg_resources.resource_filename(
                __name__, '/Data/mappings/oLCA/v2.1.1/all_stressors.xlsx'))

            # split comps and subcomps in two columns for matching with db
            olca_flows['Compartment'] = [i.split('/')[1] for i in olca_flows['comp']]
            olca_flows['Sub-compartment'] = [i.split('/')[2] for i in olca_flows['comp']]
            # only keep relevant columns
            olca_flows = olca_flows.loc[:, ['flow_id', 'flow_name', 'unit', 'Compartment', 'Sub-compartment']]
            # merge with olca_iw, it basically adds the uuids of oLCA
            db = olca_flows.merge(db, left_on=['flow_name', 'unit', 'Compartment', 'Sub-compartment'],
                                  right_on=['Elem flow name', 'Elem flow unit', 'Compartment', 'Sub-compartment'],
                                  how='left')
            # remove flows of IW+ with no link to oLCA flows
            db = db.drop(db.loc[db.loc[:, 'CF value'].isna()].index)
            # remove irrelevant columns
            db = db.drop(['flow_name', 'unit'], axis=1)

            return db

        self.olca_iw = linking(self.master_db)
        self.olca_iw_carbon_neutrality = linking(self.master_db_carbon_neutrality)

    def link_to_exiobase(self):
        """
        This method creates an openLCA method with the IW+ characterization factors.
        :return:
        """

        EXIO_IW_concordance = pd.read_excel(pkg_resources.resource_filename(
            __name__, 'Data/mappings/exiobase/EXIO_IW_concordance.xlsx'))
        EXIO_IW_concordance.set_index('EXIOBASE', inplace=True)

        self.exio_iw = pd.DataFrame(0, EXIO_IW_concordance.index, list(set(list(zip(self.master_db_carbon_neutrality.loc[:, 'Impact category'],
                                                                         self.master_db_carbon_neutrality.loc[:, 'CF unit'])))))
        self.exio_iw.columns = pd.MultiIndex.from_tuples(self.exio_iw.columns, names=['Impact category', 'CF unit'])
        self.exio_iw = self.exio_iw.T.sort_index().T

        for flow in EXIO_IW_concordance.index:
            if not EXIO_IW_concordance.loc[flow].isna().iloc[0]:
                # identify all entries (any impact category, compartment, etc.) for given flow
                CF_flow = self.master_db_carbon_neutrality.loc[self.master_db_carbon_neutrality['Elem flow name'] ==
                                                               EXIO_IW_concordance.loc[flow, 'IW']].loc[:,
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
        # cannot use +/-1 carbon biogenic approach with IO, so put back CFs from carbon neutrality approach
        self.exio_iw.loc[:, 'CO2 - waste - biogenic - air'] = 0

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

        df = self.master_db_carbon_neutrality.copy()
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
            CF_flow = self.master_db_carbon_neutrality.loc[self.master_db_carbon_neutrality['Elem flow name'] ==
                                                           other_metal_concordance.loc[flow].iloc[0]].loc[:,
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
        self.exio_iw.loc['Mineral resources use (kg deprived)',
                          'Unused Domestic Extraction - Non-Metallic Minerals - Other minerals'] = new_CF * 1000000

    def get_simplified_versions(self):

        # SimaPro
        self.simplified_version_sp = clean_up_dataframe(produce_simplified_version(self.iw_sp_carbon_neutrality).reindex(
            self.iw_sp_carbon_neutrality.columns, axis=1))

        # openLCA
        self.simplified_version_olca = clean_up_dataframe(produce_simplified_version(self.olca_iw_carbon_neutrality).reindex(
            self.olca_iw_carbon_neutrality.columns, axis=1))

        # ecoinvent
        self.simplified_version_ei38 = clean_up_dataframe(produce_simplified_version(self.ei38_iw_carbon_neutrality).reindex(
            self.ei38_iw_carbon_neutrality.columns, axis=1))

        self.simplified_version_ei39 = clean_up_dataframe(produce_simplified_version(self.ei39_iw_carbon_neutrality).reindex(
            self.ei39_iw_carbon_neutrality.columns, axis=1))

        self.simplified_version_ei310 = clean_up_dataframe(produce_simplified_version(self.ei310_iw_carbon_neutrality).reindex(
            self.ei310_iw_carbon_neutrality.columns, axis=1))

        self.simplified_version_ei311 = clean_up_dataframe(produce_simplified_version(self.ei311_iw_carbon_neutrality).reindex(
            self.ei311_iw_carbon_neutrality.columns, axis=1))

    def get_total_hh_and_eq(self):
        """
        OpenLCA doesn't allow for reliable contribution analyses for total damage categories (unlike
        SimaPro). So we create two additional impact categories "Total human health" and "Total ecosystem quality".
        :return:
        """

        total_hh = self.olca_iw.loc[self.olca_iw.loc[:, 'CF unit'] == 'DALY'].drop('Impact category', axis=1).groupby(
            by=['Elem flow name', 'Compartment', 'Sub-compartment', 'Elem flow unit', 'CF unit',
                'MP or Damage', 'flow_id']).agg({
            'CF value': sum,
            'CAS number': 'first'
        }).reset_index()
        total_hh.loc[:, 'Impact category'] = 'Total human health'

        total_eq = self.olca_iw.loc[self.olca_iw.loc[:, 'CF unit'] == 'PDF.m2.yr'].drop('Impact category',
                                                                                        axis=1).groupby(
            by=['Elem flow name', 'Compartment', 'Sub-compartment', 'Elem flow unit', 'CF unit',
                'MP or Damage', 'flow_id']).agg({
            'CF value': sum,
            'CAS number': 'first'
        }).reset_index()
        total_eq.loc[:, 'Impact category'] = 'Total ecosystem quality'

        self.olca_iw = clean_up_dataframe(pd.concat([self.olca_iw, total_hh, total_eq]))

        total_hh = self.olca_iw_carbon_neutrality.loc[self.olca_iw_carbon_neutrality.loc[:, 'CF unit'] == 'DALY'].drop(
            'Impact category', axis=1).groupby(by=['Elem flow name', 'Compartment', 'Sub-compartment',
                                                   'Elem flow unit', 'CF unit', 'MP or Damage', 'flow_id']).agg({
            'CF value': sum,
            'CAS number': 'first'
        }).reset_index()
        total_hh.loc[:, 'Impact category'] = 'Total human health'

        total_eq = self.olca_iw_carbon_neutrality.loc[self.olca_iw_carbon_neutrality.loc[:, 'CF unit'] == 'PDF.m2.yr'].drop(
            'Impact category', axis=1).groupby(by=['Elem flow name', 'Compartment', 'Sub-compartment',
                                                   'Elem flow unit', 'CF unit', 'MP or Damage', 'flow_id']).agg({
            'CF value': sum,
            'CAS number': 'first'
        }).reset_index()
        total_eq.loc[:, 'Impact category'] = 'Total ecosystem quality'

        self.olca_iw_carbon_neutrality = clean_up_dataframe(pd.concat([self.olca_iw_carbon_neutrality, total_hh, total_eq]))

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
    Method producing the simplified version of IW+ in which there are only 5 indicators:
    - Climate change, short term (=GWP100)
    - Water scarcity
    - Fossil resources
    - Remaining Human health damage
    - Remaining Ecosystem quality damage
    The damage impacts of climate change and water use on the Human health and Ecosystem quality area of protections
    are removed, hence only leaving "Remaining [...] damage". These are removed to avoid potential misuse by users
    where the impact of climate change and water use would be double-counted. Once as a midpoint category, and
    another time as a damage category. The simplified version relies on the carbon neutrality assumption.
    :return: a simplified version of the full "expert" database version
    """

    simplified_version = complete_dataframe.copy('deep')

    # midpoint categories excluded from simplified version
    midpoint_drop = ['Climate change, long term', 'Freshwater acidification', 'Freshwater ecotoxicity',
                     'Freshwater eutrophication', 'Human toxicity cancer', 'Human toxicity non-cancer',
                     'Ionizing radiations', 'Land occupation, biodiversity', 'Land transformation, biodiversity',
                     'Marine eutrophication', 'Mineral resources use', 'Ozone layer depletion',
                     'Particulate matter formation', 'Photochemical ozone formation',
                     'Plastics physical effects on biota', 'Terrestrial acidification']
    # endpoint categories excluded from simplified version
    endpoint_drop = ['Climate change, ecosystem quality, long term',
                     'Climate change, ecosystem quality, short term',
                     'Climate change, human health, long term', 'Climate change, human health, short term',
                     'Water availability, freshwater ecosystem', 'Water availability, human health',
                     'Water availability, terrestrial ecosystem', 'Marine ecotoxicity, long term',
                     'Human toxicity cancer, long term', 'Human toxicity non-cancer, long term',
                     'Freshwater ecotoxicity, long term', 'Marine acidification, long term',
                     'Terrestrial ecotoxicity, long term', 'Total human health', 'Total ecosystem quality']

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
    if 'flow_id' in simplified_version.columns: # for openLCA
        simplified_version = simplified_version.set_index(['CF unit', 'Compartment', 'Sub-compartment',
                                                       'Elem flow name', 'Elem flow unit', 'MP or Damage', 'flow_id'])
    else: # for other software
        simplified_version = simplified_version.set_index(['CF unit', 'Compartment', 'Sub-compartment',
                                                       'Elem flow name', 'Elem flow unit', 'MP or Damage'])
    # isolate and group HH CFs
    hh_simplified = simplified_version.loc['DALY'].copy()
    hh_simplified = hh_simplified.drop(['Impact category', 'CAS number',
                                        'Native geographical resolution scale'], axis=1).groupby(hh_simplified.index).sum()
    hh_simplified.index = pd.MultiIndex.from_tuples(hh_simplified.index)
    # isolate and group EQ CFs
    eq_simplified = simplified_version.loc['PDF.m2.yr'].copy()
    eq_simplified = eq_simplified.drop(['Impact category', 'CAS number',
                                        'Native geographical resolution scale'], axis=1).groupby(eq_simplified.index).sum()
    eq_simplified.index = pd.MultiIndex.from_tuples(eq_simplified.index)
    # delete HH and EQ CFs from original df
    simplified_version.drop(['DALY', 'PDF.m2.yr'], inplace=True)
    simplified_version = simplified_version.reset_index()
    # make hh_simplified respect the format of self.simplified_version_sp for concatenation
    hh_simplified = hh_simplified.reset_index()
    if len(hh_simplified.columns) == 7:  # for openLCA
        hh_simplified = hh_simplified.rename(
            columns={'level_0': 'Compartment', 'level_1': 'Sub-compartment', 'level_2': 'Elem flow name',
                     'level_3': 'Elem flow unit', 'level_4': 'MP or Damage', 'level_5': 'flow_id'})
    else:
        hh_simplified = hh_simplified.rename(
            columns={'level_0': 'Compartment', 'level_1': 'Sub-compartment', 'level_2': 'Elem flow name',
                     'level_3': 'Elem flow unit', 'level_4': 'MP or Damage'})
    hh_simplified.loc[:, 'CF unit'] = 'DALY'
    hh_simplified.loc[:, 'Impact category'] = 'Human health (residual)'
    # make eq_simplified respect the format of self.simplified_version_sp for concatenation
    eq_simplified = eq_simplified.reset_index()
    if len(eq_simplified.columns) == 7:  # for openLCA
        eq_simplified = eq_simplified.rename(
            columns={'level_0': 'Compartment', 'level_1': 'Sub-compartment', 'level_2': 'Elem flow name',
                     'level_3': 'Elem flow unit', 'level_4': 'MP or Damage', 'level_5': 'flow_id'})
    else:
        eq_simplified = eq_simplified.rename(
            columns={'level_0': 'Compartment', 'level_1': 'Sub-compartment', 'level_2': 'Elem flow name',
                     'level_3': 'Elem flow unit', 'level_4': 'MP or Damage'})
    eq_simplified.loc[:, 'CF unit'] = 'PDF.m2.yr'
    eq_simplified.loc[:, 'Impact category'] = 'Ecosystem quality (residual)'
    # concat
    simplified_version = pd.concat([simplified_version, hh_simplified, eq_simplified])
    # put back the CAS numbers
    simplified_version['CAS number'] = [cas[i] for i in simplified_version['Elem flow name']]

    simplified_version = clean_up_dataframe(simplified_version)

    simplified_version.loc[[i for i in simplified_version.index if simplified_version.loc[i, 'Impact category'] ==
                            'Climate change, short term'], 'Impact category'] = 'Carbon footprint'
    simplified_version.loc[[i for i in simplified_version.index if simplified_version.loc[i, 'Impact category'] ==
                            'Water scarcity'], 'Impact category'] = 'Water footprint - Scarcity'
    simplified_version.loc[[i for i in simplified_version.index if simplified_version.loc[i, 'Impact category'] ==
                            'Fossil and nuclear energy use'], 'Impact category'] = 'Energetic resource depletion'

    return simplified_version


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
