"""
This project aims to parse the IMPACT World+ files from the original Microsoft access database to the different
implementation in the available LCA software (SimaPro, openLCA & brightway2). A developer version of IW+ is also parsed
and is made available in an Excel format. This version is aimed for entities wishing to integrated IW+ in their
tools/databases on their own.

In the end, IMPACT World+ will span over 6 different files, in different formats:
- the core SQL database (referred to as the source version)
- the Excel version of IW+ regrouping characterization factors from the source version as well as additional
extrapolated CFs (referred to as the dev version)
- the IW+ version directly implementable in the SimaPro LCA software, in a .csv format (referred to as the SimaPro version)
- the IW+ version directly implementable in the openLCA LCA software, in a .zip format (referred to as the openLCA version)
- the IW+ version directly implementable in the brightway2 LCA software, in a .bw2package format (referred to as the
brightway2 version)
- the IW+ version directly usable with the EXIOBASE MRIO database, in a .xlsx format (referred to as the exiobase version)

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
import bw2data as bd
import bw2io as bi
from datetime import datetime
import csv
import warnings
import uuid
import logging
import sqlite3
import math
import molmass
import olca_ipc as ipc
import olca_schema as schema
from scipy.stats import gmean
from tqdm import tqdm


class Parse:
    def __init__(self, path_access_db, version, bw2_projects, bw_version):
        """
        :param path_access_db: path to the Microsoft access database (source version)
        :param version: the version of IW+ to parse
        :param bw2_projects: the name of a brightway2 project in which the database "biosphere3" is available
        :param bw_version: the version of brightway used, can be '2' or '2.5'

        Object instance variables:
        -------------------------
            - master_db : the master dataframe where basic IW CFs are stored (what is used to produce the dev.xlsx file)
                          following a +/-1 approach for biogenic carbon
            - master_db_carbon_neutrality : the master dataframe where basic IW CFs are stored (what is used to produce the dev.xlsx file)
                                            following a 0/0 approach for biogenic carbon
            - ei310_iw  : the dataframe where IW CFs linked to ecoinvent v3.10 elementary flows are stored,
                          following a +/-1 approach for biogenic carbon
            - ei310_iw_carbon_neutrality  : the dataframe where IW CFs linked to ecoinvent v3.10 elementary flows are stored
                                            following a 0/0 approach for biogenic carbon
            - simplified_version_ei310  : the dataframe where IW CFs linked to ecoinvent v3.10 elementary flows are stored
                                        for the footprint version
            - ei311_iw  : the dataframe where IW CFs linked to ecoinvent v3.11 elementary flows are stored,
                          following a +/-1 approach for biogenic carbon
            - ei311_iw_carbon_neutrality  : the dataframe where IW CFs linked to ecoinvent v3.11 elementary flows are stored
                                            following a 0/0 approach for biogenic carbon
            - simplified_version_ei311  : the dataframe where IW CFs linked to ecoinvent v3.11 elementary flows are stored
                                        for the footprint version
            - iw_sp  : the dataframe where IW CFs linked to SimaPro elementary flows are stored,
                          following a +/-1 approach for biogenic carbon
            - iw_sp_carbon_neutrality  : the dataframe where IW CFs linked to SimaPro elementary flows are stored
                                            following a 0/0 approach for biogenic carbon
            - simplified_version_sp  : the dataframe where IW CFs linked to SimaPro elementary flows are stored
                                        for the footprint version
            - olca_iw  : the dataframe where IW CFs linked to openLCA elementary flows are stored,
                          following a +/-1 approach for biogenic carbon
            - olca_iw_carbon_neutrality  : the dataframe where IW CFs linked to openLCA elementary flows are stored
                                            following a 0/0 approach for biogenic carbon
            - simplified_version_olca  : the dataframe where IW CFs linked to openLCA elementary flows are stored
                                        for the footprint version
            - exio_iw_38    : the dataframe where IW CFs linked to EXIOBASE v3.8.2 and before elementary flows are stored
            - exio_iw_39    : the dataframe where IW CFs linked to EXIOBASE v3.9.0 and after elementary flows are stored

        Object insteance methods:
        -------------------------
            - load_cfs()
            - load_basic_cfs()
            - load_climate_change_cfs()
            - load_ozone_layer_depletion_cfs()
            - load_photochemical_ozone_formation()
            - load_freshwater_acidification_cfs()
            - load_terrestrial_acidification_cfs()
            - load_marine_eutrophication_cfs()
            - load_resources_services_loss_cfs()
            - load_freshwater_eutrophication_cfs()
            - load_land_use_cfs()
            - load_particulates_cfs ()
            - load_water_scarcity_cfs()
            - load_water_availability_hh_cfs()
            - load_water_availability_fw_cfs()
            - load_water_availability_terr_cfs()
            - load_thermally_polluted_water_cfs()
            - load_physical_effects_cfs()
            - load_fisheries_cfs()
            - harmonize_regionalized_substances()
            - apply_rules()
            - create_not_regio_flows()
            - create_regio_flows_for_not_regio_ic()
            - order_things_around()
            - deal_with_biogenic_carbon()
            - deal_with_temporary_storage_of_carbon()
            - separate_ghg_indicators()
            - separate_regio_cfs()
            - link_to_ecoinvent()
            - link_to_sp()
            - link_to_exiobase()
            - link_to_olca()
            - get_simplified_versions()
            - get_total_hh_and_eq()
            - export_to_bw()
            - export_to_sp()
            - export_to_olca()
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
        self.bw_version = bw_version

        # OUTPUTs
        self.master_db = pd.DataFrame()
        self.master_db_carbon_neutrality = pd.DataFrame()
        self.master_db_not_regio = pd.DataFrame()
        self.master_db_not_regio_carbon_neutrality = pd.DataFrame()
        self.ei310_iw = pd.DataFrame()
        self.ei310_iw_carbon_neutrality = pd.DataFrame()
        self.ei311_iw = pd.DataFrame()
        self.ei311_iw_carbon_neutrality = pd.DataFrame()
        self.ei312_iw = pd.DataFrame()
        self.ei312_iw_carbon_neutrality = pd.DataFrame()
        self.simplified_version_ei310 = pd.DataFrame()
        self.simplified_version_ei311 = pd.DataFrame()
        self.simplified_version_ei312 = pd.DataFrame()
        self.iw_sp = pd.DataFrame()
        self.iw_sp_carbon_neutrality = pd.DataFrame()
        self.simplified_version_sp = pd.DataFrame()
        self.simplified_version_olca = pd.DataFrame()
        self.simplified_version_bw = pd.DataFrame()
        self.sp_data = {}
        self.olca_iw = pd.DataFrame()
        self.olca_iw_carbon_neutrality = pd.DataFrame()
        self.exio_iw_38 = pd.DataFrame()
        self.exio_iw_39 = pd.DataFrame()

        self.conn = sqlite3.connect(self.path_access_db)

        # Open openLCA. Open a database. Go to Tools/Dev tools/IPC server. Create a server with the local port 8080
        # (comes by default). Activate the port (the green olay arrow)
        # then we simply connect to that port. Now Python and openLCA can communicate via the "client" we created
        self.olca_client = ipc.Client(8080)
        # check if the opening of the openLCA client is activated or not
        try:
            self.olca_client.get_descriptor(schema.Unit, name='kg')
        except:
            print("The openLCA client. If you need it, don't forget to activate it in openLCA!")

    # -------------------------------------------- Main methods -------------------------------------------------------

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
        self.logger.info("Loading resources services loss characterization factors...")
        self.load_resources_services_loss_cfs()
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
        self.logger.info("Loading physical effects on biota characterization factors...")
        self.load_physical_effects_cfs()
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
        self.get_total_hh_and_eq_for_olca()

    def export_to_bw(self):
        """
        This method creates a brightway2 or brightway2.5 method with the IW+ characterization factors.

        :return:
        """

        self.logger.info("Exporting to brightway2...")

        for project in self.bw2_projects:
            bd.projects.set_current(project)
            # for bw2
            if 'biosphere3' in bd.databases:
                biosphere_db_name = 'biosphere3'
            # for bw2.5
            else:
                biosphere_db_name = [i for i in bd.databases if 'biosphere' in i][0]
            bio = bd.Database(biosphere_db_name)
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

            if '3.10' in project:
                ei_in_bw_normal = self.ei310_iw.merge(bw_flows_with_codes)
                ei_in_bw_carbon_neutrality = self.ei310_iw_carbon_neutrality.merge(bw_flows_with_codes)
                ei_in_bw_simple = self.simplified_version_ei310.merge(bw_flows_with_codes)
            elif '3.11' in project:
                ei_in_bw_normal = self.ei311_iw.merge(bw_flows_with_codes)
                ei_in_bw_carbon_neutrality = self.ei311_iw_carbon_neutrality.merge(bw_flows_with_codes)
                ei_in_bw_simple = self.simplified_version_ei311.merge(bw_flows_with_codes)
            elif '3.12' in project:
                ei_in_bw_normal = self.ei312_iw.merge(bw_flows_with_codes)
                ei_in_bw_carbon_neutrality = self.ei312_iw_carbon_neutrality.merge(bw_flows_with_codes)
                ei_in_bw_simple = self.simplified_version_ei312.merge(bw_flows_with_codes)
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
                total_eq.index = pd.MultiIndex.from_product(
                    [['Total ecosystem quality'], ['PDF.m2.yr'], total_eq.index])
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
                    new_method = bd.Method(name)
                    # register the new method
                    new_method.register()
                    # set its unit
                    new_method.metadata["unit"] = ic[1]

                    df = ei_in_bw.loc[[ic], ['code', 'CF value']].copy()
                    df.set_index('code', inplace=True)

                    data = []
                    for stressor in df.index:
                        data.append(((biosphere_db_name, stressor), df.loc[stressor, 'CF value']))
                    new_method.write(data)

            # -------------- For simplified version of IW+ ----------------

            ei_in_bw_simple.set_index(['Impact category', 'CF unit'], inplace=True)
            impact_categories_simple = ei_in_bw_simple.index.drop_duplicates()

            for ic in impact_categories_simple:

                name = ('IMPACT World+ Footprint ' + self.version + ' for ecoinvent v' + ei_version, ic[0])

                # initialize the "Method" method
                new_method = bd.Method(name)
                # register the new method
                new_method.register()
                # set its unit
                new_method.metadata["unit"] = ic[1]

                df = ei_in_bw_simple.loc[[ic], ['code', 'CF value']].copy()
                df.set_index('code', inplace=True)

                data = []
                for stressor in df.index:
                    data.append(((biosphere_db_name, stressor), df.loc[stressor, 'CF value']))
                new_method.write(data)

    def export_to_sp(self):
        """
        This method careates the necessary information for the csv creation in SimaPro.
        :return:
        """

        self.logger.info("Exporting to SimaPro...")

        # some elementary flows with the following region exceed the 100 characters limit of SimaPro -> drop region
        self.iw_sp_carbon_neutrality = self.iw_sp_carbon_neutrality.loc[
            ~self.iw_sp_carbon_neutrality.loc[:, 'Elem flow name'].str.contains(
                'United States of America, including overseas territories')]
        self.iw_sp = self.iw_sp.loc[~self.iw_sp.loc[:, 'Elem flow name'].str.contains(
            'United States of America, including overseas territories')]

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
                                    ['IMPACT World+ Midpoint ' + self.version + ' (incl. CO2 uptake)', '', '', '', '',
                                     ''],
                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                    ['2', '2', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                    ['IMPACT World+ ' + self.version + chr(int("007F", 16)) +
                                     'New categories:' + chr(int("007F", 16)) +
                                     '- Resources services loss' + chr(int("007F", 16)) +
                                     '- Resources services loss (adaptation)' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'Deleted categories:' + chr(int("007F", 16)) +
                                     '- Mineral resource use' + chr(int("007F", 16)) +
                                     '- Fossil and nuclear energy use' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'Renamed categories' + chr(int("007F", 16)) +
                                     '- Physical effect on biota (formerly Plastic physical effect on biota)' + chr(
                                        int("007F", 16)) +
                                     'Updated categories:' + chr(int("007F", 16)) +
                                     '- Acidification and eutrophication indicators' + chr(int("007F", 16)) +
                                     '- Physical effect on biota' + chr(int("007F", 16)) +
                                     '- Land occupation (biodiversity)' + chr(int("007F", 16)) +
                                     '- Land transformation (biodiversity)' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(
                                        int("007F", 16)) +
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
                                  ['2', '2', '', '', '', ''],
                                  ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                  ['IMPACT World+ ' + self.version + chr(int("007F", 16)) +
                                   'New categories:' + chr(int("007F", 16)) +
                                   '- Climate change, ecosystem quality, marine ecosystems (abbreviated as Climate change, EQ, mar)' + chr(
                                      int("007F", 16)) +
                                   '- Climate change, ecosystem quality, terrestrial ecosystems (abbreviated as Climate change, EQ, terr)' + chr(
                                      int("007F", 16))
                                   + chr(int("007F", 16)) +
                                   'Renamed categories' + chr(int("007F", 16)) +
                                   '- Physical effect on biota (formerly Plastic physical effect on biota)' + chr(
                                      int("007F", 16)) +
                                   'Updated categories:' + chr(int("007F", 16)) +
                                   '- Climate change damage indicators' + chr(int("007F", 16)) +
                                   '- Acidification and eutrophication indicators' + chr(
                                      int("007F", 16)) +
                                   '- Physical effect on biota' + chr(int("007F", 16)) +
                                   '- Land occupation (biodiversity)' + chr(int("007F", 16)) +
                                   '- Land transformation (biodiversity)' + chr(int("007F", 16)) +
                                   '- Marine ecotoxicity' + chr(int("007F", 16)) +
                                   '- Terrestrial ecotoxicity' + chr(int("007F", 16)) +
                                   chr(int("007F", 16)) +
                                   'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(
                                      int("007F", 16)) +
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
                                    ['2', '2', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                    ['IMPACT World+ ' + self.version + chr(int("007F", 16)) +
                                     'New categories:' + chr(int("007F", 16)) +
                                     '- Climate change, ecosystem quality, marine ecosystems (abbreviated as Climate change, EQ, mar)' + chr(
                                        int("007F", 16)) +
                                     '- Climate change, ecosystem quality, terrestrial ecosystems (abbreviated as Climate change, EQ, terr)' + chr(
                                        int("007F", 16)) +
                                     '- Resources services loss' + chr(int("007F", 16)) +
                                     '- Resources services loss (adaptation)' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'Deleted categories:' + chr(int("007F", 16)) +
                                     '- Mineral resource use' + chr(int("007F", 16)) +
                                     '- Fossil and nuclear energy use' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'Renamed categories' + chr(int("007F", 16)) +
                                     '- Physical effect on biota (formerly Plastic physical effect on biota)' + chr(
                                        int("007F", 16)) +
                                     'Updated categories:' + chr(int("007F", 16)) +
                                     '- Climate change damage indicators' + chr(int("007F", 16)) +
                                     '- Acidification and eutrophication indicators' + chr(int("007F", 16)) +
                                     '- Physical effect on biota' + chr(int("007F", 16)) +
                                     '- Land occupation (biodiversity)' + chr(int("007F", 16)) +
                                     '- Land transformation (biodiversity)' + chr(int("007F", 16)) +
                                     '- Marine ecotoxicity' + chr(int("007F", 16)) +
                                     '- Terrestrial ecotoxicity' + chr(int("007F", 16)) +
                                     chr(int("007F", 16)) +
                                     'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(
                                        int("007F", 16)) +
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
                                                    ['2', '2', '', '', '', ''],
                                                    ['IMPACT World+ ' + self.version + chr(int("007F", 16)) +
                                                     'New categories:' + chr(int("007F", 16)) +
                                                     '- Resources services loss' + chr(int("007F", 16)) +
                                                     '- Resources services loss (adaptation)' + chr(int("007F", 16)) +
                                                     chr(int("007F", 16)) +
                                                     'Deleted categories:' + chr(int("007F", 16)) +
                                                    '- Mineral resource use' + chr(int("007F", 16)) +
                                                    '- Fossil and nuclear energy use' + chr(int("007F", 16)) +
                                                     chr(int("007F", 16)) +
                                                    'Renamed categories' + chr(int("007F", 16)) +
                                                    '- Physical effect on biota (formerly Plastic physical effect on biota)' + chr(int("007F", 16)) +
                                                     'Updated categories:' + chr(int("007F", 16)) +
                                                     '- Acidification and eutrophication indicators' + chr(int("007F", 16)) +
                                                     '- Physical effect on biota' + chr(int("007F", 16)) +
                                                     '- Land occupation (biodiversity)' + chr(int("007F", 16)) +
                                                     '- Land transformation (biodiversity)' + chr(int("007F", 16)) +
                                                     chr(int("007F", 16)) +
                                                     'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(
                                                        int("007F", 16)) +
                                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes',
                                                     '', '', '', '', ''], ['', '', '', '', '', ''],
                                                    ['Category', '', '', '', '', ''], ['Others', '', '', '', '', ''],
                                                    ['', '', '', '', '', ''],
                                                    ['Use Damage Assessment', '', '', '', '', ''],
                                                    ['No', '', '', '', '', ''],
                                                    ['', '', '', '', '', ''],
                                                    ['Use Normalization', '', '', '', '', ''],
                                                    ['No', '', '', '', '', ''],
                                                    ['', '', '', '', '', ''],
                                                    ['Use Weighting', '', '', '', '', ''], ['No', '', '', '', '', ''],
                                                    ['', '', '', '', '', ''],
                                                    ['Use Addition', '', '', '', '', ''], ['No', '', '', '', '', '']]
        # metadata on the damage method
        damage_method_metadata_carboneutrality = [['Method', '', '', '', '', ''], ['', '', '', '', '', ''],
                                                  ['Name', '', '', '', '', ''],
                                                  ['IMPACT World+ Expert ' + self.version, '', '', '', '', ''],
                                                  ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                                  ['2', '2', '', '', '', ''],
                                                  ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                                  ['IMPACT World+ ' + self.version + chr(int("007F", 16)) +
                                                   'New categories:' + chr(int("007F", 16)) +
                                                   '- Climate change, ecosystem quality, marine ecosystems (abbreviated as Climate change, EQ, mar)' + chr(
                                                      int("007F", 16)) +
                                                   '- Climate change, ecosystem quality, terrestrial ecosystems (abbreviated as Climate change, EQ, terr)' + chr(
                                                      int("007F", 16))
                                                   + chr(int("007F", 16)) +
                                                   'Renamed categories' + chr(int("007F", 16)) +
                                                   '- Physical effect on biota (formerly Plastic physical effect on biota)' + chr(
                                                      int("007F", 16)) +
                                                   'Updated categories:' + chr(int("007F", 16)) +
                                                   '- Climate change damage indicators' + chr(int("007F", 16)) +
                                                   '- Acidification and eutrophication indicators' + chr(
                                                      int("007F", 16)) +
                                                   '- Physical effect on biota' + chr(int("007F", 16)) +
                                                   '- Land occupation (biodiversity)' + chr(int("007F", 16)) +
                                                   '- Land transformation (biodiversity)' + chr(int("007F", 16)) +
                                                   '- Marine ecotoxicity' + chr(int("007F", 16)) +
                                                   '- Terrestrial ecotoxicity' + chr(int("007F", 16)) +
                                                   chr(int("007F", 16)) +
                                                   'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(
                                                      int("007F", 16)) +
                                                   'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes',
                                                   '', '', '', '', ''], ['', '', '', '', '', ''],
                                                  ['Category', '', '', '', '', ''], ['Others', '', '', '', '', ''],
                                                  ['', '', '', '', '', ''],
                                                  ['Use Damage Assessment', '', '', '', '', ''],
                                                  ['Yes', '', '', '', '', ''],
                                                  ['', '', '', '', '', ''],
                                                  ['Use Normalization', '', '', '', '', ''],
                                                  ['Yes', '', '', '', '', ''],
                                                  ['', '', '', '', '', ''],
                                                  ['Use Weighting', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                                  ['', '', '', '', '', ''],
                                                  ['Use Addition', '', '', '', '', ''], ['Yes', '', '', '', '', '']]
        # metadata on the combined method
        combined_method_metadata_carboneutrality = [['Method', '', '', '', '', ''], ['', '', '', '', '', ''],
                                                    ['Name', '', '', '', '', ''],
                                                    ['IMPACT World+ ' + self.version, '', '', '', '', ''],
                                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                                    ['2', '2', '', '', '', ''],
                                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                                    ['IMPACT World+ ' + self.version + chr(int("007F", 16)) +
                                                     'New categories:' + chr(int("007F", 16)) +
                                                     '- Climate change, ecosystem quality, marine ecosystems (abbreviated as Climate change, EQ, mar)' + chr(int("007F", 16)) +
                                                     '- Climate change, ecosystem quality, terrestrial ecosystems (abbreviated as Climate change, EQ, terr)' + chr(int("007F", 16)) +
                                                     '- Resources services loss' + chr(int("007F", 16)) +
                                                     '- Resources services loss (adaptation)' + chr(int("007F", 16)) +
                                                     chr(int("007F", 16)) +
                                                     'Deleted categories:' + chr(int("007F", 16)) +
                                                    '- Mineral resource use' + chr(int("007F", 16)) +
                                                    '- Fossil and nuclear energy use' + chr(int("007F", 16)) +
                                                     chr(int("007F", 16)) +
                                                    'Renamed categories' + chr(int("007F", 16)) +
                                                    '- Physical effect on biota (formerly Plastic physical effect on biota)' + chr(int("007F", 16)) +
                                                     'Updated categories:' + chr(int("007F", 16)) +
                                                     '- Climate change damage indicators' + chr(int("007F", 16)) +
                                                     '- Acidification and eutrophication indicators' + chr(int("007F", 16)) +
                                                     '- Physical effect on biota' + chr(int("007F", 16)) +
                                                     '- Land occupation (biodiversity)' + chr(int("007F", 16)) +
                                                     '- Land transformation (biodiversity)' + chr(int("007F", 16)) +
                                                     '- Marine ecotoxicity' + chr(int("007F", 16)) +
                                                     '- Terrestrial ecotoxicity' + chr(int("007F", 16)) +
                                                     chr(int("007F", 16)) +
                                                     'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(
                                                        int("007F", 16)) +
                                                     'Full list of changes available here: https://github.com/CIRAIG/IWP_Reborn/tree/master/Report_changes',
                                                     '', '', '', '', ''], ['', '', '', '', '', ''],
                                                    ['Category', '', '', '', '', ''], ['Others', '', '', '', '', ''],
                                                    ['', '', '', '', '', ''],
                                                    ['Use Damage Assessment', '', '', '', '', ''],
                                                    ['Yes', '', '', '', '', ''],
                                                    ['', '', '', '', '', ''],
                                                    ['Use Normalization', '', '', '', '', ''],
                                                    ['Yes', '', '', '', '', ''],
                                                    ['', '', '', '', '', ''],
                                                    ['Use Weighting', '', '', '', '', ''], ['Yes', '', '', '', '', ''],
                                                    ['', '', '', '', '', ''],
                                                    ['Use Addition', '', '', '', '', ''], ['Yes', '', '', '', '', '']]
        # metadata on the simplified method
        simplified_method_metadata = [['Method', '', '', '', '', ''], ['', '', '', '', '', ''],
                                      ['Name', '', '', '', '', ''],
                                      ['IMPACT World+ Footprint ' + self.version, '', '', '', '', ''],
                                      ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                      ['2', '2', '', '', '', ''],
                                      ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                      ['IMPACT World+ Footprint ' + self.version + chr(int("007F", 16)) +
                                       'New categories:' + chr(int("007F", 16)) +
                                       '- Resources services loss' + chr(int("007F", 16)) +
                                       '- Resources services loss (adaptation)' + chr(int("007F", 16)) +
                                       chr(int("007F", 16)) +
                                       'Deleted categories:' + chr(int("007F", 16)) +
                                       '- Mineral resource use' + chr(int("007F", 16)) +
                                       '- Fossil and nuclear energy use' + chr(int("007F", 16)) +
                                       chr(int("007F", 16)) +
                                       'Updated categories:' + chr(int("007F", 16)) +
                                       '- Human health (residual)' + chr(int("007F", 16)) +
                                       '- Ecosystem quality (residual)' + chr(int("007F", 16)) +
                                       'For details on what the footprint version of IW+ entails, please consult this page: '
                                       'https://www.impactworldplus.org/version-2-0-1/' + chr(int("007F", 16)) +
                                       'For more information on IMPACT World+ and its methodology: https://www.impactworldplus.org.' + chr(
                                          int("007F", 16)) +
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
        weighting_info_damage_carboneutrality = [[df.loc[i].tolist()[0], df.loc[i].tolist()[1], '', '', '', ''] for i in
                                                 df.index]
        weighting_info_damage_carboneutrality[4] = ['Climate change, HH, LT', '1.00E+00', '', '', '', '']
        weighting_info_damage_carboneutrality[5] = ['Climate change, HH, ST', '1.00E+00', '', '', '', '']
        weighting_info_damage_carboneutrality[10] = ['Ionizing radiations, human health', '1.00E+00', '', '', '', '']
        weighting_info_damage_carboneutrality[13] = ['Photochemical ozone formation, human health', '1.00E+00', '', '',
                                                     '', '']
        weighting_info_damage_carboneutrality.insert(22, ['Fisheries impact', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality[27] = ['Ionizing radiations, ecosystem quality', '1.00E+00', '', '', '',
                                                     '']
        weighting_info_damage_carboneutrality.insert(32, ['Marine ecotoxicity, long term', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(33, ['Marine ecotoxicity, short term', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(35,
                                                     ['Photochemical ozone formation, ecosystem quality', '1.00E+00',
                                                      '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(36, ['Physical effects on biota', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(38,
                                                     ['Terrestrial ecotoxicity, long term', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(39, ['Terrestrial ecotoxicity, short term', '1.00E+00', '', '', '',
                                                          ''])
        weighting_info_damage_carboneutrality[20] = ['Climate change, EQ, terr, LT', '1.00E+00', '', '', '', '']
        weighting_info_damage_carboneutrality[21] = ['Climate change, EQ, terr, ST', '1.00E+00', '', '', '', '']
        weighting_info_damage_carboneutrality.insert(22,  ['Climate change, EQ, mar, LT', '1.00E+00', '', '', '', ''])
        weighting_info_damage_carboneutrality.insert(23, ['Climate change, EQ, mar, ST', '1.00E+00', '', '', '', ''])

        weighting_info_combined_carboneutrality = weighting_info_damage_carboneutrality.copy()
        weighting_info_combined_carboneutrality[11] = ['Ozone layer depletion (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[12] = ['Particulate matter formation (damage)', '1.00E+00', '', '', '',
                                                       '']
        weighting_info_combined_carboneutrality[25] = ['Freshwater acidification (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[28] = ['Freshwater eutrophication (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[30] = ['Land occupation, biodiversity (damage)', '1.00E+00', '', '', '',
                                                       '']
        weighting_info_combined_carboneutrality[31] = ['Land transformation, biodiversity (damage)', '1.00E+00', '', '',
                                                       '', '']
        weighting_info_combined_carboneutrality[36] = ['Marine eutrophication (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[38] = ['Physical effects on biota (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined_carboneutrality[39] = ['Terrestrial acidification (damage)', '1.00E+00', '', '', '', '']

        weighting_info_damage = [[df.loc[i].tolist()[0], df.loc[i].tolist()[1], '', '', '', ''] for i in df.index]
        weighting_info_damage[4] = ['Climate change, HH, LT, fossil', '1.00E+00', '', '', '', '']
        weighting_info_damage[5] = ['Climate change, HH, ST, fossil', '1.00E+00', '', '', '', '']
        weighting_info_damage[10] = ['Ionizing radiations, human health', '1.00E+00', '', '', '', '']
        weighting_info_damage[13] = ['Photochemical ozone formation, human health', '1.00E+00', '', '', '', '']
        weighting_info_damage[20] = ['Climate change, EQ, terr, LT, fossil', '1.00E+00', '', '', '', '']
        weighting_info_damage[21] = ['Climate change, EQ, terr, ST, fossil', '1.00E+00', '', '', '', '']
        weighting_info_damage.insert(22, ['Climate change, EQ, mar, LT, fossil', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(23, ['Climate change, EQ, mar, ST, fossil', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(24, ['Fisheries impact', '1.00E+00', '', '', '', ''])
        weighting_info_damage[29] = ['Ionizing radiations, ecosystem quality', '1.00E+00', '', '', '', '']
        weighting_info_damage.insert(34, ['Marine ecotoxicity, long term', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(35, ['Marine ecotoxicity, short term', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(37,
                                     ['Photochemical ozone formation, ecosystem quality', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(38, ['Physical effects on biota', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(40, ['Terrestrial ecotoxicity, long term', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(41, ['Terrestrial ecotoxicity, short term', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(6, ['Climate change, HH, LT, biogenic', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(7, ['Climate change, HH, ST, biogenic', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(8, ['Climate change, HH, LT, CO2 uptake', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(9, ['Climate change, HH, ST, CO2 uptake', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(10, ['Climate change, HH, LT, land transformation', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(11, ['Climate change, HH, ST, land transformation', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(28, ['Climate change, EQ, terr, LT, biogenic', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(29, ['Climate change, EQ, terr, ST, biogenic', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(30, ['Climate change, EQ, terr, LT, CO2 uptake', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(31, ['Climate change, EQ, terr, ST, CO2 uptake', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(32,
                                     ['Climate change, EQ, terr, LT, land transformation', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(33,
                                     ['Climate change, EQ, terr, ST, land transformation', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(36, ['Climate change, EQ, mar, LT, biogenic', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(37, ['Climate change, EQ, mar, ST, biogenic', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(38, ['Climate change, EQ, mar, LT, CO2 uptake', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(39, ['Climate change, EQ, mar, ST, CO2 uptake', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(40,
                                     ['Climate change, EQ, mar, LT, land transformation', '1.00E+00', '', '', '', ''])
        weighting_info_damage.insert(41,
                                     ['Climate change, EQ, mar, ST, land transformation', '1.00E+00', '', '', '', ''])

        weighting_info_combined = weighting_info_damage.copy()
        weighting_info_combined[17] = ['Ozone layer depletion (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[18] = ['Particulate matter formation (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[43] = ['Freshwater acidification (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[46] = ['Freshwater eutrophication (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[48] = ['Land occupation, biodiversity (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[49] = ['Land transformation, biodiversity (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[54] = ['Marine eutrophication (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[56] = ['Physical effects on biota (damage)', '1.00E+00', '', '', '', '']
        weighting_info_combined[57] = ['Terrestrial acidification (damage)', '1.00E+00', '', '', '', '']

        # extracting midpoint CFs
        d_ic_unit = self.iw_sp.loc[self.iw_sp['MP or Damage'] == 'Midpoint',
                                   ['Impact category', 'CF unit']].drop_duplicates().set_index('Impact category').iloc[
                    :,
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
                                                     ['Impact category', 'CF unit']].drop_duplicates().set_index(
            'Impact category').iloc[:,
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
                                   ['Impact category', 'CF unit']].drop_duplicates().set_index('Impact category').iloc[
                    :,
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
                                                     ['Impact category', 'CF unit']].drop_duplicates().set_index(
            'Impact category').iloc[:,
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
        same_names = ['Freshwater acidification', 'Freshwater eutrophication', 'Land occupation, biodiversity',
                      'Land transformation, biodiversity', 'Marine eutrophication', 'Ozone layer depletion',
                      'Particulate matter formation', 'Terrestrial acidification',
                      'Physical effects on biota']
        combined_values = []
        for j in ic_unit.index:
            combined_values.append(['', '', '', '', '', ''])
            combined_values.append(['Impact category', '', '', '', '', ''])
            if ic_unit.loc[j, 'Impact category'] in same_names:
                if ic_unit.loc[j, 'CF unit'] in ['DALY', 'PDF.m2.yr']:
                    combined_values.append([ic_unit.loc[j, 'Impact category'] + ' (damage)',
                                            ic_unit.loc[j, 'CF unit'], '', '', '', ''])
                else:
                    combined_values.append([ic_unit.loc[j, 'Impact category'] + ' (midpoint)',
                                            ic_unit.loc[j, 'CF unit'], '', '', '', ''])
            else:
                combined_values.append([ic_unit.loc[j, 'Impact category'],
                                        ic_unit.loc[j, 'CF unit'], '', '', '', ''])
            combined_values.append(['', '', '', '', '', ''])
            combined_values.append(['Substances', '', '', '', '', ''])
            df = self.iw_sp.loc[[i for i in self.iw_sp.index if (
                    self.iw_sp.loc[i, 'Impact category'] == ic_unit.loc[j, 'Impact category'] and
                    self.iw_sp.loc[i, 'CF unit'] == ic_unit.loc[j, 'CF unit'])]]
            df = df[['Compartment', 'Sub-compartment', 'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit']]
            for i in df.index:
                if type(df.loc[i, 'CAS number']) == float:
                    df.loc[i, 'CAS number'] = ''
                combined_values.append(df.loc[i].tolist())

        ic_unit = self.iw_sp_carbon_neutrality.loc[:, ['Impact category', 'CF unit']].drop_duplicates()
        same_names = ['Freshwater acidification', 'Freshwater eutrophication', 'Land occupation, biodiversity',
                      'Land transformation, biodiversity', 'Marine eutrophication', 'Ozone layer depletion',
                      'Particulate matter formation', 'Terrestrial acidification',
                      'Physical effects on biota']
        combined_values_carboneutrality = []
        for j in ic_unit.index:
            combined_values_carboneutrality.append(['', '', '', '', '', ''])
            combined_values_carboneutrality.append(['Impact category', '', '', '', '', ''])
            if ic_unit.loc[j, 'Impact category'] in same_names:
                if ic_unit.loc[j, 'CF unit'] in ['DALY', 'PDF.m2.yr']:
                    combined_values_carboneutrality.append([ic_unit.loc[j, 'Impact category'] + ' (damage)',
                                                            ic_unit.loc[j, 'CF unit'], '', '', '', ''])
                else:
                    combined_values_carboneutrality.append([ic_unit.loc[j, 'Impact category'] + ' (midpoint)',
                                                            ic_unit.loc[j, 'CF unit'], '', '', '', ''])
            else:
                combined_values_carboneutrality.append([ic_unit.loc[j, 'Impact category'],
                                                        ic_unit.loc[j, 'CF unit'], '', '', '', ''])
            combined_values_carboneutrality.append(['', '', '', '', '', ''])
            combined_values_carboneutrality.append(['Substances', '', '', '', '', ''])
            df = self.iw_sp_carbon_neutrality.loc[[i for i in self.iw_sp_carbon_neutrality.index if (
                    self.iw_sp_carbon_neutrality.loc[i, 'Impact category'] == ic_unit.loc[j, 'Impact category'] and
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

        locations = self.olca_client.get_all(schema.Location)
        locations = {i.code: i.to_ref() for i in locations}
        flows = self.olca_client.get_all(schema.Flow)
        flows = {i.id: i.to_ref() for i in flows}
        locations_undefined_in_olca = []

        def write_to_olca(olca_db, new_impact_method):
            for impact_category in tqdm(set(olca_db.loc[:, 'Impact category'])):
                # select impact category in IW+
                dff = olca_db.loc[olca_db.loc[:, 'Impact category'] == impact_category].copy()
                # create the impact category in oLCA
                new_impact_category = schema.new_impact_category(impact_category)
                new_impact_category.category = new_impact_method.category + '/' + new_impact_method.name
                new_impact_category.ref_unit = dff.loc[:, 'CF unit'].iloc[0]
                new_impact_category.impact_factors = []
                # add the category to the method
                new_impact_method.impact_categories.append(new_impact_category.to_ref())

                # loop through the IW+ CFs
                for ix in dff.index:
                    # for every row, it's a CF we need to create
                    new_impact_factor = schema.ImpactFactor()
                    # we find the flow "item" in oLCA
                    new_impact_factor.flow = flows[dff.loc[ix, 'flow_id']]
                    # dump in the CF value
                    new_impact_factor.value = dff.loc[ix, 'CF value']
                    # for regionalized impact categories, we need to add the location
                    if type(dff.loc[ix, 'Location']) == str:
                        if dff.loc[ix, 'Location'] in locations.keys():
                            # for the locations that exist in oLCA, add the location
                            new_impact_factor.location = locations[dff.loc[ix, 'Location']]
                            # finally we add the CF to the impact category
                            new_impact_category.impact_factors.append(new_impact_factor)
                            # if the location does not exist in oLCA, we do not create the CF
                        else:
                            locations_undefined_in_olca.append(dff.loc[ix, 'Location'])
                    else:
                        new_impact_category.impact_factors.append(new_impact_factor)
                self.olca_client.put(new_impact_category)

            self.olca_client.put(new_impact_method)

        df = self.olca_iw_carbon_neutrality.copy()
        same_names = ['Freshwater acidification', 'Freshwater eutrophication', 'Land occupation, biodiversity',
                      'Land transformation, biodiversity', 'Marine eutrophication', 'Ozone layer depletion',
                      'Particulate matter formation', 'Terrestrial acidification',
                      'Physical effects on biota']
        df.loc[(df.loc[:, 'Impact category'].isin(same_names) & (
                    df.loc[:, 'MP or Damage'] == 'Midpoint')), 'Impact category'] += ' (midpoint)'
        df.loc[(df.loc[:, 'Impact category'].isin(same_names) & (
                    df.loc[:, 'MP or Damage'] == 'Damage')), 'Impact category'] += ' (damage)'
        new_method = schema.ImpactMethod()
        new_method.category = 'IMPACT World+ LCIA Methods'
        new_method.version = '2.2'
        new_method.name = 'IMPACT World+ v2.2'
        new_method.description = 'For more information on IMPACT World+, check our website: https://www.impactworldplus.org/'
        new_method.impact_categories = []
        write_to_olca(df, new_method)

        df = self.olca_iw.loc[self.olca_iw.loc[:, 'MP or Damage'] == 'Damage'].copy()
        new_method = schema.ImpactMethod()
        new_method.category = 'IMPACT World+ LCIA Methods'
        new_method.version = '2.2'
        new_method.name = 'IMPACT World+ v2.2 - Expert - incl. CO2 uptake'
        new_method.description = 'For more information on IMPACT World+, check our website: https://www.impactworldplus.org/'
        new_method.impact_categories = []
        write_to_olca(df, new_method)

        df = self.olca_iw.loc[self.olca_iw.loc[:, 'MP or Damage'] == 'Midpoint'].copy()
        new_method = schema.ImpactMethod()
        new_method.category = 'IMPACT World+ LCIA Methods'
        new_method.version = '2.2'
        new_method.name = 'IMPACT World+ v2.2 - Midpoint - incl. CO2 uptake'
        new_method.description = 'For more information on IMPACT World+, check our website: https://www.impactworldplus.org/'
        new_method.impact_categories = []
        write_to_olca(df, new_method)

        df = self.olca_iw_carbon_neutrality.loc[
            self.olca_iw_carbon_neutrality.loc[:, 'MP or Damage'] == 'Damage'].copy()
        new_method = schema.ImpactMethod()
        new_method.category = 'IMPACT World+ LCIA Methods'
        new_method.version = '2.2'
        new_method.name = 'IMPACT World+ v2.2 - Expert'
        new_method.description = 'For more information on IMPACT World+, check our website: https://www.impactworldplus.org/'
        new_method.impact_categories = []
        write_to_olca(df, new_method)

        df = self.olca_iw_carbon_neutrality.loc[
            self.olca_iw_carbon_neutrality.loc[:, 'MP or Damage'] == 'Midpoint'].copy()
        new_method = schema.ImpactMethod()
        new_method.category = 'IMPACT World+ LCIA Methods'
        new_method.version = '2.2'
        new_method.name = 'IMPACT World+ v2.2 - Midpoint'
        new_method.description = 'For more information on IMPACT World+, check our website: https://www.impactworldplus.org/'
        new_method.impact_categories = []
        write_to_olca(df, new_method)

        df = self.simplified_version_olca.copy()
        new_method = schema.ImpactMethod()
        new_method.category = 'IMPACT World+ LCIA Methods'
        new_method.version = '2.2'
        new_method.name = 'IMPACT World+ v2.2 - Footprint'
        new_method.description = 'For more information on IMPACT World+, check our website: https://www.impactworldplus.org/'
        new_method.impact_categories = []
        write_to_olca(df, new_method)

    def produce_files(self):
        """
        Function producing the different IW+ files for the different versions.
        :return: the IW+ files
        """

        # Note, for openLCA we just export the files directly from the software

        self.logger.info("Creating all the files...")

        path = pkg_resources.resource_filename(__name__, '/Databases/Impact_world_' + self.version)

        # if the folders are not there yet, create them
        if not os.path.exists(path + '/Dev/'):
            os.makedirs(path + '/Dev/')
        if not os.path.exists(path + '/ecoinvent/'):
            os.makedirs(path + '/ecoinvent/')
        if not os.path.exists(path + '/exiobase/'):
            os.makedirs(path + '/exiobase/')
        if not os.path.exists(path + '/bw' + self.bw_version.replace('.', '')):
            os.makedirs(path + '/bw' + self.bw_version.replace('.', ''))
        if not os.path.exists(path + '/SimaPro/'):
            os.makedirs(path + '/SimaPro/')
        if not os.path.exists(path + '/openLCA/'):
            os.makedirs(path + '/openLCA/')

        # Dev version
        self.master_db.to_excel(path + '/Dev/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_dev.xlsx')
        self.master_db_carbon_neutrality.to_excel(path + '/Dev/impact_world_plus_' + self.version + '_dev.xlsx')

        # ecoinvent versions in Excel format
        self.ei310_iw.to_excel(
            path + '/ecoinvent/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_expert_version_ecoinvent_v310.xlsx')
        self.ei311_iw.to_excel(
            path + '/ecoinvent/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_expert_version_ecoinvent_v311.xlsx')
        self.ei312_iw.to_excel(
            path + '/ecoinvent/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_expert_version_ecoinvent_v312.xlsx')
        self.ei310_iw_carbon_neutrality.to_excel(
            path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v310.xlsx')
        self.ei311_iw_carbon_neutrality.to_excel(
            path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v311.xlsx')
        self.ei312_iw_carbon_neutrality.to_excel(
            path + '/ecoinvent/impact_world_plus_' + self.version + '_expert_version_ecoinvent_v312.xlsx')

        # ecoinvent version in DataFrame format
        self.simplified_version_ei310.to_excel(
            path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v310.xlsx')
        self.simplified_version_ei311.to_excel(
            path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v311.xlsx')
        self.simplified_version_ei312.to_excel(
            path + '/ecoinvent/impact_world_plus_' + self.version + '_footprint_version_ecoinvent_v312.xlsx')

        # exiobase version in DataFrame format
        self.exio_iw_38.to_excel(path + '/exiobase/impact_world_plus_' + self.version +
                                 '_expert_version_exiobase_3.8.2_and_before.xlsx')
        self.exio_iw_39.to_excel(path + '/exiobase/impact_world_plus_' + self.version +
                                 '_expert_version_exiobase_3.9_and_after.xlsx')

        # brightway2 versions in bw2package format
        for project in self.bw2_projects:
            bd.projects.set_current(project)
            if '3.10' in project:
                IW_ic = [bd.Method(ic) for ic in list(bd.methods) if
                         ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                          self.version in ic[0] and "for ecoinvent" in ic[0] and
                          ' (incl. CO2 uptake)' in ic[0] and 'regionalized' not in ic[0])]
                bi.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                  ' (incl. CO2 uptake)_brightway' + self.bw_version +
                                                                  '_expert_version_ei310',
                                                  folder=path + '/bw' + self.bw_version.replace('.', '') + '/')
                IW_ic = [bd.Method(ic) for ic in list(bd.methods) if
                         ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                          self.version in ic[0] and "for ecoinvent" in ic[0] and
                          ' (incl. CO2 uptake)' not in ic[0] and 'regionalized' not in ic[0])]
                bi.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                  '_brightway' + self.bw_version +
                                                                  '_expert_version_ei310',
                                                  folder=path + '/bw' + self.bw_version.replace('.', '') + '/')
                IW_ic = [bd.Method(ic) for ic in list(bd.methods) if
                         ('IMPACT World+' in ic[0] and 'Footprint' in ic[0] and
                          self.version in ic[0] and "for ecoinvent" in ic[0] and 'regionalized' not in ic[0])]
                bi.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                  '_brightway' + self.bw_version +
                                                                  '_footprint_version_ei310',
                                                  folder=path + '/bw' + self.bw_version.replace('.', '') + '/')
            elif '3.11' in project:
                IW_ic = [bd.Method(ic) for ic in list(bd.methods) if
                         ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                          self.version in ic[0] and "for ecoinvent" in ic[0] and
                          ' (incl. CO2 uptake)' in ic[0] and 'regionalized' not in ic[0])]
                bi.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                  ' (incl. CO2 uptake)_brightway' + self.bw_version +
                                                                  '_expert_version_ei311',
                                                  folder=path + '/bw' + self.bw_version.replace('.', '') + '/')
                IW_ic = [bd.Method(ic) for ic in list(bd.methods) if
                         ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                          self.version in ic[0] and "for ecoinvent" in ic[0] and
                          ' (incl. CO2 uptake)' not in ic[0] and 'regionalized' not in ic[0])]
                bi.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                  '_brightway' + self.bw_version +
                                                                  '_expert_version_ei311',
                                                  folder=path + '/bw' + self.bw_version.replace('.', '') + '/')
                IW_ic = [bd.Method(ic) for ic in list(bd.methods) if
                         ('IMPACT World+' in ic[0] and 'Footprint' in ic[0] and
                          self.version in ic[0] and "for ecoinvent" in ic[0] and 'regionalized' not in ic[0])]
                bi.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                  '_brightway' + self.bw_version +
                                                                  '_footprint_version_ei311',
                                                  folder=path + '/bw' + self.bw_version.replace('.', '') + '/')
            elif '3.12' in project:
                IW_ic = [bd.Method(ic) for ic in list(bd.methods) if
                         ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                          self.version in ic[0] and "for ecoinvent" in ic[0] and
                          ' (incl. CO2 uptake)' in ic[0] and 'regionalized' not in ic[0])]
                bi.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                  ' (incl. CO2 uptake)_brightway' + self.bw_version +
                                                                  '_expert_version_ei312',
                                                  folder=path + '/bw' + self.bw_version.replace('.', '') + '/')
                IW_ic = [bd.Method(ic) for ic in list(bd.methods) if
                         ('IMPACT World+' in ic[0] and 'Footprint' not in ic[0] and
                          self.version in ic[0] and "for ecoinvent" in ic[0] and
                          ' (incl. CO2 uptake)' not in ic[0] and 'regionalized' not in ic[0])]
                bi.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                  '_brightway' + self.bw_version +
                                                                  '_expert_version_ei312',
                                                  folder=path + '/bw' + self.bw_version.replace('.', '') + '/')
                IW_ic = [bd.Method(ic) for ic in list(bd.methods) if
                         ('IMPACT World+' in ic[0] and 'Footprint' in ic[0] and
                          self.version in ic[0] and "for ecoinvent" in ic[0] and 'regionalized' not in ic[0])]
                bi.package.BW2Package.export_objs(IW_ic, filename='impact_world_plus_' + self.version +
                                                                  '_brightway' + self.bw_version +
                                                                  '_footprint_version_ei312',
                                                  folder=path + '/bw' + self.bw_version.replace('.', '') + '/')

        # SimaPro version in csv format
        with open(
                path + '/SimaPro/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_midpoint_version_simapro.csv',
                'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['midpoint_method_metadata'] +
                self.sp_data['midpoint_values'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(
                path + '/SimaPro/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_expert_version_simapro.csv',
                'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['damage_method_metadata'] +
                self.sp_data['damage_values'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_damage'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path + '/SimaPro/impact_world_plus_' + self.version + ' (incl. CO2 uptake)_simapro.csv', 'w',
                  newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['combined_method_metadata'] +
                self.sp_data['combined_values'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_combined'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path + '/SimaPro/impact_world_plus_' + self.version + '_footprint_version_simapro.csv', 'w',
                  newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['simplified_method_metadata'] +
                self.sp_data['simplified_values'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path + '/SimaPro/impact_world_plus_' + self.version + '_midpoint_version_simapro.csv', 'w',
                  newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data[
                    'midpoint_method_metadata_carboneutrality'] +
                self.sp_data['midpoint_values_carboneutrality'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path + '/SimaPro/impact_world_plus_' + self.version + '_expert_version_simapro.csv', 'w',
                  newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data[
                    'damage_method_metadata_carboneutrality'] +
                self.sp_data['damage_values_carboneutrality'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_damage_carboneutrality'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path + '/SimaPro/impact_world_plus_' + self.version + '_simapro.csv', 'w',
                  newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data[
                    'combined_method_metadata_carboneutrality'] +
                self.sp_data['combined_values_carboneutrality'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_combined_carboneutrality'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])

    # ----------------------------------------- Secondary methods -----------------------------------------------------

    def load_basic_cfs(self):
        """
        Loading the basic CFs. By basic we mean that these CFs do not require further treatment.

        Concerned impact categories:
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

        :return: updated master_db
        """

        self.logger.info("Loading ionizing radiations characterization factors...")
        ionizing = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - IonizingRadiations]', con=self.conn)

        self.logger.info("Loading marine acidification characterization factors...")
        mar_acid = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - MarineAcidification]', con=self.conn)

        elem_flow_list = pd.read_sql(sql='SELECT * FROM [SI - Mapping with elementary flows]', con=self.conn)

        self.logger.info("Loading toxicity characterization factors...")
        toxicity = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - HumanTox]', con=self.conn)
        toxicity = toxicity.merge(elem_flow_list.loc[:, ['Name IW+', 'CAS-Usetox2_FW']], left_on=['CAS number'],
                                  right_on=['CAS-Usetox2_FW'], how='left').drop_duplicates()
        toxicity = toxicity.drop(['Elem flow name', 'CAS number'], axis=1)
        toxicity = toxicity.rename(columns={'Name IW+': 'Elem flow name', 'CAS-Usetox2_FW': 'CAS number'})

        self.logger.info("Loading freshwater ecotoxicity characterization factors...")
        fw_ecotoxicity = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - FreshwaterEcotox]', con=self.conn)
        fw_ecotoxicity = fw_ecotoxicity.merge(elem_flow_list.loc[:, ['Name IW+', 'CAS-Usetox2_FW']],
                                              left_on=['CAS number'],
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

        self.master_db = pd.concat([ionizing, mar_acid, toxicity, fw_ecotoxicity, mar_ecotoxicity, terr_ecotoxicity])

        self.master_db = clean_up_dataframe(self.master_db)

    def load_climate_change_cfs(self):
        """
        Loading the CFs for the climate change impact categories.

        Concerned impact categories:
            - Climate change, short term
            - Climate change, long term
            - Climate change, human health, short term
            - Climate change, human health, long term
            - Climate change, ecosystem quality, terrestrial ecosystem, short term
            - Climate change, ecosystem quality, terrestrial ecosystem, long term
            - Climate change, ecosystem quality, marine ecosystem, short term
            - Climate change, ecosystem quality, marine ecosystem, long term

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

        # get fate factors
        fate_factors = pd.read_sql('SELECT * FROM "SI - Climate change - fate factors (K/kg)"',
                                   con=self.conn).set_index(['Name IW+', 'CAS IW+'])

        # get effect factors
        effect_factors = pd.read_sql('SELECT * FROM "SI - Climate change - effect factors"', con=self.conn).set_index(
            'index')
        HH_effect_factor = effect_factors.loc['Total', 'Human health (DALY/K/yr)']
        EQ_terr_effect_factor = effect_factors.loc['Total', 'Ecosystem quality - terrestrial species (PDF.m2/K/yr)']
        EQ_mar_effect_factor = effect_factors.loc['Total', 'Ecosystem quality - marine species (PDF.m2/K/yr)']

        # Climate change, human health, short term
        GWP_damage_HH_short = (fate_factors.iloc[:, :101].sum(1) * HH_effect_factor).reset_index()
        GWP_damage_HH_short.columns = ['Elem flow name', 'CAS number', 'CF value']
        GWP_damage_HH_short.loc[:, 'Impact category'] = 'Climate change, human health, short term'
        GWP_damage_HH_short.loc[:, 'CF unit'] = 'DALY'
        GWP_damage_HH_short.loc[:, 'Compartment'] = 'Air'
        GWP_damage_HH_short.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_damage_HH_short.loc[:, 'Elem flow unit'] = 'kg'
        GWP_damage_HH_short.loc[:, 'MP or Damage'] = 'Damage'
        GWP_damage_HH_short.loc[:, 'Native geographical resolution scale'] = 'Global'

        # Climate change, human health, long term
        GWP_damage_HH_long = (fate_factors.iloc[:, 101:].sum(1) * HH_effect_factor).reset_index()
        GWP_damage_HH_long.columns = ['Elem flow name', 'CAS number', 'CF value']
        GWP_damage_HH_long.loc[:, 'Impact category'] = 'Climate change, human health, long term'
        GWP_damage_HH_long.loc[:, 'CF unit'] = 'DALY'
        GWP_damage_HH_long.loc[:, 'Compartment'] = 'Air'
        GWP_damage_HH_long.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_damage_HH_long.loc[:, 'Elem flow unit'] = 'kg'
        GWP_damage_HH_long.loc[:, 'MP or Damage'] = 'Damage'
        GWP_damage_HH_long.loc[:, 'Native geographical resolution scale'] = 'Global'

        # Climate change, ecosystem quality, terrestrial ecosystem, short term
        GWP_damage_EQ_terr_short = (fate_factors.iloc[:, :101].sum(1) * EQ_terr_effect_factor).reset_index()
        GWP_damage_EQ_terr_short.columns = ['Elem flow name', 'CAS number', 'CF value']
        GWP_damage_EQ_terr_short.loc[:, 'Impact category'] = 'Climate change, ecosystem quality, terrestrial ecosystem, short term'
        GWP_damage_EQ_terr_short.loc[:, 'CF unit'] = 'PDF.m2.yr'
        GWP_damage_EQ_terr_short.loc[:, 'Compartment'] = 'Air'
        GWP_damage_EQ_terr_short.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_damage_EQ_terr_short.loc[:, 'Elem flow unit'] = 'kg'
        GWP_damage_EQ_terr_short.loc[:, 'MP or Damage'] = 'Damage'
        GWP_damage_EQ_terr_short.loc[:, 'Native geographical resolution scale'] = 'Global'

        # Climate change, ecosystem quality, terrestrial ecosystem, long term
        GWP_damage_EQ_terr_long = (fate_factors.iloc[:, 101:].sum(1) * EQ_terr_effect_factor).reset_index()
        GWP_damage_EQ_terr_long.columns = ['Elem flow name', 'CAS number', 'CF value']
        GWP_damage_EQ_terr_long.loc[:, 'Impact category'] = 'Climate change, ecosystem quality, terrestrial ecosystem, long term'
        GWP_damage_EQ_terr_long.loc[:, 'CF unit'] = 'PDF.m2.yr'
        GWP_damage_EQ_terr_long.loc[:, 'Compartment'] = 'Air'
        GWP_damage_EQ_terr_long.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_damage_EQ_terr_long.loc[:, 'Elem flow unit'] = 'kg'
        GWP_damage_EQ_terr_long.loc[:, 'MP or Damage'] = 'Damage'
        GWP_damage_EQ_terr_long.loc[:, 'Native geographical resolution scale'] = 'Global'

        # Climate change, ecosystem quality, marine ecosystem, short term
        GWP_damage_EQ_mar_short = (fate_factors.iloc[:, :101].sum(1) * EQ_mar_effect_factor).reset_index()
        GWP_damage_EQ_mar_short.columns = ['Elem flow name', 'CAS number', 'CF value']
        GWP_damage_EQ_mar_short.loc[:, 'Impact category'] = 'Climate change, ecosystem quality, marine ecosystem, short term'
        GWP_damage_EQ_mar_short.loc[:, 'CF unit'] = 'PDF.m2.yr'
        GWP_damage_EQ_mar_short.loc[:, 'Compartment'] = 'Air'
        GWP_damage_EQ_mar_short.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_damage_EQ_mar_short.loc[:, 'Elem flow unit'] = 'kg'
        GWP_damage_EQ_mar_short.loc[:, 'MP or Damage'] = 'Damage'
        GWP_damage_EQ_mar_short.loc[:, 'Native geographical resolution scale'] = 'Global'

        # Climate change, ecosystem quality, marine ecosystem, long term
        GWP_damage_EQ_mar_long = (fate_factors.iloc[:, 101:].sum(1) * EQ_mar_effect_factor).reset_index()
        GWP_damage_EQ_mar_long.columns = ['Elem flow name', 'CAS number', 'CF value']
        GWP_damage_EQ_mar_long.loc[:, 'Impact category'] = 'Climate change, ecosystem quality, marine ecosystem, long term'
        GWP_damage_EQ_mar_long.loc[:, 'CF unit'] = 'PDF.m2.yr'
        GWP_damage_EQ_mar_long.loc[:, 'Compartment'] = 'Air'
        GWP_damage_EQ_mar_long.loc[:, 'Sub-compartment'] = '(unspecified)'
        GWP_damage_EQ_mar_long.loc[:, 'Elem flow unit'] = 'kg'
        GWP_damage_EQ_mar_long.loc[:, 'MP or Damage'] = 'Damage'
        GWP_damage_EQ_mar_long.loc[:, 'Native geographical resolution scale'] = 'Global'

        # deal with biogenic methane
        def add_biogenic_methane(df):
            dff = df.loc[df.loc[:, 'Elem flow name'] == 'Methane, fossil'].copy()
            dff.loc[:, 'Elem flow name'] = 'Methane, biogenic'
            dff.loc[:, 'CF value'] *= (
                        GWP_midpoint.loc[GWP_midpoint.loc[:, 'Elem flow name'] == 'Methane, biogenic', 'CF value'].iloc[0] /
                        GWP_midpoint.loc[GWP_midpoint.loc[:, 'Elem flow name'] == 'Methane, fossil', 'CF value'].iloc[0])
            df = pd.concat([df, dff])
            return df

        GWP_damage_HH_short = add_biogenic_methane(GWP_damage_HH_short)
        GWP_damage_HH_long = add_biogenic_methane(GWP_damage_HH_long)
        GWP_damage_EQ_terr_short = add_biogenic_methane(GWP_damage_EQ_terr_short)
        GWP_damage_EQ_terr_long = add_biogenic_methane(GWP_damage_EQ_terr_long)
        GWP_damage_EQ_mar_short = add_biogenic_methane(GWP_damage_EQ_mar_short)
        GWP_damage_EQ_mar_long = add_biogenic_methane(GWP_damage_EQ_mar_long)

        self.master_db = pd.concat([self.master_db, GWP_midpoint, GTP_midpoint,
                                    GWP_damage_EQ_terr_short, GWP_damage_EQ_terr_long,
                                    GWP_damage_EQ_mar_short, GWP_damage_EQ_mar_long,
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
        data = pd.read_sql('SELECT * FROM [CF - regionalized - AcidFW - aggregated]', self.conn)

        # concatenating/formatting all the data in one dataframe
        concat_data = pd.DataFrame()

        for ix, row in data.iterrows():
            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Freshwater acidification', 'kg SO2 eq', 'Air', '(unspecified)',
                                                   'Ammonia, ' + row['Short name ecoinvent'], '7664-41-7',
                                                   row['CF NH3 (kg SO2 eq)'], 'kg', 'Midpoint', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Freshwater acidification', 'PDF.m2.yr', 'Air', '(unspecified)',
                                                   'Ammonia, ' + row['Short name ecoinvent'], '7664-41-7',
                                                   row['CF NH3 (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Freshwater acidification', 'kg SO2 eq', 'Air', '(unspecified)',
                                                   'Nitrogen oxides, ' + row['Short name ecoinvent'], '11104-93-1',
                                                   row['CF NOx (kg SO2 eq)'], 'kg', 'Midpoint', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Freshwater acidification', 'PDF.m2.yr', 'Air', '(unspecified)',
                                                   'Nitrogen oxides, ' + row['Short name ecoinvent'], '11104-93-1',
                                                   row['CF NOx (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Freshwater acidification', 'kg SO2 eq', 'Air', '(unspecified)',
                                                   'Sulfur dioxide, ' + row['Short name ecoinvent'], '7446-09-05',
                                                   row['CF SO2 (kg SO2 eq)'], 'kg', 'Midpoint', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Freshwater acidification', 'PDF.m2.yr', 'Air', '(unspecified)',
                                                   'Sulfur dioxide, ' + row['Short name ecoinvent'], '7446-09-05',
                                                   row['CF SO2 (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Freshwater acidification', 'kg SO2 eq', 'Air', '(unspecified)',
                                                   'Nitric acid, ' + row['Short name ecoinvent'], '7697-37-2',
                                                   row['CF HNO3 (kg SO2 eq)'], 'kg', 'Midpoint', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Freshwater acidification', 'PDF.m2.yr', 'Air', '(unspecified)',
                                                   'Nitric acid, ' + row['Short name ecoinvent'], '7697-37-2',
                                                   row['CF HNO3 (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

        concat_data = concat_data.reset_index().drop('index', axis=1)

        concat_data.loc[[i for i in concat_data.index if
                         concat_data.loc[i, 'Elem flow name'].split(', ')[-1] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF',
                                                                                  'RME', 'UN-OCEANIA']],
                        'Native geographical resolution scale'] = 'Continent'
        concat_data.loc[[i for i in concat_data.index if concat_data.loc[i, 'Elem flow name'].split(', ')[
            -1] == 'GLO'], 'Native geographical resolution scale'] = 'Global'
        concat_data.loc[[i for i in concat_data.index if concat_data.loc[i, 'Elem flow name'].split(', ')[
            -1] == 'RoW'], 'Native geographical resolution scale'] = 'Other region'

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - Stoechiometry]', self.conn)
        stoc = stoc.loc[stoc.loc[:, 'Impact category'] == 'Freshwater acidification']

        for ix in stoc.index:
            proxy = stoc.loc[ix, 'Proxy molecule']
            comp = stoc.loc[ix, 'Compartment']
            df = pd.DataFrame()
            if proxy == 'NH3':
                # don't create a duplicate
                if stoc.loc[ix, 'Elem flow name'] != 'Ammonia':
                    df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Ammonia') &
                                     ~(concat_data.loc[:, 'Elem flow name'].str.contains('Ammonia, as N'))].loc[
                        concat_data.loc[:, 'Compartment'] == comp].copy('deep')
                    df.loc[:, 'Elem flow name'] = [i.replace('Ammonia', stoc.loc[ix, 'Elem flow name']) for i in
                                                   df.loc[:, 'Elem flow name']]
            if proxy == 'HNO3':
                # don't create a duplicate
                if stoc.loc[ix, 'Elem flow name'] != 'Nitric acid':
                    df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Nitric acid')].loc[
                        concat_data.loc[:, 'Compartment'] == comp].copy('deep')
                    df.loc[:, 'Elem flow name'] = [i.replace('Nitric acid', stoc.loc[ix, 'Elem flow name']) for i in
                                                   df.loc[:, 'Elem flow name']]
            if proxy == 'SO2':
                # don't create a duplicate
                if stoc.loc[ix, 'Elem flow name'] != 'Sulfur dioxide':
                    df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Sulfur dioxide')].loc[
                        concat_data.loc[:, 'Compartment'] == comp].copy('deep')
                    df.loc[:, 'Elem flow name'] = [i.replace('Sulfur dioxide', stoc.loc[ix, 'Elem flow name']) for i in
                                                   df.loc[:, 'Elem flow name']]
            if proxy == 'NOx':
                # don't create a duplicate
                if stoc.loc[ix, 'Elem flow name'] != 'Nitrogen oxides':
                    df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Nitrogen oxides')].loc[
                        concat_data.loc[:, 'Compartment'] == comp].copy('deep')
                    df.loc[:, 'Elem flow name'] = [i.replace('Nitrogen oxides', stoc.loc[ix, 'Elem flow name']) for i in
                                                   df.loc[:, 'Elem flow name']]
            if not df.empty:
                df.loc[:, 'CAS number'] = stoc.loc[ix, 'CAS number']
                df.loc[:, 'CF value'] *= stoc.loc[ix, 'Proxy ratio']

                concat_data = clean_up_dataframe(pd.concat([concat_data, df]))

        # concat with master_db
        self.master_db = pd.concat([self.master_db, concat_data])
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
        data = pd.read_sql('SELECT * FROM [CF - regionalized - AcidTerr - aggregated]', self.conn)

        # concatenating/formatting all the data in one dataframe
        concat_data = pd.DataFrame()

        for ix, row in data.iterrows():
            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Terrestrial acidification', 'kg SO2 eq', 'Air', '(unspecified)',
                                                   'Ammonia, ' + row['Short name ecoinvent'], '7664-41-7',
                                                   row['CF NH3 (kg SO2 eq)'], 'kg', 'Midpoint', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Terrestrial acidification', 'PDF.m2.yr', 'Air', '(unspecified)',
                                                   'Ammonia, ' + row['Short name ecoinvent'], '7664-41-7',
                                                   row['CF NH3 (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Terrestrial acidification', 'kg SO2 eq', 'Air', '(unspecified)',
                                                   'Nitrogen oxides, ' + row['Short name ecoinvent'], '11104-93-1',
                                                   row['CF NOx (kg SO2 eq)'], 'kg', 'Midpoint', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Terrestrial acidification', 'PDF.m2.yr', 'Air', '(unspecified)',
                                                   'Nitrogen oxides, ' + row['Short name ecoinvent'], '11104-93-1',
                                                   row['CF NOx (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Terrestrial acidification', 'kg SO2 eq', 'Air', '(unspecified)',
                                                   'Sulfur dioxide, ' + row['Short name ecoinvent'], '7446-09-05',
                                                   row['CF SO2 (kg SO2 eq)'], 'kg', 'Midpoint', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Terrestrial acidification', 'PDF.m2.yr', 'Air', '(unspecified)',
                                                   'Sulfur dioxide, ' + row['Short name ecoinvent'], '7446-09-05',
                                                   row['CF SO2 (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

        concat_data = concat_data.reset_index().drop('index', axis=1)
        concat_data.loc[[i for i in concat_data.index if
                         concat_data.loc[i, 'Elem flow name'].split(', ')[-1] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF',
                                                                                  'RME', 'UN-OCEANIA']],
                        'Native geographical resolution scale'] = 'Continent'
        concat_data.loc[[i for i in concat_data.index if concat_data.loc[i, 'Elem flow name'].split(', ')[
            -1] == 'GLO'], 'Native geographical resolution scale'] = 'Global'
        concat_data.loc[[i for i in concat_data.index if concat_data.loc[i, 'Elem flow name'].split(', ')[
            -1] == 'RoW'], 'Native geographical resolution scale'] = 'Other region'

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - Stoechiometry]', self.conn)
        stoc = stoc.loc[stoc.loc[:, 'Impact category'] == 'Terrestrial acidification']

        for ix in stoc.index:
            proxy = stoc.loc[ix, 'Proxy molecule']
            comp = stoc.loc[ix, 'Compartment']
            df = pd.DataFrame()
            if proxy == 'NH3':
                # don't create a duplicate
                if stoc.loc[ix, 'Elem flow name'] != 'Ammonia':
                    df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Ammonia') &
                                     ~(concat_data.loc[:, 'Elem flow name'].str.contains('Ammonia, as N'))].loc[
                        concat_data.loc[:, 'Compartment'] == comp].copy('deep')
                    df.loc[:, 'Elem flow name'] = [i.replace('Ammonia', stoc.loc[ix, 'Elem flow name']) for i in
                                                   df.loc[:, 'Elem flow name']]
            if proxy == 'SO2':
                # don't create a duplicate
                if stoc.loc[ix, 'Elem flow name'] != 'Sulfur dioxide':
                    df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Sulfur dioxide')].loc[
                        concat_data.loc[:, 'Compartment'] == comp].copy('deep')
                    df.loc[:, 'Elem flow name'] = [i.replace('Sulfur dioxide', stoc.loc[ix, 'Elem flow name']) for i in
                                                   df.loc[:, 'Elem flow name']]
            if proxy == 'NOx':
                # don't create a duplicate
                if stoc.loc[ix, 'Elem flow name'] != 'Nitrogen oxides':
                    df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Nitrogen oxides')].loc[
                        concat_data.loc[:, 'Compartment'] == comp].copy('deep')
                    df.loc[:, 'Elem flow name'] = [i.replace('Nitrogen oxides', stoc.loc[ix, 'Elem flow name']) for i in
                                                   df.loc[:, 'Elem flow name']]
            if not df.empty:
                df.loc[:, 'CAS number'] = stoc.loc[ix, 'CAS number']
                df.loc[:, 'CF value'] *= stoc.loc[ix, 'Proxy ratio']

                concat_data = clean_up_dataframe(pd.concat([concat_data, df]))

        # concat with master_db
        self.master_db = pd.concat([self.master_db, concat_data])
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
        data = pd.read_sql('SELECT * FROM [CF - regionalized - MarEutro - aggregated]', self.conn)

        concat_data = pd.DataFrame()

        for ix, row in data.iterrows():
            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Marine eutrophication', 'kg N N-lim eq', 'Air', '(unspecified)',
                                                   'Ammonia, ' + row['Short name ecoinvent'], '7664-41-7',
                                                   row['CF NH3 (kg N N-lim eq)'], 'kg', 'Midpoint', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Marine eutrophication', 'PDF.m2.yr', 'Air', '(unspecified)',
                                                   'Ammonia, ' + row['Short name ecoinvent'], '7664-41-7',
                                                   row['CF NH3 (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Marine eutrophication', 'kg N N-lim eq', 'Air', '(unspecified)',
                                                   'Nitrogen oxides, ' + row['Short name ecoinvent'], '11104-93-1',
                                                   row['CF NOx (kg N N-lim eq)'], 'kg', 'Midpoint', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Marine eutrophication', 'PDF.m2.yr', 'Air', '(unspecified)',
                                                   'Nitrogen oxides, ' + row['Short name ecoinvent'], '11104-93-1',
                                                   row['CF NOx (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Marine eutrophication', 'kg N N-lim eq', 'Air', '(unspecified)',
                                                   'Nitric acid, ' + row['Short name ecoinvent'], '7697-37-2',
                                                   row['CF HNO3 (kg N N-lim eq)'], 'kg', 'Midpoint', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Marine eutrophication', 'PDF.m2.yr', 'Air', '(unspecified)',
                                                   'Nitric acid, ' + row['Short name ecoinvent'], '7697-37-2',
                                                   row['CF HNO3 (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

        concat_data = concat_data.reset_index().drop('index', axis=1)

        concat_data.loc[[i for i in concat_data.index if
                         concat_data.loc[i, 'Elem flow name'].split(', ')[-1] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF',
                                                                                  'RME', 'UN-OCEANIA']],
                        'Native geographical resolution scale'] = 'Continent'
        concat_data.loc[[i for i in concat_data.index if concat_data.loc[i, 'Elem flow name'].split(', ')[
            -1] == 'GLO'], 'Native geographical resolution scale'] = 'Global'
        concat_data.loc[[i for i in concat_data.index if concat_data.loc[i, 'Elem flow name'].split(', ')[
            -1] == 'RoW'], 'Native geographical resolution scale'] = 'Other region'

        # add non-regionalized flows (water emissions)
        concat_data = clean_up_dataframe(pd.concat(
            [concat_data, pd.read_sql('SELECT * FROM [CF - not regionalized - MarEutro]', self.conn)]))
        concat_data.loc[concat_data.Compartment == 'Water', 'Native geographical resolution scale'] = 'Not regionalized'

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - Stoechiometry]', self.conn)
        stoc = stoc.loc[stoc.loc[:, 'Impact category'] == 'Marine eutrophication']

        for ix in stoc.index:
            proxy = stoc.loc[ix, 'Proxy molecule']
            comp = stoc.loc[ix, 'Compartment']
            df = pd.DataFrame()
            if proxy == 'NH3':
                # don't create a duplicate
                if stoc.loc[ix, 'Elem flow name'] != 'Ammonia':
                    df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Ammonia') &
                                     ~(concat_data.loc[:, 'Elem flow name'].str.contains('Ammonia, as N'))].loc[
                        concat_data.loc[:, 'Compartment'] == comp].copy('deep')
                    df.loc[:, 'Elem flow name'] = [i.replace('Ammonia', stoc.loc[ix, 'Elem flow name']) for i in
                                                   df.loc[:, 'Elem flow name']]
            if proxy == 'HNO3':
                # don't create a duplicate
                if stoc.loc[ix, 'Elem flow name'] != 'Nitric acid':
                    df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Nitric acid')].loc[
                        concat_data.loc[:, 'Compartment'] == comp].copy('deep')
                    df.loc[:, 'Elem flow name'] = [i.replace('Nitric acid', stoc.loc[ix, 'Elem flow name']) for i in
                                                   df.loc[:, 'Elem flow name']]
            if proxy == 'NOx':
                # don't create a duplicate
                if stoc.loc[ix, 'Elem flow name'] != 'Nitrogen oxides':
                    df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Nitrogen oxides')].loc[
                        concat_data.loc[:, 'Compartment'] == comp].copy('deep')
                    df.loc[:, 'Elem flow name'] = [i.replace('Nitrogen oxides', stoc.loc[ix, 'Elem flow name']) for i in
                                                   df.loc[:, 'Elem flow name']]
            if not df.empty:
                df.loc[:, 'CAS number'] = stoc.loc[ix, 'CAS number']
                df.loc[:, 'CF value'] *= stoc.loc[ix, 'Proxy ratio']

                concat_data = clean_up_dataframe(pd.concat([concat_data, df]))

        # concat with master_db
        self.master_db = pd.concat([self.master_db, concat_data])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_resources_services_loss_cfs(self):
        """
        Loading the CFs for the resources services loss adaptation impact category.

        Concerned impact categories:
            - Resources services loss
            - Resources services loss adaptation

        :return: updated master_db
        """

        # ----------------------------- RESEDA ----------------------------------------

        data = pd.read_sql('SELECT * FROM [CF - not regionalized - ResourcesServicesLoss]', self.conn)
        data.columns = ['Elem flow name', 'CAS number', 'CF value', 'Status', 'Elem flow unit']

        # uncomment if you want to remove isotopes
        # data = data[data.Status != 'isotope']
        data = data.drop('Status', axis=1)

        # add metadata
        data.loc[:, 'Impact category'] = 'Resources services loss'
        data.loc[:, 'CF unit'] = 'kg deficit'
        data.loc[:, 'MP or Damage'] = 'Midpoint'
        data.loc[:, 'Native geographical resolution scale'] = 'Not regionalized'

        # the CFs must be created for all compartments and sub-compartments possible
        data.loc[:, 'Compartment'] = 'Air'
        data.loc[:, 'Sub-compartment'] = '(unspecified)'
        water_comp_data = data.copy()
        water_comp_data.loc[:, 'Compartment'] = 'Water'
        soil_comp_data = data.copy()
        soil_comp_data.loc[:, 'Compartment'] = 'Soil'

        data = pd.concat([data, water_comp_data, soil_comp_data]).reset_index().drop('index', axis=1)

        subcomps_air = ['high. pop.', 'indoor', 'low. pop.', 'low. pop., long-term', 'stratosphere + troposphere']
        subcomps_water = ['lake', 'river', 'ocean', 'groundwater', 'groundwater, long-term']
        subcomps_soil = ['agricultural', 'forestry', 'industrial']

        for subcomp in subcomps_air:
            df = data[data.Compartment == 'Air'].copy()
            df.loc[:, 'Sub-compartment'] = subcomp
            data = pd.concat([data, df])
        for subcomp in subcomps_water:
            df = data[data.Compartment == 'Water'].copy()
            df.loc[:, 'Sub-compartment'] = subcomp
            data = pd.concat([data, df])
        for subcomp in subcomps_soil:
            df = data[data.Compartment == 'Soil'].copy()
            df.loc[:, 'Sub-compartment'] = subcomp
            data = pd.concat([data, df])

        data = data.reset_index().drop('index', axis=1)

        # concat with master_db
        self.master_db = pd.concat([self.master_db, data])
        self.master_db = clean_up_dataframe(self.master_db)

        # ------------------------------ ACP ----------------------------------------

        data = pd.read_sql('SELECT * FROM [CF - not regionalized - ResourcesServicesLossAdaptation]', self.conn)
        data.columns = ['Elem flow name', 'CAS number', 'CF value', 'Status', 'Elem flow unit']

        # uncomment if you want to remove isotopes
        # data = data[data.Status != 'isotope']
        data = data.drop('Status', axis=1)

        # add metadata
        data.loc[:, 'Impact category'] = 'Resources services loss (adaptation)'
        data.loc[:, 'CF unit'] = 'MJ'
        data.loc[:, 'MP or Damage'] = 'Midpoint'
        data.loc[:, 'Native geographical resolution scale'] = 'Not regionalized'

        # the CFs must be created for all compartments and sub-compartments possible
        data.loc[:, 'Compartment'] = 'Air'
        data.loc[:, 'Sub-compartment'] = '(unspecified)'
        water_comp_data = data.copy()
        water_comp_data.loc[:, 'Compartment'] = 'Water'
        soil_comp_data = data.copy()
        soil_comp_data.loc[:, 'Compartment'] = 'Soil'

        data = pd.concat([data, water_comp_data, soil_comp_data]).reset_index().drop('index', axis=1)

        subcomps_air = ['high. pop.', 'indoor', 'low. pop.', 'low. pop., long-term', 'stratosphere + troposphere']
        subcomps_water = ['lake', 'river', 'ocean', 'groundwater', 'groundwater, long-term']
        subcomps_soil = ['agricultural', 'forestry', 'industrial']

        for subcomp in subcomps_air:
            df = data[data.Compartment == 'Air'].copy()
            df.loc[:, 'Sub-compartment'] = subcomp
            data = pd.concat([data, df])
        for subcomp in subcomps_water:
            df = data[data.Compartment == 'Water'].copy()
            df.loc[:, 'Sub-compartment'] = subcomp
            data = pd.concat([data, df])
        for subcomp in subcomps_soil:
            df = data[data.Compartment == 'Soil'].copy()
            df.loc[:, 'Sub-compartment'] = subcomp
            data = pd.concat([data, df])

        data = data.reset_index().drop('index', axis=1)

        # concat with master_db
        self.master_db = pd.concat([self.master_db, data])
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
        data = pd.read_sql('SELECT * FROM [CF - regionalized - EutroFW - aggregated]', self.conn)

        # ------------------------------ FORMAT DATA ------------------------------------
        concat_data = pd.DataFrame()

        for ix, row in data.iterrows():
            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(
                                         ['Freshwater eutrophication', 'kg PO4 P-lim eq', 'Water', '(unspecified)',
                                          'Phosphate, ' + row['Short name ecoinvent'], '14265-44-2',
                                          row['CF PO4 (kg PO4 P-lim eq)'], 'kg', 'Midpoint', 'Country'],
                                         index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                'MP or Damage', 'Native geographical resolution scale']).T])

            concat_data = pd.concat([concat_data,
                                     pd.DataFrame(['Freshwater eutrophication', 'PDF.m2.yr', 'Water', '(unspecified)',
                                                   'Phosphate, ' + row['Short name ecoinvent'], '14265-44-2',
                                                   row['CF PO4 (PDF.m2.yr)'], 'kg', 'Damage', 'Country'],
                                                  index=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment',
                                                         'Elem flow name', 'CAS number', 'CF value', 'Elem flow unit',
                                                         'MP or Damage', 'Native geographical resolution scale']).T])

        concat_data = concat_data.reset_index().drop('index', axis=1)

        concat_data.loc[[i for i in concat_data.index if
                         concat_data.loc[i, 'Elem flow name'].split(', ')[-1] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF',
                                                                                  'RME', 'UN-OCEANIA']],
                        'Native geographical resolution scale'] = 'Continent'
        concat_data.loc[[i for i in concat_data.index if concat_data.loc[i, 'Elem flow name'].split(', ')[
            -1] == 'GLO'], 'Native geographical resolution scale'] = 'Global'

        # ------------------------------ APPLYING STOECHIOMETRIC RATIOS --------------------------
        stoc = pd.read_sql('SELECT * FROM [SI - Stoechiometry]', self.conn)
        stoc = stoc.loc[stoc.loc[:, 'Impact category'] == 'Freshwater eutrophication']

        for ix in stoc.index:
            df = concat_data[concat_data.loc[:, 'Elem flow name'].str.contains('Phosphate')].loc[
                concat_data.loc[:, 'Compartment'] == stoc.loc[ix, 'Compartment']].copy('deep')

            df.loc[:, 'Elem flow name'] = [i.replace('Phosphate', stoc.loc[ix, 'Elem flow name']) for i in
                                           df.loc[:, 'Elem flow name']]
            df.loc[:, 'CAS number'] = stoc.loc[ix, 'CAS number']
            df.loc[:, 'CF value'] *= stoc.loc[ix, 'Proxy ratio']

            concat_data = clean_up_dataframe(pd.concat([concat_data, df]))

        # concat with master_db
        self.master_db = pd.concat([self.master_db, concat_data])
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
        data = pd.read_sql('SELECT * FROM [CF - regionalized - Land use - aggregated]', self.conn)

        data = data.loc[:, ['Elem flow name', 'CFs (PDF.m2.yr)']]
        data = data.rename(columns={'CFs (PDF.m2.yr)': 'CF value'})
        data.loc[data.loc[:, 'Elem flow name'].str.contains(
            'Occupation'), 'Impact category'] = 'Land occupation, biodiversity'
        data.loc[data.loc[:, 'Elem flow name'].str.contains(
            'Transformation'), 'Impact category'] = 'Land transformation, biodiversity'
        data.loc[data.loc[:, 'Elem flow name'].str.contains('Occupation'), 'Elem flow unit'] = 'm2.yr'
        data.loc[data.loc[:, 'Elem flow name'].str.contains('Transformation'), 'Elem flow unit'] = 'm2'
        data.loc[:, 'CF unit'] = 'PDF.m2.yr'
        data.loc[:, 'Compartment'] = 'Raw'
        data.loc[:, 'Sub-compartment'] = 'land'
        data.loc[:, 'MP or Damage'] = 'Damage'
        data.loc[:, 'Native geographical resolution scale'] = 'Country'

        # calculate midpoint values
        midpoint_data = data.copy()
        midpoint_data.loc[[i for i in midpoint_data.index if 'Occupation' in midpoint_data.loc[i, 'Elem flow name']],
                          'CF value'] /= (
            midpoint_data.loc[midpoint_data.loc[:, 'Elem flow name'] ==
                              'Occupation, annual crops, GLO', 'CF value'].iloc[0]
        )
        midpoint_data.loc[
            [i for i in midpoint_data.index if 'Transformation' in midpoint_data.loc[i, 'Elem flow name']],
            'CF value'] /= abs(
            midpoint_data.loc[midpoint_data.loc[:, 'Elem flow name'] ==
                              'Transformation, from annual crops, GLO', 'CF value'].iloc[0]
        )
        midpoint_data.loc[:, 'MP or Damage'] = 'Midpoint'
        midpoint_data.loc[
            midpoint_data.loc[:, 'Elem flow name'].str.contains('Occupation'), 'CF unit'] = 'm2 arable land eq .yr'
        midpoint_data.loc[
            midpoint_data.loc[:, 'Elem flow name'].str.contains('Transformation'), 'CF unit'] = 'm2 arable land eq'

        data = pd.concat([data, midpoint_data]).reset_index().drop('index', axis=1)

        data.loc[[i for i in data.index if
                  data.loc[i, 'Elem flow name'].split(', ')[-1] in ['RNA', 'RLA', 'RER', 'RAS', 'RAF',
                                                                    'RME', 'UN-OCEANIA']],
                 'Native geographical resolution scale'] = 'Continent'
        data.loc[[i for i in data.index if
                  data.loc[i, 'Elem flow name'].split(', ')[-1] == 'GLO'],
                 'Native geographical resolution scale'] = 'Global'

        data.loc[[i for i in data.index if
                  data.loc[i, 'Elem flow name'].split(', ')[-1] == 'RoW'],
                 'Native geographical resolution scale'] = 'Other region'

        forest_not_used = data.loc[data.loc[:, 'Elem flow name'].str.contains('artificial areas')].copy()
        forest_not_used.loc[:, 'Elem flow name'] = [i.replace('artificial areas', 'forest/grassland, not used') for i in
                                                    forest_not_used.loc[:, 'Elem flow name']]
        forest_not_used.loc[:, 'CF value'] = 0

        data = pd.concat([data, forest_not_used]).reset_index().drop('index', axis=1)

        # concat with master_db
        self.master_db = pd.concat([self.master_db, data])
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
                                                'United States', 'Australia', 'Brazil', 'China', 'India', 'Indonesia']

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
                            water_data.loc[str(i + j), 'Elem flow name'] = 'Water, ' + \
                                                                           conc.loc[data.ecoinvent_shortname[i]].iloc[j]
                        elif data.loc[i, 'Water type'] == 'agri':
                            water_data.loc[str(i + j), 'Elem flow name'] = 'Water, agri, ' + \
                                                                           conc.loc[data.ecoinvent_shortname[i]].iloc[j]
                        elif data.loc[i, 'Water type'] == 'non-agri':
                            water_data.loc[str(i + j), 'Elem flow name'] = 'Water, non-agri, ' + \
                                                                           conc.loc[data.ecoinvent_shortname[i]].iloc[j]
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
                                                    'Australia, including overseas territories' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'BALTSO' in all_data.loc[i, 'Elem flow name'] or
                                                    'CENTREL' in all_data.loc[i, 'Elem flow name'] or
                                                    'CUSMA/T-MEC/USMCA' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta and Quebec' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Canada without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canary Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'Central Asia' in all_data.loc[i, 'Elem flow name'] or
                                                    'China w/o Inner Mongol' in all_data.loc[i, 'Elem flow name'] or
                                                    'Crimea' in all_data.loc[i, 'Elem flow name'] or
                                                    'Cyprus No Mans Area' in all_data.loc[i, 'Elem flow name'] or
                                                    'Dhekelia Base' in all_data.loc[i, 'Elem flow name'] or
                                                    'ENTSO-E' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Austria' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without NORDEL (NCPA)' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Europe without Switzerland' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland and Austria' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Europe without Switzerland and France' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Europe, without Russia and Türkiye' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'FSU' in all_data.loc[i, 'Elem flow name'] or
                                                    'France, including overseas territories' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Guantanamo Bay' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Africa' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Asia, without China and GCC' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, EU27 & EFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Gulf Cooperation Council' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, North America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, North America, without Quebec' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, Russia & RER w/o EU27 & EFTA' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, South America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IN-Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'MRO' in all_data.loc[i, 'Elem flow name'] or
                                                    'NAFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'NORDEL' in all_data.loc[i, 'Elem flow name'] or
                                                    'NPCC' in all_data.loc[i, 'Elem flow name'] or
                                                    'North America without Quebec' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Northern Cyprus' in all_data.loc[i, 'Elem flow name'] or
                                                    'Québec, HQ distribution network' in all_data.loc[
                                                        i, 'Elem flow name'] or
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
                                                    'UCTE without Germany and France' in all_data.loc[
                                                        i, 'Elem flow name'] or
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
                                                    'United States of America, including overseas territories' in
                                                    all_data.loc[i, 'Elem flow name'] or
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
                                  'agri' not in all_data.loc[i, 'Elem flow name'] and all_data.loc[
                                      i, 'Compartment'] == 'Raw'])

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
                        water_data.loc[str(i + j), 'Elem flow name'] = 'Water, ' + \
                                                                       conc.loc[data.ecoinvent_shortname[i]].iloc[j]
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
                                                    'Australia, including overseas territories' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'BALTSO' in all_data.loc[i, 'Elem flow name'] or
                                                    'CENTREL' in all_data.loc[i, 'Elem flow name'] or
                                                    'CUSMA/T-MEC/USMCA' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta and Quebec' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Canada without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canary Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'Central Asia' in all_data.loc[i, 'Elem flow name'] or
                                                    'China w/o Inner Mongol' in all_data.loc[i, 'Elem flow name'] or
                                                    'Crimea' in all_data.loc[i, 'Elem flow name'] or
                                                    'Cyprus No Mans Area' in all_data.loc[i, 'Elem flow name'] or
                                                    'Dhekelia Base' in all_data.loc[i, 'Elem flow name'] or
                                                    'ENTSO-E' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Austria' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without NORDEL (NCPA)' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Europe without Switzerland' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland and Austria' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Europe without Switzerland and France' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Europe, without Russia and Türkiye' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'FSU' in all_data.loc[i, 'Elem flow name'] or
                                                    'France, including overseas territories' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Guantanamo Bay' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Africa' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Asia, without China and GCC' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, EU27 & EFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Gulf Cooperation Council' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, North America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, North America, without Quebec' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, Russia & RER w/o EU27 & EFTA' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, South America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IN-Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'MRO' in all_data.loc[i, 'Elem flow name'] or
                                                    'NAFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'NORDEL' in all_data.loc[i, 'Elem flow name'] or
                                                    'NPCC' in all_data.loc[i, 'Elem flow name'] or
                                                    'North America without Quebec' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Northern Cyprus' in all_data.loc[i, 'Elem flow name'] or
                                                    'Québec, HQ distribution network' in all_data.loc[
                                                        i, 'Elem flow name'] or
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
                                                    'UCTE without Germany and France' in all_data.loc[
                                                        i, 'Elem flow name'] or
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
                                                    'United States of America, including overseas territories' in
                                                    all_data.loc[i, 'Elem flow name'] or
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
                                                    'Australia, including overseas territories' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'BALTSO' in all_data.loc[i, 'Elem flow name'] or
                                                    'CENTREL' in all_data.loc[i, 'Elem flow name'] or
                                                    'CUSMA/T-MEC/USMCA' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canada without Alberta and Quebec' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Canada without Quebec' in all_data.loc[i, 'Elem flow name'] or
                                                    'Canary Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'Central Asia' in all_data.loc[i, 'Elem flow name'] or
                                                    'China w/o Inner Mongol' in all_data.loc[i, 'Elem flow name'] or
                                                    'Crimea' in all_data.loc[i, 'Elem flow name'] or
                                                    'Cyprus No Mans Area' in all_data.loc[i, 'Elem flow name'] or
                                                    'Dhekelia Base' in all_data.loc[i, 'Elem flow name'] or
                                                    'ENTSO-E' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Austria' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without NORDEL (NCPA)' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Europe without Switzerland' in all_data.loc[i, 'Elem flow name'] or
                                                    'Europe without Switzerland and Austria' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Europe without Switzerland and France' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Europe, without Russia and Türkiye' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'FSU' in all_data.loc[i, 'Elem flow name'] or
                                                    'France, including overseas territories' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Guantanamo Bay' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Africa' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Asia, without China and GCC' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, EU27 & EFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, Gulf Cooperation Council' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, North America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IAI Area, North America, without Quebec' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, Russia & RER w/o EU27 & EFTA' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'IAI Area, South America' in all_data.loc[i, 'Elem flow name'] or
                                                    'IN-Islands' in all_data.loc[i, 'Elem flow name'] or
                                                    'MRO' in all_data.loc[i, 'Elem flow name'] or
                                                    'NAFTA' in all_data.loc[i, 'Elem flow name'] or
                                                    'NORDEL' in all_data.loc[i, 'Elem flow name'] or
                                                    'NPCC' in all_data.loc[i, 'Elem flow name'] or
                                                    'North America without Quebec' in all_data.loc[
                                                        i, 'Elem flow name'] or
                                                    'Northern Cyprus' in all_data.loc[i, 'Elem flow name'] or
                                                    'Québec, HQ distribution network' in all_data.loc[
                                                        i, 'Elem flow name'] or
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
                                                    'UCTE without Germany and France' in all_data.loc[
                                                        i, 'Elem flow name'] or
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
                                                    'United States of America, including overseas territories' in
                                                    all_data.loc[i, 'Elem flow name'] or
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
                    water_data.loc[i, 'Elem flow name'] = 'Water, well, in ground, ' + conc.loc[
                        data.ecoinvent_shortname[i]]
                    water_data.loc[i, 'CF value'] = data.loc[i, 'CF (PDF.m2.yr/m3)']
                else:
                    for j in range(0, len(conc.loc[data.ecoinvent_shortname[i]])):
                        water_data.loc[str(i + j), 'Elem flow name'] = 'Water, well, in ground, ' + \
                                                                       conc.loc[data.ecoinvent_shortname[i]].iloc[j]
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
                                                        'Australia, including overseas territories' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'BALTSO' in water_data.loc[i, 'Elem flow name'] or
                                                        'CENTREL' in water_data.loc[i, 'Elem flow name'] or
                                                        'CUSMA/T-MEC/USMCA' in water_data.loc[i, 'Elem flow name'] or
                                                        'Canada without Alberta' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Canada without Alberta and Quebec' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Canada without Quebec' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Canary Islands' in water_data.loc[i, 'Elem flow name'] or
                                                        'Central Asia' in water_data.loc[i, 'Elem flow name'] or
                                                        'China w/o Inner Mongol' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Crimea' in water_data.loc[i, 'Elem flow name'] or
                                                        'Cyprus No Mans Area' in water_data.loc[i, 'Elem flow name'] or
                                                        'Dhekelia Base' in water_data.loc[i, 'Elem flow name'] or
                                                        'ENTSO-E' in water_data.loc[i, 'Elem flow name'] or
                                                        'Europe without Austria' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Europe without NORDEL (NCPA)' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Europe without Switzerland' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Europe without Switzerland and Austria' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Europe without Switzerland and France' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Europe, without Russia and Türkiye' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'FSU' in water_data.loc[i, 'Elem flow name'] or
                                                        'France, including overseas territories' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Guantanamo Bay' in water_data.loc[i, 'Elem flow name'] or
                                                        'IAI Area, Africa' in water_data.loc[i, 'Elem flow name'] or
                                                        'IAI Area, Asia, without China and GCC' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, EU27 & EFTA' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, Gulf Cooperation Council' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, North America' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, North America, without Quebec' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, Russia & RER w/o EU27 & EFTA' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, South America' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IN-Islands' in water_data.loc[i, 'Elem flow name'] or
                                                        'MRO' in water_data.loc[i, 'Elem flow name'] or
                                                        'NAFTA' in water_data.loc[i, 'Elem flow name'] or
                                                        'NORDEL' in water_data.loc[i, 'Elem flow name'] or
                                                        'NPCC' in water_data.loc[i, 'Elem flow name'] or
                                                        'North America without Quebec' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Northern Cyprus' in water_data.loc[i, 'Elem flow name'] or
                                                        'Québec, HQ distribution network' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'RER w/o AT+BE+CH+DE+FR+IT' in water_data.loc[
                                                            i, 'Elem flow name'] or
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
                                                        'UCTE without Germany and France' in water_data.loc[
                                                            i, 'Elem flow name'] or
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
                                                        'United States of America, including overseas territories' in
                                                        water_data.loc[i, 'Elem flow name'] or
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
                        water_data.loc[str(i + j), 'Elem flow name'] = 'Water, cooling, unspecified natural origin, ' + \
                                                                       conc.loc[geo].iloc[j]
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
                                                        'Australia, including overseas territories' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'BALTSO' in water_data.loc[i, 'Elem flow name'] or
                                                        'CENTREL' in water_data.loc[i, 'Elem flow name'] or
                                                        'CUSMA/T-MEC/USMCA' in water_data.loc[i, 'Elem flow name'] or
                                                        'Canada without Alberta' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Canada without Alberta and Quebec' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Canada without Quebec' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Canary Islands' in water_data.loc[i, 'Elem flow name'] or
                                                        'Central Asia' in water_data.loc[i, 'Elem flow name'] or
                                                        'China w/o Inner Mongol' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Crimea' in water_data.loc[i, 'Elem flow name'] or
                                                        'Cyprus No Mans Area' in water_data.loc[i, 'Elem flow name'] or
                                                        'Dhekelia Base' in water_data.loc[i, 'Elem flow name'] or
                                                        'ENTSO-E' in water_data.loc[i, 'Elem flow name'] or
                                                        'Europe without Austria' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Europe without NORDEL (NCPA)' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Europe without Switzerland' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Europe without Switzerland and Austria' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Europe without Switzerland and France' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Europe, without Russia and Türkiye' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'FSU' in water_data.loc[i, 'Elem flow name'] or
                                                        'France, including overseas territories' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Guantanamo Bay' in water_data.loc[i, 'Elem flow name'] or
                                                        'IAI Area, Africa' in water_data.loc[i, 'Elem flow name'] or
                                                        'IAI Area, Asia, without China and GCC' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, EU27 & EFTA' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, Gulf Cooperation Council' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, North America' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, North America, without Quebec' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, Russia & RER w/o EU27 & EFTA' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IAI Area, South America' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'IN-Islands' in water_data.loc[i, 'Elem flow name'] or
                                                        'MRO' in water_data.loc[i, 'Elem flow name'] or
                                                        'NAFTA' in water_data.loc[i, 'Elem flow name'] or
                                                        'NORDEL' in water_data.loc[i, 'Elem flow name'] or
                                                        'NPCC' in water_data.loc[i, 'Elem flow name'] or
                                                        'North America without Quebec' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'Northern Cyprus' in water_data.loc[i, 'Elem flow name'] or
                                                        'Québec, HQ distribution network' in water_data.loc[
                                                            i, 'Elem flow name'] or
                                                        'RER w/o AT+BE+CH+DE+FR+IT' in water_data.loc[
                                                            i, 'Elem flow name'] or
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
                                                        'UCTE without Germany and France' in water_data.loc[
                                                            i, 'Elem flow name'] or
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
                                                        'United States of America, including overseas territories' in
                                                        water_data.loc[i, 'Elem flow name'] or
                                                        'WECC' in water_data.loc[i, 'Elem flow name'] or
                                                        'WEU' in water_data.loc[i, 'Elem flow name'])],
                       'Native geographical resolution scale'] = 'Other region'

        # concat with master_db
        self.master_db = pd.concat([self.master_db, water_data])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_physical_effects_cfs(self):
        """
        Load CFs for plastics and naturel fibers physical effects on biota impact.

        Concerned impact categories:
            - Physical effect on biota

        :return: update master_db
        """

        original_cfs = pd.read_sql('SELECT * from "CF - not regionalized - PhysicalImpactonBiota"',
                                   con=self.conn)

        original_cfs.drop(['Geometric st.dev.', 'Lower limit 95% CI', 'Upper limit 95% CI'], axis=1, inplace=True)
        original_cfs.loc[:, 'Impact category'] = 'Physical effects on biota'
        original_cfs.loc[:, 'Native geographical resolution scale'] = 'Global'
        original_cfs.loc[:, 'MP or Damage'] = ['Midpoint' if original_cfs.loc[i, 'CF unit'] == 'CTUe' else 'Damage' for
                                               i in original_cfs.index]

        CAS = {'EPS': '9003-53-6',
               'HDPE': '25087-34-7',
               'LDPE': '9002-88-4',
               'PA/Nylon': '25038-54-4',
               'PET': '25038-59-9',
               'PHA': '117068-64-1',
               'PLA': '26100-51-6',
               'PP': '9003-07-0',
               'PS': '9003-53-6',
               'PVC': '9002-86-2',
               'TRWP': None,
               'Cotton': None,
               'Linen': None,
               'Lyocell': None,
               'Modal': None,
               'PAN/Acrylic': '9065-11-6',
               'PU/Spandex': '9009-54-5',
               'Rayon': None,
               'Viscose': None}

        original_cfs.loc[:, 'CAS number'] = [CAS[i] for i in original_cfs.loc[:, 'Polymer type']]
        original_cfs = original_cfs.rename(columns={'Recommended CF (geometric mean)': 'CF value'})
        original_cfs.loc[:, 'Elem flow name'] = [
            original_cfs.loc[i, 'Shape'] + ' - ' + original_cfs.loc[i, 'Polymer type'] + ' (' +
            str(original_cfs.loc[i, 'Size']) + ' µm diameter)' for i in original_cfs.index]

        # we create CFs for default sizes which depend on the shape of the microplastics
        default_sizes = {'Beads/spheres': 1000,
                         'Film fragments': 100,
                         'Microfibers/cylinders': 10}
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
        df = self.master_db.loc[self.master_db.loc[:, 'Impact category'] == 'Physical effects on biota'].loc[
            self.master_db.loc[:, 'Sub-compartment'] == 'lake'].copy()
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
            [data,
             pd.concat([pd.DataFrame([i[1] + ' fish, discarded, GLO' for i in glo.index], columns=['Elem flow name']),
                        pd.DataFrame(glo.loc[:, 'Discard CF (PDF.m2.yr/t discarded fish)'].values,
                                     columns=['CF value'])], axis=1)])
        data = clean_up_dataframe(data)
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

        for substance in ['Ammonia', 'Nitrogen oxides', 'Sulfur dioxide']:
            all_existing_geos = set(self.master_db.loc[(self.master_db.loc[:, 'Elem flow name'].str.contains(
                substance + ', ') & ~(self.master_db.loc[:, 'Elem flow name'].str.contains(', as N, '))),
                                                       'Elem flow name'])
            all_existing_geos = [i.split(substance + ', ')[1] for i in all_existing_geos]

            for indicator in ['Particulate matter formation', 'Freshwater acidification', 'Terrestrial acidification',
                              'Marine eutrophication']:
                existing_geos = set(self.master_db.loc[self.master_db.loc[:, 'Impact category'] == indicator].loc[
                                        self.master_db.loc[:, 'Elem flow name'].str.contains(
                                            substance + ', '), 'Elem flow name'])
                existing_geos = [i.split(substance + ', ')[1] for i in existing_geos]
                if existing_geos:
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

        # also ensure consistency for substance made from stochiometric ratios
        harmonization_stoechiometry = {'Ammonia, as N': 'Ammonia',
                                       'Ammonium carbonate': 'Ammonia',
                                       'Ammonium nitrate': 'Nitrogen oxides',
                                       'Ammonium, ion': 'Ammonia',
                                       'Nitrate': 'Nitrogen oxides',
                                       'Nitric acid': 'Nitrogen oxides',
                                       'Nitrite': 'Nitrogen oxides',
                                       'Nitrogen dioxide': 'Nitrogen oxides',
                                       'Nitric oxide': 'Nitrogen oxides',
                                       'Sulfate': 'Sulfur dioxide',
                                       'Sulfur trioxide': 'Sulfur dioxide',
                                       'Sulfuric acid': 'Sulfur dioxide'
                                       }

        for substance in harmonization_stoechiometry:
            for indicator in ['Particulate matter formation', 'Freshwater acidification',
                              'Terrestrial acidification', 'Marine eutrophication']:
                existing_geos = set(self.master_db.loc[self.master_db.loc[:, 'Impact category'] == indicator].loc[
                                        self.master_db.loc[:, 'Elem flow name'].str.contains(
                                            substance + ', '), 'Elem flow name'])
                existing_geos = [i.split(substance + ', ')[1] for i in existing_geos]
                if existing_geos:
                    existing_geos_ref = set(
                        self.master_db.loc[self.master_db.loc[:, 'Impact category'] == indicator].loc[
                            (self.master_db.loc[:, 'Elem flow name'].str.contains(
                                harmonization_stoechiometry[substance] + ', ') & ~(
                                self.master_db.loc[:, 'Elem flow name'].str.contains(
                                    ', as N, '))), 'Elem flow name'])
                    existing_geos_ref = [i.split(harmonization_stoechiometry[substance] + ', ')[1] for i in
                                         existing_geos_ref]

                    missing_geos = set(existing_geos_ref) - set(existing_geos)
                    for missing_geo in missing_geos:
                        df = self.master_db.loc[self.master_db.loc[:, 'Impact category'] == indicator].loc[
                            self.master_db.loc[:, 'Elem flow name'] == substance + ', ' + map.loc[
                                missing_geo, 'continent']].copy('deep')
                        df.loc[:, 'Elem flow name'] = substance + ', ' + missing_geo
                        if ('-' in missing_geo and missing_geo not in ['ENTSO-E', 'UN-SEASIA']) or len(
                                missing_geo) == 2:
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
        proxy.set_index(['Impact category', 'CF unit', 'Compartment', 'Sub-compartment', 'Elem flow name'],
                        inplace=True)
        proxy.update(self.master_db.set_index(['Impact category', 'CF unit', 'Compartment',
                                               'Sub-compartment', 'Elem flow name']))
        proxy = proxy.reset_index()
        self.master_db = pd.concat([self.master_db, proxy]).drop_duplicates()

        self.master_db = clean_up_dataframe(self.master_db)

        # ------------ groundwater and ocean subcomps ------------

        # for the groundwater and ocean subcomps, only in some impact categories are the values equal to unspecified
        water_comp = self.master_db.loc[
            [i for i in self.master_db.index if self.master_db.loc[i, 'Compartment'] == 'Water']]
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
                             'Ionizing radiations, human health',
                             'Ionizing radiations',
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

        long_term_cats = ['Climate change, ecosystem quality, terrestrial ecosystem', 'Climate change, human health',
                          'Climate change, ecosystem quality, marine ecosystem', 'Freshwater ecotoxicity',
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
        Some regionalized emissions (e.g., Ammonia) also impact non-regionalized impact categories. Hence, this method
        creates non-regionalized CFs for these emissions.
        :return:
        """

        regio_flows_per_ic = {}

        for ic in set(self.master_db.loc[:, 'Impact category']):
            df = self.master_db[self.master_db['Impact category'] == ic].copy()
            dff = df.loc[
                df.loc[:, 'Elem flow name'].str.contains(', FR'), ['Elem flow name', 'Compartment']].drop_duplicates()
            regio_flows_per_ic[ic] = set(zip(dff['Elem flow name'], dff['Compartment']))
        regio_flows_per_ic = {k: set([(i[0].split(', FR')[0], i[1]) for i in v]) for k, v in regio_flows_per_ic.items()}

        flows_to_create = {}

        for ic in set(self.master_db.loc[:, 'Impact category']):
            df = self.master_db[self.master_db['Impact category'] == ic].copy()
            for flow in set.union(*[set(i) for i in list(regio_flows_per_ic.values())]):
                dff = df.loc[(df.loc[:, 'Elem flow name'] == flow[0]) & (df.loc[:, 'Compartment'] == flow[1])]
                # if the CF is equal to zero we don't care
                if dff.loc[:, 'CF value'].sum() != 0:
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
                self.master_db[self.master_db['Elem flow name'].str.contains(substance[0], na=False)].index.tolist()]
            df = df[df.loc[:, 'Native geographical resolution scale'].isin(['Country', 'Continent', 'Other region'])]
            regions = set(
                [i.split(substance[0] + ', ')[-1] for i in df.loc[:, 'Elem flow name'] if substance[0] + ', ' in i])
            ic_flows_to_create = {k for k, v in flows_to_create.items() if substance in v}
            for ic in ic_flows_to_create:
                dff = self.master_db.loc[(self.master_db.loc[:, 'Elem flow name'] == substance[0]) &
                                         (self.master_db.loc[:, 'Compartment'] == substance[1])].loc[
                    self.master_db.loc[:, 'Impact category'] == ic].copy()
                list_flow_added = ([substance[0] + ', ' + i for i in regions] * len(dff))
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
            df.loc[df.loc[:, 'Impact category'] == 'Climate change, ecosystem quality, terrestrial ecosystem, long term', 'CF value'] = -(
                df.loc[df.loc[:, 'Impact category'] == 'Climate change, ecosystem quality, terrestrial ecosystem, short term', 'CF value']).iloc[0]
            df.loc[df.loc[:,'Impact category'] == 'Climate change, ecosystem quality, marine ecosystem, long term', 'CF value'] = - (
                df.loc[df.loc[:,'Impact category'] == 'Climate change, ecosystem quality, marine ecosystem, short term', 'CF value']).iloc[0]
            df.loc[df.loc[:, 'Impact category'] == 'Climate change, human health, long term', 'CF value'] = -(
                df.loc[df.loc[:, 'Impact category'] == 'Climate change, human health, short term', 'CF value']).iloc[0]
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
            i + ', biogenic' for i in self.master_db.loc[
                self.master_db['Elem flow name'].isin(biogenic) &
                self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category']]

        self.master_db.loc[
            self.master_db['Elem flow name'].isin(land_use) &
            self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category'] = [
            i + ', land transformation' for i in self.master_db.loc[
                self.master_db['Elem flow name'].isin(land_use) &
                self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category']]

        self.master_db.loc[
            self.master_db['Elem flow name'].isin(CO2_uptake) &
            self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category'] = [
            i + ', CO2 uptake' for i in self.master_db.loc[
                self.master_db['Elem flow name'].isin(CO2_uptake) &
                self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category']]

        self.master_db.loc[
            ~self.master_db['Elem flow name'].isin(CO2_uptake + biogenic + land_use) &
            self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category'] = [
            i + ', fossil' for i in self.master_db.loc[
                ~self.master_db['Elem flow name'].isin(CO2_uptake + biogenic + land_use) &
                self.master_db['Impact category'].str.contains('Climate change', na=False), 'Impact category']]

    def separate_regio_cfs(self):
        """
        Method to obtain two different versions of master_db. One with regionalized factors that will be used for
        SimaPro and openLCA versions. One with only non regionalized factors that will be used for brightway and
        ecoinvent versions.
        """

        self.master_db_not_regio = self.master_db.loc[[i for i in self.master_db.index if
                                                       self.master_db.loc[
                                                           i, 'Native geographical resolution scale'] not in [
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

        latest_ei_version = '3.12'

        elem_flow_uuid = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/mappings/ei' + latest_ei_version.replace('.', '') + '/ei_elem_flow_uuids.xlsx'))

        for db_format in ['normal', 'carbon neutrality']:
            if db_format == 'normal':
                ei_iw_db = self.master_db_not_regio.copy()
            elif db_format == 'carbon neutrality':
                ei_iw_db = self.master_db_not_regio_carbon_neutrality.copy()

            # -------------- Mapping substances --------------

            mapping = pd.read_excel(pkg_resources.resource_stream(
                __name__, '/Data/mappings/ei' + latest_ei_version.replace('.', '') + '/ei_iw_mapping.xlsx'))
            ei_mapping = mapping.loc[:, ['ecoinvent name', 'iw name']].dropna()
            not_one_for_one = ei_mapping[ei_mapping.loc[:, 'iw name'].duplicated(False)]
            one_for_one = ei_mapping[~ei_mapping.loc[:, 'iw name'].duplicated(False)]

            # for one_for_one it's easy! We just replace one for one
            ei_iw_db.loc[:, 'Elem flow name'] = ei_iw_db.loc[:, 'Elem flow name'].replace(
                one_for_one.loc[:, 'iw name'].tolist(), one_for_one.loc[:, 'ecoinvent name'].tolist())

            # for not_one_for_one it's harder, e.g., the "Zinc" substance from iw+ must be linked to multiple elementary flows in ecoinvent
            unique_not_one_for_one = set(not_one_for_one.loc[:, 'iw name'])
            for subst in tqdm(unique_not_one_for_one, leave=True):
                ei_df = not_one_for_one.loc[
                    [i for i in not_one_for_one.index if not_one_for_one.loc[i, 'iw name'] == subst]]
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
            # CFs for resources "in ground" should only be for Fossil and nuclear energy use
            minerals = [i for i in ei_iw_db.index if (', in ground' in ei_iw_db.loc[i, 'Elem flow name'] and
                                                      ei_iw_db.loc[i, 'Impact category'] not in [
                                                          'Fossil and nuclear energy use'] and
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

            with open(pkg_resources.resource_filename(
                    __name__, "Data/mappings/ei" + latest_ei_version.replace('.', '') + "/comps.json"), "r") as f:
                comps = json.load(f)
            with open(pkg_resources.resource_filename(
                    __name__, "Data/mappings/ei" + latest_ei_version.replace('.', '') + "/subcomps.json"), "r") as f:
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
                self.ei312_iw = ei_iw_db.copy('deep')

                # add ecoinvent elem flow uuids
                self.ei312_iw = self.ei312_iw.merge(
                    elem_flow_uuid.loc[:, ['Name', 'Compartment', 'Subcompartment', 'ID']],
                    right_on=['Name', 'Compartment', 'Subcompartment'],
                    left_on=['Elem flow name', 'Compartment', 'Sub-compartment'],
                    how='left')

                only_in_312 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == 3.12].dropna(
                    subset=['iw name']).loc[:, 'ecoinvent name'])

                self.ei311_iw = self.ei312_iw.drop(
                    [i for i in self.ei312_iw.index if self.ei312_iw.loc[i, 'Elem flow name'] in
                     only_in_312]).copy('deep')

                only_in_311 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == 3.11].dropna(
                    subset=['iw name']).loc[:, 'ecoinvent name'])

                self.ei310_iw = self.ei311_iw.drop(
                    [i for i in self.ei311_iw.index if self.ei311_iw.loc[i, 'Elem flow name'] in
                     only_in_311]).copy('deep')

                self.ei312_iw = self.ei312_iw.dropna(subset=['ID']).drop_duplicates()
                self.ei311_iw = self.ei311_iw.dropna(subset=['ID']).drop_duplicates()
                self.ei310_iw = self.ei310_iw.dropna(subset=['ID']).drop_duplicates()

            elif db_format == 'carbon neutrality':

                self.ei312_iw_carbon_neutrality = ei_iw_db.copy('deep')
                # add ecoinvent elem flow uuids
                self.ei312_iw_carbon_neutrality = self.ei312_iw_carbon_neutrality.merge(
                    elem_flow_uuid.loc[:, ['Name', 'Compartment', 'Subcompartment', 'ID']],
                    right_on=['Name', 'Compartment', 'Subcompartment'],
                    left_on=['Elem flow name', 'Compartment', 'Sub-compartment'],
                    how='left')

                only_in_312 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == 3.12].dropna(
                    subset=['iw name']).loc[:, 'ecoinvent name'])

                self.ei311_iw_carbon_neutrality = self.ei312_iw_carbon_neutrality.drop(
                    [i for i in self.ei312_iw_carbon_neutrality.index if
                     self.ei312_iw_carbon_neutrality.loc[i, 'Elem flow name'] in
                     only_in_312]).copy('deep')

                only_in_311 = list(mapping[mapping.loc[:, 'introduced in ei v.'] == 3.11].dropna(
                    subset=['iw name']).loc[:, 'ecoinvent name'])

                self.ei310_iw_carbon_neutrality = self.ei311_iw_carbon_neutrality.drop(
                    [i for i in self.ei311_iw_carbon_neutrality.index if
                     self.ei311_iw_carbon_neutrality.loc[i, 'Elem flow name'] in
                     only_in_311]).copy('deep')

                self.ei312_iw_carbon_neutrality = self.ei312_iw_carbon_neutrality.dropna(subset=['ID']).drop_duplicates()
                self.ei311_iw_carbon_neutrality = self.ei311_iw_carbon_neutrality.dropna(subset=['ID']).drop_duplicates()
                self.ei310_iw_carbon_neutrality = self.ei310_iw_carbon_neutrality.dropna(subset=['ID']).drop_duplicates()

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
            sp = sp.loc[:, ['Name', 'Name IW+']].dropna()
            differences = sp.loc[sp.Name != sp.loc[:, 'Name IW+']]
            double_iw_flow = sp.loc[sp.loc[:, 'Name IW+'].duplicated(), 'Name IW+'].tolist()

            m = db.merge(differences[['Name IW+', 'Name']], how='left', left_on='Elem flow name', right_on='Name IW+')

            # rows that have any mapping
            matched = m[m['Name'].notna()].copy()
            matched['Elem flow name'] = matched['Name']
            matched = matched.drop(columns=['Name IW+', 'Name'])

            # concat originals + all mapped copies
            db = pd.concat([db, matched], ignore_index=True)

            # remove duplicates names from IW+ (we just added everything together, i.e., SP names AND IW names)
            db = db.loc[~((db.loc[:, "Elem flow name"].isin(differences.loc[:, 'Name IW+'])) & ~(
                db.loc[:, "Elem flow name"].isin(double_iw_flow)))]
            # remove other duplicated flows, basically all coming from land categories
            db = db.loc[~(db.loc[:, 'Elem flow name'].str.contains('annual crops')) &
                        ~(db.loc[:, 'Elem flow name'].str.contains('permanent crops')) &
                        ~(db.loc[:, 'Elem flow name'].str.contains('forest, used, extensive')) &
                        ~(db.loc[:, 'Elem flow name'].str.contains('forest, used, intensive')) &
                        ~(db.loc[:, 'Elem flow name'].str.contains('pasture/meadow')) &
                        ~(db.loc[:, 'Elem flow name'].str.contains('without Russia and Türkiye'))]

            # ------------------------------- ODDITIES ----------------------------------

            # need to change the unit of land occupation flows to match SP nomenclature
            db.loc[db.loc[:, 'Elem flow unit'] == 'm2.yr', 'Elem flow unit'] = 'm2a'

            # Some water names are reserved names in SimaPro, so we modify it
            db.loc[db['Elem flow name'] == 'Water, agri', 'Elem flow name'] = 'Water/m3, agri'
            db.loc[db['Elem flow name'] == 'Water, non-agri', 'Elem flow name'] = 'Water/m3, non-agri'

            # fix comp/subcomp of pesky biogenic carbon elementary flows
            db = db.drop(
                db.loc[db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                         'Carbon dioxide, in air'])].loc[
                    db.loc[:, 'Sub-compartment'] != '(unspecified)'].index)
            db = db.drop(
                db.loc[db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                         'Carbon dioxide, in air'])].loc[
                    (db.loc[:, 'Compartment'] != 'Air')].index)
            db.loc[db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                     'Carbon dioxide, in air']), 'Compartment'] = 'Raw'
            db.loc[db.loc[:, 'Elem flow name'].isin(['Carbon dioxide, non-fossil, resource correction',
                                                     'Carbon dioxide, in air']), 'Sub-compartment'] = 'in air'
            db = db.drop(db.loc[db.loc[:, 'Elem flow name'] == 'Carbon dioxide, biogenic, uptake'].index)

            # somehow the unit of the radioactive flow Phosphorus-32 in SimaPro is kg and not Bq/kBq. Delete ...
            db = db.drop(db.loc[db.loc[:, 'Elem flow name'] == 'Phosphorus-32'].index)

            # --------------------------- UNIT CONVERSIONS ----------------------------------

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

            # change Becquerels into kiloBecquerels
            db.loc[[i for i in db.index if db.loc[i, 'Elem flow unit'] == 'Bq'], 'CF value'] *= 1000
            db.loc[[i for i in db.index if db.loc[i, 'Elem flow unit'] == 'Bq'], 'Elem flow unit'] = 'kBq'

            # finally, SimaPro limits to 12 characters the size of the CF unit. We rename some to avoid issues
            db.loc[db.loc[:, 'CF unit'] == 'kg CFC-11 eq', 'CF unit'] = 'kg CFC11 eq'
            db.loc[db.loc[:, 'CF unit'] == 'm2 arable land eq', 'CF unit'] = 'm2 ar ld eq'
            db.loc[db.loc[:, 'CF unit'] == 'm2 arable land eq .yr', 'CF unit'] = 'm2 ar ld.yr eq'

            # ------------------------------ NAMING SHENANIGANS ----------------------------

            # the names of CC EQ categories exceed 40 characters (i.e., the limit of SimaPro), so we rename them
            db.loc[db.loc[:, 'Impact category'] ==
                   'Climate change, ecosystem quality, terrestrial ecosystem, short term', 'Impact category'] = \
                'Climate change, EQ, terr, ST'
            db.loc[db.loc[:, 'Impact category'] ==
                   'Climate change, ecosystem quality, terrestrial ecosystem, long term', 'Impact category'] = \
                'Climate change, EQ, terr, LT'
            db.loc[db.loc[:, 'Impact category'] ==
                   'Climate change, ecosystem quality, marine ecosystem, short term', 'Impact category'] = \
                'Climate change, EQ, mar, ST'
            db.loc[db.loc[:, 'Impact category'] ==
                   'Climate change, ecosystem quality, marine ecosystem, long term', 'Impact category'] = \
                'Climate change, EQ, mar, LT'
            # for consistency, also change the human health endpoint categories' names
            db.loc[db.loc[:, 'Impact category'] ==
                   'Climate change, human health, short term', 'Impact category'] = \
                'Climate change, HH, ST'
            db.loc[db.loc[:, 'Impact category'] ==
                   'Climate change, human health, long term', 'Impact category'] = \
                'Climate change, HH, LT'

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
                    'Climate change, ecosystem quality, terrestrial ecosystem, short term'), 'Impact category'] = [
                    'Climate change, EQ, terr, ST' + i.split('Climate change, ecosystem quality, terrestrial ecosystem, short term')[1] for i in
                    db.loc[db.loc[:, 'Impact category'].str.contains(
                        'Climate change, ecosystem quality, terrestrial ecosystem, short term'), 'Impact category']]

                db.loc[db.loc[:, 'Impact category'].str.contains(
                    'Climate change, ecosystem quality, terrestrial ecosystem, long term'), 'Impact category'] = [
                    'Climate change, EQ, terr, LT' + i.split('Climate change, ecosystem quality, terrestrial ecosystem, long term')[1] for i in
                    db.loc[db.loc[:, 'Impact category'].str.contains(
                        'Climate change, ecosystem quality, terrestrial ecosystem, long term'), 'Impact category']]
                db.loc[db.loc[:, 'Impact category'].str.contains(
                    'Climate change, ecosystem quality, marine ecosystem, short term'), 'Impact category'] = [
                    'Climate change, EQ, mar, ST' + i.split('Climate change, ecosystem quality, marine ecosystem, short term')[1] for i in
                    db.loc[db.loc[:, 'Impact category'].str.contains(
                        'Climate change, ecosystem quality, marine ecosystem, short term'), 'Impact category']]

                db.loc[db.loc[:, 'Impact category'].str.contains(
                    'Climate change, ecosystem quality, marine ecosystem, long term'), 'Impact category'] = [
                    'Climate change, EQ, mar, LT' + i.split('Climate change, ecosystem quality, marine ecosystem, long term')[1] for i in
                    db.loc[db.loc[:, 'Impact category'].str.contains(
                        'Climate change, ecosystem quality, marine ecosystem, long term'), 'Impact category']]

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
                __name__, '/Data/mappings/oLCA/v2.5/oLCA_mapping.xlsx'), index_col=0).loc[:,
                   ['Name', 'Name IW+']].dropna()
            differences = olca.loc[olca.Name != olca.loc[:, 'Name IW+']]
            double_iw_flow = olca.loc[olca.loc[:, 'Name IW+'].duplicated(), 'Name IW+'].tolist()

            m = db.merge(differences[['Name IW+', 'Name']], how='left', left_on='Elem flow name', right_on='Name IW+')

            # rows that have any mapping
            matched = m[m['Name'].notna()].copy()
            matched['Elem flow name'] = matched['Name']
            matched = matched.drop(columns=['Name IW+', 'Name'])

            # concat originals + all mapped copies
            db = pd.concat([db, matched], ignore_index=True)

            # remove duplicates names from IW+ (we just added everything together, i.e., oLCA names AND IW names)
            db = db.loc[~((db.loc[:, "Elem flow name"].isin(differences.loc[:, 'Name IW+'])) & ~(
                db.loc[:, "Elem flow name"].isin(double_iw_flow)))]

            # ------------------------------------ UNITS --------------------------------

            # Unit name changes and others
            db.loc[db.loc[:, 'Elem flow unit'] == 'Bq', 'CF value'] *= 1000
            db.loc[db.loc[:, 'Elem flow unit'] == 'Bq', 'Elem flow unit'] = 'kBq'
            db.loc[db.loc[:, 'Elem flow unit'] == 'm2.yr', 'Elem flow unit'] = 'm2*a'
            db.loc[db.loc[:, 'Elem flow unit'] == 'kgy', 'Elem flow unit'] = 'kg*a'

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
            with open(pkg_resources.resource_filename(__name__, '/Data/mappings/oLCA/v2.5/comps.json'), 'r') as f:
                comps = json.load(f)
            db.Compartment = [{v: k for k, v in comps.items()}[i] for i in db.Compartment]

            db.loc[db.loc[:, 'Sub-compartment'] == '(unspecified)', 'Sub-compartment'] = 'unspecified'
            db.loc[db.loc[:, 'Sub-compartment'] == 'groundwater', 'Sub-compartment'] = 'ground water'
            db.loc[
                db.loc[:, 'Sub-compartment'] == 'groundwater, long-term', 'Sub-compartment'] = 'ground water, long-term'
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
                __name__, '/Data/mappings/oLCA/v2.5/all_stressors.xlsx'), index_col=0)

            # split comps and subcomps in two columns for matching with db
            olca_flows['Compartment'] = [i.split('/')[1] for i in olca_flows['comp']]
            olca_flows['Sub-compartment'] = [i.split('/')[2] for i in olca_flows['comp']]
            # only keep relevant columns
            olca_flows = olca_flows.loc[:, ['flow_id', 'flow_name', 'unit', 'Compartment', 'Sub-compartment']]
            # merge with olca_iw, it basically adds the uuids of oLCA
            olca_db = olca_flows.merge(db, left_on=['flow_name', 'unit', 'Compartment', 'Sub-compartment'],
                                       right_on=['Elem flow name', 'Elem flow unit', 'Compartment', 'Sub-compartment'],
                                       how='left')
            # remove flows of IW+ with no link to oLCA flows
            olca_db = olca_db.drop(olca_db.loc[olca_db.loc[:, 'CF value'].isna()].index)
            # remove irrelevant columns
            olca_db = olca_db.drop(['flow_name', 'unit'], axis=1)

            with open(pkg_resources.resource_filename(
                    __name__, '/Data/mappings/oLCA/v2.5/flows_to_spatialize.json'), 'r') as f:
                spatialized_flows = json.load(f)

            # adding spatialized_flows for mapped flow names (e.g., Sulfur oxides based on Sulfur dioxide)
            for spatialized_flow in spatialized_flows:
                # take GB because it's unambiguous
                df = db.loc[db.loc[:, 'Elem flow name'].str.contains(spatialized_flow + ', GB')].copy()
                if df.empty:
                    df = db.loc[db.loc[:, 'Elem flow name'].str.contains(
                        olca.set_index('Name').loc[spatialized_flow, 'Name IW+'] + ', ')].copy()
                    df.loc[:, 'Elem flow name'] = [
                        i.replace(olca.set_index('Name').loc[spatialized_flow, 'Name IW+'], spatialized_flow)
                        for i in df.loc[:, 'Elem flow name']]
                    db = clean_up_dataframe(pd.concat([db, df]))

            # extract all existing locations within IMPACT World+
            regions_water = set([i.split('Water, cooling, unspecified natural origin, ')[1] for i in
                                 self.master_db_carbon_neutrality.loc[:, 'Elem flow name']
                                 if 'Water, cooling, unspecified natural origin, ' in i])
            regions_land = set([i.split('Occupation, annual crops, ')[1] for i in
                                self.master_db_carbon_neutrality.loc[:, 'Elem flow name']
                                if 'Occupation, annual crops, ' in i])
            regions_acid = set(
                [i.split('Nitrogen oxides, ')[1] for i in self.master_db_carbon_neutrality.loc[:, 'Elem flow name']
                 if 'Nitrogen oxides, ' in i])
            regions_pm = set([i.split('Particulates, < 2.5 um, ')[1] for i in
                              self.master_db_carbon_neutrality.loc[:, 'Elem flow name']
                              if 'Particulates, < 2.5 um, ' in i])
            iw_locations = set(list(regions_water) + list(regions_land) + list(regions_acid) + list(regions_pm))

            # separating location from name and put it in a new column / also add flow ids of openLCA from non-spatialized flow
            for spatialized_flow in spatialized_flows:
                df = db.loc[
                    db.loc[:, 'Elem flow name'].isin([spatialized_flow + ', ' + i for i in iw_locations])].copy()

                df.loc[:, 'Location'] = [i.split(spatialized_flow + ', ')[1] for i in df.loc[:, 'Elem flow name']]
                df.loc[:, 'Elem flow name'] = [spatialized_flow for i in df.loc[:, 'Elem flow name']]

                df = olca_flows.merge(df, left_on=['flow_name', 'unit', 'Compartment', 'Sub-compartment'],
                                      right_on=['Elem flow name', 'Elem flow unit', 'Compartment', 'Sub-compartment'],
                                      how='right').dropna(subset=['flow_id']).drop(['flow_name', 'unit'], axis=1)

                olca_db = pd.concat([olca_db, df])

            olca_db = clean_up_dataframe(olca_db)

            return olca_db

        self.olca_iw = linking(self.master_db)
        self.olca_iw_carbon_neutrality = linking(self.master_db_carbon_neutrality)

    def link_to_exiobase(self):
        """
        This method creates an openLCA method with the IW+ characterization factors.
        :return:
        """

        for exio_version in ['3.8', '3.9']:
            EXIO_IW_concordance = pd.read_excel(pkg_resources.resource_filename(
                __name__, 'Data/mappings/exiobase/EXIO_' + exio_version.replace('.', '_') + '_IW_concordance.xlsx'))
            EXIO_IW_concordance.set_index('EXIOBASE', inplace=True)

            C = pd.DataFrame(0, EXIO_IW_concordance.index,
                             list(set(list(zip(self.master_db_carbon_neutrality.loc[:, 'Impact category'],
                                               self.master_db_carbon_neutrality.loc[:, 'CF unit'])))))
            C.columns = pd.MultiIndex.from_tuples(C.columns, names=['Impact category', 'CF unit'])
            C = C.T.sort_index().T

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
                    else:
                        try:
                            CFs = pd.pivot(CF_flow, values='CF value', index=['Compartment', 'Sub-compartment'],
                                           columns=['Impact category', 'CF unit']).loc[('raw', '(unspecified)')].fillna(0)
                        except KeyError:
                            CFs = pd.pivot(CF_flow, values='CF value', index=['Compartment', 'Sub-compartment'],
                                           columns=['Impact category', 'CF unit']).loc[('raw', 'biotic')].fillna(0)
                    CFs.name = flow
                    # dumping the CF values in the C matrix
                    C.update(pd.DataFrame(CFs).T)

            # EXIOBASE land occupation in km2 while IW in m2, so we convert
            C.loc[:, [i for i in C.columns if 'Land' in i[0]]] *= 1000000
            # EXIOBASE water flows in Mm3 while IW in m3, so we convert
            C.loc[:, [i for i in C.columns if 'Water' in i[0]]] *= 1000000
            # more common to have impact categories as index
            C = C.T
            # keep same index format as previously
            C.index = [(i[0] + ' (' + i[1] + ')') for i in C.index]
            # HFC and PFC linked to "carbon dioxide" but should not impact marine acidification
            C.loc[['Marine acidification, long term (PDF.m2.yr)', 'Marine acidification, short term (PDF.m2.yr)',
                   'Resources services loss (adaptation) (MJ)'],
                  ['HFC - air', 'PFC - air']] = 0

            if exio_version == '3.8':
                self.exio_iw_38 = C.copy()
            elif exio_version == '3.9':
                self.exio_iw_39 = C.copy()

    def get_simplified_versions(self):

        # SimaPro
        self.simplified_version_sp = clean_up_dataframe(
            produce_simplified_version(self.iw_sp_carbon_neutrality).reindex(
                self.iw_sp_carbon_neutrality.columns, axis=1))

        # openLCA
        self.simplified_version_olca = clean_up_dataframe(
            produce_simplified_version(self.olca_iw_carbon_neutrality).reindex(
            self.olca_iw_carbon_neutrality.columns, axis=1))

        # ecoinvent
        self.simplified_version_ei310 = clean_up_dataframe(
            produce_simplified_version(self.ei310_iw_carbon_neutrality).reindex(
                self.ei310_iw_carbon_neutrality.columns, axis=1))

        self.simplified_version_ei311 = clean_up_dataframe(
            produce_simplified_version(self.ei311_iw_carbon_neutrality).reindex(
                self.ei311_iw_carbon_neutrality.columns, axis=1))

        self.simplified_version_ei312 = clean_up_dataframe(
            produce_simplified_version(self.ei312_iw_carbon_neutrality).reindex(
                self.ei312_iw_carbon_neutrality.columns, axis=1))

    def get_total_hh_and_eq_for_olca(self):
        """
        OpenLCA doesn't allow for reliable contribution analyses for total damage categories (unlike
        SimaPro). So we create two additional impact categories "Total human health" and "Total ecosystem quality".
        :return:
        """

        total_hh = pd.concat([self.olca_iw.loc[self.olca_iw.Location.isna()].loc[self.olca_iw.loc[:, 'CF unit'] == 'DALY'].drop(
            ['Impact category', 'CAS number', 'MP or Damage'], axis=1).groupby(
            by=['Elem flow name', 'Compartment', 'Sub-compartment',
                'Elem flow unit', 'flow_id', 'CF unit']).agg(
            {'CF value': sum}).reset_index(),
                              self.olca_iw.loc[self.olca_iw.loc[:, 'CF unit'] == 'DALY'].drop(
                                  ['Impact category', 'CAS number', 'MP or Damage',
                                   'Native geographical resolution scale'], axis=1).dropna(subset='Location').groupby(
                                  by=['Elem flow name', 'Compartment', 'Sub-compartment',
                                      'Elem flow unit', 'flow_id',
                                      'Location', 'CF unit']).sum().reset_index()])
        total_hh.loc[:, 'Impact category'] = 'Total human health'

        total_eq = pd.concat([self.olca_iw.loc[self.olca_iw.Location.isna()].loc[self.olca_iw.loc[:, 'CF unit'] == 'PDF.m2.yr'].drop(
            ['Impact category', 'CAS number', 'MP or Damage'], axis=1).groupby(
            by=['Elem flow name', 'Compartment', 'Sub-compartment',
                'Elem flow unit', 'flow_id', 'CF unit']).agg(
            {'CF value': sum}).reset_index(),
                              self.olca_iw.loc[self.olca_iw.loc[:, 'CF unit'] == 'DALY'].drop(
                                  ['Impact category', 'CAS number', 'MP or Damage',
                                   'Native geographical resolution scale'], axis=1).dropna(subset='Location').groupby(
                                  by=['Elem flow name', 'Compartment', 'Sub-compartment',
                                      'Elem flow unit', 'flow_id', 'CF unit',
                                      'Location']).sum().reset_index()])
        total_eq.loc[:, 'Impact category'] = 'Total ecosystem quality'

        self.olca_iw = clean_up_dataframe(pd.concat([self.olca_iw, total_hh, total_eq]))

        total_hh = pd.concat([self.olca_iw_carbon_neutrality.loc[self.olca_iw_carbon_neutrality.Location.isna()].loc[
                                  self.olca_iw_carbon_neutrality.loc[:, 'CF unit'] == 'DALY'].drop(
            ['Impact category', 'CAS number', 'MP or Damage'], axis=1).groupby(
            by=['Elem flow name', 'Compartment', 'Sub-compartment',
                'Elem flow unit', 'CF unit', 'flow_id']).agg(
            {'CF value': sum}).reset_index(),
                              self.olca_iw_carbon_neutrality.loc[self.olca_iw_carbon_neutrality.loc[:, 'CF unit'] == 'DALY'].drop(
                                  ['Impact category', 'CAS number', 'MP or Damage',
                                   'Native geographical resolution scale'], axis=1).dropna(subset='Location').groupby(
                                  by=['Elem flow name', 'Compartment', 'Sub-compartment',
                                      'Elem flow unit', 'flow_id', 'CF unit',
                                      'Location']).sum().reset_index()])
        total_hh.loc[:, 'Impact category'] = 'Total human health'

        total_eq = pd.concat([self.olca_iw_carbon_neutrality.loc[self.olca_iw_carbon_neutrality.Location.isna()].loc[
                                  self.olca_iw_carbon_neutrality.loc[:, 'CF unit'] == 'PDF.m2.yr'].drop(
            ['Impact category', 'CAS number', 'MP or Damage'], axis=1).groupby(
            by=['Elem flow name', 'Compartment', 'Sub-compartment',
                'Elem flow unit', 'CF unit', 'flow_id']).agg(
            {'CF value': sum}).reset_index(),
                              self.olca_iw_carbon_neutrality.loc[self.olca_iw_carbon_neutrality.loc[:, 'CF unit'] == 'DALY'].drop(
                                  ['Impact category', 'CAS number', 'MP or Damage',
                                   'Native geographical resolution scale'], axis=1).dropna(subset='Location').groupby(
                                  by=['Elem flow name', 'Compartment', 'Sub-compartment',
                                      'Elem flow unit', 'flow_id', 'CF unit',
                                      'Location']).sum().reset_index()])
        total_eq.loc[:, 'Impact category'] = 'Total ecosystem quality'

        self.olca_iw_carbon_neutrality = clean_up_dataframe(
            pd.concat([self.olca_iw_carbon_neutrality, total_hh, total_eq]))

# -------------- Support modules -------------------

# taken from the fair==1.6.2 Python package
def meinshausen(
        C,
        Cpi=np.array([277.15, 731.41, 273.87]),
        a1=-2.4785e-07, b1=0.00075906, c1=-0.0021492, d1=5.2488,
        a2=-0.00034197, b2=0.00025455, c2=-0.00024357, d2=0.12173,
        a3=-8.9603e-05, b3=-0.00012462, d3=0.045194,
        F2x=3.71, scale_F2x=True
):
    """Modified Etminan relationship from Meinshausen et al 2019
    https://gmd.copernicus.org/preprints/gmd-2019-222/gmd-2019-222.pdf
    table 3

    Inputs:
        C: [CO2, CH4, N2O] concentrations, [ppm, ppb, ppb]

    Keywords:
        Cpi: pre-industrial [CO2, CH4, N2O] concentrations. Should use defaults
        a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, d3: coefficients
        F2x: radiative forcing from a doubling of CO2.
        scale_F2x: boolean. Scale the calculated value to the specified F2x?

    Returns:
        3-element array of radiative forcing: [F_CO2, F_CH4, F_N2O]

    """
    # Tune the coefficient of CO2 forcing to acheive desired F2x, using
    # pre-industrial CO2 and N2O. F2x_etminan ~= 3.801.
    scaleCO2 = 1
    if scale_F2x:
        F2x_etminan = (
                              -2.4e-7 * Cpi[0] ** 2 + 7.2e-4 * Cpi[0] - 2.1e-4 * Cpi[2] + 5.36) * np.log(2)
        scaleCO2 = F2x / F2x_etminan

    F = np.zeros(3)

    # CO2
    Camax = Cpi[0] - b1 / (2 * a1)
    if Cpi[0] < C[0] <= Camax:  # the most likely case
        alphap = d1 + a1 * (C[0] - Cpi[0]) ** 2 + b1 * (C[0] - Cpi[0])
    elif C[0] <= Cpi[0]:
        alphap = d1
    else:
        alphap = d1 - b1 ** 2 / (4 * a1)
    alphaN2O = c1 * np.sqrt(C[2])
    F[0] = (alphap + alphaN2O) * np.log(C[0] / Cpi[0]) * scaleCO2

    # CH4
    F[1] = (a3 * np.sqrt(C[1]) + b3 * np.sqrt(C[2]) + d3) * (np.sqrt(C[1]) - np.sqrt(Cpi[1]))

    # N2O
    F[2] = (a2 * np.sqrt(C[0]) + b2 * np.sqrt(C[2]) + c2 * np.sqrt(C[1]) + d2) * (np.sqrt(C[2]) - np.sqrt(Cpi[2]))

    return F


M_ATMOS = 5.1352E18
M_AIR = 28.97E-3
M_CO2 = 44.01E-3
M_C = 12.0E-3
M_CH4 = 16.043E-3
M_N2O = 44.0E-3


# function from Official Working Group1 IPCC Github repo: https://github.com/IPCC-WG1/Chapter-7/tree/main/src/ar6/metrics
def co2_analytical(H, d, q, a=np.array([0.2173, 0.2240, 0.2824, 0.2763]), alpha_co2=np.array([0, 394.4, 36.54, 4.304]),
                   co2=409.85, n2o=332.091, co2_ra=0.05):
    """Calculates baseline metrics for a 1 ppm CO2 perturbation.

    Inputs:
    -------
    H : float or `np.ndarray`
        time horizon(s) of interest
    co2 : float, optional
        baseline concentrations of CO2, ppmv
    n2o : float, optional
        baseline concentrations of N2O, ppbv
    co2_ra : float, optional
        tropospheric rapid adjustment enhancement of CO2 forcing, expressed as a decimal
    d : `np.ndarray`, optional
        2-element array of fast and slow timescales to climate warming impulse response function
    q : `np.ndarray`, optional
        2-element array of fast and slow contributions to climate warming impulse response function
    a : `np.ndarray`, optional
        4-element array of partition fractions of CO2 atmospheric boxes, slow to fast
    alpha_co2 : `np.ndarray`, optional
        4-element array of time constants of CO2 atmospheric boxes, slow to fast

    Returns:
    --------
    (rf, agwp, agtp, iagtp) : tuple of float or `np.ndarray`
        rf : Effective radiative forcing from a 1 ppmv increase in CO2
        agwp : Absolute global warming potential of CO2, W m-2 yr kg-1
        agtp : Absolute global temperature change potential of CO2, K kg-1
        iagtp : Integrated absolute global temperature change potential, K kg-1
    """
    # the CH4 concentration does not affect CO2 forcing, so we hardcode an approximate 2019 value
    re = meinshausen(np.array([co2 + 1, 1866.3, n2o]), np.array([co2, 1866.3, n2o]), scale_F2x=False)[0] * (1 + co2_ra)

    ppm2kg = 1E-6 * (M_CO2 / M_AIR) * M_ATMOS
    A = re / ppm2kg  # W/m2/kg

    agtp = H * 0.
    iagtp = H * 0.
    rf = H * 0.
    agwp = H * 0.
    for j in np.arange(2):
        if (j == 0):
            rf = rf + A * a[0]
            agwp = agwp + A * a[0] * H
        agtp = agtp + A * a[0] * q[j] * (1 - np.exp(-H / d[j]))
        iagtp = iagtp + A * a[0] * q[j] * (H - d[j] * (1 - np.exp(-H / d[j])))

        for i in np.arange(1, 4):
            if (j == 0):
                rf = rf + A * a[i] * np.exp(-H / alpha_co2[i])
                agwp = agwp + A * a[i] * alpha_co2[i] * \
                       (1 - np.exp(-H / alpha_co2[i]))
            agtp = agtp + A * a[i] * alpha_co2[i] * q[j] * \
                   (np.exp(-H / alpha_co2[i]) -
                    np.exp(-H / d[j])) / (alpha_co2[i] - d[j])
            iagtp = iagtp + A * a[i] * alpha_co2[i] * q[j] * \
                    (alpha_co2[i] * (1 - np.exp(-H / alpha_co2[i])) -
                     d[j] * (1 - np.exp(-H / d[j]))) / \
                    (alpha_co2[i] - d[j])

    return rf, agwp, agtp, iagtp


# function from Official Working Group1 IPCC Github repo: https://github.com/IPCC-WG1/Chapter-7/tree/main/src/ar6/metrics
def ch4_analytical(H, d, q, co2=409.85, ch4=1866.3275, n2o=332.091, ch4_ra=-0.14, ch4_o3=1.4e-4, ch4_h2o=0.00004,
                   alpha_ch4=11.8):
    """Calculates metrics for a 1 ppb CH4 perturbation.

    Inputs:
    -------
    H : float or `np.ndarray`
        time horizon(s) of interest
    co2 : float, optional
        baseline concentrations of CO2, ppmv
    ch4: float, optional
        baseline concentrations of CH4, ppbv
    n2o : float, optional
        baseline concentrations of N2O, ppbv
    ch4_ra : float, optional
        tropospheric rapid adjustment enhancement of CH4 forcing
    ch4_o3 : float, optional
        radiative efficiency increase of CH4 emissions due to O3 formation, W m-2 (ppb CH4)-1
    ch4_h2o : float, optional
        radiative efficiency increase of CH4 emissions due to stratospheric H2O formation, W m-2 (ppb CH4)-1
    d : `np.ndarray`, optional
        2-element array of fast and slow timescales to climate warming impulse response function
    q : `np.ndarray`, optional
        2-element array of fast and slow contributions to climate warming impulse response function
    alpha_ch4 : float
        perturbation lifetime of CH4, years

    Returns:
    --------
    (rf, agwp, agtp, iagtp) : tuple of float or `np.ndarray`
        rf : Effective radiative forcing from a 1 ppbv increase in CH4
        agwp : Absolute global warming potential of CH4, W m-2 yr kg-1
        agtp : Absolute global temperature change potential of CH4, K kg-1
        iagtp : Integrated absolute global temperature change potential, K kg-1
    """
    re = meinshausen(np.array([co2, ch4 + 1, n2o]), np.array([co2, ch4, n2o]), scale_F2x=False)[1] * (1 + ch4_ra)
    ppb2kg = 1e-9 * (M_CH4 / M_AIR) * M_ATMOS
    A = (re + ch4_o3 + ch4_h2o) / ppb2kg

    agtp = H * 0.
    iagtp = H * 0.
    rf = H * 0.
    agwp = H * 0.

    rf = rf + A * np.exp(-H / (alpha_ch4))
    agwp = agwp + A * alpha_ch4 * (1 - np.exp(-H / alpha_ch4))
    for j in np.arange(2):
        agtp = agtp + A * alpha_ch4 * q[j] * \
               (np.exp(-H / (alpha_ch4)) -
                np.exp(-H / d[j])) / (alpha_ch4 - d[j])
        iagtp = iagtp + A * alpha_ch4 * q[j] * \
                (alpha_ch4 * (1 - np.exp(-H / (alpha_ch4))) -
                 d[j] * (1 - np.exp(-H / d[j]))) / \
                (alpha_ch4 - d[j])
    return rf, agwp, agtp, iagtp


# function from Official Working Group1 IPCC Github repo: https://github.com/IPCC-WG1/Chapter-7/tree/main/src/ar6/metrics
def n2o_analytical(H, d, q, co2=409.85, ch4=1866.3275, n2o=332.091, n2o_ra=0.07, n2o_o3=5.5e-4, f_n2o_ch4=-1.7,
                   ch4_ra=-0.14, ch4_o3=1.4e-4, ch4_h2o=0.00004, alpha_n2o=109):
    """Calculates metrics for a 1 ppb N2O perturbation.

    Inputs:
    -------
    H : float or `np.ndarray`
        time horizon(s) of interest
    co2 : float, optional
        baseline concentrations of CO2, ppmv
    ch4: float, optional
        baseline concentrations of CH4, ppbv
    n2o : float, optional
        baseline concentrations of N2O, ppbv
    n2o_ra : float, optional
        tropospheric rapid adjustment enhancement of N2O forcing
    n2o_o3 : float, optional
        radiative efficiency increase of N2O emissions due to O3 formation, W m-2 (ppb N2O)-1
    f_n2o_ch4 : float, optional
        feedback change in methane lifetime due to N2O emissions, (ppb CH4) (ppb N2O)-1
    ch4_ra : float, optional
        tropospheric rapid adjustment enhancement of CH4 forcing
    ch4_o3 : float, optional
        radiative efficiency increase of CH4 emissions due to O3 formation, W m-2 (ppb CH4)-1
    ch4_h2o : float, optional
        radiative efficiency increase of CH4 emissions due to stratospheric H2O formation, W m-2 (ppb CH4)-1
    d : `np.ndarray`, optional
        2-element array of fast and slow timescales to climate warming impulse response function
    q : `np.ndarray`, optional
        2-element array of fast and slow contributions to climate warming impulse response function
    alpha_n2o : float
        perturbation lifetime of N2O, years

    Returns:
    --------
    (rf, agwp, agtp, iagtp) : tuple of float or `np.ndarray`
        rf : Effective radiative forcing from a 1 ppbv increase in CH4
        agwp : Absolute global warming potential of CH4, W m-2 yr kg-1
        agtp : Absolute global temperature change potential of CH4, K kg-1
        iagtp : Integrated absolute global temperature change potential, K kg-1
    """
    re_n2o = meinshausen(np.array([co2, ch4, n2o + 1]), np.array([co2, ch4, n2o]), scale_F2x=False)[2] * (
            1 + n2o_ra) + n2o_o3
    re_ch4 = meinshausen(np.array([co2, ch4 + 1, n2o]), np.array([co2, ch4, n2o]), scale_F2x=False)[1] * (
            1 + ch4_ra) + ch4_o3 + ch4_h2o
    ppb2kg = 1e-9 * (M_N2O / M_AIR) * M_ATMOS
    # Add in a component for the destruction of methane from AR5 8.SM.11.3.3
    A = (re_n2o + f_n2o_ch4 * re_ch4) / ppb2kg

    agtp = H * 0.
    iagtp = H * 0.
    rf = H * 0.
    agwp = H * 0.
    rf = rf + A * np.exp(-H / alpha_n2o)
    agwp = agwp + A * alpha_n2o * (1 - np.exp(-H / alpha_n2o))
    for j in np.arange(2):
        agtp = agtp + A * alpha_n2o * q[j] * (np.exp(-H / alpha_n2o) -
                                              np.exp(-H / d[j])) / \
               (alpha_n2o - d[j])
        iagtp = iagtp + A * alpha_n2o * q[j] * \
                (alpha_n2o * (1 - np.exp(-H / (alpha_n2o))) -
                 d[j] * (1 - np.exp(-H / d[j]))) / \
                (alpha_n2o - d[j])

    return rf, agwp, agtp, iagtp


# function from Official Working Group1 IPCC Github repo: https://github.com/IPCC-WG1/Chapter-7/tree/main/src/ar6/metrics
def halogen_analytical(H, d, q, alpha, re, mass, halogen_ra=0):
    """Calculates metrics for a 1 ppt perturbation from halogenated gas.

    Inputs:
    -------
    H : float or `np.ndarray`
        time horizon(s) of interest
    alpha : float
        atmospheric lifetime, years
    re : float
        radiative efficiency, W m-2 ppb-1
    mass : float
        molecular mass, kg mol-1
    d : `np.ndarray`, optional
        2-element array of fast and slow timescales to climate warming impulse response function
    q : `np.ndarray`, optional
        2-element array of fast and slow contributions to climate warming impulse response function
    halogen_ra : float, optional
        tropospheric rapid adjustment enhancement of halogen forcing

    Returns:
    --------
    (rf, agwp, agtp, iagtp) : tuple of float or `np.ndarray`
        rf : Effective radiative forcing from a 1 ppbv increase in CH4
        agwp : Absolute global warming potential of CH4, W m-2 yr kg-1
        agtp : Absolute global temperature change potential of CH4, K kg-1
        iagtp : Integrated absolute global temperature change potential, K kg-1
    """
    ppb2kg = 1e-9 * (mass / M_AIR) * M_ATMOS
    A = re / ppb2kg * (1 + halogen_ra)

    agtp = H * 0.
    iagtp = H * 0.
    rf = H * 0.
    agwp = H * 0.
    rf = rf + A * np.exp(-H / alpha)
    agwp = agwp + A * alpha * (1 - np.exp(-H / alpha))
    for j in np.arange(2):
        agtp = agtp + A * alpha * q[j] * (np.exp(-H / alpha) -
                                          np.exp(-H / d[j])) / \
               (alpha - d[j])
        iagtp = iagtp + A * alpha * q[j] * \
                (alpha * (1 - np.exp(-H / (alpha))) - d[j] *
                 (1 - np.exp(-H / d[j]))) / \
                (alpha - d[j])
    return rf, agwp, agtp, iagtp


# function from Official Working Group1 IPCC Github repo: https://github.com/IPCC-WG1/Chapter-7/tree/main/src/ar6/metrics
def carbon_cycle_adjustment(H, d, q, agtp, co2=409.85, n2o=332.091, co2_ra=0.05):
    """Calculates adjustment to metrics based on carbon cycle feedback

    Inputs:
    -------
    H : `np.ndarray` of float
        reguarly spaced time horizons of interest, yr
    agtp : `np.ndarray` of float
        Unadjusted Absolute Global Temperature Change Potential evaluated at each time horizon of H
    co2 : float, optional
        baseline concentrations of CO2, ppmv
    n2o : float, optional
        baseline concentrations of N2O, ppbv
    co2_ra : float, optional
        tropospheric rapid adjustment enhancement of CO2 forcing, expressed as a decimal
    d : `np.ndarray`, optional
        2-element array of fast and slow timescales to climate warming impulse response function
    q : `np.ndarray`, optional
        2-element array of fast and slow contributions to climate warming impulse response function

    Returns:
    --------
    (rf_cc, agwp_cc, agtp_cc) : tuple of `np.ndarray`
        rf_cc : Increase in effective radiative forcing due to carbon cycle adjustment
        agwp : Increase in absolute global warming potential due to carbon cycle adjustment, W m-2 yr kg-1
        agtp : Increase in absolute global temperature change potential due to carbon cycle adjustment, K kg-1
    """
    dts = H[1]
    rf_co2, agwp_co2, agtp_co2, iagtp_co2 = co2_analytical(H, co2=co2, n2o=n2o, co2_ra=co2_ra, d=d, q=q)

    agtp_cc = H * 0.
    agwp_cc = H * 0.
    rf_cc = H * 0.
    F_CO2 = H * 0.
    a = np.array([0.6368, 0.3322, 0.0310])  # Gasser et al. 2017
    alpha = np.array([2.376, 30.14, 490.1])

    gamma = 3.015 * 1E12  # kgCO2/yr/K  Gasser et al. 2017
    r_f = H * 0.
    r_f[0] = np.sum(a) / dts
    for i in np.arange(0, 3):
        r_f = r_f - (a[i] / alpha[i]) * np.exp(-H / alpha[i])

    for j in np.arange(H.size):
        for i in np.arange(j + 1):
            F_CO2[j] = F_CO2[j] + agtp[i] * gamma * r_f[j - i] * dts
    for j in np.arange(H.size):
        for i in np.arange(j + 1):
            rf_cc[j] = rf_cc[j] + F_CO2[i] * rf_co2[j - i] * dts * \
                       (M_CO2 / M_C)
            agwp_cc[j] = agwp_cc[j] + F_CO2[i] * agwp_co2[j - i] * dts * \
                         (M_CO2 / M_C)
            agtp_cc[j] = agtp_cc[j] + F_CO2[i] * agtp_co2[j - i] * dts * \
                         (M_CO2 / M_C)
    return rf_cc, agwp_cc, agtp_cc


def produce_simplified_version(complete_dataframe):
    """
    Method producing the simplified version of IW+ in which there are only 5 indicators:
    - Climate change, short term (=GWP100)
    - Water scarcity
    - Resources services loss
    - Resources services loss (adaptation)
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
                     'Marine eutrophication', 'Ozone layer depletion', 'Particulate matter formation',
                     'Photochemical ozone formation', 'Physical effects on biota', 'Terrestrial acidification']
    # endpoint categories excluded from simplified version
    endpoint_drop = ['Climate change, ecosystem quality, terrestrial ecosystem, long term',
                     'Climate change, ecosystem quality, terrestrial ecosystem, short term',
                     'Climate change, ecosystem quality, marine ecosystem, long term',
                     'Climate change, ecosystem quality, marine ecosystem, short term',
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
    if 'flow_id' in simplified_version.columns:  # for openLCA
        simplified_version = simplified_version.set_index(['CF unit', 'Compartment', 'Sub-compartment',
                                                           'Elem flow name', 'Elem flow unit', 'MP or Damage',
                                                           'flow_id', 'Location'])
    else:  # for other software
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
    if len(hh_simplified.columns) == 8:  # for openLCA
        hh_simplified = hh_simplified.rename(
            columns={'level_0': 'Compartment', 'level_1': 'Sub-compartment', 'level_2': 'Elem flow name',
                     'level_3': 'Elem flow unit', 'level_4': 'MP or Damage', 'level_5': 'flow_id',
                     'level_6': 'Location'})
    else:
        hh_simplified = hh_simplified.rename(
            columns={'level_0': 'Compartment', 'level_1': 'Sub-compartment', 'level_2': 'Elem flow name',
                     'level_3': 'Elem flow unit', 'level_4': 'MP or Damage'})
    hh_simplified.loc[:, 'CF unit'] = 'DALY'
    hh_simplified.loc[:, 'Impact category'] = 'Human health (residual)'
    # make eq_simplified respect the format of self.simplified_version_sp for concatenation
    eq_simplified = eq_simplified.reset_index()
    if len(eq_simplified.columns) == 8:  # for openLCA
        eq_simplified = eq_simplified.rename(
            columns={'level_0': 'Compartment', 'level_1': 'Sub-compartment', 'level_2': 'Elem flow name',
                     'level_3': 'Elem flow unit', 'level_4': 'MP or Damage', 'level_5': 'flow_id',
                     'level_6': 'Location'})
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
                            'Resources services loss'], 'Impact category'] = 'Resources services loss'
    simplified_version.loc[[i for i in simplified_version.index if simplified_version.loc[i, 'Impact category'] ==
                            'Resources services loss adaptation'], 'Impact category'] = 'Resources services loss adaptation'

    return simplified_version


def clean_up_dataframe(df):
    # remove duplicates
    df = df.drop_duplicates()
    # fix index
    df = df.reset_index().drop('index', axis=1)
    return df
