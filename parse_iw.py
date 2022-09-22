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

file name: parse_iw.py
author: Maxime Agez
e-mail: maxime.agez@polymtl.ca
date created: 02-04-22
python version= 3.9
"""

import pyodbc
import pandas as pd
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

class Parse:
    def __init__(self, path_access_db, version, bw2_project=None):
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
            - iw_sp : the dataframe where IW CFs linked to SimaPro elementary flows are stored

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

        self.path_access_db = path_access_db
        self.version = str(version)
        self.bw2_project = bw2_project

        # SQL Alchemy warnings are annoying. Let's remove them from sight!
        warnings.catch_warnings()
        warnings.simplefilter("ignore")

        # OUTPUTs
        self.master_db = pd.DataFrame()
        self.master_db_not_regio = pd.DataFrame()
        self.ei35_iw = pd.DataFrame()
        self.ei36_iw = pd.DataFrame()
        self.ei371_iw = pd.DataFrame()
        self.ei38_iw = pd.DataFrame()
        self.ei35_iw_as_matrix = pd.DataFrame()
        self.ei36_iw_as_matrix = pd.DataFrame()
        self.ei371_iw_as_matrix = pd.DataFrame()
        self.ei38_iw_as_matrix = pd.DataFrame()
        self.iw_sp = pd.DataFrame()
        self.simplified_version_sp = pd.DataFrame()
        self.simplified_version_bw = pd.DataFrame()
        self.sp_data = {}
        self.olca_iw = pd.DataFrame()
        self.olca_data = {}

        # connect to the access database
        self.conn = pyodbc.connect(r'Driver={Microsoft Access Driver (*.mdb, *.accdb)};'
                                   r'DBQ='+self.path_access_db+';')

# -------------------------------------------- Main methods ------------------------------------------------------------

    def load_cfs(self):
        """
        Load the characterization factors and stored them in master_db.
        :return: updated master_db
        """

        print("Loading basic characterization factors...")
        self.load_basic_cfs()
        print("Loading acidification and eutrophication characterization factors...")
        self.load_acid_eutro_cfs()
        print("Loading land use characterization factors...")
        self.load_land_use_cfs()
        print("Loading particulate matter characterization factors...")
        self.load_particulates_cfs()
        print("Loading water scarcity characterization factors...")
        self.load_water_scarcity_cfs()
        print("Loading water availability characterization factors...")
        self.load_water_availability_eq_cfs()
        self.load_water_availability_hh_cfs()
        self.load_water_availability_terr_cfs()
        print("Loading thermally polluted water characterization factors...")
        self.load_thermally_polluted_water_cfs()

        print("Applying rules...")
        self.apply_rules()

        print("Treating regionalized factors...")
        self.create_not_regio_flows()
        self.create_regio_flows_for_not_regio_ic()
        self.order_things_around()
        self.separate_regio_cfs()

        print("Linking to ecoinvent elementary flows...")
        self.link_to_ecoinvent()

        print("Linking to SimaPro elementary flows...")
        self.link_to_sp()

        print("Linking to openLCA elementary flows...")
        self.link_to_olca()

        self.get_simplified_versions()

    def export_to_bw2(self, ei_flows_version=None):
        """
        This method creates a brightway2 method with the IW+ characterization factors.
        :param ei_flows_version: [str] Provide a specific ei version (e.g., 3.6) to be used to determine the elementary flows
                                 to be linked to iw+. Default values = eiv3.8 (in 2022)

        :return:
        """

        print("Exporting to brightway2...")

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
        else:
            ei_in_bw = self.ei38_iw.merge(bw_flows_with_codes)
        ei_in_bw_simple = self.simplified_version_bw.merge(bw_flows_with_codes)

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

            if ei_in_bw_simple.loc[[ic], 'MP or Damage'].iloc[0] == 'Midpoint':
                # create the name of the method
                name = ('IMPACT World+ ' + self.version + ' - Combined midpoint-damage profile', ic[0])
            else:
                # create the name of the method
                if ic[1] == 'DALY':
                    name = ('IMPACT World+ ' + self.version + ' - Combined midpoint-damage profile', ic[0])
                else:
                    name = ('IMPACT World+ ' + self.version + ' - Combined midpoint-damage profile', ic[0])

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

        print("Exporting to SimaPro...")

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
                                    ['', '', '', '', '', ''], ['', '', '', '', '', ''],
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
                                  ['IMPACTWorld+ Damage ' + self.version, '', '', '', '', ''],
                                  ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                  [self.version.split('.')[0],self.version.split('.')[1], '', '', '', ''],
                                  ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                  ['', '', '', '', '', ''], ['', '', '', '', '', ''],
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
                                    ['', '', '', '', '', ''], ['', '', '', '', '', ''],
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
                                    ['IMPACTWorld+ ' + self.version + ' combined midpoint-damage', '', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Version', '', '', '', '', ''],
                                    [self.version.split('.')[0],self.version.split('.')[1], '', '', '', ''],
                                    ['', '', '', '', '', ''], ['Comment', '', '', '', '', ''],
                                    ['', '', '', '', '', ''], ['', '', '', '', '', ''],
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
                      'Particulate matter formation','Photochemical oxidant formation','Terrestrial acidification']
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
                        'combined_values':combined_values, 'simplified_values':simplified_values}

    def export_to_olca(self):
        """
        This method creates the necessary information for the creation of json files in openLCA.
        :return:
        """

        # metadata method (category folder in oLCA app)
        id_category = str(uuid.uuid4())

        category_metadata = {"@context": "http://greendelta.github.io/olca-schema/context.jsonld",
                             "@type": "Category",
                             "@id": id_category,
                             "name": "CIRAIG methods",
                             "version": "0.00.000",
                             "modelType": "IMPACT_METHOD"}

        # metadata json for lcia methods
        dict_ = self.olca_iw.reset_index().loc[:, ['Impact category', 'CF unit']].drop_duplicates().to_dict('list')
        category_names = list(zip(dict_['Impact category'], dict_['CF unit']))
        category_names = [(i[0], i[1], str(uuid.uuid4())) for i in category_names]
        category_names = {(i[0], i[1]): i[2] for i in category_names}

        id_iw = str(uuid.uuid4())

        metadata_iw = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "ImpactMethod",
            "@id": id_iw,
            "name": "IMPACT World+ v2.0",
            "lastChange": "2022-09-15T17:25:43.725-05:00",
            "category": {
                "@type": "Category",
                "@id": id_category,
                "name": "CIRAIG methods",
                "categoryType": "ImpactMethod"},
            'impactCategories': []
        }

        for cat in category_names:
            metadata_iw['impactCategories'].append(
                {"@type": "ImpactCategory",
                 "@id": category_names[cat],
                 "name": cat[0],
                 "refUnit": cat[1]}
            )

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


        cf_dict = {}
        for cat in category_names:
            cf_values = {
                "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
                "@type": "ImpactCategory",
                "@id": category_names[cat],
                "name": cat[0],
                "version": "02.00.000",
                "referenceUnitName": cat[1],
                "impactFactors": []
            }

            dff = self.olca_iw.loc[cat].copy()

            for flow_id in dff.index:
                cf_values["impactFactors"].append({
                    "@type": "ImpactFactor",
                    "value": dff.loc[flow_id, 'CF value'],
                    "flow": {
                        "@type": "Flow",
                        "@id":flow_id,
                        "name": dff.loc[flow_id, 'flow_name'],
                        "categoryPath": ["Elementary flows",
                                         eval(dff.loc[flow_id, 'comp'])[0],
                                         eval(dff.loc[flow_id, 'comp'])[1]],
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
                        "name": dff.loc[flow_id, 'unit'],
                        "categoryPath": ["Technical flow properties"]
                    }
                })

            cf_dict[cat] = cf_values

        # normalization/weighting files
        norm = {}
        for i in category_names:
            if i[1] == 'DALY':
                norm[(i[0], i[1], category_names[i])] = {"normalisationFactor": 13.7, "weightingFactor": 5401.459854}
            elif i[1] == 'PDF.m2.yr':
                norm[(i[0], i[1], category_names[i])] = {"normalisationFactor": 1.01E-4, "weightingFactor": 1386.138614}

        normalization = {
            "@context": "http://greendelta.github.io/olca-schema/context.jsonld",
            "@type": "NwSet",
            "@id": str(uuid.uuid4()),
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

        self.olca_data = {'category_metadata': category_metadata, 'metadata_iw': metadata_iw, 'cf_dict': cf_dict,
                          'normalization':normalization}

    def produce_files(self):
        """
        Function producing the different IW+ files for the different versions.
        :return: the IW+ files
        """

        print("Creating all the files...")

        path = pkg_resources.resource_filename(__name__, '/Databases/Impact_world_' + self.version)

        # if the folders are not there yet, create them
        if not os.path.exists(path + '/Excel/'):
            os.makedirs(path + '/Excel/')
        if not os.path.exists(path + '/DataFrame/'):
            os.makedirs(path + '/DataFrame/')
        if not os.path.exists(path + '/bw2/'):
            os.makedirs(path + '/bw2/')
        if not os.path.exists(path + '/SimaPro/'):
            os.makedirs(path + '/SimaPro/')
        if not os.path.exists(path + '/openLCA/'):
            os.makedirs(path + '/openLCA/')

        # Dev version
        self.master_db.to_excel(path + '/Excel/impact_world_plus_' + self.version + '_dev.xlsx')

        # ecoinvent versions in Excel format
        self.ei35_iw.to_excel(path + '/Excel/impact_world_plus_' + self.version + '_ecoinvent_v35.xlsx')
        self.ei36_iw.to_excel(path + '/Excel/impact_world_plus_' + self.version + '_ecoinvent_v36.xlsx')
        self.ei371_iw.to_excel(path + '/Excel/impact_world_plus_' + self.version + '_ecoinvent_v371.xlsx')
        self.ei38_iw.to_excel(path + '/Excel/impact_world_plus_' + self.version + '_ecoinvent_v38.xlsx')

        # ecoinvent version in DataFrame format
        self.ei35_iw_as_matrix.to_excel(path + '/DataFrame/impact_world_plus_' + self.version + '_ecoinvent_v35.xlsx')
        self.ei36_iw_as_matrix.to_excel(path + '/DataFrame/impact_world_plus_' + self.version + '_ecoinvent_v36.xlsx')
        self.ei371_iw_as_matrix.to_excel(path + '/DataFrame/impact_world_plus_' + self.version + '_ecoinvent_v371.xlsx')
        self.ei38_iw_as_matrix.to_excel(path + '/DataFrame/impact_world_plus_' + self.version + '_ecoinvent_v38.xlsx')

        # brightway2 versions in bw2package format
        IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Combined' not in ic[0])]
        bw2io.package.BW2Package.export_objs(IW_ic, filename='IMPACT_World+_'+self.version, folder=path+'/bw2/')
        # bw2 combined version
        IW_ic = [bw2.Method(ic) for ic in list(bw2.methods) if ('IMPACT World+' in ic[0] and 'Combined' in ic[0])]
        bw2io.package.BW2Package.export_objs(IW_ic, filename='IMPACT_World+_'+self.version+'_combined_version',
                                             folder=path+'/bw2/')

        # SimaPro version in csv format
        with open(path+'/SimaPro/IMPACT_World+_'+self.version+'_Midpoint.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['midpoint_method_metadata'] +
                self.sp_data['midpoint_values'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/IMPACT_World+_'+self.version+'_Damage.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['damage_method_metadata'] +
                self.sp_data['damage_values'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_damage'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/IMPACT_World+_'+self.version+'.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['combined_method_metadata'] +
                self.sp_data['combined_values'] + [['', '', '', '', '', '']])
            writer.writerows(self.sp_data['weighting_info_combined'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])
        with open(path+'/SimaPro/IMPACT_World+_'+self.version+'_combined_midpoint-damage.csv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerows(
                self.sp_data['metadata'] + [['', '', '', '', '', '']] + self.sp_data['simplified_method_metadata'] +
                self.sp_data['simplified_values'] + [['', '', '', '', '', '']])
            writer.writerows([['End', '', '', '', '', '']])

        # create the openLCA version (zip file)
        zipObj = zipfile.ZipFile(path + '/openLCA/CIRAIG_methods.zip', 'w')
        if not os.path.exists(path + '/openLCA/oLCA_folders/categories/'):
            os.makedirs(path + '/openLCA/oLCA_folders/categories/')
        with open(path + '/openLCA/oLCA_folders/categories/' + self.olca_data['category_metadata']['@id'] + '.json',
                  'w') as f:
            json.dump(self.olca_data['category_metadata'], f)
        zipObj.write(path + '/openLCA/oLCA_folders/categories/' + self.olca_data['category_metadata']['@id'] + '.json')
        if not os.path.exists(path + '/openLCA/oLCA_folders/lcia_methods/'):
            os.makedirs(path + '/openLCA/oLCA_folders/lcia_methods/')
        with open(path + '/openLCA/oLCA_folders/lcia_methods/' + self.olca_data['metadata_iw']['@id'] + '.json', 'w') as f:
            json.dump(self.olca_data['metadata_iw'], f)
        zipObj.write(path + '/openLCA/oLCA_folders/lcia_methods/' + self.olca_data['metadata_iw']['@id'] + '.json')
        if not os.path.exists(path + '/openLCA/oLCA_folders/lcia_categories/'):
            os.makedirs(path + '/openLCA/oLCA_folders/lcia_categories/')
        for cat in self.olca_data['cf_dict'].keys():
            with open(path + '/openLCA/oLCA_folders/lcia_categories/' + self.olca_data['cf_dict'][cat]['@id'] + '.json',
                      'w') as f:
                json.dump(self.olca_data['cf_dict'][cat], f)
            zipObj.write(path + '/openLCA/oLCA_folders/lcia_categories/' + self.olca_data['cf_dict'][cat]['@id'] + '.json')
        if not os.path.exists(path + '/openLCA/oLCA_folders/nw_sets/'):
            os.makedirs(path + '/openLCA/oLCA_folders/nw_sets/')
        with open(path + '/openLCA/oLCA_folders/nw_sets/' + self.olca_data['normalization']['@id'] + '.json',
                  'w') as f:
            json.dump(self.olca_data['normalization'], f)
        zipObj.write(path + '/openLCA/oLCA_folders/nw_sets/' + self.olca_data['normalization']['@id'] + '.json')
        zipObj.close()
        # use shutil to simplify the folder structure within the zip file
        shutil.make_archive(path + '/openLCA/CIRAIG_methods', 'zip', path + '/openLCA/oLCA_folders/')

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

        scipy.sparse.save_npz(path+'/for_hybrid_ecoinvent/ei35/Ecoinvent_not_regionalized.npz',
                              scipy.sparse.csr_matrix(self.ei35_iw_as_matrix))
        scipy.sparse.save_npz(path+'/for_hybrid_ecoinvent/ei36/Ecoinvent_not_regionalized.npz',
                              scipy.sparse.csr_matrix(self.ei36_iw_as_matrix))
        scipy.sparse.save_npz(path+'/for_hybrid_ecoinvent/ei371/Ecoinvent_not_regionalized.npz',
                              scipy.sparse.csr_matrix(self.ei371_iw_as_matrix))
        scipy.sparse.save_npz(path+'/for_hybrid_ecoinvent/ei38/Ecoinvent_not_regionalized.npz',
                              scipy.sparse.csr_matrix(self.ei38_iw_as_matrix))

# ----------------------------------------- Secondary methods ----------------------------------------------------------

    def load_basic_cfs(self):
        """
        Loading the basic CFs. By basic we mean that these CFs do not require further treatment.

        Concerned impact categories:
            - Climate change, ecosystem quality, long term
            - Climate change, ecosystem quality, short term
            - Climate change, human health, long term
            - Climate change, human health, short term
            - Climate change, long term
            - Climate change, short term
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
        self.master_db = pd.read_sql(sql='SELECT * FROM [CF - not regionalized - All other impact categories]',
                                     con=self.conn).drop('ID', axis=1)

        # for some reason, some queries are not happening correctly, e.g., "Perfluorooctane, PFC-7-1-18\r\n\r\nPerfluorooctane, PFC-7-1-18"
        # we correct it by splitting at the right place and replacing the value
        read_sql_issues = [i for i in self.master_db.index if '\n' in self.master_db.loc[i, 'Elem flow name']]
        self.master_db.loc[read_sql_issues, 'Elem flow name'] = [i.split("\r")[0] for i in
                                                               self.master_db.loc[read_sql_issues, 'Elem flow name']]

    def load_acid_eutro_cfs(self):
        """
        Loading the CFs for the acidification and eutrophication impact categories. This includes CFs coming from the
        original article of IW+, as well as other CFs extrapolated from stoechiometric ratios.

        Concerned impact categories:
            - Freshwater acidification
            - Terrestrial acidification
            - Freshwater eutrophication
            - Marine eutrophication

        :return: updated master_db
        """

        concerned_cat = [('CF - regionalized - AcidFW - aggregated', 'Freshwater acidification'),
                         ('CF - regionalized - AcidTerr - aggregated', 'Terrestrial acidification'),
                         ('CF - regionalized - EutroFW - aggregated', 'Freshwater eutrophication'),
                         ('CF - regionalized - EutroMar - aggregated', 'Marine eutrophication')]

        stoec_ratios = pd.read_sql(sql='SELECT * FROM [SI - Acid*,Eutro* - Elem flow mapping IW Simapro]',
                                   con=self.conn)

        for cat in concerned_cat:
            # load the corresponding access table
            db = pd.read_sql(sql='SELECT * FROM [' + cat[0] + ']', con=self.conn)
            # converting 3-letters ISO codes to 2-letters ISO codes
            db = convert_country_codes(db)
            # select the impact category within the dataframe
            stoec_ratios_cat = stoec_ratios.loc[
                [i for i in stoec_ratios.index if stoec_ratios.loc[i, 'Impact category'] == cat[1]]]
            # the dataframe where we will reconstruct the database
            df = pd.DataFrame(None, columns=self.master_db.columns)
            for i in stoec_ratios_cat.index:
                # a simple counter to know which id we are at
                id_count = len(self.master_db)
                base_flow = stoec_ratios_cat.loc[i, 'Corresponding elem flow in IW']
                base_values = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == base_flow]]
                new_values = base_values.loc[:, 'Weighted Average'] * stoec_ratios_cat.loc[i, 'Stoechiometric ratio']
                for j in base_values.index:
                    # populating the metadata and values for the added CFs
                    df.loc[id_count, 'Impact category'] = cat[1]
                    df.loc[id_count, 'CF unit'] = base_values.loc[j, 'Unit'].split('[')[1].split(']')[0].split('/')[0]
                    df.loc[id_count, 'Elem flow unit'] = 'kg'
                    df.loc[id_count, 'Compartment'] = base_values.loc[j, 'Compartment']
                    df.loc[id_count, 'Sub-compartment'] = base_values.loc[j, 'Sub-compartment']
                    df.loc[id_count, 'Elem flow name'] = stoec_ratios_cat.loc[i, 'Elem flow name in Simapro'] + ', ' + \
                                                         base_values.loc[j, 'Region code']
                    df.loc[id_count, 'CAS number'] = stoec_ratios_cat.loc[i, 'CAS number']
                    df.loc[id_count, 'MP or Damage'] = base_values.loc[j, 'MP or Damage']
                    df.loc[id_count, 'Native geographical resolution scale'] = base_values.loc[j, 'Resolution']
                    df.loc[id_count, 'CF value'] = new_values.loc[j]
                    id_count += 1
                self.master_db = pd.concat([self.master_db, df])

        self.master_db = clean_up_dataframe(self.master_db)

    def load_land_use_cfs(self):
        """
        Loading the CFs for the land use impact categories.

        Concerned impact categories:
            - Land occupation, biodiversity
            - Land transformation, biodiversity

        :return: updated master_db
        """

        # --------- Land occupation -----------

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - LandOcc - aggregated]', con=self.conn)

        # converting 3-letters ISO codes to 2-letters ISO codes
        db = convert_country_codes(db)

        # the dataframe where we will reconstruct the database
        df = pd.DataFrame(None, columns=self.master_db.columns)

        id_count = len(self.master_db)

        for i in db.index:
            df.loc[id_count, 'Impact category'] = db.loc[i, 'Impact category']
            df.loc[id_count, 'CF unit'] = db.loc[i, 'Unit'].split('[')[1].split(']')[0].split('/')[0]
            df.loc[id_count, 'Elem flow unit'] = db.loc[i, 'Unit'].split('[')[1].split(']')[0].split('/')[1]
            # hardcoded comp and subcomp
            df.loc[id_count, 'Compartment'] = 'Raw'
            df.loc[id_count, 'Sub-compartment'] = 'land'
            df.loc[id_count, 'Elem flow name'] = 'Occupation, ' + db.loc[i, 'Elem flow'].lower() + ', ' + db.loc[i, 'Region code']
            # careful, different name
            df.loc[id_count, 'MP or Damage'] = db.loc[i, 'MP or Damage']
            df.loc[id_count, 'Native geographical resolution scale'] = db.loc[i, 'Resolution']
            df.loc[id_count, 'CF value'] = db.loc[i, 'Weighted Average']
            id_count += 1

        self.master_db = pd.concat([self.master_db, df])

        self.master_db = clean_up_dataframe(self.master_db)

        # --------- Land transformation -----------

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - LandTrans - aggregated]', con=self.conn)

        # converting 3-letters ISO codes to 2-letters ISO codes
        db = convert_country_codes(db)

        # the dataframe where we will reconstruct the database
        df = pd.DataFrame(None, columns=self.master_db.columns)

        id_count = len(self.master_db)

        for i in db.index:
            df.loc[id_count, 'Impact category'] = db.loc[i, 'Impact category']
            df.loc[id_count, 'CF unit'] = db.loc[i, 'Unit'].split('[')[1].split(']')[0].split('/')[0]
            df.loc[id_count, 'Elem flow unit'] = db.loc[i, 'Unit'].split('[')[1].split(']')[0].split('/')[1]
            # hardcoded comp and subcomp
            df.loc[id_count, 'Compartment'] = 'Raw'
            df.loc[id_count, 'Sub-compartment'] = 'land'
            df.loc[id_count, 'Elem flow name'] = 'Transformation, from '+ db.loc[i, 'Elem flow'].lower() + ', ' + \
                                                 db.loc[i, 'Region code']
            # careful, different name
            df.loc[id_count, 'MP or Damage'] = db.loc[i, 'MP or Damage']
            df.loc[id_count, 'Native geographical resolution scale'] = db.loc[i, 'Resolution']
            df.loc[id_count, 'CF value'] = - db.loc[i, 'Weighted Average']
            id_count += 1
            df.loc[id_count, 'Impact category'] = db.loc[i, 'Impact category']
            df.loc[id_count, 'CF unit'] = db.loc[i, 'Unit'].split('[')[1].split(']')[0].split('/')[0]
            df.loc[id_count, 'Elem flow unit'] = db.loc[i, 'Unit'].split('[')[1].split(']')[0].split('/')[1]
            # hardcoded comp and subcomp
            df.loc[id_count, 'Compartment'] = 'Raw'
            df.loc[id_count, 'Sub-compartment'] = 'land'
            df.loc[id_count, 'Elem flow name'] = 'Transformation, to '+ db.loc[i, 'Elem flow'].lower() + ', ' + \
                                                 db.loc[i, 'Region code']
            # careful, different name
            df.loc[id_count, 'MP or Damage'] = db.loc[i, 'MP or Damage']
            df.loc[id_count, 'Native geographical resolution scale'] = db.loc[i, 'Resolution']
            df.loc[id_count, 'CF value'] = db.loc[i, 'Weighted Average']
            id_count += 1

        self.master_db = pd.concat([self.master_db, df])

        self.master_db = clean_up_dataframe(self.master_db)

    def load_particulates_cfs(self):
        """
        Load CFs for particulate mater formation.

        Concerned impact categories:
            - Particulate matter formation

        :return: updated master_db
        """

        # ------- Particulate matter formation ---------

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - PartMatterForm - native]', con=self.conn)

        subcomps = {'Urban': 'high. pop.',
                    'Rural': 'low. pop.',
                    'Remote': 'remote',
                    'Population-weighted average': '(unspecified)'}

        flow_specificity = {'High-stack': 'high-stack',
                            'Low-stack': 'low-stack',
                            'Ground-level': 'ground-level',
                            'Emission-weighted average': 'em. av.'}

        ratio = {'Particulates': 0.6,
                 'Particulates, < 10 um': 0.6,
                 'Particulates, > 10 um': 0,
                 'Particulates, > 2.5 um, and < 10um': 0,
                 'Particulates, < 2.5 um': 1}

        substances = {'SO2': ('Sulfur dioxide', '007446-09-5'),
                      'NH3': ('Ammonia', '007664-41-7'),
                      'NOx': ('Nitrogen oxides', '011104-93-1')}

        def add_generic_pm_intel(master_db, id_count, midpoint_or_damage):
            master_db.loc[id_count, 'Impact category'] = 'Particulate matter formation'
            master_db.loc[id_count, 'Elem flow unit'] = 'kg'
            master_db.loc[id_count, 'Native geographical resolution scale'] = 'Global'
            master_db.loc[id_count, 'Compartment'] = 'Air'
            master_db.loc[id_count, 'MP or Damage'] = midpoint_or_damage
            if midpoint_or_damage == 'Damage':
                master_db.loc[id_count, 'CF unit'] = 'DALY'
            elif midpoint_or_damage == 'Midpoint':
                master_db.loc[id_count, 'CF unit'] = 'kg PM2.5 eq'

        # for particulates
        for flow in ratio:
            for subcomp in subcomps:
                for spec in flow_specificity:
                    # Midpoint
                    id_count = len(self.master_db)
                    add_generic_pm_intel(self.master_db, id_count, 'Midpoint')
                    self.master_db.loc[id_count, 'Elem flow name'] = flow + ', ' + flow_specificity[spec]
                    self.master_db.loc[id_count, 'Sub-compartment'] = subcomps[subcomp]
                    self.master_db.loc[id_count, 'CF value'] = (
                            db.loc[[i for i in db.index if (db.loc[i, 'Elem flow'] == 'Primary PM2.5' and
                                                            db.loc[i, 'Archetype 1'] == subcomp and
                                                            db.loc[i, 'Archetype 2'] == spec and
                                                            db.loc[i, 'MP or Damage'] == 'Midpoint'
                                                            )], 'CF value'].iloc[0] * ratio[flow])
                    # Endpoint
                    id_count += 1
                    add_generic_pm_intel(self.master_db, id_count, 'Damage')
                    self.master_db.loc[id_count, 'Elem flow name'] = flow + ', ' + flow_specificity[spec]
                    self.master_db.loc[id_count, 'Sub-compartment'] = subcomps[subcomp]
                    self.master_db.loc[id_count, 'CF value'] = (
                            db.loc[[i for i in db.index if (db.loc[i, 'Elem flow'] == 'Primary PM2.5' and
                                                            db.loc[i, 'Archetype 1'] == subcomp and
                                                            db.loc[i, 'Archetype 2'] == spec and
                                                            db.loc[i, 'MP or Damage'] == 'Damage'
                                                            )], 'CF value'].iloc[0] * ratio[flow])
                # add the flow itself, i.e., add "Particulates" and not "Particulates, high-stack"
                id_count = len(self.master_db)
                add_generic_pm_intel(self.master_db, id_count, 'Midpoint')
                self.master_db.loc[id_count, 'Elem flow name'] = flow
                self.master_db.loc[id_count, 'Sub-compartment'] = subcomps[subcomp]
                # default value is taken with "Emission-weighted average"
                self.master_db.loc[id_count, 'CF value'] = (
                        db.loc[[i for i in db.index if (db.loc[i, 'Elem flow'] == 'Primary PM2.5' and
                                                        db.loc[i, 'Archetype 1'] == subcomp and
                                                        db.loc[i, 'Archetype 2'] == 'Emission-weighted average' and
                                                        db.loc[i, 'MP or Damage'] == 'Midpoint'
                                                        )], 'CF value'].iloc[0] * ratio[flow])
                id_count += 1
                add_generic_pm_intel(self.master_db, id_count, 'Damage')
                self.master_db.loc[id_count, 'Elem flow name'] = flow
                self.master_db.loc[id_count, 'Sub-compartment'] = subcomps[subcomp]
                # default value is taken with "Emission-weighted average"
                self.master_db.loc[id_count, 'CF value'] = (
                        db.loc[[i for i in db.index if (db.loc[i, 'Elem flow'] == 'Primary PM2.5' and
                                                        db.loc[i, 'Archetype 1'] == subcomp and
                                                        db.loc[i, 'Archetype 2'] == 'Emission-weighted average' and
                                                        db.loc[i, 'MP or Damage'] == 'Damage'
                                                        )], 'CF value'].iloc[0] * ratio[flow])
        # for SO2, NH3 et NOx
        for subs in substances:
            for subcomp in subcomps:
                id_count = len(self.master_db)
                add_generic_pm_intel(self.master_db, id_count, 'Midpoint')
                self.master_db.loc[id_count, 'Elem flow name'] = substances[subs][0]
                self.master_db.loc[id_count, 'CAS number'] = substances[subs][1]
                self.master_db.loc[id_count, 'Sub-compartment'] = subcomps[subcomp]
                self.master_db.loc[id_count, 'CF value'] = (
                        db.loc[[i for i in db.index if (db.loc[i, 'Elem flow'] == subs and
                                                        db.loc[i, 'Archetype 1'] == subcomp and
                                                        db.loc[i, 'MP or Damage'] == 'Midpoint'
                                                        )], 'CF value'].iloc[0])
                id_count += 1
                add_generic_pm_intel(self.master_db, id_count, 'Damage')
                self.master_db.loc[id_count, 'Elem flow name'] = substances[subs][0]
                self.master_db.loc[id_count, 'CAS number'] = substances[subs][1]
                self.master_db.loc[id_count, 'Sub-compartment'] = subcomps[subcomp]
                self.master_db.loc[id_count, 'CF value'] = (
                    db.loc[[i for i in db.index if (db.loc[i, 'Elem flow'] == subs and
                                                    db.loc[i, 'Archetype 1'] == subcomp and
                                                    db.loc[i, 'MP or Damage'] == 'Damage'
                                                    )], 'CF value'].iloc[0])

        # add "unspecified" emissions, i.e., not high/low stack, not em. av. and not ground-level
        df = self.master_db.loc[[i for i in self.master_db.index if
                                 'em. av.' in self.master_db.loc[i, 'Elem flow name']]].copy()
        df.loc[:, 'Elem flow name'] = [i.split(', em. av.')[0] for i in df.loc[:, 'Elem flow name']]
        self.master_db = pd.concat([self.master_db, df])

        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_scarcity_cfs(self):
        """
        Load CFs for water scarcity

        Concerned impact categories:
            - Water scarcity

        :return: updated master_db
        """

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - WaterScarc - aggregated]', con=self.conn)

        # converting 3-letters ISO codes to 2-letters ISO codes
        db = convert_country_codes(db)

        def add_generic_scarcity_intel(master_db, id_count):
            master_db.loc[id_count, 'Impact category'] = 'Water scarcity'
            master_db.loc[id_count, 'Elem flow unit'] = 'm3'
            master_db.loc[id_count, 'MP or Damage'] = 'Midpoint'
            master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
            master_db.loc[id_count, 'CF unit'] = 'm3 world-eq'
            master_db.loc[id_count, 'CAS number'] = '007732-18-5'

        for flow in ['UNKNOWN', 'AGRI', 'NON-AGRI']:
            data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == flow]]
            for j in data.index:
                # Water comp
                id_count = len(self.master_db)
                add_generic_scarcity_intel(self.master_db, id_count)
                self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[j, 'Resolution']
                self.master_db.loc[id_count, 'Compartment'] = 'Water'
                self.master_db.loc[id_count, 'CF value'] = - data.loc[j, 'Weighted Average']
                if flow == 'UNKNOWN':
                    self.master_db.loc[id_count, 'Elem flow name'] = 'Water' + ', ' + data.loc[j, 'Region code']
                else:
                    self.master_db.loc[id_count, 'Elem flow name'] = 'Water, ' + flow.lower() + ', ' + data.loc[
                        j, 'Region code']
                # Raw comp
                id_count += 1
                add_generic_scarcity_intel(self.master_db, id_count)
                self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[j, 'Resolution']
                self.master_db.loc[id_count, 'Compartment'] = 'Raw'
                self.master_db.loc[id_count, 'CF value'] = data.loc[j, 'Weighted Average']
                if flow == 'UNKNOWN':
                    self.master_db.loc[id_count, 'Elem flow name'] = 'Water' + ', ' + data.loc[j, 'Region code']
                else:
                    self.master_db.loc[id_count, 'Elem flow name'] = 'Water, ' + flow.lower() + ', ' + data.loc[
                        j, 'Region code']

        self.master_db = clean_up_dataframe(self.master_db)

        # creating the other water flows from the default water flow
        other_water = ['Water, lake', 'Water, river', 'Water, unspecified natural origin',
                       'Water, well, in ground', 'Water, cooling, unspecified natural origin']
        for water in other_water:
            df = self.master_db.loc[
                [i for i in self.master_db.index if (self.master_db.loc[i, 'Impact category'] == 'Water scarcity' and
                                                     len(self.master_db.loc[i, 'Elem flow name'].split(',')) == 2 and
                                                     'agri' not in self.master_db.loc[i, 'Elem flow name'] and
                                                     self.master_db.loc[i,'Compartment'] == 'Raw')]]
            df.loc[:, 'Elem flow name'] = [water + ', ' + i.split(', ')[1] for i in df.loc[:, 'Elem flow name']]
            self.master_db = pd.concat([self.master_db, df])
            self.master_db = clean_up_dataframe(self.master_db)

        self.master_db = clean_up_dataframe(self.master_db)

        # add the RoW geography based on the Global value
        df = self.master_db.loc[
            [i for i in self.master_db.index if (self.master_db.loc[i, 'Impact category'] == 'Water scarcity' and
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

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - WaterAvailab_EcosysFW - native]', con=self.conn)

        # Proper aggregation to determine the global value was not yet performed. Median is used instead.
        def add_generic_water_avail_eq_intel(master_db, id_count):
            master_db.loc[id_count, 'Impact category'] = 'Water availability, freshwater ecosystem'
            master_db.loc[id_count, 'Native geographical resolution scale'] = 'Country'
            master_db.loc[id_count, 'MP or Damage'] = 'Damage'
            master_db.loc[id_count, 'CAS number'] = '007732-18-5'
            master_db.loc[id_count, 'CF unit'] = 'PDF.m2.yr'
            master_db.loc[id_count, 'Elem flow unit'] = 'm3'
            master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'

        with open(pkg_resources.resource_filename(__name__, '/Data/geographies_water.json'), 'r') as f:
            regions = json.load(f)

        # creating the other water flows from the default water flow
        other_water = ['Water, lake', 'Water, river', 'Water, unspecified natural origin',
                       'Water, well, in ground', 'Water, cooling, unspecified natural origin']

        for region in regions:
            # Water comp
            id_count = len(self.master_db)
            add_generic_water_avail_eq_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water, ' + region
            self.master_db.loc[id_count, 'Compartment'] = 'Water'
            self.master_db.loc[id_count, 'CF value'] = - db.loc[:, 'CF value'].median()
            for water in other_water:
                # Raw comp
                id_count = len(self.master_db)
                add_generic_water_avail_eq_intel(self.master_db, id_count)
                self.master_db.loc[id_count, 'Elem flow name'] = water + ', ' + region
                self.master_db.loc[id_count, 'Compartment'] = 'Raw'
                self.master_db.loc[id_count, 'CF value'] = db.loc[:, 'CF value'].median()

        # changes value of resolution scale for the global flow
        self.master_db.loc[[i for i in self.master_db.index if ', GLO' in self.master_db.loc[i,'Elem flow name']],
                           'Native geographical resolution scale'] = 'Global'

        self.master_db = clean_up_dataframe(self.master_db)

        # add the RoW geography based on the Global value
        df = self.master_db.loc[
            [i for i in self.master_db.index if (self.master_db.loc[i, 'Impact category'] ==
                                                 'Water availability, freshwater ecosystem' and
                                                 'GLO' in self.master_db.loc[i, 'Elem flow name'])]]
        df['Elem flow name'] = [i.split(', GLO')[0] + ', RoW' for i in df['Elem flow name']]
        df['Native geographical resolution scale'] = 'Other region'
        self.master_db = pd.concat([self.master_db, df])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_availability_hh_cfs(self):
        """
        Load CFs for water availability human health.

        Concerned impact categories:
            - Water availability, human health

        :return: update master_db
        """

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - WaterAvailab_HH - aggregated]', con=self.conn)

        # converting 3-letters ISO codes to 2-letters ISO codes
        db = convert_country_codes(db)

        # remove the Iraq-Saudi Arabia Neutral Zone geography which creates issues with country_converter because it yields a list and not a string
        db = clean_up_dataframe(db.drop([i for i in db.index if type(db.loc[i, 'Region code']) != str]))

        def add_generic_water_avai_hh_intel(master_db, id_count):
            master_db.loc[id_count, 'Impact category'] = 'Water availability, human health'
            master_db.loc[id_count, 'Elem flow unit'] = 'm3'
            master_db.loc[id_count, 'MP or Damage'] = 'Damage'
            master_db.loc[id_count, 'CF unit'] = 'DALY'

        # Water comp / unspecified subcomp
        data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == 'Unknown']]
        for i in data.index:
            id_count = len(self.master_db)
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[i, 'Resolution']
            self.master_db.loc[id_count, 'Compartment'] = 'Water'
            self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water' + ', ' + data.loc[i, 'Region code']
            self.master_db.loc[id_count, 'CF value'] = -data.loc[i, 'Weighted Average']

        # Water comp / lake subcomp
        data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == 'Surface']]
        for i in data.index:
            id_count = len(self.master_db)
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[i, 'Resolution']
            self.master_db.loc[id_count, 'Compartment'] = 'Water'
            self.master_db.loc[id_count, 'Sub-compartment'] = 'lake'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water' + ', ' + data.loc[i, 'Region code']
            self.master_db.loc[id_count, 'CF value'] = -data.loc[i, 'Weighted Average']

        # Water comp / river subcomp
        data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == 'Surface']]
        for i in data.index:
            id_count = len(self.master_db)
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[i, 'Resolution']
            self.master_db.loc[id_count, 'Compartment'] = 'Water'
            self.master_db.loc[id_count, 'Sub-compartment'] = 'river'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water' + ', ' + data.loc[i, 'Region code']
            self.master_db.loc[id_count, 'CF value'] = -data.loc[i, 'Weighted Average']

        # Water comp / groundwater subcomp
        data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == 'Ground']]
        for i in data.index:
            id_count = len(self.master_db)
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[i, 'Resolution']
            self.master_db.loc[id_count, 'Compartment'] = 'Water'
            self.master_db.loc[id_count, 'Sub-compartment'] = 'groundwater'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water' + ', ' + data.loc[i, 'Region code']
            self.master_db.loc[id_count, 'CF value'] = -data.loc[i, 'Weighted Average']

        # Raw comp / unspecified water
        data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == 'Unknown']]
        for i in data.index:
            id_count = len(self.master_db)
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[i, 'Resolution']
            self.master_db.loc[id_count, 'Compartment'] = 'Raw'
            self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water, unspecified natural origin' + ', ' + data.loc[
                i, 'Region code']
            self.master_db.loc[id_count, 'CF value'] = data.loc[i, 'Weighted Average']

        # Raw comp / lake water
        data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == 'Surface']]
        for i in data.index:
            id_count = len(self.master_db)
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[i, 'Resolution']
            self.master_db.loc[id_count, 'Compartment'] = 'Raw'
            self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water, lake' + ', ' + data.loc[i, 'Region code']
            self.master_db.loc[id_count, 'CF value'] = data.loc[i, 'Weighted Average']

        # Raw comp / river water
        data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == 'Surface']]
        for i in data.index:
            id_count = len(self.master_db)
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[i, 'Resolution']
            self.master_db.loc[id_count, 'Compartment'] = 'Raw'
            self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water, river' + ', ' + data.loc[i, 'Region code']
            self.master_db.loc[id_count, 'CF value'] = data.loc[i, 'Weighted Average']

        # Raw comp / cooling, unspecified natural origin water
        data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == 'Unknown']]
        for i in data.index:
            id_count = len(self.master_db)
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[i, 'Resolution']
            self.master_db.loc[id_count, 'Compartment'] = 'Raw'
            self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water, cooling, unspecified natural origin' + ', ' + data.loc[i, 'Region code']
            self.master_db.loc[id_count, 'CF value'] = data.loc[i, 'Weighted Average']

        # Raw comp / well water
        data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow'] == 'Ground']]
        for i in data.index:
            id_count = len(self.master_db)
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Native geographical resolution scale'] = data.loc[i, 'Resolution']
            self.master_db.loc[id_count, 'Compartment'] = 'Raw'
            self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water, well, in ground' + ', ' + data.loc[
                i, 'Region code']
            self.master_db.loc[id_count, 'CF value'] = data.loc[i, 'Weighted Average']

        # drop NaN values
        self.master_db.dropna(subset=['CF value'], inplace=True)

        # add the RoW geography based on the Global value
        df = self.master_db.loc[
            [i for i in self.master_db.index if (self.master_db.loc[i, 'Impact category'] ==
                                                 'Water availability, human health' and
                                                 'GLO' in self.master_db.loc[i, 'Elem flow name'])]]
        df['Elem flow name'] = [i.split(', GLO')[0] + ', RoW' for i in df['Elem flow name']]
        df['Native geographical resolution scale'] = 'Other region'
        self.master_db = pd.concat([self.master_db, df])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_availability_terr_cfs(self):
        """
        Load CFs for water availability terrestrial ecosystem.

        Concerned impact categories:
            - Water availability, terrestrial ecosystem

        :return: update master_db
        """

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - WaterAvailab_Terr - aggregated]', con=self.conn)

        self.master_db = pd.concat([self.master_db, db.drop('ID', axis=1)])

        self.master_db = clean_up_dataframe(self.master_db)

        # add the RoW geography based on the Global value
        df = self.master_db.loc[
            [i for i in self.master_db.index if (self.master_db.loc[i, 'Impact category'] ==
                                                 'Water availability, terrestrial ecosystem' and
                                                 'GLO' in self.master_db.loc[i, 'Elem flow name'])]]
        df['Elem flow name'] = [i.split(', GLO')[0] + ', RoW' for i in df['Elem flow name']]
        df['Native geographical resolution scale'] = 'Other region'
        self.master_db = pd.concat([self.master_db, df])
        self.master_db = clean_up_dataframe(self.master_db)

    def load_thermally_polluted_water_cfs(self):
        """
        Load CFs for thermally polluted water.

        Concerned impact categories:
            - Thermally polluted water

        :return: update master_db
        """

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - ThermallyPollutedWater - aggregated]', con=self.conn)

        self.master_db = pd.concat([self.master_db, db.drop('ID', axis=1)])

        self.master_db = clean_up_dataframe(self.master_db)

        # add the RoW geography based on the Global value
        df = self.master_db.loc[
            [i for i in self.master_db.index if (self.master_db.loc[i, 'Impact category'] ==
                                                 'Thermally polluted water' and
                                                 'GLO' in self.master_db.loc[i, 'Elem flow name'])]]
        df['Elem flow name'] = [i.split(', GLO')[0] + ', RoW' for i in df['Elem flow name']]
        df['Native geographical resolution scale'] = 'Other region'
        self.master_db = pd.concat([self.master_db, df])
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
                    'Soil': ['industrial','agricultural'],
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
                                            'Freshwater ecotoxicity','Human toxicity cancer',
                                            'Human toxicity non-cancer','Marine acidification']

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

        regio_flows = set([', '.join(i.split(', ')[:-1]) for i
                           in self.master_db[self.master_db['Native geographical resolution scale'] == 'Country'].loc[
                                                                  :, 'Elem flow name']])

        regio_regions = set([i.split(', ')[-1] for i
                             in self.master_db[self.master_db['Native geographical resolution scale'] == 'Country'].loc[:,
                                                        'Elem flow name']])

        regio_ic = set(
            self.master_db[self.master_db['Native geographical resolution scale'] == 'Country'].loc[:, 'Impact category'])

        # identify which regionalized flows need to be characterized for non-regionalized imapct categories

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

        versions_ei = ['3.5', '3.6', '3.7.1', '3.8']

        for version_ei in versions_ei:
            print("Linking to ecoinvent"+str(version_ei)+" elementary flows ...")
            ei_iw_db = self.master_db_not_regio.copy()

            # -------------- Mapping substances --------------

            ei_mapping = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/mappings/ei'+
                                                                     version_ei.replace('.','')+
                                                                     '/ei_iw_mapping.xlsx'))
            ei_mapping = ei_mapping.loc[:, ['ecoinvent name', 'iw name']].dropna()
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
                                                                     version_ei.replace('.','')+
                                                                     "/comps.json"), "r") as f:
                comps = json.load(f)
            with open(pkg_resources.resource_filename(__name__, "Data/mappings/ei"+
                                                                     version_ei.replace('.','')+
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

            if version_ei == '3.5':
                self.ei35_iw = ei_iw_db
            elif version_ei == '3.6':
                self.ei36_iw = ei_iw_db
            elif version_ei == '3.7.1':
                self.ei371_iw = ei_iw_db
            elif version_ei == '3.8':
                self.ei38_iw = ei_iw_db

            # introducing UUID for stressors of ecoinvent
            stressors_ei = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/metadata/ei'+
                                                                     version_ei.replace('.','')+
                                                                     '/stressors.xlsx'))
            # matching with uuids
            df_ei = stressors_ei.set_index(['name', 'unit', 'comp', 'subcomp']).drop('cas', axis=1)
            df_ei.index.names = (None, None, None, None)
            df_iw = ei_iw_db.set_index(['Elem flow name', 'Elem flow unit', 'Compartment', 'Sub-compartment'])
            df_iw.index.names = (None, None, None, None)
            ei_iw_db = df_ei.join(df_iw).dropna(subset=['Impact category', 'CF value'])
            ei_iw_db = ei_iw_db.set_index('id').drop(['CAS number', 'MP or Damage', 'Native geographical resolution scale'],
                                               axis=1)
            # changing into a dataframe format readily available for matrix calculations
            ei_iw_db = ei_iw_db.pivot_table(values='CF value', index='id', columns=['Impact category', 'CF unit']).fillna(0)
            ei_iw_db = ei_iw_db.reindex(stressors_ei.id).fillna(0)

            if version_ei == '3.5':
                self.ei35_iw_as_matrix = ei_iw_db.T
            elif version_ei == '3.6':
                self.ei36_iw_as_matrix = ei_iw_db.T
            elif version_ei == '3.7.1':
                self.ei371_iw_as_matrix = ei_iw_db.T
            elif version_ei == '3.8':
                self.ei38_iw_as_matrix = ei_iw_db.T

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
        sp = pd.read_excel(pkg_resources.resource_filename(__name__, '/Data/mappings/SP/sp_mapping.xlsx')).dropna()
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

        olca_mapping = pd.read_excel(pkg_resources.resource_filename(__name__, '/Data/mappings/oLCA/oLCA_mapping.xlsx'))
        olca_flows = pd.read_excel(pkg_resources.resource_filename(__name__, '/Data/mappings/oLCA/all_stressors.xlsx'))

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
                                      right_on=['oLCA name'], how='left').dropna(subset='iw name')
        # split comp and subcomp
        olca_flows['Sub-compartment'] = [eval(i)[1] for i in olca_flows['comp']]
        olca_flows['Compartment'] = [eval(i)[0] for i in olca_flows['comp']]

        # map compartments between oLCA and IW
        with open(pkg_resources.resource_filename(__name__, '/Data/mappings/oLCA/comps.json'), 'r') as f:
            comps = json.load(f)
        olca_flows['Compartment'] = [comps[i] for i in olca_flows['Compartment']]

        # map sub-compartments between oLCA and IW
        with open(pkg_resources.resource_filename(__name__, '/Data/mappings/oLCA/subcomps.json'), 'r') as f:
            subcomps = json.load(f)

        # remove weird LT subcomps that only exist in openLCA flows
        olca_flows = olca_flows.drop([i for i in olca_flows.index if eval(olca_flows.loc[i, 'comp'])[1] in [
            'river, long-term', 'fresh water, long-term']])
        olca_flows['Sub-compartment'] = [subcomps[i] for i in olca_flows['Sub-compartment']]

        # switch iw names for olca names
        self.olca_iw = olca_flows.merge(self.olca_iw, how='inner',
                                        left_on=['iw name', 'Sub-compartment', 'Compartment'],
                                        right_on=['Elem flow name', 'Sub-compartment', 'Compartment'])

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
                          'Ionizing radiation, human health'], 'CF value'] /= 1000

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
            df = df.loc[df['comp'] == converting.loc[i, 'comp']]
            df = df.reset_index()
            df['flow_id'] = converting.loc[i, 'flow_id']
            df['flow_name'] = converting.loc[i, 'flow_name']
            df = df.set_index(['Impact category', 'CF unit', 'flow_id'])
            self.olca_iw = pd.concat([self.olca_iw, df])

    def get_simplified_versions(self, ei_flows_version=None):

        self.simplified_version_sp = clean_up_dataframe(produce_simplified_version(self.iw_sp).reindex(
            self.iw_sp.columns, axis=1))

        if ei_flows_version == '3.5':
            self.simplified_version_bw = clean_up_dataframe(produce_simplified_version(self.ei35_iw).reindex(
                self.ei35_iw.columns, axis=1))
        elif ei_flows_version == '3.6':
            self.simplified_version_bw = clean_up_dataframe(produce_simplified_version(self.ei36_iw).reindex(
                self.ei36_iw.columns, axis=1))
        elif ei_flows_version == '3.7.1':
            self.simplified_version_bw = clean_up_dataframe(produce_simplified_version(self.ei371_iw).reindex(
                self.ei371_iw.columns, axis=1))
        elif ei_flows_version == '3.8':
            self.simplified_version_bw = clean_up_dataframe(produce_simplified_version(self.ei38_iw).reindex(
                self.ei38_iw.columns, axis=1))
        else:
            self.simplified_version_bw = clean_up_dataframe(produce_simplified_version(self.ei38_iw).reindex(
                self.ei38_iw.columns, axis=1))

# -------------- Support modules -------------------

def produce_simplified_version(complete_dataframe):
    """
    Method producing the simplified version of IW+ in which there are only 6 indicators:
    - Climate change, short term (=GWP100)
    - Water scarcity
    - Fossil resources
    - Mineral resources
    - Remaining Human health damage
    - Remaining Ecosystem quality damage
    The damage impacts of climate change and water use on the Human health and Ecosystem quality area of protections
    are removed, hence only leaving "Remaining [...] damage". These are removed to avoid potential misuse by users
    where the impact of climate change and water use would be double-counted. Once as a midpoint category, and
    another time as a damage category.
    :return:
    """

    simplified_version = complete_dataframe.copy()

    # midpoint categories excluded from simplified version
    midpoint_drop = ['Climate change, long term', 'Freshwater acidification', 'Freshwater ecotoxicity',
                     'Freshwater eutrophication',
                     'Human toxicity cancer', 'Human toxicity non-cancer', 'Ionizing radiations',
                     'Land occupation, biodiversity',
                     'Land transformation, biodiversity', 'Marine eutrophication', 'Ozone layer depletion',
                     'Particulate matter formation', 'Photochemical oxidant formation', 'Terrestrial acidification']
    # endpoint categories excluded from simplified version
    endpoint_drop = ['Climate change, ecosystem quality, long term',
                     'Climate change, ecosystem quality, short term',
                     'Climate change, human health, long term', 'Climate change, human health, short term',
                     'Water availability, freshwater ecosystem', 'Water availability, human health',
                     'Water availability, terrestrial ecosystem']
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

    return simplified_version


def clean_up_dataframe(df):
    # remove duplicates
    df = df.drop_duplicates()
    # fix index
    df = df.reset_index().drop('index',axis=1)
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


def extract_olca_flows():
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
