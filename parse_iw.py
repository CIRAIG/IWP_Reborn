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
import scipy.sparse

class Parse:
    def __init__(self, path_access_db, version):
        """
        :param path_access_db: path to the Microsoft access database (source version)
        :param version: the version of IW+ to parse

        Object instance variables:
        -------------------------
            - master_db : the master dataframe where CFs will be stored
            - conn : the connector to the access database
        """

        self.path_access_db = path_access_db
        self.version = str(version)

        # OUTPUTs
        self.master_db = pd.DataFrame()
        self.ei35_iw = pd.DataFrame()
        self.ei36_iw = pd.DataFrame()
        self.ei371_iw = pd.DataFrame()
        self.ei38_iw = pd.DataFrame()
        self.ei35_iw_as_matrix = pd.DataFrame()
        self.ei36_iw_as_matrix = pd.DataFrame()
        self.ei371_iw_as_matrix = pd.DataFrame()
        self.ei38_iw_as_matrix = pd.DataFrame()

        # connect to the access database
        self.conn = pyodbc.connect(r'Driver={Microsoft Access Driver (*.mdb, *.accdb)};'
                                   r'DBQ='+self.path_access_db+';')

    def load_cfs(self):
        """
        Load the characterization factors and stored them in master_db.
        :return: updated master_db
        """

        self.load_basic_cfs()
        self.load_acid_eutro_cfs()
        self.load_land_use_cfs()
        self.load_particulates_cfs()
        self.load_water_scarcity_cfs()
        self.load_water_availability_eq_cfs()
        self.load_water_availability_hh_cfs()
        self.load_water_availability_terr_cfs()
        self.load_thermally_polluted_water_cfs()

        self.apply_rules()

        self.master_db = self.master_db.sort_values(by=['Impact category', 'Elem flow name'])
        self.master_db = self.master_db.reset_index().drop('index', axis=1)

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
            # only keep global average values
            db = db.loc[[i for i in db.index if db.loc[i, 'Resolution'] in ['Global', 'Not regionalized']]]
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
                    df.loc[id_count, 'Elem flow name'] = stoec_ratios_cat.loc[i, 'Elem flow name in Simapro']
                    df.loc[id_count, 'CAS number'] = stoec_ratios_cat.loc[i, 'CAS number']
                    df.loc[id_count, 'MP or Damage'] = base_values.loc[j, 'MP or Damage']
                    df.loc[id_count, 'Native geographical resolution scale'] = 'Global'
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

        # only selecting global average CFs
        db = db.loc[[i for i in db.index if db.loc[i, 'Resolution'] in 'Global']]

        # the dataframe where we will reconstruct the database
        df = pd.DataFrame(None, columns=self.master_db.columns)

        for i in db.index:
            id_count = len(self.master_db)
            df.loc[id_count, 'Impact category'] = db.loc[i, 'Impact category']
            df.loc[id_count, 'CF unit'] = db.loc[i, 'Unit'].split('[')[1].split(']')[0].split('/')[0]
            df.loc[id_count, 'Elem flow unit'] = db.loc[i, 'Unit'].split('[')[1].split(']')[0].split('/')[1]
            # hardcoded comp and subcomp
            df.loc[id_count, 'Compartment'] = 'Raw'
            df.loc[id_count, 'Sub-compartment'] = 'land'
            df.loc[id_count, 'Elem flow name'] = db.loc[i, 'Elem flow']
            # careful, different name
            df.loc[id_count, 'MP or Damage'] = db.loc[i, 'MP or Damage']
            df.loc[id_count, 'Native geographical resolution scale'] = 'Global'
            df.loc[id_count, 'CF value'] = db.loc[i, 'Weighted Average']
            id_count += 1
            self.master_db = pd.concat([self.master_db, df])

        self.master_db = clean_up_dataframe(self.master_db)

        # --------- Land transformation -----------

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - LandTrans - aggregated]', con=self.conn)

        # only selecting global average CFs
        db = db.loc[[i for i in db.index if db.loc[i, 'Resolution'] in 'Global']]

        # the dataframe where we will reconstruct the database
        df = pd.DataFrame(None, columns=self.master_db.columns)

        for i in db.index:
            id_count = len(self.master_db)
            df.loc[id_count, 'Impact category'] = db.loc[i, 'Impact category']
            df.loc[id_count, 'CF unit'] = db.loc[i, 'Unit'].split('[')[1].split(']')[0].split('/')[0]
            df.loc[id_count, 'Elem flow unit'] = db.loc[i, 'Unit'].split('[')[1].split(']')[0].split('/')[1]
            # hardcoded comp and subcomp
            df.loc[id_count, 'Compartment'] = 'Raw'
            df.loc[id_count, 'Sub-compartment'] = 'land'
            df.loc[id_count, 'Elem flow name'] = db.loc[i, 'Elem flow']
            # careful, different name
            df.loc[id_count, 'MP or Damage'] = db.loc[i, 'MP or Damage']
            df.loc[id_count, 'Native geographical resolution scale'] = 'Global'
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

        def add_generic_scarcity_intel(master_db, id_count):
            master_db.loc[id_count, 'Impact category'] = 'Water scarcity'
            master_db.loc[id_count, 'Elem flow unit'] = 'm3'
            master_db.loc[id_count, 'Native geographical resolution scale'] = 'Global'
            master_db.loc[id_count, 'MP or Damage'] = 'Midpoint'
            master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
            master_db.loc[id_count, 'CF unit'] = 'm3 world-eq'
            master_db.loc[id_count, 'CAS number'] = '007732-18-5'

        for flow in ['UNKNOWN', 'AGRI', 'NON-AGRI']:
            data = db.loc[
                [i for i in db.index if (db.loc[i, 'Resolution'] in 'Global' and db.loc[i, 'Elem flow'] == flow)]]
            # Water comp
            id_count = len(self.master_db)
            add_generic_scarcity_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Compartment'] = 'Water'
            self.master_db.loc[id_count, 'CF value'] = - data.loc[:, 'Weighted Average'].iloc[0]
            if flow == 'UNKNOWN':
                self.master_db.loc[id_count, 'Elem flow name'] = 'Water'
            else:
                self.master_db.loc[id_count, 'Elem flow name'] = 'Water, ' + flow.lower()
            # Raw comp
            id_count += 1
            add_generic_scarcity_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Compartment'] = 'Raw'
            self.master_db.loc[id_count, 'CF value'] = data.loc[:, 'Weighted Average'].iloc[0]
            if flow == 'UNKNOWN':
                self.master_db.loc[id_count, 'Elem flow name'] = 'Water'
            else:
                self.master_db.loc[id_count, 'Elem flow name'] = 'Water, ' + flow.lower()

        self.master_db = clean_up_dataframe(self.master_db)

        # creating the other water flows from the default water flow
        other_water = ['Water, lake','Water, river','Water, unspecified natural origin',
                       'Water, well, in ground','Water, cooling, unspecified natural origin']
        for water in other_water:
            df = self.master_db.loc[
                [i for i in self.master_db.index if (self.master_db.loc[i, 'Impact category'] == 'Water scarcity' and
                                                   self.master_db.loc[i, 'Elem flow name'] == 'Water' and
                                                   self.master_db.loc[i, 'Compartment'] == 'Raw')]]
            df.loc[:,'Elem flow name'] = water
            self.master_db = pd.concat([self.master_db,df])
            self.master_db = clean_up_dataframe(self.master_db)

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
            master_db.loc[id_count, 'Native geographical resolution scale'] = 'Global'
            master_db.loc[id_count, 'MP or Damage'] = 'Damage'
            master_db.loc[id_count, 'CF unit'] = 'PDF.m2.yr'
            master_db.loc[id_count, 'Elem flow unit'] = 'm3'
            master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'

        # Water comp
        id_count = len(self.master_db)
        add_generic_water_avail_eq_intel(self.master_db, id_count)
        self.master_db.loc[id_count, 'Elem flow name'] = 'Water'
        self.master_db.loc[id_count, 'Compartment'] = 'Water'
        self.master_db.loc[id_count, 'CF value'] = - db.loc[:, 'CF value'].median()
        # Raw comp
        id_count += 1
        add_generic_water_avail_eq_intel(self.master_db, id_count)
        self.master_db.loc[id_count, 'Elem flow name'] = 'Water'
        self.master_db.loc[id_count, 'Compartment'] = 'Raw'
        self.master_db.loc[id_count, 'CF value'] = db.loc[:, 'CF value'].median()

        self.master_db = clean_up_dataframe(self.master_db)

        # creating the other water flows from the default water flow
        other_water = ['Water, lake','Water, river','Water, unspecified natural origin',
                       'Water, well, in ground','Water, cooling, unspecified natural origin']
        for water in other_water:
            df = self.master_db.loc[
                [i for i in self.master_db.index if (self.master_db.loc[i, 'Impact category'] == 'Water availability, freshwater ecosystem' and
                                                   self.master_db.loc[i, 'Elem flow name'] == 'Water' and
                                                   self.master_db.loc[i, 'Compartment'] == 'Raw')]]
            df.loc[:,'Elem flow name'] = water
            self.master_db = pd.concat([self.master_db,df])
            self.master_db = clean_up_dataframe(self.master_db)

        self.master_db = clean_up_dataframe(self.master_db)

    def load_water_availability_hh_cfs(self):
        """
        Load CFs for water availability human health.

        Concerned impact categories:
            - Water availability, human health

        :return: update master_db
        """

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - WaterAvailab_HH - aggregated]', con=self.conn)

        def add_generic_water_avai_hh_intel(master_db, id_count):
            master_db.loc[id_count, 'Impact category'] = 'Water availability, human health'
            master_db.loc[id_count, 'Elem flow unit'] = 'm3'
            master_db.loc[id_count, 'Native geographical resolution scale'] = 'Global'
            master_db.loc[id_count, 'MP or Damage'] = 'Damage'
            master_db.loc[id_count, 'CF unit'] = 'DALY'

        # Water comp / unspecified subcomp
        id_count = len(self.master_db)
        add_generic_water_avai_hh_intel(self.master_db, id_count)
        self.master_db.loc[id_count, 'Compartment'] = 'Water'
        self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
        self.master_db.loc[id_count, 'Elem flow name'] = 'Water'
        self.master_db.loc[id_count, 'CF value'] = - db.loc[[i for i in db.index if (db.loc[i, 'Resolution'] in 'Global' and
                                                                                db.loc[i, 'Elem flow'] == 'Unknown')],
                                                       'Weighted Average'].iloc[0]
        # Water comp / lake subcomp
        id_count += 1
        add_generic_water_avai_hh_intel(self.master_db, id_count)
        self.master_db.loc[id_count, 'Compartment'] = 'Water'
        self.master_db.loc[id_count, 'Sub-compartment'] = 'lake'
        self.master_db.loc[id_count, 'Elem flow name'] = 'Water'
        self.master_db.loc[id_count, 'CF value'] = - db.loc[[i for i in db.index if (db.loc[i, 'Resolution'] in 'Global' and
                                                                                db.loc[i, 'Elem flow'] == 'Surface')],
                                                       'Weighted Average'].iloc[0]
        # Water comp / river subcomp
        id_count += 1
        add_generic_water_avai_hh_intel(self.master_db, id_count)
        self.master_db.loc[id_count, 'Compartment'] = 'Water'
        self.master_db.loc[id_count, 'Sub-compartment'] = 'river'
        self.master_db.loc[id_count, 'Elem flow name'] = 'Water'
        self.master_db.loc[id_count, 'CF value'] = - db.loc[[i for i in db.index if (db.loc[i, 'Resolution'] in 'Global' and
                                                                                db.loc[i, 'Elem flow'] == 'Surface')],
                                                       'Weighted Average'].iloc[0]
        # Water comp / groundwater subcomp
        id_count += 1
        add_generic_water_avai_hh_intel(self.master_db, id_count)
        self.master_db.loc[id_count, 'Compartment'] = 'Water'
        self.master_db.loc[id_count, 'Sub-compartment'] = 'groundwater'
        self.master_db.loc[id_count, 'Elem flow name'] = 'Water'
        self.master_db.loc[id_count, 'CF value'] = - db.loc[[i for i in db.index if (db.loc[i, 'Resolution'] in 'Global' and
                                                                                db.loc[i, 'Elem flow'] == 'Ground')],
                                                       'Weighted Average'].iloc[0]
        # Raw comp / unspecified water
        id_count += 1
        add_generic_water_avai_hh_intel(self.master_db, id_count)
        self.master_db.loc[id_count, 'Compartment'] = 'Raw'
        self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
        self.master_db.loc[id_count, 'Elem flow name'] = 'Water, unspecified natural origin'
        self.master_db.loc[id_count, 'CF value'] = db.loc[[i for i in db.index if (db.loc[i, 'Resolution'] in 'Global' and
                                                                              db.loc[i, 'Elem flow'] == 'Unknown')],
                                                     'Weighted Average'].iloc[0]
        # Raw comp / lake water
        id_count += 1
        add_generic_water_avai_hh_intel(self.master_db, id_count)
        self.master_db.loc[id_count, 'Compartment'] = 'Raw'
        self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
        self.master_db.loc[id_count, 'Elem flow name'] = 'Water, lake'
        self.master_db.loc[id_count, 'CF value'] = db.loc[[i for i in db.index if (db.loc[i, 'Resolution'] in 'Global' and
                                                                              db.loc[i, 'Elem flow'] == 'Surface')],
                                                     'Weighted Average'].iloc[0]
        # Raw comp / river water
        id_count += 1
        add_generic_water_avai_hh_intel(self.master_db, id_count)
        self.master_db.loc[id_count, 'Compartment'] = 'Raw'
        self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
        self.master_db.loc[id_count, 'Elem flow name'] = 'Water, river'
        self.master_db.loc[id_count, 'CF value'] = db.loc[[i for i in db.index if (db.loc[i, 'Resolution'] in 'Global' and
                                                                              db.loc[i, 'Elem flow'] == 'Surface')],
                                                     'Weighted Average'].iloc[0]
        # Raw comp / well water
        id_count += 1
        add_generic_water_avai_hh_intel(self.master_db, id_count)
        self.master_db.loc[id_count, 'Compartment'] = 'Raw'
        self.master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'
        self.master_db.loc[id_count, 'Elem flow name'] = 'Water, well, in ground'
        self.master_db.loc[id_count, 'CF value'] = db.loc[[i for i in db.index if (db.loc[i, 'Resolution'] in 'Global' and
                                                                              db.loc[i, 'Elem flow'] == 'Ground')],
                                                     'Weighted Average'].iloc[0]

    def load_water_availability_terr_cfs(self):
        """
        Load CFs for water availability terrestrial ecosystem.

        Concerned impact categories:
            - Water availability, terrestrial ecosystem

        :return: update master_db
        """

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - WaterAvailab_Terr - aggregated]', con=self.conn)

        def add_generic_water_avai_terr_intel(master_db, id_count):
            master_db.loc[id_count, 'Impact category'] = 'Water availability, terrestrial ecosystem'
            master_db.loc[id_count, 'Elem flow unit'] = 'm3'
            master_db.loc[id_count, 'Native geographical resolution scale'] = 'Global'
            master_db.loc[id_count, 'MP or Damage'] = 'Damage'
            master_db.loc[id_count, 'CF unit'] = 'PDF.m2.yr'
            master_db.loc[id_count, 'Sub-compartment'] = 'in water'
            master_db.loc[id_count, 'Compartment'] = 'Raw'

        data = db.loc[[i for i in db.index if db.loc[i, 'Elem flow name'] in ['Water, shallow well, in ground',
                                                                              'Water, well, in ground']]]
        for water_flow in data.index:
            id_count = len(self.master_db)
            add_generic_water_avai_terr_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Elem flow name'] = data.loc[water_flow, 'Elem flow name']
            self.master_db.loc[id_count, 'CF value'] = data.loc[water_flow, 'CF value']

        self.master_db = clean_up_dataframe(self.master_db)

    def load_thermally_polluted_water_cfs(self):
        """
        Load CFs for thermally polluted water.

        Concerned impact categories:
            - Thermally polluted water

        :return: update master_db
        """

        db = pd.read_sql(sql='SELECT * FROM [CF - regionalized - ThermallyPollutedWater - aggregated]', con=self.conn)

        def add_generic_thermally_intel(master_db, id_count):
            master_db.loc[id_count, 'Impact category'] = 'Thermally polluted water'
            master_db.loc[id_count, 'Elem flow unit'] = 'm3'
            master_db.loc[id_count, 'Native geographical resolution scale'] = 'Global'
            master_db.loc[id_count, 'MP or Damage'] = 'Damage'
            master_db.loc[id_count, 'CF unit'] = 'PDF.m2.yr'
            master_db.loc[id_count, 'Sub-compartment'] = 'in water'
            master_db.loc[id_count, 'Compartment'] = 'Raw'

        data = db.loc[
            [i for i in db.index if db.loc[i, 'Elem flow name'] in ['Water, cooling, unspecified natural origin']]]
        for water_flow in data.index:
            id_count = len(self.master_db)
            add_generic_thermally_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Elem flow name'] = data.loc[water_flow, 'Elem flow name']
            self.master_db.loc[id_count, 'CF value'] = data.loc[water_flow, 'CF value']

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
                    'Soil': ['industrial'],
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
                          'ocean': ['Ionizing radiation, ecosystem quality',
                                    'Ionizing radiation, human health',
                                    'Ionizing radiations',
                                    'Marine eutrophication']}
        to_zero = {'groundwater': ['Freshwater ecotoxicity',
                                   'Freshwater ecotoxicity, long term',
                                   'Freshwater ecotoxicity, short term',
                                   'Freshwater eutrophication',
                                   'Human toxicity cancer',
                                   'Human toxicity cancer, long term',
                                   'Human toxicity cancer, short term',
                                   'Human toxicity non cancer',
                                   'Human toxicity non cancer, long term',
                                   'Human toxicity non cancer, short term',
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
                                              'Human toxicity non cancer',
                                              'Human toxicity non cancer, long term',
                                              'Human toxicity non cancer, short term',
                                              'Ionizing radiation, ecosystem quality',
                                              'Ionizing radiation, human health',
                                              'Ionizing radiations',
                                              'Marine eutrophication'],
                   'ocean': ['Freshwater ecotoxicity',
                             'Freshwater ecotoxicity, long term',
                             'Freshwater ecotoxicity, short term',
                             'Freshwater eutrophication',
                             'Human toxicity cancer',
                             'Human toxicity cancer, long term',
                             'Human toxicity cancer, short term',
                             'Human toxicity non cancer',
                             'Human toxicity non cancer, long term',
                             'Human toxicity non cancer, short term',
                             'Water availability, freshwater ecosystem',
                             'Water availability, human health',
                             'Water scarcity']}
        for subcomp in to_unspecified:
            for cat in to_unspecified[subcomp]:
                data = water_comp.loc[[i for i in water_comp.index if (water_comp.loc[i, 'Impact category'] == cat and
                                                                       water_comp.loc[
                                                                           i, 'Sub-compartment'] == '(unspecified)')]]
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
            # slice dataframe to only keep the corresponding subcomp
            data = self.master_db.loc[[i for i in self.master_db.index if (self.master_db.loc[i, 'Sub-compartment']
                                                                           == long_term_subcomps[subcomp])]].copy()
            for cat in long_term_cats:
                # slice dataframe to only keep corresponding impact category
                df = data.loc[[i for i in data.index if (cat in data.loc[i, 'Impact category'] and
                                                         cat != data.loc[i, 'Impact category'])]]
                # remove the "short term" and "long term" from impact category name
                df.loc[:, 'Impact category'] = [','.join(i.split(',')[:-1]) for i in df.loc[:, 'Impact category']]
                # now that they have the same category name, we can add their CF value by merging dataframes
                df = df.groupby(by=['Impact category', 'CF unit', 'Compartment', 'Sub-compartment', 'Elem flow name',
                                    'CAS number', 'Elem flow unit', 'MP or Damage',
                                    'Native geographical resolution scale']).sum().reset_index()
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

        # special case climate change, short term, i.e., the only short term impact category at midpoint
        self.master_db.loc[[i for i in self.master_db.index if
                            (self.master_db.loc[i, 'Impact category'] == 'Climate change, short term' and
                             self.master_db.loc[i, 'Sub-compartment'] == 'low. pop., long-term')], 'CF value'] = 0

    def link_to_ecoinvent(self):
        """
        Function that links names of substance from IW+ to the names of ecoinvent.
        :return: self.ei35_iw, self.ei36_iw, self.ei371_iw, self.ei38_iw
        """

        versions_ei = ['3.5', '3.6', '3.7.1', '3.8']

        for version_ei in versions_ei:
            ei_iw_db = self.master_db.copy()

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
                        ['Fossil and nuclear energy use', 'MJ deprived', 'natural resource', subcomp, energy, None, 1.0,
                         'MJ', 'Midpoint', 'Global'],
                        ['Impact category', 'CF unit', 'Compartment', 'Sub-compartment', 'Elem flow name', 'CAS number',
                         'CF value',
                         'Elem flow unit', 'MP or Damage', 'Native geographical resolution scale']).T
                    ei_iw_db = pd.concat([ei_iw_db, add])
                    ei_iw_db = clean_up_dataframe(ei_iw_db)

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

    def produce_files(self):
        """
        Function producing the different IW+ files for the different versions.
        :return: the IW+ files
        """

        path = pkg_resources.resource_filename(__name__, '/Databases/Impact_world_' + self.version)

        if not os.path.exists(path + '/Excel/'):
            os.makedirs(path + '/Excel/')
        if not os.path.exists(path + '/DataFrame/'):
            os.makedirs(path + '/DataFrame/')

        self.master_db.to_excel(path + '/Excel/impact_world_plus_' + self.version + '_dev.xlsx')
        self.ei35_iw.to_excel(path + '/Excel/impact_world_plus_' + self.version + '_ecoinvent_v35.xlsx')
        self.ei36_iw.to_excel(path + '/Excel/impact_world_plus_' + self.version + '_ecoinvent_v36.xlsx')
        self.ei371_iw.to_excel(path + '/Excel/impact_world_plus_' + self.version + '_ecoinvent_v371.xlsx')
        self.ei38_iw.to_excel(path + '/Excel/impact_world_plus_' + self.version + '_ecoinvent_v38.xlsx')
        self.ei35_iw_as_matrix.to_excel(path + '/DataFrame/impact_world_plus_' + self.version + '_ecoinvent_v35.xlsx')
        self.ei36_iw_as_matrix.to_excel(path + '/DataFrame/impact_world_plus_' + self.version + '_ecoinvent_v36.xlsx')
        self.ei371_iw_as_matrix.to_excel(path + '/DataFrame/impact_world_plus_' + self.version + '_ecoinvent_v371.xlsx')
        self.ei38_iw_as_matrix.to_excel(path + '/DataFrame/impact_world_plus_' + self.version + '_ecoinvent_v38.xlsx')

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

# -------------- Support modules -------------------

def clean_up_dataframe(df):
    # remove duplicates
    df = df.drop_duplicates()
    # fix index
    df = df.reset_index().drop('index',axis=1)
    return df


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
