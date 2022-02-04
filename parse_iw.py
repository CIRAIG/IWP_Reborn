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
        self.master_db = None

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

        clean_up_dataframe(self.master_db)

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

        clean_up_dataframe(self.master_db)

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

        clean_up_dataframe(self.master_db)

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

        clean_up_dataframe(self.master_db)

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
            master_db.loc[id_count, 'CAS number'] == '007732-18-5'

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

        clean_up_dataframe(self.master_db)

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

        clean_up_dataframe(self.master_db)

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
            master_db.loc[id_count, 'Sub-compartment'] = '(unspecified)'

        data = db.loc[[i for i in db.index if db.loc[i, 'Resolution'] in 'Global']]
        for water_flow in data.index:
            # Water comp
            id_count = len(self.master_db)
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Compartment'] = 'Water'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water, ' + data.loc[water_flow, 'Elem flow'].lower()
            self.master_db.loc[id_count, 'CF value'] = - data.loc[water_flow, 'Weighted Average']
            # Raw comp
            id_count += 1
            add_generic_water_avai_hh_intel(self.master_db, id_count)
            self.master_db.loc[id_count, 'Compartment'] = 'Raw'
            self.master_db.loc[id_count, 'Elem flow name'] = 'Water, ' + data.loc[water_flow, 'Elem flow'].lower()
            self.master_db.loc[id_count, 'CF value'] = data.loc[water_flow, 'Weighted Average']

        clean_up_dataframe(self.master_db)

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

        clean_up_dataframe(self.master_db)

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

        clean_up_dataframe(self.master_db)

    def produce_files(self, path):
        """
        Function producing the different IW+ files for the different versions.
        :param path: path of the folder where to store the files
        :return: the IW+ files
        """

        self.master_db.to_excel(path+'/impact_world_plus_dev_v'+self.version+'.xlsx')

# -------------- Support modules -------------------

def clean_up_dataframe(df):
    # remove duplicates
    df = df.drop_duplicates()
    # fix index
    df = df.reset_index().drop('index',axis=1)
