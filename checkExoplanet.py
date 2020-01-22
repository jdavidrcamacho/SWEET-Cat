#!/usr/bin/env python
# -*- coding: utf8 -*-
from urllib.request import urlopen
import pandas as pd
import numpy as np
import time
from astropy.io import votable
import warnings
from astropy import coordinates as coord
from astropy import units as u
# For fun, but still useful
from clint.textui import puts, colored
warnings.simplefilter("ignore")
from io import open
from six import indexbytes


def writeFile(fname, data):
    """ Write data to a file """
    with open(fname, 'w') as f:
        f.write(data)


class Update:
    """ Check for updates to SWEET-Cat comparing with exoplanet.eu """
    def __init__(self, controversial, download = False):
        self.controversial = controversial
        self.download = download
        self.fname = 'exo.csv'
        self.blacklist = []
        # Kapteyn's can't be added with the ' in the website
        self.readSC()
        self.downloadExoplanet()


    def downloadExoplanet(self):
        """
        Download the table from exoplanetEU and save it to a file (exo.csv).
        Return a pandas DataFrame sorted in 'update'.
        """
        if self.download:
            response = urlopen("http://exoplanet.eu/catalog/votable")
            table = response.read()
            with open("exo.xml", "wb") as f:
                f.write(table)
            self.xml2csv()
        df = pd.read_csv(self.fname)
        df = df[(df.detection_type == 'Radial Velocity') \
                | (df.detection_type == 'Primary Transit') \
                | (df.detection_type == 'Astrometry')]
        self.exoplanet = df
        names=map(lambda x: self.remove_planet(x),self.exoplanet['name'])
        self.exo_names = map(lambda x: x.strip(), names)


    def xml2csv(self):
        """ Convert the saved xml file to csv and read with pandas """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            vo = votable.parse('exo.xml', invalid='mask', pedantic=False)
            vo = vo.get_first_table().to_table(use_names_over_ids=True)
            df = vo.to_pandas()
            #Divide the data in Confirmed and not.
            df[df.planet_status == 'Confirmed'].to_csv('exo.csv', index=False)
            df[(df.planet_status == 'Unconfirmed') \
               | (df.planet_status == 'Candidate')].to_csv('exo_cont.csv', index=False)


    def remove_planet(self, name):
        """ Remove the trailing b, c, d, etc in the stellar name """
        planets = 'bcdefghijB'
        for planet in planets:
            if name.endswith(' %s' % planet):
                
                return name[:-2]
        #some exoplanets have .01 or .02 in the name 
        if name.endswith('.01') or name.endswith('.02') or name.endswith('.2'):
            return name[:-3]                
        return name


    def readSC(self):
        names_ = ['name', 'hd', 'ra', 'dec', 'V', 'Verr', 'p', 'perr',
                  'pflag', 'Teff', 'Tefferr', 'logg', 'logger',
                  'n1', 'n2', 'vt', 'vterr', 'feh', 'feherr', 'M', 'Merr',
                  'author', 'link', 'source', 'update', 'comment', 'n3']
        SC = pd.read_csv('WEBSITE_online.rdb', delimiter='\t', names=names_)
        self.sc_names = map(lambda x: x.lower().replace(' ', '').replace('-', ''), SC.name)
        self.sc_names = map(str.strip, self.sc_names)
        self.sc_names_orig = map(lambda x: x.strip(), SC.name)
        self.coordinates = SC.loc[:, ['ra', 'dec']]


    def _sccoordinates(self, idx):
        """
        Transform the coordinates to degrees
        
        Parameters
        ----------
        idx : int
            The index of the SWEET-Cat catalogue
          
        Return
        ------
        RAsc : float
             RA in degrees
        DEsc : float
             DEC in degrees
        """
        aux = self.coordinates['ra'].values[idx]
        RAsc = (float(aux[0:2])+float(aux[3:5])/60.+float(aux[6:])/3600.)*15.
        aux = self.coordinates['dec'].values[idx]
        if aux[0] == '-':
            DEsc = float(aux[0:3])-float(aux[4:6])/60.-float(aux[7:])/3600.
        else:
            DEsc = float(aux[0:3])+float(aux[4:6])/60.+float(aux[7:])/3600.
        return RAsc, DEsc


    def update(self):
        #We have this already, but without the ' in the name.
        print('\n*** Matching data base ***')
        NewStars = []
        #from exoplanet.eu
        coordExo = coord.SkyCoord(ra = self.exoplanet['ra'].values, 
                                 dec = self.exoplanet['dec'].values, 
                                 unit = (u.deg,u.deg), frame = 'icrs')
        #from sweet-cat
        coordSC = coord.SkyCoord(ra = self.coordinates['ra'].values, 
                                dec = self.coordinates['dec'].values, 
                                unit = (u.hourangle,u.deg), frame = 'icrs')
        for i, exo_name in enumerate(self.exo_names):
            new = exo_name
            tmp = new.lower().replace(' ', '').replace('-', '') 
            sep = coordExo[i].separation(coordSC).arcsecond
            ind = np.where(sep<5.)[0]
            if len(ind)==0:
                try:
                    #it didn't find by position but it finds by name
                    position = self.sc_names.index(tmp)
                except:
                    position = -1
                    #it didn't find by position and neither by name
                    if (tmp not in self.blacklist):
                        NewStars.append(new)
        NewStars = sorted(list(set(NewStars)))
        Nstars = len(NewStars)
        if Nstars:
            puts(' -> ' + colored.green(str(Nstars) + " new exoplanet available!"))
            writeFile('names.txt', '\n'.join(NewStars))
            updated=False
        else:
            puts(colored.clean('*** No new updates available ***'))
            updated=True
        #removing planets that are not in Exoplanet.eu anymore
        NewStars = []
        for i, scname in enumerate(self.sc_names_orig):
            sep = coordSC[i].separation(coordExo).arcsecond
            ind = np.where(sep<5.)[0]
            if len(ind)==0:
                try:
                    #it didn't find by position but it finds by name
                    position=self.exo_names.index(scname)
                except:
                    position = -1
                    #it didn't find by position and neither by name
                    if (tmp not in self.blacklist):
                        NewStars.append(scname)
        NewStars = sorted(list(set(NewStars)))
        Nstars = len(NewStars)
        if Nstars:
            puts(' -> ' + colored.green(str(Nstars) + " exoplanet has to be removed!\n"))
            print('\n'.join(NewStars))
        else:
            puts(colored.clean(' -> No planet to remove!'))
            if updated:        
                puts(colored.clean('*** SWEET-Cat is up to date ***'))
                puts(colored.green('    Good job '))


if __name__ == '__main__':
    with open('starnotfoundinsimbad.list', 'a') as f:
        f.write(str(time.strftime("%d-%m-%Y"))+'\n')
    new = Update(controversial=False, download=True)
    new.update()
