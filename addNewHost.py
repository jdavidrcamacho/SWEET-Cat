#!/usr/bin/env python
# -*- coding: utf8 -*-
import numpy as np
from astropy import coordinates as coord
from astropy import units as u
import pandas as pd
from clint.textui import puts, colored
import time
from ParallaxSpec import parallax
from astroquery.simbad import Simbad
import warnings
warnings.filterwarnings('ignore')
from astroquery.irsa_dust import IrsaDust
from astroquery.vizier import Vizier


def GAIAplx(ra,de):
    v = Vizier(columns=["*", "+_r"], catalog='I/345/gaia2')
    pos = coord.SkyCoord(ra=ra, dec=de, unit=(u.hourangle,u.deg), frame='icrs', 
                         obstime='J2000')
    result = v.query_region(pos, radius="10s", catalog='I/345/gaia2')
    #Moving the positions to 2000
    try:
        deltat = -15.5
        sep = []
        for ig, _ in enumerate(result[0]['Source']):
            raold=result[0]['RA_ICRS'].data[ig]+(result[0]['pmRA'].data[ig]*deltat)/3600000.
            deold=result[0]['DE_ICRS'].data[ig]+(result[0]['pmDE'].data[ig]*deltat)/3600000.
            posold = coord.ICRS(ra=raold * u.deg, dec=deold * u.deg)
            sep.append(pos.separation(posold).arcsecond)
        indG = np.argmin(sep)
        if sep[indG]<1.5 and result[0]['Plx'].data[indG]>0:
            return str(round(result[0]['Plx'].data[indG],2)), \
                    str(round(result[0]['e_Plx'].data[indG],2))
    except:
        return 'NULL','NULL'
    return 'NULL','NULL'


def torres(name, teff=False, logg=False, feh=False):
    """
    Calculates the mass and error from Torres. See source for more information
    """
    from TorresMass import massTorres
    T, Terr = teff
    L, Lerr = logg
    F, Ferr = feh
    try:
        Terr, Lerr, Ferr = float(Terr), float(Lerr), float(Ferr)
        T, L, F = float(T), float(L), float(F)
    except ValueError:
        puts(colored.red('No mass derived for this star...'))
        return 'NULL', 'NULL'
    M, Merr = massTorres(T, Terr, L, Lerr, F, Ferr)
    puts(colored.green('Done'))
    return round(M, 2), round(Merr, 2)


def variable_assignment(digits):
    """ Giving values to our stellar parameters """
    x = input(digits)
    if len(x) == 0:
        x = 'NULL'
    return x


if __name__ == '__main__':
    with open('names.txt') as f:
        stars = f.readlines()
    f.close()
    var = 'Y'
    #Read the data from exoplanet.eu
    fields = ['star_name', 'ra', 'dec', 'mag_v', 'star_metallicity', 
              'star_metallicity_error_min','star_metallicity_error_max',
              'star_teff','star_teff_error_min','star_teff_error_max']
    exo_all = pd.read_csv('exo.csv', skipinitialspace=True, usecols=fields)
    #Remove trailing whitespaces
    exo_all.star_name = exo_all.star_name.str.strip()
    output = 'WEBSITE_online_ADD.rdb'
    for i, star in enumerate(stars):
        star = star.strip('\n')
        exo = exo_all[exo_all.star_name == star]
        next = True
        print('')
        print('Star: ' + colored.green(star))
        try:
            name = exo.star_name.values[0]
        except IndexError as e:
            print('')
            puts(colored.red(star) + ' not found. Star added in the file manual.list.')
            print('')
            manual = open('manual.list', "a")
            manual.write(star+'\n')
            manual.close()
            next = False
            #Update the list of new hosts
            with open('names.txt', 'w') as names:
                #if the last star was added so no star is updated
                if i+1 == len(stars):
                    names.write('')
                else:
                    for j in stars[i+1:]:
                        names.write(j)
            names.close()
            print('-------------------------------')
        #if the star is found in the exoplanet.eu
        if next:
            print('')
            var = input('Continue? [Y/N]: ')
            if var.upper().strip()=='Y':
                #Get RA and dec
                ra, dec = float(exo.ra.values[0]), float(exo.dec.values[0])
                c = coord.SkyCoord(ra, dec, unit=(u.degree, u.degree), frame='icrs')
                RA = list(c.ra.hms)
                RA[0] = str(int(RA[0])).zfill(2)
                RA[1] = str(int(RA[1])).zfill(2)
                RA[2] = str(round(RA[2], 2)).zfill(4)
                if len(RA[2]) == 4:
                    RA[2] += '0'
                RA = "{0} {1} {2}".format(*RA)
                DEC = list(c.dec.dms)
                DEC[0] = str(int(DEC[0])).zfill(2)
                DEC[1] = str(abs(int(DEC[1]))).zfill(2)
                DEC[2] = str(abs(round(DEC[2], 2))).zfill(4)
                if int(DEC[0]) > 0:
                    DEC[0] = '+'+DEC[0]
                if len(DEC[2]) == 4:
                    DEC[2] += '0'
                DEC = "{0} {1} {2}".format(*DEC)
                #search in Simbad the parallax, Vmag and spectral type
                customSimbad = Simbad()
                customSimbad.add_votable_fields('plx', 'plx_error', 'flux(V)', 
                                                'flux_error(V)', 'sptype', 
                                                'otype', 'ids')
                result = customSimbad.query_region(coord.SkyCoord(ra = c.ra, 
                                                                  dec = c.dec,
                                                                  frame = 'icrs'),radius='15s')
                empty = 'NULL'
                #Here comes the user interface part...
                puts(colored.black('\nStandard parameters\n'))
                #The metallicity
                if ~np.isnan(exo.star_metallicity_error_min.values[0]) \
                    and ~np.isnan(exo.star_metallicity_error_max.values[0]):
                    errFeH_exo = (exo.star_metallicity_error_min.values[0] \
                                  +exo.star_metallicity_error_max.values[0])/2
                elif ~np.isnan(exo.star_metallicity_error_min.values[0]):
                    errFeH_exo = exo.star_metallicity_error_min.values[0]
                elif ~np.isnan(exo.star_metallicity_error_max.values[0]):
                    errFeH_exo = exo.star_metallicity_error_max.values[0]
                else:
                    errFeH_exo = np.nan 
                FeH_exo = exo.star_metallicity.values[0]
                if np.isnan(FeH_exo):
                    puts('The ' + colored.yellow('[Fe/H]'))
                    FeH = variable_assignment(2)
                    puts('The error on ' + colored.yellow('[Fe/H]'))
                    Ferr = variable_assignment(2)
                else:
                    FeH = round(float(FeH_exo), 2)
                    if np.isnan(errFeH_exo):
                        puts('The error on ' + colored.yellow('[Fe/H]'))
                        Ferr = variable_assignment(2)
                    else:
                        Ferr=round(errFeH_exo, 2)
                #The effective temperature
                if ~np.isnan(exo.star_teff_error_min.values[0]) \
                    and ~np.isnan(exo.star_teff_error_max.values[0]):
                    errTeff_exo = (exo.star_teff_error_min.values[0] \
                                   +exo.star_teff_error_max.values[0])/2
                elif ~np.isnan(exo.star_teff_error_min.values[0]):
                    errTeff_exo = exo.star_teff_error_min.values[0]
                elif ~np.isnan(exo.star_teff_error_max.values[0]):
                    errTeff_exo = exo.star_teff_error_max.values[0]
                else:
                    errTeff_exo = np.nan 
                Teff_exo = exo.star_teff.values[0]
                if np.isnan(Teff_exo):
                    puts('The ' + colored.yellow('Teff'))
                    Teff = variable_assignment(0)
                    puts('The error on ' + colored.yellow('Teff'))
                    Tefferr = variable_assignment(0)
                else:
                    #the Teff is not float
                    Teff = int(Teff_exo)
                    if ~np.isnan(errTeff_exo):
                        Tefferr = int(errTeff_exo)
                    else:
                        puts('The error on ' + colored.yellow('Teff'))
                        Tefferr = variable_assignment(0)
                #The log g
                puts('The ' + colored.yellow('logg'))
                logg = variable_assignment("> ")
                puts('The error on ' + colored.yellow('logg'))
                loggerr = variable_assignment("> ")
                #The mass
                puts(colored.magenta('Calculating the mass...'))
                M, Merr = torres(name, [Teff, Tefferr], [logg, loggerr], feh=[FeH, Ferr])
                #The microturbulence number
                puts('The '+colored.yellow('microturbulence'))
                vt = variable_assignment("> ")
                puts('The error on '+colored.yellow('microturbulence'))
                vterr = variable_assignment("> ")
                #Author and link to ADS
                puts('Who is the '+colored.yellow('author?'))
                author = input('> ').strip()
                if author == '':
                    author = empty                
                puts('Link to article ('+colored.yellow('ADS')+')')
                link = input('> ').strip()
                if link == '':
                    link = empty
                #Source flag
                puts(colored.yellow('Source flag'))
                source = input('(0/1) > ')
                if source == '':
                    source = '0'
                V_exo=exo.mag_v.values[0]
                try:
                    #select the star and not the planet, they have the same coordinates
                    if len(result)>1:
                        indr = np.where((result['OTYPE']!='Planet') \
                                        &(result['OTYPE']!='Planet?') \
                                        &(result['OTYPE'][1]!='brownD*'))[0][0]
                    else:
                        indr = 0
                    RA = str(result['RA'][indr])[:11]
                    DEC = str(result['DEC'][indr])[:12]
                    #The HD number
                    HD = empty
                    for iname in result['IDS'][indr].split('|'):
                        if iname[:2] == 'HD':
                            HD = iname.replace('HD ','')
                    #The V magnitude
                    if type(result['FLUX_V'][indr])!=np.ma.core.MaskedConstant:
                        V = round(float(result['FLUX_V'][indr]), 2)
                        if type(result['FLUX_ERROR_V'][indr])!=np.ma.core.MaskedConstant:
                            Verr = round(float(result['FLUX_ERROR_V'][indr]), 2)
                        else:
                            print('\nV magnitude = '+str(V))
                            puts('The error on ' + colored.yellow('V magnitude'))
                            Verr = variable_assignment("> ")
                            if Verr == '':
                                Verr = 'NULL'
                    else:
                        if ~np.isnan(V_exo):
                            V = round(float(V_exo), 2)
                        else:    
                            puts('The ' + colored.yellow('V magnitude'))
                            V = variable_assignment("> ")
                            if V == '':
                                V = 'NULL'
                        print('\nV magnitude = '+str(V))
                        puts('The error on ' + colored.yellow('V magnitude'))
                        Verr = variable_assignment("> ")
                        if Verr == '':
                            Verr = 'NULL'
                    # The parallax
                    plx,eplx=GAIAplx(RA, DEC)
                    if plx!='NULL':
                        p = plx
                        perr = eplx
                        pflag = 'GAIADR2' 
                    elif type(result['PLX_VALUE'][indr])!=np.ma.core.MaskedConstant:
                        p = round(float(result['PLX_VALUE'][indr]),2)
                        if type(result['PLX_VALUE'][indr])!=np.ma.core.MaskedConstant:
                            perr = round(float(result['PLX_ERROR'][indr]),2)
                        else:
                            perr=empty
                        pflag = 'Simbad'
                    else:
                        try:
                            pos = coord.SkyCoord(ra=ra, dec=dec, 
                                                 unit=(u.hourangle,u.deg), 
                                                 frame='icrs')
                            #AvSF = Schlafly & Finkbeiner 2011 (ApJ 737, 103)
                            tableAv = IrsaDust.get_query_table(pos, radius='02d',
                                                               section='ebv',
                                                               timeout=60)
                            Av = tableAv['ext SandF mean'].data[0]
                            Averr = tableAv['ext SandF std'].data[0]
                        except:
                            Av = 0
                            Averr = 0
                        try:    
                            p , perr = map(lambda x: round(x,2), 
                                           parallax(Teff,Tefferr, float(logg), 
                                                    float(loggerr), V, Verr, 
                                                    M, Merr, Av, Averr))
                            pflag = 'Spec'
                        except:
                            p = 'NULL'
                            perr = 'NULL'
                            pflag = 'NULL' 
                    #Comments
                    if result['SP_TYPE'][indr]!='' and result['SP_TYPE'][indr][0]=='M':
                        comment = result['SP_TYPE'][indr]
                    else:                    
                        puts('Any '+colored.yellow('comments'))
                        puts('E.g. if we have a M dwarf...')
                        comment = input('> ')
                        if comment == '':
                            comment = 'NULL'   
                except:
                    #The HD number
                    puts('The '+colored.yellow('HD number'))
                    HD = input('> ')
                    if HD == '':
                        HD = 'NULL' 
                    #The V magnitude
                    if ~np.isnan(V_exo):
                        V = round(float(V_exo), 2)
                    else:    
                        puts('The ' + colored.yellow('V magnitude'))
                        V = variable_assignment("> ")
                    print('\nV magnitude = '+str(V))
                    puts('The error on ' + colored.yellow('V magnitude'))
                    Verr = variable_assignment("> ")
                    # The parallax
                    plx,eplx=GAIAplx(RA,DEC)
                    if plx!='NULL':
                        p = plx
                        perr = eplx
                        pflag = 'GAIADR2' 
                    else:
                        try:
                            pos=coord.SkyCoord(ra=RA, dec=DEC, 
                                               unit=(u.hourangle,u.deg), 
                                               frame='icrs')
                            #AvSF = Schlafly & Finkbeiner 2011 (ApJ 737, 103)
                            tableAv = IrsaDust.get_query_table(pos, radius='02d',
                                                               section='ebv',
                                                               timeout=60)
                            Av=tableAv['ext SandF mean'].data[0]
                            Averr=tableAv['ext SandF std'].data[0]
                        except:
                            Av=0
                            Averr=0
                        try:
                            p,perr = map(lambda x: round(x,2), 
                                         parallax(Teff, Tefferr, logg, loggerr,
                                                  V, Verr, M, Merr, Av, Averr))
                            pflag = 'Spec'
                            #print p,perr
                        except:
                            p = 'NULL'
                            perr = 'NULL'
                            pflag = 'NULL' 
                    #Comments
                    puts('Any '+colored.yellow('comments'))
                    puts('E.g. if we have a M dwarf...')
                    comment = input('> ')
                    if comment == '':
                        comment = 'NULL'
                #Last update
                update = str(time.strftime("%Y-%m-%d"))
                params = [name, HD, RA, DEC, V, Verr, p, perr, pflag, Teff, 
                          Tefferr,logg, loggerr, 'NULL', 'NULL', vt, vterr, 
                          FeH, Ferr, M, Merr, author, link, source, update, 
                          comment]
                params = map(str, params)
                #New host information
                with open(output, 'a') as f:
                    f.write('\n'+'\t'.join(params) + '\tNULL')
                f.close()
                #Update the list of new hosts
                with open('names.txt', 'w') as names:
                    #if the last star was added so no star is updated
                    if i+1==len(stars):
                        names.write('')
                    else:
                        for j in stars[i+1:]:
                            names.write(j)
                names.close()
                print('')
                print('-------------------------------')
            else:
                print('Bye then (¬_¬)')
                break
