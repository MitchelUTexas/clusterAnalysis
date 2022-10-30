#Finds Distance estimates for a specified cluster
#Error propigation via taylor approximation (assuming symettric gaussian errors)
    #TODO: add better cluster membership filtering with algorithms
    #TODO: add isochrone fit, also get age

#Imports
from ssl import PEM_cert_to_DER_cert
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import spatial as spat
from scipy import stats as stats_sci
import statistics as stat
import sys

#Stopping printing of warnings
pd.set_option("mode.chained_assignment", None)

#Main
def main():
#----------SET USER PARAMETERS HERE---------
    #Set z-score for pm exclusion
    z=1
    #Choose GAIA CSV input file name
    GAIAname = "M48Gaia.csv"
    SIMBADname = "M48Sim.csv"
    #Set to true if using a very nearby cluster for testing purposes
    TestData = False
#----------SET USER PARAMETERS HERE---------


    #Import Query Results, see GAIA and SIMBAD documentation for csv file variable names
    gaia = pd.read_csv(GAIAname)
    simbad = pd.read_csv(SIMBADname,sep=";")
    simbad = simbad.rename(columns=lambda x: x.strip())
    simbad = simbad.rename(columns=lambda x: "spec_type" if x=="spec. type" else x)

    #Parse Simbad data
    simbad["spec_type"] = simbad["spec_type"].str.strip()
        #Keep only main sequence stars (luminosity class 5) with spectral types
    simbad = simbad[simbad["spec_type"].str.match("^[OBAFGKM]\d+(\.\d*)?V$",na=False)]
        #remove duplicates by Ra, Dec
    simbad.drop_duplicates(subset=["coord1 (ICRS,J2000/2000)"])
        #Ignoring 'A value is trying to be set on a copy of a slice from a DataFrame.' warnings, this is what I am trying to do
        #Conver Ra,Dec to decimal degrees
    simbad["ra"] = simbad["coord1 (ICRS,J2000/2000)"].map(lambda x: float(x.split()[0])*(15) + float(x.split()[1])*(15/60) + float(x.split()[2])*(15/3600))
    simbad["Dec"] = simbad["coord1 (ICRS,J2000/2000)"].map(lambda x: np.sign(float(x.split()[3]))*(np.abs(float(x.split()[3])) + float(x.split()[4])*(1/60) + float(x.split()[5])*(1/3600)) )
        #Drop old coordinate format
    simbad = simbad.drop(columns=["coord1 (ICRS,J2000/2000)"])
        #parse simbad B,V magnitudes
    simbad['Mag B'] = simbad['Mag B'].map(lambda x: x.strip())
    simbad['Mag V'] = simbad['Mag V'].map(lambda x: x.strip())
    simbad = simbad[(simbad['Mag B'] != '~') & (simbad['Mag V'] != '~')]
    simbad['Mag B'] = simbad['Mag B'].astype(float)
    simbad['Mag V'] = simbad['Mag V'].astype(float)
    simbad = simbad.rename(columns={'Mag B':"B_app",'Mag V':'V_app'})
        #get simbad color
    simbad["(B-V)_app"] = simbad['B_app'].astype(float) - simbad["V_app"].astype(float)
            #setting .03 mag errors
    simbad['err_V_app'] = .03
    simbad['err_B_app'] = .03
    simbad['err_(B-V)_app'] = .03 * np.sqrt(2) #by error propigation

    #Filtering gaia data for nonvalues and similar apparant magnitudes
    maxMag = 10
    gaia = gaia.dropna(subset=["phot_g_mean_mag","phot_bp_mean_mag","phot_rp_mean_mag","parallax"])
    gaia = gaia[(gaia["phot_g_mean_mag"] < maxMag) & (gaia["phot_bp_mean_mag"] < maxMag) & (gaia["phot_rp_mean_mag"] < maxMag) & (gaia["parallax"] > 0)]
        #Get trig parallax distance to each star
    gaia["dist_trig_parallax"] = 1/(gaia["parallax"]/1000)
    gaia["err_dist_trig_parallax"] = np.abs(-1000/(gaia["parallax"]**2)*gaia["parallax_error"])
        #get gaia color
    gaia["bp-rp"] = gaia["phot_bp_mean_mag"] - gaia["phot_rp_mean_mag"]

    #print different databases' ra, dec
    plt.figure()
    plt.title("Ra, Dec of Database Stars")
    plt.scatter(gaia["ra"],gaia["dec"],label="gaia",color='blue')
    plt.scatter(simbad["ra"],simbad["Dec"],label='simbad',color='red')
    plt.xlabel("Ra [degrees]")
    plt.ylabel("Dec [degrees]")
    plt.legend()
    plt.show()
    
    #Plot proper motions with z std circled
    plt.figure()
    plt.title("Proper Motions")
    plt.xlabel("PM RA [mas/yr]")
    plt.ylabel("PM Dec [mas/yr]")
    plt.scatter(gaia["pmra"],gaia["pmdec"])
    ra_mean = stat.mean(gaia["pmra"])
    dec_mean = stat.mean(gaia["pmdec"])
    ra_std = stat.stdev(gaia["pmra"])
    dec_std = stat.stdev(gaia["pmdec"])
    plt.xlim((ra_mean - 3*ra_std*z),(ra_mean + 5*ra_std*z))
    plt.ylim((dec_mean - 3*dec_std*z),(dec_mean + 5*dec_std*z))
    plt.gca().add_patch(Ellipse(xy=(ra_mean,dec_mean),width=(2*ra_std*z),height=(2*dec_std*z),fill=False,color="red"))
    plt.show()

    #Join gaia and simbad data
    data = coordMatch(gaia,simbad)
    #Get cluster members by removing pm and parallax z-std outliers
    data = getMems(data,z)
    gaia_mems = getMems(gaia,z)
    if TestData:
        gaia_mems_real = gaia_mems.copy()
        gaia_mems = data.copy()
        #drop cluster members without needed gspphot data
    gaia_mems = gaia_mems.dropna(subset=['ebpminrp_gspphot','ebpminrp_gspphot_upper','ebpminrp_gspphot_lower'])
    #Get intrinsic magnitudes and color excesses from spec types
    data = getSpecParams(data)

    #get cluster metallicity
    FeH = np.nanmean( gaia_mems["fem_gspspec"] - gaia_mems["mh_gspspec"] )
    FeH_lower = np.nanmean( gaia_mems["fem_gspspec_lower"] - gaia_mems["mh_gspspec_lower"] )
    FeH_upper = np.nanmean( gaia_mems["fem_gspspec_upper"] - gaia_mems["mh_gspspec_upper"] )
        #drop intermediary values
    gaia_mems = gaia_mems.drop(columns=["fem_gspspec","fem_gspspec_lower","fem_gspspec_upper","mh_gspspec","mh_gspspec_lower","mh_gspspec_upper"])
    print("Fe/H =", figRound(FeH), ", lowerConf =", figRound(FeH_lower), ", upperConf = ", figRound(FeH_upper))

    #Plot color-color
        #gaia colors
    plt.figure()
    plt.title("Color-Color (not dereddened)")
    plt.scatter(gaia_mems["phot_g_mean_mag"] - gaia_mems["phot_bp_mean_mag"],gaia_mems["phot_bp_mean_mag"] - gaia_mems["phot_rp_mean_mag"])
    plt.xlabel("(G-BP)_app [mag]")
    plt.ylabel("(BP-RP)_app [mag]")
    plt.show()
        #drop gaia magnitudes
    data = data.drop(columns=["phot_g_mean_mag","phot_bp_mean_mag","phot_rp_mean_mag","bp-rp"])

    #Format Data For Ease of Use
    data["dMod_app"] = data["V_app"] - data["V_intr"]
    data["err_dMod_app"] = np.sqrt(data["err_V_app"]**2 + data["err_V_intr"]**2)
    data["colorExcess"] = data["(B-V)_app"] - data["(B-V)_intr"]
    data["err_colorExcess"] = np.sqrt(data["err_(B-V)_app"]**2 + data["err_(B-V)_intr"]**2)

    gaia_mems["dMod_app"] = gaia_mems["phot_g_mean_mag"] - gaia_mems["mg_gspphot"]
    gaia_mems["err_dMod_app"] = (gaia_mems['mg_gspphot_upper'] - gaia_mems['mg_gspphot_lower'])/2 #no errors given for magnitudes in gaia
    gaia_mems["colorExcess"] = gaia_mems['ebpminrp_gspphot']
    gaia_mems["err_colorExcess"] = (gaia_mems["ebpminrp_gspphot_upper"] - gaia_mems['ebpminrp_gspphot_lower'])/2
    gaia_mems["(BP-RP)_intr"] = gaia_mems['bp-rp'] - gaia_mems["colorExcess"]
    gaia_mems["(BP-RP)_intr_err"] = gaia_mems["err_colorExcess"] #no errors given for magnitudes in gaia

    #Plot Cluster Trig Parallax distances
    plt.figure()
    plt.title("Trig Parallax Distances (Cluster Only)")
    plt.xlabel("Distance (pc)")
    plt.ylabel("Freq")
    plt.hist((gaia_mems["dist_trig_parallax"]), bins=100, range=(0,1000))
    plt.show()

    #Plot cluster ra, dec
    plt.figure()
    plt.title("Ra, Dec of Cluster Stars")
    plt.scatter(gaia_mems["ra"],gaia_mems["dec"])
    plt.xlabel("Ra [degrees]")
    plt.ylabel("Dec [degrees]")
    plt.show()

    #Plot Observed CMD
    plt.figure()
    plt.gca().invert_yaxis()
    plt.title("Observed CMD")
    plt.xlabel("(BP-RP)_app [Mag]")
    plt.ylabel("G_app [Mag]")
    plt.scatter(gaia_mems["bp-rp"],gaia_mems["phot_g_mean_mag"], color="red", label = "=Cluster Stars")
    plt.legend()
    plt.show()

    #Apply Variable Extinction Method to whole cluster
    if not TestData:
        R, R_err, dMod_real,dMod_real_err = varExt(gaia_mems)
    else:
        R, R_err, dMod_real,dMod_real_err = varExt(data)
    #Print Results
    print("R = " + str(sciRound(R,R_err)))
        #Use R=3.1 if regression R values doesn't make sense
    if R<0:
        R,R_err = 3.1,.1
        print("Assuming R=3.1 +/- .1")
        #Calculate distance estimate with error
    dist_ext = 10**((dMod_real + 5) / 5)
    err_dist_ext = dist_ext * 1/5 * np.log(10) * dMod_real_err

    if not TestData:
        #Create variable Extinciton graph
        plt.figure()
        plt.title("Variable Extinction (Non-Reddened Removed from R calculation)")
        plt.xlabel("E(BP-RP) [Mag]")
        plt.ylabel("dMod_app [Mag]")
        plt.scatter(gaia_mems["colorExcess"],gaia_mems["dMod_app"], color="red", label = "Cluster Stars")
        es = np.linspace(np.min(gaia_mems["colorExcess"]),np.max(gaia_mems["colorExcess"]),10000)
        plt.plot(es,(dMod_real + es*R), color="blue", label = "Variable Extinction Trend")
        plt.legend()
        plt.show()
    else:
        #Create variable Extinciton graph
        plt.figure()
        plt.title("Variable Extinction (Non-Reddened Removed from R calculation)")
        plt.xlabel("E(B-V) [Mag]")
        plt.ylabel("dMod_app [Mag]")
        plt.scatter(data["colorExcess"],data["dMod_app"], color="red", label = "Cluster Stars")
        es = np.linspace(np.min(data["colorExcess"]),np.max(data["colorExcess"]),10000)
        plt.plot(es,(dMod_real + es*R), color="blue", label = "Variable Extinction Trend")
        plt.legend()
        plt.show()

    #correct for extinction
    data = ext_correct(data,'both',R,R_err)
    gaia_mems = ext_correct(gaia_mems,'gaia',R,R_err)

    #Plot corrected CMD
    plt.figure()
    plt.gca().invert_yaxis()
    plt.title("Corrected CMD")
    plt.xlabel("(BP-RP)_intr [Mag]")
    plt.ylabel("G_rel [Mag]")
    plt.scatter(gaia_mems["(BP-RP)_intr"],gaia_mems["G_rel"], color="red", label = "Cluster Stars")
    plt.legend()
    plt.show()

    #plot both extinctions
    plt.figure()
    plt.title("Extinction Comparison")
    plt.xlabel("A_V variable extinction [Mag]")
    plt.ylabel("A_V gspphot [Mag]")
    plt.scatter(R*gaia_mems["colorExcess"],gaia_mems['ag_gspphot'])
    lp = np.linspace(min(R*gaia_mems['colorExcess']),max(R*gaia_mems['colorExcess']),1000)
    plt.plot(lp,lp)
    plt.show()

    #Get individual spectroscopic parallax distance
    data["dist_spec_parallax"] = 10**((data['V_rel'] - data["V_intr"] + 5)/5)
    data["err_dist_spec_parallax"] = .25*data['dist_spec_parallax'] #per Tsvetkov et al. paper, propigated to average
    #Plot Cluster Spec Parallax distances
    plt.figure()
    plt.title("Spec Parallax Distances (Cluster Only)")
    plt.xlabel("Distance (pc)")
    plt.ylabel("Freq")
    plt.hist(data["dist_spec_parallax"], bins=100, range=(0,max(data['dist_spec_parallax'])))
    plt.show()

    #Get Distance Estimates and print them
        #cephid distnce
    cep = gaia.dropna(subset=['pf','ag_gspphot'])
    cep['absMag_cepheid']= -1*(2.76*(np.log10(cep["pf"]))-1.0)-4.16
    cep['dist_cep'] = 10**((cep["int_average_g"] - cep['ag_gspphot'] - cep['absMag_cepheid'] +5)/5)
    cep['absMag_cepheid_err'] = 2.76*(np.log10(cep["pf_error"]))
    cep['dist_cep_err'] = cep['dist_cep'] * np.log(10) * (np.sqrt(cep["int_average_g_error"]**2 + cep['absMag_cepheid_err']**2)/5)
    dist_cep = np.mean(cep['dist_cep'])
    dist_cep_err = np.sqrt(np.nansum((cep["dist_cep_err"]**2))) / cep.shape[0]
        #Trig parallax
    dist_trigParallax = np.nanmean(gaia_mems["dist_trig_parallax"])
    err_dist_trigParallax = np.sqrt(np.nansum((gaia_mems["err_dist_trig_parallax"]**2))) / gaia_mems.shape[0]
        #spec parallax
    dist_specParallax = stat.mean(data["dist_spec_parallax"])
    err_dist_specParallax = np.sqrt(np.sum((data["err_dist_spec_parallax"]**2))) / data.shape[0]
        #MS Fit
    dMod_MS, dMod_MS_err = MS_fit(gaia_mems)
    dist_MS = 10**((dMod_MS + 5)/5)
    dist_MS_err = dist_MS * 1/5 * np.log(10) * dMod_MS_err
        #print
    print("Trig Parallax: " + str(sciRound(dist_trigParallax,err_dist_trigParallax)) + " pc")
    print("Extinction: " + str(sciRound(dist_ext,err_dist_ext)) + " pc")
    print("Spec Parallax " + str(sciRound(dist_specParallax,err_dist_specParallax)) + " pc")
    if not TestData:
        print("MS: " + str(sciRound(dist_MS,dist_MS_err)) + ' pc')
    else:
        a,b,c,d = getZAMS(gaia_mems)
        print('MS Line: ' + str(figRound(a)) + 'x^3 + ' + str(figRound(b)) + 'x^2 + ' + str(figRound(c)) + 'x + ' + str(figRound(a)))
        #Plot corrected CMD
        plt.figure()
        plt.gca().invert_yaxis()
        plt.title("Corrected CMD")
        plt.xlabel("(BP-RP)_intr [Mag]")
        plt.ylabel("G_rel [Mag]")
        plt.scatter(gaia_mems["(BP-RP)_intr"],gaia_mems["G_rel"], color="red", label = "Cluster Stars")
        cs = np.linspace(np.min(gaia_mems["(BP-RP)_intr"]),np.max(gaia_mems["(BP-RP)_intr"]),1000)
        plt.plot(cs, a*cs**3 + b*cs**2 + c*cs + d, color='blue', label='MS Line')
        plt.legend()
        plt.show()
    if cep.shape[0] > 0:
        print('Cepheid PL: ' + str(sciRound(dist_cep,dist_cep_err)) + ' pc, ' + cep.shape[0] + ' cepheids')
    else:
        print('No cepheids found in cluster')
    #virial thm
    if not TestData:
        mass_gaia, radius_gaia = virial(gaia_mems)
        print('Mass =', figRound(mass_gaia),'solar masses')
        print('Radius =', figRound(radius_gaia), 'pc')
    else:
        mass_gaia, radius_gaia = virial(gaia_mems_real)
        print('Mass =', figRound(mass_gaia),'solar masses')
        print('Radius =', figRound(radius_gaia), 'pc')

#Returns pandas dataframe by removing z-std outliers (for proper motion and parallax)
def getMems(stars,z):
    stars = stars[(np.abs(stats_sci.zscore(stars["pmra"]))<z) & (np.abs(stats_sci.zscore(stars["pmdec"]))<z) & (np.abs(stats_sci.zscore(stars["dist_trig_parallax"]))<z)]
    return stars

#Applies Variable Extinction Method and returns R value and distance estimate
def varExt(stars): 
         #perform fit
    popt,pcov = curve_fit(ext_fit, stars["colorExcess"], stars["dMod_app"],sigma=np.sqrt(stars["err_dMod_app"]**2 + (3.1)**2 * stars["err_colorExcess"]**2))
        #get fit parameters
    perr = np.sqrt(np.diag(pcov))
    r = popt[0]
    r_err = perr[0]
    dMod_real = popt[1]
    dMod_real_err = perr[1]
        #Return R and corrected distance modulus for whole cluster
    return r, r_err, dMod_real, dMod_real_err


#Returns pandas dataframe to correct for extinction and reddening
def ext_correct(stars,type,R,R_err):
    #Use intrinsic color from now on to account for reddening
    if type=='both':
        stars["V_rel"] = stars["V_app"] - R*stars["colorExcess"]
        stars["err_V_rel"] = np.sqrt(stars["err_V_app"]**2 + R**2 * stars["err_colorExcess"]**2 + stars['colorExcess']**2 *R_err**2)
    else:
        stars["G_rel"] = stars["phot_g_mean_mag"] - R*stars["colorExcess"]
            #gaia gives no errors for magnitudes
        stars["err_G_rel"] = np.sqrt(stars["err_colorExcess"]**2 + stars['colorExcess']**2 *R_err**2)
    return stars
    
#Equation for variable extinction regression
def ext_fit(e, r, dMod_real):
    dMod_app = dMod_real + r*e
    return dMod_app

#round num to 2 sig figs
def figRound(num):
    return float("%.2g" % num)

#get number decimals
def numDec(num):
    if num == 0:
        return 2
    else:
        return -int(np.floor(np.log10(abs(num)))) + 1

#do value and error rounding (to 2 sig figs in error)
def sciRound(num,err):
    err = figRound(err)
    num = round(num, numDec(err))
    return [num,err]

#match gaia and simbad data by RA, Dec
def coordMatch(gaia, simbad):
        #merge_asof needs sorted dataframes
    gaia = gaia.sort_values(by="ra")
    simbad = simbad.sort_values(by="ra")
        #get join tolerance from average ra,dec error in gaia
        #sqrt 2 factor approximates additon of error from simbad
    tolRa = np.abs( np.nanmean(gaia["ra_error"])) / 2 * np.sqrt(2)
    tolDec = np.abs( np.nanmean(gaia["dec_error"])) / 2 * np.sqrt(2)
        #merge_asof only allows merge on one columns, so drop collumns that are above error for dec
    matches = pd.merge_asof(gaia,simbad,on='ra',tolerance=tolRa,direction='nearest')
    matches = matches[np.abs(matches["dec"] - matches["Dec"]) <= tolDec]
        #renaming merged columns
    matches = matches.rename(columns={'ra_x':'ra'})
    if matches.shape[0] == 0:
        print("ERROR: No Coordinate Matches")
        sys.exit()
    return matches

#Get intrinsic color and magnitude from spectral type
def getSpecParams(frame):
            #create list of ordered spec types
    order = "OBAFGKM_"
    order2 = [("O" + str(i)) for i in range(5,10)]
    order2.extend([(let + str(i)) for let in "BAFGK" for i in range(0,10)])
    order2.extend([("M" + str(i)) for i in range(0,9)])
    #Used tabuled M_V values from Tsvetkov et al. (https://link.springer.com/content/pdf/10.1134/S1063773708010039.pdf)
        #Johnson system
    tabV = {"O5":-5.6,"O9":-4.5,"B0":-4.0,"B5":-4.2,"A0":0.6,"A5":1.9,"F0":2.7,"F5":3.5,"G0":4.4,"G5":5.1,"K0":5.9,"K5":7.3,"M0":8.8,"M5":12.3}
    #Used tabuled M_V values from:
        #Loktin and Beshenov for O4-B8 (http://ezproxy.lib.utexas.edu/login?url=https://search.ebscohost.com/login.aspx?direct=true&db=a9h&AN=7289091&site=ehost-live)
            #B9 = avg(A0,B8)
        #Grenier et al. for A0-F5 (https://adsabs.harvard.edu/pdf/1985A%26A...145..331G)
        #Mikami and Heck for F6-K6 (https://adsabs.harvard.edu/pdf/1982PASJ...34..529M)
            #K2-4 adjusted by error for high and low in range (i.e., K2 = K0-3 - error)
            #K5 = K6 - error; K7 = K6 + error; K8 = K6 +2*error; K9 = avg(M0,K8)
        #Mikami for M0-M4 (https://ui.adsabs.harvard.edu/abs/1978PASJ...30..191M/abstract)
            #adjusted by error for high and low in range (i.e., M2 = M0-1 - error)
        #Tsvetkov et al. for M5 and O5 (https://link.springer.com/content/pdf/10.1134/S1063773708010039.pdf)
            #M6-8 via linear continutation of M4 to M5 difference
    tabV = {"O5":-5.6,'O6':-5.38,'O7':-4.8,'O8':-4.66,'O9':-4.18,
        'B0':-3.55,'B1':-2.84,'B2':-2.11,'B3':-1.56,'B4':-1.38,'B5':-1.2,'B6':-.86,'B7':-.58,'B8':-.12,'B9':-.4,
        'A0':.4,'A1':.5,'A2':.6,'A3':.8,'A4':1,'A5':1.4,'A6':1.6,'A7':1.7,'A8':1.8,'A9':1.9,
        'F0':2.3,'F1':2.5,'F2':2.9,'F3':3.1,'F4':3.2,'F5':3.4, 'F6':3.43,'F7':3.43,'F8':3.53,'F9':3.53,
        'G0':3.83,'G1':3.83,'G2':4.52,'G3':4.52,'G4':4.64,'G5':4.64,'G6':4.11,'G7':4.11,'G8':5.03,'G9':5.03,
        'K0':5.65,'K1':5.65,'K2':6.05,'K3':6.44,'K4':6.83,'K5':6.9,'K6':7.22,'K7':7.54,'K8':7.86,'K9':8.08,
        'M0':8.3,'M1':8.9,'M2':8.1,'M3':8.7,'M4':9.3,'M5':12.3,'M6':15.3,"M7":18.3,"M8":21.3
        }
    #Used tabulated B-V values from FitzGerald (https://ui.adsabs.harvard.edu/abs/1970A%26A.....4..234F/abstract)
        #missing B-V tabulated values were taken as an average of the surrounding bins (A6,K6,K8)
        #Johnson system        
    tabB_V={"O5":-0.32,'O6':-.32,'O7':-.32,'O8':-.31,'O9':-.31,
        'B0':-.3,'B1':-.26,'B2':-.24,'B3':-.2,'B4':-.18,'B5':-.16,'B6':-.14,'B7':-.13,'B8':-.11,'B9':-.07,
        'A0':-.01,'A1':.02,'A2':.05,'A3':.08,'A4':.12,'A5':.15,'A6':.175,'A7':.2,'A8':.27,"A9":.3,
        'F0':.32,'F1':.34,"F2":.35,'F3':.41,'F4':.42,'F5':.45,'F6':.48,'F7':.5,'F8':.53,'F9':.56,
        'G0':.6,'G1':.62,'G2':.63,'G3':.65,'G4':.66,'G5':.68,'G6':.72,'G7':.73,'G8':.74,'G9':.76,
        'K0':.81,'K1':.86,'K2':.92,'K3':.95,'K4':1,'K5':1.15,'K6':1.24,'K7':1.33,'K8':1.35,'K9':1.37,
        'M0':1.37,'M1':1.47,'M2':1.47,'M3':1.47,'M4':1.52,"M5":1.61,'M6':1.64,'M7':1.68,'M8':1.77
        }

        #spec types round to nearest tabulated type
    frame["spec_type"] = frame["spec_type"].map(lambda x: x[:1] + str(round(float(x[1:-1]))) if (round(float(x[1:-1])) < 10) else order[order.index(x[:1])+1] + "0")
        #deal with last listed types (drop if after M8)
    frame["spec_type"] = frame["spec_type"].map(lambda x: "drop" if (x=="_0" or (float(x[1:])>8 and x[1:]=='M')) else x)
    frame = frame[frame['spec_type'] != 'drop']
        #get intrinsic/tabulated values
            #also deal with edge of arrays: I use np.sign to check if index=0 or if index>0 because lambda doesn't support elif
    frame["V_intr"] = frame["spec_type"].map(lambda x: tabV[x])
    frame["err_V_intr"] = frame["spec_type"].map(lambda x: .5) #per Tsvetkov et al. paper 
    frame["(B-V)_intr"] = frame["spec_type"].map(lambda x: tabB_V[x])
    frame["err_(B-V)_intr"] = frame["spec_type"].map(lambda x: .03) #approx difference between Tsvetkov et al. and FitxGerald colors

    return frame

#get radius and virial mass
def virial(df):
    m_pc = 3.086*10**(16) #meters per pc
    G = 6.6743*10**(-11) #gravitational constant
    m_sun = 1.98847*10**(30) #solar mass
        #drop NaNs
    df = df.dropna(subset=['radial_velocity','ra','dec','dist_trig_parallax'])
        #convert to m/s
    df['radial_velocity'] = df['radial_velocity'] * 10**3
        #get center of cluster   
    ra = np.mean(df['ra']) * np.pi/180
    dec = np.mean(df['dec']) * np.pi/180
    dist = np.mean(df['dist_trig_parallax']) * m_pc
            #--> cartesian
    x = dist*np.cos(dec)*np.cos(ra)
    y = dist*np.cos(dec)*np.sin(ra)
    z = dist * np.sin(dec)
        #get cartisian coords of each star
    df['x'] = df['dist_trig_parallax'] * m_pc * np.cos(df['dec']*np.pi/180) * np.cos(df['ra']*np.pi/180)
    df['y'] = df['dist_trig_parallax'] * m_pc * np.cos(df['dec']*np.pi/180) * np.sin(df['ra']*np.pi/180)
    df['z'] = df['dist_trig_parallax'] * m_pc * np.sin(df['dec']*np.pi/180)
    pts = list(zip(df.x,df.y,df.z))
        #get radial dist to each star
    df['radial'] = np.sqrt((df['x']-x)**2 + (df['y']-y)**2 + (df['z']-z)**2)
        #estimate radius of cluster with convex hull
    vol = spat.ConvexHull(pts).volume
    radius = (3*vol/(4*np.pi))**(1/3)
        #get velocity dispersion at radius
    temp = df[df['radial'] <= radius]
    disp = stat.stdev(temp['radial_velocity'])
        #estimate cluster mass
    m = radius*disp**2/G * (3/2)
    #convert radius to pc, mass to solar masses
    radius /= m_pc
    m /= m_sun
    #return mass and radius
    return m, radius

#MS fitting
def MS_fit(df):
    popt, pcov = curve_fit(zams,df["(BP-RP)_intr"],df['G_rel'],sigma=df['err_G_rel'])
    dMod = popt[0]
    dMod_err = np.sqrt(np.diag(pcov)[0])
    return dMod, dMod_err

#zams line (+distance modulus) for MS fit
def zams(color,dMod):
    absMag = 9.8*color**3 -28*color**2 + 30*color + 9.8 #getZAMS() result for hyades
    return absMag + dMod

def getZAMS(df):
    a,b,c,d = np.polyfit(df['(BP-RP)_intr'],df['G_rel'],3)
    return a,b,c,d

#Run Main
main()
