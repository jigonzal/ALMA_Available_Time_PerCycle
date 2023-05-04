import numpy as np
import matplotlib.pyplot as plt
import zipfile
import glob
import shutil 
import xml.etree.ElementTree as ET
import argparse
import astropy.units as u

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io import ascii
from astropy.modeling import models

from astroplan import Observer, FixedTarget,AltitudeConstraint
from astroplan.utils import time_grid_from_range
from astroplan.plots import plot_sky,plot_altitude
from matplotlib.backends.backend_pdf import PdfPages

cc = plt.rcParams['axes.prop_cycle'].by_key()['color']

def GetParametersSB(SB):
    Output = {'SB':SB}
    tree = ET.parse(SB)
    root = tree.getroot()
    for child in root:
        if child.tag=='{Alma/ObsPrep/ObsProject}ObsUnitControl':
            
            Output['arrayRequested'] = child.attrib['arrayRequested']
            if child.attrib['arrayRequested'] != 'TWELVE-M':
                print('Worng Array:',child.attrib['arrayRequested'])
                exit()
            for subchild in child:
                if subchild.tag=='{Alma/ObsPrep/ObsProject}estimatedExecutionTime':
                    ExecutionTimeXML = float(subchild.text)
                    if subchild.attrib['unit']=='min':
                        ExecutionTimeXML = ExecutionTimeXML/60.0
                    elif subchild.attrib['unit']=='s':
                        ExecutionTimeXML = ExecutionTimeXML/3600.0
                    Output['ExecutionTime']= ExecutionTimeXML

        if child.tag=='{Alma/ObsPrep/SchedBlock}SchedulingConstraints':
            aux = child.attrib['representativeReceiverBand']
            aux = aux[-2:]
            Output['Band']=aux
            for subchild in child:
                if subchild.tag=='{Alma/ObsPrep/SchedBlock}nominalConfiguration':
                    aux = subchild.text
                    aux = aux[0]+aux[3:5]
                    Output['NominalConfiguration']=aux
                if subchild.tag=='{Alma/ObsPrep/SchedBlock}representativeCoordinates':
                    for j in subchild:
                        if j.tag=='{Alma/ValueTypes}longitude':
                            # print(j.text,j.attrib['unit'])
                            Output['ra']=j.text
                        if j.tag=='{Alma/ValueTypes}latitude':
                            Output['dec']=j.text

    return Output

def LoadTimeAvailable(path):
    ALMA_Times = {}
    ConfigurationsFiles = glob.glob(path+'/C*.txt')
    ConfigurationsFiles.sort()
    for cf in ConfigurationsFiles:
        tt = ascii.read(cf)
        # print(tt)
        tt['Prop_B1to7'] = np.zeros_like(tt['col1'])
        tt['Prop_B8'] = np.zeros_like(tt['col1'])
        tt['Prop_B9'] = np.zeros_like(tt['col1'])
        tt['Prop_B10'] = np.zeros_like(tt['col1'])
        ALMA_Times[cf[-7:-4]] = tt
    return ALMA_Times


def UpdateALMA_Times(ALMA_Times,NomilalConfiguration,TimeDistribution,Band):
    Band = float(Band)
    if Band<=7:
        ALMA_Times[NomilalConfiguration]['Prop_B1to7'] = ALMA_Times[NomilalConfiguration]['Prop_B1to7'] + TimeDistribution
    elif Band==8:
        ALMA_Times[NomilalConfiguration]['Prop_B8'] = ALMA_Times[NomilalConfiguration]['Prop_B8'] + TimeDistribution
    elif Band==9:
        ALMA_Times[NomilalConfiguration]['Prop_B9'] = ALMA_Times[NomilalConfiguration]['Prop_B9'] + TimeDistribution
    elif Band==10:
        ALMA_Times[NomilalConfiguration]['Prop_B10'] = ALMA_Times[NomilalConfiguration]['Prop_B10'] + TimeDistribution                

    return ALMA_Times


def TimeFromTable(TablePath,MinimumElevation,SigmaDistributionTime):
    ALMA_Times = LoadTimeAvailable('Alma_cycle10_tables')

    try:
        table = ascii.read(TablePath)
    except:
        print('Problems reading table.. exiting')
        exit()

    FirstTarget = True
    for i in range(len(table['Name'])):
        Output = {'ExecutionTime':table['TotalTime'][i],
                  'Band':table['Band'][i],
                  'NominalConfiguration':table['Configuration'][i],
                  'ra':table['RA'][i],
                  'dec':table['DEC'][i]}

        Alma = Observer.at_site("Alma")
        constraints = [AltitudeConstraint(MinimumElevation*u.deg, 90*u.deg)]
        start_time = Time('2024-07-01 23:48:00')
        end_time = Time('2024-07-02 23:48:00')
        time_resolution = 1 * u.hour
        time_grid = time_grid_from_range([start_time, end_time],time_resolution=time_resolution)

        LST = Alma.local_sidereal_time(time_grid)
        lst = np.array([round(lst.hour,0) for lst in LST])

        Target = FixedTarget(coord=SkyCoord(ra=float(Output['ra'])*u.deg, dec=float(Output['dec'])*u.deg), name='Target')

        observability_grid = np.zeros((len(constraints), len(time_grid)))

        for i, constraint in enumerate(constraints):
            observability_grid[i, :] = constraint(Alma, Target, times=time_grid)
        observability_grid = observability_grid[0]
        ExecutionTime = Output['ExecutionTime']

        
        g1 = models.Gaussian1D(amplitude=1 / (SigmaDistributionTime * np.sqrt(2 * np.pi)), mean=Target.ra.hour, stddev=SigmaDistributionTime)
        if Target.ra.hour<=8.0:
            lst2 = lst*1
            lst2[lst2>16] = lst2[lst2>16]-24
        elif Target.ra.hour>=16.0:
            lst2 = lst*1
            lst2[lst2<8] = lst2[lst2<8]+24            
        else:
            lst2 = lst*1
        TimeDistribution = g1(lst2)
        TimeDistribution[observability_grid==0] = 0
        TimeDistribution = TimeDistribution/np.sum(TimeDistribution)

        TimeDistribution = TimeDistribution*ExecutionTime
        TimeDistribution = TimeDistribution[np.argsort(lst)]
        observability_grid = observability_grid[np.argsort(lst)]
        lst = lst[np.argsort(lst)]

        if FirstTarget:
            plt.figure()
            plot_altitude(Target, Alma, time_grid)
            plt.tight_layout()
            plt.savefig('Elevation_First_Target_Table.pdf')
            plt.close('all')

            plt.figure()
            plot_sky(Target, Alma, time_grid)
            plt.tight_layout()
            plt.savefig('Sky_First_Target_Table.pdf')   
            plt.close('all')

            plt.figure()
            ax1 = plt.gca()
            ax1.set_ylabel('observability')
            ax1.plot(lst,TimeDistribution,label='Time Distribution',ds='steps-mid',color=cc[0])
            ax1.fill_between(lst, observability_grid*np.max(TimeDistribution)*1.1, y2=0,label='Above Elevation Limit of '+str(MinimumElevation)+' deg',color=cc[1],alpha=0.5)

            ax1.set_ylabel('Time Distribution [Hours]')
            ax1.set_xlabel('LST [Hours]')
            ax1.legend(loc=0)
            plt.savefig('TimeDistribution_Obs_First_Target_Table.pdf')
            FirstTarget = False



        ALMA_Times = UpdateALMA_Times(ALMA_Times,Output['NominalConfiguration'],TimeDistribution,Output['Band'])

    TotalProposed = 0

    with PdfPages('TimePressurePerConfiguration_Table.pdf') as pdf:
        for key in ALMA_Times.keys():
            table = ALMA_Times[key]

            plt.figure()
            plt.title(key)
            plt.bar(table['col1'],table['col2'],label='Total',width=0.9)
            plt.bar(table['col1'],table['col3'],label='LP',width=0.9)
            plt.bar(table['col1'],table['col4'],label='Band8',width=0.9)
            plt.bar(table['col1'],table['col5'],label='Band9',width=0.9)
            plt.bar(table['col1'],table['col6'],label='Band10',width=0.9)
            plt.bar(table['col1'],table['Prop_B1to7'],label='Prop_B1to7',width=0.6,edgecolor='black')
            plt.bar(table['col1'],table['Prop_B8'],label='Prop_B8',width=0.5,edgecolor='black')
            plt.bar(table['col1'],table['Prop_B9'],label='Prop_B9',width=0.4,edgecolor='black')
            plt.bar(table['col1'],table['Prop_B10'],label='Prop_B10',width=0.3,edgecolor='black')
            TotalProposed = TotalProposed + np.sum(table['Prop_B1to7']) + np.sum(table['Prop_B8']) + np.sum(table['Prop_B9']) + np.sum(table['Prop_B10']) 
            plt.legend(loc=0,ncol=3)
            plt.xlabel('LST [Hours]')
            plt.ylabel('Hours')
            pdf.savefig()

            table.round(1)
            table.rename_column('col1', 'LST') 
            table.rename_column('col2', 'Total') 
            table.rename_column('col3', 'LP') 
            table.rename_column('col4', 'Band8') 
            table.rename_column('col5', 'Band9') 
            table.rename_column('col6', 'Band10') 
            table.write(key+'_Table.html', format='html',overwrite=True)
        print('Total proposed time:',round(TotalProposed,1),'hours')
    return

def TimeFromAOT(AOTPath,MinimumElevation,SigmaDistributionTime):
    ALMA_Times = LoadTimeAvailable('Alma_cycle10_tables')

    shutil.rmtree('tmp')

    try:
        print('Processing AOT file:',AOTPath)
        with zipfile.ZipFile(AOTPath, 'r') as zip_ref:
            zip_ref.extractall('tmp')
    except:
        print('Problems extracting file:',AOTPath)

    SB = glob.glob('tmp/SchedBlock*.xml')
    if len(SB)==0:
        print('No SchedBlock in AOT file. Try submitting the aot file.. ')
        return
    

    FirstTarget = True
    for sb in SB:
        Output = GetParametersSB(sb)

        Alma = Observer.at_site("Alma")
        constraints = [AltitudeConstraint(MinimumElevation*u.deg, 90*u.deg)]
        start_time = Time('2024-07-01 23:48:00')
        end_time = Time('2024-07-02 23:48:00')
        time_resolution = 1 * u.hour
        time_grid = time_grid_from_range([start_time, end_time],time_resolution=time_resolution)

        LST = Alma.local_sidereal_time(time_grid)
        lst = np.array([round(lst.hour,0) for lst in LST])

        Target = FixedTarget(coord=SkyCoord(ra=float(Output['ra'])*u.deg, dec=float(Output['dec'])*u.deg), name='Target')

        observability_grid = np.zeros((len(constraints), len(time_grid)))

        for i, constraint in enumerate(constraints):
            observability_grid[i, :] = constraint(Alma, Target, times=time_grid)
        observability_grid = observability_grid[0]
        ExecutionTime = Output['ExecutionTime']

        
        g1 = models.Gaussian1D(amplitude=1 / (SigmaDistributionTime * np.sqrt(2 * np.pi)), mean=Target.ra.hour, stddev=SigmaDistributionTime)
        if Target.ra.hour<=8.0:
            lst2 = lst*1
            lst2[lst2>16] = lst2[lst2>16]-24
        elif Target.ra.hour>=16.0:
            lst2 = lst*1
            lst2[lst2<8] = lst2[lst2<8]+24            
        else:
            lst2 = lst*1
        TimeDistribution = g1(lst2)
        TimeDistribution[observability_grid==0] = 0
        TimeDistribution = TimeDistribution/np.sum(TimeDistribution)

        TimeDistribution = TimeDistribution*ExecutionTime
        TimeDistribution = TimeDistribution[np.argsort(lst)]
        observability_grid = observability_grid[np.argsort(lst)]
        lst = lst[np.argsort(lst)]

        if FirstTarget:
            plt.figure()
            plot_altitude(Target, Alma, time_grid)
            plt.tight_layout()
            plt.savefig('Elevation_First_Target_AOT.pdf')
            plt.close('all')

            plt.figure()
            plot_sky(Target, Alma, time_grid)
            plt.tight_layout()
            plt.savefig('Sky_First_Target_AOT.pdf')   
            plt.close('all')

            plt.figure()
            ax1 = plt.gca()
            ax1.set_ylabel('Observability')
            ax1.plot(lst,TimeDistribution,label='Time Distribution',ds='steps-mid',color=cc[0])
            ax1.fill_between(lst, observability_grid*np.max(TimeDistribution)*1.1, y2=0,label='Above Elevation Limit of '+str(MinimumElevation)+' deg',color=cc[1],alpha=0.5)

            ax1.set_ylabel('Time Distribution [Hours]')
            ax1.set_xlabel('LST [Hours]')
            ax1.legend(loc=0)
            plt.savefig('TimeDistribution_Obs_First_Target_AOT.pdf')
            FirstTarget = False



        ALMA_Times = UpdateALMA_Times(ALMA_Times,Output['NominalConfiguration'],TimeDistribution,Output['Band'])

    TotalProposed = 0

    with PdfPages('TimePressurePerConfiguration_AOT.pdf') as pdf:
        for key in ALMA_Times.keys():
            table = ALMA_Times[key]

            plt.figure()
            plt.title(key)
            plt.bar(table['col1'],table['col2'],label='Total',width=0.9)
            plt.bar(table['col1'],table['col3'],label='LP',width=0.9)
            plt.bar(table['col1'],table['col4'],label='Band8',width=0.9)
            plt.bar(table['col1'],table['col5'],label='Band9',width=0.9)
            plt.bar(table['col1'],table['col6'],label='Band10',width=0.9)
            plt.bar(table['col1'],table['Prop_B1to7'],label='Prop_B1to7',width=0.6,edgecolor='black')
            plt.bar(table['col1'],table['Prop_B8'],label='Prop_B8',width=0.5,edgecolor='black')
            plt.bar(table['col1'],table['Prop_B9'],label='Prop_B9',width=0.4,edgecolor='black')
            plt.bar(table['col1'],table['Prop_B10'],label='Prop_B10',width=0.3,edgecolor='black')
            TotalProposed = TotalProposed + np.sum(table['Prop_B1to7']) + np.sum(table['Prop_B8']) + np.sum(table['Prop_B9']) + np.sum(table['Prop_B10']) 
            plt.legend(loc=0,ncol=3)
            plt.xlabel('LST [Hours]')
            plt.ylabel('Hours')
            pdf.savefig()

            table.round(1)
            table.rename_column('col1', 'LST') 
            table.rename_column('col2', 'Total') 
            table.rename_column('col3', 'LP') 
            table.rename_column('col4', 'Band8') 
            table.rename_column('col5', 'Band9') 
            table.rename_column('col6', 'Band10') 
            table.write(key+'_AOT.html', format='html',overwrite=True)
            
        print('Total proposed time:',round(TotalProposed,1),'hours')
    return
   
def main():

	#Parse the input arguments
    parser = argparse.ArgumentParser(description="Python script that estimate the time used per configuration and compares it with the available time")
    parser.add_argument('-AOTFile', type=str,default='default.aot', required=False,help = 'Path to the .aot file to explore. [Default:default.aot]')
    parser.add_argument('-TableFile', type=str, default='table.txt', required=False , help = 'Path to the ascii table to get the time used per configuration and band in the proposal. [Default:table.txt]')
    parser.add_argument('-MinimumElevation', type=float, default = 30.0, required=False,help = 'Minimum elevation for the observations to be done. OT default is 20 degree [Default:30.0]')
    parser.add_argument('-SigmaDistributionTime', type=float, default = 2.0, required=False,help = 'Sigma value for the Gaussian distribution modelling the window of osbervability, in hours.   [Default:2.0]')
    parser.add_argument('-IgnoreWarning', type=str, default = 'True',choices=['True','False'], required=False,help = 'Option to limit the verbose of waning messages, specially in astropy coords system. ')
    args = parser.parse_args()

    #Checking input arguments

    if args.IgnoreWarning == 'True':
        import warnings
        warnings.filterwarnings("ignore")

    if args.AOTFile!='default.aot':
        print(50*'*')
        print(50*'*')
        print('Using the AOT file:'+args.AOTFile)
        TimeFromAOT(args.AOTFile,args.MinimumElevation,args.SigmaDistributionTime)
        print(50*'*')
        print(50*'*')

    if args.TableFile!='table.txt':
        print(50*'*')
        print(50*'*')        
        print('Using the ascii table:'+args.TableFile)    
        TimeFromTable(args.TableFile,args.MinimumElevation,args.SigmaDistributionTime)
        print(50*'*')
        print(50*'*')
main()
