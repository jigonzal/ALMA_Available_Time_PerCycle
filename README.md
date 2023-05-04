# ALMA_Available_Time_PerCycle
Quick and dirty Python script to compare the time requested per band and configuration for the ALMA observatory with the time available for that cycle

First you will need to install the modules numpy, matplotlib, astropy and astroplan. 

How to use ir? You have to run the script as follows:

python EstimateAlmaTimePerConfig.py -AOTFile CRISTAL_PROGRAM.aot -TableFile Example_Table.dat -MinimumElevation 30 -SigmaDistributionTime 2

If you run the script as "python EstimateAlmaTimePerConfig.py -h" you will get a short scription of the parameters. 

In short. The script can use a .aot file that has been already submitted and the Program tab has been activated. For those proposals Scheduling Blocks (SB) are generated and can be read by the script. I am not sure if this option will work for all proposal in this cycle but we can check. 

Another option is to provide an ascii table with the representative coordinates of the science goal or SB (you can have multiple SB per science goal), the band, configuration and time requested per SB. 

MinimumElevation parameter is to control how close to the horizon the observations will be allowed. The OT has a intrinsic limit of 20 degree but you can go as high as you want if you want a less elongated beam. The unit is in degrees

SigmaDistributionTime is the sigma value for a Gaussian distribution around the LST of the target. This is the way I am modeling the window of observability of the targets. A smaller SigmaDistributionTime means that all the observations will be done when the source is transiting. A larger SigmaDistributionTime means that some flexibility is added to the observations. 

The LST where the targets can be obvserved correspond to those where the eleveation of the source is above the elevation limit and within the window of observability given by the Gaussian model. 

What are the outputs? You get some plots of elevation and positions on the sky for the first target in the table or in the AOT file. These plots should give an idea of how the time is distributed across LSTs. 

TimePressurePerConfiguration.pdf file has the distribution of time avaliable and requested per configuration and per band. It is a visual aid to check how possible is that your observations will be completed in cycle 10. Additionally, html tables are created for each configuration with the estimated time distribution requested. The table will tell you if you are requested too much time for a given configuration and band. 

This is just a quick and dirty script to estimate the time for ALMA cycle 10. If you get Nans values try checkin that the inputs make sense or if the sources can be above the elevation limit more than 1 hour a day. You can modify the script or send me an email to jgonzalez@carnegiescience.edu for feedback and comments. 

Jorge

