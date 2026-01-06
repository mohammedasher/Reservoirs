import pybamm as pb;import pandas as pd;import numpy as np;
import os, json,openpyxl,traceback,multiprocessing,scipy.optimize,sys
import matplotlib.pyplot as plt;
import pickle,imageio,timeit,random,time, signal
from scipy.io import savemat,loadmat;
from pybamm import constants,exp;import matplotlib as mpl
fs=17; font = {'family' : 'DejaVu Sans','size'   : fs}
mpl.rc('font', **font)

########################     Global settings!!!
rows_per_file = 1;  Scan_end_end = 675
purpose_i = "Full_8Exps_AddLi_Pore" 
# define options:
On_HPC =  True;         Runshort="Reservoir";    Add_Rest = False
Plot_Exp=True;          Timeout=True;     Return_Sol=True;   
Check_Small_Time=True;  R_from_GITT = False
fs = 13; dpi = 100; Re_No =2
Options = [ 
    On_HPC,Runshort,Add_Rest,
    Plot_Exp,Timeout,Return_Sol,
    Check_Small_Time,R_from_GITT,
    dpi,fs]
Timelimit = int(3600*70) # give 48 hours!

Picks_special = [
    1
]



if On_HPC:
    pbs_array_index = os.environ.get("PBS_ARRAY_INDEX", "1")
    i_bundle = Picks_special[int(os.environ["PBS_ARRAY_INDEX"])-1]
    print(f"Running for i_bundle = {i_bundle}")
else:
    i_bundle = 1 # Picks_special[-1]; # manually specify
Scan_start = (i_bundle-1)*rows_per_file+1;    
Scan_end   = min(Scan_start + rows_per_file-1, Scan_end_end)    
purpose = f"{purpose_i}_Case_{Scan_start}_{Scan_end}"
Target  = f'/{purpose}/'
# interpetation: Simnon suggested, with cracking activation, heat transfer
para_csv = f"Bundle_{i_bundle}.csv"  # name of the random file to get parameters

# Path setting:
if On_HPC:                          # Run on HPC
    Path_csv = f"InputData/{purpose_i}/" 
    Path_to_ExpData = "InputData/" 
    BasicPath=os.getcwd() 
    Para_file = Path_csv +  para_csv
else:
    # Add path to system to ensure Fun_P2 can be used
    import sys  
    str_path_0 = os.path.abspath(os.path.join(pb.__path__[0],'..'))
    str_path_1 = os.path.abspath(
        os.path.join(str_path_0,"wip/Rio_Code/Fun_P2"))
    sys.path.append(str_path_1) 
    Path_to_ExpData = os.path.expanduser(
        "~/SimSave/InputData/") # for Linux
    BasicPath =  os.path.expanduser(
        "~/SimSave/Reservoir_Paper")
    Para_file = BasicPath+f'/Get_Random_sets/{purpose_i}/'+para_csv

# import all functions 
from Fun_P2 import * 
from Custom_Para_Func import * # import customized functions used in parameters

# Load input file
Para_dict_list = load_combinations_from_csv(Para_file)
pool_no = len(Para_dict_list) # do parallel computing if needed
# Update 240624: specify ageing sets based on previous calculation
try:
    df_ageset = pd.read_csv(BasicPath+f'/Get_Random_sets/{purpose_i}/'+"Age_set.csv")
    age_set_pre = int(df_ageset.loc[df_ageset['Pick_i'] == i_bundle, 'age_set'].iloc[0])
    # update 240701: Reset Re_No here:
    Re_No = int(df_ageset.loc[df_ageset['Pick_i'] == i_bundle, 'Re_run'].iloc[0])
except:
    age_set_pre = []
else:
    pass

Options.append(age_set_pre)

midc_merge_all = [];  Sol_RPT_all = [];  Sol_AGE_all = []
Path_List = [BasicPath, Path_to_ExpData,Target,purpose] 
# Run the model
if Re_No == 0:
    midc_merge,Sol_RPT,Sol_AGE,DeBug_Lists = Run_P2_Excel (
        Para_dict_list[0], Path_List, 
        Re_No, Timelimit, Options) 
elif Re_No > 0:
    Scan_i = int(Para_dict_list[0][ 'Scan No'] )
    with open(
        BasicPath + Target+"Mats/" + str(Scan_i)
        + f'_Re_{Re_No-1}-Save_for_Reload.pkl', 'rb') as file:
        Save_for_Reload = pickle.load(file)
    midc_merge,Sol_RPT,Sol_AGE,DeBug_Lists = Run_Restart (
        Save_for_Reload,  Para_dict_list[0],  Path_List, 
        Re_No,      Timelimit, Options)