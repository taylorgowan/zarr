import zarr
import iris
import iris_grib
import xarray as xr
from numcodecs import Blosc, Zstd
import requests
import copy
import gc
from datetime import datetime, timedelta
from cf_units import Unit
import boto3
import s3fs
import numpy as np


def get_idx_analysis(day,filetype):
    '''
    The Iris package doesn't get all the info we need
    from the HRRR grib2 files. This function queries
    pando to get the .idx files and create a list of
    variable names for the zarr groups we will create
    
    day - Date we are interested in.  This ensures that
          we are getting an accurate variable list
          
    filetype - Either 'sfc' or 'prs' to determine whether 
          we only get the surface data or for all pressure
          levels
    '''
    if filetype == 'sfc':
        fileidx_anl = 'https://pando-rgw01.chpc.utah.edu/hrrr/sfc/'+day+'/hrrr.t00z.wrfsfcf00.grib2.idx'
       
        
    if filetype == 'prs':
        fileidx_anl = 'https://pando-rgw01.chpc.utah.edu/hrrr/prs/'+day+'/hrrr.t00z.wrfprsf00.grib2.idx'
        
    
    # Clean up the text a bit
    # We don't want spaces in our path names
    
    #
    
    idxraw_anl = requests.get(fileidx_anl).text.replace(' ', '_').replace('-','_')
    idxraw_anl = idxraw_anl.replace('_mb','mb').replace('_m_','m_')
    idxraw_anl = idxraw_anl.replace('0_0_day','')
    idxraw_anl = idxraw_anl.replace('_(considered_as_a_single_layer)','_single_layer').replace('_K','K')
    idxraw_anl = idxraw_anl.replace('var_discipline=0_center=7_local_table=1_parmcat=7_parm=204','EFHL')
    idxraw_anl = idxraw_anl.replace('var_discipline=0_center=7_local_table=1_parmcat=7_parm=206','CANGLE')
    idxraw_anl = idxraw_anl.replace('var_discipline=0_center=7_local_table=1_parmcat=7_parm=205','ESP')
    idxraw_anl = idxraw_anl.replace('var_discipline=0_center=7_local_table=1_parmcat=7_parm=200','MNUPHL')
    idxraw_anl = idxraw_anl.replace('var_discipline=0_center=7_local_table=1_parmcat=16_parm=201','unknown')
    idx_anl = idxraw_anl.split(':')
    
    n = len(idx_anl)

    zarr_group = []
    
    for i in range(3,n,6):
        timestep = idx_anl[i+2]
        if timestep == 'anl':
            varname = idx_anl[i+1]+'/'+idx_anl[i]
        else:
            varname = idx_anl[i+1]+'/'+idx_anl[i]+timestep

        zarr_group.append(varname)
    
    return zarr_group

########################################################

def awsAnlIDX(day,filetype):
    '''
    The Iris package doesn't get all the info we need
    from the HRRR grib2 files. This function queries
    pando to get the .idx files and create a list of
    variable names for the zarr groups we will create
    
    day - Date we are interested in.  This ensures that
          we are getting an accurate variable list
          
    filetype - Either 'sfc' or 'prs' to determine whether 
          we only get the surface data or for all pressure
          levels
    '''
    KEY = 'AKIATQDTW4IE5ILHHKTS'
    SECRET = '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'
    REGION = 'us-west-1'

    s3 = boto3.resource(
        service_name='s3',
        region_name=REGION,
        aws_access_key_id=KEY,
        aws_secret_access_key=SECRET
    )
    
    if filetype == 'sfc':
        fileobj_anl = s3.Bucket('hrrr').Object('sfc/'+day+'/hrrr.t00z.wrfsfcf00.grib2.idx').get()['Body']
        
    elif filetype == 'prs':
        fileobj_anl = s3.Bucket('hrrr').Object('prs/'+day+'/hrrr.t00z.wrfprsf00.grib2.idx').get()['Body']
        
        
    # Clean up the text a bit
    # We don't want spaces in our path names
    
    idxraw_anl = fileobj_anl.read().decode('utf-8')
    idxraw_anl = idxraw_anl.replace(' ', '_').replace('-','_')
    idxraw_anl = idxraw_anl.replace('_mb','mb').replace('_m_','m_')
    idxraw_anl = idxraw_anl.replace('0_0_day','')
    idxraw_anl = idxraw_anl.replace('_(considered_as_a_single_layer)','_single_layer').replace('_K','K')
    idxraw_anl = idxraw_anl.replace('var_discipline=0_center=7_master_table=2_parmcat=17_parm=1','LTPINX')
    idxraw_anl = idxraw_anl.replace('var_discipline=0_center=7_local_table=1_parmcat=','')
    idxraw_anl = idxraw_anl.replace('7_parm=204','EFHL').replace('7_parm=206','CANGLE')
    idxraw_anl = idxraw_anl.replace('7_parm=205','ESP').replace('7_parm=200','MNUPHL')
    idxraw_anl = idxraw_anl.replace('var_discipline=0_center=7_local_table=1_parmcat=16_parm=201','unknown')
    idx_anl = idxraw_anl.split(':')
    
    n = len(idx_anl)

    zarr_group = []
    
    for i in range(3,n,6):
        timestep = idx_anl[i+2]
        if timestep == 'anl':
            varname = idx_anl[i+1]+'/'+idx_anl[i]
        else:
            varname = idx_anl[i+1]+'/'+idx_anl[i]+timestep

        zarr_group.append(varname)
    
    return zarr_group


########################################################

def get_idx_forecast(day,hrs):
    
    zarr_group = []
    end = hrs+1
    
    for hr in range(1,end):
    
        fileidx = 'https://pando-rgw01.chpc.utah.edu/hrrr/sfc/'+day+'/hrrr.t00z.wrfsfcf{:02d}.grib2.idx'.format(hr)
        
        h = '{}'.format(hr)
        h_1 = '{}_'.format(hr-1)
    
        idxraw = requests.get(fileidx).text.replace(' ', '_').replace('-','_')
        idxraw = idxraw.replace('_mb','mb').replace('_m_','m_')
        idxraw = idxraw.replace(h_1+h+'_hour','_1hr').replace('0_'+h+'_hour','')
        idxraw = idxraw.replace('_(considered_as_a_single_layer)','_single_layer').replace('_K','K')
        idxraw = idxraw.replace('var_discipline=0_center=7_master_table=2_parmcat=17_parm=1','LTPINX')
        idxraw = idxraw.replace('var_discipline=0_center=7_local_table=1_parmcat=','')
        idxraw = idxraw.replace('7_parm=204','EFHL').replace('7_parm=206','CANGLE')
        idxraw = idxraw.replace('7_parm=205','ESP').replace('7_parm=200','MNUPHL')
        idxraw = idxraw.replace('var_discipline=0_center=7_local_table=1_parmcat=16_parm=201','unknown')
        idx_fcst = idxraw.split(':')
    
        for i in range(3,len(idx_fcst),6):
            timestep = idx_fcst[i+2]
            
            if timestep == h+'_hour_fcst':
                varname = idx_fcst[i+1]+'/'+idx_fcst[i]
            else:
                varname = idx_fcst[i+1]+'/'+idx_fcst[i]+timestep

            zarr_group.append(varname)
    
    return zarr_group

########################################################

def awsFcstIDX(day,hrs):
    
    KEY = 'AKIATQDTW4IE5ILHHKTS'
    SECRET = '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'
    REGION = 'us-west-1'

    s3 = boto3.resource(
        service_name='s3',
        region_name=REGION,
        aws_access_key_id=KEY,
        aws_secret_access_key=SECRET
    )
    
    zarr_group = []
    end = hrs+1
    
    for hr in range(1,end):
        
        fileName = 'sfc/'+day+'/hrrr.t00z.wrfsfcf{:02d}.grib2.idx'.format(hr)
        
        fileobj = s3.Bucket('hrrr').Object(fileName).get()['Body']
        
        h = '{}'.format(hr)
        h_1 = '{}_'.format(hr-1)
        
        idxraw = requests.get(fileidx).text.replace(' ', '_').replace('-','_')
        idxraw = idxraw.replace('_mb','mb').replace('_m_','m_')
        idxraw = idxraw.replace(h_1+h+'_hour','_1hr').replace('0_'+h+'_hour','')
        idxraw = idxraw.replace('_(considered_as_a_single_layer)','_single_layer').replace('_K','K')
        idxraw = idxraw.replace('var_discipline=0_center=7_master_table=2_parmcat=17_parm=1','LTPINX')
        idxraw = idxraw.replace('var_discipline=0_center=7_local_table=1_parmcat=','')
        idxraw = idxraw.replace('7_parm=204','EFHL').replace('7_parm=206','CANGLE')
        idxraw = idxraw.replace('7_parm=205','ESP').replace('7_parm=200','MNUPHL')
        idxraw = idxraw.replace('var_discipline=0_center=7_local_table=1_parmcat=16_parm=201','unknown')
        idx_fcst = idxraw.split(':')
    
        for i in range(3,len(idx_fcst),6):
            timestep = idx_fcst[i+2]
            
            if timestep == h+'_hour_fcst':
                varname = idx_fcst[i+1]+'/'+idx_fcst[i]
            else:
                varname = idx_fcst[i+1]+'/'+idx_fcst[i]+timestep

            zarr_group.append(varname)
    
    return zarr_group

########################################################


def rename_cubes(cubes, names):
    '''
    Function for renmaing cubes using
    the info from the .idx files
    
    names - list of zarr group paths
    cubes - list of grib messages or "cubes"
    '''
    
    for v in range(0, len(names)):
        cubes[v].rename(names[v])
        

    return cubes





def update_coords(cubes):
    '''
    Function to fix the dimensional coordinates.
    Occasionally, after merging, the time coordinate
    will be missing from the dim_coords, so we need
    to create it
    '''
    timecoord = cubes[1].coord('time').points
    fcstcoord = cubes[1].coord('forecast_period').points
    
    time = iris.coords.DimCoord(timecoord, long_name='time', var_name='time',
                                units=Unit('hours since 1970-01-01 00:00:00', calendar='gregorian'))
    
    newtime = timecoord[1:]
    fcstperiod = fcstcoord[1:]
    
    time_1 = iris.coords.DimCoord(newtime, long_name='time_1', var_name='time_1',
                                  units=Unit('hours since 1970-01-01 00:00:00', calendar='gregorian'))
    fcst_1 = iris.coords.AuxCoord(fcstperiod, long_name='forecast_period_1',
                                   var_name='forecast_period_1',
                                  units=Unit('hours since 1970-01-01 00:00:00', calendar='gregorian'))
    
    # This is for the accumulated precipitation variables
    # We skip hours 1, 24, and 48 as they have their own variables
    newtime1 = timecoord[np.r_[1:24,25:47]]
    fcstperiod1 = fcstcoord[np.r_[1:24,25:47]]
    
    time_2 = iris.coords.DimCoord(newtime1, long_name='time_2', var_name='time_2',
                                  units=Unit('hours since 1970-01-01 00:00:00', calendar='gregorian'))
    fcst_2 = iris.coords.AuxCoord(fcstperiod1, long_name='forecast_period_2', var_name='forecast_period_2',
                                 units=Unit('hours since 1970-01-01 00:00:00', calendar='gregorian'))
    
    for cube in cubes:
        
        try:
            cube.remove_coord('height')
        except:
            pass
        
        try:
            cube.remove_coord('pressure')
            
        except:
            pass
        
        

        if cube.shape[0] == 17:
            cube.remove_coord('time')
            cube.add_dim_coord(time_1, 0)
            cube.remove_coord('forecast_period')
            cube.add_aux_coord(fcst_1,0)
            
        elif cube.shape[0] == 45:
            cube.remove_coord('time')
            cube.add_dim_coord(time_2,0)
            cube.remove_coord('forecast_period')
            cube.add_aux_coord(fcst_2,0)
            
        elif cube.shape[0] == 1059:
            try:
                cube.remove_coord('time')
                cube.remove_coord('forecast_period')
                
            except:
                pass
            
        elif cube.shape[0] < 45:
            diffVar = cube.shape[0]
            
            newtime1 = timecoord[0:diffVar]
            fcstperiod1 = fcstcoord[0:diffVar]
    
            time_3 = iris.coords.DimCoord(newtime1, long_name='time_2', var_name='time_2',
                                  units=Unit('hours since 1970-01-01 00:00:00', calendar='gregorian'))
            fcst_3 = iris.coords.AuxCoord(fcstperiod1, long_name='forecast_period_3', var_name='forecast_period_3',
            units=Unit('hours since 1970-01-01 00:00:00', calendar='gregorian'))
            
            cube.remove_coord('time')
            cube.add_dim_coord(time_3,0)
            cube.remove_coord('forecast_period')
            cube.add_aux_coord(fcst_3,0)
            
        else:
            cube.remove_coord('time')
            cube.add_dim_coord(time, 0)
            
            diffVar = 1
            
            
            
        
    gc.collect()
        
    return cubes, diffVar


def cube_to_ds(cube):
    # Convert Iris cube to Xarray Dataset
    return xr.DataArray.from_iris(cube).to_dataset()


def cube_to_da(cube):
    # Convert Iris cube to Xarray DataArray
    return xr.DataArray.from_iris(cube)

def cubelist_to_dalist(cubelist):
    # Convert Iris Cubelist to list of Xarray DataArrays
    dalist = []
    for cube in cubelist:
        dalist.append(cube_to_da(cube))
    return dalist




