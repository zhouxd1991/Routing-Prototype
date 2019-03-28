# The main programe for creating the routing map (with HTU)
# The module has to be divided into several parts 
# 1. divide the basin to HTUs
# 2. define the dams and river HTUs
# 3. run the routing model if runoff/drainage data are provided. 
#
from run_def import *
import os
import sys
sys.path.append('/home/xzhou/ROUTING_ISO/prototype_thesis/')
sys.path.append('/home/xzhou/ROUTING_ISO/prototype_thesis/subroutines')
sys.path.append('/home/xzhou/ROUTING_ISO/prototype_thesis/sub_routing')
import findbasins
from routines import *       # all definitions are defined in routines.py

print 'Processing running: ', basin_name + '-'+ mode
if mode.find('Zun') > -1 and styear <1982:
    styear=1982

start_time=datetime.datetime.now()
cwd = os.getcwd()

surfmap = '/home/xzhou/ROUTING_ISO/prototype/surfmap/'
codefolder = '/home/xzhou/ROUTING_ISO/prototype/sub_routing/'
# ============================================================================
# The folder and files for a few inputs.  
'''
irrigfile = '/homedata/xzhou/irrig/irrig_CN_1km.nc'         # the irrigation location map 
pathfolder = '/bdd/MEDI/workspaces/xzhou/Pathway/CN_FINAL/'                           # the Pathway folder
pathfile=pathfolder + 'path.nc'                                                     # the path file
grid_area_file = "/homedata/xzhou/stations/gridarea.nc"         # the area for each 0.008333degree grids
file_rt = '/bdd/MEDI/workspaces/xzhou/Pathway/CN_CN/routing.nc'     # routing file at 1km, with topoinx for each 1km
'''
irrig_file = surfmap+'irrig_CN_1km.nc'
trip_file = surfmap+'trip.nc'
fac_file = surfmap+'fac.nc'
orog_file = surfmap+'orog.nc'
path_file = surfmap+'path.nc'
grid_area_file = surfmap+'gridarea.nc'
file_rt = surfmap + 'routing.nc'
distance_file = surfmap + 'route_distance.nc'
celer_file = surfmap + 'celer_asia_r.nc'
# ============================================================================
# create the output folder
#
# the restart file for HTUs
restart_file = outfolder+ 'restart-'+basin_name+'.nc'

command = 'mkdir '+ outfolder + '/' + basin_name + '-' + mode[:-1]
os.system(command)
outfolder = outfolder + '/' + basin_name + '-' + mode[:-1] + '/'
command = 'cp ' + cwd + '/*.f90 ' + outfolder
os.system(command)
command = 'cp ' + cwd + '/*.py ' + outfolder
os.system(command)
command = 'cp ' + codefolder + '/*.f90 ' + outfolder
os.system
# 
figfolder = outfolder + '/figure/'
os.system('mkdir '+ figfolder)
#
pdf = PdfPages(figfolder+'figures.pdf')
#
# ============================================================================
# package clip : to extract the area identified by the boundaries of lat and lon
# file_rt
lat, lon = clip.CHECK(SOUTH_NORTH, WEST_EAST, file_rf, full=2 )
lat_forcing, lon_forcing = clip.CHECK(SOUTH_NORTH, WEST_EAST, file_rf, full=2 )
res_rf = lon[1] - lon[0]
nbpt = len(lat) * len(lon)
iim = len(lat)
jjm = len(lon)
#
restart_file = outfolder+ '/restart-'+basin_name+'.nc'
if restart_id ==1 :
    try:
        o = Dataset(restart_file,'r')
        lat_restart = o.variables['lat'][:]
        lon_restart = o.variables['lon'][:]
        # 
        nbasmax=len(o.dimensions['nbasmax'])
        nbpt=len(o.dimensions['nbpt'])
        if ( lat_restart[0] == lat[0] and lon_restart[0]==lon[0] and nbpt_restart == nbpt and nbasmax_given==nbasmax):
            restart_id = 1
            print 'Restart file existed.'
        else:
            print 'The restart file is not corresponding to your setting'
            print 'Restart file setting', lat_restart[0], lon_restart[0] , nbpt_restart, nbasmax_given
            print 'Your new setting',lat[0], lon[0], nbpt, nbasmax
            print 'The model will run, but the variables will be saved in file ', restart_file
            restart_id = 0
        o.close()
    except RuntimeError:
        print 'The restart file does not exist'
        print 'The model will run, but the variables will be saved in file ', restart_file
        restart_file = restart_file
        restart_id=0
        nbasmax = nbasmax_given
else:
    nbasmax = nbasmax_given
#
# =================================================================
# Part 1. prepare the coordinates and the index 
# 
# adjust the boundary of the regions 
WEST_EAST = [lon[0]-res_rf/2. , lon[-1]+res_rf/2.]
SOUTH_NORTH = [lat[0]-res_rf/2. , lat[-1]+ res_rf/2.]

# The ranges for high-resolution routing map (1km)
lat_rt, lon_rt = clip.CHECK(SOUTH_NORTH, WEST_EAST, path_file,  full=2)

# The ranges for high-resolution map (1km) for plotting. 
lat_rt_plot = np.zeros(len(lat_rt))
for i in range(len(lat_rt)):
    lat_rt_plot[i] = round(lat_rt[i],4)
lon_rt_plot = np.zeros(len(lon_rt))
for i in range(len(lon_rt)):
    lon_rt_plot[i] = round(lon_rt[i],4)

# the ratio of two resolution
res_rt= lon_rt[1]-lon_rt[0]
r_ratio = int( res_rf / res_rt + 0.1)
# 
# VERY IMPORTANT: create the index
latmin = np.zeros(nbpt,'f')
latmax = np.zeros(nbpt,'f')
lonmin = np.zeros(nbpt,'f')
lonmax = np.zeros(nbpt,'f')
index = np.zeros((iim,jjm),'d')
ib = 0
for i in range(iim):
    for j in range(jjm):
        latmin[ib] = r_ratio*i
        latmax[ib] = r_ratio*(i+1)
        lonmin[ib] = r_ratio*j
        lonmax[ib] = r_ratio*(j+1)
        index[i,j] = ib
        ib += 1

# findig: find the number of grid according to index of (lat,lon)
findig=np.zeros((len(lat_rt),len(lon_rt)),'d')
for it in range(len(lat_rt)):
    for jt in range(len(lon_rt)):
        i = int(it/r_ratio)
        j = int(jt/r_ratio)
        findig[it,jt] = index[i,j]

# End of part 1
# =================================================================
#
# Part 2: Deal with the HTUs, divide the HTUs according to the water velocity

basins_area = clip.CLIPFILE(lat_rt, lon_rt, grid_area_file, 'area') * 1000000.  # m^2 at 0.008333 degree (1km) resolution

if restart_id == 0 :
    # Read the necessary data
    trip = clip.CLIPFILE(lat_rt, lon_rt, trip_file, 'trip')
    trip[trip>10000]=-1
    t1, t2,fac = clip.CLIPFILE(lat_rt, lon_rt, fac_file, 'fac', 3)
    t1, t2, orog = clip.CLIPFILE(lat_rt, lon_rt, orog_file , 'orog', 3)
    t1, t2, distance = clip.CLIPFILE(lat_rt, lon_rt, distance_file ,'distance', 3)

    t1, t2, celer = clip.CLIPFILE_bias(lat_rt, lon_rt, celer_file, 'celer', 3, lat_bias=2)
    # check the celerity with trip
    '''
    pdf1 = PdfPages(figfolder+'check_celer.pdf')
    a1 = 0; a2 = 50
    b1 = 0; b2= 50
    lat_rt_p = lat_rt[a1:a2]
    lon_rt_p = lon_rt[b1:b2]
    mapplot0408.MAPPLOT_pdf(pdf1, lat_rt_p, lon_rt_p, trip[a1:a2,b1:b2])
    mapplot0408.MAPPLOT_pdf(pdf1, lat_rt_p, lon_rt_p, fac[a1:a2,b1:b2])
    mapplot0408.MAPPLOT_pdf(pdf1, lat_rt_p, lon_rt_p, celer[a1:a2,b1:b2])

    pdf1.close()
    #sys.exit()
    '''

    # There is a shift (2 grid cells, 2km) in the lat coordinates in the celerity data
    # which results in the inconsistence with the direction data. Not found the reason!
    celer[celer.recordmask]=-99.9
    topoindex = clip.CLIPFILE(lat_rt, lon_rt, file_rt, 'topoind')
    #trip[((fac==1) & (trip==98))] = -1
    #fac[((fac==1) & (trip==98))] = 0
    #
    # =================================================================
    #
    # STEP 1: prepare the river and dams
    #
    # ==================================================================
    # Read the dams first, because those HTUs with dams (might less than the criterion of river) should be considered as rivers
    # only the dams grid need to search upstream basins
    # read the dams locations, record the characteristics of the dams.  
    dam_lat = []; dam_lon=[]; 
    if do_dams == 1:
        o = Dataset(dam_file, 'r')
        try :
            nbdams=len(o.dimensions['dams'])
        except:
            nbdams=len(o.dimensions['stations'])

        nb = 0
        dam_no=[]
        dam_latitude = []; dam_longitude=[]
        tmp_dam_vref=[] ; tmp_dam_area=[]; tmp_dam_year=[]
        dam_height = []

        for i in range(nbdams):
            tmplat = o.variables['LAT_DD'][i]
            tmplon = o.variables['LONG_DD'][i]
            if (tmplat >= lat_rt[0] and tmplat<= lat_rt[-1] and tmplon >= lon_rt[0] and tmplon<= lon_rt[-1]):
                tmp = np.abs(lat_rt - tmplat)
                ilat = np.where(tmp == np.min(tmp))[0]
                tmp = np.abs(lon_rt - tmplon)
                ilon = np.where(tmp == np.min(tmp))[0]
                dam_lat.append(ilat)
                dam_lon.append(ilon)
                dam_latitude.append(tmplat)
                dam_longitude.append(tmplon)

                #dam_no.append(o.variables['GRAND_ID'][i])
                tmp_dam_vref.append(o.variables['CAP_REP'][i] * 1000000* 1000)  # 10^6 m^3 => kg
                tmp_dam_area.append(o.variables['AREA_REP'][i])
                tmp_dam_year.append(o.variables['YEAR'][i])
                try:
                    dam_height.append(o.variables['DAM_HGT_M'][i])
                except:
                    dam_height.append(o.variables['HGT_DAM_M'][i])
                nb = nb+1
        o.close()
        print "The number of dams in the domain", nb
    #
    # If the dam is not activated or there are no dams in the selected basin
    # we assume their is one dam at (0,0), which makes the variable dam_lat and dam_lon can be passed to the fortran subroutine.
    if len(dam_lat) == 0:
        do_dams = 0 
        nbdams = 1
        dam_lat.append(0)
        dam_lon.append(0)
    #
    # =============================================================
    # The main codes
    # 
    # Find and clip HTUs 
    # debug=1 #
    # tr: this is an additional shift of cumulative time when dividing the HTUs
    # Because the division condition is that grid cells within an increment of time unit (1h) are grouped into the same HTU.
    # While because the cumulative time is not integer, for example, the previous cumulative time is 0.68h, then the upstream
    # cell can be 1.13h (0.68+0.45). In this case, the cumulative time in all HTUs are shorter than 1h, 
    # which will delay the river routing
    # If tr is larger than 1, for instance 1.2, the grid cells with cumulative time less than 1.2 will be grouped,
    # and that larger than 1.2 are grouped to the other HTU. Then the average cumulative time could be higher, and closer
    # to the given routing time unit (1h).

    # It is set as 1.0 now, so the function is deactivated. 
    tr = 1.
    #
    # output files:
    # SIG_TIME: the flow time in a grid cell (1km)
    # CUM_TIME: the flow time to the boundaries of a grid (0.5dg or 0.25dg or 0.1dg)
    # attention: here we can export sig_time : the time to travel in a single HTU
    # the cum_time here calculates the cumulative time out of the gridbox,rather than the HTU!!

    print "print start doing findbasins"
    celer, basins, lsoutlet, SIG_TIME,  CUM_TIME, route_togrid, route_tobasin, HTU_lsriver, lsdam= \
            findbasins.findbasins_main( nbasmax, routing_time_step, debug, do_dams, fac_criterion, tr, do_pfa, maxpercent, trip, fac,\
            basins_area, celer, latmin, latmax, lonmin, lonmax, findig, dam_lat, dam_lon)
    #celer, basins, lsoutlet, SIG_TIME,  CUM_TIME, route_togrid, route_tobasin, HTU_lsriver, lsdam= \
    #        findbasins.findbasins_main( nbasmax, routing_time_step, debug, do_dams, fac_criterion, tr, trip, fac,\
    #        basins_area, celer, latmin, latmax, lonmin, lonmax, findig, dam_lat, dam_lon)

    '''
    tmpbasins = basins - 1
    tmp = tmpbasins * 1.0
    tmp[tmp<=-1]=np.nan
    pdf1 = PdfPages(figfolder+'figure_basin.pdf')
    MAPPLOT_pdf(pdf1, lat_rt, lon_rt, tmp, lsoutlet, trip, title='basins', lat_int = 0.25, lon_int=0.25, SN=[26,26.5], WE=[103.5, 104])
    tmp = fac
    tmp[tmp<1000] = np.nan
    tmp1 = lsoutlet
    tmp1[:] = -1
    MAPPLOT_pdf(pdf1, lat_rt, lon_rt, tmp, tmp1, trip, title='basins', lat_int = 0.25, lon_int=0.25, SN=[26,26.5], WE=[103.5,104])
    print "print figure_basin done"
    pdf1.close()
    #sys.exit()
    '''
    # If the maximum number of HTUs in all grids are less than the given number. 
    nbasmax = np.nanmax(basins)
    route_togrid = route_togrid[:,:nbasmax]
    route_tobasin = route_tobasin[:,:nbasmax]
    HTU_lsriver = HTU_lsriver[:,:nbasmax]
    
    # Calculatie the HTU diagnostics. 
    # HTU_cumtime >> the cumulative time of HTU, the maximum time for a grid cell flowing out of the HTU
    HTU_fac, HTU_area, HTU_distance, HTU_cumtime, topo_resid, topoindex_h= \
            findbasins.htu_cal( nbasmax, debug, latmin, latmax, lonmin, lonmax, \
            basins, lsoutlet, fac, trip, basins_area, distance, SIG_TIME, CUM_TIME, topoindex) 
    
    #
    print 'total number of grids', nbpt
    print 'Finished finding basins, maximum basins', np.nanmax(basins)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # if needs to plot the maps (for checking) 
    if do_plotmap:
        tmp = basins * 1.0
        tmp[tmp<=-1]=np.nan
        #CUM_TIME[np.isnan(tmp)] = np.nan
        tmp1 = np.copy(lsoutlet)
        tmp1[:] = 0

        if basins.shape[0]<400:
            #MAPPLOT(figfolder + 'basin-new_2', lat_rt, lon_rt, tmp, lsoutlet, trip)
            #MAPPLOT(figfolder + 'cum-time-new_2', lat_rt, lon_rt, CUM_TIME, lsoutlet, trip)
            #MAPPLOT(figfolder + 'sig-time-new_2', lat_rt, lon_rt, SIG_TIME, lsoutlet, trip)
            #MAPPLOT(figfolder + 'celer', lat_rt, lon_rt, celer, lsoutlet, trip)
            MAPPLOT_pdf(pdf, lat_rt, lon_rt, tmp, lsoutlet, trip, title='basins')
            MAPPLOT_pdf(pdf, lat_rt, lon_rt, CUM_TIME, lsoutlet, trip, title='cum_time')
            MAPPLOT_pdf(pdf, lat_rt, lon_rt, SIG_TIME, lsoutlet, trip, title='sig_time')
            MAPPLOT_pdf(pdf, lat_rt, lon_rt, celer, lsoutlet, trip, title='celer')
        else:
            MAPPLOT_pdf(pdf, lat_rt, lon_rt, tmp, tmp1, trip, title='basins')
            MAPPLOT_pdf(pdf, lat_rt, lon_rt, CUM_TIME, tmp1, trip, title='cum_time')
            MAPPLOT_pdf(pdf, lat_rt, lon_rt, SIG_TIME, tmp1, trip, title='sig_time')
            MAPPLOT_pdf(pdf, lat_rt, lon_rt, celer, tmp1, trip, title='celer')
            MAPPLOT_pdf(pdf, lat_rt, lon_rt, fac, tmp1, trip, title='fac')
            #MAPPLOT(figfolder + 'basin-new_2', lat_rt, lon_rt, tmp, tmp1, trip)
            # CUM_TIME >> time to grid boundary
            #MAPPLOT(figfolder + 'cum-time-new_2', lat_rt, lon_rt, CUM_TIME, tmp1, trip)
            #MAPPLOT(figfolder + 'pixel-time-new_2', lat_rt, lon_rt, pixel_time, tmp1, trip)
            #MAPPLOT(figfolder + 'sig-time-new_2', lat_rt, lon_rt, SIG_TIME, tmp1, trip)
            #MAPPLOT(figfolder + 'celer', lat_rt, lon_rt, celer, tmp1, trip)

        #MAPPLOT(figfolder + 'basin-new_2', lat_rt, lon_rt, tmp, tmp1, trip)
        #MAPPLOT(figfolder + 'cum-time-new_2', lat_rt, lon_rt, CUM_TIME, tmp1, trip)
        #tmp = plothtu1(figfolder+'HTU_cumtime', lat_rt, lon_rt, basins-1, findig, HTU_cumtime)

        print 'test plotting finished'
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # 
    # If do_clip, the coordinates in outlet_lat, outlet_lon will be the outlets that identify the basins
    # However, the outlet_lat, outlet_lon are not necessary the points where the discharge should write out.
    # 
    if (do_clip and len(outlet_lat) >= 1) : 
        basin_clip, grid_clip, route_togrid, route_tobasin = findbasins.clipbasin(0, outlet_lat, outlet_lon, lat_rt, lon_rt, \
                findig, basins, route_togrid, route_tobasin,HTU_fac)
        HTU_lsriver = HTU_lsriver * basin_clip
        HTU_area = HTU_area * basin_clip
        HTU_fac = HTU_fac * basin_clip
        #
    else:
        basin_clip = np.zeros((nbpt, nbasmax),'f')
        basin_clip[HTU_area>0]=1
        grid_clip = np.zeros((len(lat_rt),len(lon_rt)),'f')
        grid_clip[:] = 1

    # finished the HTU numbering and then transfer it from Fortran to Python, minus 1 to fit the python number system
    basins = basins -1 
    lsoutlet = lsoutlet - 1
    route_togrid = route_togrid[:,:nbasmax]-1
    route_tobasin = route_tobasin[:,:nbasmax]-1
    HTU_lsriver = HTU_lsriver[:,:nbasmax]
    #
    # If dam regulation is considered, the dams that are not located in the basin are excluded. 
    # If there are dams, the rivers which connect to dams should be adjusted 
    #
    plot_dam_lat = []
    plot_dam_lon = []
    if do_dams:
        lsdam_ig = np.zeros(nbdams,'i')
        lsdam_ib = np.zeros(nbdams,'i')
        lsdam_id = np.zeros(nbdams,'i')
        #
        lsdam_ig=[]
        lsdam_ib = []
        lsdam_id = []
        dam_vref = []
        dam_year = []
        dam_area = []
        nb = 0
        for i in range(lsdam.shape[0]):
            for j in range(lsdam.shape[1]):
                if lsdam[i,j]>0:    
                    ig = findig[i,j]
                    ib = basins[i,j]
                    id = lsdam[i,j]-1
                    if HTU_area[ig,ib] > 0 :
                        lsdam_id.append(nb)
                        lsdam_ig.append(findig[i,j])
                        lsdam_ib.append(basins[i,j])
                        dam_vref.append(tmp_dam_vref[id])
                        dam_year.append(tmp_dam_year[id])
                        dam_area.append(tmp_dam_area[id])
                        plot_dam_lat.append(dam_latitude[id])
                        plot_dam_lon.append(dam_longitude[id])
                        nb += 1
        nbdams = nb 

        print "dams in region", nbdams

        # 2. Check if the dam is on the river, if not, shift the dams to a river HTU or out of the grid
        #
        for i in range(nbdams):
            ig = lsdam_ig[i]
            ib = lsdam_ib[i]
            if (HTU_lsriver[ig,ib] <1):  
                HTU_lsriver[ig,ib] = 1
                igt = route_togrid[ig,ib]
                ibt = route_tobasin[ig,ib]
                while (igt >= 0 and HTU_lsriver[igt,ibt] < 1):
                    # search until finding a river
                    HTU_lsriver[igt,ibt]=1
                    ig = igt; ib=ibt
                    igt = route_togrid[ig,ib]
                    ibt = route_tobasin[ig,ib]

        #print "sum of lsriver", np.sum(HTU_lsriver)
        #print 'max in HTU_fac' , np.max(HTU_fac)

    # find the upstream HTU linkage
    # which is used to propagate the downstream water demand
    nbne = 4
    upstream_htu = findbasins.updown( nbne, 0,  HTU_lsriver, route_togrid+1, route_tobasin+1, HTU_fac) 
    upstream_htu = upstream_htu - 1
    #
    # ==================================================================
    # The rivers connect to dams should be a separate HTU
    if do_dams: 
        # for dams in the source, add the upstream_htu of the dams (to get the inflow)
        for i in range(nbdams):
            ig = lsdam_ig[i]
            ib = lsdam_ib[i]
            tmp = 0

            nb = 0
            #if np.max(upstream_htu[ig,ib,:,:]) < 1: # no upstream HTU
            for igt in range(nbpt):
                for ibt in range(nbasmax):
                    if route_togrid[igt,ibt]==ig and route_tobasin[igt,ibt]==ib :
                        if HTU_fac[igt,ibt] > tmp:
                            tmp = HTU_fac[igt,ibt]
                            upstream_htu[ig,ib,nb,0]=igt
                            upstream_htu[ig,ib,nb,1]=ibt
                            nb = nb + 1
                            HTU_lsriver[igt,ibt] = 1
            #print "dam ", i, ig,ib, HTU_lsriver[ig,ib],HTU_fac[ig,ib],upstream_htu[ig,ib,:,:]

    # ===============================================
    # save the variables into netcdf file, in case they will be used repectly
    print "Writing the basins file"
    o = Dataset(restart_file,'w')
    o.history = "Temporary save of the files"
    o.createDimension('lat', iim)
    o.createDimension('lon', jjm)
    o.createDimension('n_lat_rt', len(lat_rt))
    o.createDimension('n_lon_rt', len(lon_rt))
    o.createDimension('nbpt', nbpt)
    o.createDimension('nbasmax', nbasmax)
    o.createDimension('nbne', nbne)
    o.createDimension('layer2', 2)
    if do_dams:
        o.createDimension('nbdams', nbdams)

    #
    vlat=o.createVariable('lat','f', ('lat',))
    vlat.units='degree_north'
    vlat.longname='latitude'
    vlat[:]=lat[:]
    # 
    vlon=o.createVariable('lon','f', ('lon',))
    vlon.units='degree_east'
    vlon.longname='longitude'
    vlon[:]=lon[:]
    #
    vlat_rt=o.createVariable('lat_rt','f', ('n_lat_rt',))
    vlat_rt.units='degree_north'
    vlat_rt.longname='latitude'
    vlat_rt[:]=lat_rt[:]
    # 
    vlon_rt=o.createVariable('lon_rt','f', ('n_lon_rt',))
    vlon_rt.units='degree_east'
    vlon_rt.longname='longitude'
    vlon_rt[:]=lon_rt[:]
    #
    vbasins = o.createVariable('basins', 'f', ('n_lat_rt','n_lon_rt',))
    vbasins.longname = 'Number of the basins'
    vbasins[:] = basins[:]    
    #
    vbasins = o.createVariable('findig', 'f', ('n_lat_rt','n_lon_rt',))
    vbasins.longname = 'IG of the the basins'
    vbasins[:] = findig[:]    
    #
    vbasins = o.createVariable('celer', 'f', ('n_lat_rt','n_lon_rt',))
    vbasins.longname = 'celerity of the grids'
    vbasins[:] = celer[:]    
    #
    vbasins = o.createVariable('fac', 'f', ('n_lat_rt','n_lon_rt',))
    vbasins.longname = 'fac of the grids'
    vbasins[:] = fac[:]    
    #
    vbasins = o.createVariable('orog', 'i', ('n_lat_rt','n_lon_rt',))
    vbasins.longname = 'elevation of the grids(DEM)'
    vbasins[:] = orog[:]    
    #
    #vbasins_orig = o.createVariable('basins_orig', 'f', ('n_lat_rt','lon_rt',))
    #vbasins_orig.longname = 'Number of the basins'
    #vbasins_orig[:] = basins_orig[:]    
    #
    voutlet = o.createVariable('outlet', 'f', ('n_lat_rt','n_lon_rt',))
    voutlet.longname = 'outlet of the basins'
    voutlet[:] = lsoutlet[:]    
    #
    vgrid_clip = o.createVariable('grid_clip', 'f', ('n_lat_rt','n_lon_rt',))
    vgrid_clip.longname = 'proportion in the domain'
    vgrid_clip[:] = grid_clip[:]    
    #
    #topo_resid_lalo = np.zeros((nbasmax,iim,jjm),'f')
    #HTU_fac_lalo = np.zeros((nbasmax,iim,jjm),'f')
    #HTU_distance_lalo = np.zeros((nbasmax,iim,jjm),'f')
    #HTU_orog_lalo = np.zeros((nbasmax,iim,jjm),'f')
    HTU_area_lalo = np.zeros((nbasmax,iim,jjm),'f')
    HTU_lsriver_lalo = np.zeros((nbasmax,iim,jjm),'d')
    #route_togrid_lalo = np.zeros((nbasmax,iim,jjm),'d')
    #route_tobasin_lalo = np.zeros((nbasmax,iim,jjm),'d')
    for i in range(nbasmax):
        #topo_resid_lalo[i,:,:] = grid_lalo(topo_resid[:,i])
        #HTU_fac_lalo[i,:,:] = grid_lalo(HTU_fac[:,i])
        #HTU_distance_lalo[i,:,:] = grid_lalo(HTU_distance[:,i])
        #HTU_orog_lalo[i,:,:] = grid_lalo(HTU_orog[:,i])
        HTU_area_lalo[i,:,:] = grid_lalo(HTU_area[:,i], iim, jjm)
        HTU_lsriver_lalo[i,:,:] = grid_lalo(HTU_lsriver[:,i], iim, jjm)
        #route_togrid_lalo[i,:,:] = grid_lalo(route_togrid[:,i])
        #route_tobasin_lalo[i,:,:] = grid_lalo(route_tobasin[:,i])

    v = o.createVariable('topo_resid','f',('nbpt','nbasmax',))
    v.units='-'
    v.longname='Basin topographic index for each HTU'
    v[:] = topo_resid[:]
    #
    #v = o.createVariable('topo_resid_lalo','f',('nbasmax','lat','lon'))
    #v.units='-'
    #v.longname='Basin topographic index for each lalo grid'
    #v[:] = topo_resid_lalo[:]
    #
    v = o.createVariable('HTU_lsriver','f',('nbpt','nbasmax',))
    v.units='-'
    v.longname='If the HTU is a river'
    v[:] = HTU_lsriver[:]
    #
    v = o.createVariable('HTU_lsriver_lalo','f',('nbasmax','lat','lon'))
    v.units='-'
    v.longname='If the grid is a river'
    v[:] = HTU_lsriver_lalo[:]
    #
    v = o.createVariable('HTU_fac','f',('nbpt','nbasmax',))
    v.units='-'
    v.longname='Largest fac value in a HTU'
    v[:] = HTU_fac[:]
    #
    v = o.createVariable('HTU_cumtime','f',('nbpt','nbasmax',))
    v.units='-'
    v.longname='Cumulative time of the HTU'
    v[:] = HTU_cumtime[:]
    #
    #v = o.createVariable('HTU_fac_lalo','f',('nbasmax','lat','lon'))
    #v.units='-'
    #v.longname='Largest fac value in a HTU at a grid'
    #v[:] = HTU_fac_lalo[:]
    #
    v = o.createVariable('HTU_distance','f',('nbpt','nbasmax',))
    v.units='km'
    v.longname='Distance to the oulet in a HTU'
    v[:] = HTU_distance[:]
    #
    #v = o.createVariable('HTU_distance_lalo','f',('nbasmax','lat','lon'))
    #v.units='km'
    #v.longname='Distance to the outlet HTU at a grid'
    #v[:] = HTU_distance_lalo[:]
    #
    #v = o.createVariable('HTU_orog','f',('nbpt','nbasmax',))
    #v.units='m'
    #v.longname='Average elevation of the HTU'
    #v[:] = HTU_orog[:]
    #
    #v = o.createVariable('HTU_orog_lalo','f',('nbasmax','lat','lon'))
    #v.units='m'
    #v.longname='Average elevation of the HTU at a grid'
    #v[:] = HTU_orog_lalo[:]
    #
    v = o.createVariable('HTU_area','f',('nbpt','nbasmax',))
    v.units='m^2'
    v.longname='Area for each HTU'
    v[:] = HTU_area[:]
    #
    v = o.createVariable('HTU_area_lalo','f',('nbasmax','lat','lon'))
    v.units='m^2'
    v.longname='Area for each HTU at a grid'
    v[:] = HTU_area_lalo[:]
    #
    v = o.createVariable('route_togrid','d',('nbpt','nbasmax',))
    v.longname='routing to next grid'
    v[:]=route_togrid[:]
    #
    #v = o.createVariable('route_togrid_lalo','d',('nbasmax','lat','lon'))
    #v.longname='routing to next grid'
    #v[:]=route_togrid_lalo[:]
    #
    v = o.createVariable('basin_clip','d',('nbpt','nbasmax',))
    v.longname='basins belong to the defined outlets'
    v[:]=basin_clip[:]
    #
    v = o.createVariable('route_tobasin','d',('nbpt','nbasmax',))
    v.longname='routing to next basin'
    v[:]=route_tobasin[:]
    #
    #v = o.createVariable('route_tobasin_lalo','d',('nbasmax','lat','lon'))
    #v.longname='routing to next basin'
    #v[:]=route_tobasin_lalo[:]
    #
    v = o.createVariable('upstream_htu','d',('nbpt','nbasmax','nbne','layer2'))
    v.longname='upstream htu of the current htu'
    v[:]=upstream_htu[:]
    #
    if do_dams:
        v = o.createVariable('lsdam_id','i', ('nbdams')) 
        v[:] = lsdam_id[:]

        v = o.createVariable('lsdam_ig','i', ('nbdams')) 
        v[:] = lsdam_ig[:]

        v = o.createVariable('lsdam_ib','i', ('nbdams')) 
        v[:] = lsdam_ib[:]

        v = o.createVariable('dam_vref','f', ('nbdams')) 
        v[:] = dam_vref[:]

        v = o.createVariable('dam_year','f', ('nbdams')) 
        v[:] = dam_year[:]

        v = o.createVariable('dam_area','f', ('nbdams')) 
        v[:] = dam_area[:]

    o.close()
else:
    print "   >> reading data from restart file"
    o = Dataset(restart_file, 'r')
    basins = o.variables['basins'][:]
    #basins_orig = o.variables['basins_orig'][:]
    lsoutlet = o.variables['outlet'][:]
    grid_clip = o.variables['grid_clip'][:]
    topo_resid = o.variables['topo_resid'][:]
    HTU_fac = o.variables['HTU_fac'][:]
    #HTU_orog = o.variables['HTU_orog'][:]
    HTU_area = o.variables['HTU_area'][:]
    route_togrid = o.variables['route_togrid'][:]
    route_tobasin = o.variables['route_tobasin'][:]
    basin_clip  = o.variables['basin_clip'][:]
    HTU_lsriver  = o.variables['HTU_lsriver'][:]
    upstream_htu = o.variables['upstream_htu'][:]
    if do_dams:
        nbdams = len(o.dimensions['nbdams'])
        print ">>>>>>>>nbdams",nbdams
        lsdam_id = o.variables['lsdam_id'][:]
        lsdam_ig = o.variables['lsdam_ig'][:]
        lsdam_ib = o.variables['lsdam_ib'][:]
        dam_vref = o.variables['dam_vref'][:]
        dam_year = o.variables['dam_year'][:]
        dam_area = o.variables['dam_area'][:]

    o.close()

# calculate the area of the grid, which is used to sum the total amount of water as runoff and drainage.  
total_area_grid = np.sum(HTU_area,axis=1)
total_area_lalo = grid_lalo( np.sum(HTU_area,axis=1), iim, jjm) 
total_area_lalo[total_area_lalo==0]= 1.e-10

'''
if do_plotmap:
    tmp = plothtu1(figfolder+'HTU-2', lat_rt, lon_rt, basins, findig, HTU_fac, 1)
    tmp = plothtu1(figfolder+'HTU_area', lat_rt, lon_rt, basins, findig, HTU_area, 1)
    tmp = plothtu1(figfolder+'topoind', lat_rt, lon_rt, basins, findig, topo_resid, 0)
    rivermask = plothtu1(figfolder+'rivers', lat_rt, lon_rt, basins, findig, HTU_lsriver, 1)
    #MAPPLOT1(figfolder+"river-dam", lat_rt_plot, lon_rt_plot, rivermask)
    #basin_clip1 = plothtu(figfolder+'basin_clip', lat_rt, lon_rt, basin_clip, 1)
    #MAPPLOT1(figfolder+"river-dam", lat_rt_plot, lon_rt_plot, rivermask)
'''

pdf.close()
#sys.exit()
# ===================================================================
#
# STEP 2 link the irrigation grids to the river HTUs
# 
# ===================================================================
# Ir_points: the location of irrigation points 
a1,a2,Ir_points = clip.CLIPFILE([lat_rt[0]-0.008333, lat_rt[-1]], lon_rt, irrig_file, 'irrig',full=3)

a1,a2,Ir_points = clip.CLIPFILE(lat_rt, lon_rt, irrig_file, 'irrig',full=3)
Ir_points[np.isnan(Ir_points)]=0
#
# Ir_togrid and Ir_tobasin show the linkage to the HTUs (ig,ib)
Ir_togrid = np.zeros((len(lat_rt), len(lon_rt)),'d')
Ir_togrid[:]=-1
Ir_tobasin = np.zeros((len(lat_rt), len(lon_rt)),'d')

print "Ir_points size", Ir_points.shape
print "Ir_tobasin size", Ir_togrid.shape, Ir_tobasin.shape

# read the cost, linkage of the grid cells from independent file 
a1,a2,COST = clip.CLIPFILE(lat_rt, lon_rt, path_file, 'cost',full=3)
lat_index = clip.CLIPFILE(lat_rt, lon_rt, path_file, 'lat_index')
lon_index = clip.CLIPFILE(lat_rt, lon_rt, path_file, 'lon_index')
ab_lat = clip.CLIPFILE(lat_rt, lon_rt, path_file, 'ab_lat')
ab_lon = clip.CLIPFILE(lat_rt, lon_rt, path_file, 'ab_lon')
#
tmp = np.zeros((len(lat_rt), len(lon_rt)), 'd')
tmploc= np.zeros((len(lat_rt),len(lon_rt),2),'d')
for i in range(len(lat_rt)):
    for j in range(len(lon_rt)):
        if Ir_points[i,j] == 1:   
            ig = findig[i,j]
            ib = basins[i,j]
            if basin_clip[int(ig),int(ib)] == 1:
                if (np.isnan(ab_lat[i,j]) == 0):
                    loc0 = np.where(lat_index[:,0] == ab_lat[i,j])[0]
                    loc1 = np.where(lon_index[0,:] == ab_lon[i,j])[0]

                    if (len(loc0) == 1 and len(loc1) == 1):
                        ig = int(findig[loc0,loc1])
                        ib = int(basins[loc0,loc1])

                        if (basin_clip[ig, ib] == 1 and HTU_lsriver[ig,ib] == 1):
                            Ir_togrid[i,j]=ig
                            Ir_tobasin[i,j]=ib
                        if HTU_lsriver[ig,ib] == 1:
                            tmp[i,j] = 2
                            tmploc[i,j,0] = loc0
                            tmploc[i,j,1] = loc1
                        else:
                            tmp[i,j] = 1
                            tmploc[i,j,0] = loc0
                            tmploc[i,j,1] = loc1

tmp1 = np.copy(tmp)

tmp[np.isnan(tmp)]=0
tmp_area = np.nansum(tmp * basins_area)
tmp[tmp==2] = 1
#tmp[Ir_togrid >=0 ]=1
irrigation_area = np.nansum(tmp * basins_area)
print "total irrigation area over the domain" ,np.nansum( tmp*basins_area )/1000000,"km^2"
print "area which connected to the rivers", (tmp_area - irrigation_area)/1000000 , "km^2"
print "total basin area over the domain" ,np.sum( HTU_area ),"m^2"

if do_plotmap:
    pdf = PdfPages( figfolder + basin_name + '-pathway.pdf')
    IRRIGMAP(pdf, lat_rt_plot, lon_rt_plot, lat_rt, lon_rt,basins, findig,HTU_lsriver, COST,plot_dam_lon, plot_dam_lat,  tmp1, tmploc, fac, lat_int=0.5, lon_int=0.5)
    pdf.close()
irrigation_average = 24.37 * 1e8 / irrigation_area * 1000. / 365./ 86400.     # mm/s 
print "the average irrigation rate" , irrigation_average * 86400 * 365 , 'mm/yr'

#sys.exit()
# ==================================================================
#
# Part 3. Routing the runoff, the time axis is added
#
# ==================================================================

def ROUTING(year): 
    global restart_basin_id
    global nbdams
    kt = dt_write/ dt_routing
    # =========================================
    # 
    # try to find the restart file for reservoirs 
    restart_basin_file = outfolder+ 'restart_basin-' + basin_name + '-' + mode + str(year-1) + '.nc'  

    if restart_basin_id == 1:
        try:
            o = Dataset(restart_basin_file,'r')
            lat_restart = o.variables['lat'][:]
            lon_restart = o.variables['lon'][:]
            nbdams = len(o.dimensions['nbdams'])
            print 'number of nbdams ',nbdams
            # 
            if ( lat_restart[0] == lat[0] and lon_restart[0]==lon[0] and len(lat_restart) == len(lat) and len(lon_restart)==len(lon)):
                restart_basin_id = 1
                print 'Reservoirs restart file existed.'
            else:
                print 'The restart basin file is not corresponding to your setting'
                print 'The model will run, but the variables will be saved in file ', restart_basin_file
                restart_basin_id = 0
            o.close()
        except RuntimeError:
            print 'The restart file does not exist'
            print 'The model will run, but the variables will be saved in file ', restart_basin_file
            restart_basin_id=0
    # 
    # ttm: find whether the year is leap year, to determine the total time steps
    ttm = leap(year)

    # ================================
    # Setting of the hedging parameters for dam regulation
    # dams, eco_flows, flooding?
    # 
    alpha = np.zeros((nbpt,nbasmax,nbkt),'float64')     # If the release should be higher than the demand
    hedging = np.zeros((nbpt,nbasmax,nbkt),'float64')   # hedging factor
    TVref = np.zeros((ttm,nbdams),'float64')            # reference storage (vary in time)
    Vmax = np.zeros(nbdams,'float64')                   # the maximum capacity of a dam
    Vmin = np.zeros(nbdams,'float64')                   # the minimum capacity of a dam 
    #Vmin = np.zeros((nbpt,nbasmax),'float64')
    seconds = np.zeros(6,'f')                           # the seconds that divide different periods (in different periods, the regulation rules can be different
    lsdam = np.zeros((nbpt,nbasmax),'i')                # whether the HTU is a dam.
    lsdam[:]=-1
    # 
    if do_dams == 1: 
        print "there are", nbdams , "dams in the region"
        #for i in range(nbdams):
        for i in np.arange(0,nbdams):
            ig = lsdam_ig[i]
            ib = lsdam_ib[i]
            if dam_year[i] > year : 
                #print "dam", i, "is not considered since the built year", dam_year[i], "is after this year", year
                lsdam[ig,ib]=-1    
            else:
                lsdam[ig,ib]= i
                #
                Vmax[i] = dam_vref[i]
                # 
                if do_vmin == 1:
                    # the way to determine the minimum dam storage (dead capacity)
                    Vmin[i] = 0.356 * dam_vref[i]   # R^2 = 0.902
                else:
                    Vmin[i] = 0
                # 
                if do_vmax == 1:
                    kkt = kt_vmax
                else:
                    kkt = 1

                # flexiable maximum dam storage
                TVref[0,i] = dam_vref[i] 
                seconds[0] = 0
               
                t = (datetime.date(1900,flood_month[0],1) - datetime.date(1900,1,1)).days
                TVref[t,i] = dam_vref[i] 
                seconds[1] = t * 86400
     
                t = (datetime.date(1900,flood_month[1],1) - datetime.date(1900,1,1)).days
                TVref[t,i] = dam_vref[i] * kkt 
                seconds[2] = t * 86400
                
                t = (datetime.date(1900,flood_month[2],1) - datetime.date(1900,1,1)).days
                TVref[t,i] = dam_vref[i] * kkt
                seconds[3] = t * 86400
            
                t = (datetime.date(1900,flood_month[3],1) - datetime.date(1900,1,1)).days
                TVref[t,i] = dam_vref[i]  
                seconds[4] = t * 86400
        
                TVref[-1,i] = dam_vref[i]  
                seconds[-1] = t * 86400

                hedging[ig,ib,1] = hedging_eco
                hedging[ig,ib,3] = hedging_irrig
                alpha[ig,ib,1] = alpha_eco
                alpha[ig,ib,3] = alpha_irrig

                #print "dam", i, "is considered since the built year", dam_year[i], "before this year", year

    else:
        lsdam[:] = -1
        nbdams = 1
        Vmax = [0]          # to ensure that the variable can be passed to fortran subroutine
        Vmin = [0]
    print ">>>>>>>>> reading dams finished"
    print datetime.datetime.now()-start_time
    # ================================================
    #
    # Read the runoff / drainage files, restart files and define the output files 
    #  
    # ================================================
    # The file structure of the runoff/drainage are not the same in different forcings
    # For example: in Zun's simulation, all the variables are stored in a same file for the whole time period.
    # the name of the variables, and the units  are also different, which need to be specified here. 

    # Qs: runoff
    # Qsb: drainage 
    # Qirrig_req: irrigation requirement
    # Qirrig: irrigation estimated by ORCHIDEE

    if mode.find('CMA')>=0 or mode.find('E2O')>= 0 or mode.find('ITPCAS')>=0 or mode.find(EXPNAME)>=0 : 
        runoff_file_IN = runoff_IN_name + str(year)+'.nc'    
        runoff_file_IY = runoff_IY_name + str(year)+'.nc'    

        a1 = 0
        a2 = (datetime.date(year+1,01,01) - datetime.date(year,01,01)).days

        Qs_var = 'Qs'
        Qsb_var = 'Qsb'
        Qirrig_req_var = 'Qirrig_req'
        Qirrig_var = 'Qirrig'

        # units: kg/m2/s >> mm/s
        Qs_file= runoff_file_IN
        Qsb_file = runoff_file_IN
        Qs_in = clip.CLIPFILE_3(lat, lon, [a1,a2], Qs_file, Qs_var)
        Qsb_in = clip.CLIPFILE_3(lat, lon, [a1,a2], Qsb_file, Qsb_var)
        #DIS_in = clip.CLIPFILE_POINT (outlet_lat, outlet_lon, runoff_file_IN, 'Dis')

        Qirrig_file = runoff_file_IY
        Qirrig_req = clip.CLIPFILE_3(lat, lon, [a1,a2], Qirrig_file, Qirrig_req_var)
        Qirrig = clip.CLIPFILE_3(lat, lon, [a1,a2], Qirrig_file, Qirrig_var)
    
    elif mode.find('Zun') >=0 :
        yr_start = 1982
        a1 = (datetime.date(year,01,01) - datetime.date(yr_start,01,01)).days
        a2 = (datetime.date(year+1,01,01) - datetime.date(yr_start,01,01)).days

        runoff_file_IN = '/bdd/MEDI/workspaces/xzhou/ORC/Zun/NCRR4XD_GSWP3_ni_05deg_CN_D.nc'
        runoff_file_IY = '/bdd/MEDI/workspaces/xzhou/ORC/Zun/NCRR4XD_GSWP3_ir_05deg_CN_D.nc'
        #runoff_file_IN = '/bdd/MEDI/workspaces/xzhou/ORC/Zun/new/NCRR_GSWP3_ni_05deg_8214Y_CN_M.nc'
        #runoff_file_IY = '/bdd/MEDI/workspaces/xzhou/ORC/Zun/new/NCRR_GSWP3_irdp_05deg_8214Y_CN_M.nc'

        Qs_var = 'runoff_soil'
        Qsb_var = 'drainage_soil'
        Qirrig_req_var = 'transpot'
        Qirrig_var = 'irrig_fin'

        # units: mm/day >> mm/s
        Qs_file= runoff_file_IN
        Qsb_file = runoff_file_IN
        Qs_in_layer = clip.CLIPFILE_4(lat, lon, [a1,a2], Qs_file, Qs_var, full=1)
        lat1, lon1, Qs_in_layer = clip.CLIPFILE_4(lat, lon, [a1,a2], Qs_file, Qs_var, full=3)
        if len(lat1) != len(lat) or len(lon1)!=len(lon):
            print "the size of lat-lon doesn't correspond", len(lat1), len(lat), len(lon1), len(lon)
            print np.nanmin(lat1), np.nanmax(lat1), np.nanmin(lat), np.nanmax(lat)
            print np.nanmin(lon1), np.nanmax(lon1), np.nanmin(lon), np.nanmax(lon)
            sys.exit()
        Qsb_in_layer = clip.CLIPFILE_4(lat, lon, [a1,a2], Qsb_file, Qsb_var, full=1)
        #Qs_in = np.nansum(Qs_in_layer, axis=1)/86400
        #Qsb_in = np.nansum(Qsb_in_layer, axis=1)/86400
        Qs_in = np.nanmean(Qs_in_layer, axis=1)/86400
        Qsb_in = np.nanmean(Qsb_in_layer, axis=1)/86400
        #DIS_in = clip.CLIPFILE_POINT (outlet_lat, outlet_lon, runoff_file_IN, 'Dis')   

        if do_irrigation:
            Qirrig_file = runoff_file_IY
            Qirrig_req_layer = clip.CLIPFILE_4(lat, lon, [a1,a2], Qirrig_file, Qirrig_req_var,full=1)
            Qirrig_layer = clip.CLIPFILE_4(lat, lon, [a1,a2], Qirrig_file, Qirrig_var,full=1)
            #Qirrig_req = np.nansum(Qirrig_req_layer,axis=1)/86400
            #Qirrig = np.nansum(Qirrig_layer,axis=1)/86400
            Qirrig_req = np.nanmean(Qirrig_req_layer,axis=1)/86400
            Qirrig = np.nanmean(Qirrig_layer,axis=1)/86400
        else:
            Qirrig = np.copy(Qs_in)
            Qirrig_req = np.copy(Qs_in)
            Qirrig[:] =0
            Qirrig_req[:] =0

    Qs_in[Qs_in>100]=0
    Qsb_in[Qsb_in>100] =0
    Qs_in[Qs_in<0]=0
    Qsb_in[Qsb_in<0]=0
    Qs = daytorouting(Qs_in, dt_routing,  'mean')
    Qsb = daytorouting(Qsb_in, dt_routing, 'mean')

    # the irrigation deficit. Now it takes the value Qirrig_req as the deficit 
    # so that the model can estimate how much the requirement can be satisfied. 
    Qirrig_def_in = Qirrig_req
    if do_irrigation:
        #Qirrig_def_in = Qirrig_req-Qirrig
        Qirrig_def_in[Qirrig_def_in > 100] = 0
    else:
        Qirrig_def_in[:]=0

    # modify the values in daily to dt_routing 
    # if dt_routing < 1 day. unit :  mm/s 
    Qirrig_def = daytorouting(Qirrig_def_in, dt_routing, 'mean')
    Qirrig_def[Qirrig_def < 0.0] = 0

    # if we can use observed data to correct the Qs, Qsb and irrigation (assimilation)
    # default kt_discharge = 1, kt_irrig = 1
    if do_assimilation_discharge:
        Qs = Qs * kt_discharge
        Qsb = Qsb * kt_discharge
    if do_assimilation_irrigation:
        Qirrig_def = Qirrig_def * kt_irrig

    print ">>>>>>>>>finished reading runoff, irrigation "
    print datetime.datetime.now()-start_time
    #
    # ==========================================================================================
    #
    # variables stored for dams
    # 
    # ==========================================================================================
    if do_dams:
        dam_streamr = np.zeros((ttm * kt, nbdams, 2), 'float64')
        dam_discharge = np.zeros((ttm * kt, nbdams, 2), 'float64')
        dam_d = np.zeros((ttm * kt, nbdams, nbkt), 'float64')
        dam_Re = np.zeros((ttm * kt, nbdams, nbkt), 'float64')

        # possible maximum Vtarget, if the built year is later than the routing year, the value is 0. 
        TTVref = np.zeros((ttm*kt,nbdams),'float64')
        x1 = np.arange(ttm*kt) / float(kt)
        for i in range(nbdams):
            x = []
            y = []
            if dam_year[i] <= year : 
                for t in range(ttm):
                    if TVref[t,i] > 0:
                        x.append(t)
                        y.append(TVref[t,i])
                TTVref[:,i] = np.interp(x1, x, y)

    # Find the ecoflow settings, percentile of the historical flow in same HTU
    # setting the parameters for ecoflow 

    ecoflow = np.zeros((nbpt,nbasmax),'float64')
    if do_ecoflow: 
        o = Dataset( ecofile, 'r')
        try:
            streamr = o.variables['ecoflow'][:]     # kg/m2/s
        except KeyError:
            streamr = o.variables['streamr'][:]     # kg/m2/s
    
        for ib in range(nbasmax):
            if ecoflow_percentile > 0 :
                tmp = np.percentile(streamr[:,ib,:,:], 100-ecoflow_percentile, axis=0)
            else:
                tmp = np.nanmean(streamr[:,ib,:,:],axis=0)
            ecoflow[:,ib] = lalo_grid( tmp[:,:])  * total_area_grid * dt_routing
    
    # The exceedflow represents that the discharge/water level can not exceed a level.
    # It is in general a percentile of the historial flow 
    exceedflow = np.zeros((nbpt,nbasmax),'float64')
    if do_flood: 
        o = Dataset( floodfile, 'r')
        tmp = o.variables['discharge'][:]     # m3/s
        for ib in range(nbasmax):
            tmp_flood = np.percentile(tmp[:,ib,:,:], 100-flood_percentile, axis=0)
            exceedflow[:,ib] = lalo_grid( tmp_flood[:,:])    # m3/s
    #print "average of ecoflow criterion: " , np.nanmean(ecoflow)

    # =======================================================================
    # The other necessary model variables 
    # stream_in : the value in stream reservoir. It has five types: 0-natural, 1-ecoflow, 2-irrigation, 3-domestic, 4-hydropower
    # fast_in: value in fast reservoir, it has two types: 0-natural, 1-anthropogenic
    # slow_in: value in slwo reservoir, it has two types
    # e: extraction from stream reservoir
    # Fs: the outflow from stream reservoir

    stream_in = np.zeros((nbpt,nbasmax,nbkt),'float64')
    fast_in = np.zeros((nbpt,nbasmax,2),'float64')
    slow_in = np.zeros((nbpt,nbasmax,2),'float64')
    e = np.zeros((nbpt,nbasmax,nbkt),'float64')
    Fs = np.zeros((nbpt,nbasmax,nbkt),'float64')
    #
    if restart_basin_id == 1:
        # Read the water stored in routing reservoirs and water fluxes at the last time step. 
        print "Restarted file for reservoirs exit"
        print "    >> read the restarted routing reservoirs"
        o = Dataset(restart_basin_file,'r')
        fast_tmp=o.variables['fastres'][:]
        slow_tmp=o.variables['slowres'][:]
        stream_tmp=o.variables['streamres'][:]
        e_tmp = o.variables['e'][:]
        fs_tmp = o.variables['Fs'][:]
        o.close()
        #
        # change the shape
        for ib in range(nbasmax):
            for ik in range(2):
                fast_in[:,ib,ik] = lalo_grid(fast_tmp[ib,ik,:,:]) * total_area_grid 
                slow_in[:,ib,ik] = lalo_grid(slow_tmp[ib,ik,:,:]) * total_area_grid 
            for ik in range(nbkt):
                stream_in[:,ib,ik] = lalo_grid(stream_tmp[ib,ik,:,:]) * total_area_grid 
                e[:,ib,ik] = lalo_grid(e_tmp[ib,ik,:,:]) * total_area_grid * dt_routing
                Fs[:,ib,ik] = lalo_grid(fs_tmp[ib,ik,:,:]) * total_area_grid * dt_routing
        # print "using restarted reservoirs" , np.mean(stream_in)
    
    system_in_n = np.nansum(fast_in[:,:,0]) + np.nansum(slow_in[:,:,0])+ np.nansum(stream_in[:,:,0])
    system_in_g = np.nansum(fast_in[:,:,1]) + np.nansum(slow_in[:,:,1])+ np.nansum(stream_in[:,:,1:])
    print "Routing reservoir input",  system_in_n , system_in_g
    print datetime.datetime.now()-start_time
    
    ## ======================================================================
    # variables set for routing processes
    # igout, ibout, is the grid where the water flows out of the selected basins. 
    # if the basin has been cliped, the outlet of the cliped basin is recorded. 
    igout=[]
    ibout=[]
    for ig in range(nbpt):
        for ib in range(nbasmax):
            if basin_clip[ig,ib] >0 and route_togrid[ig,ib] <0 :
                igout.append(ig)
                ibout.append(ib)
    #
    #print "upstream ", upstream_htu[igout,ibout,0,:]
    #print "upstream ", upstream_htu[igout,ibout,1,:]
    #
    
    # d: the water demand for different types 
    # du: the unsatisfied demand 
    # Re: released water from dams (it is 0 for other HTUs) 
    d = np.zeros((nbpt,nbasmax,nbkt),'float64')
    du = np.zeros((nbpt,nbasmax,nbkt),'float64')
    Re = np.zeros((nbpt,nbasmax,nbkt),'float64') 

    # V: the water stored in the dam HTUs
    # lslake: whether if a HTU is a lake 
    # Vref_lake: the maximum volume for a lake HTU
    V = np.zeros((nbpt,nbasmax),'float64')
    lslake = np.zeros((nbpt,nbasmax),'float64')
    Vref_lake = np.zeros((nbpt,nbasmax),'float64')
    # the lake is not activated now
    lslake[:]= -1

    # The ratio of the time for propagation and routing
    TT = dt_propagate / dt_routing

    # du_propagate: the unsatisfied demand for propagation
    # flood_propagate: the space left for flooding downstream
    du_propagate = np.zeros((TT, nbpt, nbasmax, nbkt),'float64')
    flood_propagate = np.zeros((TT, nbpt, nbasmax,2),'float64')

    # total_flux_in: the total amount of water fluxes into the system
    # total_irrigation: the total irrigation requirement
    # total_irrigation_dl: the total satisfied water after regulation
    # total_irrigation_du: the total unsatisfied water after regulation.
    total_flux_in = np.zeros(ttm*kt, 'float64')
    #fs_out= np.zeros((ttm*kt,nbkt), 'float64')
    total_irrigation = np.zeros(ttm*kt, 'float64')
    total_irrigation_deficit = np.zeros(ttm*kt, 'float64')
    total_irrigation_actual = np.zeros(ttm*kt, 'float64')
    total_irrigation_dl = np.zeros(ttm*kt, 'float64')
    total_irrigation_du = np.zeros(ttm*kt, 'float64')
    T_irrig_demand = np.zeros(( ttm*kt, nbpt),'f')
    T_irrig_deficit = np.zeros(( ttm*kt, nbpt),'f')
    T_irrig_actual = np.zeros(( ttm*kt, nbpt),'f')
    
    QIRRIG_PIXEL_demand = np.zeros((ttm*kt, len(lat_rt), len(lon_rt)), 'f')
    QIRRIG_PIXEL_deficit = np.zeros((ttm*kt, len(lat_rt), len(lon_rt)), 'f')
    QIRRIG_PIXEL_add = np.zeros((ttm*kt, len(lat_rt), len(lon_rt)), 'f')
    # T_dl: the total water remand at a grid 
    # T_du: the total unsatisfied water at a grid
    # T_d: the total water demand at all downstream grid
    # T_e: the total water outflow from a grid
    T_dl= np.zeros(( ttm*kt, nbpt, nbkt),'f')
    T_du= np.zeros(( ttm*kt, nbpt, nbkt),'f')
    T_d= np.zeros(( ttm*kt, nbpt, nbkt),'f')
    T_e= np.zeros(( ttm*kt, nbpt, nbkt),'f')

    # ##### 
    # dis: the extra output at the points the users gave
    # the 10 values are: 
    # natural: 0-value at HTU, 2- maximum at grid, 4- minimum at grid, 6-std at grid, 
    #       8-sum of discharge flowing out the grid
    # anthropogenic : 1 , 3, 5, 7, 9
    dis = np.zeros((ttm*kt, nbpt,10),'float64')

    # T_streamr: the water in streamr reservoir (which is used for estimating the ecoflow)
    # T_discharge: the discharge at each HTU (which is used for estimating the ecoflow)
    if mode.find('nat')>=0:
        T_streamr = np.zeros((ttm*kt, nbpt, nbasmax, nbkt),'f')
        T_discharge = np.zeros((ttm*kt, nbpt, nbasmax, 2),'f')

    # output settings. Select the output cells which are used to identify different discharge of the grid
    dismask = np.zeros((nbpt,nbasmax),'f')
    dismask[:] = np.nan
    dismask_out = np.zeros((nbpt,nbasmax),'f')
    dismask_out[:] = np.nan
    tmp1 = np.zeros(nbasmax,'f')
    tmp2 = np.zeros(nbasmax,'f')
    for ig in range(nbpt):
        tmp = HTU_lsriver[ig,:]
        tmp1[:] = np.nan
        tmp2[:] = np.nan
        if np.nansum(tmp) > 0:
            # if there are rivers, take all the river HTU
            tmp1[tmp > 0] =1
            for ib in range(nbasmax):
                if tmp[ib] >0 and route_togrid[ig,ib] != ig :
                    tmp2[ib] = 1
        else:
            # if there are no rivers, take the largest 
            tmp3 = HTU_fac[ig,:]
            tmp1[tmp3==np.nanmax(tmp3)] = 1
            tmp2 = np.copy(tmp1)
        dismask[ig,:] = tmp1
        dismask_out[ig,:] = tmp2
    #
    # display the discharge at selected points
    # the grid (ig), basin (ib), latitude (lat), longitude(lon), and name (name) for the points that need to write.
    display_ig=[]
    display_ib=[]
    display_lat = []
    display_lon = []
    display_name = []
    loop = 0
    
    for i in range(len(display_coordinates_lat)):
        d_lat = display_coordinates_lat[i]
        d_lon = display_coordinates_lon[i]
        if (d_lat >= lat_rt[0] and d_lat <= lat_rt[-1] and d_lon >= lon_rt[0] and d_lon<= lon_rt[-1]):
            tmp1 = np.abs(lat_rt - d_lat)
            i_lat = np.where(tmp1 == np.min(tmp1))
            tmp1 = np.abs(lon_rt - d_lon)
            i_lon = np.where(tmp1 == np.min(tmp1))
            #
            ig = int(findig[i_lat, i_lon])
            ib = int(basins[i_lat, i_lon])
            
            while (loop < 10 and HTU_fac[ig,ib] < 1000):
                loop += 1
                igt = route_togrid[ig,ib]
                ibt = route_tobasin[ig,ib]
                ig = igt
                ib = ibt

            display_ig.append(ig)
            display_ib.append(ib)
            display_lat.append(d_lat)
            display_lon.append(d_lon)
            display_name.append(names[i])
    #
    print "The number of points to display", len(display_ig)
    print datetime.datetime.now()-start_time
    if len(display_ig) > 0:
        dis_display = np.zeros((ttm*kt, len(display_ig), 2),'f')
    #
    # ========================================
    #
    # !!! start the routing process
    #
    # ========================================
    print "!!! start the routing process >>>>>>>>>>>>>>>>>>> "
    print datetime.datetime.now()-start_time
    # 
    #for it in range(ttm*kt):
    # it_propagate: used to count the time steps for propagation. 
    # For example, if it reaches 23, the porpagation starts, and then it goes back to 0.
    it_propagate = 0

    DT = ttm * kt       # total time steps. Users can set a smaller value for testing the model.  
    #if do_testrouting:
    #    DT = 2000
    for it in range(DT):
        # The unit of runoff and drainage should be converted from kg/dt => kg
        runoff = lalo_grid(Qs[it,:,:]) * dt_routing
        drainage = lalo_grid(Qsb[it,:,:]) * dt_routing
        #
        if do_dams:
            Vtarget=TTVref[it,:]
        else:
            Vtarget = [0]
        #
        # step 1. calculate the irrigation netereq 
        # subroutine routing.distribute: because the irrigation deficit is given by grids (coarse resolution)
        # the subroutine distributes the irrigation demand to the irrigation points
        Qirrig_pixel_flux = routing.distribute(Qirrig_def[it,:,:], Ir_points)          # kg/m2/s, now the Qirrig_def equals the Qirrig_req
        Qirrig_pixel_demand = Qirrig_pixel_flux * basins_area * dt_routing * grid_clip       # convert to kg. 

        # subroutine routing.collect_htu: collect the irrigation demands to HTU
        irrig_needs = routing.collect_htu(nbpt, nbasmax, Qirrig_pixel_demand, findig, basins)          # kg
       
        # step 1. we update the reservoirs and outflows
        # Attention: the possible ourflow is estimated at the end of each time step
        # and then in the next time stpe, the outflow is modified with the remained water in reservoris and inflow. 
        # The modified outflow is stored for output. 
        stream_reservoir, fast_reservoir, slow_reservoir, discharge, flux_in, irrig_actual, tmp_fs  = routing.routing_flow_main(0, route_togrid, route_tobasin, \
                Vmax, Vmin, topo_resid , lsdam, lslake, stream_tcst, fast_tcst, slow_tcst, runoff, drainage, HTU_area, dt_routing, igout, ibout ,\
                e, stream_in, fast_in, slow_in, Fs, irrig_needs)
        #
        discharge = discharge / 1000. / dt_routing      # kg => m3/s
        #
        total_flux_in[it] = flux_in / dt_routing        # kg/s
        # 
        irrig_deficit = irrig_needs - irrig_actual      # kg

        # step 2. water demand and extraction 
        # step 2.1 collect the demand information for irrigation
        # distribute the irrig_deficit to points in proportion
        # This step estimates the unsatisfied irrigation demand which is used to propagate 
        # Qirrig_pixel_deficit(iim, jjm) : the deficit at pixels 
        Qirrig_pixel_deficit = routing.redistribute(findig, basins, irrig_deficit, irrig_needs, Qirrig_pixel_demand) 
        
        dl = np.zeros((nbpt,nbasmax,nbkt),'f')      # it is important to renew the dl to zeros at each time step. 
        dl[:,:,3] = routing.collect_demand( Qirrig_pixel_deficit, Ir_togrid, Ir_tobasin, HTU_lsriver)        # kg

        # water demand for ecoflow
        tmp = ecoflow - np.sum(stream_in[:,:,1:], axis=2)
        tmp[tmp<0]=0
        dl[:,:,1] = tmp 
        if do_ecoflow == 0:
            dl[:,:,1] = 0
        #
        #!!!! important, irrigation demand can use local reservoirs
        # You should know at each time, how much water is needed for each small grids
        #
        # Qirrig_def is at resolution 0.1 degree, while the model is at 1km
        # so we should create the 1km irrigation demand according to the 1km irrigation map


        # step 2.2 collect information on the river grids
        # Estimate the extration and unsatisfied water demand at HTUs.
        # The water demand is multiplied by TT, why?
        # because the propogation only conduct with a time interval of TT time steps. 
        # but here, is there any differences for ecoflow water and the others?
        du, e = routing.routing_waterdemand( dt_routing, HTU_lsriver, dl, TT * d, stream_reservoir)
        
        # the river extraction can be added to irrigation points 
        # the reverse process as the collect_demand. 
        Qirrig_pixel_add = routing.distribute_extract( dl[:,:,3], e[:,:,3], Ir_togrid, Ir_tobasin, HTU_lsriver, Qirrig_pixel_deficit  )

        QIRRIG_PIXEL_demand[it,:,:] = Qirrig_pixel_demand
        QIRRIG_PIXEL_deficit[it,:,:] = Qirrig_pixel_deficit
        QIRRIG_PIXEL_add[it,:,:] = Qirrig_pixel_add
        
        #total_irrigation[it] = np.sum(Qirrig_pixel_deficit)   # kg/s for write into a file
        total_irrigation[it] = np.sum(irrig_needs)/dt_routing   # kg/s for write into a file
        total_irrigation_deficit[it] = np.sum(irrig_deficit)/dt_routing   # kg/s for write into a file
        total_irrigation_actual[it] = np.sum(irrig_actual)/dt_routing   # kg/s for write into a file

        T_irrig_demand[it,:] = np.sum(irrig_needs,axis=1)/dt_routing
        T_irrig_deficit[it,:] = np.sum(irrig_deficit,axis=1)/dt_routing
        T_irrig_actual[it,:] = np.sum(irrig_actual,axis=1)/dt_routing

        #total_irrigation_dl[it] = np.sum(dl[:,:,3]) / dt_routing
        #total_irrigation_du[it] = np.sum(du[:,:,3]) / dt_routing

        T_dl[it,:,:] = np.sum(dl,axis=1)/dt_routing
        T_du[it,:,:] = np.sum(du,axis=1)/dt_routing
        T_d[it,:,:] = np.sum(d,axis=1)/dt_routing
        T_e[it,:,:] = np.sum(e,axis=1)/dt_routing

        du_propagate[it_propagate,:,:,:] = du
        flood_propagate[it_propagate,:,:,:] = discharge[:,:,:]
 
        # to overcome the difference between dt_routing and dt_propagate, the propagating process is isolated. 
        # simulate the outflow from dams. 
        it_propagate += 1
        if it_propagate == TT :  
            # The total unsatisfied water demand is the sum of the pervious TT time steps.
            du = np.sum(du_propagate,axis=0)
            #
            # define the regulation periods (0, 1, 2, 3)
            # 
            if do_dams:
                for i in range(4):
                    if it * dt_routing >= seconds[i] and it * dt_routing < seconds[i+1]: 
                       it_period = i + 1
            else:
                it_period = 0
            #
            # The current discharge is estimated as the mean values.
            dis_current = np.mean(flood_propagate, axis=0)  
            flood = dis_current[:,:,1] - exceedflow     # check if the discharge has exceeded the flooding control flow
            flood[flood<500.]=0.                            # flooding discharge with less than a criterion will be negelected. 
            #if (np.nanmax(flood)>0):  print "maximum of flood" , np.nanmax(flood)
            if do_flood ==0 :
                flood[:] = 0.
            #
            # subroutine: to propagate the water demands.
            tmp_d, tmp_Re =  routing.routing_propagate( 0, it_period, propagation_option , dt_propagate, upstream_htu, Vmax, Vmin, Vtarget,\
                     hedging, alpha, HTU_lsriver, lsdam,\
                     igout, ibout, stream_reservoir, du, route_togrid, route_tobasin, HTU_fac, dis_current , flood)
            
            # Attention: the water demand is divided to TT time steps
            # Attention: the estimated dam release is released at next TT time steps. 
            d = tmp_d / TT
            Re = tmp_Re / TT
            it_propagate = 0
            du_propagate[:]=0
            flood_propagate[:]=0
            #
        # Then we estimate the outflow from normal stream reservoir and lakes. 
        Fs =  routing.routing_smooth(stream_tcst, lake_tcst, dt_routing, topo_resid, lsdam, lslake, Vref_lake, stream_reservoir - e, Re) 
        # 
        # do another loop. 
        stream_in = np.copy(stream_reservoir)
        slow_in = np.copy(slow_reservoir)
        fast_in = np.copy(fast_reservoir)

        # be careful, the discharge at single grid is not necessary the outlet from the grid,
        # so there might be some problem when assessing the water balance using the discharge.
        # Thus there are many choices to select a value from the subbasins to represent the discharge value in this grid
        # As a old way, the largest value is chosen, 
        #  dis[it,:,:] = np.max(discharge, axis=1)
        # Now we can change the codes and store a few other values into the values (maximum, total flow out, average)
        # dismask / dismask_out

        for ik in range(2):
            tmp = discharge[:,:,ik] * dismask
            for ig in range(nbpt):
                tmp1 = tmp[ig,:]
                tmp1[tmp1<0.1*np.nanmax(tmp1)] = np.nan
                tmp[ig,:] = tmp1 
            dis[it,:,ik] = np.nanmax(tmp,axis=1)
            dis[it,:,ik+2] = np.nanmean(tmp,axis=1)
            dis[it,:,ik+4] = np.nanmin(tmp,axis=1)
            dis[it,:,ik+6] = np.nanstd(tmp,axis=1)
            tmp2 = discharge[:,:,ik] * dismask_out
            dis[it,:,ik+8] = np.nansum(tmp2,axis=1)
        #
        # Extract the discharge at given points
        if len(display_ig) > 0:
            for ik in range(len(display_ig)):
                ig = display_ig[ik]
                ib = display_ib[ik]
                dis_display[it,ik,:] = discharge[ig,ib,:]
        # 
        # prepare the variables for dams as output
        if do_dams : 
            for idam in range(nbdams):
                ig = lsdam_ig[idam]
                ib = lsdam_ib[idam]
                dam_streamr[it, idam, 0] = stream_reservoir[ig,ib,0]        # kg
                dam_streamr[it, idam, 1] = stream_reservoir[ig,ib,1]
                dam_discharge[it,idam,0] = discharge[ig,ib,0]  # m3/s
                dam_discharge[it,idam,1] = discharge[ig,ib,1] 
                dam_d[it,idam,:] = d[ig,ib,:]/1000./dt_routing           # kg 
                dam_Re[it,idam,:] = Re[ig,ib,:]/1000./dt_routing         # kg

        if mode.find('nat')>=0:
            T_discharge[it,:,:,:] = discharge 
            T_streamr[it,:,:,:] = stream_reservoir[:,:,:] 

        if do_testrouting or do_process:
            update_progress.update_progress( float(it)/(ttm*kt-1))

    # END of all calculating. 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #
    # Write the files 
    # 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # 
    krt = dt_write / dt_routing # the frequence of the model output
    dis_tmp = routingtowrite(dis, krt)

    dis_w = np.zeros((ttm, iim, jjm, 10),'float64')
    for it in range(ttm):
        for ik in range(10):
            dis_w[it,:,:,ik] = grid_lalo(dis_tmp[it,:,ik], iim,jjm)   

    if len(display_ig) > 0:
        dis_display_w = routingtowrite(dis_display, krt)

    if mode.find('nat')>=0:
        streamr_w = np.zeros((ttm, iim, jjm, nbasmax,nbkt), 'float64')
        tmp = routingtowrite(T_streamr, krt)
        T_dis_w = np.zeros((ttm, iim, jjm, nbasmax, 2),'float64')
        T_dis_tmp = routingtowrite(T_discharge, krt)
        for it in range(ttm):
            for ib in range(nbasmax):
                for ik in range(nbkt):
                    streamr_w[it,:,:,ib,ik] = grid_lalo(tmp[it,:,ib,ik], iim, jjm ) / total_area_lalo / dt_routing 
                for ik in range(2):
                    T_dis_w[it,:,:,ib,ik] = grid_lalo(T_dis_tmp[it,:,ib,ik], iim, jjm)    
    #
    #tmp_ecoflow = np.percentile(ecoflow, 100-ecoflow_percentile, axis=0)
    #ecoflow = np.zeros((iim,jjm,nbasmax),'f')
    #for ib in range(nbasmax):
    #    ecoflow[:,:,ib] = grid_lalo(tmp_ecoflow[:,ib]) / total_area_lalo / dt_routing 

    total_flux_in_w = routingtowrite(total_flux_in, krt)
    total_irrigation_w = routingtowrite(total_irrigation, krt)
    total_irrigation_deficit_w = routingtowrite(total_irrigation_deficit, krt)
    total_irrigation_actual_w = routingtowrite(total_irrigation_actual, krt)
    #total_irrigation_dl_w = routingtowrite(total_irrigation_dl, krt)
    #total_irrigation_du_w = routingtowrite(total_irrigation_du, krt)

    streamr = np.zeros((iim,jjm,nbasmax,nbkt),'float64')
    fastr=np.zeros((iim,jjm, nbasmax,2),'float64')
    slowr=np.zeros((iim,jjm, nbasmax,2),'float64')
    tmp_e = np.copy(e)
    e = np.zeros((iim,jjm,nbasmax,nbkt), 'float64')
    tmp_fs = np.copy(Fs)
    Fs = np.zeros((iim,jjm,nbasmax,nbkt), 'float64')
    for ib in range(nbasmax):
        for ik in range(2):
            fastr[:,:,ib,ik] = grid_lalo(fast_reservoir[:,ib,ik], iim, jjm) / total_area_lalo 
            slowr[:,:,ib,ik] = grid_lalo(slow_reservoir[:,ib,ik], iim, jjm) / total_area_lalo 
        for ik in range(nbkt):
            streamr[:,:,ib,ik] = grid_lalo(stream_reservoir[:,ib,ik], iim, jjm) / total_area_lalo
            e[:,:,ib,ik] = grid_lalo(tmp_e[:,ib,ik], iim, jjm) / total_area_lalo / dt_routing
            Fs[:,:,ib,ik] = grid_lalo(tmp_fs[:,ib,ik], iim, jjm) / total_area_lalo / dt_routing

    dl = np.zeros((ttm,iim,jjm,nbkt),'f')
    du = np.zeros((ttm,iim,jjm,nbkt),'f')
    d = np.zeros((ttm,iim,jjm,nbkt),'f')
    e_grid = np.zeros((ttm,iim,jjm,nbkt),'f')
    tmp0 = routingtowrite(T_dl, krt)
    tmp1 = routingtowrite(T_du, krt)
    tmp3 = routingtowrite(T_d, krt)
    tmp2 = routingtowrite(T_e, krt)
    for it in range(ttm):
        for ik in range(nbkt):
            dl[it,:,:,ik] = grid_lalo(tmp0[it,:,ik], iim, jjm) / total_area_lalo          # kg/m2/s
            du[it,:,:,ik] = grid_lalo(tmp1[it,:,ik], iim, jjm) / total_area_lalo
            d[it,:,:,ik] = grid_lalo(tmp3[it,:,ik], iim, jjm) / total_area_lalo
            e_grid[it,:,:,ik] = grid_lalo(tmp2[it,:,ik], iim, jjm) / total_area_lalo

    irrig_demand_grid = np.zeros((ttm,iim,jjm),'f')
    irrig_deficit_grid = np.zeros((ttm,iim,jjm),'f')
    irrig_actual_grid = np.zeros((ttm,iim,jjm),'f')
    tmp0 = routingtowrite(T_irrig_demand, krt)
    tmp1 = routingtowrite(T_irrig_deficit, krt)
    tmp2 = routingtowrite(T_irrig_actual, krt)
    for it in range(ttm):
        irrig_demand_grid[it,:,:] = grid_lalo(tmp0[it,:], iim, jjm) / total_area_lalo          # kg/m2/s
        irrig_deficit_grid[it,:,:] = grid_lalo(tmp1[it,:], iim, jjm) / total_area_lalo          # kg/m2/s
        irrig_actual_grid[it,:,:] = grid_lalo(tmp2[it,:], iim, jjm) / total_area_lalo          # kg/m2/s



    # irrigation at pixel > kg/m2/s
    qirrig_pixel_demand_w = np.zeros((ttm, len(lat_rt), len(lon_rt)),'f') 
    qirrig_pixel_deficit_w = np.zeros((ttm, len(lat_rt), len(lon_rt)),'f') 
    qirrig_pixel_add_w = np.zeros((ttm, len(lat_rt), len(lon_rt)),'f') 
    tmp0 = routingtowrite(QIRRIG_PIXEL_demand, krt) 
    tmp1 = routingtowrite(QIRRIG_PIXEL_deficit, krt) 
    tmp2 = routingtowrite(QIRRIG_PIXEL_add, krt)  
    for it in range(ttm):
        qirrig_pixel_demand_w[it, :,: ] = tmp0[it,:,:] / basins_area 
        qirrig_pixel_deficit_w[it, :,: ] = tmp1[it,:,:] / basins_area 
        qirrig_pixel_add_w[it, :,: ] = tmp2[it,:,:] / basins_area 
    # ===========================
    # store the variables 
    # the current unit for the variables are volumn (kg)
    # convert to kg/m2/dt and m3/s

    restart_basin_file_out = outfolder+ 'restart_basin-' + basin_name + '-' + mode + str(year) + '.nc' 
    o = Dataset(restart_basin_file_out, 'w')
    restart_basin_id = 1
    if write_history:
        o.history=history           # history:defined in the run_def, general description of the model settings.  
    o.createDimension('lat',iim)
    o.createDimension('lon',jjm)
    o.createDimension('lat_rt',len(lat_rt))
    o.createDimension('lon_rt',len(lon_rt))
    o.createDimension('time',ttm)
    o.createDimension('layer2',2)
    o.createDimension('layer10',10)
    o.createDimension('nbasmax',nbasmax)
    o.createDimension('nbkt',nbkt)
    if len(display_ig) > 0:
        o.createDimension('nb_display', len(display_ig))
        strlen = 25
        o.createDimension('str',strlen)
    
    o.createDimension('nbdams',nbdams)
    #
    vlat=o.createVariable('lat','f', ('lat',))
    vlat.units='degree_north'
    vlat.longname='latitude'
    vlat[:]=lat[:]
    # 
    vlon=o.createVariable('lon','f', ('lon',))
    vlon.units='degree_east'
    vlon.longname='longitude'
    vlon[:]=lon[:]
    #
    vlat=o.createVariable('lat_rt','f', ('lat_rt',))
    vlat.units='degree_north'
    vlat.longname='latitude'
    vlat[:]=lat_rt[:]
    # 
    vlon=o.createVariable('lon_rt','f', ('lon_rt',))
    vlon.units='degree_east'
    vlon.longname='longitude'
    vlon[:]=lon_rt[:]
    #
    vtime = o.createVariable('time','f',('time',))
    vtime.units= "seconds since "+ str(year) + "-01-01 00:00:00"
    for i in range(ttm):
        vtime[i] = 86400 * i
    #
    # Variables for restart file.
    v = o.createVariable('streamres','float64',('nbasmax','nbkt','lat','lon',))
    v.units='kg/m2'
    v.longname = 'stream reservoir'
    for ib in range(nbasmax):
        for ik in range(nbkt):
            v[ib,ik,:,:] = streamr[:,:,ib,ik]
    #
    v = o.createVariable('fastres','float64',('nbasmax','layer2','lat','lon',))
    v.units='kg/m2'
    v.longname='fast reservoir'
    for ib in range(nbasmax):
        for ik in range(2):
            v[ib,ik,:,:] = fastr[:,:,ib,ik]
    #
    v = o.createVariable('slowres','float64',('nbasmax','layer2','lat','lon',))
    v.units='kg/m2'
    v.longname='slow reservoir'
    v[:] = slowr[:]
    for ib in range(nbasmax):
        for ik in range(2):
            v[ib,ik,:,:] = slowr[:,:,ib,ik]
    #
    v = o.createVariable('e','float64',('nbasmax','nbkt','lat','lon',))
    v.units='kg/m2/s'
    v.longname = 'extraction water'
    for ib in range(nbasmax):
        for ik in range(nbkt):
            v[ib,ik,:,:] = e[:,:,ib,ik]
    #
    v = o.createVariable('Fs','float64',('nbasmax','nbkt','lat','lon',))
    v.units='kg/m2/s'
    v.longname = 'flow from the stream reservoir of last time step'
    for ib in range(nbasmax):
        for ik in range(nbkt):
            v[ib,ik,:,:] = Fs[:,:,ib,ik]
    #
    # Variables for the whole period. 
    v = o.createVariable('flux_in','float64',('time',))
    v.units='kg/s'
    v.longname = 'water flux in (runoff+drainage)'
    v[:] = total_flux_in_w[:]
    #
    v = o.createVariable('irrig','float64',('time',))
    v.units='kg/s'
    v.longname = 'local irrigation demand'
    v[:] = total_irrigation_w[:]
    #
    v = o.createVariable('irrig_deficit','float64',('time',))
    v.units='kg/s'
    v.longname = 'local irrigation deficit'
    v[:] = total_irrigation_deficit_w[:]
    #
    v = o.createVariable('irrig_actual','float64',('time',))
    v.units='kg/s'
    v.longname = 'local irrigation actual'
    v[:] = total_irrigation_actual_w[:]
    #
    #v = o.createVariable('irrig_dl','float64',('time',))
    #v.units='kg/s'
    #v.longname = 'irrigation demand'
    #v[:] = total_irrigation_dl_w[:]
    #
    #v = o.createVariable('irrig_du','float64',('time',))
    #v.units='kg/s'
    #v.longname = 'Unsatisfied irrigation demand'
    #v[:] = total_irrigation_du_w[:]
    #
    if mode.find('nat')>=0:
        v = o.createVariable('streamr', 'f',('time','nbasmax','nbkt','lat','lon',))
        v.units='kg/m2/s'
        v.longname = 'Stream routing water '
        for ib in range(nbasmax):
            for ik in range(nbkt):
                v[:,ib,ik,:,:] = streamr_w[:,:,:,ib,ik]
        #
        v = o.createVariable('discharge', 'f',('time','nbasmax','layer2','lat','lon',))
        v.units='m3/s'
        v.longname = 'Discharge'
        for ib in range(nbasmax):
            for ik in range(2):
                v[:,ib,ik,:,:] = T_dis_w[:,:,:,ib,ik]
    #
    v = o.createVariable('irrig_demand_pixel', 'f',('time','lat_rt','lon_rt',))
    v.units='kg/m2/s'
    v.longname = 'Irrigation water demand at pixels'
    v[:] = qirrig_pixel_demand_w[:,:,:]
    #
    v = o.createVariable('irrig_deficit_pixel', 'f',('time','lat_rt','lon_rt',))
    v.units='kg/m2/s'
    v.longname = 'Irrigation water deficit at pixels'
    v[:] = qirrig_pixel_deficit_w[:,:,:]
    #
    v = o.createVariable('irrig_add_pixel', 'f',('time','lat_rt','lon_rt',))
    v.units='kg/m2/s'
    v.longname = 'Irrigation extraction for rivers to pixels'
    v[:] = qirrig_pixel_add_w[:,:,:]
    #
    v = o.createVariable('irrig_demand_grid', 'f',('time','lat','lon',))
    v.units='kg/m2/s'
    v.longname = 'Irrigation water demand at grids'
    v[:] = irrig_demand_grid[:,:,:]
    #
    v = o.createVariable('irrig_deficit_grid', 'f',('time','lat','lon',))
    v.units='kg/m2/s'
    v.longname = 'Irrigation water deficit at grids'
    v[:] = irrig_deficit_grid[:,:,:]
    #
    v = o.createVariable('irrig_actual_grid', 'f',('time','lat','lon',))
    v.units='kg/m2/s'
    v.longname = 'Actual irrigation water at grids'
    v[:] = irrig_actual_grid[:,:,:]
    #
    v = o.createVariable('dl', 'f',('time','nbkt','lat','lon',))
    v.units='kg/m2/s'
    v.longname = 'Water requirement'
    for ik in range(nbkt):
        v[:,ik,:,:] = dl[:,:,:,ik]
    #
    v = o.createVariable('d', 'f',('time','nbkt','lat','lon',))
    v.units='kg/m2/s'
    v.longname = 'Unsatisfied downstream water requirement'
    for ik in range(nbkt):
        v[:,ik,:,:] = d[:,:,:,ik]
    #
    v = o.createVariable('du', 'f',('time','nbkt','lat','lon',))
    v.units='kg/m2/s'
    v.longname = 'Unsatisfied water requirement'
    for ik in range(nbkt):
        v[:,ik,:,:] = du[:,:,:,ik]
    #
    v = o.createVariable('e_grid', 'f',('time','nbkt','lat','lon',))
    v.units='kg/m2/s'
    v.longname = 'Extraction for water requirement'
    for ik in range(nbkt):
        v[:,ik,:,:] = e_grid[:,:,:,ik]
    #
    v = o.createVariable('dis','f',('time','layer10','lat','lon',))
    v.units='m3/s'
    v.longname = 'discharge, at max, mean, min, std, sum of outlets, 2*5'
    for ik in range(10):
        v[:,ik,:,:] = dis_w[:,:,:,ik]
    #
    #
    if len(display_ig) > 0:
        v = o.createVariable('display_dis','f',('time','nb_display','layer2'))
        v.units='m3/s'
        v.longname ='discharge at the exact points that given by users'
        v[:] = dis_display_w
        #
        v = o.createVariable('display_lat', 'f', 'nb_display')
        v.units='-'
        v.longname = 'latitude for the points that displayed'
        v[:] = display_lat[:]
        #
        v = o.createVariable('display_lon', 'f', 'nb_display')
        v.units='-'
        v.longname = 'longitude for the points that displayed'
        v[:] = display_lon[:]
    #
        v = o.createVariable('name','c',('nb_display','str'))
        v.long_name = 'station name'
        for i in range(len(display_ig)):
            v[i,:] = np.asarray(list(display_name[i].ljust(strlen)))
    #
    if do_dams:
        v = o.createVariable('dam_streamr', 'f' , ('time', 'nbdams', 'layer2',))
        v.units='kg'
        v.longname = 'Stream reservoir for dam HTU'
        v[:] = routingtowrite(dam_streamr,krt)[:]
        #
        v = o.createVariable('dam_discharge', 'f' , ('time', 'nbdams', 'layer2',))
        v.units='m^3/s'
        v.longname = 'Discharge from the dam HTU'
        v[:] = routingtowrite(dam_discharge,krt)[:]

        v = o.createVariable('dam_d', 'f' , ('time', 'nbdams', 'nbkt',))
        v.units='m^3/s'
        v.longname = 'Downstream unsatisfied demand'
        v[:] = routingtowrite(dam_d,krt)[:]

        v = o.createVariable('dam_Re', 'f' , ('time', 'nbdams', 'nbkt',))
        v.units='m^3/s'
        v.longname = 'Water release from the dam'
        v[:] = routingtowrite(dam_Re,krt)[:]

    o.close()

#if do_plotmap:
#    sys.exit()
print "Starting calculating routing, be patient......"
print datetime.datetime.now()-start_time
for year in range(styear, enyear+1):
    print "Doing the routing for year", year
    ROUTING(year)
    #sys.exit()



