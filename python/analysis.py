# Python 2.7 script to be exectued from the Python shell in a GRASS GIS
# mapset set to the coordinate system WGS84/UTM zone 32N (EPSG:32632). The Python
# module 'numpy' has to be installed and the 'v.centerpoint' add-on has to be
# installed to GRASS GIS.

# In the Python shell in GRASS GIS, set the working directory to the location of 
# this file and execute the script. This can be done by inserting the local file 
# path in the os.chdir() function that is called in the commented out code below, 
# and running the code in the shell:

# import os
# os.chdir([insert local filepath]/analysis.py')
# execfile('analysis.py')

import os
import sys
import math
import bisect
import numpy as np
from numpy import genfromtxt
from datetime import datetime

# Used to time the script.
startTime = datetime.now()

# Map layers required (these are imported below the definition of functions).
studyDtm = "studyarea" 
regionDtm = "region"
sites = "sites"
lakes =  "lakes"
sediments = "sediments" 

# Polygons splitting region along the border between bamble and porsgrunn. 
# These are used in the sampleByPhase() and emergence() functions, defined below.
splitSouth = "reg_split_south"
splitNorth = "reg_split_north"

# Path to precompiled shoreline displacement curves formatted for r.recode.
# Used with emergence().
recodeBamble = os.path.join(os.getcwd(),
                            '../gis_data/shoreline/recode_bamble.txt')
recodeGunnarsrod = os.path.join(os.getcwd(),
                                '../gis_data/shoreline/recode_gunnarsrod.txt')

# Specifies in how many directions fetch should be estimated.
# Used in computeFetch().
nInterval = 48

# Sets the r.viewshed parameters 'maximum viewing distance'
# and 'observerer elevation' for computeViewshed(), defined below.
maxDistance = '10000'
obsElev = '1.6'

# The size of the buffer to create around sites for sampling assumed non-sites.
bufferDist = '500'

# The number of random samples generated for each phase for
# each of the two sampling frames.
nSample = '1000'

########################### Defintion of functions ############################

# Function to create random samples. Samples are created seperately for each of
# the two areas covered by the two shoreline displacement curves and then
# combined. Elevation and slope range of sites also constrains the samples.
def sampleByPhase(phase, siteVect, bambleVect, gunnarsrodVect, dtmMap, nPts):

    # Makes sure the computational region is correct.
    grass.run_command('g.region', raster = dtmMap)

    # Use lakes as mask for constraining random samples.
    grass.run_command('r.mask', overwrite = True, flags = 'i', vector = lakes)

    # Find maximum slope value among sites in Bamble.
    slopeData = grass.read_command('v.db.univar', flags = 'g', map = siteVect,
                                   column = 'slope_average',
                                   where = " phase = '" + phase +
                                   "' AND municipality = 'bamb'")
    slopeDataByRow = slopeData.split('\n')
    slopeMax = slopeDataByRow[2].split('=')[1]

    # Create raster based on max slope in Bamble.
    slopeRangeBamble = 'slope_range_bamble' 
    grass.mapcalc(slopeRangeBamble + ' = if(' + slopeMap + ' <= ' + slopeMax +
                  ', 1, null())', overwrite = True)

    # Find range of elevation values among sites in Bamble.
    elevData = grass.read_command('v.db.univar', flags = 'g', map = siteVect,
                                  column = 'elev_minimum', 
                                  where = "phase = '" +
                                  phase + "' AND municipality = 'bamb'")
    elevDataByRow = elevData.split('\n')
    minElev = elevDataByRow[1].split('=')[1]
    maxElev = elevDataByRow[2].split('=')[1]

    # Limit elevation map by the bamble vector.
    grass.run_command('r.mask', overwrite = True, vector = bambleVect)
    
    # Create raster based on elevation range.
    elevRangeBamble = 'elev_range_bamble'
    grass.mapcalc(elevRangeBamble + ' = if(' + dtmMap + ' >= ' + minElev + 
                ' && ' + dtmMap + ' <= ' + maxElev + ', 1, null())', 
                overwrite = True)

    # Combine the null values defined by the elevation and slope ranges. 
    conRastBamble = 'con_rast_bamble'
    grass.mapcalc(conRastBamble + ' = ' + slopeRangeBamble + '/' +
                  elevRangeBamble, overwrite = True)

    # Create 500 m buffer around site polygons for constraining sampling.
    grass.run_command('v.buffer', input = siteVect, overwrite = True,
                      output = 'temp_buffer', distance = bufferDist,
                      where = "phase = '" + phase + "' AND municipality = 'bamb'")

    # Create convex hull around 500 m buffer for constraining sampling.
    tempHull = 'temp_hull_bamb_' + phase
    grass.run_command('v.hull', overwrite = True, input = 'temp_buffer',
                      output = tempHull)

    # Create sampling vector within constraining geometries using slope- and
    # elevation range rasters.
    conBufferBamble = 'temp_con_buffer_bamble'
    grass.run_command('r.mask', overwrite = True, vector = 'temp_buffer')
    grass.run_command('r.to.vect', overwrite = True, input = conRastBamble,
                      output = conBufferBamble, type = 'area')

    conHullBamble = 'temp_con_hull_bamble'
    grass.run_command('r.mask', overwrite = True, vector = tempHull)
    grass.run_command('r.to.vect', overwrite = True, input = conRastBamble,
                      output = conHullBamble, type = 'area')
    grass.run_command('r.mask', flags = 'r')

    # Draw samples within constrained areas.
    bufferSampleBamble = 'buffer_sample_bamble'
    grass.run_command('v.random', overwrite = True, output = bufferSampleBamble,
                      npoints = nPts, restrict = conBufferBamble)
    hullSampleBamble = 'hull_sample_bamble'
    grass.run_command('v.random', overwrite = True, output = hullSampleBamble,
                      npoints = nPts, restrict = conHullBamble)
    
    # Add tables to hold category number and name of sample type. 
    grass.run_command('v.db.addtable', map = bufferSampleBamble)
    grass.run_command('v.db.addtable', map = hullSampleBamble)
    
    # Add categories to each sample.
    grass.run_command('v.category', overwrite = True,
                      input = bufferSampleBamble,
                      output = 'temp_buffer_sample', option = 'add')
    grass.run_command('v.category', overwrite = True,
                      input = hullSampleBamble,
                      output = 'temp_hull_sample', option = 'add')
    
    # Add column to hold sample type ('buffer' or 'hull').
    grass.run_command('v.db.addcolumn', map = bufferSampleBamble,
                      columns = 'sample_type character(20)')
    grass.run_command('v.db.addcolumn', map = hullSampleBamble,
                      columns = 'sample_type character(20)')
    
    # Insert sample type.
    grass.run_command('v.db.update', map = bufferSampleBamble, layer = '1',
                      column = 'sample_type', value = 'buffer')
    grass.run_command('v.db.update', map = hullSampleBamble, layer = '1',
                      column = 'sample_type', value = 'hull')
    
    # Combine the two sample types.
    grass.run_command('v.patch', overwrite = True, flags = 'e',
                      input = bufferSampleBamble + ',' + hullSampleBamble,
                      output = 'bamble_sample')

    # Repeat for Gunnarsrod area.
    grass.run_command('r.mask', overwrite = True, flags = 'i', vector = lakes)

    slopeData = grass.read_command('v.db.univar', flags = 'g', map = siteVect,
                                   column = 'slope_average',
                                   where = " phase = '" + phase +
                                   "' AND municipality != 'bamb'")
    slopeDataByRow = slopeData.split('\n')
    slopeMin = slopeDataByRow[1].split('=')[1]
    slopeMax = slopeDataByRow[2].split('=')[1]

    slopeRangeGunnarsrod = 'slope_range_gunnarsrod' 
    grass.mapcalc(slopeRangeGunnarsrod + ' = if(' + slopeMap + ' <= ' +
                  slopeMax + ', 1, null())', overwrite = True)

    elevData = grass.read_command('v.db.univar', flags = 'g', map = siteVect,
                                  column = 'elev_minimum', where = "phase = '" +
                                  phase + "' AND municipality != 'bamb'")
    elevDataByRow = elevData.split('\n')
    minElev = elevDataByRow[1].split('=')[1]
    maxElev = elevDataByRow[2].split('=')[1]

    grass.run_command('r.mask', overwrite = True, vector = gunnarsrodVect)
    elevRangeGunnarsrod = 'elev_range_gunnarsrod'
    grass.mapcalc(elevRangeGunnarsrod + ' = if(' + dtmMap + ' >= ' + minElev +
                  ' && ' + dtmMap + ' <= ' + maxElev + ', 1, null())',
                  overwrite = True)

    conRastGunnarsrod = 'con_rast_gunnarsrod'
    grass.mapcalc(conRastGunnarsrod + ' = ' + slopeRangeGunnarsrod + '/' +
                  elevRangeGunnarsrod, overwrite = True)

    grass.run_command('v.buffer', input = siteVect, overwrite = True,
                      output = 'temp_buffer', distance = bufferDist,
                      where = " phase = '" + phase +
                      "' AND municipality != 'bamb'")
    tempHull = 'temp_hull_gunn_' + phase
    grass.run_command('v.hull', overwrite = True, input = 'temp_buffer',
                      output = 'tempHull')

    conBufferGunnarsrod = 'temp_con_buffer_gunnarsrod'
    grass.run_command('r.mask', overwrite = True, vector = 'temp_buffer')
    grass.run_command('r.to.vect', overwrite = True, input = conRastGunnarsrod,
                      output = conBufferGunnarsrod, type = 'area')

    conHullGunnarsrod = 'temp_constr_hull_gunnarsrod'
    grass.run_command('r.mask', overwrite = True, vector = 'tempHull')
    grass.run_command('r.to.vect', overwrite = True, input = conRastGunnarsrod,
                      output = conHullGunnarsrod, type = 'area')
    grass.run_command('r.mask', flags = 'r')

    bufferSampleGunnarsrod = 'buffer_sample_gunnarsrod'
    grass.run_command('v.random', overwrite = True,
                      output = bufferSampleGunnarsrod,
                      npoints = nPts, restrict = conBufferGunnarsrod)
    
    hullSampleGunnarsrod = 'hull_sample_gunnarsrod'
    grass.run_command('v.random', overwrite = True,
                      output = hullSampleGunnarsrod, npoints = nPts,
                      restrict = conHullGunnarsrod)
    
    grass.run_command('v.db.addtable', map = bufferSampleGunnarsrod)
    grass.run_command('v.db.addtable', map = hullSampleGunnarsrod)
    
    grass.run_command('v.category', overwrite = True,
                      input = bufferSampleGunnarsrod,
                      output = 'temp_buffer_sample', option = 'add')
    grass.run_command('v.category', overwrite = True,
                      input = hullSampleGunnarsrod,
                      output = 'temp_hull_sample', option = 'add')
    
    grass.run_command('v.db.addcolumn', map = bufferSampleGunnarsrod,
                      columns = 'sample_type character(20)')
    grass.run_command('v.db.addcolumn', map = hullSampleGunnarsrod,
                      columns = 'sample_type character(20)')
    
    grass.run_command('v.db.update', map = bufferSampleGunnarsrod, layer = '1',
                      column = 'sample_type', value = 'buffer')
    grass.run_command('v.db.update', map = hullSampleGunnarsrod, layer = '1',
                      column = 'sample_type', value = 'hull')
    
    grass.run_command('v.patch', overwrite = True, flags = 'e',
                      input = bufferSampleGunnarsrod + ',' +
                      hullSampleGunnarsrod, output = 'gunnarsrod_sample')

    # Combine the restricting polygons (used for visualisation).
    grass.run_command('v.patch', overwrite = True, flags = 'e',
                      input = conBufferBamble + ',' + conBufferGunnarsrod,
                      output = 'constr_buffer_' + phase)
    grass.run_command('v.patch', overwrite = True, flags = 'e',
                      input = conHullBamble + ',' + conHullGunnarsrod,
                      output = 'constr_hull_' + phase)

    # Combine the two sample vectors.
    sample = 'sample_' + phase
    grass.run_command('v.patch', overwrite = True, flags = 'e',
                      input = 'bamble_sample' + ',' + 'gunnarsrod_sample',
                      output = sample)

    # Add column holding what phase the sample belongs to.
    grass.run_command('v.db.addcolumn', map = sample,
                      columns = 'phase character(20)')
    grass.run_command('v.db.update', map = sample, layer = '1',
                      column = 'phase', value = phase)

    # Return name of sample vector.
    return(sample)

# Function to retrieve various data on the point features.
# To be used for fetch and viewshed.
def getPtData(pointFeature, minElev, centElevCol, xCoord, yCoord):
    pointData = grass.read_command("v.report", map = pointFeature,
                                   option = "coor")
    pointDataByLine = pointData.split('\n')
    nPoints = len(pointDataByLine) - 2
    fields = pointDataByLine[0].split('|')

    # Identify placement of columns of interest (cat is hardcoded as this is 
    # always present in a GRASS vector).
    i = 0
    for field in fields:
        if field == "cat":
            iCat = i
        if field == minElev:
            iMinElev = i
        if field == centElevCol:
            iCentElev = i
        if field == xCoord:
            iEasting = i
        if field == yCoord:
            iNorthing = i
        i = i + 1
    
    # Retrieve the data for each point for the columns identified above.
    pointData = []
    for point in range(1, nPoints + 1):
        fields = pointDataByLine[point].split('|')
        cat = fields[iCat]
        minElev = fields[iMinElev]
        centElevation = fields[iCentElev]
        easting = fields[iEasting]
        northing = fields[iNorthing]
        pointData.append([cat, minElev, centElevation, easting, northing])
        
    return(pointData)

# Function to retrieve data on polygon features. Follows the same structure as 
# getPtData() above, except minimum elevation is only retrieved if it is given
# as an argument. This is because the polygon cats need to retrieved before 
# min elevation has been found for the buffers of the random samples.
def getPolyData(polyFeature, minElev = None):
    polyData = grass.read_command("v.db.select", map = polyFeature)
    polyDataByLine = polyData.split('\n')
    nPoly = len(polyDataByLine) - 2
    fields = polyDataByLine[0].split('|')

    i = 0
    for field in fields:
        if field == "cat":
            iCat = i
        if minElev != None:
            if field == minElev:
                iMinElev = i
        i = i + 1
    
    polyData = []
    for poly in range(1, nPoly + 1):
        fields = polyDataByLine[poly].split('|')
        cat = fields[iCat]
        if minElev != None:
            minElev = fields[iMinElev]
            polyData.append([cat, minElev])
        else:
            polyData.append(cat)

    return(polyData)

# Utility function to identify sites consisting of multiple polygons, 
# create a centerpoint from these, delete the polygons in the old layer,
# and reinsert the point in their place. See: 
# https://gis.stackexchange.com/questions/330184/replacing-multiple-polygon
# -features-by-their-shared-center-point-in-grass-gis 
# Here this only pertains to Askeladden ID 236024. This was defined by the 
# archaeologists in the field as one site broken up by a quarry, hence the 
# multiple polygons.
def cleanMultiFeatures(vector, outputVect):
    cats1 = []
    cats2 = []

    # Retrieve cats by coordinates.
    # This does not report for cats with several polygons.
    cats = grass.read_command('v.report', map = vector, option = 'coor')
    catsByRow = cats.split('\n')
    for row in catsByRow:
        cats1.append(row.split('|')[0])

    # Retrieve cats by area.
    # This returns cats for multi-polygon features as well. 
    cats = grass.read_command('v.report', map = vector, option = 'area')
    catsByRow = cats.split('\n')
    for row in catsByRow:
        cats2.append(row.split('|')[0])

    # Identify what cats exist in cats2 but not in cats1
    s = set(cats1)
    multiCats = ([x for x in cats2 if x not in s])   

    # Loop over list of cats and create centerpoint of the centerpoints of the
    # polygons
    if len(multiCats) > 0:
        for i in multiCats:
            grass.run_command('v.extract', overwrite = True, input = vector,
                                cats = multiCats, output = 'multipoly')
            grass.run_command('v.centerpoint', overwrite = True, 
                                input = 'multipoly', output = 'polycenter')
            grass.run_command('v.centerpoint', overwrite = True, 
                                input = 'polycenter', 
                                output = 'centerpoint_' + str(i))
            grass.run_command('v.db.addtable', map = 'centerpoint_' + str(i))
            grass.run_command('db.copy', overwrite = True, from_table = vector,
                                to_table = 'centerpoint_' + str(i), 
                                where = "cat = " + str(i))
            grass.run_command('v.edit', map = vector, tool = 'delete',
                                 cats = i)
            grass.run_command('v.patch', overwrite = True, flags = 'e', 
                                input = vector + ',' + 'centerpoint_' + str(i),
                                 output = outputVect)
    return(outputVect)

# Function to retrieve the mode of sediment values beneath polygons. 
# Unfortunately, v.rast.stats does not provide the mode, so this was
# done by rasterising the polygons and using r.mode.
def sedimentMode(sedimentRast, attribute, polyFeature, polyData):

    # Set resolution to 1 for rasterisation of polygons
    grass.run_command('g.region', raster = sedimentRast, res= '1')  

    # Rasterise the polygons
    grass.run_command('v.to.rast', overwrite = True, input = polyFeature,
                      type = 'area', output = 'poly_rast', use = 'cat')

    # Convert the polygon raster to the mode of the underlying sediment raster
    grass.run_command('r.mode', overwrite = True, base = 'poly_rast',
                      cover = sedimentRast, output = 'temp_mode')

    # Retrieve the mode with v.rast.stats ('maximum' is abritrary as the
    # raster it is reading only consists of the mode beneath each site).
    grass.run_command('v.rast.stats', flags = 'c', map = polyFeature,
                      raster = 'temp_mode', column_prefix = attribute,
                      method = 'maximum')

    # Rename the column to reflect that the value is the mode
    grass.run_command('v.db.renamecolumn', map = polyFeature,
                  column = attribute + '_maximum,' + attribute + '_mode')

    # Here overlapping polygons are extracted and treated individually. 
    # See the main script for more on this.

    # Identify overlapping polygons (I did not find a straightforward
    # function for this)
    catData = grass.read_command('v.category', input = polyFeature,
                                 option = 'print')
    catDataByRow = catData.split('\n')

    # This for-loop tries to split each line in two. If that works then there
    # are overlapping polygons at that line. This follows from how v.category
    # returns categories above. 
    overlapCats = []
    for row in catDataByRow:
        try: 
            row.split('/')[1]
            overlapCats.append(row.split('/')[0])
            overlapCats.append(row.split('/')[1])
        except IndexError:
            continue

    # If there are overlapping polygons, loop over the category numbers        
    if overlapCats != []:
        for cat in overlapCats:
            
            # Extract overlapping polygon, creating a temporary feature
            grass.run_command('v.extract', flags = 'd', overwrite = True,
                            input = polyFeature, cats = cat,
                            output = 'temp_overlap')

            # Repeat the above steps for the overlapping polygon
            grass.run_command('v.to.rast', overwrite = True, 
                            input = 'temp_overlap', type = 'area', 
                            output = 'temp_poly_rast', use = 'cat')
            
            grass.run_command('r.mode', overwrite = True, 
                            base = 'temp_poly_rast', cover = sedimentRast,
                            output = 'temp_mode')

            grass.run_command('v.rast.stats', flags = 'c', 
                            map = 'temp_overlap', raster = 'temp_mode',
                             column_prefix = attribute, method = 'maximum')

            # Read the value from the temporary feature
            sedimentValue = grass.read_command('v.db.select', 
                                            map = 'temp_overlap',
                                             columns = attribute + '_maximum')

            # And upload it to the original polygon feature   
            grass.run_command('v.db.update', map = polyFeature, layer = '1',
                            column = attribute + '_mode',
                            value = str(sedimentValue.split('\n')[1]),
                            where = 'cat = ' + cat)

    # Reset resolution
    grass.run_command('g.region', raster = studyDtm)

# Function to compute viewshed size at a defined max and 5 km distance
# with sealevel raised to minimum elevation of the polygon (buffer/site polygon)
# associated with the analysed point.
def computeViewshed(pointFeature, pointData, dtmMap, maxDist, obsElev): 
    # Makes sure the computational region is correct
    grass.run_command('g.region', raster = dtmMap)

    # Assumes the input pointData to be looped has been retrieved with the 
    # getPtData function defined above
    i = 0
    for point in pointData:
        pointCat = str(point[0])
        pointCoords = str(point[3]) + ',' + str(point[4])
        
        # Raise sea level to match the lowest elevation of the point
        grass.mapcalc('temp_sealevel' + ' = if(' + dtmMap +' < ' +  
            str(point[1]) + ', 0, ' + studyDtm + ')', overwrite = True)
        
        #viewMap = pointCat + "_view" # Uncomment to keep each unique viewshed 
        viewMap = "temp_view"
        
        # Compute viewshed (c flag indicates that the curvature of the earth is
        # taken into account)
        grass.run_command("r.viewshed", overwrite = True, flags = 'bc',
                            input = 'temp_sealevel', output = viewMap,
                            coordinate = pointCoords,  max_dist = maxDist, 
                            observer_elevation = obsElev)
        
        # Rerieve number of visible cells from the viewshed map (at maxDist)
        viewcount = grass.read_command("r.stats", flags = "c", input = viewMap)
        viewcountByLine = viewcount.split("\n")
        viewsizeMax = viewcountByLine[1].split(" ")[1]

        # Create 5 km buffer to retrieve the view at 5 km
        grass.run_command('v.buffer', flags = 't', overwrite = True,
                            input = pointFeature, cats = pointCat, 
                            output = 'temp_view_buffer', distance = '5000')
        grass.run_command('v.rast.stats', map = 'temp_view_buffer',
                             raster = viewMap, column_prefix = 'view_5k',
                             method = 'sum')
        viewReport = grass.read_command('v.db.select', 
                                            map = 'temp_view_buffer',  
                                            columns = 'view_5k_sum',  
                                            where = 'cat = ' + pointCat)
        viewsize5k = viewReport.split('\n')[1]
        
        # Update the columns of the input pointFeature with the viewshed values
        # (these columns have to already exist in the attribute table)
        grass.run_command('v.db.update', map = pointFeature, layer = '1', 
                            column = 'viewsize_5k', value = viewsize5k, 
                            where = 'cat = ' + pointCat)
        grass.run_command('v.db.update', map = pointFeature, layer = '1',
                             column = 'viewsize_10k', value = viewsizeMax,
                             where = 'cat = ' + pointCat)
        i = i + 1

# Function to find the diagonal of a raster. Used to find the max extent for
# fetchlines to avoid edge-effects.
def diagonal(regionRaster):
    # Find extent of study region
    rasterInfo = grass.read_command('r.info', flags = 'g', map = regionRaster)
    infoByLine = rasterInfo.split('\n')

    # Retrieves the number of rows and columns, and the resolution of raster
    for line in infoByLine:
     if line.split('=')[0] == 'rows':
      rows = float(line.split('=')[1])
     if line.split('=')[0] == 'cols':
      cols = float(line.split('=')[1])
     if line.split('=')[0] == 'nsres':
      res = float(line.split('=')[1])
      
    # Finds the length of the diagonal in map units.
    diag = math.sqrt(rows**2 + cols**2)*res
    return(diag)

# Function to estimate average fetch from a point, 
# using two rasters (studyarea and larger region).       
def computeFetch(pointFeature, pointData, mapMax,
                 dtmMap, regionVect, nIntervals):
    
    # Finds degrees per interval to achieve desired number of lines.
    degInterval = 360.0 / nIntervals

    elev = []
    for point in pointData:
        pointCat = str(point[0])
        pointX = str(point[3])
        pointY = str(point[4])

        # Make sure that the resolution and comp region are correct
        # for sea level change.
        grass.run_command('g.region', vector = regionVect)
        
        # Raise the sealevel to the point, rounded up to make sure the point
        # ends up in the sea. This is necessary for the subsequent clipping
        # to work. If there already exists a sea vector at a point's elevation,
        # this is reused.
        seaVect = 'seavect_' + str(int(math.ceil(float(point[2]))))
        if math.ceil(float(point[2])) not in elev:
            # Create binary representation of sealevel.
            seaBin = 'sea_bin_' + str(int(math.ceil(float(point[2]))))
            grass.mapcalc(seaBin + ' = if(' + dtmMap +' < ' +
                             str(math.ceil(float(point[2]))) +
                             ', 1, null())', overwrite = True)

            # Vectorise the binary sea raster.
            grass.run_command('r.to.vect', overwrite = True, input = seaBin,
                                 output = seaVect, type = 'area')
            elev.append(math.ceil(float(point[2])))
        
        # Create buffer around point. 
        polyBuffer = "temp_buffer"
        lineBuffer = "temp_buffer_line"
        grass.run_command('v.buffer', overwrite = True, input = pointFeature,
                            cats = pointCat, output = polyBuffer,
                            distance = mapMax)
        
        # Adds table to the vector database and establishes
        # connection (required for next step).
        grass.run_command('v.db.addtable', map = polyBuffer)
        
        # Converts the border of the polygon buffer to a line.
        grass.run_command('v.to.lines', overwrite = True, input = polyBuffer,
                             output = lineBuffer)
     
        # Finds the length of the line buffer.
        bufferReport = grass.read_command('v.report', map = lineBuffer,
                                             option = 'length')
        bReportByLine = bufferReport.split('\n')
        bufferLength = float(bReportByLine[1].split('|')[1])

        # Finds the required distance between end-points of each radian.
        deg = bufferLength / 360
        distInterval = deg * degInterval

        # Creates a text document to hold the rules to be used in
        # v.segment below.
        dist = 0.0
        segRules = 'segment_rules.txt'
        f = open(segRules, 'w+')
        for i in range(1, nInterval + 1):
            f.write("P %d 1 %f\n" % (i, dist))
            dist = dist + distInterval
        f.close()

        # Create end-point every n degree along line buffer.
        circPts = "temp_segment"
        grass.run_command('v.segment', overwrite = True, input = lineBuffer,
                             output = circPts, rules = segRules)
     
        # Retrieve the coordinates for end-points.
        radialPtsReport = grass.read_command('v.report', map = circPts,
                                                option = 'coor')
        radialReportByLine = radialPtsReport.split('\n')
        nPts = len(radialReportByLine) - 2

        # Set up ascii text document to create radial lines
        # between central point and end-points.
        radialsTxt = 'radials_rules.txt'
        f = open(radialsTxt, 'w+')
        for pt in range(1, nPts + 1):
            xCoord = radialReportByLine[pt].split('|')[1]
            yCoord = radialReportByLine[pt].split('|')[2]
            f.write("L 2 1\n %s %s\n %s %s\n 1 %d\n" %
                 (pointX, pointY, xCoord, yCoord, pt))
        f.close()

        # Create radials from central point to end-points.
        radials = "temp_radials"
        grass.run_command('v.edit', overwrite = True, map = radials,
                            tool = 'create')
        grass.run_command('v.edit', flags = 'n', map = radials, tool = 'add',
                            input = radialsTxt)

        # Clip radials by the vector representing the 
        # sea in the study region.
        radialsClipped1 = 'temp_radials_clipped1'
        grass.run_command('v.clip', overwrite = True, input = radials,
                            clip = seaVect, output = radialsClipped1)
        
        # Clip radials by the vector representing the the modern day
        # coastline in the larger region.
        radialsClipped = "temp_radials_clipped"
        grass.run_command('v.clip', overwrite = True, input = radialsClipped1,
                            clip = regionVect, output = radialsClipped)
        
        # Retrieve lines that are intersecting the point after being clipped.
        # v.select seems to have an issue with rounding error, causing it to
        # miss some of the lines intersecting the point, see 
        # https://gis.stackexchange.com/questions/325823/retrieving-lines
        # -radiating-from-point-with-v-select-in-grass-gis-7-6
        fetchLines = 'fetch_lines_' + pointCat
        radialsIds = grass.read_command('v.edit', map = radialsClipped,
                                          type = 'line', tool = 'select', 
                                          coords = pointX + ',' + pointY)
        grass.run_command('v.edit', overwrite = True, map = fetchLines,
                            tool = 'create')
        grass.run_command('v.edit', map = fetchLines, type = 'line',
                            tool = 'copy', ids = radialsIds,
                            bgmap = radialsClipped)

        # Find the mean length of the fetch lines.
        fetchReport = grass.read_command('v.report', map = fetchLines,
                        option = 'length')
        fetchByLine = fetchReport.split('\n')
        nFetch = len(fetchByLine) - 2
        fetchLengths = []
        for line in range(1, nFetch + 1):
            fetchLengths.append(float(fetchByLine[line].split('|')[1]))
        avgFetch = sum(fetchLengths)/len(fetchLengths)

        # Update the point feature layer (column avg_fetch has to already exist).
        grass.run_command('v.db.update', map = pointFeature, layer = '1', 
                            column = 'avg_fetch', value = avgFetch, 
                            where = 'cat = ' + pointCat)
    
    # Reset computational region.
    grass.run_command('g.region', raster = studyDtm)

# Function to identify if point is on island (y/n), and if so of what size.
def location(pointFeature, pointData, dtmMap):

    # Create column to hold island size.
    grass.run_command('v.db.addcolumn', map = pointFeature,
                          columns = ['location character(20)',
                                     'island_size double precision'])

    # Make sure computational region is correct.
    grass.run_command('g.region', raster = dtmMap)
    
    # Create a vector outline of the computational region, to be used below.
    grass.run_command('v.in.region', output = 'dtm_edge', type = 'line',
                      overwrite = True)

    # Loop over points.
    for point in pointData:
        pointCat = str(point[0])

        # As v.what.vect cannot be restricted by SQL query or cats,
        # each point is extracted here.
        grass.run_command('v.extract', flags = 'd', overwrite = True,
                  input = pointFeature, cats = pointCat, output = 'temp_pt')
        
        # Raise the sealevel to the point.
        grass.mapcalc('temp_sealevel' + ' = if(' + dtmMap +' < ' +  
            str(math.floor(float(point[2]))) + ', null(), 1)', overwrite = True)

        # Vectorise the elevation raster.
        seaVect = 'temp_seavect_' + pointCat
        grass.run_command('r.to.vect', overwrite = True,
                          input = 'temp_sealevel',
                          output = seaVect, type = 'area')
        
        # Add a column to hold area and find area of each polygon in the vector.
        grass.run_command('v.db.addcolumn', map = seaVect,
                          columns = 'area double precision')
        grass.run_command('v.to.db', map = seaVect, column = 'area',
                          option = 'area')

        # Add column to hold polygon type and assign 'island' to all polygons.
        grass.run_command('v.db.addcolumn', map = seaVect,
                          columns = 'class character(20)')
        grass.run_command('v.db.update', map = seaVect,
                          column = 'class', value = 'island')
  
        # Retrieve polygons touching the edge of the computational region and
        # assign them the class 'mainland'. While this might wrongly classify
        # some islands as mainland, none of these stretch into the sampling
        # regions.
        grass.run_command('v.select', overwrite = True, ainput = seaVect,
                          binput= 'dtm_edge', output = 'temp_mainland',
                          operator = 'intersects')
        mainData = grass.read_command("v.db.select", map = 'temp_mainland')
        mainDataByLine = mainData.split('\n')

        # Retrieve categories for each mainland polygon.
        mainCats = []
        for poly in mainDataByLine[1:-1]:
            cat = poly.split('|')[0]
            mainCats.append(cat)

        # Code all polygons with these category numbers as mainland.
        grass.run_command('v.db.update', map = seaVect, column = 'class',
                          value = 'mainland', where = 'cat in ' + '(' +
                          str(mainCats).strip('[').strip(']') + ')')

        # Manual inspection at various sea levels indicates that no island is
        # ever larger than around 13.000 hectares, and that the smallest part 
        # of the mainland ever reaching into the sampling region is always 
        # smaller than around 50.000 hectares. All polygons larger than 20.000
        # were consequently coded 'mainland'.
        grass.run_command('v.db.update', map = seaVect, column = 'class',
                          where = 'area > 200000000', value = 'mainland')

        # Find whether point is on island or mainland.
        grass.run_command('v.what.vect', map = 'temp_pt',
                          column = 'location', query_map = seaVect,
                          query_column = 'class')

        # Retrieve location category from temporary point.
        islandLocation = grass.read_command('v.db.select', map = 'temp_pt',  
                                            columns = 'location')

        # Update point with the size of the polygon they are on.
        grass.run_command('v.what.vect', map = 'temp_pt',
                          column = 'island_size', query_map = seaVect,
                          query_column = 'area')

        # Retrieve polygon size from temporary point.
        islandSize = grass.read_command('v.db.select', map = 'temp_pt',  
                                            columns = 'island_size')

        # Update the original point feature with the values.
        grass.run_command('v.db.update', map = pointFeature, layer = '1', 
                          column = 'location',
                          value = islandLocation.split('\n')[1], 
                          where = 'cat = ' + pointCat)

        grass.run_command('v.db.update', map = pointFeature, layer = '1', 
                          column = 'island_size',
                          value = islandSize.split('\n')[1],
                          where = 'cat = ' + pointCat)
        
    # After the for-loop, replace island size with null for all points on
    # mainland.
    grass.run_command('v.db.update', map = pointFeature, column = 'island_size',
                          where = "location = 'mainland'", value = 'NULL')


# Function used to retrieve standard deviation in year of emergence of 
# shoreline at 3 different radii around site polygons, exluding areas 
# below sea level at the assumed time of occupation.
def emergence(polyFeature, polyData, emergeRaster, dtmMap):
        # Loop over polygons.
        for poly in polyData:
            # Skip polygons above 100 masl.
            if float(poly[1]) < 100:
                
                polyCat = poly[0]
                # In case of overlapping polygons in polyFeature, each feature
                # is extracted and treated individually.
                grass.run_command('v.extract', overwrite = True,
                  input = polyFeature, cats = polyCat, output = 'temp_feature')
                
                # Mask area below sea level.
                grass.mapcalc('temp_sealevel' + ' = if(' + dtmMap +' < ' +  
                                str(poly[1]) + ', 1, null())',
                                overwrite = True)
                grass.run_command('r.mask', overwrite = True, flags = 'i',
                                  raster = 'temp_sealevel')
                
                # Create buffer with distance 50 m.
                grass.run_command('v.buffer', flags = 't', overwrite = True,
                                  input = 'temp_feature',  
                                  output = 'temp_shore_buffer', distance = '50')
                
                # Retrieve standard deviation of year of emergence within buffer.
                grass.run_command('v.rast.stats', flags = 'c',
                                  map = 'temp_shore_buffer',
                                  raster = emergeRaster,
                                  column_prefix = 'emergence_50',
                                  method = 'stddev')

                # Retrieve the value.
                emerge50 = grass.read_command('v.db.select', 
                                              map = 'temp_shore_buffer',  
                                              columns = 'emergence_50_stddev',  
                                              where = 'cat = ' + polyCat)

                # Repeat for buffer with distance 500 m.
                grass.run_command('v.buffer', flags = 't', overwrite = True,
                                 input = 'temp_feature',  
                                 output = 'temp_shore_buffer', distance = '500')
            
                grass.run_command('v.rast.stats', flags = 'c',
                                  map = 'temp_shore_buffer',
                                  raster = emergeRaster,
                                  column_prefix = 'emergence_500',
                                  method = 'stddev')
                
                emerge500 = grass.read_command('v.db.select', 
                                               map = 'temp_shore_buffer',  
                                               columns = 'emergence_500_stddev',  
                                               where = 'cat = ' + polyCat)
                
                # Repeat for buffer with distance 1000 m.
                grass.run_command('v.buffer', flags = 't', overwrite = True,
                                input = 'temp_feature', 
                                output = 'temp_shore_buffer', distance = '1000')
            
                grass.run_command('v.rast.stats', flags = 'c',
                                  map = 'temp_shore_buffer',
                                  raster = emergeRaster,
                                  column_prefix = 'emergence_1k',
                                  method = 'stddev')

                emerge1k = grass.read_command('v.db.select', 
                                              map = 'temp_shore_buffer',  
                                              columns = 'emergence_1k_stddev',  
                                              where = 'cat = ' + polyCat)

                # Upload the above retrieved values to the attribute table of
                # the polygon feature. The columns need to have been predefined. 
                emergence50 = emerge50.split('\n')[1]
                grass.run_command('v.db.update', map = polyFeature, layer = '1',
                                  column = 'emergence_50',
                                  value = str(emergence50),
                                  where = 'cat = ' + polyCat)
                
                emergence500 = emerge500.split('\n')[1]
                grass.run_command('v.db.update', map = polyFeature, layer = '1',
                                  column = 'emergence_500',
                                  value = str(emergence500),
                                  where = 'cat = ' + polyCat)
                
                emergence1k = emerge1k.split('\n')[1]
                grass.run_command('v.db.update', map = polyFeature, layer = '1',
                                  column = 'emergence_1k',
                                  value = str(emergence1k),
                                  where = 'cat = ' + polyCat)

                # Remove mask, both for the repeated loop and for any
                # analysis to follow.
                grass.run_command('r.mask', flags = 'r')

################################# Main script ##################################
# Import the necessary layers to the mapset.
grass.run_command('v.import',
                  input = os.path.join(os.path.dirname(os.getcwd()),
                                    'gis_data', sites + '.gpkg'),
                                    overwrite = True)
grass.run_command('v.import',
                  input = os.path.join(os.path.dirname(os.getcwd()),
                                    'gis_data', lakes + '.gpkg'),
                                    overwrite = True)
grass.run_command('v.import',
                  input = os.path.join(os.path.dirname(os.getcwd()),
                                    'gis_data', sediments + '.gpkg'),
                                    overwrite = True)
grass.run_command('v.import',
                  input = os.path.join(os.path.dirname(os.getcwd()),
                                    'gis_data', splitSouth + '.gpkg'),
                                    overwrite = True)
grass.run_command('v.import',
                  input = os.path.join(os.path.dirname(os.getcwd()),
                                    'gis_data', splitNorth + '.gpkg'),
                                    overwrite = True)

# The digital terrain models (DTMs) are provided as multiple tiles for storage
# purposes. These therefore first have to be loaded in and then merged.

# Retrieve name of all studyarea DTM tiles (skip files not
# ending in .tif).
study_dtms = []
for file in os.listdir(os.path.join(os.path.dirname(os.getcwd()),
                                    'gis_data\studyarea')):
    if file.endswith(".tif"):
        study_dtms.append(file)

# Loop over and import each of these.
# study_names is to hold the name of each raster tile.
study_names = []
for i, n in enumerate(study_dtms, start=1):
    grass.run_command('r.import',
                      input = os.path.join(os.path.dirname(os.getcwd()),
                                    'gis_data\studyarea', n),
                      output = 'studyarea_' + str(i), overwrite = True)
    study_names.append('studyarea_' + str(i))

# Set computational region to all of the raster tiles.
grass.run_command('g.region', raster = ','.join(study_names))

# Combine these into one with r.patch, setting the output to
# studyDtm as defined at the start of the script.
grass.run_command('r.patch', input = ','.join(study_names),
                  output = studyDtm, overwrite = True)

# Repeat the above steps for the raster tiles of the larger
# region.
region_dtms = []
for file in os.listdir(os.path.join(os.path.dirname(os.getcwd()),
                                    'gis_data\\region')): # Note double escape
    if file.endswith(".tif"):
        region_dtms.append(file)

region_names = []
for i, n in enumerate(region_dtms, start=1):
    grass.run_command('r.import',
                      input = os.path.join(os.path.dirname(os.getcwd()),
                                    'gis_data\\region', n),
                      output = 'region_' + str(i), overwrite = True)
    region_names.append('region_' + str(i))

grass.run_command('g.region', raster = ','.join(region_names))
grass.run_command('r.patch', input = ','.join(region_names),
                  output = regionDtm, overwrite = True)

# Set computational region back to study area.
grass.run_command('g.region', raster = studyDtm)

# Finds the lowest elevation for each site.
grass.run_command('v.rast.stats', flags = 'c', map = sites, raster = studyDtm,
                    column_prefix = 'elev', method = 'minimum')

# A couple of the site polygons are small enough to not hit raster cell 
# centers, and are consequently not given a value with v.rast.stats. These
# are treated here, where the 'i' flag means that the value is interpolated
# from the four closest raster cells.
grass.run_command('v.what.rast', flags = 'i', map = sites, raster = studyDtm,
                  type = 'centroid', column = 'elev_minimum', 
                  where = 'elev_minimum ISNULL')

# Extract sites with adequate quality level for further 
# analysis (see the paper for more on this).
sitesTmp = 'sites_tmp'
quality = 3
grass.run_command('v.extract', overwrite = True,
                  input = sites, where = 'quality <= ' +
                  str(quality), output = sitesTmp)

# Add column to hold chronological phase.
grass.run_command('v.db.addcolumn', map = sitesTmp,
                    columns = ['phase varchar(15)'])

# Retrieve shoreline displacement curve for Bamble. This consists of two curves
# representing the confidence interval. Here this uncertainty is simply ignored
# and the mean of the two lines is used.
bamble = genfromtxt('../gis_data/shoreline/bamble_curve.csv',
                    delimiter = ',', names = True)
bamble['upperx'] = (bamble['upperx'] - 1950) * -1
bamble['lowerx'] = (bamble['lowerx'] - 1950) * -1

# Interpolate elevation values corresponding to each calendar year between
# -9500 and 1950 for each of the two curves (upper and lower).
calYears = range(-9500, 1950)
uppery = np.interp(calYears, bamble['upperx'], bamble['uppery'])
lowery = np.interp(calYears, bamble['lowerx'], bamble['lowery'])

# Retrieve the mean of the two curves.
meany = np.mean(np.array([[uppery], [lowery]]), axis = 0)

# Define look-up to ascribe sites to chronological phases.
# early mesolithic (em), middle mesolithic (mm), late mesolithic (lm),
# outside range (NULL).
phases = [(-8200, 'em'), (-6300, 'mm'), (-3600, 'lm'), (0, 'NULL')]

# Loop over sites from Bamble, find the mean shoreline date and ascribe site
# to corresponding phase.
elevDat = grass.read_command('v.db.select', map = sitesTmp,
                             where = "municipality == 'bamb'")
elevByRow = elevDat.split('\n')
for row in elevByRow[1:-1]:

    # Find shoreline date.
    sdate = np.interp(row.split('|')[4],
                      list(reversed(meany.reshape(-1).tolist())),
                      list(reversed(calYears)))
    
    # Check what phase this corresponds to in the look-up.
    phase = phases[bisect.bisect_right(phases, (sdate,))][1]

    # And update the record.
    grass.run_command('v.db.update', map = sitesTmp, layer = '1', 
                    column = 'phase', value = phase, 
                    where = 'cat = ' + row.split('|')[0])
    
# Repeat above steps for sites in Gunnarsrod area
# (Porsgrunn and Larvik municipalities).
gunnarsrod = genfromtxt('../gis_data/shoreline/gunnarsrod_curve.csv',
                        delimiter = ',', names = True)
gunnarsrod['upperx'] = (gunnarsrod['upperx'] - 1950) * -1
gunnarsrod['lowerx'] = (gunnarsrod['lowerx'] - 1950) * -1

calYears = range(-9500, 1950)
uppery = np.interp(calYears, gunnarsrod['upperx'], gunnarsrod['uppery'])
lowery = np.interp(calYears, gunnarsrod['lowerx'], gunnarsrod['lowery'])
meany = np.mean(np.array([[uppery], [lowery]]), axis = 0)

sitDat = grass.read_command('v.db.select', map = sitesTmp,
                            where = "municipality != 'bamb'")
sitByRow = sitDat.split('\n')
for row in sitByRow[1:-1]:
    sdate = np.interp(row.split('|')[4],
                      list(reversed(meany.reshape(-1).tolist())),
                      list(reversed(calYears)))
    phase = phases[bisect.bisect_right(phases, (sdate,))][1]
    grass.run_command('v.db.update', map = sitesTmp, layer = '1', 
                    column = 'phase', value = phase, 
                    where = 'cat = ' + row.split('|')[0])
    
# Exclude sites with a shoreline date outside the desired range
# (ascribed NULL to phase).
sitesQu = 'sites_q' + str(quality)
grass.run_command('v.extract', input = sitesTmp, overwrite = True,
                  where = "phase != 'NULL'", output = sitesQu)

# Create slope and aspect maps.
slopeMap = 'dtm_slope'
aspectMap = 'dtm_aspect'
grass.run_command('r.slope.aspect', flags = 'n', overwrite = True,
                  elevation = studyDtm, slope = slopeMap, aspect = aspectMap)

# Find average slope on sites.
grass.run_command('v.rast.stats', flags = 'c', map = sitesQu, raster = slopeMap,
                  column_prefix = 'slope', method = 'average')

# Find average slope for the remaining sites (site polygons not hitting center 
# of raster cells, see the step for finding elevation above).
grass.run_command('v.what.rast', flags = 'i', map = sitesQu, raster = slopeMap,
                  type = 'centroid', column = 'slope_average', 
                  where = 'slope_average ISNULL')

# Find average aspect on sites.
grass.run_command('v.rast.stats', flags = 'c', map = sitesQu,
                  raster = aspectMap, column_prefix = 'aspect',
                  method = 'average')

# Find average aspect for the remaining sites.
grass.run_command('v.what.rast', flags = 'i', map = sitesQu, raster = aspectMap,
                  type = 'centroid', column = 'aspect_average', 
                  where = 'aspect_average ISNULL')

# Rasterise the sediment polygons by infiltration level. Resolution of 1 seems good
# enough for both sediment and site polygons (see definition of sedimentMode()
# for more on this).
grass.run_command('g.region', vector = sediments, res= '1')
sedInfilRast = 'sed_infil_rast'
grass.run_command('v.to.rast', overwrite = True, input = sediments,
                 type = 'area', output = sedInfilRast, use = 'attr',
                  attribute_column = 'infilt')

# Reset computational region and resolution.
grass.run_command('g.region', raster = studyDtm)

# Add columns to sites that are to hold various data.
grass.run_command('v.db.addcolumn', map = sitesQu,
                    columns = ['infiltration_mode int',
                               'emergence_50 double precision',
                               'emergence_500 double precision',
                               'emergence_1k double precision'])

# Retrieve polygon data to be used for sedimentMode.
sitesPolyData = getPolyData(sitesQu, 'elev_minimum')

# Find the mode of the infiltration abilities of sediment types beneath sites.
sedimentMode(sedInfilRast, 'infiltration', sitesQu, sitesPolyData)

# Create shoreline emergence raster.
# First remove the present day sea from the raster.
dtmSea = 'dtm_sea'
grass.mapcalc(dtmSea + ' = if(' + studyDtm + '<= 0, null(),' + studyDtm + ')',
             overwrite = True)

# Create emergence map for the southern part of the study area,
# following the Bamble curve.
grass.run_command('r.mask', vector = splitSouth)
grass.run_command('r.recode', input = dtmSea, output = 'temp_south',
                 rules = recodeBamble, overwrite = True)

# Create emergence map for the northern part of the study area,
# following the Gunnarsrod curve.
grass.run_command('r.mask', overwrite = True, vector = splitNorth)
grass.run_command('r.recode', input = dtmSea, output = 'temp_north',
                 rules = recodeGunnarsrod, overwrite = True)

# Mask by present day lakes (where there is no depth data), and patch the two
# rasters. The i flag indicates inverse mask.
grass.run_command('r.mask', overwrite = True, flags = 'i', vector = lakes)
emergeMap = 'emerg'
grass.run_command('r.patch', overwrite = True, input = 'temp_north,temp_south',
                 output = emergeMap)

# Remove mask.
grass.run_command('r.mask', flags = 'r')

# Find stdv year of emergence for sites.
emergence(sitesQu, sitesPolyData, emergeMap, studyDtm)

# Create random sample of points for each phase.
# Number of samples is set at the top of the script.
sampleEm = sampleByPhase('em', sitesQu, splitSouth,
                         splitNorth,studyDtm, nSample)
sampleMm = sampleByPhase('mm', sitesQu, splitSouth,
                         splitNorth, studyDtm, nSample)
sampleLm = sampleByPhase('lm', sitesQu, splitSouth,
                         splitNorth, studyDtm, nSample)

# Combine these into one vector feature.
samplePts = 'sample_pts'
grass.run_command('v.patch', overwrite = True, flags = 'e', input = sampleEm +
                ',' + sampleMm + ',' + sampleLm, output = samplePts)

# Find the size of sites.
grass.run_command('v.db.addcolumn', map = sitesQu,
                  columns = ['site_size double precision'])
grass.run_command('v.to.db', map = sitesQu, option = 'area',
                 columns = 'site_size')

# Use median site size to determine the buffer size for random points.
sizeData = grass.read_command('v.univar', flags = 'ge', map = sitesQu,
                            column = 'site_size')
sizeDataByRow = sizeData.split('\n')
medianSize = sizeDataByRow[17].split('=')[1]
sampleDist = math.sqrt(float(medianSize)/math.pi)

sampleBuffers = samplePts + '_buffer'
grass.run_command('v.buffer', flags = 't', input = samplePts,
                 overwrite = True, output = sampleBuffers,
                 distance = sampleDist)
                            
# Add columns to sample buffers to hold data to be retrieved.
grass.run_command('v.db.addcolumn', map = sampleBuffers,
                   columns = ['elev_minimum double precision', 
                              'aspect_average double precision',
                              'infiltration_mode int',
                              'emergence_50 double precision',
                              'emergence_500 double precision',
                              'emergence_1k double precision'])

# Get sample buffer categories to be used in the below loop.
samplesPolyData = getPolyData(sampleBuffers)

# Loop over sample buffers to retrieve elevation and aspect data. 
# The topological vector format of GRASS (see v.import) means raster statistics 
# with overlapping buffers is less straightforward. As this can happen with 
# random points, all buffers are looped over and treated individually (the
# site data does not overlap).
for polyCat in samplesPolyData:
                            
      # Extract sample buffer.
      grass.run_command('v.extract', flags = 'd', overwrite = True,
                input = sampleBuffers, cats = polyCat,
                        output = 'temp_feature')
                            
      # Find the lowest elevation for sample buffer (to be used with viewshed).
      grass.run_command('v.rast.stats', flags = 'c', map = 'temp_feature',
                        raster = studyDtm, column_prefix = 'elev',
                        method = 'minimum')
                            
      # Find average aspect for random sample buffers.
      grass.run_command('v.rast.stats', flags = 'c', map = 'temp_feature',
                        raster = aspectMap, column_prefix = 'aspect',
                        method = 'average')
                            
      # Extract the values.                    
      values = grass.read_command('v.db.select', 
                         map = 'temp_feature',  
                         columns = 'elev_minimum,aspect_average',  
                         where = 'cat = ' + polyCat)
                            
      # Update the sampleBuffers feature with the elevation value.                 
      grass.run_command('v.db.update', map = sampleBuffers, layer = '1',
                                column = 'elev_minimum',
                                value = values.split('\n')[1].split('|')[0],
                                where = 'cat = ' + polyCat)

      # Update the sampleBuffers feature with the aspect value.                    
      grass.run_command('v.db.update', map = sampleBuffers, layer = '1',
                                column = 'aspect_average',
                                value = values.split('\n')[1].split('|')[1],
                                where = 'cat = ' + polyCat)

# Get sample buffer data, now with elevation data.
samplesPolyData = getPolyData(sampleBuffers, 'elev_minimum')

# As this was much slower, the sedimentMode() function identifies and only 
# treats overlapping polygons individually, unlike above where each one was
# extracted. Sediment infiltration levels for sample buffers.
sedimentMode(sedInfilRast, 'infiltration', sampleBuffers, samplesPolyData)

# Stdv year of emergence around sample buffers (this function is set to 
# extract each individual polygon).
emergence(sampleBuffers, samplesPolyData, emergeMap, studyDtm)

# Copy the data retrieved with the sample buffers back to the sample points.
grass.run_command('db.copy', overwrite = True, from_table = sampleBuffers, 
                to_table = samplePts)

# Finds the elevation of each sample point (to be used with functions below).
grass.run_command('v.what.rast', map = samplePts, raster = studyDtm,
                  type = 'point', column = 'cent_elev')

# Identify centerpoint of site constituted by multiple polygons
# (site with ID 236024, see definition of function) and replace 
# these with the single point.
sitesClean = cleanMultiFeatures(sitesQu, "sites_clean")

# Finds the average coordinate point for the rest of the sites.
siteCentPts = 'site_cent_pts'
grass.run_command('v.extract', overwrite = True, input = sitesClean,
                 output = 'temp_centroids')
grass.run_command('v.type', overwrite = True, input = 'temp_centroids',
                 output = siteCentPts, from_type = 'centroid', 
                 to_type = 'point')

# Finds the elevation of each site center point. 
grass.run_command('v.what.rast', map = siteCentPts, raster = studyDtm,
                  type = 'point', column = 'cent_elev')

# Create the different maps needed for fetch analysis.
# First set computational region to larger region, but retain resolution of 10 
# to not loose resolution in manipulation of study region raster.
grass.run_command('g.region', raster = regionDtm, res = '10')

# Create raster of larger region that holds null in the smaller study area.
grass.mapcalc('temp_zero = ' + studyDtm + '* 0', overwrite = True)
grass.mapcalc('temp_zero_null = if(isnull(temp_zero), 1, 0)', overwrite = True)
grass.mapcalc('temp_region_zero = ' + regionDtm + '/ temp_zero_null',
             overwrite = True)

# Creates an expanded version of the study area raster that holds zero in the
# larger region to allow fetch lines to extend beyond the study area before
# being clipped by the larger region. This is continously manipulated and 
# vectorised in computeFetch().
dtmExp = studyDtm + "_exp"
grass.mapcalc(dtmExp + ' = if(isnull(temp_region_zero),' + studyDtm + ', 0)',
             overwrite = True)
             
# Sets up the larger regional raster and subsequent regional sea vector to be
# empty in the study area to only clip fetch lines that escape the studyarea.
lregion = 'lregion'
grass.mapcalc(lregion + ' = if(isnull(temp_region_zero), 0,' + regionDtm + ')',
             overwrite = True)
grass.mapcalc('regional_sea = if(' + lregion +' <= 0, 1, null())',
             overwrite = True)
regSeaVect = 'reg_sea_vect'
grass.run_command('r.to.vect', overwrite = True, input = 'regional_sea',
                  output = regSeaVect, type = 'area')

# Reset the computational region.
grass.run_command('g.region', raster = studyDtm)

# Retrive the point data to be used for viewhsed and island analyses.
siteData = getPtData(siteCentPts, 'elev_minimum', 'cent_elev', 'x', 'y')
sampleData = getPtData(samplePts, 'elev_minimum', 'cent_elev', 'x', 'y')

# Find location (island/mainland) of site and sample points.
location(siteCentPts, siteData, studyDtm)
location(samplePts, sampleData, studyDtm)

# Add columns to hold fetch and viewshed data.
grass.run_command('v.db.addcolumn', map = siteCentPts,
                 columns = ['avg_fetch double precision',
                            'viewsize_5k double precision',
                            'viewsize_10k double precision'])
grass.run_command('v.db.addcolumn', map = samplePts,
                 columns = ['avg_fetch double precision',
                            'viewsize_5k double precision',
                            'viewsize_10k double precision'])

# Compute viewshed for site and sample points.
computeViewshed(siteCentPts, siteData, studyDtm, maxDistance, obsElev)
computeViewshed(samplePts, sampleData, studyDtm, maxDistance, obsElev)

# Find diagonal of regional dtm to find max fetch-line length.
maxDist = diagonal(regionDtm)

# Compute average fetch for site and sample points.
computeFetch(siteCentPts, siteData, maxDist, dtmExp, regSeaVect, nInterval)
computeFetch(samplePts, sampleData, maxDist, dtmExp, regSeaVect, nInterval)

# Export attribute tables as csv files.
grass.run_command('db.out.ogr', overwrite = True, input = siteCentPts,
                  output = '../gis_output/site_data.csv')
grass.run_command('db.out.ogr', overwrite = True, input = samplePts,
                  output = '../gis_output/sample_data.csv')

# Print the time it took to execute the script.
print datetime.now() - startTime